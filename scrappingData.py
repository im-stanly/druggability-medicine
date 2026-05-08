"""Search, filter, and download human protein structures from RCSB PDB.

Criteria:
- Human proteins (taxonomy id 9606)
- New releases (default: last 4 years)
- Resolution between 2-3 A (lower is better)
- R-free between 0.15-0.25 (lower is better)
- AlphaFold global pLDDT > 70 (higher is better)

Downloaded files:
- RCSB: .cif.gz
- AlphaFold: .cif (downloaded) -> .cif.gz (compressed locally)

Output layout:
- <output_dir>/<UNIPROT_ID>/PDB-<ENTRY_ID>.cif.gz
- <output_dir>/<UNIPROT_ID>/AF-*.cif.gz

Note: AlphaFold pLDDT is fetched from the AlphaFold DB API using UniProt IDs.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import time
from dataclasses import dataclass
from datetime import date, timedelta
from typing import Dict, Iterable, List, Optional, Tuple

import requests
from rcsbapi.search import Sort
from rcsbapi.search import search_attributes as attrs


RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{}"
RCSB_POLYMER_ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{}/{}"
RCSB_CIF_GZ_URL = "https://files.rcsb.org/download/{}.cif.gz"
ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{}"


@dataclass
class EntryRecord:
	entry_id: str
	release_date: Optional[str]
	resolution: Optional[float]
	r_free: Optional[float]
	uniprot_ids: List[str]
	alphafold_uniprot_id: Optional[str]
	alphafold_plddt: Optional[float]
	rcsb_cif_gz_path: Optional[str] = None
	alphafold_cif_gz_path: Optional[str] = None


def build_search_query(
	start_date: str,
	resolution_min: float,
	resolution_max: float,
	rfree_min: float,
	rfree_max: float,
) -> Tuple[object, bool]:
	filters: List[object] = []
	rfree_in_search = True

	# Experimental X-ray entries only; required for resolution and R-free metrics.
	filters.append(attrs.exptl.method == "X-RAY DIFFRACTION")
	filters.append(attrs.rcsb_accession_info.initial_release_date >= start_date)
	filters.append(attrs.rcsb_entry_info.resolution_combined >= resolution_min)
	filters.append(attrs.rcsb_entry_info.resolution_combined <= resolution_max)
	filters.append(attrs.refine.ls_R_factor_R_free >= rfree_min)
	filters.append(attrs.refine.ls_R_factor_R_free <= rfree_max)
	filters.append(
		attrs.rcsb_entity_source_organism.scientific_name == "Homo sapiens"
	)

	query = filters[0]
	for term in filters[1:]:
		query &= term
	return query, rfree_in_search


def search_entry_ids(
	start_date: str,
	resolution_min: float,
	resolution_max: float,
	rfree_min: float,
	rfree_max: float,
	limit: Optional[int] = None,
) -> Tuple[List[str], bool]:
	query, rfree_in_search = build_search_query(
		start_date, resolution_min, resolution_max, rfree_min, rfree_max
	)

	entry_ids = list(
		query(
			return_type="entry",
			sort=Sort(
				sort_by="rcsb_accession_info.initial_release_date",
				direction="desc",
			),
		)
	)

	if limit:
		entry_ids = entry_ids[:limit]

	return entry_ids, rfree_in_search


def safe_get(url: str, session: requests.Session, timeout: int = 30) -> Optional[dict]:
	response = session.get(url, timeout=timeout)
	if response.status_code == 404:
		return None
	response.raise_for_status()
	return response.json()


def extract_resolution(entry_json: dict) -> Optional[float]:
	values = (
		entry_json.get("rcsb_entry_info", {}).get("resolution_combined", [])
		if entry_json
		else []
	)
	if not values:
		return None
	try:
		return float(min(values))
	except (TypeError, ValueError):
		return None


def extract_r_free(entry_json: dict) -> Optional[float]:
	refine_blocks = entry_json.get("refine", []) if entry_json else []
	rfree_values = []
	for block in refine_blocks:
		if "ls_R_factor_R_free" in block and block["ls_R_factor_R_free"] is not None:
			rfree_values.append(block["ls_R_factor_R_free"])
	if not rfree_values:
		return None
	try:
		return float(min(rfree_values))
	except (TypeError, ValueError):
		return None


def extract_release_date(entry_json: dict) -> Optional[str]:
	return (
		entry_json.get("rcsb_accession_info", {}).get("initial_release_date")
		if entry_json
		else None
	)


def extract_polymer_entity_ids(entry_json: dict) -> List[str]:
	identifiers = entry_json.get("rcsb_entry_container_identifiers", {})
	return identifiers.get("polymer_entity_ids", []) if identifiers else []


def fetch_uniprot_ids(
	entry_id: str, entity_ids: Iterable[str], session: requests.Session
) -> List[str]:
	uniprot_ids: List[str] = []
	for entity_id in entity_ids:
		url = RCSB_POLYMER_ENTITY_URL.format(entry_id, entity_id)
		entity_json = safe_get(url, session)
		if not entity_json:
			continue
		if not is_human_entity(entity_json):
			continue
		identifiers = entity_json.get("rcsb_polymer_entity_container_identifiers", {})
		uniprot_ids.extend(identifiers.get("uniprot_ids", []) or [])
	# Deduplicate while preserving order.
	seen = set()
	deduped = []
	for uid in uniprot_ids:
		if uid not in seen:
			seen.add(uid)
			deduped.append(uid)
	return deduped


def is_human_entity(entity_json: dict) -> bool:
	def has_human_organism(items: Optional[List[dict]]) -> bool:
		if not items:
			return False
		for item in items:
			if item.get("ncbi_taxonomy_id") == 9606:
				return True
			if item.get("scientific_name") == "Homo sapiens":
				return True
		return False

	return has_human_organism(entity_json.get("rcsb_entity_source_organism")) or has_human_organism(
		entity_json.get("rcsb_entity_host_organism")
	)


def fetch_alphafold_info(
	uniprot_id: str,
	session: requests.Session,
	cache: Dict[str, Tuple[Optional[float], Optional[str]]],
) -> Tuple[Optional[float], Optional[str]]:
	if uniprot_id in cache:
		return cache[uniprot_id]

	url = ALPHAFOLD_API_URL.format(uniprot_id)
	response = session.get(url, timeout=30)
	if response.status_code == 404:
		cache[uniprot_id] = (None, None)
		return cache[uniprot_id]
	response.raise_for_status()
	payload = response.json()
	if not payload:
		cache[uniprot_id] = (None, None)
		return cache[uniprot_id]
	# AlphaFold DB returns a list; use the first item.
	value = payload[0].get("globalMetricValue")
	cif_url = payload[0].get("cifUrl")
	plddt = float(value) if value is not None else None
	cache[uniprot_id] = (plddt, cif_url)
	return cache[uniprot_id]


def gzip_file(path: str) -> str:
	if path.endswith(".gz"):
		return path
	gz_path = f"{path}.gz"
	if os.path.exists(gz_path):
		return gz_path
	with open(path, "rb") as source, gzip.open(gz_path, "wb") as target:
		shutil.copyfileobj(source, target)
	os.remove(path)
	return gz_path


def download_rcsb_cif_gz(entry_id: str, output_dir: str, session: requests.Session) -> str:
	os.makedirs(output_dir, exist_ok=True)
	url = RCSB_CIF_GZ_URL.format(entry_id)
	path = os.path.join(output_dir, f"PDB-{entry_id}.cif.gz")
	if os.path.exists(path):
		return path
	with session.get(url, stream=True, timeout=60) as response:
		response.raise_for_status()
		with open(path, "wb") as handle:
			for chunk in response.iter_content(chunk_size=1024 * 64):
				if chunk:
					handle.write(chunk)
	return path


def download_alphafold_cif_gz(
	uniprot_id: str, cif_url: Optional[str], output_dir: str, session: requests.Session
) -> Optional[str]:
	if not cif_url:
		return None
	os.makedirs(output_dir, exist_ok=True)
	filename = os.path.basename(cif_url)
	if not filename.endswith(".cif"):
		filename = f"AF-{uniprot_id}.cif"
	path = os.path.join(output_dir, filename)
	gz_path = f"{path}.gz"
	if os.path.exists(gz_path):
		return gz_path
	with session.get(cif_url, stream=True, timeout=60) as response:
		response.raise_for_status()
		with open(path, "wb") as handle:
			for chunk in response.iter_content(chunk_size=1024 * 64):
				if chunk:
					handle.write(chunk)
	return gzip_file(path)


def collect_records(
	entry_ids: Iterable[str],
	rfree_min: float,
	rfree_max: float,
	plddt_min: float,
	download: bool,
	output_dir: str,
	polite_delay: float,
	debug: bool,
) -> List[EntryRecord]:
	session = requests.Session()
	alphafold_cache: Dict[str, Tuple[Optional[float], Optional[str]]] = {}
	records: List[EntryRecord] = []
	stats = {
		"entries_total": 0,
		"entry_json_missing": 0,
		"rfree_missing": 0,
		"rfree_out_of_range": 0,
		"no_entity_ids": 0,
		"no_human_uniprot": 0,
		"no_alphafold_model": 0,
		"plddt_below_threshold": 0,
		"records_saved": 0,
	}

	for entry_id in entry_ids:
		stats["entries_total"] += 1
		entry_json = safe_get(RCSB_ENTRY_URL.format(entry_id), session)
		if not entry_json:
			stats["entry_json_missing"] += 1
			continue

		r_free = extract_r_free(entry_json)
		if r_free is None:
			stats["rfree_missing"] += 1
			continue
		if not (rfree_min <= r_free <= rfree_max):
			stats["rfree_out_of_range"] += 1
			continue

		entity_ids = extract_polymer_entity_ids(entry_json)
		if not entity_ids:
			stats["no_entity_ids"] += 1
			continue
		uniprot_ids = fetch_uniprot_ids(entry_id, entity_ids, session)
		if not uniprot_ids:
			stats["no_human_uniprot"] += 1
			continue

		best_plddt: Optional[float] = None
		best_uid: Optional[str] = None
		best_cif_url: Optional[str] = None
		for uid in uniprot_ids:
			plddt, cif_url = fetch_alphafold_info(uid, session, alphafold_cache)
			if plddt is None or cif_url is None:
				continue
			if best_plddt is None or plddt > best_plddt:
				best_plddt = plddt
				best_uid = uid
				best_cif_url = cif_url

		if best_plddt is None:
			stats["no_alphafold_model"] += 1
			continue
		if best_plddt < plddt_min:
			stats["plddt_below_threshold"] += 1
			continue

		record = EntryRecord(
			entry_id=entry_id,
			release_date=extract_release_date(entry_json),
			resolution=extract_resolution(entry_json),
			r_free=r_free,
			uniprot_ids=uniprot_ids,
			alphafold_uniprot_id=best_uid,
			alphafold_plddt=best_plddt,
		)

		if download:
			protein_dir = os.path.join(output_dir, best_uid or entry_id)
			record.rcsb_cif_gz_path = download_rcsb_cif_gz(entry_id, protein_dir, session)
			record.alphafold_cif_gz_path = download_alphafold_cif_gz(
				best_uid or "", best_cif_url, protein_dir, session
			)

		records.append(record)
		stats["records_saved"] += 1
		if polite_delay > 0:
			time.sleep(polite_delay)

	if debug:
		print("Debug stats:")
		for key, value in stats.items():
			print(f"- {key}: {value}")

	return records


def write_output(records: List[EntryRecord], output_path: str) -> None:
	payload = [record.__dict__ for record in records]
	with open(output_path, "w", encoding="utf-8") as handle:
		json.dump(payload, handle, indent=2)


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description="Search and download recent human protein structures from RCSB and AlphaFold."
	)
	parser.add_argument("--years", type=int, default=4)
	parser.add_argument("--resolution-min", type=float, default=2.0)
	parser.add_argument("--resolution-max", type=float, default=3.0)
	parser.add_argument("--rfree-min", type=float, default=0.15)
	parser.add_argument("--rfree-max", type=float, default=0.25)
	parser.add_argument("--plddt-min", type=float, default=70.0)
	parser.add_argument("--limit", type=int, default=500)
	parser.add_argument("--no-download", dest="download", action="store_false")
	parser.add_argument("--output-dir", default="data/01_raw/proteins")
	parser.add_argument("--output-json", default="data/01_raw/rcsb_hits.json")
	parser.add_argument("--polite-delay", type=float, default=0.25)
	parser.add_argument("--debug", action="store_true")
	parser.set_defaults(download=True)
	return parser.parse_args()


def main() -> None:
	args = parse_args()
	start_date = (date.today() - timedelta(days=args.years * 365)).isoformat()

	entry_ids, _ = search_entry_ids(
		start_date,
		args.resolution_min,
		args.resolution_max,
		args.rfree_min,
		args.rfree_max,
		args.limit,
	)
	if args.debug:
		print(f"Search returned {len(entry_ids)} entry IDs")

	if not entry_ids:
		print("No entries found by RCSB search criteria.")
		return

	records = collect_records(
		entry_ids=entry_ids,
		rfree_min=args.rfree_min,
		rfree_max=args.rfree_max,
		plddt_min=args.plddt_min,
		download=args.download,
		output_dir=args.output_dir,
		polite_delay=args.polite_delay,
		debug=args.debug,
	)

	# Rank: resolution asc, R-free asc, AlphaFold pLDDT desc
	records.sort(
		key=lambda r: (
			r.resolution if r.resolution is not None else 99.0,
			r.r_free if r.r_free is not None else 99.0,
			-(r.alphafold_plddt if r.alphafold_plddt is not None else 0.0),
		)
	)

	os.makedirs(os.path.dirname(args.output_json), exist_ok=True)
	write_output(records, args.output_json)

	print(f"Saved {len(records)} hits to {args.output_json}")


if __name__ == "__main__":
	main()
