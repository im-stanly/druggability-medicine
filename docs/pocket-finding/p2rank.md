# P2RANK — zasada działania (z claude xd)

## Czym jest P2RANK?
[Krivák, R., Hoksza, D. P2Rank: machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure. J Cheminform 10, 39 (2018).](https://doi.org/10.1186/s13321-018-0285-8)

P2RANK (Protein to Rank) to narzędzie bioinformatyczne oparte na uczeniu maszynowym, służące do **przewidywania miejsc wiązania ligandów** w strukturach białek. Jego głównym zastosowaniem jest identyfikacja potencjalnych kieszeni wiążących ligandy na powierzchni białka bez konieczności korzystania z zewnętrznych baz danych, szablonów strukturalnych ani oprogramowania trzeciego.

---

## Etapy działania algorytmu

### 1. Generowanie punktów powierzchniowych

P2RANK jako dane wejściowe przyjmuje trójwymiarową strukturę białka (plik PDB). Na jej podstawie wyznaczana jest **powierzchnia dostępna dla rozpuszczalnika** (ang. *solvent accessible surface*, SAS). Na tej powierzchni rozmieszczane są równomiernie punkty próbkowania — tzw. *near-surface points* — które stanowią podstawowe jednostki analizy.

### 2. Reprezentacja lokalnego otoczenia chemicznego

Dla każdego punktu powierzchniowego algorytm analizuje jego **lokalne chemiczne sąsiedztwo** — zbiór atomów białka leżących w promieniu punktu. Cechy fizykochemiczne tych atomów (np. hydrofobowość, polarność, ładunek, objętość van der Waalsa) są rzutowane na punkt i agregowane w jeden wektor cech opisujący dane mikrośrodowisko.

### 3. Klasyfikacja przy użyciu Random Forest

Wektor cech każdego punktu jest podawany na wejście **klasyfikatora Random Forest** — zespołu drzew decyzyjnych. Model, wytrenowany wcześniej na znanych strukturach białko-ligand, przypisuje każdemu punktowi wartość **ligandability** — prawdopodobieństwo, że dany punkt leży w pobliżu rzeczywistego miejsca wiązania liganda.

### 4. Klastrowanie punktów

Punkty o wysokim przewidywanym *ligandability* są grupowane metodą **klastrowania** (np. hierarchicznego lub opartego na gęstości). Każdy klaster odpowiada potencjalnemu miejscu wiązania na powierzchni białka.

### 5. Ranking przewidywanych kieszeni

Zidentyfikowane klastry są **rangowane** według zagregowanej wartości *ligandability* ich punktów. W wyniku działania algorytmu użytkownik otrzymuje listę kandydatów miejsc wiążących, uszeregowaną od najbardziej do najmniej prawdopodobnych.

---

## Kluczowe cechy P2RANK

| Cecha | Opis |
|---|---|
| **Brak szablonów** | Nie korzysta z homologicznych struktur ani baz sekwencji |
| **Szybkość** | Analiza pojedynczego białka trwa poniżej 1 sekundy |
| **Niezależność** | Nie wymaga zewnętrznych narzędzi (np. NACCESS, Fpocket) |
| **Skalowalność** | Obsługuje wielowątkowe przetwarzanie dużych zbiorów danych |
| **Otwarte oprogramowanie** | Dostępny jako narzędzie linii poleceń i biblioteka Java |

---

## Porównanie z innymi metodami

P2RANK osiąga wyższą dokładność niż popularne narzędzia takie jak **Fpocket** czy **SiteHound**, a także dorównuje lub przewyższa metody oparte na głębokim uczeniu, takie jak **DeepSite**. Jego przewaga wynika z połączenia wydajności obliczeniowej z bogatą reprezentacją cech fizykochemicznych.

---

## Zastosowania

- Wirtualny przesiew leków (*structure-based virtual screening*)
- Przewidywanie funkcji białek
- Wyznaczanie miejsc alosterycznych
- Zautomatyzowane potoki bioinformatyczne (*pipelines*)

---

## Podsumowanie

```
Struktura białka (PDB)
        │
        ▼
Generowanie punktów powierzchniowych (SAS)
        │
        ▼
Ekstrakcja cech fizykochemicznych lokalnego sąsiedztwa
        │
        ▼
Klasyfikacja Random Forest → wartość ligandability
        │
        ▼
Klastrowanie punktów o wysokim ligandability
        │
        ▼
Ranking kieszeni wiążących
```

P2RANK jest szczególnie przydatny w zautomatyzowanych badaniach przesiewowych, gdzie kluczowa jest zarówno prędkość działania, jak i brak zależności od zewnętrznych zasobów strukturalnych.

## Format wyjścia

### Prediction output 

   For each structure file `{struct_file}` in the dataset, P2Rank generates several output files:
   * `{struct_file}_predictions.csv`: lists **predicted pockets** in order of score, including each pocket's score, center coordinates, adjacent residues, adjacent protein surface atoms, and a calibrated probability of being a ligand-binding site.
   * `{struct_file}_residues.csv`: lists **all residues** from the input protein along with their scores, mapping to predicted pockets, and a calibrated probability of being a ligand-binding residue.
   * **PyMol and ChimeraX visualizations**: `.pml` and `.cxc` scripts in `visualizations/` directory  with additional files in `data/`.
     * Optional settings:
       * Use `-visualizations 0` to disable visualization generation.
       * Use `-vis_renderers 'pymol,chimerax'` to toggle specific renderers on/off.
       * Use `-vis_copy_proteins 0` to prevent copying protein structures to the visualizations directory (faster, but visualizations won't be portable). 
   * **SAS points data**: coordinates and ligandability scores for solvent-accessible surface (SAS) points are saved in `visualizations/data/{struct_file}_points.pdb.gz`. Here:
     * Residue sequence number (position 23-26) represents the pocket rank (0 indicates no pocket).
     * B-factor column contains predicted ligandability score.