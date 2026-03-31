init:
	uv init && uv sync --all-groups

sync:
	uv sync --all-groups

# can be used instead of uv or for deployment purposes
generate_requirements:
	uv export --all-groups --format requirements-txt > requirements.txt

format:
	uv run ruff format

test:
	uv run pytest