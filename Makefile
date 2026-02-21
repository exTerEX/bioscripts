.PHONY: help install lint format clean

help:
	@echo "Available targets:"
	@echo "  install   - Install dependencies with uv"
	@echo "  test      - Run all tests"
	@echo "  coverage  - Run tests with coverage report"
	@echo "  docs      - Run sphinx-build for generating documentation"
	@echo "  lint      - Run linting checks"
	@echo "  format    - Format code with ruff"
	@echo "  clean     - Remove build artifacts"

install:
	uv sync --dev

lint: install
	uv run ruff check src/dupfinder/
	uv run pyright src/dupfinder/

format: install
	uv run ruff format src/dupfinder/ tests/
	uv run ruff check --fix src/dupfinder/ tests/

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf .ruff_cache/
	find . -type d -name __pycache__ -exec rm -rf {} +
