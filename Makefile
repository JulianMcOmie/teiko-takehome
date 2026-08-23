# Entry points for the grader. `make setup` once, then `make pipeline`, then
# `make dashboard`.

PYTHON ?= python3

.PHONY: setup pipeline dashboard clean

setup:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -r requirements.txt

# Runs end to end with no manual intervention: builds the database from the CSV,
# then writes every table and figure into outputs/.
pipeline:
	@echo "=== Part 1: building the database ==="
	$(PYTHON) load_data.py
	@echo ""
	@echo "=== Part 2: cell population frequencies ==="
	$(PYTHON) -m src.summary
	@echo ""
	@echo "=== Part 3: responders vs non-responders ==="
	$(PYTHON) -m src.stats
	@echo ""
	@echo "=== Part 4: baseline subset ==="
	$(PYTHON) -m src.subsets
	@echo ""
	@echo "Pipeline complete. Tables and figures are in outputs/."

dashboard:
	$(PYTHON) -m streamlit run dashboard.py

clean:
	rm -f cell_counts.db
	rm -rf outputs __pycache__ src/__pycache__
