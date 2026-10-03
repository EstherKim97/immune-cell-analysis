VENV ?= .venv
SYSTEM_PYTHON ?= python3
PYTHON ?= $(VENV)/bin/python

.PHONY: setup pipeline dashboard

setup:
	$(SYSTEM_PYTHON) -m venv $(VENV)
	$(PYTHON) -m pip install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) frequency_table.py
	$(PYTHON) response_analysis.py
	$(PYTHON) subset_analysis.py

dashboard:
	$(PYTHON) -m streamlit run app.py --server.headless true --server.port 8501
