# Disease2Target Agent

An object-oriented Python mini project for retrieving and ranking
disease-associated therapeutic target candidates.

## Features
- Disease-name search through Open Targets
- Disease and target classes using Python OOP
- Disease-target association retrieval
- Ranking by database association score
- CSV report generation
- Automated tests using pytest

## Installation
    python3 -m venv .venv
    source .venv/bin/activate
    python -m pip install -r requirements.txt

## Run
    python main.py --disease "Alzheimer disease" --top 10

## Tests
    pytest -v

## Data source
https://platform.opentargets.org/

## Limitation
Association scores indicate database evidence, not clinical efficacy.
This project is intended for educational and research purposes.
