# Disease2Target Agent

A bioinformatics project for identifying and analyzing disease associated therapeutic target candidates using biomedical databases and Python.

## Features

- Disease name search using Open Targets
- Retrieval of disease associated targets
- Ranking of targets using database association scores
- Automatic selection of a candidate target
- Protein information retrieval from UniProt
- Retrieval of recent matching PubMed publications with titles, IDs, publication dates, and links
- Biological pathway retrieval from Reactome
- Retrieval of clinical candidates and available mechanism of action information
- CSV report generation
- Detailed JSON report generation
- Automated testing using pytest

## Installation

Create a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install the required dependencies:

```bash
python -m pip install -r requirements.txt
```

## Usage

Run the project with a disease name:

```bash
python main.py --disease "Alzheimer disease" --top 10
```

Other examples:

```bash
python main.py --disease "Parkinson disease" --top 10
python main.py --disease "Type 2 diabetes mellitus" --top 10
python main.py --disease "Familial hypercholesterolemia" --top 10
```

## Output Reports

The project generates the following reports:

- `results/target_report.csv` contains ranked disease associated targets, gene identifiers, association scores, and candidate selection status.
- `results/target_report_details.json` contains disease information, selected target details, UniProt entries, PubMed publications, Reactome pathways, clinical candidates, and scientific interpretation.

## Testing

Run the automated tests using:

```bash
python -m pytest -v
```

## Data Sources

- Open Targets: https://platform.opentargets.org/
- UniProt: https://www.uniprot.org/
- PubMed: https://pubmed.ncbi.nlm.nih.gov/
- Reactome: https://reactome.org/

## Limitations

Association scores represent database evidence and are not probabilities of clinical efficacy. Biological pathways provide context but do not necessarily represent a drug mechanism of action. Clinical candidates associated with a target are not necessarily approved treatments for the selected disease.

This project is intended for educational and research purposes and does not provide clinical recommendations.
