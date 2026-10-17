# Disease2Target Agent

A Python bioinformatics project for identifying, prioritizing, and analyzing disease-associated therapeutic target candidates using biomedical databases.

## Features

- Disease search using Open Targets
- Retrieval and ranking of disease-associated targets
- Automatic selection of the highest-ranked target
- Protein information retrieval from UniProt
- Recent matching PubMed publications with titles, PMIDs, dates, and links
- Biological pathway retrieval from Reactome
- Target-associated clinical candidates and available mechanism information
- Curated disease-level treatment options for familial hypercholesterolemia
- LDL receptor dependence information for curated FH treatments
- CSV reports and detailed JSON reports
- Automated testing with pytest

## Workflow

```mermaid
flowchart TD
    A["Disease name input"] --> B["Open Targets disease search"]
    B --> C["Retrieve disease-associated targets"]
    C --> D["Rank targets by association score"]
    D --> E["Select top-ranked target"]
    E --> F["Retrieve UniProt protein information"]
    E --> G["Retrieve matching PubMed publications"]
    E --> H["Retrieve Reactome pathways"]
    E --> I["Retrieve target-associated clinical candidates"]
    B --> J["Retrieve curated disease-level treatments when supported"]
    F --> K["Generate CSV and JSON reports"]
    G --> K
    H --> K
    I --> K
    J --> K
