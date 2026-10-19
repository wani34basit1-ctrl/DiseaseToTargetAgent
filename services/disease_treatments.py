"""Curated disease-level treatment information.

These records are separate from target-associated clinical candidates.
They are informational and are not personalized treatment recommendations.
"""

DISEASE_TREATMENTS = {
    "familial hypercholesterolemia": [
        {
            "drug": "Atorvastatin",
            "class": "Statin",
            "mechanism": (
                "Inhibits HMG-CoA reductase, reducing hepatic cholesterol "
                "synthesis and increasing LDL receptor expression."
            ),
            "ldlr_dependence": "Primarily LDLR-dependent",
            "scope": "Commonly used in FH management",
        },
        {
            "drug": "Ezetimibe",
            "class": "Cholesterol absorption inhibitor",
            "mechanism": (
                "Inhibits NPC1L1-mediated intestinal cholesterol absorption, "
                "helping lower circulating LDL cholesterol."
            ),
            "ldlr_dependence": "LDL lowering partly depends on LDLR activity",
            "scope": "Often combined with a statin",
        },
        {
            "drug": "Evolocumab",
            "class": "PCSK9 monoclonal antibody",
            "mechanism": (
                "Inhibits PCSK9, helping preserve LDL receptors and increase "
                "LDL clearance."
            ),
            "ldlr_dependence": "LDLR-dependent",
            "scope": "Used in eligible patients with FH",
        },
        {
            "drug": "Inclisiran",
            "class": "PCSK9-targeting siRNA",
            "mechanism": (
                "Reduces hepatic PCSK9 production, supporting LDL receptor "
                "recycling and LDL cholesterol clearance."
            ),
            "ldlr_dependence": "LDLR-dependent",
            "scope": "Eligibility depends on indication and local approval",
        },
        {
            "drug": "Evinacumab",
            "class": "ANGPTL3 monoclonal antibody",
            "mechanism": (
                "Inhibits ANGPTL3 and lowers LDL cholesterol through a "
                "mechanism that does not require functional LDL receptors."
            ),
            "ldlr_dependence": "LDLR-independent",
            "scope": "Specifically relevant to eligible patients with HoFH",
        },
    ]
}


def get_disease_treatments(disease_name):
    """Return curated treatments for a supported disease."""
    if not disease_name:
        return []

    normalized = " ".join(disease_name.lower().replace("_", " ").split())

    aliases = {
        "fh": "familial hypercholesterolemia",
        "familial hypercholesterolaemia": "familial hypercholesterolemia",
        "homozygous familial hypercholesterolemia": "familial hypercholesterolemia",
        "heterozygous familial hypercholesterolemia": "familial hypercholesterolemia",
    }

    normalized = aliases.get(normalized, normalized)
    return [
        {
            **item,
            "evidence_level": "disease-level",
            "evidence_source": "curated_disease_treatment_catalogue",
        }
        for item in DISEASE_TREATMENTS.get(normalized, [])
    ]
