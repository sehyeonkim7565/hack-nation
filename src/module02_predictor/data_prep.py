"""
Module 02 data prep: build the real, labeled n=25 genome x antibiotic
phenotype table from kpneu_amr_phenotypes.tsv.

IMPORTANT SCOPE NOTE (read before using train_loocv.py):
Module 01's audit found that only 1 of these 25 labeled genomes
(1226680.5) has real, species-matched nucleotide sequence available in
kpneu_genomes.fna. A supervised model trained on AMRFinderPlus/BLAST gene
presence-absence features, as originally scoped, therefore cannot be
trained or honestly evaluated on this cohort -- there is only one
genomic-feature row.

To still deliver a real, non-fabricated statistical baseline that can be
evaluated with LOOCV across all n=25 labeled genomes, train_loocv.py uses
each genome's OTHER real, lab-observed antibiotic phenotypes as features
to predict a held-out target antibiotic ("cross-antibiotic co-resistance"
baseline). This is legitimate, unmodified real data -- but it is a
retrospective-cohort statistical exercise, NOT a genome-derived model. It
cannot be applied to a freshly uploaded genome with no prior susceptibility
testing, which is exactly what Module 03 needs to acknowledge and does
(see predict.py / app.py docstrings).
"""
import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PHENOTYPES_TSV = REPO_ROOT / "kpneu_amr_phenotypes.tsv"

TARGET_ANTIBIOTICS = [
    "trimethoprim/sulfamethoxazole",
    "gentamicin",
    "cefepime",
]

PHENOTYPE_TO_LABEL = {"Resistant": 1, "Susceptible": 0, "Intermediate": None}


def load_labeled_cohort() -> dict:
    """Returns {genome_id: {antibiotic: 'Resistant'|'Susceptible'|'Intermediate'}}
    restricted to genomes with at least one non-blank phenotype call.
    The known-corrupt trailing line in the source TSV (a non-UTF8,
    single-field garbage row -- see project data-quality notes) is
    skipped automatically because it doesn't parse as 4 tab-separated
    fields."""
    table: dict[str, dict[str, str]] = {}
    with open(PHENOTYPES_TSV, encoding="utf-8", errors="ignore") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        for row in reader:
            if len(row) != 4:
                continue
            genome_id, _, antibiotic, phenotype = row
            if not phenotype:
                continue
            table.setdefault(genome_id, {})[antibiotic] = phenotype
    return table


def get_all_antibiotics(table: dict) -> list:
    seen = set()
    for genome_labels in table.values():
        seen.update(genome_labels.keys())
    return sorted(seen)


if __name__ == "__main__":
    table = load_labeled_cohort()
    print(f"Labeled genomes: {len(table)}")
    print(f"Distinct antibiotics with >=1 real label: {len(get_all_antibiotics(table))}")
    for ab in TARGET_ANTIBIOTICS:
        counts = {"Resistant": 0, "Susceptible": 0, "Intermediate": 0}
        for genome_labels in table.values():
            p = genome_labels.get(ab)
            if p in counts:
                counts[p] += 1
        print(f"  {ab}: R={counts['Resistant']} S={counts['Susceptible']} I={counts['Intermediate']}")
