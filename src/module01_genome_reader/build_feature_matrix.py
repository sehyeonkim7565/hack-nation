"""
Module 01 step C: turn results/amr_gene_hits.tsv into a reusable
genome x gene presence/absence binary matrix.

This is genome-count-agnostic: point it at any set of genome FASTA files
run through split_genomes.py + detect_amr_genes.py and it will produce a
matrix with one row per genome and one column per distinct gene family
observed across the cohort. Today that cohort has 2 real genomes; the
script scales unchanged if more sequenced genomes are added later.
"""
import csv
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
HITS_TSV = REPO_ROOT / "results" / "amr_gene_hits.tsv"
OUT_TSV = REPO_ROOT / "results" / "gene_presence_absence_matrix.tsv"


def main():
    genome_genes = defaultdict(set)
    with open(HITS_TSV, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genome_genes[row["genome_id"]].add(row["gene_name"])

    all_genes = sorted({g for genes in genome_genes.values() for g in genes})
    genome_ids = sorted(genome_genes.keys())

    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["genome_id"] + all_genes)
        for gid in genome_ids:
            row = [gid] + [1 if g in genome_genes[gid] else 0 for g in all_genes]
            writer.writerow(row)

    print(f"{len(genome_ids)} genomes x {len(all_genes)} gene features -> {OUT_TSV}")


if __name__ == "__main__":
    main()
