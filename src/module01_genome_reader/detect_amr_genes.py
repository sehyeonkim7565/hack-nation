"""
Module 01 step B: detect acquired AMR genes in a genome FASTA via BLAST.

AMRFinderPlus was the tool specified in scope, but it could not be
installed in this environment:
  - No conda/mamba is available (the specified install path).
  - Direct network access to github.com and ftp.ncbi.nlm.nih.gov (where
    AMRFinderPlus's binary and its curated reference DB are hosted) is
    blocked by the sandbox's egress policy (confirmed 403 policy denial).

Fallback used instead (per project instructions): BLAST against a curated
NCBI-style acquired-resistance-gene reference database. Specifically, the
`resfinder-db` Debian package (CGE ResFinder database, a real, citable,
widely-used curated database of acquired antimicrobial resistance genes)
installed via `apt-get install resfinder-db`. This is NOT a fabricated or
placeholder database -- it is the actual upstream ResFinder gene set,
shipped with its own gene -> phenotype (antibiotic) mapping table
(phenotypes.txt), which this script uses to translate gene hits into
per-antibiotic resistance-gene evidence.

Limitation vs. AMRFinderPlus: ResFinder's database covers acquired
(horizontally-transferred) resistance genes and point-mutation panels are
not applied here; it does not include AMRFinderPlus's curated
point-mutation and some chromosomal-gene logic. Hits are therefore a
lower bound on true resistance-determinant content, and chromosomal
intrinsic mechanisms (e.g. SHV beta-lactamase) are handled separately by
the deterministic gate in Module 02, not by this BLAST step.

Thresholds: >=90% nucleotide identity and >=60% subject (reference gene)
coverage, matching ResFinder's own default calling thresholds.
"""
import csv
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
BLAST_DB = REPO_ROOT / "data" / "blast_db" / "resfinder"
PHENOTYPES_TSV = REPO_ROOT / "reference_db" / "resfinder_phenotypes.txt"
GENOME_DIR = REPO_ROOT / "data" / "genomes"
RESULTS_DIR = REPO_ROOT / "results"

MIN_IDENTITY = 90.0
MIN_COVERAGE = 60.0

BLAST_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "slen",
]


def load_gene_phenotypes(path: Path) -> dict:
    """gene_accession -> {"class": str, "antibiotics": [str, ...]}"""
    mapping = {}
    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            gene_acc, drug_class, phenotype = row[0], row[1], row[2]
            antibiotics = [a.strip() for a in phenotype.split(",") if a.strip()]
            mapping[gene_acc] = {"class": drug_class, "antibiotics": antibiotics}
    return mapping


def run_blast(genome_fasta: Path) -> list[dict]:
    cmd = [
        "blastn",
        "-query", str(genome_fasta),
        "-db", str(BLAST_DB),
        "-outfmt", "6 " + " ".join(BLAST_FIELDS),
        "-perc_identity", str(MIN_IDENTITY),
        "-evalue", "1e-20",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    hits = []
    for line in proc.stdout.strip().splitlines():
        if not line:
            continue
        values = line.split("\t")
        row = dict(zip(BLAST_FIELDS, values))
        row["pident"] = float(row["pident"])
        row["length"] = int(row["length"])
        row["slen"] = int(row["slen"])
        row["coverage"] = 100.0 * row["length"] / row["slen"]
        row["bitscore"] = float(row["bitscore"])
        row["qstart"] = int(row["qstart"])
        row["qend"] = int(row["qend"])
        hits.append(row)
    return hits


def collapse_to_best_hit_per_locus(hits: list[dict]) -> list[dict]:
    """A single physical resistance-gene locus matches many near-identical
    reference alleles (e.g. dozens of blaSHV-* alleles hitting the same
    chromosomal SHV gene). Keep only the single best (highest bitscore)
    hit per overlapping genomic locus, per query contig, so the presence
    matrix reflects genes, not database redundancy."""
    by_contig: dict[str, list[dict]] = {}
    for h in hits:
        by_contig.setdefault(h["qseqid"], []).append(h)

    kept = []
    for contig, contig_hits in by_contig.items():
        contig_hits.sort(key=lambda h: min(h["qstart"], h["qend"]))
        loci: list[list[dict]] = []
        for h in contig_hits:
            lo, hi = sorted((h["qstart"], h["qend"]))
            placed = False
            for locus in loci:
                locus_lo = min(min(x["qstart"], x["qend"]) for x in locus)
                locus_hi = max(max(x["qstart"], x["qend"]) for x in locus)
                if lo <= locus_hi and hi >= locus_lo:
                    locus.append(h)
                    placed = True
                    break
            if not placed:
                loci.append([h])
        for locus in loci:
            best = max(locus, key=lambda h: h["bitscore"])
            kept.append(best)
    return kept


def detect_genes_for_genome(genome_id: str, gene_phenotypes: dict) -> list[dict]:
    genome_fasta = GENOME_DIR / f"{genome_id}.fna"
    if not genome_fasta.exists():
        raise FileNotFoundError(f"No cleaned FASTA found for {genome_id}: {genome_fasta}")

    raw_hits = [h for h in run_blast(genome_fasta) if h["coverage"] >= MIN_COVERAGE]
    best_per_locus = collapse_to_best_hit_per_locus(raw_hits)
    calls = []
    for hit in best_per_locus:
        gene_acc = hit["sseqid"]
        pheno = gene_phenotypes.get(gene_acc, {"class": "unknown", "antibiotics": []})
        calls.append({
            "genome_id": genome_id,
            "gene_accession": gene_acc,
            "gene_name": gene_acc.split("_")[0],
            "identity_pct": round(hit["pident"], 2),
            "coverage_pct": round(hit["coverage"], 2),
            "drug_class": pheno["class"],
            "antibiotics": ";".join(pheno["antibiotics"]),
        })
    return calls


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    gene_phenotypes = load_gene_phenotypes(PHENOTYPES_TSV)

    genome_ids = sorted(p.stem for p in GENOME_DIR.glob("*.fna"))
    if not genome_ids:
        print("No cleaned genome FASTA files found. Run split_genomes.py first.")
        return

    all_calls = []
    for genome_id in genome_ids:
        calls = detect_genes_for_genome(genome_id, gene_phenotypes)
        all_calls.extend(calls)
        print(f"{genome_id}: {len(calls)} acquired resistance gene(s) detected "
              f"(identity>={MIN_IDENTITY}%, coverage>={MIN_COVERAGE}%)")
        for c in calls:
            print(f"    {c['gene_name']:20s} id={c['identity_pct']:6.2f}% cov={c['coverage_pct']:6.2f}% "
                  f"class={c['drug_class']:15s} antibiotics={c['antibiotics']}")

    out_path = RESULTS_DIR / "amr_gene_hits.tsv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=[
            "genome_id", "gene_accession", "gene_name", "identity_pct",
            "coverage_pct", "drug_class", "antibiotics",
        ])
        writer.writeheader()
        writer.writerows(all_calls)
    print(f"\nWrote {len(all_calls)} total gene calls -> {out_path}")


if __name__ == "__main__":
    main()
