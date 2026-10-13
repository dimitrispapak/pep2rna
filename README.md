cat > README.md <<'EOF'
# pep2rna

Detect peptides encoded by annotated ORFs in RNA-seq reads.

Given a list of peptides (e.g. from a peptide-array experiment), this
script scans RNA-seq BAM files for reads whose translated sequence
contains one of the peptides, using annotated nuORFdb and MetamORF loci
to restrict the search.

## Install (linux)

```bash
uv venv
source .venv/bin/activate
uv pip install -r requirements.txt
```

## Usage
```bash
python pep2rna.py \
    --peptides  peptides.tsv \
    --nuorf     /path/PA_nuORFdb_v1.2_protein.fasta \
    --metamorf  /path/MetamORF_Hsapiens_aa.fasta \
    --bam-dir   /path/to/BAM/ \
    --workers   28 \
    --output    hits.tsv
```
