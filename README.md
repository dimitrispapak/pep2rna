# pep2rna

Detect peptides encoded by annotated ORFs in RNA-seq reads.

`pep2rna` scans RNA-seq BAM files for reads whose translated nucleotide
sequence contains any of a user-supplied set of peptides. Candidate
genomic loci for each peptide are restricted to entries in two public
ORF catalogs — **nuORFdb** and **MetamORF** — so the search is targeted
rather than genome-wide. This makes the tool useful for asking
questions like *"is there any RNA-level evidence that peptide X is
expressed?"* when X has been observed at the protein level (e.g. by
mass spectrometry or a peptide array) but its transcript of origin is
unclear.

---

## 1. Goal

Given:

- a list of target **peptides**,
- two **ORF reference catalogs** (nuORFdb, MetamORF) in FASTA format,
- a directory of **RNA-seq BAM files**,

determine, for every (BAM, peptide) pair, how many reads whose
translation contains the peptide map to a locus annotated as encoding
that peptide.

The intended use is exploratory: a hit is evidence, not proof. Absence
of a hit is not evidence of absence (coverage, editing, introns,
sequencing errors, and the incompleteness of the ORF catalogs all
matter).

---

## 2. Pipeline structure

The script runs in five logical stages:

```
  peptides.tsv ──┐
                 │   (1) build pattern → id map
                 ▼
       pattern_to_ids {peptide: {id, ...}}
                 │
                 │   (2) build Aho-Corasick automaton
                 ▼
       Aho-Corasick automaton over peptide patterns
                 │
                 │   (3) scan MetamORF + nuORFdb FASTA headers
                 ▼
       {peptide: {ids, loci[(chrom, start, end, strand)]}}
                 │
                 │   (4) merge sources, collapse overlapping loci
                 ▼
       final peptide → loci table
                 │
                 │   (5) for each BAM, for each peptide, for each locus:
                 │       fetch reads → strand filter → mRNA → translate
                 ▼
       hits.tsv  (bam, peptide, chrom, start, end, strand, reads, hits)
```

### Stage details

**(1) Peptide → identifier map.**
Read `peptides.tsv`. Each row associates a peptide with an identifier.
Only rows whose identifier is in the requested set (`--pattern-ids`,
default: a built-in set) are kept. The output is
`{peptide: {id1, id2, ...}}` — a peptide may map to several ids.

**(2) Aho-Corasick automaton.**
Build a single automaton over all peptide patterns so that all peptides
can be found in a sequence in one pass — O(sequence length) rather than
O(number of peptides × sequence length).

**(3) Scan ORF catalogs.**
For each FASTA in nuORFdb and MetamORF, stream protein sequences. For
each sequence, run the automaton. Every peptide occurrence yields one
locus, derived from the FASTA header. Loci are deduplicated per peptide
so repeated headers don't inflate counts.

The two catalogs use different header conventions, so they are parsed
separately:

- **MetamORF** header example (synthetic):
  `>... chr=chr<C> strand=<+|- > start_pos=<int> stop_pos=<int> ...`
- **nuORFdb** header example (synthetic):
  `>chr<C>_<start>-<end>:<strand>|<...>`

**(4) Merge and collapse.**
Loci from both catalogs are merged per peptide. Overlapping intervals
within the same `(chrom, strand)` are collapsed into a single locus
spanning the union. Strand is part of the grouping key: a `+` and a `-`
locus that overlap are **not** merged, since they encode different
strands of the same region.

**(5) BAM scan.**
For each BAM (in parallel across `--workers` processes), for each
peptide, for each locus:

1. `bam.fetch(chrom, start-1, end)` — 1-based inclusive → 0-based half-open.
2. Drop reads below `--min-mapq` (default `0`).
3. Keep only reads on the locus strand.
4. Reconstruct the read's mRNA-oriented nucleotide sequence using
   `get_aligned_pairs()`, complementing and reversing as needed.
5. Translate the sequence in all three reading frames using the
   standard genetic code (table 1).
6. If the peptide appears in any frame, record a hit.

The result is one row per (BAM, peptide, locus) with ≥1 hit.

---

## 3. Inputs

### 3.1 `peptides.tsv`

A tab-separated file with a header row. Each subsequent row pairs a
peptide sequence with an identifier that groups it with other peptides
from the same source.

| Column | Name (suggested) | Type   | Required | Description |
|--------|------------------|--------|----------|-------------|
| 1      | `peptide`        | string | yes      | Amino-acid sequence to search for. Uppercase, one-letter codes, no whitespace. Must be a substring of some translated ORF to yield a hit. |
| 2      | `identifier`     | string | yes      | Group label. Only identifiers in the requested set (`--pattern-ids`) are kept. |

Additional columns are ignored. Blank lines and rows with fewer than
two fields are skipped.

**Example (synthetic):**

```
peptide	identifier
PEPTIDE1	GROUP_A
PEPTIDE2	GROUP_B
PEPTIDE3	GROUP_C
```

**Notes:**

- The header line is required — the script skips the first line
  unconditionally.
- Peptide sequences are matched **exactly** (case-sensitive). Stick to
  uppercase.
- A peptide appearing in multiple rows is fine; ids are unioned.
- A single `identifier` may be associated with multiple peptides; they
  are treated as separate patterns.
- Any peptide whose length is `L` requires reads of at least `3L`
  aligned nucleotides to be considered — the script silently drops
  shorter reads before translation.
- By default, only a preset list of identifiers is kept. To use a
  different set, pass `--pattern-ids id1,id2,...`.

### 3.2 `--nuorf`: nuORFdb protein FASTA

Protein sequences from nuORFdb. Headers must be of the form:

```
>chr<chrom>_<start>-<end>:<strand>|<...>
```

Synthetic example:

```
>chr<C>_<start>-<end>:<strand>|<transcript_id>|...
```

The script parses `chrom`, `start`, `end`, and `strand` from the part
before the first `|`. Coordinates are treated as 1-based inclusive.

### 3.3 `--metamorf`: MetamORF protein FASTA

Protein sequences from MetamORF. Headers are expected to contain the
following whitespace-separated `key=value` pairs:

```
>... chr=chr<C> strand=<+|- > start_pos=<int> stop_pos=<int> ...
```

Only these four keys are read; all other content is ignored.

### 3.4 `--bam-dir`: BAM directory

A directory containing one or more `.bam` files. Files must be:

- **coordinate-sorted** (required by `pysam.AlignmentFile.fetch`),
- **indexed** — a `.bai` file alongside each `.bam` (or `.csi` for
  large references).

The script does **not** sort or index for you. If your BAMs come
straight out of an aligner, run:

```bash
for b in *.bam; do
    samtools sort -o "${b%.bam}.sorted.bam" "$b"
    samtools index "${b%.bam}.sorted.bam"
done
```

Chromosome naming must match the catalogs. Some catalogs use `chr6`,
others use bare `6`; the script prepends `chr` to whatever the locus
parser yields, so both catalogs end up producing `chr6` for fetching.

---

## 4. Outputs

### 4.1 `hits.tsv` (or `--output`)

A TSV with a header row and one row per (BAM, peptide, locus) with at
least one hit:

| Column     | Type    | Description |
|------------|---------|-------------|
| `bam_file` | string  | BAM filename (not full path). |
| `peptide`  | string  | Matched peptide sequence. |
| `chrom`    | string  | Chromosome with `chr` prefix. |
| `start`    | integer | Locus start, 1-based inclusive. |
| `end`      | integer | Locus end, 1-based inclusive. |
| `strand`   | `+`/`-` | Strand of the annotated ORF. |
| `reads`    | integer | Number of reads examined at the locus after filtering. |
| `hits`     | integer | Of those reads, how many contained the peptide in some frame. |

**Interpretation.**

- `reads` is the denominator: how many reads actually made it through
  the strand and length filters at that locus.
- `hits` is the numerator: how many of those reads, once translated,
  contained the peptide.
- A row with `reads=3, hits=1` and a row with `reads=300, hits=1` mean
  very different things. Always report both.
- A peptide may appear in multiple rows (multiple loci, multiple BAMs).

### 4.2 stderr summary

Information-level messages go to stderr so stdout/TSV stays clean:

```
[info] N distinct peptide patterns for M ids
[info] K / N peptides have loci
[info] B BAM files to scan
Total Hits: H
[info] wrote R rows to hits.tsv
```

`Total Hits` is the sum of the `hits` column — **not** the number of
rows.

---

## 5. Usage

```bash
python pep2rna.py \
    --peptides  peptides.tsv \
    --nuorf     /path/to/nuORFdb_protein.fasta \
    --metamorf  /path/to/MetamORF_Hsapiens_aa.fasta \
    --bam-dir   /path/to/BAM/ \
    --workers   28 \
    --output    hits.tsv
```

Run `python pep2rna.py --help` for the full option list.

### Key options

| Option | Default | Effect |
|---|---|---|
| `--peptides` | required | Path to `peptides.tsv`. |
| `--nuorf` | required | nuORFdb FASTA. |
| `--metamorf` | required | MetamORF FASTA. |
| `--bam-dir` | required | Directory of indexed BAMs. |
| `--pattern-ids` | built-in set | Comma-separated identifiers to keep. |
| `--output` | `hits.tsv` | Output TSV path. |
| `--workers` | `28` | Worker processes for the BAM scan. |
| `--min-mapq` | `0` | Minimum MAPQ for a read to count. |
| `--quiet` | off | Suppress report on stdout (still writes TSV). |

---

## 6. Install

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

`requirements.txt`:

```
pysam>=0.22
pyahocorasick>=2.0
biopython>=1.81
```

---

## 7. Caveats

- **Matching is peptide-in-frame, not exact ORF reconstruction.** A
  read counts as a hit if the peptide appears anywhere in a translated
  frame of the mapped region — including in a frame not annotated as
  the ORF's frame. This is deliberate (it tolerates
  annotation/coordinate disagreement), but it means hits are
  permissive evidence.
- **No intron handling.** The mRNA is reconstructed from aligned
  pairs, so spliced reads are effectively concatenated in reference
  order. Peptides spanning a splice junction may or may not be
  detected depending on the aligner's representation.
- **Deletions are dropped, insertions are retained** in the
  reconstructed mRNA. Frame shifts from indels are not modeled.
- **No error correction.** Sequencing errors that introduce a stop
  codon or a frameshift will suppress a hit.
- **Strand is enforced by the annotated ORF, not inferred** from the
  read. Overlapping opposite-strand genes are treated separately.
- **Catalogs are incomplete.** Absence of a hit does not mean the
  peptide is not expressed; it may simply not be in nuORFdb or
  MetamORF.

---

## 8. License

MIT — see `LICENSE`.
