# rtk2 mosaic mode

`mosaic` reports pairs of sequence-related gene representatives from the same MGS
whose abundance profiles support the interpretation that they are variants of one
locus. It does not merge or rewrite the input matrix.

## Inputs

```text
rtk2 mosaic \
  -i gene_abundance.tsv.gz \
  -reference SB.clusters.core \
  -paf catalog.self.paf.gz \
  [-groups gene_to_nog.tsv] \
  -o output_prefix
```

- The abundance matrix has gene IDs in column 1, a header, and one non-negative
  abundance per sample. Gzip input is supported.
- The MGS/core file uses MGS in column 1 and either one gene or a comma-separated
  gene list in column 2. In six-column `SB.clusters.core`, columns 3–6 are the
  observation count, multicopy count, weighted MGS-assignment score, and marker
  designation. Baseline genes reproduce the MATAFILER core criteria: assignment
  score between 0.9 and 3 and multicopy fraction below 0.1, relaxed to 0.2 for
  marker genes.
- The optional group file restricts candidates to the same NOG or other functional
  group. Supported layouts include `gene<TAB>group`, `group<TAB>gene`,
  `MGS<TAB>gene<TAB>group`, and `group<TAB>count<TAB>gene1,gene2,...`.
- The PAF must describe catalog self-alignments and include both query directions
  for a pair. Generate base-level alignments with minimap2 `-c` if
  `-mosaicRequireCigar` will be used. Retain secondary alignments because the mode
  uses the runner-up hit to reject ambiguous partners.

## Evidence model

1. Keep whole-gene PAF relationships passing identity, query coverage, target
   coverage, and length-ratio thresholds.
2. Require a reciprocal best partner. Hits are ranked by matching bases divided by
   the longer gene length; the best hit must beat the runner-up by the configured
   relative score margin.
3. Build a per-sample MGS baseline as the median abundance of other eligible core
   genes. Markers are used when at least `-mosaicMinBaselineGenes` remain;
   otherwise all eligible core genes are used. Candidate genes are excluded from
   their own baseline. `-mosaicAllCoreBaseline` disables marker preference.
4. Divide both candidate abundances by the MGS baseline. The combined normalized
   copy number should be close to one and fit one copy materially better than
   either candidate alone.
5. Reject pairs whose detected variants overlap too much. Sample overlap is joint
   support divided by the smaller variant support. Abundance overlap is the sum of
   per-sample minimum copy numbers divided by the smaller total variant abundance.

The default overlap limits are both 0.15, preserving the intent of the Perl
presence/absence safeguard while using abundance as additional evidence. Every
sequence relationship receives an explicit confirmation or rejection reason.

## Mosaic versus decluter

The modes answer related but different questions:

- `decluter` starts from CD-HIT clusters, tests abundance co-exclusion within
  those supplied clusters, merges accepted rows, and writes a replacement
  matrix plus `concat.list`.
- `decluter2` performs the same kind of matrix consolidation with MMseqs2-style
  cluster input.
- `mosaic` starts from genes assigned to the same MGS, requires reciprocal
  whole-gene PAF evidence (and optionally a shared functional group), and tests
  whether their MGS-normalized abundances complement one another. It writes a
  decision report and confirmed-pair list but does not alter the matrix.

Use `decluter` or `decluter2` when the desired result is a broadly consolidated
matrix driven by an existing sequence-cluster definition. Use `mosaic` when the
question is specifically whether representatives split by catalog clustering are
alternative forms of one MGS locus and the evidence must remain auditable.

## Outputs

- `<prefix>.mosaic.tsv`: all candidate relationships, sequence metrics, normalized
  abundance metrics, status, and rejection reason.
- `<prefix>.concat.list`: only confirmed pairs as `left,right`, suitable for the
  existing downstream concatenation workflow.
- `<prefix>.mosaic.summary.tsv`: input and decision counts.

The first version is intentionally pairwise and report-only. Outgroup selection,
alignment construction, and tree building remain in the Perl workflow; this keeps
those phylogenetic steps separate from evidence-based candidate discovery.

## MATAFILER4 integration

[`integration/MATAFILER4_mosaic_rtk.patch`](../integration/MATAFILER4_mosaic_rtk.patch)
updates `prepare_mosaic_loci.pl` and its workflow regression. Apply it from a
MATAFILER4 checkout with:

```text
git apply /path/to/rtk/integration/MATAFILER4_mosaic_rtk.patch
```

The integrated order is:

1. Select genes with shared NOGs and write the minimap2 query FASTA.
2. Write or reuse `<output>.minimap2.paf` using minimap2 `-c`.
3. Run configured program `rare` (rtk2) in `mosaic` mode against
   `$GCd/Matrix.mat.gz`, the MGS core file, PAF, and temporary gene-to-NOG map.
4. Publish rtk-confirmed pairs in the existing MATAFILER catalogue format.
5. Reuse the PAF in Perl for outgroup selection and leave later phylogenetic
   processing unchanged.
