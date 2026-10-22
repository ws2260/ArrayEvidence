# Local DRAGEN Array: IDAT to VCF

## Scope

This vignette documents the local DRAGEN Array workflow used
for AIDA GSA processing in October 2026.

Pipeline: IDAT -> GTC -> GRCh38 VCF -> ArrayEvidence.

ArrayEvidence starts downstream of genotype calling.
DRAGEN Array and Illumina assay resources are not distributed
with this repository.

## Tested environment

- Platform: Linux executable under WSL
- DRAGEN Array CLI: 1.4.0
- Runtime: conda environment with compatible ICU libraries

Configure the executable and runtime library paths:

    DRAGENA="/path/to/dragena"
    DRAGENA_LIB="/path/to/conda/envs/dragen-array/lib"
    export LD_LIBRARY_PATH="$DRAGENA_LIB:${LD_LIBRARY_PATH:-}"

The ICU configuration addresses the error:
"Couldn't find a valid ICU package installed".

In the tested CLI, --help may return exit code 1 even when
help text is displayed correctly.

## Required resources

- Paired red/green IDAT files
- Compatible BPM assay manifest
- Compatible EGT cluster model
- CSV assay manifest containing SourceSeq
- GRCh38 reference FASTA and corresponding .fai

The AIDA workflow used GSA-24v3-0_A2 BPM and CSV files,
and GSA-24v3-0_A1 EGT. Verify resource compatibility
before applying these resources to another cohort.

## Configure paths

    IDAT="/path/to/idat"
    BPM="/path/to/manifest.bpm"
    EGT="/path/to/clusters.egt"
    CSV="/path/to/manifest.csv"
    REF="/path/to/GRCh38.fa"
    GTC="/path/to/output/gtc"
    VCF="/path/to/output/vcf"

## Stage 1: IDAT to GTC

    "$DRAGENA" genotype call --idat-folder "$IDAT" --bpm-manifest "$BPM" --cluster-file "$EGT" --output-folder "$GTC" --num-threads 4 --gencall-cutoff 0.15

The thread count and cutoff reproduce the recorded
AIDA command; they are not universal recommendations.

Review calling logs, sample completeness, and call-rate
QC before conversion.

## Stage 2: GTC to VCF

    "$DRAGENA" genotype gtc-to-vcf --gtc-folder "$GTC" --bpm-manifest "$BPM" --csv-manifest "$CSV" --genome-fasta-file "$REF" --output-folder "$VCF"

This command interface was verified against the local
DRAGEN Array CLI 1.4.0 help output.

The AIDA workflow produced 48 readable compressed VCFs.
This does not establish genotype accuracy.

## Stage 3: ArrayEvidence

For each single-sample VCF:

    arrayevidence build --manifest /path/to/assay_manifest.csv --vcf /path/to/sample.vcf.gz --reference "$REF" --sample-id PSEUDONYM --output /path/to/new_bundle

Verify the bundle:

    arrayevidence verify /path/to/new_bundle

Generate aggregate QC:

    arrayevidence qc /path/to/new_bundle

## Reproducibility

Record software versions, input checksums, assay resource
versions, genome reference, command lines, exit codes,
calling QC, output checksums, and ArrayEvidence commit.

Keep participant-derived data, proprietary resources,
and private filesystem paths outside the public repository.

## Interpretation

VCF absence is not evidence of homozygous reference.
ArrayEvidence allele MATCH does not establish genotype
concordance or clinical validity.

For independent cohorts, freeze the upstream calling
protocol and ArrayEvidence version before processing data.
