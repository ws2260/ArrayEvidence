# ArrayEvidence Synthetic Performance Benchmark

## Scope

This benchmark characterizes the runtime and peak memory usage of
the ArrayEvidence `build` command using deterministic synthetic SNV
inputs. It does not evaluate genotype accuracy.

## Implementation

- Git commit: e1eed5c370352c9ec8dff92552fc63e15844f1c2
- Development state: unreleased v0.1.1 work
- Installed package metadata: 0.1.0 (editable installation)

## Environment

- OS: WSL2, Linux 6.6.87.2-microsoft-standard-WSL2
- CPU: Intel Core Ultra 9 185H (22 logical CPUs available)
- Memory: 29 GiB available to WSL
- Python: 3.12.14
- Measurement tool: GNU time

## Workload

A fixed-seed generator produced synthetic reference sequences,
Illumina-style assay manifests, and single-sample VCFs.

- Assay counts: 1,000; 10,000; 100,000
- VCF coordinate coverage: 80%
- Variant class: biallelic SNVs
- Genotypes: 0/0, 0/1, 1/1
- Random seed: 20261023 (Python random.Random)
- Generator dependency: pysam (including FASTA indexing)
- Reference-aware allele validation
- Three independent builds per workload
- Separate output directory for each run

## Results

| Assays | VCF records | Runtime runs (s) | Mean runtime (s) | Mean peak RSS (MiB) |
|---:|---:|---|---:|---:|
| 1,000 | 800 | 0.11, 0.11, 0.10 | 0.107 | 33.75 |
| 10,000 | 8,000 | 0.44, 0.43, 0.42 | 0.430 | 58.68 |
| 100,000 | 80,000 | 4.08, 4.07, 3.91 | 4.020 | 297.27 |

All nine builds succeeded. One representative bundle per workload
was independently verified and its QC classifications matched the
expected synthetic assay and candidate counts.

## Interpretation

Runtime increased approximately 9.35-fold between 10,000 and
100,000 assays, consistent with approximately linear scaling
over this workload interval.

## Limitations

These results describe one synthetic SNV workload on one WSL2
environment. They do not establish performance for real GSA
manifests, indels, multiallelic variants, different coordinate
coverage, other hardware, or concurrent workloads.

The benchmark does not establish genotype accuracy or biological
validity.

## Reproducing the benchmark

Run the following commands from the repository root in an
environment with ArrayEvidence installed.

Generate synthetic inputs:

    python benchmarks/generate.py --assays 10000 \
      --output /tmp/arrayevidence-benchmark-reproduction

Measure a build:

    /usr/bin/time -f 'Elapsed: %e s | Peak RSS: %M KB' \
      arrayevidence build \
      --manifest /tmp/arrayevidence-benchmark-reproduction/synthetic_manifest.csv \
      --vcf /tmp/arrayevidence-benchmark-reproduction/synthetic.vcf \
      --reference /tmp/arrayevidence-benchmark-reproduction/synthetic.fa \
      --sample-id SYNTHETIC_BENCHMARK \
      --output /tmp/arrayevidence-benchmark-reproduction/bundle

Verify the result:

    arrayevidence verify /tmp/arrayevidence-benchmark-reproduction/bundle
    arrayevidence qc /tmp/arrayevidence-benchmark-reproduction/bundle

Generate each workload once, then reuse its synthetic inputs
for three independent builds. Each build must use a new bundle
output path because ArrayEvidence refuses to overwrite existing
bundles. Repeat with 1,000, 10,000, and 100,000 assays.

GNU time reports peak RSS in KiB; divide by 1,024 to obtain MiB.
Results may vary with hardware, software versions, and system load.
