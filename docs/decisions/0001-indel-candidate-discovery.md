# ADR-0001: Reference-normalized indel candidate discovery

Date: 2026-10-21
Status: Accepted for candidate discovery; reconciliation integration pending

## Context

Illumina GSA manifest indel coordinates and VCF indel coordinates
do not necessarily identify equivalent alleles using coordinate
equality alone.

The ArrayEvidence publication-track implementation must establish
allele equivalence independently of genotype interpretation.

Development data: authorized AIDA GSA dataset.
Reference assembly: GRCh38.

No participant-derived records are included in this document.

## Decision

Use reference-normalized sequence events to identify equivalent
simple insertion and deletion alleles.

Candidate discovery must:

1. Validate manifest-derived events against the reference.
2. Validate VCF REF alleles against the reference.
3. Normalize simple indels using reference-aware left normalization.
4. Match by normalized sequence-event identity.
5. Preserve all equivalent VCF records.
6. Distinguish unique, multiple, absent, and unresolved evidence.
7. Never infer homozygous-reference genotypes from absent records.
8. Keep allele equivalence separate from genotype interpretation.

Do not apply a fixed coordinate offset.

## Implementation checkpoint

Implemented components:

- indel_events.py
- indel_normalize.py
- manifest_indel.py
- indel_match.py
- indel_index.py

The normalized candidate index is not yet integrated into
the main reconcile() function.

Full test suite at checkpoint: 132 passing tests.

## AIDA development audit

Manifest indel assays: 10,118
VCF records examined: 646,941

Normalized candidate discovery:

- Unique equivalent VCF record: 9,433
- Multiple equivalent VCF records: 102
- Unresolved manifest reconstruction: 577
- Unmapped assay: 6

Unique-match coordinate relationships:

- Same reported coordinate: 2,375
- VCF coordinate one base upstream: 7,058

Unique-match genotype statuses:

- HOM_REF: 9,370
- CALLED_NONREF: 49
- NO_CALL: 14

All 102 multiple-match cases contained two equivalent VCF
records, both classified as HOM_REF.

Manifest reconstruction failures:

- UNSUPPORTED_FLANKS: 449
- NO_FLANK_ALIGNMENT: 99
- AMBIGUOUS_FLANK_ALIGNMENT: 26
- REFERENCE_UNAVAILABLE: 2
- REFERENCE_ALLELE_NOT_IN_MANIFEST: 1

These are development-cohort observations, not independent
accuracy estimates or fixed performance acceptance targets.

## Interpretation

Coordinate-only candidate selection is insufficient for indels.

Reference-normalized sequence-event matching resolves the
previously observed 387 apparent allele mismatches among
coordinate-selected indel comparisons by identifying equivalent
VCF representations at nearby positions.

Multiple equivalent records remain distinct provenance sources.

A unique equivalent record does not, by itself, establish
genotype accuracy.

## Known limitations

1. Candidate discovery currently indexes simple insertions and
   deletions, not all complex sequence variants.
2. The current index returns VCF records rather than explicit
   (record, ALT index) allele-level candidates.
3. Multiple equivalent VCF records are not yet reconciled at
   genotype level.
4. Manifest reconstruction remains unresolved for some assays.
5. VCF absence must not be interpreted as HOM_REF.
6. AIDA is development data and cannot provide independent
   held-out validation.
7. The main reconciliation pipeline has not yet incorporated
   normalized indel candidate discovery.

## Next milestone: Allele-aware reconciliation

Follow specification-driven red-green-refactor TDD.

Planned sequence:

1. Define allele-level candidate identity and invariants.
2. Write failing synthetic tests.
3. Implement allele-aware candidate indexing.
4. Implement allele-specific genotype projection.
5. Integrate with reconcile() without changing established
   SNV behavior.
6. Add end-to-end synthetic integration tests.
7. Run a predefined AIDA regression audit.
8. Freeze the method before independent held-out validation.

Do not tune matching rules solely to improve AIDA concordance.

## Data governance

The publication repository must not contain participant-derived
genotypes, sequencing data, raw array data, or participant-derived
test fixtures.

Synthetic fixtures must be generated independently of real
participant records.

Authorized datasets remain outside the source repository.
