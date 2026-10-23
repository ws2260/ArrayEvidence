# ArrayEvidence — AI Development Instructions

## Mandatory session recovery

Before proposing implementation, analysis, or a numbered step:

1. Read `CURRENT_STATE.md`.
2. Read `PROJECT_CHARTER.md` if present.
3. Read `docs/ARCHITECTURE.md` and `docs/ROADMAP.md`.
4. Read relevant architecture decisions and validation documents.
5. Inspect Git status, recent commits, and relevant source/tests.
6. Verify the active environment before claiming test results.

If the current state is unavailable, outdated, or contradictory,
ask the user for the latest handoff and execution logs.

Never reconstruct missing project facts from conversational memory.

## Architecture and release discipline

- Preserve the agreed Module 1, 2A, 2B, and 3 boundaries.
- Treat `docs/ARCHITECTURE.md` as authoritative.
- Distinguish released features from proposed functionality.
- Keep v0.1.0 validation separate from subsequent development.
- Do not silently change scientific semantics or output contracts.
- Require evidence, regression tests, and documented decisions
  for scientific-method changes.

## Scientific integrity and data governance

- Never reconstruct participant-level results or genotype calls
  from conversational memory.
- Verify scientific claims against source artifacts and logs.
- Preserve technical ambiguity, provenance, and failure states.
- Do not equate allele reconciliation with genotype accuracy.
- Keep participant-derived data and restricted resources
  outside public Git.
- Use synthetic fixtures for repository tests.
- Do not use generative AI for production clinical adjudication.

## Working protocol

- Work one command or reviewable step at a time.
- Wait for the user's output before proceeding.
- Never repeat completed work without a verified reason.
- Do not discard, overwrite, stage, or commit existing changes
  without explicit user approval.
- Do not claim unexecuted tests or workflows have passed.
- Prioritize the current grant progress report when applicable.

## Handoff maintenance

After significant milestones, update `CURRENT_STATE.md` with:

- Date and verified Git state.
- Completed operations and actual test results.
- Files changed.
- Release and module status.
- Unresolved questions and pending decisions.
- Exact next action.

If the session loses context, explicitly say:

"I cannot reliably reconstruct the current ArrayEvidence state.
Please provide the latest CURRENT_STATE.md and relevant logs
before we continue."

Do not guess.
