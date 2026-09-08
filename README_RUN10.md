# Run 10 — V60 referee close-out

This workload addresses the V60 cross-check referee report dated 8 September 2026.

## What it closes

1. **FaIR programme-feedback attribution correction**
   - Reproduces the archived Run 8 strict pair as a regression control.
   - Sets stored-carbon reversal to zero in the removal-off world.
   - Makes permafrost release configuration-specific and responsive to each arm's own FaIR surface-temperature trajectory.
   - Runs the two corrections separately and together.

2. **Hector programme-feedback sensitivity**
   - Runs the previously prepared Run 9 factorial correction in the same workflow.

3. **Leakage-corrected inverse solve**
   - Additional post-2183 CDR carries additional stored-carbon reversal at the canonical cumulative reversal/CDR fraction.
   - Solves 2200, 2300 and 2400 independently.

4. **p90 restoration criterion**
   - Solves both ensemble-median and 90th-percentile criteria for 280 ppm.
   - Also reports the equivalent member-relative-preindustrial-baseline criteria.

## Files to put in repository root

- `external_validation_net_co2_trajectory.csv`
- `requirements.txt`
- `run_fair_referee_closeout.py`
- `verify_fair_feedback_sensitivity.py`
- `verify_fair_policy_inverse.py`
- `run_hector_feedback_sensitivity.R`
- `verify_hector_feedback_sensitivity.py`

Put `science-gates-run10.yml` at:

`.github/workflows/science-gates-run10.yml`

Then run **Actions → V60 referee closeout — feedback + policy inverse → Run workflow**.

## Expected jobs

Five jobs should be visible: one FaIR feedback job, three independent FaIR inverse jobs (2200/2300/2400), and one Hector feedback job. Successful jobs upload separate artifacts so one failure does not destroy results from the others.

## Deliberately not claimed by this workload

A cross-model identical neutral non-CO2 forcing prescription is not attempted here because FaIR and Hector have different species/input architectures. The manuscript should state that the residual FaIR/Hector difference is of the same order as the already measured FaIR background sensitivity and remains physically unresolved unless a separately justified common-forcing protocol is developed.
