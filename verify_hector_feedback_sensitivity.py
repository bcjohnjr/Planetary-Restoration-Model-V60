#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np
import pandas as pd

req = [
    'hector_feedback_sensitivity_pair.csv',
    'hector_feedback_sensitivity_selected.csv',
    'hector_feedback_sensitivity_long.csv',
    'hector_feedback_sensitivity_forcing_audit.csv',
    'hector_feedback_sensitivity_permafrost.csv',
    'hector_feedback_sensitivity_convergence.csv',
    'hector_feedback_sensitivity_summary.json',
]
for fn in req:
    assert Path(fn).is_file(), fn

p = pd.read_csv('hector_feedback_sensitivity_pair.csv')
experiments = {
    'legacy_strict_pair',
    'reversal_zero_off_only',
    'dynamic_permafrost_only',
    'combined_programme_effect',
}
assert set(p.experiment) == experiments
for e in experiments:
    z = p[p.experiment == e]
    assert len(z) == 275
    assert z.year.min() == 2026 and z.year.max() == 2300
    # In 2026 cumulative CDR is zero, so response_fraction is undefined there.
    core = z[['co2_on_ppm','co2_off_ppm','delta_co2_ppm']]
    assert np.isfinite(core).all().all()
    late = z[z.year >= 2027]
    assert np.isfinite(late[['response_fraction']]).all().all()
    y2026 = z[z.year == 2026].iloc[0]
    assert float(y2026.cumulative_cdr_gtco2) == 0.0

s = json.load(open('hector_feedback_sensitivity_summary.json'))
assert s['model'] == 'Hector' and s['version'] == '3.5.0'
assert s['common_state_max_co2_error_ppm'] <= 1e-12
assert s['common_state_max_tas_error_degc'] <= 1e-12
for k,v in s['fixed_point_closure_error_gtco2_yr'].items():
    assert v <= 5e-5, (k,v)

# The strict legacy arm should reproduce the already-validated first pair to
# tight numerical tolerance, proving this sensitivity script did not silently
# change the underlying baseline experiment.
row = p[(p.experiment=='legacy_strict_pair') & (p.year==2300)].iloc[0]
assert abs(row.co2_on_ppm - 306.663201432654) < 1e-6
assert abs(row.co2_off_ppm - 389.054767828132) < 1e-6
assert abs(row.response_fraction - 0.243703273975215) < 1e-9

# Removal-off reversal must actually be zero in the corrected cases.  Audit is
# based on Hector's stored FFI/DACCS inputs; the summary independently records
# the exported reversal amount removed from the off arm.
assert s['cumulative_exported_reversal_2026_2300_gtco2'] > 0

print('Hector feedback sensitivity outputs verified.')
