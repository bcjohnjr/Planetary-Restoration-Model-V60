#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np
import pandas as pd

required = [
    'fair_feedback_sensitivity_pair.csv',
    'fair_feedback_sensitivity_selected.csv',
    'fair_feedback_sensitivity_members.npz',
    'fair_feedback_sensitivity_permafrost.csv',
    'fair_feedback_sensitivity_convergence.csv',
    'fair_feedback_sensitivity_forcing_audit.csv',
    'fair_feedback_sensitivity_summary.json',
]
for fn in required:
    assert Path(fn).is_file(), fn

p = pd.read_csv('fair_feedback_sensitivity_pair.csv')
experiments = {
    'legacy_strict_pair',
    'reversal_zero_off_only',
    'dynamic_permafrost_only',
    'combined_programme_effect',
}
assert set(p.experiment) == experiments
for e in experiments:
    z = p[p.experiment == e]
    assert len(z) == 651, (e, len(z))
    assert float(z.timebound_year.min()) == 1750.0
    assert float(z.timebound_year.max()) == 2400.0
    # 2026 and 2027 have zero cumulative programme CDR; response fractions
    # are therefore undefined until timebound 2028. CO2 deltas remain auditable.
    core = z[['delta_co2_p05_ppm','delta_co2_p50_ppm','delta_co2_p95_ppm']]
    assert np.isfinite(core).all().all()
    late = z[z.timebound_year >= 2028]
    assert np.isfinite(late[['fraction_p05','fraction_p50','fraction_p95']]).all().all()

# Exact regression lock to the independently archived Run 8 pair.
def row(exp, year):
    z = p[p.experiment == exp]
    return z.iloc[int(np.argmin(np.abs(z.timebound_year.to_numpy(float)-year)))]

r = row('legacy_strict_pair', 2300)
# Run-8 regression lock. Allow only sub-micro-ppm platform floating-point noise.
expected_delta_2300 = 104.03176804841752
expected_frac_2300 = 0.3083303157940670
expected_frac_2184 = 0.3698648683773295
assert abs(float(r.delta_co2_p50_ppm) - expected_delta_2300) < 5e-7, (float(r.delta_co2_p50_ppm), expected_delta_2300)
assert abs(float(r.fraction_p50) - expected_frac_2300) < 5e-9, (float(r.fraction_p50), expected_frac_2300)
r2184 = row('legacy_strict_pair', 2184)
assert abs(float(r2184.fraction_p50) - expected_frac_2184) < 5e-9, (float(r2184.fraction_p50), expected_frac_2184)

# Removing programme-contingent reversal from the off arm must lower the
# attributed benefit relative to the strict legacy construction.
rz = row('reversal_zero_off_only', 2300)
assert float(rz.delta_co2_p50_ppm) < float(r.delta_co2_p50_ppm)
assert float(rz.fraction_p50) < float(r.fraction_p50)

s = json.load(open('fair_feedback_sensitivity_summary.json'))
assert s['model'] == 'FaIR 2.2.4'
assert s['configs'] == 841
assert s['common_state_co2_maxabs_ppm_at_2026'] <= 1e-12
assert 0.03 < s['canonical_reversal_fraction_of_cdr'] < 0.06
for k, v in s['fixed_point_closure_error_gtco2_yr'].items():
    assert float(v) <= 5e-4, (k, v)

z = np.load('fair_feedback_sensitivity_members.npz', allow_pickle=False)
assert z['selected_years'].shape == (8,)
assert z['config_ids'].shape == (841,)
for k in z.files:
    if k.endswith('_co2') or k.endswith('_tas'):
        assert z[k].shape == (8, 841), (k, z[k].shape)
        assert np.isfinite(z[k]).all(), k

# Explicit forcing audit: programme-contingent reversal is zero only where it
# should be zero, rather than silently altered in the legacy controls.
a = pd.read_csv('fair_feedback_sensitivity_forcing_audit.csv')
for name in ['reversal_zero_off', 'dynamic_off_no_reversal']:
    assert np.max(np.abs(a[a['case']==name].reversal_gtco2_yr.to_numpy(float))) == 0.0
assert np.max(a[a['case']=='legacy_off'].reversal_gtco2_yr.to_numpy(float)) > 0

print('FaIR feedback sensitivity outputs verified.')
