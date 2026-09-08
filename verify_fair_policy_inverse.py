#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--target', type=int, choices=[2200,2300,2400], required=True)
args = ap.parse_args()
y = args.target

summary_csv = Path(f'fair_policy_inverse_{y}_summary.csv')
members_csv = Path(f'fair_policy_inverse_{y}_members.csv')
summary_json = Path(f'fair_policy_inverse_{y}_summary.json')
for p in [summary_csv, members_csv, summary_json]:
    assert p.is_file(), p

s = pd.read_csv(summary_csv)
expected = {
    'absolute_280_median',
    'member_relative_baseline_median',
    'absolute_280_p90',
    'member_relative_baseline_p90',
}
assert set(s.criterion) == expected
assert set(s.target_year.astype(int)) == {y}
assert set(s.status).issubset({'solved','already_met','not_bracketed'})
assert np.isfinite(s.canonical_extra_reversal_fraction).all()
assert np.all((s.canonical_extra_reversal_fraction > 0.03) &
              (s.canonical_extra_reversal_fraction < 0.06))

m = pd.read_csv(members_csv)
solved = s[s.status.isin(['solved','already_met'])]
assert len(m) == 841 * len(solved), (len(m), len(solved))

for _, r in solved.iterrows():
    z = m[m.criterion == r.criterion]
    assert len(z) == 841
    rate = float(r.extra_cdr_rate_gtco2_yr)
    assert rate >= 0
    # Leak correction must actually be carried whenever extra CDR is positive.
    if rate > 0:
        assert float(r.extra_reversal_rate_gtco2_yr) > 0
        assert float(r.cumulative_extra_reversal_gtco2) > 0
        assert float(r.net_extra_removal_after_extra_reversal_gtco2) < float(r.cumulative_extra_cdr_gtco2)

    if r.criterion == 'absolute_280_median':
        assert abs(np.median(z.co2_ppm.to_numpy(float)) - float(r.co2_p50)) < 1e-8
        assert float(r.co2_p50) <= 280.05
    elif r.criterion == 'member_relative_baseline_median':
        d = z.co2_minus_baseline_ppm.to_numpy(float)
        assert np.median(d) <= 0.05
    elif r.criterion == 'absolute_280_p90':
        q90 = float(np.quantile(z.co2_ppm.to_numpy(float), .90))
        assert abs(q90 - float(r.co2_p90)) < 1e-8
        assert q90 <= 280.05
        assert float(r.fraction_le_280) >= 0.895
    elif r.criterion == 'member_relative_baseline_p90':
        d = z.co2_minus_baseline_ppm.to_numpy(float)
        q90 = float(np.quantile(d, .90))
        assert abs(q90 - float(r.co2_minus_baseline_p90)) < 1e-8
        assert q90 <= 0.05
        assert float(r.fraction_le_own_baseline) >= 0.895

# p90 policy criterion cannot require less removal than the median criterion.
by = s.set_index('criterion')
if by.loc['absolute_280_p90','status'] == 'solved' and by.loc['absolute_280_median','status'] == 'solved':
    assert float(by.loc['absolute_280_p90','extra_cdr_rate_gtco2_yr']) >= float(by.loc['absolute_280_median','extra_cdr_rate_gtco2_yr'])
if by.loc['member_relative_baseline_p90','status'] == 'solved' and by.loc['member_relative_baseline_median','status'] == 'solved':
    assert float(by.loc['member_relative_baseline_p90','extra_cdr_rate_gtco2_yr']) >= float(by.loc['member_relative_baseline_median','extra_cdr_rate_gtco2_yr'])

j = json.load(open(summary_json))
assert j['model'] == 'FaIR 2.2.4' and j['configs'] == 841 and j['target_year'] == y
assert 'Additional post-2183 CDR carries' in j['extra_leakage_treatment']

print(f'FaIR policy inverse {y} outputs verified.')
