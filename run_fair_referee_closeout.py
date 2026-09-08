#!/usr/bin/env python3
"""FaIR 2.2.4 referee close-out experiments for V60/V61.

This script addresses the remaining scientific objections in the V60 cross-check:

1) Paired attribution feedback sensitivity
   * legacy strict pair (regression control)
   * stored-carbon reversal zero in removal-off only
   * dynamic permafrost only
   * combined programme-effect pair: reversal=0 in removal-off and permafrost
     responds to each arm/configuration's own temperature trajectory.

2) Policy inverse solve
   * carries programme-contingent leakage on the *additional* post-2183 CDR
     at the canonical cumulative reversal / cumulative CDR fraction;
   * solves both median and p90 criteria for absolute 280 ppm and each
     configuration's own calibrated preindustrial baseline.

The dynamic-permafrost formulation is the same one used in the Hector Run 9
sensitivity: 18 PgC per degC above 1.20 C, running-peak commitment, 50-year
release e-folding, no reinjection of historical commitment already implicit in
the common state.  FaIR's 841 calibrated configurations are coupled separately:
permafrost is configuration-specific because each member has its own TAS path.

Important indexing note: programme emissions for year y occupy FaIR timepoint
y+0.5 and affect timebound y+1.  The pre-programme common state is therefore
timebound 2026; annual permafrost release in year y is solved simultaneously
with end-of-year TAS at timebound y+1 by fixed-point iteration.
"""

from pathlib import Path
import argparse
import json
import math
import os
import time

import numpy as np
import pandas as pd
from fair import FAIR
from fair.interface import fill, initialise
from fair.io import read_properties

ROOT = Path('.')
FREP = Path(os.environ.get('FAIR_SOURCE_DIR', 'FAIR-v2.2.4'))
D = FREP / 'examples/data/calibrated_constrained_ensemble'
PAR = D / 'calibrated_constrained_parameters_calibration1.4.1.csv'
SP = D / 'species_configs_properties_calibration1.4.1.csv'
EM = D / 'extensions_1750-2500.csv'
FO = D / 'volcanic_solar.csv'
TR = ROOT / 'external_validation_net_co2_trajectory.csv'
SCENARIO = 'medium-extension'

cfg = pd.read_csv(PAR, index_col=0)
assert len(cfg) == 841
ids = np.asarray(cfg.index.astype(str))
base = np.asarray(cfg['baseline_concentration[CO2]'], float)
NCONFIG = len(cfg)
GTCO2_PER_PPM = 2.124 * (44.009 / 12.011)
GTC_TO_GTCO2 = 44.009 / 12.011

can = pd.read_csv(TR).sort_values('year').reset_index(drop=True)
required = [
    'year', 'gross_co2_gtco2', 'permafrost_co2_gtco2',
    'stored_carbon_reversal_gtco2', 'cdr_gtco2', 'net_co2_gtco2'
]
assert all(c in can.columns for c in required)
assert int(can.year.iloc[0]) == 2026 and int(can.year.iloc[-1]) == 2183
assert not can.year.duplicated().any()
mb = (can.gross_co2_gtco2 + can.permafrost_co2_gtco2
      + can.stored_carbon_reversal_gtco2 - can.cdr_gtco2)
assert float(np.max(np.abs(mb - can.net_co2_gtco2))) < 1e-9

last = can.iloc[-1]
CAN_CDR = float(can.cdr_gtco2.sum())
PEAK_CDR = float(can.cdr_gtco2.max())
LC = float(last.cdr_gtco2)
LG = float(last.gross_co2_gtco2)
LP = float(last.permafrost_co2_gtco2)
LR = float(last.stored_carbon_reversal_gtco2)
REVERSAL_FRACTION = float(can.stored_carbon_reversal_gtco2.sum() / CAN_CDR)

# Same project formulation as Hector Run 9.
PF_CARBON_PGC_PER_C = 18.0
PF_ACTIVATION_C = 1.20
PF_TAU_YR = 50.0
PF_DAMPING = 0.45
PF_TOL_GTCO2_YR = 1e-4
PF_MAX_ITER = 36


def quantiles(x, probs=(.05, .5, .95)):
    return [float(v) for v in np.quantile(np.asarray(x, float), probs)]


def make(end):
    f = FAIR(ch4_method='Thornhill2021')
    f.define_time(1750, int(end), 1)
    f.define_scenarios([SCENARIO])
    f.define_configs(cfg.index)
    species, props = read_properties(filename=SP)
    f.define_species(species, props)
    f.allocate()
    return f


def initialise_fair(f):
    f.fill_from_csv(emissions_file=EM, forcing_file=FO)
    fill(
        f.forcing,
        f.forcing.sel(specie='Volcanic')
        * cfg['forcing_scale[Volcanic]'].values.squeeze(),
        specie='Volcanic',
    )
    fill(
        f.forcing,
        f.forcing.sel(specie='Solar')
        * cfg['forcing_scale[Solar]'].values.squeeze(),
        specie='Solar',
    )
    f.fill_species_configs(SP)
    f.override_defaults(PAR)
    initialise(f.concentration, f.species_configs['baseline_concentration'])
    initialise(f.forcing, 0)
    initialise(f.temperature, 0)
    initialise(f.cumulative_emissions, 0)
    initialise(f.airborne_emissions, 0)
    initialise(f.ocean_heat_content_change, 0)


def co2(f):
    return np.asarray(
        f.concentration.loc[dict(scenario=SCENARIO, specie='CO2')], float
    )


def temp(f):
    return np.asarray(
        f.temperature.loc[dict(scenario=SCENARIO, layer=0)], float
    )


def tb_index(f, year):
    return int(np.argmin(np.abs(np.asarray(f.timebounds, float) - float(year))))


def annual_years(end):
    # FaIR define_time(..., end, 1) has timepoints through end-0.5, so the
    # final emissions year is end-1.  Canonical years 2026..2183 therefore
    # produce the concentration reported at timebound 2184.
    return np.arange(2026, int(end), dtype=int)


def component_vectors(end):
    years = annual_years(end)
    gross = np.empty(len(years), float)
    pf_exported = np.empty(len(years), float)
    reversal = np.empty(len(years), float)
    cdr = np.empty(len(years), float)
    by_year = can.set_index('year')
    for i, y in enumerate(years):
        if y <= 2183:
            r = by_year.loc[y]
            gross[i] = float(r.gross_co2_gtco2)
            pf_exported[i] = float(r.permafrost_co2_gtco2)
            reversal[i] = float(r.stored_carbon_reversal_gtco2)
            cdr[i] = float(r.cdr_gtco2)
        else:
            gross[i] = LG
            pf_exported[i] = LP
            reversal[i] = LR
            cdr[i] = LC
    return years, gross, pf_exported, reversal, cdr


def override_co2_matrix(f, years, net_gtco2):
    """Replace CO2 FFI with a year x config matrix; set AFOLU to zero."""
    years = np.asarray(years, int)
    arr = np.asarray(net_gtco2, float)
    if arr.ndim == 1:
        arr = np.repeat(arr[:, None], NCONFIG, axis=1)
    assert arr.shape == (len(years), NCONFIG), (arr.shape, len(years), NCONFIG)

    ffi = np.asarray(
        f.emissions.loc[dict(scenario=SCENARIO, specie='CO2 FFI')], float
    ).copy()
    af = np.asarray(
        f.emissions.loc[dict(scenario=SCENARIO, specie='CO2 AFOLU')], float
    ).copy()
    idx = {int(y): i for i, y in enumerate(years)}
    for ti, tp in enumerate(np.asarray(f.timepoints, float)):
        y = int(np.floor(tp))
        if y in idx:
            ffi[ti, :] = arr[idx[y], :]
            af[ti, :] = 0.0
    f.emissions.loc[dict(scenario=SCENARIO, specie='CO2 FFI')] = ffi
    f.emissions.loc[dict(scenario=SCENARIO, specie='CO2 AFOLU')] = af


def run_net(net_gtco2, end):
    years = annual_years(end)
    f = make(end)
    initialise_fair(f)
    override_co2_matrix(f, years, net_gtco2)
    f.run(progress=False)
    return f


def annual_end_tas(f, years):
    tt = temp(f)
    out = np.empty((len(years), NCONFIG), float)
    for i, y in enumerate(years):
        out[i, :] = tt[tb_index(f, y + 1), :]
    return out


def permafrost_from_tas(end_tas, common_tas):
    """Return annual GtCO2/yr release, member by member."""
    end_tas = np.asarray(end_tas, float)
    common_tas = np.asarray(common_tas, float)
    assert end_tas.ndim == 2 and end_tas.shape[1] == NCONFIG
    assert common_tas.shape == (NCONFIG,)
    running_peak = common_tas.copy()
    prev_commit = PF_CARBON_PGC_PER_C * np.maximum(0.0, common_tas - PF_ACTIVATION_C)
    pool = np.zeros(NCONFIG, float)
    rel = 1.0 - math.exp(-1.0 / PF_TAU_YR)
    out = np.zeros_like(end_tas)
    for i in range(end_tas.shape[0]):
        running_peak = np.maximum(running_peak, end_tas[i, :])
        target = PF_CARBON_PGC_PER_C * np.maximum(0.0, running_peak - PF_ACTIVATION_C)
        new_commit = np.maximum(0.0, target - prev_commit)
        pool += new_commit
        release_pgc = pool * rel
        pool -= release_pgc
        out[i, :] = release_pgc * GTC_TO_GTCO2
        prev_commit = target
    return out


def build_net(end, removal_on, reversal_off_zero, pf_matrix, extra=0.0,
              carry_extra_leakage=False):
    years, gross, _, reversal, cdr = component_vectors(end)
    pf = np.asarray(pf_matrix, float)
    if pf.ndim == 1:
        pf = np.repeat(pf[:, None], NCONFIG, axis=1)
    assert pf.shape == (len(years), NCONFIG)

    rev = reversal.copy()
    if (not removal_on) and reversal_off_zero:
        rev[:] = 0.0

    # The extra inverse CDR is applied only after the canonical trajectory,
    # beginning in 2184.  Carry its associated reversal contemporaneously at
    # the canonical cumulative reversal / cumulative CDR ratio.  This is a
    # conservative/simple correction requested by the referee and is reported
    # explicitly rather than hidden inside the canonical reversal series.
    extra_vec = np.where(years >= 2184, float(extra), 0.0)
    extra_reversal = (extra_vec * REVERSAL_FRACTION
                      if carry_extra_leakage else np.zeros_like(extra_vec))

    source = gross[:, None] + pf + rev[:, None] + extra_reversal[:, None]
    uptake = cdr[:, None] + extra_vec[:, None] if removal_on else 0.0
    return years, source - uptake, cdr + (extra_vec if removal_on else 0.0), extra_reversal


def common_state_tas(end=2400):
    """Get member-specific common pre-programme TAS at timebound 2026."""
    years, gross, pf_exp, reversal, cdr = component_vectors(end)
    # Legacy on pathway is enough because both paired arms are identical at 2026.
    net = gross + pf_exp + reversal - cdr
    f = run_net(net, end)
    return temp(f)[tb_index(f, 2026), :].copy()


def solve_dynamic(end, removal_on, reversal_off_zero, label, common_tas, seed_pf):
    years = annual_years(end)
    pf = np.asarray(seed_pf, float)
    if pf.ndim == 1:
        pf = np.repeat(pf[:, None], NCONFIG, axis=1)
    conv_rows = []
    final_f = None
    final_target = None

    for it in range(1, PF_MAX_ITER + 1):
        _, net, _, _ = build_net(
            end, removal_on, reversal_off_zero, pf,
            extra=0.0, carry_extra_leakage=False,
        )
        f = run_net(net, end)
        target_pf = permafrost_from_tas(annual_end_tas(f, years), common_tas)
        err = float(np.max(np.abs(target_pf - pf)))
        conv_rows.append({
            'case': label,
            'iteration': it,
            'max_abs_flux_change_gtco2_yr': err,
            'cumulative_pf_p05_gtco2': quantiles(np.sum(target_pf, axis=0))[0],
            'cumulative_pf_p50_gtco2': quantiles(np.sum(target_pf, axis=0))[1],
            'cumulative_pf_p95_gtco2': quantiles(np.sum(target_pf, axis=0))[2],
            'peak_pf_max_gtco2_yr': float(np.max(target_pf)),
        })
        print('PF_ITER', label, it, 'max_change', err,
              'cum_pf_p50', conv_rows[-1]['cumulative_pf_p50_gtco2'], flush=True)
        if err < PF_TOL_GTCO2_YR:
            pf = target_pf
            _, net, _, _ = build_net(
                end, removal_on, reversal_off_zero, pf,
                extra=0.0, carry_extra_leakage=False,
            )
            final_f = run_net(net, end)
            final_target = permafrost_from_tas(
                annual_end_tas(final_f, years), common_tas
            )
            break
        pf = PF_DAMPING * target_pf + (1.0 - PF_DAMPING) * pf

    if final_f is None:
        raise RuntimeError(f'dynamic permafrost did not converge: {label}')
    closure = float(np.max(np.abs(final_target - pf)))
    if closure > 5 * PF_TOL_GTCO2_YR:
        raise RuntimeError(f'final PF closure failed {label}: {closure}')
    return final_f, final_target, pd.DataFrame(conv_rows), closure


def paired_table(fon, foff, cdr_schedule, label):
    rows = []
    aa, bb = co2(fon), co2(foff)
    # cdr_schedule has annual years 2026..end-1.
    years = annual_years(int(round(float(np.asarray(fon.timebounds)[-1]))))
    cum = np.cumsum(np.asarray(cdr_schedule, float))
    cumulative_by_tb = {int(y + 1): float(cum[i]) for i, y in enumerate(years)}
    for i, t in enumerate(np.asarray(fon.timebounds, float)):
        yy = int(round(float(t)))
        cumulative = cumulative_by_tb.get(yy, 0.0)
        delta = bb[i, :] - aa[i, :]
        frac = (delta * GTCO2_PER_PPM / cumulative
                if cumulative > 0 else np.full(NCONFIG, np.nan))
        dq = quantiles(delta)
        fq = quantiles(frac) if cumulative > 0 else [math.nan] * 3
        rows.append({
            'experiment': label,
            'timebound_year': float(t),
            'cumulative_cdr_gtco2': cumulative,
            'delta_co2_p05_ppm': dq[0],
            'delta_co2_p50_ppm': dq[1],
            'delta_co2_p95_ppm': dq[2],
            'fraction_p05': fq[0],
            'fraction_p50': fq[1],
            'fraction_p95': fq[2],
        })
    return pd.DataFrame(rows)


def selected_rows(df):
    out = []
    for y in [2040, 2100, 2156, 2184, 2200, 2300, 2400]:
        i = int(np.argmin(np.abs(df.timebound_year.to_numpy(float) - y)))
        out.append(df.iloc[i])
    return pd.DataFrame(out)


def run_feedback():
    end = 2400
    years, gross, pf_exp, reversal, cdr = component_vectors(end)
    seed = np.repeat(pf_exp[:, None], NCONFIG, axis=1)

    # Legacy strict pair: exact Run 8 construction, used as regression lock.
    legacy_on_net = gross + pf_exp + reversal - cdr
    legacy_off_net = gross + pf_exp + reversal
    print('RUN legacy strict pair', flush=True)
    legacy_on = run_net(legacy_on_net, end)
    legacy_off = run_net(legacy_off_net, end)

    common_tas = temp(legacy_on)[tb_index(legacy_on, 2026), :].copy()
    if float(np.max(np.abs(common_tas - temp(legacy_off)[tb_index(legacy_off, 2026), :]))) > 1e-12:
        raise RuntimeError('legacy common-state TAS mismatch')
    common_co2_err = float(np.max(np.abs(
        co2(legacy_on)[tb_index(legacy_on, 2026), :]
        - co2(legacy_off)[tb_index(legacy_off, 2026), :]
    )))

    # Reversal correction only, with prescribed permafrost retained.
    print('RUN reversal-zero-off only', flush=True)
    no_rev_off_net = gross + pf_exp  # no stored reversal, no CDR
    no_rev_off = run_net(no_rev_off_net, end)

    # Dynamic PF cases.  Dynamic-on is shared by both dynamic comparisons.
    print('RUN dynamic permafrost on arm', flush=True)
    dyn_on, pf_on, conv_on, close_on = solve_dynamic(
        end, True, False, 'dynamic_on', common_tas, seed
    )
    print('RUN dynamic permafrost off arm with legacy reversal', flush=True)
    dyn_off_same_rev, pf_off_same, conv_off_same, close_off_same = solve_dynamic(
        end, False, False, 'dynamic_off_same_reversal', common_tas, seed
    )
    print('RUN dynamic permafrost off arm with reversal zero', flush=True)
    dyn_off_no_rev, pf_off_zero, conv_off_zero, close_off_zero = solve_dynamic(
        end, False, True, 'dynamic_off_no_reversal', common_tas, seed
    )

    tables = []
    tables.append(paired_table(legacy_on, legacy_off, cdr, 'legacy_strict_pair'))
    tables.append(paired_table(legacy_on, no_rev_off, cdr, 'reversal_zero_off_only'))
    tables.append(paired_table(dyn_on, dyn_off_same_rev, cdr, 'dynamic_permafrost_only'))
    tables.append(paired_table(dyn_on, dyn_off_no_rev, cdr, 'combined_programme_effect'))
    pair = pd.concat(tables, ignore_index=True)
    pair.to_csv('fair_feedback_sensitivity_pair.csv', index=False)
    pd.concat([selected_rows(x) for x in tables], ignore_index=True).to_csv(
        'fair_feedback_sensitivity_selected.csv', index=False
    )

    # Store member arrays at selected years to allow independent quantile audit.
    sel_years = np.array([2026, 2040, 2100, 2156, 2184, 2200, 2300, 2400])
    cases = {
        'legacy_on': legacy_on, 'legacy_off': legacy_off,
        'reversal_zero_off': no_rev_off,
        'dynamic_on': dyn_on,
        'dynamic_off_same_reversal': dyn_off_same_rev,
        'dynamic_off_no_reversal': dyn_off_no_rev,
    }
    npz = {'selected_years': sel_years, 'config_ids': ids.astype('U')}
    for name, f in cases.items():
        npz[f'{name}_co2'] = np.vstack([co2(f)[tb_index(f, y), :] for y in sel_years])
        npz[f'{name}_tas'] = np.vstack([temp(f)[tb_index(f, y), :] for y in sel_years])
    np.savez_compressed('fair_feedback_sensitivity_members.npz', **npz)

    # Permafrost audit is long-form only for selected years plus member quantiles,
    # avoiding a huge CSV while retaining the full selected-member arrays above.
    pf_rows = []
    for label, pf in [
        ('dynamic_on', pf_on),
        ('dynamic_off_same_reversal', pf_off_same),
        ('dynamic_off_no_reversal', pf_off_zero),
    ]:
        for y in [2026, 2040, 2100, 2156, 2183, 2200, 2300, 2399]:
            j = int(np.where(years == y)[0][0])
            qq = quantiles(pf[j, :])
            pf_rows.append({
                'case': label, 'emissions_year': y,
                'pf_p05_gtco2_yr': qq[0], 'pf_p50_gtco2_yr': qq[1],
                'pf_p95_gtco2_yr': qq[2],
            })
    pd.DataFrame(pf_rows).to_csv('fair_feedback_sensitivity_permafrost.csv', index=False)
    pd.concat([conv_on, conv_off_same, conv_off_zero], ignore_index=True).to_csv(
        'fair_feedback_sensitivity_convergence.csv', index=False
    )

    # Explicit forcing/path audit at selected programme years.
    audit_rows = []
    for y in [2026, 2040, 2100, 2156, 2183, 2200, 2300, 2399]:
        j = int(np.where(years == y)[0][0])
        for label, removal_on, rev_zero, pf in [
            ('legacy_on', True, False, np.repeat(pf_exp[:,None], NCONFIG, axis=1)),
            ('legacy_off', False, False, np.repeat(pf_exp[:,None], NCONFIG, axis=1)),
            ('reversal_zero_off', False, True, np.repeat(pf_exp[:,None], NCONFIG, axis=1)),
            ('dynamic_on', True, False, pf_on),
            ('dynamic_off_same_reversal', False, False, pf_off_same),
            ('dynamic_off_no_reversal', False, True, pf_off_zero),
        ]:
            rev = 0.0 if ((not removal_on) and rev_zero) else reversal[j]
            qq = quantiles(pf[j, :])
            audit_rows.append({
                'case': label, 'emissions_year': y,
                'gross_gtco2_yr': gross[j], 'reversal_gtco2_yr': rev,
                'cdr_gtco2_yr': cdr[j] if removal_on else 0.0,
                'permafrost_p05_gtco2_yr': qq[0],
                'permafrost_p50_gtco2_yr': qq[1],
                'permafrost_p95_gtco2_yr': qq[2],
            })
    pd.DataFrame(audit_rows).to_csv('fair_feedback_sensitivity_forcing_audit.csv', index=False)

    # Regression and headline summary.
    def get(df, exp, year):
        z = df[df.experiment == exp]
        return z.iloc[int(np.argmin(np.abs(z.timebound_year.to_numpy(float) - year)))]

    legacy_2300 = get(pair, 'legacy_strict_pair', 2300)
    combined = {str(y): get(pair, 'combined_programme_effect', y).to_dict()
                for y in [2156, 2184, 2200, 2300, 2400]}
    summary = {
        'model': 'FaIR 2.2.4',
        'calibration': 'fair-calibrate 1.4.1',
        'configs': 841,
        'gtco2_per_ppm': GTCO2_PER_PPM,
        'canonical_cdr_2026_2183_gtco2': CAN_CDR,
        'canonical_reversal_2026_2183_gtco2': float(can.stored_carbon_reversal_gtco2.sum()),
        'canonical_reversal_fraction_of_cdr': REVERSAL_FRACTION,
        'common_state_co2_maxabs_ppm_at_2026': common_co2_err,
        'common_state_tas_definition': 'member-specific FaIR timebound 2026 before programme emissions take effect',
        'permafrost_formulation': {
            'carbon_pgc_per_degree_c': PF_CARBON_PGC_PER_C,
            'activation_degc': PF_ACTIVATION_C,
            'tau_years': PF_TAU_YR,
            'running_peak_commitment': True,
            'historical_commitment_reinjected': False,
            'fixed_point_damping': PF_DAMPING,
            'tolerance_gtco2_per_year': PF_TOL_GTCO2_YR,
        },
        'fixed_point_closure_error_gtco2_yr': {
            'dynamic_on': close_on,
            'dynamic_off_same_reversal': close_off_same,
            'dynamic_off_no_reversal': close_off_zero,
        },
        'legacy_2300': legacy_2300.to_dict(),
        'combined_selected': combined,
        'interpretation': (
            'combined_programme_effect is the corrected paired attribution: stored-carbon '
            'reversal is programme-contingent and therefore zero in removal-off; permafrost '
            'is generated from each arm/configuration own TAS path.'
        ),
    }
    Path('fair_feedback_sensitivity_summary.json').write_text(json.dumps(summary, indent=2))
    print('FAIR_FEEDBACK_COMPLETE', flush=True)


def endpoint_stats(f, target):
    c = co2(f)[tb_index(f, target), :]
    t = temp(f)[tb_index(f, target), :]
    d = c - base
    return {
        'c': c, 't': t, 'd': d,
        'co2_q': quantiles(c),
        'delta_baseline_q': quantiles(d),
        'co2_p90': float(np.quantile(c, .90)),
        'delta_baseline_p90': float(np.quantile(d, .90)),
        'fraction_le_280': float(np.mean(c <= 280.0)),
        'fraction_le_own_baseline': float(np.mean(c <= base)),
        'temp_q': quantiles(t),
    }


def run_inverse(target):
    target = int(target)
    if target not in (2200, 2300, 2400):
        raise ValueError('target must be 2200, 2300, or 2400')
    years, gross, pf_exp, reversal, cdr = component_vectors(target)
    pf = np.repeat(pf_exp[:, None], NCONFIG, axis=1)
    cache = {}

    def evaluate(rate):
        key = round(float(rate), 10)
        if key not in cache:
            t0 = time.time()
            _, net, _, extra_rev = build_net(
                target, True, False, pf,
                extra=float(rate), carry_extra_leakage=True,
            )
            f = run_net(net, target)
            st = endpoint_stats(f, target)
            if np.any(~np.isfinite(st['c'])) or np.any(st['c'] <= 0):
                raise RuntimeError(f'nonphysical endpoint target={target}, rate={rate}')
            st['extra_reversal_rate'] = float(rate) * REVERSAL_FRACTION
            st['extra_reversal_cumulative'] = float(np.sum(extra_rev))
            cache[key] = st
            print('POLICY_EVAL', target, rate,
                  'p50', st['co2_q'][1], 'p90', st['co2_p90'],
                  'f280', st['fraction_le_280'],
                  'seconds', round(time.time()-t0, 1), flush=True)
        return cache[key]

    criteria = {
        'absolute_280_median': lambda s: float(np.median(s['c']) - 280.0),
        'member_relative_baseline_median': lambda s: float(np.median(s['d'])),
        'absolute_280_p90': lambda s: float(np.quantile(s['c'], .90) - 280.0),
        'member_relative_baseline_p90': lambda s: float(np.quantile(s['d'], .90)),
    }

    def solve(name, metric):
        lo = 0.0
        flo = metric(evaluate(lo))
        if flo <= 0:
            return 0.0, 'already_met', 0.0, 0.0
        hi = 1.0
        fhi = metric(evaluate(hi))
        while fhi > 0 and hi < 1024:
            lo, flo = hi, fhi
            hi *= 2.0
            fhi = metric(evaluate(hi))
        if fhi > 0:
            return None, 'not_bracketed', lo, hi
        # 10 rounds gives <= ~0.03 Gt/yr resolution even in a 32 Gt bracket.
        for _ in range(10):
            mid = (lo + hi) / 2.0
            fm = metric(evaluate(mid))
            if fm > 0:
                lo, flo = mid, fm
            else:
                hi, fhi = mid, fm
        return hi, 'solved', lo, hi

    rows = []
    members = []
    n_years = target - 2184
    for name, metric in criteria.items():
        rate, status, lo, hi = solve(name, metric)
        row = {
            'target_year': target, 'criterion': name, 'status': status,
            'extra_cdr_rate_gtco2_yr': rate, 'bracket_low': lo, 'bracket_high': hi,
            'canonical_extra_reversal_fraction': REVERSAL_FRACTION,
        }
        if rate is not None:
            st = evaluate(rate)
            row.update({
                'number_of_extra_cdr_years': n_years,
                'cumulative_extra_cdr_gtco2': float(rate) * n_years,
                'extra_reversal_rate_gtco2_yr': st['extra_reversal_rate'],
                'cumulative_extra_reversal_gtco2': st['extra_reversal_cumulative'],
                'net_extra_removal_after_extra_reversal_gtco2':
                    float(rate) * n_years - st['extra_reversal_cumulative'],
                'total_post2183_cdr_rate_gtco2_yr': LC + float(rate),
                'co2_p05': st['co2_q'][0], 'co2_p50': st['co2_q'][1],
                'co2_p90': st['co2_p90'], 'co2_p95': st['co2_q'][2],
                'co2_minus_baseline_p05': st['delta_baseline_q'][0],
                'co2_minus_baseline_p50': st['delta_baseline_q'][1],
                'co2_minus_baseline_p90': st['delta_baseline_p90'],
                'co2_minus_baseline_p95': st['delta_baseline_q'][2],
                'fraction_le_280': st['fraction_le_280'],
                'fraction_le_own_baseline': st['fraction_le_own_baseline'],
                'temp_p05': st['temp_q'][0], 'temp_p50': st['temp_q'][1],
                'temp_p95': st['temp_q'][2],
                'historical_program_peak_cdr_gtco2_yr': PEAK_CDR,
                'exceeds_historical_program_peak': bool(LC + float(rate) > PEAK_CDR),
            })
            for cid, b, c, d, t in zip(ids, base, st['c'], st['d'], st['t']):
                members.append({
                    'target_year': target, 'criterion': name,
                    'extra_cdr_rate_gtco2_yr': rate,
                    'config': cid, 'baseline_co2_ppm': b,
                    'co2_ppm': c, 'co2_minus_baseline_ppm': d,
                    'temperature_degc': t,
                    'le_280': bool(c <= 280.0),
                    'le_own_baseline': bool(c <= b),
                })
        rows.append(row)
        print('POLICY_SOLUTION', target, name, status, rate, lo, hi, flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(f'fair_policy_inverse_{target}_summary.csv', index=False)
    pd.DataFrame(members).to_csv(f'fair_policy_inverse_{target}_members.csv', index=False)
    summary = {
        'model': 'FaIR 2.2.4', 'calibration': 'fair-calibrate 1.4.1',
        'configs': 841, 'target_year': target,
        'canonical_cdr_2026_2183_gtco2': CAN_CDR,
        'canonical_terminal_cdr_gtco2_yr': LC,
        'canonical_reversal_fraction_of_cdr': REVERSAL_FRACTION,
        'extra_leakage_treatment': (
            'Additional post-2183 CDR carries contemporaneous stored-carbon reversal '
            'equal to the canonical cumulative stored-reversal / cumulative CDR fraction.'
        ),
        'criteria': list(criteria),
        'inverse': rows,
    }
    Path(f'fair_policy_inverse_{target}_summary.json').write_text(json.dumps(summary, indent=2))
    print('FAIR_POLICY_INVERSE_COMPLETE', target, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=['feedback', 'inverse'], required=True)
    ap.add_argument('--target', type=int, choices=[2200, 2300, 2400])
    args = ap.parse_args()
    if args.mode == 'feedback':
        if args.target is not None:
            ap.error('--target is only valid with --mode inverse')
        run_feedback()
    else:
        if args.target is None:
            ap.error('--target is required with --mode inverse')
        run_inverse(args.target)


if __name__ == '__main__':
    main()
