# Referee Report — V60.0 Crosscheck Package and Run 9 Sensitivity Design

**Packages:** `Planetary_Restoration_Model_V60_0_COMPLETE_CROSSCHECK` (226 files) and `Hector_Feedback_Sensitivity_Run9_FLAT` (3 files, scripts only)
**Recommendation:** **Revision, and for the first time a short one.** Three of the four blocking gates are closed and independently verified. One correction stands between this and a defensible headline number.

---

## 1. Independent verification of the new results

I recomputed the reported quantiles from the raw member arrays rather than accepting the summary files.

**The common-state audit holds exactly.** From `fair_paired_members.npz` directly: `max|CO₂_on − CO₂_off|` at 2025 and 2026 is 0.000e+00 ppm, and the same for temperature at 2026. Not "small" — zero.

**Every reported quantile reproduces to four decimals.** Recomputing member-wise `(off − on) × 7.782459 / cumulative_CDR` and then taking percentiles across the 841 configurations:

| Year | p05 reported / recomputed | p50 reported / recomputed | p95 reported / recomputed |
|---|---|---|---|
| 2040 | 0.6552 / 0.6552 | 0.7127 / 0.7127 | 0.7522 / 0.7522 |
| 2100 | 0.3965 / 0.3965 | 0.4720 / 0.4720 | 0.5257 / 0.5257 |
| 2156 | 0.3351 / 0.3351 | 0.4075 / 0.4075 | 0.4560 / 0.4560 |
| 2184 | 0.2984 / 0.2984 | 0.3699 / 0.3699 | 0.4200 / 0.4200 |
| 2200 | 0.2865 / 0.2865 | 0.3573 / 0.3573 | 0.4095 / 0.4095 |
| 2300 | 0.2545 / 0.2545 | 0.3083 / 0.3083 | 0.3690 / 0.3690 |
| 2400 | 0.2457 / 0.2457 | 0.2835 / 0.2835 | 0.3410 / 0.3410 |

**The quantile-differencing error is fixed and the difference is measurable.** At 2300, differencing the p05/p50/p95 concentration curves gives [0.2684, 0.3087, 0.3566], width 0.088. The correct member-wise calculation gives [0.2545, 0.3083, 0.3690], width 0.115. The medians agree to 0.0004; the old method compresses the tails by about 30%. The retirement of the [0.5678, 0.5872] envelope is correct and the replacement is defensible.

**The abstract's numbers check out against the archive.** 95.71 ppm at 2184, 744.9 GtCO₂ (95.71058764829235 × 7.782459079177421 = 744.85), fraction 0.370 with p05–p95 of 0.298–0.420, and the six inverse rates 31.46 / 9.46 / 6.99 and 33.21 / 9.87 / 7.24 GtCO₂ yr⁻¹ all match the summary JSONs exactly.

**The inverse arithmetic is self-consistent.** 26.1875 × 16 = 419.0; 4.1875 × 116 = 485.75; 1.7109375 × 216 = 369.5625; and 26.1875 + 5.275424296 = 31.4629 for the total post-2183 rate.

**Provenance is auditable.** Raw GitHub Actions ZIPs are archived unmodified alongside the processed outputs, with run ID 34255177367 and commit 4aed7706 named in the guide.

---

## 2. What the package now establishes

**The non-CO₂ background is not the explanation.** The neutral-future-non-CO₂ sensitivity moves the median response fraction by only −0.0093 to −0.0121 across 2100–2400 (0.4720 → 0.4626 at 2100; 0.3083 → 0.2962 at 2300). That is a genuinely useful negative result: the attribution is robust to a strongly altered background, and the remaining Hector/FaIR gap cannot simply be attributed to forcing.

**The Hector/FaIR comparison is now quantitative.** Hector's paired fractions (0.2808 at 2200, 0.2437 at 2300) sit just below FaIR's p05 (0.2865, 0.2545) — a gap of 0.006–0.011, which is roughly the size of the neutral-background effect itself. That is worth stating explicitly: the divergence is small and of the same order as a known systematic, which makes it plausible but unconfirmed that a neutral-background Hector rerun would close it. The package correctly declines to claim the divergence has been physically explained.

**The inverse solve is the paper's contribution, and it delivers.** Pulling the 280 ppm median date from 2400 to 2200 raises the required post-2183 removal rate from 6.99 to 31.46 GtCO₂ yr⁻¹ — a 4.5× increase, and more than double the programme's own 15.2 GtCO₂ yr⁻¹ peak. That is a quantified deadline-cost curve for atmospheric restoration, it is exactly what `NOVELTY_BENCHMARK_V59_1.md` identified as the novelty discriminator, and it is a genuinely publishable result.

**The manuscript is now the right paper.** 3,789 words, titled "Financing planetary restoration under state-dependent carbon-sink reversal," with the stablecoin gone entirely from the text and an introduction that states plainly that sink reversal is not a new discovery and identifies the contribution as the pathway-specific removal burden. Compared with the V53.1 manuscript this is a different and far better document.

---

## 3. The one remaining scientific objection, and it applies to the headline number

The Run 9 sensitivity script is well designed — better than what I asked for. It factorially separates the two corrections into `reversal_zero_off_only`, `dynamic_permafrost_only` and `combined_programme_effect`, keeps `legacy_strict_pair` as a reference, solves the temperature–permafrost coupling by damped fixed-point iteration with an independent final closure check, and the verifier includes a regression lock proving the legacy arm reproduces the first experiment to 1e-9. That last control is the right instinct.

**But Run 9 has not been executed, and it exists only for Hector.** The headline number in the abstract is FaIR's 0.370, and `run_fair_science_experiments.py` builds the removal-off arm at lines 64–67 as `gross_co2_gtco2 + permafrost_co2_gtco2 + stored_carbon_reversal_gtco2`, with the extension holding `off[y] = LG + LP + LR`. That is the identical construction I objected to in Hector:

- **Stored-carbon reversal is charged to a world that never ran the programme.** Cumulative reversal to 2400 is roughly 271 GtCO₂ (83.4 through 2183 plus 217 years at the terminal 0.866 GtCO₂ yr⁻¹). Removing it from the counterfactual would lower removal-off CO₂ by order 7–10 ppm against deltas of 104–115 ppm — a 6–9% overstatement of the attributed benefit.
- **Permafrost is prescribed identically across arms with different warming.** Run 9's own formulation suggests this is second-order: with a running-peak commitment at 18 PgC °C⁻¹ above 1.20 °C, the between-arm difference in cumulative permafrost is of order 8 GtCO₂, about 1% of the delta.

Net expected effect: the response fraction should fall by roughly 5–9%, so 0.308 at 2300 would land near 0.285. That is a prediction you can check against the run, and if the combined case lands far outside it, something in the implementation needs a second look.

**Two things follow.** First, Run 9 must be ported to FaIR — correcting Hector alone leaves the abstract's number uncorrected. Second, and importantly: **the inverse solve is immune to this objection**, because it is a forward solve on the removal-on arm with no counterfactual. The paper's novel contribution does not move.

---

## 4. Secondary points a referee will raise

**The inverse solve does not carry its own leakage.** `paths(extra=...)` holds reversal at the terminal `LR` regardless of how much additional CDR is imposed. Adding 419 GtCO₂ (2200) or 370 GtCO₂ (2400) of extra storage should generate additional leakage at the programme's own 4.14% reversal rate, roughly 15–17 GtCO₂. This biases the required rates downward — restoration is slightly harder than reported. Second-order, but it should be either corrected or stated.

**The cumulative extra CDR is non-monotonic in the deadline: 419 GtCO₂ for 2200, 485.75 for 2300, 369.56 for 2400.** A referee will stop at this. The maximum near 2300 presumably reflects a trade-off between the length of the window over which sinks re-release and the time available for natural drawdown, but the paper does not explain it and must. It is also a more interesting result than the monotone rate curve, because it implies an intermediate deadline is the most expensive one in cumulative terms.

**Solving for the ensemble median means half the ensemble misses.** `fraction_le_280` is ≈0.50 at every solved rate, by construction, and `fraction_le_own_baseline` is 0.05 at 2200 under the absolute-280 criterion. The abstract's wording ("place the ensemble median at or below 280 ppm") is accurate, but for a restoration target a coin-flip criterion is weak. Report the p90 solve as well; it is the policy-relevant number and it will be substantially larger.

**The affiliations block still contains an unresolved editorial note.** The front matter reads "author-supplied educational associations" and "Formal institutional affiliations should be confirmed by each author before submission," with two affiliations marked historical. This was flagged several packages ago. It is a to-do left in submitted text and will draw an editorial query.

**`__pycache__` is back.** Nine `.pyc` files (cpython-313), hashed into `SHA256SUMS.txt` as deliverable content. This was closed in V61.2 and has reopened. Everything else verifies: 225/225 hashes, tests exit 0.

---

## 5. Assessment

Six packages in a row moved the container and not the contents. The last two moved the science, and this one moved it decisively: the paired FaIR attribution across 841 configurations, member-wise uncertainty replacing the quantile-differencing artifact, the neutral-background control, and three inverse solves are collectively the work that was blocking submission from the first review onward. The results were verified here from raw member arrays, not taken on trust, and they hold.

The path to submission:

1. Port Run 9 to FaIR and execute both. Report the corrected response fraction as the headline; expect roughly 5–9% lower.
2. Either carry leakage on the inverse-solved extra CDR or state that it is neglected and bounded.
3. Explain the non-monotonic cumulative requirement, and add a p90 inverse solve alongside the median.
4. Run the neutral-background Hector comparator, or state in the manuscript that the residual FaIR/Hector gap is of the same order as the background sensitivity and is therefore not yet resolved.
5. Resolve the affiliations, and drop `__pycache__` from the build.

None of that is new research. It is one script port, one rerun, and two paragraphs. The paper is close.
