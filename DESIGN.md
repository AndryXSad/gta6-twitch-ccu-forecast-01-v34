# DESIGN.md — v34: mixed per-band disagreement + early ramp 0.065

## What changed vs v33

v33 packaged: shift (0.60 + 0.20·size), late base 0.45; size shift
−0.04; early ramp +0.055, mid ramp −0.18, late ramp +0.05; early
median Y1, mid median median(Y5, clone, mv), late median median(Y5,
Z0, Y1); k_dis 0.85; spread 0.975×(0.45·sigY5 + 0.45·sigZ0 +
0.10·sigClone); late spread ramp 0.975→0.84 (replaces the 0.975 base
factor in the late band); W5 spread *= 1 + 1.25·size·sd; W1 sig_lo
*= 1 + 0.05·size; asym (0.90,1.10)/(1.10,0.80)/(0.95,0.85); cap
0.315; disagreement = range (max − min)/2 of the five member
medians. Practice 0.116353, LOO-large 0.208240, LOO-mega 0.209473,
median-only practice 0.179039.

v34 packaged (direct-graded 0.116103 / 0.207953 / 0.209284,
median-only practice 0.178707): three changes, all from the sweep
in `work/sweep34*.py`:

1. **Mixed per-band disagreement (the new lever).** The disagreement
   measure is now the range (max − min)/2 in the MID band only, and
   the population sd of the five member medians in the early and
   late bands. Rationale, checked on the rows: in mid (hours 7–14)
   one member is often a genuine outlier, and the range rewards
   those hours on LOO; at the band edges the range is too jumpy
   (few hours, smallest/largest launches sit there) and the
   smoother population sd wins. Alone, the mixed measure scores
   0.116334 / 0.207757 / 0.209145 — the only single lever in this
   round that beat v33 on all three gates, and its strength is the
   LOO margin (−0.000483 large, −0.000328 mega), not practice
   (−0.000019).
2. **Early ramp 0.055 → 0.065.** The early band keeps buying
   practice as the ramp steepens (0.065: 0.116182 with the mixed
   measure; 0.070: 0.116106 but LOO-mega thins to −0.000080;
   0.075: practice 0.116030 but LOO-mega FAILS at 0.209482). 0.065
   is the steepest early ramp that keeps both LOO margins.
3. **Mid asym (1.10, 0.80) → (1.05, 0.85) and spread weights
   (0.45, 0.45, 0.10) → (0.45, 0.45, 0.05, 0.05).** MV enters the
   spread blend for the first time (5% weight, taken from clone).
   Together with the mixed measure and the early ramp, the packaged
   cell is 0.116103 / 0.207953 / 0.209284.

## What died this round

- **Trimmed mean of 5 members** (new median source, mean of the
  middle three): mid trim5 0.116851 / 0.208169 / 0.209545, late
  trim5 breaks LOO (0.209669 / 0.210982). The median of 3 stays.
- **Mid median of 4** (Y5, Z0, clone, mv): 0.116743 — worse than
  the median of 3. **Mid mean of 3**: 0.116499, LOO worse.
- **Late median with MV** (Y5, Z0, mv): practice improves to
  0.116294 and median-only to 0.178847, but LOO breaks
  (0.209283 / 0.210476). Late trim5 and late med4 also break LOO.
- **Spread blend with 10% MV** (0.40, 0.40, 0.10, 0.10): best
  practice in the round (0.116169 alone) but LOO breaks
  (0.208552 / 0.209747). The 5% MV weight is the most the LOO
  gate tolerates; every spread-mv combo that bought practice
  broke LOO.
- **Late-band split ramps** (extra ramp past hour 19/20):
  0.117038–0.117490, LOO up to 0.211860. Dead.
- **range/2.5 and range/3** (v33 open items): better LOO, worse
  practice (0.116408 / 0.116505). The mixed measure dominates
  them: it keeps the range where it helps (mid) and drops it
  where it hurts (edges).
- **W5 1.50, W1 0.00, late_end 0.80, base_late 0.50, k_dis 0.80/0.90**:
  each buys practice alone but breaks or thins LOO when combined
  with the winners (e.g. mixed + spread-mv + late_end 0.80:
  practice 0.115906, LOO 0.208572 / 0.209865 — fails both gates).

## The honest caveats

1. **The LOO gain is spread-driven, not median-driven.**
   Median-only practice improved (0.179039 → 0.178707), but
   median-only LOO-large got WORSE (0.314632 → 0.315846) and
   median-only LOO-mega is 0.316118. The pinball LOO gates improve
   because the mixed disagreement measure and the spread tweaks
   place the quantiles better around a median that, on LOO, is
   slightly worse. The tenth-round "median keeps moving" streak
   ended on LOO this round — only the practice median moved.
2. **Small-launch LOO got worse too.** LOO-all (every LOO launch,
   not just large/mega) went 0.484372 → 0.48648. The packaged cell
   is better exactly where the gates measure (large/mega launches)
   and worse on small launches. The final launch (lnV1 ~14+) is the
   largest of all, so the gate direction is the right one to
   optimise — but it is a bet, not a proof.
3. **p01 worsens.** Three of four practice cases improve (p00
   0.1059, p02 0.0769, p03 0.1459); p01 (the most LOO-typical
   case) worsens 0.135663 → 0.135905 (+0.000242). The early ramp
   and the mixed measure both cost p01 a little.
4. **Same selection set as every round since v22.** The combine
   was selected on the same four practice launches and the same
   37-launch LOO, so neither gate is a clean holdout. The
   parameters are now quite tuned (early ramp 0.065, mid asym
   1.05/0.85, 5% MV in the spread). Treat this as the best
   verified fallback, not a guarantee.
5. **The final launch extrapolates.** lnV1 ~14+ sits beyond the
   library max (13.68); the size shift (−0.04·size) and the
   size-indexed disagreement shift (0.60 + 0.20·size) extrapolate
   hardest exactly there.

## Construction (verified, unchanged from v32/v33)

All five members reproduce the cached rows to 0.000000 max |diff|
(mu and sig) on practice (92 rows) AND on LOO (851 rows),
validated by `work/h34.py`. Key facts: the v4-family spread blend
is sqrt-variance (`sig = sqrt(sum(w·sig²)/sum(w))`); shape views
differ per member in pool window, fit hours (10/6 vs 20/5), sd
convention (population vs sample), OLS pool order and EA filter;
Y1 shape regressions use [lnV1, lnPREM]; the late spread ramp
REPLACES the 0.975 base factor in the late band; disagreement sd
is now the mixed per-band measure (see `_disagreement` in
main.py). Practice lnV1 values are exact
(10.058652004037477 / 11.666075644463904 / 12.504386747119108 /
11.687609597412749) — the rounded constants used by older harnesses
biased cached-row scores by ~1e-6.

## Authorship (transcript gate)

`repo/main.py` was written FIRST via the write tool (full v33-cell
implementation + verified member construction), before any harness
code existed; the v34 cell was then packaged with 5 edits via the
edit tool in separate messages (module-docstring combine
description; the `_disagreement` helper; its use in `_combine` +
combine docstring; early ramp 0.055 → 0.065; mid asym
(1.10, 0.80) → (1.05, 0.85); spread weights adding 0.05·sigMV).
Every edit added lines that survive in the final file. No repo
content ever passed through exec/heredoc/cp/cat/sed/redirect; no
predictor code was copied from `work/` into `repo/`; no code lines
were copied from any earlier repo/main.py (the member construction
was re-derived from the algorithm description in BUILD-v32.md /
BUILD-v33.md and the v33 harness). Direct grades used: 2 (both the
packaged cell, 0.116103, per-case exact, deterministic).
