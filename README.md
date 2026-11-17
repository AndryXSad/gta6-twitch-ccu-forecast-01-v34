# GTA VI launch-day Twitch CCU forecaster — v34

[![make.fun · gta6-twitch-ccu-forecast-01](https://make.fun/badge/gta6-twitch-ccu-forecast-01.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)
[![rank](https://make.fun/badge/gta6-twitch-ccu-forecast-01/AndryXSad/gta6-twitch-ccu-forecast-01-v34/rank.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01/leaderboard?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)
[![reward](https://make.fun/badge/gta6-twitch-ccu-forecast-01/AndryXSad/gta6-twitch-ccu-forecast-01-v34/reward.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)

A make.fun `gta6-twitch-ccu-forecast-01` entry. It forecasts, for each of
hours 2–24 of Grand Theft Auto VI's release day (2026-11-19 UTC), five
quantiles (5%, 25%, 50%, 75%, 95%) of the concurrent Twitch viewers in the
"Grand Theft Auto VI" category, from the first hour's print and the
launch-day curves of past big game launches shipped inside each case.

## How it works

`main.py` (Python 3.13, standard library only, no network) builds five
member predictors in log space from the 37 comparable past launches:

1. **Y5** — v4-family blend, tier EA filter (drop tier == ea), shape
   regressions on `[lnV1, g12]`. `g12` is the "unlock momentum" feature:
   lnV1 minus log viewers 12 hours before release day.
2. **Z0** — v4-family blend, naive EA filter (drop launches whose unlock
   note contains "early access"), shape regressions on `[lnV1, g12]`.
3. **Y1** — v4-family blend, naive EA filter, shape regressions on
   `[lnV1, lnPREM]`.
4. **Clone** — v9-clone baseline: whole-library shape views, sample-sd
   spreads, weighted-mean blend (0.6 OLS / 0.2 console-shape / 0.2
   all-shape) with spread 0.81× the weighted mean of the view spreads.
5. **MV** — per-hour median across four views (per-hour OLS, console
   shape, all shape, clone forecast); spread = mean of the four views'
   spreads. It is the mid-band median source together with Y5 and
   clone, and the fifth member in the disagreement measure.

A v4-family member mixes three views per hour (0.5 OLS, 0.2
console-only shape, 0.3 all-launch shape): the median is the weighted
mean, the spread is the **sqrt-variance blend**
`sqrt(sum(w * sig^2) / sum(w))`. The shape views differ per member
(fit_hours 10 vs 20, min_fits 6 vs 5, windowed vs whole-library
pools, population vs sample sd, OLS pool counted before vs after the
hour filter). All five members reproduce the cached rows to 0.000000
on practice AND on LOO.

The **log-median is chosen per hour band** — the median is the bottleneck:

| band | hours | median source | ramp |
| ---- | ----- | ------------- | ---- |
| early | 2–6 | Y1 log-median | 0.00 → +0.065 |
| mid | 7–14 | median of Y5 / clone / MV log-medians | 0.00 → −0.18 |
| late | 15–24 | median of Y5 / Z0 / Y1 log-medians | 0.00 → +0.05 |

On top of that, a **size-dependent shift**: above lnV1 11.5 the median is
pulled down by 0.04 per log unit of size, in every band, and a **deeper,
size-indexed disagreement shift** `mu += (base + 0.20 * size) * sd`
(base 0.60 early/mid, 0.45 late). The disagreement measure is **mixed
per band** — new in v34: in the mid band it is the range
`(max − min) / 2` of the five member medians (hours where one member is
an outlier get a bigger shift, which the LOO gate rewards there); in
the early and late bands it is the **population sd** of the five member
medians (the range is too jumpy at the band edges). A pure range
measure everywhere scores 0.116353 (v33); the mixed measure alone
scores 0.116334 / 0.207757 / 0.209145 — the LOO gain is where it earns
its keep.

The **log-spread** is 0.975× the weighted mean of the member spreads
(45% Y5 + 45% Z0 + 5% clone + 5% MV — MV enters the spread blend for
the first time in v34), widened by a **disagreement factor**
`1 + 0.85 * sd(member medians)`, shaped per band with separate
lower/upper scales: early (0.90, 1.10), mid (1.05, 0.85), late
(0.95, 0.85), narrowed across the late band (spread scale 0.975 → 0.84
from hour 15 to hour 24 — this REPLACES the 0.975 base factor in the
late band), and finally **capped at 0.315**. Two size levers from v31:
**W1** grows only the lower spread with size (`slo *= 1 + 0.05 * size`)
and **W5** lets size amplify the disagreement effect on both sides
(`*= 1 + 1.25 * size * sd`). Quantiles are lognormal.

## Results

Direct-graded with the exact RULES scorer on the four practice cases
("Average viewers" answers): **practice pinball loss 0.116103**
(p00 0.1057 / p01 0.1359 / p02 0.0769 / p03 0.1459),
LOO-large **0.207953**, LOO-mega **0.209284**. That beats v33
(0.116353 / 0.208240 / 0.209473) on practice by −0.000250, with
LOO-large −0.000287 and LOO-mega −0.000189 — three of the four
practice cases improve individually (p01 worsens slightly, +0.000242).
Median-only practice is **0.178707** (from 0.179039), but median-only
LOO-large is **0.315846** (from 0.314632) — the LOO gain this round is
spread-driven, not a better median on LOO; see `DESIGN.md` for the
honest caveats. Not
submitted: the owner's rule is to submit only a candidate that reaches
the top 10.

## Run it

```sh
MAKEFUN_INPUT_DIR=<dir with *.in> MAKEFUN_OUTPUT_DIR=<dir> python3 main.py
```

Each `*.in` case produces a matching `*.out`:
`{"hours": [{"hour": h, "q": [q05, q25, q50, q75, q95]}]}` for h = 2..24.
