"""Forecast GTA VI's launch-day Twitch viewers, hour by hour.

Round v34: starts from the verified v33 cell and searches for a new
combine / median lever that beats it on practice AND both LOO gates.

Members (each {hour: (log-median, log-spread)} for hours 2..24):
  y5    - v4-family blend: OLS view + console/all shape views, tier
          early-access filter, shape regressors [lnV1, g12]
  z0    - v4-family blend: naive early-access filter, [lnV1, g12]
  y1    - v4-family blend: naive early-access filter, [lnV1, lnPREM]
  clone - v9-style blend: whole-library shape pools, sample sd,
          hyphen early-access filter, [lnV1, lnPREM, g12]
  mv    - per-hour median across four views (OLS, console shape,
          all shape, clone); spread = mean of the views' spreads

Two construction facts recovered in v32 and re-used here:
  * the v4-family spread blend is sqrt-variance, i.e. the spread is
    sqrt(sum(w * sig^2) / sum(w)), not a linear weighted mean;
  * the shape views differ per member in pool window, fit hours,
    sd convention, OLS pool order and early-access filter.

Combine (v34 cell, packaged this round — see DESIGN.md): per-band
median source (early = Y1, mid = median(Y5, clone, mv), late =
median(Y5, Z0, Y1)), per-band ramps (early +0.065), size shift,
size-indexed disagreement shift where the disagreement measure is
MIXED per band — the range (max - min) / 2 of the five member
medians in the mid band, the population sd of the five member
medians in the early and late bands (new in v34: the range measure
buys LOO in mid, where one member is often an outlier, while the
population sd is the more robust measure at the band edges);
3-way spread blend now 0.975*(0.45*sigY5 + 0.45*sigZ0 +
0.05*sigClone + 0.05*sigMV), disagreement factor k_dis 0.85, late
spread ramp, asymmetric scales (mid now (1.05, 0.85)), spread cap
0.315, W5 size*disagreement and W1 lower-side factors.

Honest caveat: every combine parameter was selected on the same four
practice launches and the same 37-launch LOO as every round since
v22, so neither gate is a clean holdout; the final launch (lnV1 ~14+)
sits beyond the library max (13.68), so the size-indexed terms
extrapolate hardest exactly there. This is the best verified
fallback, not a guarantee. See DESIGN.md for the full caveat.

Program contract (RULES.txt): reads $MAKEFUN_INPUT_DIR (else /input),
one *.in file per case; writes $MAKEFUN_OUTPUT_DIR (else /output),
one *.out file per case: {"hours": [{"hour": h, "q": [q05, q25, q50,
q75, q95]}]} for hours 2..24. Stdlib only, no network, deterministic.
"""

import json
import math
import os
import sys
from datetime import datetime, timedelta

LEVELS = (0.05, 0.25, 0.5, 0.75, 0.95)
Z_VALUES = (-1.6448536269514722, -0.6744897501960817, 0.0,
            0.6744897501960817, 1.6448536269514722)
HOURS = tuple(range(2, 25))
KINDS = ("y5", "z0", "y1", "clone", "mv")


# ------------------------------------------------------------------ numerics

def _mid(vals):
    s = sorted(vals)
    n = len(s)
    if n == 0:
        return 0.0
    half = n // 2
    return s[half] if n % 2 else 0.5 * (s[half - 1] + s[half])


def _avg(vals):
    return sum(vals) / len(vals) if vals else 0.0


def _sd_pop(vals):
    """Population sd (ddof=0) — the v4 family's convention."""
    if len(vals) < 2:
        return 0.0
    m = _avg(vals)
    return math.sqrt(sum((v - m) ** 2 for v in vals) / len(vals))


def _sd_samp(vals):
    """Sample sd (ddof=1) — the clone's convention."""
    if len(vals) < 2:
        return 0.0
    m = _avg(vals)
    return math.sqrt(sum((v - m) ** 2 for v in vals) / (len(vals) - 1))


def _clamp(x, lo, hi):
    return lo if x < lo else hi if x > hi else x


def _gauss(mat, rhs):
    """Solve a linear system by Gaussian elimination with partial
    pivoting; returns None when the matrix is (near-)singular."""
    n = len(rhs)
    aug = [list(mat[i]) + [rhs[i]] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[piv][col]) < 1e-11:
            return None
        aug[col], aug[piv] = aug[piv], aug[col]
        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col] / aug[col][col]
            if factor == 0.0:
                continue
            for k in range(col, n + 1):
                aug[row][k] -= factor * aug[col][k]
    return [aug[i][n] / aug[i][i] for i in range(n)]


def _fit_ls(design, ys, ridge=1e-9):
    """Least squares through the normal equations; None if singular."""
    p = len(design[0])
    gram = [[sum(r[i] * r[j] for r in design)
             + (ridge if i == j else 0.0)
             for j in range(p)] for i in range(p)]
    vec = [sum(r[i] * y for r, y in zip(design, ys)) for i in range(p)]
    return _gauss(gram, vec)


def _lin(beta, row):
    return sum(b * x for b, x in zip(beta, row))


# ---------------------------------------------------------------------- data

def _when(text):
    return datetime.strptime(text, "%Y-%m-%dT%H:%M:%SZ")


def _hourly_maps(points):
    """({datetime: viewers}, {datetime: channels}) for positive points."""
    viewers, channels = {}, {}
    for pt in points:
        val = pt.get("viewers")
        if val is None or val <= 0:
            continue
        t = _when(pt["from"])
        viewers[t] = float(val)
        ch = pt.get("channels")
        if ch is not None and ch > 0:
            channels[t] = float(ch)
    return viewers, channels


def _roughness(vmap, day_start):
    """Sd of hour-to-hour log-diffs over the 7 days before release."""
    start = _when(day_start)
    pts = sorted((t, v) for t, v in vmap.items()
                 if start - timedelta(days=7) <= t < start)
    diffs = []
    for (t0, v0), (t1, v1) in zip(pts, pts[1:]):
        if (t1 - t0).total_seconds() == 3600:
            diffs.append(math.log(v1) - math.log(v0))
    return _sd_pop(diffs) if len(diffs) >= 10 else None


def _launch_feats(vmap, day_start, hour1):
    """lnV1, lnPREM, g12, pre_rough for a launch, or None."""
    start = _when(day_start)
    prior = [v for t, v in vmap.items()
             if start - timedelta(days=7) <= t < start]
    if not prior:
        return None
    prem = math.log(_avg(prior))
    want = start - timedelta(hours=12)
    best = None
    for t, v in vmap.items():
        gap = abs((t - want).total_seconds())
        if gap <= 3 * 3600 and (best is None or gap < best[0]):
            best = (gap, v)
    g12 = None if best is None else math.log(hour1) - math.log(best[1])
    return {"lnV1": math.log(hour1), "lnPREM": prem, "g12": g12,
            "pre_rough": _roughness(vmap, day_start)}


def _classify(note):
    t = (note or "").lower()
    if "early-access-contaminated" in t or "early-access contaminated" in t:
        return "ea"
    if "tier: weak" in t:
        return "weak"
    return "comparable"


def _read_comparable(comp):
    """Feature dict + release-day curves for one comparable, or None."""
    vmap, _cmap = _hourly_maps(comp["hourly"])
    start = _when(comp["day_start"])
    v1 = vmap.get(start)
    if v1 is None:
        return None
    curve = {}
    for h in range(1, 25):
        t = start + timedelta(hours=h - 1)
        if t in vmap:
            curve[h] = vmap[t]
    if sum(1 for h in curve if h >= 2) < 20:
        return None
    feats = _launch_feats(vmap, comp["day_start"], v1)
    if feats is None:
        return None
    feats["curve"] = curve
    feats["name"] = comp["name"]
    feats["console"] = "PC" not in comp.get("platforms", "")
    note = comp.get("unlock", "") or ""
    feats["tier"] = _classify(note)
    feats["space_ea"] = "early access" in note.lower()
    feats["hyphen_ea"] = "early-access" in note.lower()
    return feats


def _read_target(case):
    """Target features for the case's own launch."""
    vmap, _cmap = _hourly_maps(case["history_hourly"])
    v1 = max(1.0, float(next(k for k in case["known"]
                             if k["hour"] == 1)["viewers"]))
    feats = _launch_feats(vmap, case["day_start"], v1)
    if feats is None:
        feats = {"lnV1": math.log(v1), "lnPREM": math.log(v1),
                 "g12": 0.0, "pre_rough": None}
    feats["curve"] = {}
    return feats


def _load_case(case):
    """(target feats, usable library) for one case dict."""
    tgt = _read_target(case)
    lib = []
    for comp in case["comparables"]:
        if "before release" in comp["name"]:
            continue
        parsed = _read_comparable(comp)
        if parsed is not None:
            lib.append(parsed)
    for key in ("g12", "pre_rough"):
        vals = sorted(f[key] for f in lib if f[key] is not None)
        if not vals:
            continue
        fallback = vals[len(vals) // 2]
        if tgt[key] is None:
            tgt[key] = fallback
        for f in lib:
            if f[key] is None:
                f[key] = fallback
    return tgt, lib


# ---------------------------------------------------------------------- views

def _size_span(lib):
    return min(f["lnV1"] for f in lib), max(f["lnV1"] for f in lib)


def _row_vec(feats, lo, hi, cols):
    """Design row; size clamped to the library range + 0.5 headroom."""
    size = _clamp(feats["lnV1"], lo, hi + 0.5)
    vals = {"lnV1": size, "lnPREM": feats["lnPREM"],
            "g12": feats["g12"]}
    return [1.0] + [vals[c] for c in cols]


def _windowed(lib, lnv1, half=2.5):
    """Launches within +/-half of the target size; whole library when
    fewer than 10 qualify."""
    close = [f for f in lib if abs(f["lnV1"] - lnv1) <= half]
    return close if len(close) >= 10 else list(lib)


def _ols_hour(lib, tgt, hour, window=2.5,
              cols=("lnV1", "lnPREM", "g12"), mult=0.8,
              population=True, count_after=False):
    """Per-hour OLS view -> (log-median, log-spread) or None.

    Pool = window-filtered launches with this hour in their curve;
    centre recentered on the median residual; spread = residual sd x
    mult. population=True uses population sd (v4 family); False uses
    sample sd (clone). count_after=True counts the pool AFTER the
    hour filter (the clone's rule); False counts the window first
    (v4's).
    """
    if count_after:
        pool = [f for f in lib if hour in f["curve"]
                and abs(f["lnV1"] - tgt["lnV1"]) <= window]
        if len(pool) < 10:
            pool = [f for f in lib if hour in f["curve"]]
    else:
        pool = [f for f in _windowed(lib, tgt["lnV1"], window)
                if hour in f["curve"]]
    if len(pool) < 8:
        return None
    lo, hi = _size_span(lib)
    design = [_row_vec(f, lo, hi, cols) for f in pool]
    ys = [math.log(f["curve"][hour]) for f in pool]
    beta = _fit_ls(design, ys)
    if beta is None:
        return None
    resid = [y - _lin(beta, r) for r, y in zip(design, ys)]
    sd = _sd_pop(resid) if population else _sd_samp(resid)
    mu = _lin(beta, _row_vec(tgt, lo, hi, cols)) + _mid(resid)
    return mu, sd * mult


def _ea_pool(base, kind):
    """EA filter on a pool: 'tier' (tier != ea), 'naive' (no SPACE
    'early access' in the raw text), 'hyphen' (no 'early-access')."""
    if kind == "tier":
        return [f for f in base if f["tier"] != "ea"]
    if kind == "naive":
        return [f for f in base if not f["space_ea"]]
    return [f for f in base if not f["hyphen_ea"]]


def _shape_hour(lib, tgt, console_only, feats=("lnV1", "lnPREM", "g12"),
                mult=1.0, population=True, windowed=False,
                ea_kind="hyphen", fit_hours=20, min_fits=5):
    """One-factor shape view -> {hour: (log-median, log-spread)}.

    Median log-ratio template over the pool; per-launch level c +
    amplitude s regressed on `feats`; amplitude clamped [0, 2];
    spread = per-hour residual sd x mult.
    """
    base = _windowed(lib, tgt["lnV1"]) if windowed else list(lib)
    pool = _ea_pool(base, ea_kind)
    if console_only:
        pool = [f for f in pool if f["console"]]
    if len(pool) < 6:
        pool = _ea_pool(base, ea_kind)
    if len(pool) < 6:
        return None
    template = {}
    for h in HOURS:
        vals = [math.log(f["curve"][h] / f["curve"][1]) for f in pool
                if h in f["curve"] and 1 in f["curve"]]
        template[h] = _mid(vals) if vals else 0.0
    fits = []
    for f in pool:
        xs = [template[h] for h in HOURS if h in f["curve"]]
        ys = [math.log(f["curve"][h] / f["curve"][1]) for h in HOURS
              if h in f["curve"]]
        if len(xs) < fit_hours:
            continue
        mx, my = _avg(xs), _avg(ys)
        den = sum((x - mx) ** 2 for x in xs)
        slope = (sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den
                 if den else 0.0)
        fits.append((f, my - slope * mx, slope))
    if len(fits) < min_fits:
        return None
    lo, hi = _size_span(lib)
    design = [_row_vec(f, lo, hi, feats) for f, _, _ in fits]
    beta_c = _fit_ls(design, [c for _, c, _ in fits])
    beta_s = _fit_ls(design, [s for _, _, s in fits])
    if beta_c is None or beta_s is None:
        return None
    xt = _row_vec(tgt, lo, hi, feats)
    c_star = _lin(beta_c, xt)
    s_star = _clamp(_lin(beta_s, xt), 0.0, 2.0)
    out = {}
    for h in HOURS:
        res = [math.log(f["curve"][h] / f["curve"][1])
               - (c + s * template[h])
               for f, c, s in fits if h in f["curve"]]
        sd = _sd_pop(res) if population else _sd_samp(res)
        out[h] = (tgt["lnV1"] + c_star + s_star * template[h], sd * mult)
    return out


# ------------------------------------------------------------------ predictors

def _member_v4(lib, tgt, ea_kind, feats, mix=(0.5, 0.2, 0.3)):
    """Per-hour (log-median, log-spread) for one v4-family member,
    using the sqrt-variance spread blend."""
    con = _shape_hour(lib, tgt, True, feats=feats, mult=1.0,
                      population=True, windowed=True,
                      ea_kind=ea_kind, fit_hours=10, min_fits=6)
    all_ = _shape_hour(lib, tgt, False, feats=feats, mult=1.0,
                       population=True, windowed=True,
                       ea_kind=ea_kind, fit_hours=10, min_fits=6)
    out = {}
    for h in HOURS:
        views = []
        ov = _ols_hour(lib, tgt, h, window=2.5, mult=0.8,
                       population=True)
        if ov is not None:
            views.append((mix[0], ov))
        else:
            views.append((mix[0], (tgt["lnV1"], 0.5)))
        if con is not None and h in con:
            views.append((mix[1], con[h]))
        if all_ is not None and h in all_:
            views.append((mix[2], all_[h]))
        tot = sum(w for w, _ in views)
        mu = sum(w * v[0] for w, v in views) / tot
        sig = math.sqrt(max(0.02, sum(w * v[1] ** 2 for w, v in views)
                            / tot))
        out[h] = (mu, sig)
    return out


def _member_clone(lib, tgt, weights=(0.6, 0.2, 0.2), spread_scale=0.81):
    """Per-hour (log-median, log-spread) for the v9-clone baseline."""
    con = _shape_hour(lib, tgt, True,
                      feats=("lnV1", "lnPREM", "g12"),
                      population=False, ea_kind="hyphen",
                      fit_hours=20, min_fits=5)
    all_ = _shape_hour(lib, tgt, False,
                       feats=("lnV1", "lnPREM", "g12"),
                       population=False, ea_kind="hyphen",
                       fit_hours=20, min_fits=5)
    out = {}
    for h in HOURS:
        views = []
        ov = _ols_hour(lib, tgt, h, window=2.5, population=False,
                       count_after=True)
        if ov is not None:
            views.append((weights[0], ov))
        if con is not None and h in con:
            views.append((weights[1], con[h]))
        if all_ is not None and h in all_:
            views.append((weights[2], all_[h]))
        if not views:
            out[h] = (tgt["lnV1"], 0.5)
            continue
        tot = sum(w for w, _ in views)
        mu = sum(w * v[0] for w, v in views) / tot
        sig = spread_scale * sum(w * v[1] for w, v in views) / tot
        out[h] = (mu, sig)
    return out


def _member_mv(lib, tgt):
    """Per-hour median across four views (the 'views' member):
    per-hour OLS (v4-family), console shape (v4, 3 feats), all shape
    (v4, 3 feats), and the clone forecast. Spread = mean of the four
    views' spreads. Validated to 0.000000 against the cached rows
    from the round that first built this member."""
    con = _shape_hour(lib, tgt, True,
                      feats=("lnV1", "lnPREM", "g12"),
                      population=True, windowed=True,
                      ea_kind="hyphen", fit_hours=10, min_fits=6)
    all_ = _shape_hour(lib, tgt, False,
                       feats=("lnV1", "lnPREM", "g12"),
                       population=True, windowed=True,
                       ea_kind="hyphen", fit_hours=10, min_fits=6)
    cl = _member_clone(lib, tgt)
    out = {}
    for h in HOURS:
        mus, sigs = [], []
        ov = _ols_hour(lib, tgt, h, window=2.5, mult=0.8,
                       population=True)
        if ov is not None:
            mus.append(ov[0])
            sigs.append(ov[1])
        if con is not None and h in con:
            mus.append(con[h][0])
            sigs.append(con[h][1])
        if all_ is not None and h in all_:
            mus.append(all_[h][0])
            sigs.append(all_[h][1])
        mus.append(cl[h][0])
        sigs.append(cl[h][1])
        out[h] = (_mid(mus), _avg(sigs))
    return out


def _member(kind, lib, tgt):
    """{hour: (log-median, log-spread)} for one member predictor."""
    if kind == "y5":
        return _member_v4(lib, tgt, "tier", ("lnV1", "g12"))
    if kind == "z0":
        return _member_v4(lib, tgt, "naive", ("lnV1", "g12"))
    if kind == "y1":
        return _member_v4(lib, tgt, "naive", ("lnV1", "lnPREM"))
    if kind == "clone":
        return _member_clone(lib, tgt)
    if kind == "mv":
        return _member_mv(lib, tgt)
    raise ValueError(kind)


# ------------------------------------------------------------------- combine

def _band_of(hour):
    return "early" if hour <= 6 else "mid" if hour <= 14 else "late"


def _band_pos(hour):
    bounds = {"early": (2, 6), "mid": (7, 14), "late": (15, 24)}
    lo, hi = bounds[_band_of(hour)]
    return (hour - lo) / (hi - lo) if hi > lo else 1.0


def _disagreement(mus, band):
    """Disagreement measure over the five member medians — MIXED per
    band (new in v34). In the mid band the range (max - min) / 2:
    hours where one member is an outlier get a bigger shift, which
    the LOO gate rewards there. In the early and late bands the
    population sd of the five medians: the range is too jumpy at the
    band edges (only a few hours, and the smallest/largest launches
    sit there), so the smoother sd wins on LOO."""
    vals = list(mus.values())
    if band == "mid":
        return (max(vals) - min(vals)) / 2.0
    return _sd_pop(vals)


def _combine(members, ln_v1, hour):
    """(log-median, lower log-spread, upper log-spread) for one hour.

    members: {kind: (mu, sig)} for this hour. v34 cell: shift
    (0.60 + 0.20*size), late base 0.45; size shift -0.04; ramps
    early +0.065 / mid -0.18 / late +0.05; early median Y1, mid
    median median(Y5, clone, mv), late median median(Y5, Z0, Y1);
    k_dis 0.85; spread 0.975*(0.45*sigY5 + 0.45*sigZ0 +
    0.05*sigClone + 0.05*sigMV); late spread ramp 0.975 -> 0.84
    (replaces the 0.975 base factor in the late band); W5 spread *=
    1 + 1.25*size*sd; W1 sig_lo *= 1 + 0.05*size; asym
    (0.90,1.10)/(1.05,0.85)/(0.95,0.85); cap 0.315; disagreement
    sd is MIXED per band — range/2 in mid, population sd in early
    and late (see _disagreement).
    """
    mus = {k: members[k][0] for k in KINDS}
    sigs = {k: members[k][1] for k in KINDS}
    band = _band_of(hour)
    sdv = _disagreement(mus, band)
    size = max(0.0, ln_v1 - 11.5)
    frac = _band_pos(hour)
    if band == "early":
        mu = mus["y1"]
        ramp = 0.065 * frac
        base = 0.60
    elif band == "mid":
        mu = _mid([mus["y5"], mus["clone"], mus["mv"]])
        ramp = -0.18 * frac
        base = 0.60
    else:
        mu = _mid([mus["y5"], mus["z0"], mus["y1"]])
        ramp = 0.05 * frac
        base = 0.45
    mu = mu + ramp + (-0.04) * size
    mu = mu + (base + 0.20 * size) * sdv
    spread = 0.975 * (0.45 * sigs["y5"] + 0.45 * sigs["z0"]
                      + 0.05 * sigs["clone"] + 0.05 * sigs["mv"])
    spread = spread * (1.0 + 0.85 * sdv)
    if band == "late":
        late_factor = 0.975 + (0.84 - 0.975) * frac
        spread = spread / 0.975 * late_factor
    asym = {"early": (0.90, 1.10), "mid": (1.05, 0.85),
            "late": (0.95, 0.85)}[band]
    slo = min(spread * asym[0], 0.315)
    shi = min(spread * asym[1], 0.315)
    slo = min(slo * (1.0 + 0.05 * size), 0.315)
    factor = 1.0 + 1.25 * size * sdv
    slo = min(slo * factor, 0.315)
    shi = min(shi * factor, 0.315)
    return (mu, slo, shi)


# --------------------------------------------------------------------- output

def _to_quantiles(mu_sig):
    """Five quantiles (viewers) from (log-median, lower, upper spread)."""
    mu, slo, shi = mu_sig
    out = []
    for z in Z_VALUES:
        sig = slo if z < 0 else shi if z > 0 else 0.0
        out.append(math.exp(mu + z * sig))
    return out


def _doc_ok(doc):
    """True iff doc is a valid forecast file per RULES.txt."""
    if not isinstance(doc, dict) or set(doc.keys()) != {"hours"}:
        return False
    rows = doc["hours"]
    if not isinstance(rows, list) or len(rows) != 23:
        return False
    seen = set()
    for row in rows:
        if not isinstance(row, dict) or set(row.keys()) != {"hour", "q"}:
            return False
        h = row["hour"]
        if not isinstance(h, int) or h < 2 or h > 24 or h in seen:
            return False
        seen.add(h)
        q = row["q"]
        if not isinstance(q, list) or len(q) != 5:
            return False
        prev = 0.0
        for val in q:
            if (not isinstance(val, (int, float))
                    or not math.isfinite(val)
                    or val <= 0.0 or val > 1e9 or val < prev):
                return False
            prev = val
    return seen == set(range(2, 25))


def _fallback_doc(hour1):
    """House baseline forecast (RULES.txt) — used only if the main
    forecast ever fails the output contract."""
    p = max(1.0, float(hour1))
    rows = []
    for h in HOURS:
        q = [p * math.exp(z * 0.15 * math.sqrt(h - 1)) for z in Z_VALUES]
        rows.append({"hour": h, "q": q})
    return {"hours": rows}


def _forecast_case(case):
    """{"hours": [...]} for one case dict."""
    tgt, lib = _load_case(case)
    members = {k: _member(k, lib, tgt) for k in KINDS}
    ln_v1 = tgt["lnV1"]
    rows = []
    for h in HOURS:
        cell = {k: members[k][h] for k in KINDS}
        q = _to_quantiles(_combine(cell, ln_v1, h))
        q = [min(max(v, 1.0), 1e9) for v in q]
        for i in range(1, 5):
            if q[i] < q[i - 1]:
                q[i] = q[i - 1]
        rows.append({"hour": h, "q": q})
    doc = {"hours": rows}
    if _doc_ok(doc):
        return doc
    hour1 = max(1.0, math.exp(ln_v1))
    return _fallback_doc(hour1)


def main():
    in_dir = os.environ.get("MAKEFUN_INPUT_DIR", "/input")
    out_dir = os.environ.get("MAKEFUN_OUTPUT_DIR", "/output")
    os.makedirs(out_dir, exist_ok=True)
    for name in sorted(os.listdir(in_dir)):
        if not name.endswith(".in"):
            continue
        with open(os.path.join(in_dir, name), "r", encoding="utf-8") as fh:
            case = json.load(fh)
        doc = _forecast_case(case)
        out_name = name[:-3] + ".out"
        with open(os.path.join(out_dir, out_name), "w",
                  encoding="utf-8") as fh:
            json.dump(doc, fh)
    return 0


if __name__ == "__main__":
    sys.exit(main())
