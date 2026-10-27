"""조절효과 분석 기준값을 tests/reference/moderation.json 에 저장한다.

실행: py -3.12 tests/gen_reference_moderation.py
- numpy 의 최소제곱(lstsq)과 공분산 공식으로 독립적으로 계산한다 (이 도구의 Householder QR 과 다른 방법).
  계수, 표준오차, t, p, 신뢰구간(scipy t 분포), R², F, 상호작용 R² 변화량과 F 변화량,
  단순기울기(공분산으로 구한 표준오차), Johnson-Neyman 경계(이차방정식), 그림용 예측값.
"""
import json
import os

import numpy as np
import scipy
from scipy import stats

rng = np.random.default_rng(20261027)


def ols(X, y):
    """절편을 포함한 OLS. X: (n, k) 설계 행렬(절편 제외). 반환: B, cov(B), r2, adjr2, F, pF, dfr, mse, tcoefs..."""
    n, k = X.shape
    A = np.column_stack([np.ones(n), X])
    B, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ B
    sse = float(resid @ resid)
    sst = float(((y - y.mean()) ** 2).sum())
    dfr = n - (k + 1)
    mse = sse / dfr
    cov = mse * np.linalg.inv(A.T @ A)
    r2 = 1 - sse / sst
    adj = 1 - (1 - r2) * (n - 1) / dfr
    F = (r2 / k) / ((1 - r2) / dfr)
    pF = float(stats.f.sf(F, k, dfr))
    return dict(B=B, cov=cov, r2=float(r2), adj=float(adj), F=float(F), pF=pF, dfr=dfr, mse=mse, n=n, k=k)


def design(y, x, w, covs, center, wtype, codes, ncat, ref):
    xm = x.mean() if center else 0.0
    xc = x - xm
    cols_main = [xc]
    inter = []
    info = {"xm": float(x.mean())}
    if wtype == "continuous":
        wm = w.mean() if center else 0.0
        wc = w - wm
        cols_main.append(wc)
        inter.append(xc * wc)
        info.update(wmean=float(w.mean()), wsd=float(w.std(ddof=1)), woff=float(wm))
    else:
        dummies = [c for c in range(ncat) if c != ref]
        for c in dummies:
            d = (codes == c).astype(float)
            cols_main.append(d)
            inter.append(xc * d)
        info.update(dummies=dummies)
    full = np.column_stack(cols_main + inter + list(covs))
    reduced = np.column_stack(cols_main + list(covs))
    return full, reduced, len(cols_main), len(inter), info


def case(y, x, w, covs=(), center=True, wtype="continuous", codes=None, ncat=None, ref=0, conf=95):
    full, reduced, nmain, q, info = design(y, x, w, covs, center, wtype, codes, ncat, ref)
    f = ols(full, y)
    r = ols(reduced, y)
    n = len(y)
    dr2 = f["r2"] - r["r2"]
    Fchg = (dr2 / q) / ((1 - f["r2"]) / f["dfr"])
    pchg = float(stats.f.sf(Fchg, q, f["dfr"]))
    B, cov = f["B"], f["cov"]
    se = np.sqrt(np.diag(cov))
    tcrit = float(stats.t.ppf(1 - (1 - conf / 100) / 2, f["dfr"]))
    t = B / se
    p = 2 * stats.t.sf(np.abs(t), f["dfr"])
    coef = [dict(B=float(B[i]), se=float(se[i]), t=float(t[i]), p=float(p[i]), lo=float(B[i] - tcrit * se[i]), hi=float(B[i] + tcrit * se[i])) for i in range(len(B))]
    # 공분산 행렬에서 절편을 뺀 부분: V[j][l] 은 설계 열 j, l 의 기울기 공분산
    V = cov[1:, 1:]
    slopes = []

    def srow(b, v, label, wval):
        s = float(np.sqrt(v))
        tt = b / s
        return dict(label=label, w=wval, b=float(b), se=s, t=float(tt), p=float(2 * stats.t.sf(abs(tt), f["dfr"])),
                    lo=float(b - tcrit * s), hi=float(b + tcrit * s))

    jn = None
    plot = []
    xmin, xmax = float(x.min()), float(x.max())
    covm = [float(c.mean()) for c in covs]

    def predict(xr, wcoord=None, dummy=None):
        xmod = xr - (x.mean() if center else 0.0)
        val = B[0] + B[1] * xmod
        if wtype == "continuous":
            val += B[2] * wcoord + B[3] * xmod * wcoord
        elif dummy is not None:
            val += B[1 + 1 + dummy] + B[1 + nmain + dummy] * xmod
        for c, m in enumerate(covm):
            val += B[1 + nmain + q + c] * m
        return float(val)

    if wtype == "continuous":
        sd = info["wsd"]
        for lab, k in (("평균 − 1SD", -1), ("평균", 0), ("평균 + 1SD", 1)):
            wm = (0.0 if center else info["wmean"]) + k * sd
            b = B[1] + B[3] * wm
            v = V[0, 0] + 2 * wm * V[0, 2] + wm * wm * V[2, 2]
            slopes.append(srow(b, v, lab, info["wmean"] + k * sd))
            plot.append([[xmin, predict(xmin, wm)], [xmax, predict(xmax, wm)]])
        # Johnson-Neyman (모형 좌표 = 원래 값 − shift)
        shift = info["wmean"] if center else 0.0
        lo, hi = float(w.min()) - shift, float(w.max()) - shift
        b1, b3, v11, v13, v33 = B[1], B[3], V[0, 0], V[0, 2], V[2, 2]
        A_ = b3 ** 2 - tcrit ** 2 * v33
        B_ = 2 * (b1 * b3 - tcrit ** 2 * v13)
        C_ = b1 ** 2 - tcrit ** 2 * v11
        roots = [r.real for r in np.roots([A_, B_, C_]) if abs(r.imag) < 1e-12] if abs(A_) > 1e-12 * max(abs(A_), abs(B_), abs(C_)) else [-C_ / B_]
        roots = sorted(r for r in roots if lo < r < hi)
        cuts = [lo] + roots + [hi]
        regions = []
        for a, b in zip(cuts[:-1], cuts[1:]):
            mid = (a + b) / 2
            sl = b1 + b3 * mid
            vv = v11 + 2 * mid * v13 + mid * mid * v33
            regions.append(dict(lo=float(a + shift), hi=float(b + shift), sig=bool(abs(sl) / np.sqrt(vv) > tcrit), sign=1 if sl >= 0 else -1))
        jn = dict(roots=[float(r + shift) for r in roots], regions=regions)
    else:
        dummies = info["dummies"]
        for c in range(ncat):
            if c == ref:
                slopes.append(srow(B[1], V[0, 0], f"cat{c}", None))
                plot.append([[xmin, predict(xmin)], [xmax, predict(xmax)]])
            else:
                pos = nmain + dummies.index(c)
                slopes.append(srow(B[1] + B[1 + pos], V[0, 0] + 2 * V[0, pos] + V[pos, pos], f"cat{c}", None))
                d = dummies.index(c)
                plot.append([[xmin, predict(xmin, dummy=d)], [xmax, predict(xmax, dummy=d)]])
    return dict(
        input=dict(y=[float(v) for v in y], x=[float(v) for v in x],
                   w=[float(v) for v in w] if wtype == "continuous" else [int(v) for v in codes],
                   covs=[[float(v) for v in c] for c in covs], center=center, wtype=wtype, ncat=ncat, ref=ref, conf=conf),
        expected=dict(n=n, coef=coef, r2=f["r2"], adj=f["adj"], F=f["F"], pF=f["pF"], dfr=f["dfr"], dfm=f["k"],
                      dr2=float(dr2), Fchg=float(Fchg), pchg=pchg, q=q, slopes=slopes, jn=jn, plot=plot,
                      tcrit=tcrit),
    )


cases = {}

n = 80
x = rng.normal(10, 3, n)
w = rng.normal(50, 8, n)
y = 5 + 0.4 * x + 0.1 * w + 0.05 * (x - 10) * (w - 50) + rng.normal(0, 1.5, n)
cases["연속 W, 평균 중심화, 상호작용 있음"] = case(y, x, w)

n = 70
x = rng.normal(0, 1, n) * 4 + 20
w = rng.normal(0, 1, n) * 2 + 5
c1 = rng.normal(0, 1, n)
c2 = rng.integers(20, 60, n).astype(float)
y = 3 + 0.2 * x - 0.3 * w + 0.12 * x * w + 0.5 * c1 + 0.02 * c2 + rng.normal(0, 2, n)
cases["연속 W, 중심화 없음, 공변량 2개"] = case(y, x, w, covs=(c1, c2), center=False, conf=90)

n = 60
x = rng.normal(0, 1, n)
w = rng.normal(0, 1, n)
y = 1 + 0.5 * x + 0.3 * w + rng.normal(0, 1, n)
cases["연속 W, 상호작용 없음"] = case(y, x, w)

n = 90
x = rng.normal(30, 6, n)
codes = rng.integers(0, 3, n)
slope = np.array([0.1, 0.5, 0.9])[codes]
y = 10 + slope * (x - 30) + np.array([0, 2, -1])[codes] + rng.normal(0, 3, n)
cases["범주형 W 3범주 (기준 첫 범주)"] = case(y, x, np.zeros(n), wtype="categorical", codes=codes, ncat=3, ref=0)

n = 64
x = rng.normal(0, 1, n) * 5 + 40
codes = rng.integers(0, 2, n)
c1 = rng.normal(0, 1, n)
y = 20 + np.where(codes == 1, 0.8, -0.2) * (x - 40) + 3 * codes + 1.5 * c1 + rng.normal(0, 4, n)
cases["범주형 W 2범주 (기준 마지막), 공변량"] = case(y, x, np.zeros(n), covs=(c1,), wtype="categorical", codes=codes, ncat=2, ref=1)

n = 120
x = rng.normal(0, 1, n)
w = rng.uniform(-2, 2, n)
y = 0.3 * x + 0.8 * x * w + rng.normal(0, 0.8, n)
cases["Johnson-Neyman 경계가 범위 안 (X의 효과가 W에 따라 부호가 바뀜)"] = case(y, x, w)

n = 50
x = rng.normal(0, 1, n)
w = rng.normal(0, 1, n)
y = 2 * x + 0.05 * x * w + rng.normal(0, 0.3, n)
cases["Johnson-Neyman 경계 없음 (어디서나 유의)"] = case(y, x, w)

path = os.path.join(os.path.dirname(__file__), "reference", "moderation.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "cases": [dict(name=k, **v) for k, v in cases.items()]}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "moderation cases to", path)
