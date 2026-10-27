"""매개효과 분석 기준값을 tests/reference/mediation.json 에 저장한다.

실행: py -3.12 tests/gen_reference_mediation.py
- 경로 계수(a, b, c', c)와 표준오차, t, p, 신뢰구간, R², F: numpy 최소제곱(lstsq)과 공분산 공식으로 독립 계산.
- Sobel: SE = sqrt(a² SE(b)² + b² SE(a)²), Z = ab/SE, 양측 p (scipy 정규분포).
- 부트스트랩: 이 도구와 같은 난수(mulberry32)와 같은 추출 방식(표본마다 n번 floor(rnd·n))으로 표본을 뽑고,
  표본마다 numpy lstsq 로 a_i, b_i 를 새로 구한다. 퍼센타일은 numpy 선형 보간(np.percentile 기본), 부트스트랩 SE 는 표본표준편차(ddof=1).
"""
import json
import os

import numpy as np
import scipy
from scipy import stats

rng = np.random.default_rng(20261029)
MASK = 0xFFFFFFFF


def imul(a, b):
    return (a * b) & MASK


def mulberry32(seed):
    state = {"a": seed & MASK}

    def rnd():
        a = (state["a"] + 0x6D2B79F5) & MASK
        state["a"] = a
        t = imul(a ^ (a >> 15), 1 | a)
        t = ((t + imul(t ^ (t >> 7), 61 | t)) & MASK) ^ t
        return ((t ^ (t >> 14)) & MASK) / 4294967296

    return rnd


def ols(X, y):
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
    F = (r2 / k) / ((1 - r2) / dfr)
    return dict(B=B, se=np.sqrt(np.diag(cov)), r2=float(r2), F=float(F), pF=float(stats.f.sf(F, k, dfr)), dfr=dfr, k=k)


def coef_dict(f, i, conf):
    tcrit = stats.t.ppf(1 - (1 - conf / 100) / 2, f["dfr"])
    b, se = f["B"][i], f["se"][i]
    t = b / se
    return dict(B=float(b), se=float(se), t=float(t), p=float(2 * stats.t.sf(abs(t), f["dfr"])), lo=float(b - tcrit * se), hi=float(b + tcrit * se))


def case(y, x, ms, covs=(), conf=95, B=2000, seed=1):
    n = len(y)
    m = len(ms)
    cov_cols = list(covs)
    Xa = np.column_stack([x] + cov_cols)
    amods = [ols(Xa, mm) for mm in ms]
    ymod = ols(np.column_stack([x] + list(ms) + cov_cols), y)
    tmod = ols(Xa, y)
    a = [coef_dict(f, 1, conf) for f in amods]
    b = [coef_dict(ymod, 2 + i, conf) for i in range(m)]
    cd, ct = coef_dict(ymod, 1, conf), coef_dict(tmod, 1, conf)
    sdx, sdy = float(x.std(ddof=1)), float(y.std(ddof=1))
    effects, sobel = [], []
    for i in range(m):
        ab = a[i]["B"] * b[i]["B"]
        se = float(np.sqrt(a[i]["B"] ** 2 * b[i]["se"] ** 2 + b[i]["B"] ** 2 * a[i]["se"] ** 2))
        z = ab / se
        effects.append(dict(ab=float(ab), std=float(ab * sdx / sdy)))
        sobel.append(dict(ab=float(ab), se=se, z=float(z), p=float(2 * stats.norm.sf(abs(z)))))
    total = float(sum(e["ab"] for e in effects))

    rnd = mulberry32(seed)
    ind = [[] for _ in range(m)]
    tot = []
    for _ in range(B):
        idx = np.array([int(rnd() * n) for _ in range(n)])
        xs, ys = x[idx], y[idx]
        mss = [mm[idx] for mm in ms]
        cvs = [c[idx] for c in cov_cols]
        Xa_b = np.column_stack([xs] + cvs)
        ab_sum = 0.0
        a_b = [ols_coefs(Xa_b, mm) [1] for mm in mss]
        y_b = ols_coefs(np.column_stack([xs] + mss + cvs), ys)
        for i in range(m):
            v = a_b[i] * y_b[2 + i]
            ind[i].append(v)
            ab_sum += v
        tot.append(ab_sum)
    alpha = 1 - conf / 100

    def ci(vals):
        arr = np.array(vals)
        return dict(lo=float(np.percentile(arr, 100 * alpha / 2)), hi=float(np.percentile(arr, 100 * (1 - alpha / 2))), se=float(arr.std(ddof=1)))

    return dict(
        input=dict(y=[float(v) for v in y], x=[float(v) for v in x], ms=[[float(v) for v in mm] for mm in ms],
                   covs=[[float(v) for v in c] for c in cov_cols], conf=conf, boot=B, seed=seed),
        expected=dict(n=n, a=a, b=b, cDirect=cd, cTotal=ct, effects=effects, totalIndirect=total, totalStd=float(total * sdx / sdy),
                      sobel=sobel, bootInd=[ci(v) for v in ind], bootTotal=ci(tot),
                      models=dict(a=[dict(r2=f["r2"], F=f["F"], pF=f["pF"]) for f in amods], y=dict(r2=ymod["r2"], F=ymod["F"], pF=ymod["pF"]), t=dict(r2=tmod["r2"], F=tmod["F"], pF=tmod["pF"])),
                      identity=float(cd["B"] + total - ct["B"])),
    )


def ols_coefs(X, y):
    A = np.column_stack([np.ones(X.shape[0]), X])
    B, *_ = np.linalg.lstsq(A, y, rcond=None)
    return B


cases = {}
n = 80
x = rng.normal(0, 1, n)
m1 = 0.6 * x + rng.normal(0, 1, n)
y = 0.3 * x + 0.5 * m1 + rng.normal(0, 1, n)
cases["매개변수 1개, 부분 매개"] = case(y, x, [m1], B=3000, seed=777)

n = 100
x = rng.normal(10, 2, n)
c1 = rng.normal(0, 1, n)
m1 = 0.5 * x + 0.3 * c1 + rng.normal(0, 2, n)
m2 = 0.2 * x + rng.normal(0, 2, n)
y = 0.1 * x + 0.4 * m1 + 0.05 * m2 + 0.5 * c1 + rng.normal(0, 2, n)
cases["매개변수 2개 + 공변량 1개"] = case(y, x, [m1, m2], covs=[c1], B=3000, seed=2026)

n = 60
x = rng.normal(50, 10, n)
c1 = rng.normal(0, 1, n)
c2 = rng.integers(20, 50, n).astype(float)
m1 = 0.4 * x + rng.normal(0, 6, n)
m2 = -0.2 * x + rng.normal(0, 6, n)
m3 = 0.1 * x + rng.normal(0, 6, n)
y = 0.05 * x + 0.3 * m1 - 0.2 * m2 + 0.1 * m3 + 0.3 * c1 + 0.02 * c2 + rng.normal(0, 5, n)
cases["매개변수 3개 + 공변량 2개, 90% 신뢰구간"] = case(y, x, [m1, m2, m3], covs=[c1, c2], conf=90, B=2000, seed=5)

n = 90
x = rng.normal(0, 1, n)
m1 = 0.7 * x + rng.normal(0, 1, n)
y = 0.5 * x - 0.6 * m1 + rng.normal(0, 1, n)
cases["억제 효과 (직접효과와 간접효과의 부호가 다름)"] = case(y, x, [m1], B=2000, seed=99)

n = 70
x = rng.normal(0, 1, n)
m1 = rng.normal(0, 1, n)
y = 0.4 * x + rng.normal(0, 1, n)
cases["매개 없음 (a, b 모두 0)"] = case(y, x, [m1], B=2000, seed=31)

path = os.path.join(os.path.dirname(__file__), "reference", "mediation.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "cases": [dict(name=k, **v) for k, v in cases.items()]}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "mediation cases to", path)
