"""효과크기(Cohen d, Hedges g)와 비중심 t 분포의 기준값을 tests/reference/effect.json 에 저장한다.

실행: py -3.12 tests/gen_reference_effect.py
- 비중심 t 누적확률: scipy.stats.nct.cdf 를 여러 (t, df, ncp) 조합에서 직접 저장
- 점추정: numpy 로 계산 (독립: 합동 표준편차, 대응: 차이의 표준편차 d_z 와 두 표준편차 평균 d_av, 일표본: 표본 표준편차)
- Hedges 보정 J = 1 − 3/(4·df − 1)
- 정확 신뢰구간: t = d·√n_eff 를 비중심 t 의 관측값으로 보고 scipy.optimize.brentq 로
  nct.cdf(t, df, δ_L) = 1 − α/2, nct.cdf(t, df, δ_U) = α/2 를 풀어 d = δ/√n_eff  (Hedges g 의 구간은 d 의 구간 × J)
  n_eff: 독립 n1·n2/(n1+n2), 대응·일표본 n
"""
import json
import os

import numpy as np
import scipy
from scipy import stats
from scipy.optimize import brentq

rng = np.random.default_rng(20261036)


def ci_ncp(t, df, conf):
    a = (1 - conf / 100) / 2

    def solve(p):
        f = lambda dl: stats.nct.cdf(t, df, dl) - p
        lo, hi = t - 8.0, t + 8.0
        while f(lo) < 0:
            lo -= 8
        while f(hi) > 0:
            hi += 8
        return brentq(f, lo, hi, xtol=1e-13, rtol=1e-14)

    return solve(1 - a), solve(a)          # (δ_L, δ_U)


def effect(d, n_eff, df, conf=95):
    t = d * np.sqrt(n_eff)
    lo, hi = ci_ncp(t, df, conf)
    J = 1 - 3 / (4 * df - 1)
    return dict(d=float(d), d_lo=float(lo / np.sqrt(n_eff)), d_hi=float(hi / np.sqrt(n_eff)),
                J=float(J), g=float(d * J), g_lo=float(lo / np.sqrt(n_eff) * J), g_hi=float(hi / np.sqrt(n_eff) * J), t=float(t), df=int(df), n_eff=float(n_eff))


# ---- 비중심 t 누적확률 격자 ----
grid = []
for df in [1, 2, 3, 5, 10, 30, 100, 1000, 50000]:
    for t in [-6.0, -2.5, -0.7, 0.0, 0.4, 1.9, 4.0, 9.0]:
        for ncp in [-4.0, -1.2, 0.0, 0.8, 2.5, 6.0]:
            p = float(stats.nct.cdf(t, df, ncp))
            if np.isfinite(p):                      # scipy 가 NaN 을 돌려주는 조합은 기준값에서 뺀다
                grid.append(dict(t=t, df=df, ncp=ncp, p=p))
# 효과크기에서 실제로 쓰는 모양 (큰 ncp)
for t, df, ncp in [(12.0, 58, 11.0), (-9.0, 40, -8.2), (25.0, 998, 24.0), (3.3, 9, 2.1), (0.2, 4, 0.0), (30.0, 100, 28.0)]:
    grid.append(dict(t=t, df=df, ncp=ncp, p=float(stats.nct.cdf(t, df, ncp))))

cases = []

# 독립표본
for name, n1, n2, m1, m2, s1, s2, conf in [("보통 효과", 24, 30, 52.0, 47.0, 9.0, 11.0, 95), ("큰 효과, 90% 구간", 15, 15, 70.0, 55.0, 8.0, 8.5, 90),
                                          ("작은 표본, 음수 효과", 6, 9, 20.0, 26.0, 5.0, 6.0, 95), ("큰 표본, 작은 효과", 300, 280, 100.5, 99.0, 15.0, 14.0, 95)]:
    a = np.round(rng.normal(m1, s1, n1), 2)
    b = np.round(rng.normal(m2, s2, n2), 2)
    df = n1 + n2 - 2
    sp = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1)) / df)
    d = (a.mean() - b.mean()) / sp
    cases.append(dict(kind="independent", name=name, conf=conf, a=a.tolist(), b=b.tolist(), sp=float(sp), **effect(d, n1 * n2 / (n1 + n2), df, conf)))

# 대응표본
for name, n, mu, sd, rho, shift, conf in [("사전·사후 향상", 28, 60.0, 10.0, 0.6, 5.0, 95), ("거의 변화 없음", 40, 3.0, 0.8, 0.5, 0.05, 95), ("작은 표본", 8, 50.0, 7.0, 0.7, -4.0, 99)]:
    x = rng.normal(mu, sd, n)
    y = rho * (x - mu) + np.sqrt(1 - rho ** 2) * rng.normal(0, sd, n) + mu + shift
    x, y = np.round(x, 2), np.round(y, 2)
    diff = x - y
    dz = diff.mean() / diff.std(ddof=1)
    dav = diff.mean() / ((x.std(ddof=1) + y.std(ddof=1)) / 2)
    cases.append(dict(kind="paired", name=name, conf=conf, x=x.tolist(), y=y.tolist(), d_av=float(dav), sd_diff=float(diff.std(ddof=1)), **effect(dz, n, n - 1, conf)))

# 일표본
for name, n, mu0, m, s, conf in [("기준값보다 높음", 35, 70.0, 74.5, 9.0, 95), ("기준값보다 낮음", 12, 3.5, 3.0, 0.6, 95)]:
    v = np.round(rng.normal(m, s, n), 2)
    d = (v.mean() - mu0) / v.std(ddof=1)
    cases.append(dict(kind="one", name=name, conf=conf, values=v.tolist(), mu=mu0, sd=float(v.std(ddof=1)), **effect(d, n, n - 1, conf)))

path = os.path.join(os.path.dirname(__file__), "reference", "effect.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "grid": grid, "cases": cases}, fh, ensure_ascii=False, indent=0)
print("wrote", len(grid), "nct points and", len(cases), "effect cases;", [(c["kind"], round(c["d"], 3), round(c["d_lo"], 3), round(c["d_hi"], 3)) for c in cases])
