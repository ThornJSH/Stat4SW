"""scipy로 정규성(Jarque-Bera, Shapiro-Wilk), 등분산(Levene) 기준값을 tests/reference/assumptions.json 에 저장한다.

실행: py -3.12 tests/gen_reference_assumptions.py
- jarque_bera: 정확 산술(분수)로 계산한 값을 기준으로 쓰고 scipy와 교차 확인한다.
  scipy는 모집단식 왜도·첨도(g1, g2)와 χ²(2) 꼬리확률을 쓰며, 이 도구와 같은 정의다.
- levene: center='mean'이 원래의 Levene, center='median'이 Brown-Forsythe.
"""
import json
import os
from decimal import Decimal, getcontext
from fractions import Fraction

import numpy as np
import scipy
from scipy import stats

getcontext().prec = 60
rng = np.random.default_rng(20261021)

jb_sets = {
    "정규 n=30": rng.normal(0, 1, 30),
    "정규 n=200": rng.normal(10, 3, 200),
    "로그정규 n=200": rng.lognormal(0, 1, 200),
    "균등 n=500": rng.uniform(0, 1, 500),
    "정수 1~5 n=60": rng.integers(1, 6, 60).astype(float),
    "큰 평균 n=100": 1e6 + rng.normal(0, 1, 100),
    "n=3": np.array([1.0, 2.0, 7.0]),
}


def exact_jb(values):
    """입력 double 값을 분수로 놓고 정확하게 계산한 g1, g2, JB, p (p = exp(-JB/2), χ² 자유도 2의 꼬리확률)."""
    xs = [Fraction(v) for v in values]
    n = len(xs)
    mean = sum(xs) / n
    d = [x - mean for x in xs]
    m2 = sum(v ** 2 for v in d) / n
    m3 = sum(v ** 3 for v in d) / n
    m4 = sum(v ** 4 for v in d) / n
    dm2 = Decimal(m2.numerator) / Decimal(m2.denominator)
    g1 = (Decimal(m3.numerator) / Decimal(m3.denominator)) / dm2 ** Decimal("1.5")
    g2 = (Decimal(m4.numerator) / Decimal(m4.denominator)) / dm2 ** 2 - 3
    jb = Decimal(n) / 6 * (g1 ** 2 + g2 ** 2 / 4)
    return float(g1), float(g2), float(jb), float((-jb / 2).exp())


jb = []
for name, x in jb_sets.items():
    vals = [float(v) for v in x]
    g1, g2, jbv, p = exact_jb(vals)
    # 정확 산술 값을 기준으로 쓰되 scipy와 1e-7 안에서 일치하는지 교차 확인한다.
    # (평균이 1e6처럼 큰 자료에서는 scipy 자체가 1e-9 수준으로 부정확하다)
    res = stats.jarque_bera(x)
    for label, mine, sc in (("jb", jbv, float(res.statistic)), ("p", p, float(res.pvalue))):
        assert abs(mine - sc) <= 1e-7 * max(abs(mine), 1e-300), f"{name} {label}: exact {mine!r} vs scipy {sc!r}"
    jb.append({"name": name, "values": vals, "jb": jbv, "p": p, "g1": g1, "g2": g2})

lev_sets = {
    "정규 3집단 (분산 다름)": [rng.normal(0, 1, 20), rng.normal(0, 2, 25), rng.normal(0, 4, 15)],
    "정규 2집단 (분산 같음)": [rng.normal(5, 2, 40), rng.normal(6, 2, 35)],
    "로그정규 4집단": [rng.lognormal(0, s, n) for s, n in [(0.5, 12), (1.0, 18), (0.8, 9), (1.5, 30)]],
    "정수 3집단 (동점 많음)": [rng.integers(1, 6, n).astype(float) for n in (15, 22, 18)],
    "큰 평균 2집단": [1e6 + rng.normal(0, 1, 30), 1e6 + rng.normal(0, 3, 30)],
    "크기 2짜리 집단 포함": [np.array([1.0, 4.0]), rng.normal(0, 1, 10), rng.normal(0, 2, 7)],
}
lev = []
for name, groups in lev_sets.items():
    for center in ("mean", "median"):
        res = stats.levene(*groups, center=center)
        lev.append({
            "name": name, "center": center,
            "groups": [[float(v) for v in g] for g in groups],
            "W": float(res.statistic), "p": float(res.pvalue),
        })

# Shapiro-Wilk: scipy(Royston AS R94)가 기준. 기존 기준값이 바뀌지 않도록 별도 난수열을 쓴다.
rng_sw = np.random.default_rng(20261022)
sw_sets = {
    "n=3": np.array([1.0, 2.0, 4.0]),
    "n=3 (거의 균등)": np.array([1.0, 2.0, 3.1]),
    "n=4": np.array([2.0, 3.0, 5.0, 9.0]),
    "n=5": np.array([1.2, 1.9, 2.4, 2.6, 7.0]),
    "정규 n=6": rng_sw.normal(0, 1, 6),
    "정규 n=11": rng_sw.normal(0, 1, 11),
    "정규 n=12": rng_sw.normal(0, 1, 12),
    "정규 n=30": rng_sw.normal(50, 10, 30),
    "정규 n=100": rng_sw.normal(0, 1, 100),
    "정규 n=1000": rng_sw.normal(5, 2, 1000),
    "로그정규 n=40": rng_sw.lognormal(0, 1, 40),
    "로그정규 n=500": rng_sw.lognormal(0, 1, 500),
    "균등 n=50": rng_sw.uniform(0, 1, 50),
    "정수 1~5 n=60 (동점 많음)": rng_sw.integers(1, 6, 60).astype(float),
    "이봉 n=80": np.concatenate([rng_sw.normal(-3, 1, 40), rng_sw.normal(3, 1, 40)]),
    "큰 평균 n=40": 1e6 + rng_sw.normal(0, 1, 40),
    "정규 n=5000": rng_sw.normal(0, 1, 5000),
}
sw = []
for name, x in sw_sets.items():
    res = stats.shapiro(x)
    sw.append({"name": name, "values": [float(v) for v in x], "W": float(res.statistic), "p": float(res.pvalue)})

path = os.path.join(os.path.dirname(__file__), "reference", "assumptions.json")
os.makedirs(os.path.dirname(path), exist_ok=True)
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "jb": jb, "levene": lev, "sw": sw}, fh, ensure_ascii=False, indent=0)
print("wrote", len(jb), "JB sets,", len(lev), "Levene cases and", len(sw), "Shapiro-Wilk sets to", path)
