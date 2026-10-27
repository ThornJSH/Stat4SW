"""scipy로 t검정 기준값을 tests/reference/compare.json 에 저장한다.

실행: py -3.12 tests/gen_reference_compare.py
- ttest_1samp / ttest_ind(equal_var=True/False) / ttest_rel 의 통계량, 자유도, 양측 p값, 평균차이의 신뢰구간
- 피어슨 상관과 p값(pearsonr)
scipy의 confidence_interval: 독립·대응은 평균차이의 구간, 일표본은 평균 자체의 구간(여기서는 검정값을 빼서 평균차이 구간으로 바꿔 기록).
"""
import json
import os
from decimal import Decimal, getcontext
from fractions import Fraction

import numpy as np
import scipy
from scipy import stats

getcontext().prec = 60
rng = np.random.default_rng(20261022)
f = lambda v: float(v)



def exact_one_sample(values, mu, conf):
    """입력 double 값을 분수로 놓고 정확하게 계산한 일표본 t검정 (p, 신뢰구간은 정확한 t로 scipy 분포 함수 사용)."""
    xs = [Fraction(v) for v in values]
    n = len(xs)
    mean = sum(xs) / n
    var = sum((x - mean) ** 2 for x in xs) / (n - 1)
    dv = Decimal(var.numerator) / Decimal(var.denominator)
    se = (dv / n).sqrt()
    diff = Decimal((mean - Fraction(mu)).numerator) / Decimal((mean - Fraction(mu)).denominator)
    t = diff / se
    df = n - 1
    tc = Decimal(float(stats.t.isf((1 - conf / 100) / 2, df)))
    return {"t": float(t), "df": df, "p": float(2 * stats.t.sf(abs(float(t)), df)),
            "diff": float(diff), "lo": float(diff - tc * se), "hi": float(diff + tc * se)}


one = []
for name, x, mu in [
    ("정규 n=20", rng.normal(5, 2, 20), 4.0),
    ("정규 n=200", rng.normal(100, 15, 200), 100.0),
    ("정수 n=40", rng.integers(1, 8, 40).astype(float), 3.5),
    ("큰 평균 n=50", 1e6 + rng.normal(0, 1, 50), 1e6),
    ("n=2", np.array([3.0, 5.0]), 0.0),
]:
    for conf in (0.95, 0.99):
        vals = [f(v) for v in x]
        ex = exact_one_sample(vals, mu, conf * 100)
        # 정확 산술 값을 기준으로 쓰고 scipy와 1e-7 안에서 교차 확인한다.
        # (평균이 1e6처럼 큰 자료에서는 scipy의 t 통계량이 1e-9 수준으로 부정확하다)
        res = stats.ttest_1samp(x, mu)
        ci = res.confidence_interval(conf)
        for label, mine, sc in (("t", ex["t"], f(res.statistic)), ("p", ex["p"], f(res.pvalue)),
                                ("lo", ex["lo"], f(ci.low) - mu), ("hi", ex["hi"], f(ci.high) - mu)):
            assert abs(mine - sc) <= 1e-7 * max(abs(mine), 1e-6), f"{name} {label}: exact {mine!r} vs scipy {sc!r}"
        one.append({"name": name, "conf": conf * 100, "mu": mu, "values": vals, **ex})

ind = []
for name, a, b in [
    ("정규 분산 같음", rng.normal(10, 2, 30), rng.normal(11, 2, 25)),
    ("정규 분산 다름", rng.normal(10, 1, 15), rng.normal(12, 5, 40)),
    ("로그정규", rng.lognormal(0, 0.5, 22), rng.lognormal(0.3, 1.0, 18)),
    ("정수 점수", rng.integers(1, 6, 35).astype(float), rng.integers(2, 7, 28).astype(float)),
    ("크기 2짜리 집단", np.array([1.0, 3.0]), rng.normal(0, 1, 9)),
]:
    for conf in (0.95, 0.9):
        eq = stats.ttest_ind(a, b, equal_var=True)
        we = stats.ttest_ind(a, b, equal_var=False)
        cie, ciw = eq.confidence_interval(conf), we.confidence_interval(conf)
        lev = stats.levene(a, b, center="mean")
        ind.append({
            "name": name, "conf": conf * 100,
            "a": [f(v) for v in a], "b": [f(v) for v in b],
            "diff": f(np.mean(a) - np.mean(b)),
            "equal": {"t": f(eq.statistic), "df": f(eq.df), "p": f(eq.pvalue), "lo": f(cie.low), "hi": f(cie.high)},
            "welch": {"t": f(we.statistic), "df": f(we.df), "p": f(we.pvalue), "lo": f(ciw.low), "hi": f(ciw.high)},
            "levene": {"W": f(lev.statistic), "p": f(lev.pvalue)},
        })

paired = []
for name, x, noise in [
    ("정규 상관 높음", rng.normal(50, 10, 30), 3.0),
    ("정규 상관 낮음", rng.normal(50, 10, 45), 15.0),
    ("정수", rng.integers(1, 10, 25).astype(float), 2.0),
]:
    y = x + rng.normal(1.0, noise, len(x))
    if name == "정수":
        y = np.round(y)
    for conf in (0.95,):
        res = stats.ttest_rel(x, y)
        ci = res.confidence_interval(conf)
        r, p = stats.pearsonr(x, y)
        paired.append({
            "name": name, "conf": conf * 100, "x": [f(v) for v in x], "y": [f(v) for v in y],
            "diff": f(np.mean(x - y)),
            "t": f(res.statistic), "df": int(res.df), "p": f(res.pvalue), "lo": f(ci.low), "hi": f(ci.high),
            "r": f(r), "rp": f(p),
        })

path = os.path.join(os.path.dirname(__file__), "reference", "compare.json")
os.makedirs(os.path.dirname(path), exist_ok=True)
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "one": one, "ind": ind, "paired": paired}, fh, ensure_ascii=False, indent=0)
print("wrote", len(one), "one-sample,", len(ind), "independent,", len(paired), "paired cases to", path)
