"""기술통계 기준값을 tests/reference/descriptive.json 에 저장한다.

실행: py -3.12 tests/gen_reference_descriptive.py

기준값은 입력 double 값을 그대로 분수(Fraction)로 놓고 정확하게 계산한 값이다.
scipy(bias=False 는 G1, G2 공식)와 1e-7 안에서 일치하는지 교차 확인하고, 어긋나면 중단한다.
(큰 평균 사례에서 scipy 자체가 1e-8 수준으로 부정확해서 scipy 값을 그대로 기준으로 쓰지 않는다.)
왜도·첨도의 표준오차는 JS 테스트에서 공식과 알려진 표 값으로 따로 확인한다.
"""
import json
import os
from decimal import Decimal, getcontext
from fractions import Fraction

import numpy as np
import scipy
from scipy import stats

getcontext().prec = 60
rng = np.random.default_rng(12345)

datasets = {
    "정규 n=5": rng.normal(50, 10, 5),
    "정규 n=30": rng.normal(0, 1, 30),
    "로그정규 n=200": rng.lognormal(0, 1, 200),
    "큰 평균 n=100": 1e6 + rng.normal(0, 1, 100),
    "정수 1~5 n=60": rng.integers(1, 6, 60).astype(float),
    "작은 값 n=1000": rng.normal(0, 1e-6, 1000),
    "n=4": np.array([2.0, 4.0, 4.0, 10.0]),
}


def dec(fr):
    return Decimal(fr.numerator) / Decimal(fr.denominator)


def exact(values):
    xs = [Fraction(v) for v in values]
    n = len(xs)
    total = sum(xs)
    mean = total / n
    d = [x - mean for x in xs]
    var = sum(v * v for v in d) / (n - 1)
    s = dec(var).sqrt()
    out = {
        "n": n,
        "min": float(min(xs)),
        "max": float(max(xs)),
        "range": float(max(xs) - min(xs)),
        "sum": float(total),
        "mean": float(mean),
        "variance": float(var),
        "sd": float(s),
        "seMean": float(s / Decimal(n).sqrt()),
        "skewness": None,
        "kurtosis": None,
    }
    if n > 2:
        d3 = dec(sum(v ** 3 for v in d))
        out["skewness"] = float(Decimal(n) / Decimal((n - 1) * (n - 2)) * d3 / s ** 3)
    if n > 3:
        d4 = dec(sum(v ** 4 for v in d))
        out["kurtosis"] = float(
            Decimal(n * (n + 1)) / Decimal((n - 1) * (n - 2) * (n - 3)) * d4 / s ** 4
            - Decimal(3 * (n - 1) ** 2) / Decimal((n - 2) * (n - 3))
        )
    return out


def check_against_scipy(name, ex, x):
    pairs = {
        "mean": np.mean(x),
        "variance": np.var(x, ddof=1),
        "sd": np.std(x, ddof=1),
        "seMean": stats.sem(x, ddof=1),
    }
    if ex["skewness"] is not None:
        pairs["skewness"] = stats.skew(x, bias=False)
    if ex["kurtosis"] is not None:
        pairs["kurtosis"] = stats.kurtosis(x, fisher=True, bias=False)
    for k, sc in pairs.items():
        rel = abs(float(sc) - ex[k]) / max(abs(ex[k]), 1e-300)
        assert rel < 1e-7, f"{name} {k}: scipy {float(sc)!r} vs exact {ex[k]!r} (rel {rel:.2e})"


out = []
for name, x in datasets.items():
    ex = exact([float(v) for v in x])
    check_against_scipy(name, ex, x)
    out.append({"name": name, "values": [float(v) for v in x], **ex})

path = os.path.join(os.path.dirname(__file__), "reference", "descriptive.json")
os.makedirs(os.path.dirname(path), exist_ok=True)
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"reference": "exact rational arithmetic, cross-checked with scipy " + scipy.__version__, "cases": out},
              fh, ensure_ascii=False, indent=0)
print("wrote", len(out), "datasets to", path)
