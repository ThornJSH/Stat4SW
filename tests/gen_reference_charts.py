"""상자그림 통계의 기준값을 tests/reference/charts.json 에 저장한다.

실행: py -3.12 tests/gen_reference_charts.py

튜키의 힌지는 R의 fivenum 과 같은 공식(위치 d = 1, n4, (n+1)/2, n+1-n4, n 의 평균)으로 계산한다.
이 도구의 JS 구현은 "중앙값으로 나눈 절반의 중앙값"이라 계산 방식이 다르므로 서로 독립적인 확인이 된다.
울타리(1.5, 3 × 힌지 간격), 수염, 이상치, 극단값은 정의대로 계산한다.
"""
import json
import math
import os

import numpy as np

rng = np.random.default_rng(20261026)


def fivenum(x):
    x = sorted(x)
    n = len(x)
    n4 = math.floor((n + 3) / 2) / 2
    d = [1, n4, (n + 1) / 2, n + 1 - n4, n]
    return [0.5 * (x[math.floor(v) - 1] + x[math.ceil(v) - 1]) for v in d]


def box(values):
    x = sorted(values)
    mn, q1, med, q3, mx = fivenum(x)
    iqr = q3 - q1
    inner = (q1 - 1.5 * iqr, q3 + 1.5 * iqr)
    outer = (q1 - 3 * iqr, q3 + 3 * iqr)
    within = [v for v in x if inner[0] <= v <= inner[1]]
    outliers = [v for v in x if (v < inner[0] or v > inner[1]) and not (v < outer[0] or v > outer[1])]
    extremes = [v for v in x if v < outer[0] or v > outer[1]]
    return {"n": len(x), "q1": q1, "median": med, "q3": q3,
            "lo": min(within) if within else q1, "hi": max(within) if within else q3,
            "outliers": outliers, "extremes": extremes}


sets = {
    "정규 n=20": rng.normal(50, 10, 20),
    "정규 n=21 (홀수)": rng.normal(0, 1, 21),
    "이상치 포함": np.concatenate([rng.normal(10, 1, 30), [16.5, 19.0, -3.0]]),
    "정수 점수 n=40": rng.integers(1, 10, 40).astype(float),
    "n=5": np.array([1.0, 2.0, 4.0, 8.0, 100.0]),
    "n=4": np.array([3.0, 7.0, 7.0, 20.0]),
    "n=1": np.array([5.0]),
    "n=2": np.array([2.0, 9.0]),
}
cases = [{"name": k, "values": [float(v) for v in vals], "expected": box([float(v) for v in vals])} for k, vals in sets.items()]

path = os.path.join(os.path.dirname(__file__), "reference", "charts.json")
os.makedirs(os.path.dirname(path), exist_ok=True)
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"cases": cases}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "box plot cases to", path)
