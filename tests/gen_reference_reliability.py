"""신뢰도 분석(크론바흐 α) 기준값을 tests/reference/reliability.json 에 저장한다.

실행: py -3.12 tests/gen_reference_reliability.py
- 입력 double 값을 분수(Fraction)로 놓고 정확 산술로 계산한다 (공분산 행렬, α, 문항 제외 시 α, 수정된 문항-전체 상관).
  제곱근이 필요한 값(표준편차, 상관)만 Decimal(60자리)로 계산해 float 으로 저장한다.
- 목록별 결측 제거: 모든 문항이 유효한 사례만 쓴다.
- Feldt 신뢰구간: (1-α̂)/(1-α) ~ F(n-1, (n-1)(k-1)) 이므로 α ∈ [1-(1-α̂)/F(γ/2), 1-(1-α̂)/F(1-γ/2)]. F 분위수는 scipy.
- scipy 에는 크론바흐 α 가 없어서, 같은 값을 numpy 공분산 공식으로도 계산해 교차 확인한다.
"""
import json
import os
from decimal import Decimal, getcontext
from fractions import Fraction

import numpy as np
import scipy
from scipy import stats

getcontext().prec = 60
rng = np.random.default_rng(20261026)


def dsqrt(fr):
    return (Decimal(fr.numerator) / Decimal(fr.denominator)).sqrt()


def exact(columns, conf=95):
    """columns: 문항별 값 목록 (결측은 None). 목록별 결측 제거 후 정확 산술로 계산."""
    k = len(columns)
    n_all = len(columns[0])
    keep = [i for i in range(n_all) if all(c[i] is not None for c in columns)]
    n = len(keep)
    xs = [[Fraction(c[i]) for i in keep] for c in columns]
    means = [sum(x) / n for x in xs]
    C = [[sum((xs[i][r] - means[i]) * (xs[j][r] - means[j]) for r in range(n)) / (n - 1) for j in range(k)] for i in range(k)]
    tr = sum(C[i][i] for i in range(k))
    S = sum(C[i][j] for i in range(k) for j in range(k))
    out = {"k": k, "n": n, "nAll": n_all, "scaleMean": float(sum(means)), "scaleVar": float(S)}
    alpha = None
    if S > 0:
        alpha = Fraction(k, k - 1) * (1 - tr / S)
        out["alpha"] = float(alpha)
    else:
        out["alpha"] = None

    # 표준화 α: 문항 간 상관의 평균
    corr_ok = all(C[i][i] > 0 for i in range(k))
    if corr_ok and k >= 2:
        rs = []
        for i in range(k):
            for j in range(i + 1, k):
                rs.append(float(Decimal(C[i][j].numerator) / Decimal(C[i][j].denominator) / (dsqrt(C[i][i]) * dsqrt(C[j][j]))))
        rbar = sum(rs) / len(rs)
        out["interitem"] = {"mean": rbar, "min": min(rs), "max": max(rs)}
        out["alphaStd"] = k * rbar / (1 + (k - 1) * rbar)
        # 문항 간 상관 행렬 (상삼각)
        out["r"] = [[(1.0 if i == j else float(Decimal(C[i][j].numerator) / Decimal(C[i][j].denominator) / (dsqrt(C[i][i]) * dsqrt(C[j][j])))) for j in range(k)] for i in range(k)]
    else:
        out["interitem"] = None
        out["alphaStd"] = None
        out["r"] = None

    items = []
    for i in range(k):
        rowsum = sum(C[i][j] for j in range(k))
        var_rest = S - 2 * rowsum + C[i][i]
        cov_rest = rowsum - C[i][i]
        it = {"mean": float(means[i]), "sd": float(dsqrt(C[i][i])), "n": n,
              "scaleMeanIfDeleted": float(sum(means) - means[i]), "scaleVarIfDeleted": float(var_rest)}
        if C[i][i] > 0 and var_rest > 0:
            it["corrected"] = float(Decimal(cov_rest.numerator) / Decimal(cov_rest.denominator) / (dsqrt(C[i][i]) * dsqrt(var_rest)))
        else:
            it["corrected"] = None
        if k >= 3 and var_rest > 0:
            it["alphaIfDeleted"] = float(Fraction(k - 1, k - 2) * (1 - (tr - C[i][i]) / var_rest))
        else:
            it["alphaIfDeleted"] = None
        items.append(it)
    out["items"] = items

    # Feldt 신뢰구간
    if alpha is not None and n >= 3 and k >= 2:
        g = 1 - conf / 100
        df1, df2 = n - 1, (n - 1) * (k - 1)
        a = float(alpha)
        lo = 1 - (1 - a) / stats.f.ppf(g / 2, df1, df2)
        hi = 1 - (1 - a) / stats.f.ppf(1 - g / 2, df1, df2)
        out["ci"] = {"conf": conf, "lo": float(lo), "hi": float(hi)}
    else:
        out["ci"] = None

    # numpy 로 교차 확인
    arr = np.array([[float(c[i]) for c in columns] for i in keep])
    cv = np.cov(arr, rowvar=False, ddof=1)
    if out["alpha"] is not None:
        chk = k / (k - 1) * (1 - np.trace(cv) / cv.sum())
        assert abs(chk - out["alpha"]) <= 1e-9 * max(1, abs(chk)), (chk, out["alpha"])
    return out


def likert(n, k, load=0.8, noise=0.8, reverse=()):
    trait = rng.normal(0, 1, n)
    cols = []
    for j in range(k):
        v = np.clip(np.round(3 + load * trait + noise * rng.normal(0, 1, n)), 1, 5)
        if j in reverse:
            v = 6 - v
        cols.append([float(x) for x in v])
    return cols


def with_missing(cols, rate):
    out = []
    for c in cols:
        out.append([None if rng.random() < rate else v for v in c])
    return out


cases = {}
cases["리커트 5문항 n=60"] = likert(60, 5)
cases["리커트 5문항 결측 포함 n=80"] = with_missing(likert(80, 5), 0.06)
cases["역문항 포함 (문항3·5)"] = likert(60, 5, reverse=(2, 4))
cases["2문항"] = likert(50, 2)
cases["3문항 (문항 제외 시 2문항)"] = likert(40, 3, load=0.6, noise=1.0)
cases["작은 표본 4문항 n=6"] = likert(6, 4)
cases["이분 문항 6개 (KR-20)"] = [[float(x) for x in (rng.random(40) < p).astype(int)] for p in (0.3, 0.5, 0.6, 0.4, 0.7, 0.5)]
cases["큰 평균 4문항"] = [[1e6 + v for v in c] for c in likert(30, 4)]
const = likert(40, 4)
const[1] = [3.0] * 40
cases["분산이 0인 문항 포함"] = const
# 연속 척도: 같은 잠재 요인에서 3문항
t = rng.normal(0, 1, 70)
cases["연속 척도 3문항"] = [[float(round(x, 2)) for x in (50 + 10 * t + 6 * rng.normal(0, 1, 70))] for _ in range(3)]

ref = []
for name, cols in cases.items():
    ref.append({"name": name, "columns": cols, "expected": exact(cols)})

path = os.path.join(os.path.dirname(__file__), "reference", "reliability.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "cases": ref}, fh, ensure_ascii=False, indent=0)
print("wrote", len(ref), "reliability cases to", path)
