"""OLS 회귀의 기준값을 tests/reference/regression.json 에 저장한다.

실행: py -3.12 tests/gen_reference_regression.py

기준값은 입력 double 값을 분수(Fraction)로 놓고 정규방정식을 정확하게 풀어 만든다 (Decimal로 제곱근).
p값과 신뢰구간의 임계값만 scipy 분포 함수를 쓴다. 다음으로 교차 확인한다.
- 단순회귀: scipy.stats.linregress (기울기, 절편, 표준오차, r, p)
- 모든 사례: numpy 최소제곱(lstsq)과 1e-6 안에서 일치
- 공선성 진단(고유값, 조건지수, 분산 비율): numpy.linalg.eigh 로 계산 (정확 산술 기준이 없음)
"""
import json
import os
from decimal import Decimal, getcontext
from fractions import Fraction

import numpy as np
import scipy
from scipy import stats

getcontext().prec = 60
rng = np.random.default_rng(20261025)
f = lambda v: float(v)


def solve_exact(A, b):
    """Fraction 행렬 A (p×p), b 로 A x = b 를 풀고 A^-1 도 돌려준다 (Gauss-Jordan)."""
    p = len(A)
    M = [row[:] + [Fraction(int(i == j)) for j in range(p)] + [b[i]] for i, row in enumerate(A)]
    for c in range(p):
        piv = max(range(c, p), key=lambda r: abs(M[r][c]))
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        M[c] = [v / pv for v in M[c]]
        for r in range(p):
            if r != c and M[r][c] != 0:
                fac = M[r][c]
                M[r] = [a - fac * bb for a, bb in zip(M[r], M[c])]
    inv = [row[p:2 * p] for row in M]
    x = [row[2 * p] for row in M]
    return x, inv


def dec(fr):
    return Decimal(fr.numerator) / Decimal(fr.denominator)


def exact_ols(names, X_cols, y, conf):
    """names: 예측변수 이름, X_cols: 예측변수 열(double 목록), y: double 목록."""
    n, k = len(y), len(X_cols)
    p = k + 1
    cols = [[Fraction(1)] * n] + [[Fraction(v) for v in c] for c in X_cols]
    yy = [Fraction(v) for v in y]
    A = [[sum(cols[i][t] * cols[j][t] for t in range(n)) for j in range(p)] for i in range(p)]
    rhs = [sum(cols[i][t] * yy[t] for t in range(n)) for i in range(p)]
    B, inv = solve_exact(A, rhs)
    resid = [yy[t] - sum(B[i] * cols[i][t] for i in range(p)) for t in range(n)]
    sse = sum(e * e for e in resid)
    ybar = sum(yy) / n
    sst = sum((v - ybar) ** 2 for v in yy)
    ssr = sst - sse
    dfr = n - p
    mse = sse / dfr
    s = dec(mse).sqrt()
    r2 = ssr / sst
    adj = 1 - (1 - r2) * Fraction(n - 1, dfr)
    F = (ssr / k) / mse
    # 중심화한 예측변수의 (Xc'Xc)^-1 로 VIF
    cen = []
    c2 = []
    for c in cols[1:]:
        m = sum(c) / n
        cen.append([v - m for v in c])
        c2.append(sum((v - m) ** 2 for v in c))
    Ac = [[sum(cen[i][t] * cen[j][t] for t in range(n)) for j in range(k)] for i in range(k)]
    _, invc = solve_exact(Ac, [Fraction(0)] * k)
    dw_num = sum((resid[t] - resid[t - 1]) ** 2 for t in range(1, n))
    tcrit = Decimal(float(stats.t.isf((1 - conf / 100) / 2, dfr)))

    out = {"names": ["(상수)"] + names, "B": [], "se": [], "t": [], "p": [], "lo": [], "hi": [], "beta": [None],
           "tol": [None], "vif": [None]}
    for i in range(p):
        se = s * dec(inv[i][i]).sqrt()
        b = dec(B[i])
        t = b / se
        out["B"].append(f(b)); out["se"].append(f(se)); out["t"].append(f(t))
        out["p"].append(f(2 * stats.t.sf(abs(f(t)), dfr)))
        out["lo"].append(f(b - tcrit * se)); out["hi"].append(f(b + tcrit * se))
        if i > 0:
            beta = b * (dec(c2[i - 1]) / dec(sst)).sqrt()
            vif = dec(c2[i - 1] * invc[i - 1][i - 1])
            out["beta"].append(f(beta)); out["vif"].append(f(vif)); out["tol"].append(f(1 / vif))
    out.update({
        "r2": f(r2), "r": f(dec(r2).sqrt()), "adjR2": f(adj), "s": f(s), "F": f(F),
        "pF": f(stats.f.sf(f(F), k, dfr)), "ssr": f(ssr), "sse": f(sse), "sst": f(sst),
        "dfm": k, "dfr": dfr, "dw": f(dec(dw_num) / dec(sse)),
    })
    return out


def numpy_check(X_cols, y, expected):
    X = np.column_stack([np.ones(len(y))] + [np.array(c) for c in X_cols])
    beta = np.linalg.lstsq(X, np.array(y), rcond=None)[0]
    for got, exp in zip(beta, expected["B"]):
        assert abs(got - exp) <= 1e-6 * max(abs(exp), 1e-3), (got, exp)


def diagnostics(X_cols):
    n = len(X_cols[0])
    X = np.column_stack([np.ones(n)] + [np.array(c) for c in X_cols])
    Xs = X / np.sqrt((X ** 2).sum(axis=0))
    S = Xs.T @ Xs
    w, V = np.linalg.eigh(S)
    order = np.argsort(w)[::-1]
    w, V = w[order], V[:, order]
    phi = (V ** 2) / w          # phi[i, j]: 계수 i, 차원 j
    prop = phi / phi.sum(axis=1, keepdims=True)
    return {"eig": [f(v) for v in w], "cond": [f(v) for v in np.sqrt(w[0] / w)],
            "prop": [[f(prop[i, j]) for i in range(len(w))] for j in range(len(w))]}   # prop[차원][계수]


cases = []

# 1) 단순회귀
x = rng.normal(10, 3, 25)
y = 2 + 1.5 * x + rng.normal(0, 2, 25)
lr = stats.linregress(x, y)
ex = exact_ols(["x1"], [list(x)], list(y), 95)
for got, exp in ((lr.slope, ex["B"][1]), (lr.intercept, ex["B"][0]), (lr.stderr, ex["se"][1]),
                 (lr.intercept_stderr, ex["se"][0]), (lr.rvalue, ex["r"]), (lr.pvalue, ex["p"][1])):
    assert abs(got - exp) <= 1e-9 * max(abs(exp), 1e-9), (got, exp)
numpy_check([list(x)], list(y), ex)
cases.append({"name": "단순회귀 (scipy linregress 교차 확인)", "y": [f(v) for v in y],
              "covariates": {"x1": [f(v) for v in x]}, "factor": None, "lastRef": False, "expected": ex})

# 2) 다중회귀 (공변량 3개) + 공선성 진단
n = 40
x1, x2, x3 = rng.normal(0, 1, n), rng.normal(5, 2, n), rng.uniform(0, 10, n)
y = 1 + 2 * x1 - 0.5 * x2 + 0.3 * x3 + rng.normal(0, 1.5, n)
cols = [list(x1), list(x2), list(x3)]
ex = exact_ols(["x1", "x2", "x3"], cols, list(y), 95)
numpy_check(cols, list(y), ex)
cases.append({"name": "다중회귀 (공변량 3개)", "y": [f(v) for v in y],
              "covariates": {"x1": [f(v) for v in x1], "x2": [f(v) for v in x2], "x3": [f(v) for v in x3]},
              "factor": None, "lastRef": False, "expected": ex, "diag": diagnostics(cols)})

# 3) 더미 변수 (범주 3개): 첫 범주 / 마지막 범주를 기준으로
n = 45
xc = rng.normal(50, 10, n)
g = rng.integers(1, 4, n)
y = 5 + 0.4 * xc + np.where(g == 2, 3.0, np.where(g == 3, -2.0, 0.0)) + rng.normal(0, 2, n)
for last in (False, True):
    ref = 3 if last else 1
    levels = [c for c in (1, 2, 3) if c != ref]
    dummies = [[1.0 if gv == lv else 0.0 for gv in g] for lv in levels]
    names = ["x1"] + ["g = %d" % lv for lv in levels]
    cols = [list(xc)] + dummies
    ex = exact_ols(names, cols, list(y), 99)
    numpy_check(cols, list(y), ex)
    cases.append({"name": "더미 변수 (기준 범주 %s)" % ("마지막" if last else "첫"), "y": [f(v) for v in y],
                  "covariates": {"x1": [f(v) for v in xc]},
                  "factor": {"name": "g", "values": [int(v) for v in g]}, "lastRef": last, "conf": 99, "expected": ex})

# 4) 큰 평균 (수치 안정성)
n = 30
z1 = 1e6 + rng.normal(0, 1, n)
z2 = 5e5 + rng.normal(0, 2, n)
y = 3 + 0.5 * (z1 - 1e6) - 0.2 * (z2 - 5e5) + rng.normal(0, 1, n)
cols = [list(z1), list(z2)]
ex = exact_ols(["x1", "x2"], cols, list(y), 95)
cases.append({"name": "큰 평균 (1e6, 5e5)", "y": [f(v) for v in y],
              "covariates": {"x1": [f(v) for v in z1], "x2": [f(v) for v in z2]},
              "factor": None, "lastRef": False, "expected": ex})

# 5) 강한 공선성 (VIF 큼) + 진단
n = 50
a = rng.normal(0, 1, n)
b = a + rng.normal(0, 0.05, n)
c = rng.normal(0, 1, n)
y = 1 + a + b + 0.5 * c + rng.normal(0, 1, n)
cols = [list(a), list(b), list(c)]
ex = exact_ols(["x1", "x2", "x3"], cols, list(y), 95)
numpy_check(cols, list(y), ex)
cases.append({"name": "강한 공선성", "y": [f(v) for v in y],
              "covariates": {"x1": [f(v) for v in a], "x2": [f(v) for v in b], "x3": [f(v) for v in c]},
              "factor": None, "lastRef": False, "expected": ex, "diag": diagnostics(cols)})

path = os.path.join(os.path.dirname(__file__), "reference", "regression.json")
os.makedirs(os.path.dirname(path), exist_ok=True)
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "cases": cases}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "regression cases to", path)
