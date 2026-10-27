"""상관분석과 교차분석의 기준값을 tests/reference/relation.json 에 저장한다.

실행: py -3.12 tests/gen_reference_relation.py

상관분석: scipy.stats.pearsonr(통계량, p, confidence_interval), scipy.stats.spearmanr
교차분석: scipy.stats.chi2_contingency(Pearson, 우도비, Yates 연속성 수정), fisher_exact,
          scipy.stats.contingency.association(Cramér V, 분할계수), 선형 대 선형은 확장한 자료의 피어슨 r 로 (N-1)·r²
"""
import json
import os

import numpy as np
import scipy
from scipy import stats
from scipy.stats.contingency import association

rng = np.random.default_rng(20261024)
f = lambda v: float(v)

# ---------- 상관 ----------
pairs = []
for name, x, y in [
    ("정규 양의 상관 n=30", *(lambda a: (a, a * 0.7 + rng.normal(0, 1, 30)))(rng.normal(0, 1, 30))),
    ("정수 점수 동순위 많음 n=50", *(lambda a: (a, np.clip(a + rng.integers(-2, 3, 50), 1, 7)))(rng.integers(1, 8, 50).astype(float))),
    ("음의 상관 n=12", *(lambda a: (a, -a + rng.normal(0, 0.5, 12)))(rng.normal(5, 2, 12))),
    ("큰 평균 n=40", *(lambda a: (a + 1e6, a * 0.5 + rng.normal(0, 1, 40) + 1e6))(rng.normal(0, 1, 40))),
    ("n=3", np.array([1.0, 2.0, 4.0]), np.array([2.0, 1.0, 5.0])),
]:
    entry = {"name": name, "x": [f(v) for v in x], "y": [f(v) for v in y]}
    pr = stats.pearsonr(x, y)
    sp = stats.spearmanr(x, y)
    entry["pearson"] = {"r": f(pr.statistic), "p": f(pr.pvalue)}
    entry["spearman"] = {"r": f(sp.statistic), "p": f(sp.pvalue)}
    entry["ci"] = {}
    if len(x) > 3:
        for conf in (95, 99):
            ci = pr.confidence_interval(conf / 100)
            entry["ci"][str(conf)] = {"lo": f(ci.low), "hi": f(ci.high)}
    pairs.append(entry)

# 결측이 있는 4변수 자료로 쌍별/목록별 행렬 검증
n = 60
base = rng.normal(0, 1, n)
cols = {
    "a": base + rng.normal(0, 0.8, n),
    "b": base * 0.6 + rng.normal(0, 1, n),
    "c": rng.integers(1, 6, n).astype(float),
    "d": -base + rng.normal(0, 1.2, n),
}
mask = {k: rng.random(n) < 0.12 for k in cols}   # 변수별 결측
data = {k: [None if mask[k][i] else f(cols[k][i]) for i in range(n)] for k in cols}
names = list(cols)


def matrix(listwise):
    out = {"pearson": {}, "spearman": {}}
    if listwise:
        ok = [all(not mask[k][i] for k in names) for i in range(n)]
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if j <= i:
                continue
            sel = ok if listwise else [(not mask[a][t]) and (not mask[b][t]) for t in range(n)]
            xs = np.array([cols[a][t] for t in range(n) if sel[t]])
            ys = np.array([cols[b][t] for t in range(n) if sel[t]])
            pr, sp = stats.pearsonr(xs, ys), stats.spearmanr(xs, ys)
            out["pearson"][f"{a},{b}"] = {"n": len(xs), "r": f(pr.statistic), "p": f(pr.pvalue)}
            out["spearman"][f"{a},{b}"] = {"n": len(xs), "r": f(sp.statistic), "p": f(sp.pvalue)}
    return out


matrix_case = {"names": names, "data": data, "pairwise": matrix(False), "listwise": matrix(True)}

# ---------- 교차분석 ----------
tables = []
for name, t in [
    ("2x2", [[12, 5], [7, 16]]),
    ("2x2 영 셀", [[3, 0], [1, 9]]),
    ("2x3", [[10, 20, 15], [25, 15, 5]]),
    ("3x4", [[5, 8, 9, 2], [6, 7, 4, 10], [12, 3, 5, 9]]),
    ("5x2 영 셀", [[0, 4], [3, 3], [8, 2], [10, 0], [6, 6]]),
]:
    obs = np.array(t, float)
    N = obs.sum()
    r, c = obs.shape
    chi = stats.chi2_contingency(obs, correction=False)
    lr = stats.chi2_contingency(obs, correction=False, lambda_="log-likelihood")
    entry = {
        "name": name, "table": t, "N": int(N),
        "chi2": f(chi.statistic), "p": f(chi.pvalue), "dof": int(chi.dof),
        "expected": [[f(v) for v in row] for row in chi.expected_freq],
        "lr": f(lr.statistic), "lr_p": f(lr.pvalue),
        "cramer": f(association(np.array(t), method="cramer")),
        "contingency": f(association(np.array(t), method="pearson")),
    }
    # 조정된 표준화 잔차
    R = obs.sum(axis=1, keepdims=True)
    C = obs.sum(axis=0, keepdims=True)
    E = chi.expected_freq
    adj = (obs - E) / np.sqrt(E * (1 - R / N) * (1 - C / N))
    entry["adj_resid"] = [[f(v) for v in row] for row in adj]
    # 선형 대 선형 (행·열 점수 = 1, 2, 3, ...)
    xs, ys = [], []
    for i in range(r):
        for j in range(c):
            xs += [i + 1] * int(obs[i, j])
            ys += [j + 1] * int(obs[i, j])
    rr = np.corrcoef(xs, ys)[0, 1]
    m2 = (N - 1) * rr ** 2
    entry["lbl"] = {"value": f(m2), "p": f(stats.chi2.sf(m2, 1))}
    if (r, c) == (2, 2):
        yates = stats.chi2_contingency(obs, correction=True)
        entry["yates"] = {"value": f(yates.statistic), "p": f(yates.pvalue)}
        a, b, cc, d = obs[0, 0], obs[0, 1], obs[1, 0], obs[1, 1]
        phi = np.sqrt(chi.statistic / N)
        entry["phi"] = f(phi if a * d - b * cc >= 0 else -phi)
        two = stats.fisher_exact(obs, alternative="two-sided").pvalue
        less = stats.fisher_exact(obs, alternative="less").pvalue
        greater = stats.fisher_exact(obs, alternative="greater").pvalue
        entry["fisher"] = {"two": f(two), "one": f(min(less, greater))}
    tables.append(entry)

path = os.path.join(os.path.dirname(__file__), "reference", "relation.json")
os.makedirs(os.path.dirname(path), exist_ok=True)
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "pairs": pairs, "matrix": matrix_case, "tables": tables},
              fh, ensure_ascii=False, indent=0)
print("wrote", len(pairs), "correlation pairs,", 1, "matrix case,", len(tables), "tables to", path)
