"""scipy로 비모수 검정(Mann-Whitney U, Wilcoxon 부호순위, Kruskal-Wallis, Dunn) 기준값을 tests/reference/nonparam.json 에 저장한다.

실행: py -3.12 tests/gen_reference_nonparam.py
- Mann-Whitney: method='asymptotic', use_continuity=False (동순위 보정 포함), 동순위가 없으면 method='exact'도 저장
- Wilcoxon: zero_method='wilcox'(0 차이는 뺌), correction=False, method='approx'와 (동순위·0이 없으면) 'exact'
- Kruskal-Wallis: scipy 기본(동순위 보정 H, χ² 근사)
- Dunn: 동순위 보정 분산과 Bonferroni 보정을 numpy로 직접 계산 (scipy에 없음. 독립 구현으로 교차 확인)
"""
import json
import os

import numpy as np
import scipy
from scipy import stats

rng = np.random.default_rng(20261023)
f = lambda a: [float(v) for v in a]

# ---------- Mann-Whitney ----------
mw_sets = {
    "작은 표본 동순위 없음 (5, 7)": (rng.normal(0, 1, 5), rng.normal(1, 1, 7)),
    "정규 (20, 25) 동순위 없음": (rng.normal(0, 1, 20), rng.normal(0.5, 1, 25)),
    "50 대 50 동순위 없음": (rng.normal(0, 1, 50), rng.normal(0.3, 1.2, 50)),
    "크기 1 대 9": (np.array([3.3]), rng.normal(0, 1, 9)),
    "정수 1~5 동순위 많음 (30, 40)": (rng.integers(1, 6, 30).astype(float), rng.integers(1, 6, 40).astype(float)),
    "로그정규 (15, 60)": (rng.lognormal(0, 1, 15), rng.lognormal(0.4, 1, 60)),
    "큰 표본 (200, 300)": (rng.normal(0, 1, 200), rng.normal(0.1, 1, 300)),
    "큰 평균": (1e6 + rng.normal(0, 1, 12), 1e6 + rng.normal(0.5, 1, 14)),
    "완전 분리": (np.array([1.0, 2.0, 3.0, 4.0]), np.array([10.0, 11.0, 12.0])),
}
mw = []
for name, (a, b) in mw_sets.items():
    asym = stats.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic", use_continuity=False)
    ties = len(np.unique(np.concatenate([a, b]))) < len(a) + len(b)
    case = {"name": name, "a": f(a), "b": f(b), "u1": float(asym.statistic), "pAsym": float(asym.pvalue), "ties": bool(ties), "pExact": None}
    if not ties and len(a) <= 50 and len(b) <= 50:
        ex = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact")
        case["pExact"] = float(ex.pvalue)
    mw.append(case)

# ---------- Wilcoxon ----------
def wil(diffs):
    d = np.asarray(diffs, dtype=float)
    nz = d[d != 0]
    res = stats.wilcoxon(d, zero_method="wilcox", correction=False, method="approx")
    case = {"diffs": f(d), "stat": float(res.statistic), "pAsym": float(res.pvalue), "z": float(res.zstatistic), "pExact": None}
    ties = len(np.unique(np.abs(nz))) < len(nz)
    if len(nz) and not ties and len(nz) <= 50:
        case["pExact"] = float(stats.wilcoxon(d, zero_method="wilcox", method="exact").pvalue)
    return case

wil_sets = {
    "작은 표본 (n=8) 동순위·0 없음": rng.normal(0.3, 1, 8),
    "정규 n=30": rng.normal(0.2, 1, 30),
    "n=50 동순위 없음": rng.normal(0.4, 1, 50),
    "정수 차이 (동순위·0 많음)": rng.integers(-3, 5, 40).astype(float),
    "0 차이 포함, 나머지 동순위 없음": np.concatenate([np.zeros(4), rng.normal(0.5, 1, 15)]),
    "큰 표본 n=300": rng.normal(0.05, 1, 300),
    "모두 양수": np.array([0.5, 1.1, 1.7, 2.4, 3.0, 3.9, 5.1]),
}
wilcox = []
for name, d in wil_sets.items():
    c = wil(d)
    c["name"] = name
    wilcox.append(c)

# ---------- Kruskal-Wallis, Dunn ----------
def dunn(groups):
    allv = np.concatenate(groups)
    N = len(allv)
    ranks = stats.rankdata(allv)
    _, counts = np.unique(allv, return_counts=True)
    tie = float(np.sum(counts ** 3 - counts))
    pos, mean_r, ns = 0, [], []
    for g in groups:
        mean_r.append(float(ranks[pos:pos + len(g)].mean()))
        ns.append(len(g))
        pos += len(g)
    k = len(groups)
    m = k * (k - 1) // 2
    base = N * (N + 1) / 12 - tie / (12 * (N - 1))
    out = []
    for i in range(k):
        for j in range(i + 1, k):
            diff = mean_r[i] - mean_r[j]
            z = diff / np.sqrt(base * (1 / ns[i] + 1 / ns[j]))
            p = float(min(1.0, 2 * stats.norm.sf(abs(z))))
            out.append({"i": i, "j": j, "diff": float(diff), "z": float(z), "p": p, "pAdj": float(min(1.0, p * m))})
    return out

kw_sets = {
    "정규 3집단 (10, 12, 15)": [rng.normal(0, 1, 10), rng.normal(0.6, 1, 12), rng.normal(1.2, 1, 15)],
    "2집단 (20, 25)": [rng.normal(0, 1, 20), rng.normal(0.4, 1, 25)],
    "정수 1~5 동순위 많음 4집단": [rng.integers(1, 6, n).astype(float) for n in (18, 22, 15, 30)],
    "로그정규 3집단": [rng.lognormal(s, 1, n) for s, n in [(0, 30), (0.5, 40), (1.0, 20)]],
    "집단 크기 1 포함": [np.array([2.5]), rng.normal(0, 1, 9), rng.normal(1, 1, 8)],
    "큰 평균 3집단": [1e6 + rng.normal(m, 1, 20) for m in (0, 0.5, 1.5)],
}
kw = []
for name, groups in kw_sets.items():
    res = stats.kruskal(*groups)
    kw.append({"name": name, "groups": [f(g) for g in groups], "H": float(res.statistic), "p": float(res.pvalue), "dunn": dunn(groups)})

path = os.path.join(os.path.dirname(__file__), "reference", "nonparam.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "mw": mw, "wilcoxon": wilcox, "kw": kw}, fh, ensure_ascii=False, indent=0)
print("wrote", len(mw), "Mann-Whitney,", len(wilcox), "Wilcoxon,", len(kw), "Kruskal-Wallis cases to", path)
