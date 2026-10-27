"""Friedman 검정 기준값을 tests/reference/friedman.json 에 저장한다.

실행: py -3.12 tests/gen_reference_friedman.py
- 검정 통계량과 유의확률: scipy.stats.friedmanchisquare (동순위 보정 포함, χ² 근사, 자유도 k − 1)
- 평균 순위: scipy.stats.rankdata(method='average') 를 사람(행)마다 적용
- Kendall W = χ² / (n (k − 1))
- Cochran Q: 0/1 자료에서 공식 Q = (k − 1)[k ΣC_j² − T²] / (k T − ΣR_i²) 와 Friedman χ² 가 같음을 확인해 저장
- 사후비교(Dunn–Bonferroni): z = (평균 순위 I − 평균 순위 J) / √( k(k + 1)/(6n) · C ), C = 1 − Σ_i Σ_t (t³ − t) / (n k (k² − 1)) (동순위 보정),
  p = 2·(1 − Φ(|z|)), 보정 p = min(1, p × k(k − 1)/2)
"""
import json
import os
from itertools import combinations

import numpy as np
import scipy
from scipy import stats

rng = np.random.default_rng(20261037)


def tie_sum(row):
    _, counts = np.unique(row, return_counts=True)
    return float(((counts ** 3) - counts).sum())


def analyse(name, X):
    X = np.asarray(X, float)
    n, k = X.shape
    ranks = np.array([stats.rankdata(r) for r in X])
    mean_ranks = ranks.mean(0)
    tie_total = sum(tie_sum(r) for r in X)
    C = 1 - tie_total / (n * k * (k * k - 1))
    if C > 0:
        res = stats.friedmanchisquare(*[X[:, j] for j in range(k)])
        chi2, p = float(res.statistic), float(res.pvalue)
        w = chi2 / (n * (k - 1))
    else:
        chi2 = p = w = None
    pairs = []
    m = k * (k - 1) / 2
    if C > 0:
        se = np.sqrt(k * (k + 1) / (6 * n) * C)
        for i, j in combinations(range(k), 2):
            diff = mean_ranks[i] - mean_ranks[j]
            z = diff / se
            pv = float(2 * stats.norm.sf(abs(z)))
            pairs.append(dict(i=i, j=j, diff=float(diff), z=float(z), p=pv, p_adj=min(1.0, pv * m)))
    out = dict(name=name, n=n, k=k, data=X.tolist(), medians=[float(np.median(X[:, j])) for j in range(k)],
               mean_ranks=mean_ranks.tolist(), rank_sums=ranks.sum(0).tolist(), tie_total=tie_total, chi2=chi2, p=p, w=w, pairs=pairs)
    if np.isin(X, [0, 1]).all():
        T = X.sum()
        Cj = X.sum(0)
        Ri = X.sum(1)
        denom = k * T - (Ri ** 2).sum()
        out["cochran_q"] = float((k - 1) * (k * (Cj ** 2).sum() - T ** 2) / denom) if denom > 0 else None
    return out


cases = []
# 동순위가 없는 연속 자료: 20명, 4회 측정 (시간이 갈수록 조금씩 오름)
base = rng.normal(50, 8, (20, 1))
X1 = np.round(base + np.array([0, 2.0, 4.5, 5.0]) + rng.normal(0, 4, (20, 4)), 3)
cases.append(analyse("연속 자료 20명 × 4회 (동순위 없음)", X1))
# 리커트 정수 자료(동순위 많음): 30명, 3회
base = rng.integers(2, 5, (30, 1))
X2 = np.clip(base + rng.integers(-1, 2, (30, 3)) + np.array([0, 0, 1]), 1, 5)
cases.append(analyse("5점 척도 30명 × 3회 (동순위 많음)", X2))
# 5회 측정 12명
X3 = np.round(rng.normal(0, 1, (12, 1)) + np.array([0, 0.5, 1.0, 1.0, 2.0]) + rng.normal(0, 1, (12, 5)), 2)
cases.append(analyse("12명 × 5회", X3))
# 한 사람은 모든 값이 같음 (그 행의 순위는 모두 같다)
X4 = X1[:10, :3].copy()
X4[3, :] = X4[3, 0]
cases.append(analyse("10명 × 3회, 한 사람은 모든 값이 같음", X4))
# 0/1 자료: Cochran Q 와 같은 값
X5 = (rng.random((24, 1)) * 0.6 + rng.random((24, 3)) * 0.6 + np.array([0, 0.15, 0.3]) > 0.8).astype(int)
cases.append(analyse("0/1 자료 24명 × 3회 (Cochran Q와 동일)", X5))
# 모든 값이 같음: 검정 불가
cases.append(analyse("모든 값이 같음", np.full((6, 3), 4.0)))

path = os.path.join(os.path.dirname(__file__), "reference", "friedman.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "cases": cases}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "friedman cases;", [(c["n"], c["k"], None if c["chi2"] is None else round(c["chi2"], 3)) for c in cases], "cochran:", [c.get("cochran_q") for c in cases][4])
