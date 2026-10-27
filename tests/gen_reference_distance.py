"""거리·유사도 기준값을 tests/reference/distance.json 에 저장한다.

실행: py -3.12 tests/gen_reference_distance.py
- 거리: scipy.spatial.distance (pdist/cdist) — 이 도구의 직접 구현과 다른 계산 경로
  euclidean, sqeuclidean, cityblock, chebyshev, minkowski(p), cosine(1 − 유사도), correlation(1 − r),
  mahalanobis(공분산 역행렬 VI 를 numpy 로 직접 계산해 넘김), 이분 자료 hamming/jaccard/dice
- 표준화: numpy 로 직접 (표준편차는 n−1)
- 마할라노비스 이상치: numpy 로 D² = (x−x̄)ᵀ S⁻¹ (x−x̄), 유의확률은 scipy chi2.sf, 임계값은 chi2.isf
"""
import json
import os

import numpy as np
import scipy
from scipy import stats
from scipy.spatial import distance as sd

rng = np.random.default_rng(20261031)


def standardize(X, kind):
    X = np.asarray(X, float)
    if kind == "none":
        return X
    if kind == "z":
        return (X - X.mean(0)) / X.std(0, ddof=1)
    if kind == "range01":
        return (X - X.min(0)) / (X.max(0) - X.min(0))
    if kind == "range11":
        return 2 * (X - X.min(0)) / (X.max(0) - X.min(0)) - 1
    if kind == "max1":
        return X / X.max(0)
    if kind == "mean1":
        return X / X.mean(0)
    if kind == "sd1":
        return X / X.std(0, ddof=1)
    raise ValueError(kind)


def full(v):
    return [[float(x) for x in row] for row in sd.squareform(v)]


def matrix(measure, V, p=3, VI=None):
    """V: (m, q) 벡터들. 거리 또는 유사도 행렬 (대각선 포함)"""
    if measure == "euclid":
        return full(sd.pdist(V, "euclidean"))
    if measure == "sqeuclid":
        return full(sd.pdist(V, "sqeuclidean"))
    if measure == "cityblock":
        return full(sd.pdist(V, "cityblock"))
    if measure == "chebyshev":
        return full(sd.pdist(V, "chebyshev"))
    if measure == "minkowski":
        return full(sd.pdist(V, "minkowski", p=p))
    if measure == "mahalanobis":
        return full(sd.pdist(V, "mahalanobis", VI=VI))
    if measure == "cosineDist":
        return full(sd.pdist(V, "cosine"))
    if measure == "cosine":                       # 유사도 (대각선 1)
        m = 1 - sd.squareform(sd.pdist(V, "cosine"))
        np.fill_diagonal(m, 1.0)
        return [[float(x) for x in r] for r in m]
    if measure == "corrDist":
        return full(sd.pdist(V, "correlation"))
    if measure == "correlation":
        m = 1 - sd.squareform(sd.pdist(V, "correlation"))
        np.fill_diagonal(m, 1.0)
        return [[float(x) for x in r] for r in m]
    if measure in ("matching", "jaccard", "dice"):
        name = {"matching": "hamming", "jaccard": "jaccard", "dice": "dice"}[measure]
        m = 1 - sd.squareform(sd.pdist(V.astype(bool), name))
        np.fill_diagonal(m, 1.0)
        return [[float(x) for x in r] for r in m]
    if measure == "hamming":
        return full(sd.pdist(V.astype(bool), "hamming"))
    raise ValueError(measure)


pairs = []


def add(name, target, measure, X, std, p=3, binary=False):
    Xs = standardize(X, std)
    V = Xs if target == "case" else Xs.T
    VI = None
    if measure == "mahalanobis":
        VI = np.linalg.inv(np.cov(Xs, rowvar=False))
    pairs.append(dict(name=name, target=target, measure=measure, standardize=std, p=p, data=[[float(x) for x in r] for r in X],
                      expected=matrix(measure, V, p=p, VI=VI), binary=binary))


# ---- 연속 자료: 단위가 서로 다른 4개 변수, 12사례 (모두 양수) ----
n = 12
X = np.column_stack([
    np.round(rng.uniform(20, 60, n)),                   # 연령
    np.round(rng.uniform(200, 900, n)),                 # 소득
    np.round(rng.normal(75, 10, n), 1),                 # 점수
    rng.integers(1, 6, n).astype(float),                # 만족도
])
for m in ["euclid", "sqeuclid", "cityblock", "chebyshev", "minkowski", "mahalanobis", "cosineDist", "cosine", "corrDist", "correlation"]:
    add(f"사례 간, 표준화 없음, {m}", "case", m, X, "none")
    add(f"사례 간, z점수, {m}", "case", m, X, "z")
for std in ["range01", "range11", "max1", "mean1", "sd1"]:
    add(f"사례 간, {std}, euclid", "case", "euclid", X, std)
for m in ["euclid", "sqeuclid", "cityblock", "chebyshev", "minkowski", "cosineDist", "cosine", "corrDist", "correlation"]:
    add(f"변수 간, z점수, {m}", "var", m, X, "z")
add("사례 간, 민코프스키 p=1.5", "case", "minkowski", X, "z", p=1.5)

# ---- 이분 자료: 10사례 × 8변수 (모든 행에 1이 하나 이상, 같은 행이 겹치지 않게) ----
while True:
    B = (rng.random((10, 8)) < 0.45).astype(float)
    if (B.sum(1) > 0).all() and len({tuple(r) for r in B}) == 10:
        break
for m in ["matching", "jaccard", "dice", "hamming"]:
    add(f"이분 자료, 사례 간, {m}", "case", m, B, "none", binary=True)

# ---- 마할라노비스 이상치: 40사례, 3변수, 상관이 있는 정규 자료 + 한 사례(각 변수는 정상 범위, 조합이 이상) ----
n = 40
cov = np.array([[1.0, 0.8, 0.5], [0.8, 1.0, 0.4], [0.5, 0.4, 1.0]])
Z = rng.multivariate_normal([50, 60, 70], cov * np.array([[100, 90, 40], [90, 100, 30], [40, 30, 64]]) / 1.0 * 0.1, n)
Z[7] = [Z[:, 0].mean() + 8.0, Z[:, 1].mean() - 8.0, Z[:, 2].mean()]
Z = np.round(Z, 2)
S = np.cov(Z, rowvar=False)
Sinv = np.linalg.inv(S)
dev = Z - Z.mean(0)
D2 = np.einsum("ij,jk,ik->i", dev, Sinv, dev)
k = Z.shape[1]
maha = dict(data=[[float(x) for x in r] for r in Z], k=k, n=n, d2=[float(x) for x in D2],
            p=[float(x) for x in stats.chi2.sf(D2, k)], critical=float(stats.chi2.isf(0.001, k)),
            flagged=[int(i) for i in np.where(stats.chi2.sf(D2, k) < 0.001)[0]], top=int(np.argmax(D2)))

path = os.path.join(os.path.dirname(__file__), "reference", "distance.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "pairs": pairs, "maha": maha}, fh, ensure_ascii=False, indent=0)
print("wrote", len(pairs), "matrices to", path, "| maha top", maha["top"], "D2", round(max(D2), 3), "flagged", maha["flagged"])
