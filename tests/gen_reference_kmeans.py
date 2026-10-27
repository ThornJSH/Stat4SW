"""k-평균 군집분석 기준값을 tests/reference/kmeans.json 에 저장한다.

실행: py -3.12 tests/gen_reference_kmeans.py
- 난수: mulberry32 를 32비트 정수 연산으로 똑같이 구현 (시드가 같으면 JS 와 같은 수열)
- 초기 중심: k-means++ (첫 중심은 floor(u·n), 다음은 D² 에 비례해 누적합 > u·합계인 첫 사례), 중심마다 난수 1개를 쓴다
- 반복: 순수 파이썬 Lloyd (가장 가까운 중심, 동점은 번호가 작은 쪽, 빈 군집은 자기 중심에서 가장 먼 사례로 채움). 같은 초기 중심으로
  scipy.cluster.vq.kmeans2(minit='matrix')를 돌려 최종 중심이 같은지 확인한다 (빈 군집이 없는 경우).
- 여러 번 시작해 군집 내 제곱합이 가장 작은 결과를 쓰고, 군집 번호는 크기 내림차순(동점은 원래 번호 순)으로 다시 매긴다.
- 변수별 F: scipy.stats.f_oneway, 평균 실루엣: numpy 거리 행렬로 직접 계산 (scikit-learn 은 설치돼 있지 않다)
"""
import json
import os

import numpy as np
import scipy
from scipy import stats
from scipy.cluster.vq import kmeans2
from scipy.spatial.distance import cdist

M32 = 0xFFFFFFFF


def mulberry32(seed):
    a = [seed & M32]

    def nxt():
        a[0] = (a[0] + 0x6D2B79F5) & M32
        t = a[0]
        t = ((t ^ (t >> 15)) * (1 | t)) & M32
        t = ((t + (((t ^ (t >> 7)) * (61 | t)) & M32)) & M32) ^ t
        return ((t ^ (t >> 14)) & M32) / 4294967296

    return nxt


def kpp(X, K, rnd):
    n = len(X)
    chosen = [int(np.floor(rnd() * n))]
    d2 = ((X - X[chosen[0]]) ** 2).sum(1)
    for _ in range(1, K):
        total = float(d2.sum())
        if total == 0:
            pick = int(np.floor(rnd() * n))
        else:
            u = rnd() * total
            cum = 0.0
            pick = n - 1
            for i in range(n):
                cum += float(d2[i])
                if cum > u:
                    pick = i
                    break
        chosen.append(pick)
        d2 = np.minimum(d2, ((X - X[pick]) ** 2).sum(1))
    return chosen


def lloyd(X, centers, max_iter):
    n, k = X.shape
    K = len(centers)
    centers = centers.copy()
    labels = np.full(n, -1)
    conv = False
    it = 0
    for it in range(1, max_iter + 1):
        d = ((X[:, None, :] - centers[None, :, :]) ** 2).sum(2)
        new = np.array([int(np.argmin(d[i])) for i in range(n)])      # argmin: 동점은 번호가 작은 쪽
        changed = bool((new != labels).any())
        labels = new
        dist = d[np.arange(n), labels].copy()
        counts = np.bincount(labels, minlength=K)
        for j in range(K):
            if counts[j] == 0:
                cand = [i for i in range(n) if counts[labels[i]] > 1]
                i = max(cand, key=lambda q: (dist[q], -q))
                counts[labels[i]] -= 1
                labels[i] = j
                counts[j] = 1
                dist[i] = 0.0
                changed = True
        for j in range(K):
            centers[j] = X[labels == j].mean(0)
        if not changed:
            conv = True
            break
    wss = float(sum(((X[labels == j] - centers[j]) ** 2).sum() for j in range(K)))
    return labels, centers, wss, it, conv


def run(X, K, seed, starts, max_iter=100):
    rnd = mulberry32(seed)
    best = None
    for s in range(starts):
        idx = kpp(X, K, rnd)
        lab, cen, wss, it, conv = lloyd(X, X[idx].copy(), max_iter)
        if best is None or wss < best["wss"]:
            best = dict(start=s + 1, init=idx, labels=lab, centers=cen, wss=wss, iters=it, converged=conv)
    # 크기 내림차순(동점은 원래 번호 순)으로 번호를 다시 매긴다
    sizes = np.bincount(best["labels"], minlength=K)
    order = sorted(range(K), key=lambda j: (-sizes[j], j))
    remap = {old: new for new, old in enumerate(order)}
    best["labels"] = np.array([remap[l] for l in best["labels"]])
    best["centers"] = best["centers"][order]
    return best


def silhouette(X, labels):
    n = len(X)
    D = cdist(X, X)
    K = labels.max() + 1
    s = np.zeros(n)
    for i in range(n):
        own = labels == labels[i]
        if own.sum() <= 1:
            s[i] = 0.0
            continue
        a = D[i, own].sum() / (own.sum() - 1)
        b = min(D[i, labels == j].mean() for j in range(K) if j != labels[i])
        s[i] = (b - a) / max(a, b)
    return float(s.mean())


def zscore(X):
    return (X - X.mean(0)) / X.std(0, ddof=1)


rng = np.random.default_rng(20261033)
n = 90
blob = rng.integers(0, 3, n)
means = np.array([[30.0, 300.0, 60.0], [45.0, 650.0, 80.0], [38.0, 420.0, 95.0]])
sds = np.array([4.0, 60.0, 6.0])
X = np.round(means[blob] + rng.normal(0, 1, (n, 3)) * sds, 2)

cases = []


def add(name, K, seed, starts, std, elbow=False):
    W = zscore(X) if std == "z" else X.copy()
    r = run(W, K, seed, starts)
    # scipy 와 같은 초기 중심에서 같은 최종 중심인지 (빈 군집이 없는 한 같다) — 한 번 시작한 경우만 확인
    if starts == 1:
        init = W[r["init"]].copy()
        c2, l2 = kmeans2(W, init, iter=100, minit="matrix")
        order = sorted(range(K), key=lambda j: (-np.bincount(l2, minlength=K)[j], j))
        assert np.allclose(c2[order], r["centers"], atol=1e-9), "scipy 와 최종 중심이 다름"
    tss = float(((W - W.mean(0)) ** 2).sum())
    sizes = [int(v) for v in np.bincount(r["labels"], minlength=K)]
    groups = [X[r["labels"] == j] for j in range(K)]
    fs = []
    for v in range(X.shape[1]):
        f, p = stats.f_oneway(*[g[:, v] for g in groups])
        fs.append(dict(F=float(f), p=float(p), df1=K - 1, df2=n - K))
    dists = np.sqrt(((W - r["centers"][r["labels"]]) ** 2).sum(1))
    cdm = cdist(r["centers"], r["centers"])
    item = dict(name=name, K=K, seed=seed, starts=starts, standardize=std, init=[int(i) for i in r["init"]], best_start=r["start"],
                labels=[int(v) for v in r["labels"]], centers=[[float(v) for v in row] for row in r["centers"]],
                centers_orig=[[float(v) for v in X[r["labels"] == j].mean(0)] for j in range(K)],
                wss=r["wss"], tss=tss, iters=r["iters"], converged=r["converged"], sizes=sizes,
                wss_by_cluster=[float(((W[r["labels"] == j] - r["centers"][j]) ** 2).sum()) for j in range(K)],
                dist_to_center=[float(v) for v in dists], center_dist=[[float(v) for v in row] for row in cdm], f=fs,
                silhouette=silhouette(W, r["labels"]))
    if elbow:
        item["elbow"] = []
        for k2 in range(2, 11):
            r2 = run(W, k2, seed, starts)
            item["elbow"].append(dict(k=k2, wss=r2["wss"], silhouette=silhouette(W, r2["labels"])))
    cases.append(item)


add("z점수, k=3, 10회 시작, 엘보 포함", 3, 12345, 10, "z", elbow=True)
add("z점수, k=4, 한 번 시작(scipy와 대조)", 4, 7, 1, "z")
add("표준화 없음, k=3, 5회 시작", 3, 2026, 5, "none")

path = os.path.join(os.path.dirname(__file__), "reference", "kmeans.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "data": [[float(v) for v in r] for r in X], "cases": cases}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "kmeans cases; sizes:", [c["sizes"] for c in cases], "iters:", [c["iters"] for c in cases], "conv:", [c["converged"] for c in cases],
      "explained:", [round(1 - c["wss"] / c["tss"], 3) for c in cases])
