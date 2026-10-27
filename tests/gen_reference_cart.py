"""의사결정나무(CART) 기준값을 tests/reference/cart.json 에 저장한다.

실행: py -3.12 tests/gen_reference_cart.py
scikit-learn 이 설치돼 있지 않아 numpy·순수 파이썬으로 독립 구현한다. 이 도구의 JS 구현(접두합으로 증분 계산)과 달리
모든 분할 후보를 처음부터 다시 계산한다 (같은 규칙, 다른 계산 경로). 규칙:

- 불순도: 분류는 지니 1 − Σp², 회귀는 SSE/n. 노드의 위험 R(t): 분류는 오분류 건수, 회귀는 SSE. 분할 이득 Δ = n·i(t) − nL·i(L) − nR·i(R).
- 후보: 숫자 변수는 노드 안 서로 다른 값의 중간값(≤ 왼쪽), 범주 변수는 (회귀, 2분류) 응답 통계(평균, 둘째 범주 비율)로 범주를 정렬해 앞에서부터 자른 집합,
  3범주 이상 분류는 첫 범주를 포함하는 모든 부분집합(마스크 오름차순). 양쪽이 최소 잎 크기 이상이어야 한다.
- 최고 이득: 변수 순서(숫자 변수들, 범주 변수들), 후보 순서대로 보며 a − b > 1e-10·max(1, |b|) 일 때만 바꾼다 (동점은 앞선 후보).
- 중지: n < minSplit, 깊이 ≥ maxDepth, 불순도 ≤ 1e-12, 유효한 분할 없음 또는 이득 ≤ 1e-12. 노드 번호는 전위 순회(루트 1).
- 가지치기: 약한 고리 g(t) = (R(t) − R(T_t)) / (잎 수 − 1)가 최소인 노드들을 번호 순으로 한꺼번에 접는다 (이미 접힌 조상 아래는 건너뜀).
  부분나무 S_0(전체), S_1, …, S_m(루트), α_k = S_{k−1} → S_k 로 갈 때의 g.
- 교차검증: 시드 mulberry32 로 인덱스를 Fisher–Yates 섞고 위치 % 폴드 수로 폴드를 정한다. 폴드마다 훈련 사례로 나무를 키우고 자기 부분나무열을 만든 뒤,
  β_0 = 0, β_k = √(α_k·α_{k+1}), β_m = ∞ 에서 α'_i ≤ β 인 i의 수 j 로 S'_j 를 골라 검증 사례를 예측한다 (범주를 못 본 경우 사례가 많은 쪽 자식, 같으면 왼쪽).
  xerror_k = 총 손실/n (분류 0/1, 회귀 제곱오차), se_k = 사례별 손실의 표본 표준편차/√n. 선택: cv = 최소(동점은 작은 나무), cv1se = 최소+se 이하 중 가장 작은 나무.
"""
import json
import math
import os

import numpy as np

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


def better(a, b):
    return a - b > 1e-10 * max(1.0, abs(b))


class Tree:
    def __init__(self, X, types, y, mode, C, params):
        self.X, self.types, self.y, self.mode, self.C, self.p = X, types, y, mode, C, params
        self.nodes = []

    # ---- 노드 통계 ----
    def stats(self, idx):
        n = len(idx)
        if self.mode == "class":
            counts = np.bincount(self.y[idx].astype(int), minlength=self.C)
            imp = 1 - float(((counts / n) ** 2).sum())
            pred = int(np.argmax(counts))
            err = float(n - counts.max())
            return dict(n=n, counts=[int(c) for c in counts], imp=imp, pred=pred, err=err, nimp=n * imp)
        yy = self.y[idx]
        m = float(yy.mean())
        sse = float(((yy - m) ** 2).sum())
        return dict(n=n, mean=m, imp=sse / n, pred=m, err=sse, nimp=sse)

    def nimp(self, idx):
        return self.stats(idx)["nimp"]

    def candidates(self, idx, j):
        x = self.X[j][idx]
        if self.types[j] == "num":
            vals = sorted(set(float(v) for v in x))
            for a, b in zip(vals[:-1], vals[1:]):
                thr = (a + b) / 2
                yield ("num", thr, idx[x <= thr], idx[x > thr])
        else:
            cats = sorted(set(int(v) for v in x))
            m = len(cats)
            if m < 2:
                return
            if self.mode == "reg" or self.C == 2:
                stat = {}
                for c in cats:
                    sel = self.y[idx[x == c]]
                    stat[c] = float(sel.mean()) if self.mode == "reg" else float((sel == 1).mean())
                order = sorted(cats, key=lambda c: (stat[c], c))
                subsets = [order[:t] for t in range(1, m)]
            else:
                subsets = []
                for mask in range(2 ** (m - 1) - 1):
                    s = [cats[0]] + [cats[i + 1] for i in range(m - 1) if (mask >> i) & 1]
                    subsets.append(s)
            for s in subsets:
                inl = np.isin(x, s)
                yield ("cat", sorted(s), idx[inl], idx[~inl])

    def best_split(self, idx, parent_nimp):
        best = None
        for j in range(len(self.X)):
            for kind, arg, left, right in self.candidates(idx, j):
                if len(left) < self.p["minLeaf"] or len(right) < self.p["minLeaf"]:
                    continue
                gain = parent_nimp - self.nimp(left) - self.nimp(right)
                if best is None or better(gain, best["gain"]):
                    best = dict(j=j, kind=kind, arg=arg, gain=gain, left=left, right=right)
        return best

    def grow(self, idx, depth):
        st = self.stats(idx)
        node = dict(id=len(self.nodes) + 1, depth=depth, left=0, right=0, split=None, **st)
        self.nodes.append(node)
        if st["n"] < self.p["minSplit"] or depth >= self.p["maxDepth"] or st["imp"] <= 1e-12:
            return node
        sp = self.best_split(idx, st["nimp"])
        if sp is None or sp["gain"] <= 1e-12:
            return node
        node["split"] = dict(j=sp["j"], kind=sp["kind"], arg=sp["arg"], gain=sp["gain"], major_left=len(sp["left"]) >= len(sp["right"]))
        node["left"] = self.grow(sp["left"], depth + 1)["id"]
        node["right"] = self.grow(sp["right"], depth + 1)["id"]
        return node

    # ---- 예측 ----
    def goes_left(self, node, row):
        s = node["split"]
        v = row[s["j"]]
        if s["kind"] == "num":
            return v <= s["arg"]
        if int(v) in s["arg"]:
            return True
        known = s.get("known")
        if known is not None and int(v) not in known:
            return s["major_left"]
        return False

    def predict_leaf(self, row, collapsed):
        node = self.nodes[0]
        while node["split"] is not None and node["id"] not in collapsed:
            node = self.nodes[(node["left"] if self.goes_left(node, row) else node["right"]) - 1]
        return node


def mark_known(tree, idx_all_by_node=None):
    """학습 때 노드에 있던 범주를 기록해, 못 본 범주는 사례가 많은 쪽으로 보낸다 (분기 노드마다)."""
    def rec(node, idx):
        s = node["split"]
        if s is None:
            return
        x = tree.X[s["j"]][idx]
        if s["kind"] == "cat":
            s["known"] = set(int(v) for v in x)
            inl = np.isin(x, s["arg"])
        else:
            inl = x <= s["arg"]
        rec(tree.nodes[node["left"] - 1], idx[inl])
        rec(tree.nodes[node["right"] - 1], idx[~inl])
    rec(tree.nodes[0], np.arange(len(tree.y)) if idx_all_by_node is None else idx_all_by_node)


def build(X, types, y, mode, C, params, idx):
    t = Tree(X, types, y, mode, C, params)
    t.grow(idx, 0)
    mark_known(t, idx)
    return t


def reach(tree, collapsed):
    seen = []
    def rec(node):
        seen.append(node["id"])
        if node["split"] is not None and node["id"] not in collapsed:
            rec(tree.nodes[node["left"] - 1])
            rec(tree.nodes[node["right"] - 1])
    rec(tree.nodes[0])
    return seen


def subtree_stats(tree, node, collapsed):
    if node["split"] is None or node["id"] in collapsed:
        return node["err"], 1
    r1, l1 = subtree_stats(tree, tree.nodes[node["left"] - 1], collapsed)
    r2, l2 = subtree_stats(tree, tree.nodes[node["right"] - 1], collapsed)
    return r1 + r2, l1 + l2


def prune_sequence(tree):
    collapsed = set()
    alphas, sets, leaves = [], [set()], [subtree_stats(tree, tree.nodes[0], collapsed)[1]]
    while True:
        ids = reach(tree, collapsed)
        internal = [i for i in sorted(ids) if tree.nodes[i - 1]["split"] is not None and i not in collapsed]
        if not internal:
            break
        g = {}
        for i in internal:
            node = tree.nodes[i - 1]
            rs, ls = subtree_stats(tree, node, collapsed)
            g[i] = max(0.0, (node["err"] - rs) / (ls - 1))
        gmin = min(g.values())
        for i in internal:
            if g[i] <= gmin + 1e-12 * max(1.0, abs(gmin)) and i in reach_set(tree, collapsed):
                collapsed.add(i)
        alphas.append(gmin)
        sets.append(set(collapsed))
        leaves.append(subtree_stats(tree, tree.nodes[0], collapsed)[1])
    return alphas, sets, leaves


def reach_set(tree, collapsed):
    return set(reach(tree, collapsed))


def betas(alphas):
    m = len(alphas)
    out = [0.0]
    for k in range(1, m):
        a, b = alphas[k - 1], alphas[k]
        out.append(math.sqrt(a * b) if a > 0 else 0.0)
    out.append(math.inf)
    return out[: m + 1]


def loss(tree, collapsed, row, yv, mode):
    leaf = tree.predict_leaf(row, collapsed)
    return (0.0 if int(leaf["pred"]) == int(yv) else 1.0) if mode == "class" else float((leaf["pred"] - yv) ** 2)


def cross_validate(X, types, y, mode, C, params, full_alphas, seed, nfold):
    n = len(y)
    rnd = mulberry32(seed)
    order = list(range(n))
    for i in range(n - 1, 0, -1):
        j = int(math.floor(rnd() * (i + 1)))
        order[i], order[j] = order[j], order[i]
    fold = [0] * n
    for pos, i in enumerate(order):
        fold[i] = pos % nfold
    bt = betas(full_alphas)
    losses = np.zeros((len(bt), n))
    for f in range(nfold):
        test = [i for i in range(n) if fold[i] == f]
        train = np.array([i for i in range(n) if fold[i] != f])
        tr = build(X, types, y, mode, C, params, train)
        al, sets, _ = prune_sequence(tr)
        for k, beta in enumerate(bt):
            j = sum(1 for a in al if a <= beta)
            for i in test:
                row = [X[c][i] for c in range(len(X))]
                losses[k, i] = loss(tr, sets[j], row, y[i], mode)
    xerr = losses.sum(1) / n
    se = losses.std(1, ddof=1) / math.sqrt(n)
    return xerr, se, bt


def choose(xerr, se, rule):
    m = len(xerr) - 1
    kmin = min(range(m + 1), key=lambda k: (xerr[k], -k))          # 동점은 작은 나무(큰 k)
    if rule == "cv":
        return kmin
    thr = xerr[kmin] + se[kmin]
    return max(k for k in range(m + 1) if xerr[k] <= thr + 1e-12)


def export_nodes(tree, collapsed):
    """보이는 나무를 전위 순회로 다시 번호 매긴다"""
    out, newid = [], {}

    def rec(node, parent):
        nid = len(out) + 1
        newid[node["id"]] = nid
        rec_node = dict(id=nid, old=node["id"], depth=node["depth"], n=node["n"], pred=node["pred"], imp=node["imp"], err=node["err"], left=0, right=0, split=None)
        if "counts" in node:
            rec_node["counts"] = node["counts"]
        out.append(rec_node)
        if node["split"] is not None and node["id"] not in collapsed:
            s = node["split"]
            rec_node["split"] = dict(j=s["j"], kind=s["kind"], arg=s["arg"] if s["kind"] == "num" else [int(v) for v in s["arg"]], gain=s["gain"])
            rec_node["left"] = rec(tree.nodes[node["left"] - 1], nid)
            rec_node["right"] = rec(tree.nodes[node["right"] - 1], nid)
        return nid

    rec(tree.nodes[0], 0)
    return out


def run_case(name, cols, types, y, mode, C, params, rule, seed, nfold, var_names):
    X = [np.asarray(c, float) for c in cols]
    y = np.asarray(y, float)
    n = len(y)
    idx = np.arange(n)
    tree = build(X, types, y, mode, C, params, idx)
    alphas, sets, leaves = prune_sequence(tree)
    xerr, se, bt = cross_validate(X, types, y, mode, C, params, alphas, seed, min(nfold, n))
    k = 0 if rule == "none" else choose(xerr, se, rule)
    final_sets = sets[k]
    final = export_nodes(tree, final_sets)
    # 재대입 지표와 변수 중요도 (최종 나무)
    preds = np.array([tree.predict_leaf([X[c][i] for c in range(len(X))], final_sets)["pred"] for i in range(n)])
    res = dict()
    if mode == "class":
        conf = np.zeros((C, C), int)
        for a, b in zip(y.astype(int), preds.astype(int)):
            conf[a, b] += 1
        res["confusion"] = conf.tolist()
        res["accuracy"] = float(np.trace(conf) / n)
    else:
        sse = float(((y - preds) ** 2).sum())
        sst = float(((y - y.mean()) ** 2).sum())
        res["r2"] = 1 - sse / sst
        res["rmse"] = math.sqrt(sse / n)
    imp = {}
    for nd in final:
        if nd["split"] is not None:
            imp[var_names[nd["split"]["j"]]] = imp.get(var_names[nd["split"]["j"]], 0.0) + nd["split"]["gain"]
    return dict(
        name=name, mode=mode, C=C, params=params, rule=rule, seed=seed, nfold=nfold,
        full=[dict(id=nd["id"], depth=nd["depth"], n=nd["n"], pred=nd["pred"], err=nd["err"], left=nd["left"], right=nd["right"],
                   split=None if nd["split"] is None else dict(j=nd["split"]["j"], kind=nd["split"]["kind"], arg=nd["split"]["arg"] if nd["split"]["kind"] == "num" else [int(v) for v in nd["split"]["arg"]], gain=nd["split"]["gain"]))
              for nd in tree.nodes],
        alphas=alphas, leaves=leaves, xerror=[float(v) for v in xerr], se=[float(v) for v in se], selected=k,
        final=final, importance=imp, **res,
        train_err=[subtree_stats(tree, tree.nodes[0], s)[0] for s in sets],
    )


rng = np.random.default_rng(20261034)
cases = []
data = {}

# A: 이진 분류. 숫자 2개 + 범주 1개(3범주) + 잡음 변수
n = 200
x1 = np.round(rng.normal(0, 1, n), 3)
x2 = np.round(rng.uniform(0, 10, n), 2)
c = rng.integers(0, 3, n)
noise = np.round(rng.normal(0, 1, n), 3)
logit = -0.3 + 2.2 * (x1 > 0.3) - 1.8 * (x2 > 6) + 1.0 * (c == 2)
yA = (rng.random(n) < 1 / (1 + np.exp(-logit))).astype(int)
data["A"] = dict(x1=x1.tolist(), x2=x2.tolist(), noise=noise.tolist(), c=c.tolist(), y=yA.tolist())
P1 = dict(maxDepth=5, minSplit=20, minLeaf=7)
cases.append(run_case("이진 분류, 가지치기 1-SE", [x1, x2, noise, c], ["num", "num", "num", "cat"], yA, "class", 2, P1, "cv1se", 2026, 10, ["x1", "x2", "noise", "c"]))
cases.append(run_case("이진 분류, 가지치기 최소 CV 오차", [x1, x2, noise, c], ["num", "num", "num", "cat"], yA, "class", 2, P1, "cv", 7, 5, ["x1", "x2", "noise", "c"]))
cases.append(run_case("이진 분류, 가지치기 없음", [x1, x2, noise, c], ["num", "num", "num", "cat"], yA, "class", 2, dict(maxDepth=3, minSplit=20, minLeaf=7), "none", 1, 10, ["x1", "x2", "noise", "c"]))

# B: 3분류. 숫자 2개 + 범주 1개(4범주, 모든 부분집합 탐색)
n = 150
x1 = np.round(rng.normal(0, 1, n), 3)
x2 = np.round(rng.uniform(0, 10, n), 2)
c = rng.integers(0, 4, n)
yB = np.where(x1 > 0.5, 2, np.where((c == 1) | (c == 3), 1, 0))
flip = rng.random(n) < 0.12
yB = np.where(flip, rng.integers(0, 3, n), yB)
data["B"] = dict(x1=x1.tolist(), x2=x2.tolist(), c=c.tolist(), y=yB.tolist())
cases.append(run_case("3분류, 범주 4개(부분집합 탐색)", [x1, x2, c], ["num", "num", "cat"], yB, "class", 3, P1, "cv1se", 11, 10, ["x1", "x2", "c"]))

# C: 회귀
n = 150
x1 = np.round(rng.normal(0, 1, n), 3)
x2 = np.round(rng.uniform(0, 10, n), 2)
c = rng.integers(0, 3, n)
yC = np.round(3 * (x1 > 0) + 0.5 * x2 + 2 * (c == 1) + rng.normal(0, 1, n), 3)
data["C"] = dict(x1=x1.tolist(), x2=x2.tolist(), c=c.tolist(), y=yC.tolist())
cases.append(run_case("회귀, 가지치기 1-SE", [x1, x2, c], ["num", "num", "cat"], yC, "reg", 0, dict(maxDepth=5, minSplit=20, minLeaf=7), "cv1se", 5, 10, ["x1", "x2", "c"]))

path = os.path.join(os.path.dirname(__file__), "reference", "cart.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"data": data, "cases": cases}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "cart cases")
for cs in cases:
    print(" ", cs["name"], "| full nodes", len(cs["full"]), "| leaves seq", cs["leaves"], "| selected k", cs["selected"], "| final nodes", len(cs["final"]),
          "| acc/r2", round(cs.get("accuracy", cs.get("r2", 0)), 3))
