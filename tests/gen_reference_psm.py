"""성향점수 매칭 기준값을 tests/reference/psm.json 에 저장한다.

실행: py -3.12 tests/gen_reference_psm.py
- 성향점수: numpy 로 독립 구현한 Newton-Raphson 로지스틱 회귀 (수렴 뒤 기울기 0 확인), logit = X·B
- 매칭: 순수 파이썬 반복문으로 직접 구현한 탐욕적 최근접 이웃 (남은 후보를 모두 비교)
  · 처치 사례를 logit 성향점수 내림차순(동점은 데이터 순서)으로 한 번씩 처리, 1:k 는 k 라운드(첫 라운드에서 못 찾은 처치는 탈락)
  · 거리: |logit 차이| (ps) 또는 공변량 마할라노비스 거리(maha, 공분산은 전체 사례 n−1, 더미 포함). 캘리퍼 = c × SD(logit, n−1) 안의 후보만.
  · 동점은 데이터 순서가 앞선 사례. replace=True 면 대조 사례를 다시 쓸 수 있다 (같은 처치에 같은 대조는 한 번만).
  · 가중치: 짝지은 처치 1, 대조는 속한 세트마다 1/(세트의 대조 수)를 더함. 매칭 ID: 짝지은 처치를 데이터 순서대로 1, 2, …
- 균형: SMD = (처치 평균 − 대조 평균) / 처치군 SD(매칭 전, n−1). 매칭 후는 짝지은 처치 평균과 가중 대조 평균.
- 효과: 세트별 차이 d = y처치 − mean(y대조)의 평균 (ATT), SE = sd(d, n−1)/√m, t(m−1) 신뢰구간(scipy t). 매칭 전은 Welch t(scipy).
"""
import json
import os

import numpy as np
import scipy
from scipy import stats

rng = np.random.default_rng(20261032)


def newton(X, y):
    n, k = X.shape
    beta = np.zeros(k)
    beta[0] = np.log(y.sum() / (n - y.sum()))
    for _ in range(200):
        eta = X @ beta
        p = 1 / (1 + np.exp(-eta))
        grad = X.T @ (y - p)
        H = X.T @ (X * (p * (1 - p))[:, None])
        step = np.linalg.solve(H, grad)
        beta = beta + step
        if np.max(np.abs(step)) < 1e-12:
            break
    eta = X @ beta
    p = 1 / (1 + np.exp(-eta))
    assert np.max(np.abs(X.T @ (y - p))) < 1e-6 * n
    return beta


def match(eta, treat, covmat, caliper_c, ratio, replace, method):
    n = len(eta)
    sd = float(np.std(eta, ddof=1))
    width = caliper_c * sd if caliper_c > 0 else None
    VI = np.linalg.inv(np.cov(covmat, rowvar=False)) if method == "maha" else None
    treated = [i for i in range(n) if treat[i] == 1]
    controls = [i for i in range(n) if treat[i] == 0]
    order = sorted(treated, key=lambda i: (-eta[i], i))
    sets = {t: [] for t in treated}
    used = set()
    active = list(order)
    for rnd in range(ratio):
        nxt = []
        for t in active:
            best, bd = None, None
            for c in controls:
                if c in sets[t]:
                    continue
                if not replace and c in used:
                    continue
                if width is not None and abs(eta[t] - eta[c]) > width:
                    continue
                if method == "ps":
                    dd = abs(eta[t] - eta[c])
                else:
                    diff = covmat[t] - covmat[c]
                    dd = float(np.sqrt(diff @ VI @ diff))
                if bd is None or dd < bd:
                    best, bd = c, dd
            if best is None:
                continue
            sets[t].append(best)
            used.add(best)
            nxt.append(t)
        if rnd == 0:
            active = nxt          # 첫 라운드에서 못 찾은 처치는 이후 라운드에서도 제외
        else:
            active = nxt
    matched = sorted(t for t in treated if sets[t])
    w = np.zeros(n)
    for t in matched:
        w[t] += 1.0
        for c in sets[t]:
            w[c] += 1.0 / len(sets[t])
    ids = {t: k + 1 for k, t in enumerate(matched)}
    return sets, matched, w, ids, width


def smd_rows(cols, names, treat, w, matched):
    out = []
    ti = np.where(treat == 1)[0]
    ci = np.where(treat == 0)[0]
    for name, x in zip(names, cols):
        sd_t = float(np.std(x[ti], ddof=1))
        mt0, mc0 = float(x[ti].mean()), float(x[ci].mean())
        mt1 = float(x[matched].mean())
        wc = w[ci]
        mc1 = float((x[ci] * wc).sum() / wc.sum())
        out.append(dict(name=name, mt0=mt0, mc0=mc0, smd0=(mt0 - mc0) / sd_t, mt1=mt1, mc1=mc1, smd1=(mt1 - mc1) / sd_t))
    return out


def effect(y, treat, sets, matched):
    d = np.array([y[t] - np.mean([y[c] for c in sets[t]]) for t in matched])
    m = len(d)
    att = float(d.mean())
    se = float(d.std(ddof=1) / np.sqrt(m))
    tcrit = float(stats.t.ppf(0.975, m - 1))
    # 매칭 전 Welch
    a, b = y[treat == 1], y[treat == 0]
    wt = stats.ttest_ind(a, b, equal_var=False)
    return dict(att=att, se=se, t=att / se, df=m - 1, p=float(2 * stats.t.sf(abs(att / se), m - 1)), lo=att - tcrit * se, hi=att + tcrit * se,
                m_treat=float(np.mean([y[t] for t in matched])),
                m_ctrl=float(np.mean([np.mean([y[c] for c in sets[t]]) for t in matched])),
                naive_diff=float(a.mean() - b.mean()), naive_t=float(wt.statistic), naive_df=float(wt.df), naive_p=float(wt.pvalue))


# ---- 자료: 150사례, 연속 2개 + 범주 1개(3범주), 처치는 공변량에 따라 선택, 결과에는 처치효과 2 ----
n = 150
x1 = np.round(rng.normal(0, 1, n), 3)
x2 = np.round(rng.normal(50, 10, n), 2)
g = rng.integers(0, 3, n)
eta_true = -0.2 + 0.9 * x1 + 0.04 * (x2 - 50) + np.array([0.0, 0.5, -0.4])[g]
treat = (rng.random(n) < 1 / (1 + np.exp(-eta_true))).astype(int)
y = np.round(2.0 * treat + 1.2 * x1 + 0.05 * (x2 - 50) + 0.5 * (g == 1) + rng.normal(0, 1, n), 3)

G1 = (g == 1).astype(float)
G2 = (g == 2).astype(float)
X = np.column_stack([np.ones(n), x1, x2, G1, G2])
# 평균 중심화로 계산한 뒤 되돌린다 (이 도구와 같은 수치 안정화, 결과는 같다)
Xc = X.copy()
means = X.mean(0)
means[0] = 0
Xc[:, 1:] = X[:, 1:] - means[1:]
bc = newton(Xc, treat.astype(float))
A = np.eye(5)
A[0, 1:] = -means[1:]
B = A @ bc
eta = X @ B
cov_cols = [x1, x2, G1, G2]
cov_names = ["x1", "x2", "g = 1", "g = 2"]
covmat = np.column_stack(cov_cols)

cases = []


def add(name, caliper, ratio, replace, method):
    sets, matched, w, ids, width = match(eta, treat, covmat, caliper, ratio, replace, method)
    rows = smd_rows(cov_cols + [eta], cov_names + ["logit"], treat, w, matched)
    cases.append(dict(
        name=name, options=dict(caliper=caliper, ratio=ratio, replace=replace, method=method),
        width=width, sets={str(t): cs for t, cs in sets.items() if cs}, matched=matched, weights=[float(v) for v in w],
        ids={str(t): i for t, i in ids.items()}, smd=rows, effect=effect(y, treat, sets, matched),
        n_treated=int(treat.sum()), n_control=int(n - treat.sum()), n_used_controls=int(len({c for t in matched for c in sets[t]})),
    ))


add("성향점수 최근접 1:1, 캘리퍼 0.2, 비복원", 0.2, 1, False, "ps")
add("캘리퍼 안 마할라노비스 1:1, 캘리퍼 0.2, 비복원", 0.2, 1, False, "maha")
add("성향점수 1:2, 캘리퍼 없음, 비복원", 0.0, 2, False, "ps")
add("성향점수 1:1, 캘리퍼 0.1, 복원", 0.1, 1, True, "ps")

data = dict(x1=[float(v) for v in x1], x2=[float(v) for v in x2], g=[int(v) for v in g], treat=[int(v) for v in treat], y=[float(v) for v in y])
fit = dict(B=[float(v) for v in B], eta=[float(v) for v in eta], p=[float(v) for v in 1 / (1 + np.exp(-eta))], sd_eta=float(np.std(eta, ddof=1)),
           eta_min_t=float(eta[treat == 1].min()), eta_max_t=float(eta[treat == 1].max()), eta_min_c=float(eta[treat == 0].min()), eta_max_c=float(eta[treat == 0].max()))
path = os.path.join(os.path.dirname(__file__), "reference", "psm.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "data": data, "fit": fit, "cases": cases}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "psm cases;", "matched:", [len(c["matched"]) for c in cases], "treated", int(treat.sum()), "|SMD after| max:", [round(max(abs(r["smd1"]) for r in c["smd"]), 3) for c in cases])
