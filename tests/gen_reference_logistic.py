"""로지스틱 회귀 기준값을 tests/reference/logistic.json 에 저장한다.

실행: py -3.12 tests/gen_reference_logistic.py
- 계수와 표준오차: numpy 로 독립 구현한 Newton-Raphson (정보행렬의 역행렬). 수렴 뒤 기울기(score)가 0 인지 확인한다.
  (이 도구의 가중 QR 방식과 다른 계산 경로)
- 모형 적합: 로그가능도, 가능도비 카이제곱(scipy chi2), Cox-Snell, Nagelkerke, McFadden
- Hosmer-Lemeshow: 예측확률의 10분위(np.quantile 선형 보간) 경계로 (b0, b1], … 나눔, 같은 경계는 합침, df = 집단 수 − 2
- 분류표: 예측확률 > 절단값 이면 사건
- AUC: 사건·비사건 쌍을 모두 비교(동점 0.5)하는 전수 계산, DeLong 표준오차는 V10, V01 을 직접 계산
"""
import json
import os

import numpy as np
import scipy
from scipy import stats

rng = np.random.default_rng(20261030)


def newton(X, y):
    """X: (n, k) 설계 행렬(상수항 포함). 반환 B, cov, ll."""
    n, k = X.shape
    beta = np.zeros(k)
    beta[0] = np.log(y.sum() / (n - y.sum()))
    for _ in range(200):
        eta = X @ beta
        p = 1 / (1 + np.exp(-eta))
        grad = X.T @ (y - p)
        W = p * (1 - p)
        H = X.T @ (X * W[:, None])
        step = np.linalg.solve(H, grad)
        beta = beta + step
        if np.max(np.abs(step)) < 1e-12:
            break
    eta = X @ beta
    p = 1 / (1 + np.exp(-eta))
    grad = X.T @ (y - p)
    assert np.max(np.abs(grad)) < 1e-6 * max(1.0, n), grad
    W = p * (1 - p)
    cov = np.linalg.inv(X.T @ (X * W[:, None]))
    ll = float(np.sum(y * eta - np.logaddexp(0, eta)))
    return beta, cov, ll, p


def hl_test(p, y, g=10):
    sorted_p = np.sort(p)
    raw = [float(np.quantile(sorted_p, i / g)) for i in range(g + 1)]
    breaks = [raw[0]]
    for b in raw[1:]:
        if b > breaks[-1]:
            breaks.append(b)
    ng = len(breaks) - 1
    inner = np.array(breaks[1:-1])
    idx = np.searchsorted(inner, p, side="left")
    groups = []
    for gi in range(ng):
        m = idx == gi
        groups.append(dict(n=int(m.sum()), o1=float(y[m].sum()), e1=float(p[m].sum())))
    chi = None
    df = None
    pv = None
    if ng >= 3:
        chi = 0.0
        for gr in groups:
            e1, e0 = gr["e1"], gr["n"] - gr["e1"]
            o1, o0 = gr["o1"], gr["n"] - gr["o1"]
            if e1 > 1e-12:
                chi += (o1 - e1) ** 2 / e1
            if e0 > 1e-12:
                chi += (o0 - e0) ** 2 / e0
        df = ng - 2
        pv = float(stats.chi2.sf(chi, df))
    return dict(groups=groups, chi2=None if chi is None else float(chi), df=df, p=pv)


def roc(p, y, conf):
    pos = p[y == 1]
    neg = p[y == 0]
    m, q = len(pos), len(neg)
    v10 = np.array([np.mean([1.0 if a > b else (0.5 if a == b else 0.0) for b in neg]) for a in pos])
    v01 = np.array([np.mean([1.0 if a > b else (0.5 if a == b else 0.0) for a in pos]) for b in neg])
    auc = float(v10.mean())
    se = float(np.sqrt(v10.var(ddof=1) / m + v01.var(ddof=1) / q))
    z = stats.norm.ppf(1 - (1 - conf / 100) / 2)
    se0 = np.sqrt((m + q + 1) / (12 * m * q))
    pv = float(2 * stats.norm.sf(abs((auc - 0.5) / se0)))
    return dict(auc=auc, se=se, lo=float(max(0, auc - z * se)), hi=float(min(1, auc + z * se)), p=pv)


def case(y, covs=(), factor=None, ncat=None, ref=0, cutoff=0.5, conf=95, extra=None):
    n = len(y)
    cols = [np.ones(n)] + [np.asarray(c, float) for c in covs]
    if factor is not None:
        for c in range(ncat):
            if c != ref:
                cols.append((factor == c).astype(float))
    X = np.column_stack(cols)
    # 평균 중심화한 설계로 계산한 뒤 원래 척도로 되돌린다 (큰 평균에서 조건수를 낮추기 위함)
    Xc = X.copy()
    means = X.mean(axis=0)
    means[0] = 0
    Xc[:, 1:] = X[:, 1:] - means[1:]
    bc, covc, ll, p = newton(Xc, y.astype(float))
    A = np.eye(X.shape[1])
    A[0, 1:] = -means[1:]
    B = A @ bc
    cov = A @ covc @ A.T
    se = np.sqrt(np.diag(cov))
    z = stats.norm.ppf(1 - (1 - conf / 100) / 2)
    n1 = y.sum()
    ll0 = float(n1 * np.log(n1 / n) + (n - n1) * np.log((n - n1) / n))
    lr = 2 * (ll - ll0)
    k = X.shape[1] - 1
    cs = 1 - np.exp(2 * (ll0 - ll) / n)
    nk = cs / (1 - np.exp(2 * ll0 / n))
    wald = (B / se) ** 2
    pred = (p > cutoff).astype(int)
    tp = int(((y == 1) & (pred == 1)).sum())
    fn = int(((y == 1) & (pred == 0)).sum())
    fp = int(((y == 0) & (pred == 1)).sum())
    tn = int(((y == 0) & (pred == 0)).sum())
    return dict(
        input=dict(y=[int(v) for v in y], covs=[[float(v) for v in c] for c in covs],
                   factor=None if factor is None else [int(v) for v in factor], ncat=ncat, ref=ref, cutoff=cutoff, conf=conf, **(extra or {})),
        expected=dict(
            n=n, n1=int(n1), k=k,
            B=[float(v) for v in B], se=[float(v) for v in se], wald=[float(v) for v in wald],
            p=[float(v) for v in stats.chi2.sf(wald, 1)], expB=[float(v) for v in np.exp(B)],
            lo=[float(v) for v in np.exp(B - z * se)], hi=[float(v) for v in np.exp(B + z * se)],
            ll=ll, ll0=ll0, lr=float(lr), lrP=float(stats.chi2.sf(lr, k)), cs=float(cs), nk=float(nk), mf=float(1 - ll / ll0),
            hl=hl_test(p, y), cls=dict(tp=tp, fn=fn, fp=fp, tn=tn), roc=roc(p, y, conf),
            probs=[float(v) for v in p],
        ),
    )


def sim(n, coefs, intercept, xs):
    eta = intercept + sum(c * x for c, x in zip(coefs, xs))
    return (rng.random(n) < 1 / (1 + np.exp(-eta))).astype(int)


cases = {}

n = 100
x1 = rng.normal(0, 1, n)
x2 = rng.normal(50, 10, n)
y = sim(n, [1.0, 0.04], -2.0, [x1, x2 - 50])
cases["연속 독립변수 2개 n=100"] = case(y, covs=[x1, x2])

n = 120
x1 = rng.normal(0, 1, n)
f = rng.integers(0, 3, n)
y = sim(n, [0.8], -0.3, [x1]) | (f == 2) * (rng.random(n) < 0.3)
cases["공변량 1개 + 범주형 3범주 (기준 첫 범주)"] = case(y.astype(int), covs=[x1], factor=f, ncat=3, ref=0)

n = 200
x1 = rng.normal(0, 1, n)
x2 = rng.normal(0, 1, n)
y = sim(n, [1.2, -0.6], -2.4, [x1, x2])
cases["사건이 드문 경우 (약 10%), 절단값 0.2"] = case(y, covs=[x1, x2], cutoff=0.2)

n = 90
x1 = rng.normal(10000, 500, n)
x2 = rng.normal(0, 1, n)
y = sim(n, [0.004, 0.7], -40.0, [x1, x2])
cases["큰 평균의 독립변수 (조건수), 90% 신뢰구간"] = case(y, covs=[x1, x2], conf=90)

n = 150
x1 = rng.normal(0, 1, n)
f = rng.integers(0, 2, n)
y = sim(n, [3.0, 1.0], 0.0, [x1, f])
cases["강한 효과 (완전 분리는 아님), 범주형 2범주 기준 마지막"] = case(y, covs=[x1], factor=f, ncat=2, ref=1)

n = 60
x1 = rng.integers(1, 6, n).astype(float)
y = sim(n, [0.6], -1.8, [x1])
cases["정수 독립변수 (동점 예측확률 많음), n=60"] = case(y, covs=[x1])

# ---- 분석 창에서 종속변수를 바로 이분하는 경우 (시트에서 미리 0/1로 바꾸지 않은 자료) ----
n = 140
x0 = rng.normal(0, 1, n)
x1 = rng.normal(0, 1, n)
latent = 0.9 * x0 - 0.5 * x1 + rng.logistic(0, 1, n)
raw = 1 + np.digitize(latent, [-1.0, 0.3, 1.5])           # 1~4 등급
y = (raw >= 3).astype(int)
cases["종속변수가 4등급(1~4): 3·4등급을 사건으로"] = case(y, covs=[x0, x1], extra=dict(raw=[int(v) for v in raw], eventValues=["3", "4"]))

n = 130
x0 = rng.normal(0, 1, n)
x1 = rng.normal(0, 1, n)
score = np.round(70 + 8 * x0 + 8 * rng.normal(0, 1, n) + 1.5 * x1, 1)
y = (score >= 70).astype(int)
cases["연속 점수 70 이상을 사건으로"] = case(y, covs=[x0, x1], extra=dict(raw=[float(v) for v in score], eventAtLeast=70.0))

n = 120
x1 = rng.normal(0, 1, n)
f = rng.integers(0, 3, n)
y = sim(n, [0.8], -0.2, [x1]) | (f == 0) * (rng.random(n) < 0.35)
cases["범주형 3범주에서 가운데 범주를 기준으로 지정"] = case(y.astype(int), covs=[x1], factor=f, ncat=3, ref=1, extra=dict(refCategory="1"))

path = os.path.join(os.path.dirname(__file__), "reference", "logistic.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "cases": [dict(name=k, **v) for k, v in cases.items()]}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "logistic cases to", path)
