"""서열(비례오즈) 로지스틱 회귀 기준값을 tests/reference/ologit.json 에 저장한다.

실행: py -3.12 tests/gen_reference_ologit.py   (statsmodels 가 개발용으로 설치돼 있어야 한다)
기준값은 세 갈래로 얻고 서로 일치하는지 이 스크립트가 확인한다.
 1) statsmodels OrderedModel(distr='logit') 의 최대우도 추정: 계수, 문턱값(delta 방법으로 표준오차), 로그가능도, 절편만 있는 모형의 로그가능도, 예측 확률
 2) numpy 로 정의대로 직접 계산: 로그가능도 함수를 따로 구현해 추정값에서 평가하고, 중심차분 기울기가 0, 중심차분 헤세 행렬로 구한 표준오차가 statsmodels 와 같은지 확인
 3) Brant 검정: 문턱값마다 statsmodels Logit(Y > j)의 계수·공분산을 쓰고, 계수 사이의 공분산(Brant 1990 식)과 Wald 통계량은 numpy 로 계산
    — Brant 통계량은 비교할 도구(R brant 등)가 없어 같은 식을 numpy 로 따로 쓴 것이라 식 자체를 잘못 이해했다면 못 잡는다. 문턱값별 계수와 공분산은 statsmodels 로 독립 확인된다.
유사 R²(Cox & Snell, Nagelkerke, McFadden)와 분류표는 정의대로 numpy 로 계산한다.
"""
import json
import os
import warnings

import numpy as np
import scipy
import statsmodels
import statsmodels.api as sm
from scipy import stats
from statsmodels.miscmodels.ordinal_model import OrderedModel

warnings.filterwarnings("ignore")
rng = np.random.default_rng(20261102)
chk = lambda x, y, what, tol=1e-6, atol=1e-8: np.testing.assert_allclose(np.asarray(x, float), np.asarray(y, float), rtol=tol, atol=atol, err_msg=what)


def loglik(phi, X, y, J):
    k = X.shape[1]
    beta, theta = phi[:k], phi[k:]
    eta = X @ beta
    th = np.concatenate([[-np.inf], theta, [np.inf]])
    F = 1 / (1 + np.exp(-(th[y + 1] - eta)))
    G = 1 / (1 + np.exp(-(th[y] - eta)))
    F = np.where(np.isinf(th[y + 1]), 1.0, F)
    G = np.where(np.isinf(th[y]), 0.0, G)
    return float(np.log(F - G).sum())


def num_hess(f, x, h=1e-4):
    p = len(x)
    H = np.zeros((p, p))
    for i in range(p):
        for j in range(i, p):
            xs = [x.copy() for _ in range(4)]
            xs[0][i] += h; xs[0][j] += h
            xs[1][i] += h; xs[1][j] -= h
            xs[2][i] -= h; xs[2][j] += h
            xs[3][i] -= h; xs[3][j] -= h
            H[i, j] = H[j, i] = (f(xs[0]) - f(xs[1]) - f(xs[2]) + f(xs[3])) / (4 * h * h)
    return H


def analyse(name, y_codes, cols, names, cats, conf=95):
    y_codes = np.asarray(y_codes)
    J = len(cats)
    y = np.array([cats.index(v) for v in y_codes])
    X = np.column_stack(cols)
    n, k = X.shape
    df = {nm: X[:, i] for i, nm in enumerate(names)}
    mod = OrderedModel(y, X, distr="logit")
    res = mod.fit(method="bfgs", disp=False, maxiter=2000, gtol=1e-10)
    beta = np.asarray(res.params[:k])
    raw = np.asarray(res.params[k:])
    theta = np.concatenate([[raw[0]], raw[0] + np.cumsum(np.exp(raw[1:]))])
    # delta 방법: θ = T(raw)
    Jm = np.zeros((J - 1, len(res.params)))
    Jm[0, k] = 1
    for j in range(1, J - 1):
        Jm[j, k] = 1
        for q in range(1, j + 1):
            Jm[j, k + q] = np.exp(raw[q])
    cov_p = np.asarray(res.cov_params())
    T = np.zeros((len(res.params), len(res.params)))
    T[:k, :k] = np.eye(k)
    T[k:, k:] = Jm[:, k:]
    cov = T @ cov_p @ T.T
    se = np.sqrt(np.diag(cov))
    phi = np.concatenate([beta, theta])
    ll = float(res.llf)
    # numpy 로 같은 로그가능도, 기울기 0, 헤세
    chk(loglik(phi, X, y, J), ll, name + " numpy 로그가능도", 1e-9)
    grad = np.array([(loglik(phi + e, X, y, J) - loglik(phi - e, X, y, J)) / 2e-5 for e in np.eye(len(phi)) * 1e-5])
    assert np.abs(grad).max() < 1e-3, name + " 기울기 " + str(np.abs(grad).max())
    H = num_hess(lambda v: loglik(v, X, y, J), phi)
    cov_num = np.linalg.inv(-H)
    chk(np.sqrt(np.diag(cov_num)), se, name + " 표준오차 (수치 헤세)", 2e-4, 1e-6)
    counts = np.bincount(y, minlength=J)
    ll0 = float((counts * np.log(counts / n)).sum())
    chk(ll0, res.llnull, name + " 절편만 모형 로그가능도", 1e-9)
    chi2 = 2 * (ll - ll0)
    cs = 1 - np.exp(2 * (ll0 - ll) / n)
    nk = cs / (1 - np.exp(2 * ll0 / n))
    probs = np.asarray(res.predict())
    # numpy 확률
    eta = X @ beta
    th = np.concatenate([[-np.inf], theta, [np.inf]])
    P = np.array([[1 / (1 + np.exp(-(th[j + 1] - e))) - 1 / (1 + np.exp(-(th[j] - e))) for j in range(J)] for e in eta])
    chk(P, probs, name + " 예측 확률", 1e-6, 1e-8)
    pred = P.argmax(1)
    table = np.zeros((J, J), int)
    for a, b in zip(y, pred):
        table[a, b] += 1
    z = stats.norm.ppf(1 - (1 - conf / 100) / 2)
    out = dict(name=name, n=int(n), J=J, k=k, names=names, cats=cats, y=[int(v) for v in y_codes], cols=[c.tolist() for c in cols], conf=conf,
               beta=beta.tolist(), theta=theta.tolist(), se=se.tolist(), cov=cov.tolist(), ll=ll, ll0=ll0, chi2=float(chi2), df=k,
               p=float(stats.chi2.sf(chi2, k)), cs=float(cs), nk=float(nk), mf=float(1 - ll / ll0), probs=P[:6].tolist(), counts=counts.tolist(),
               table=table.tolist(), acc=float(np.trace(table) / n))

    # ---- Brant ----
    S = J - 1
    fits, pis, covs = [], [], []
    Xc = sm.add_constant(X)
    for m in range(S):
        yb = (y > m).astype(int)
        r = sm.Logit(yb, Xc).fit(disp=False, maxiter=200)
        fits.append(np.asarray(r.params)); covs.append(np.asarray(r.cov_params())); pis.append(np.asarray(r.predict(Xc)))
    q = k + 1
    V = {}
    for m in range(S):
        for l in range(m, S):
            w = pis[l] - pis[m] * pis[l]
            Sml = (Xc * w[:, None]).T @ Xc
            V[(m, l)] = covs[m] @ Sml @ covs[l]
            V[(l, m)] = V[(m, l)].T
    big = np.zeros((S * k, S * k))
    for m in range(S):
        for l in range(S):
            big[m * k:(m + 1) * k, l * k:(l + 1) * k] = V[(m, l)][1:, 1:]
    bvec = np.concatenate([fits[m][1:] for m in range(S)])

    def wald(vars_):
        rows = []
        for v in vars_:
            for m in range(S - 1):
                r = np.zeros(S * k)
                r[m * k + v] = 1; r[(m + 1) * k + v] = -1
                rows.append(r)
        D = np.array(rows)
        diff = D @ bvec
        Mx = D @ big @ D.T
        c2 = float(diff @ np.linalg.solve(Mx, diff))
        return dict(chi2=c2, df=len(rows), p=float(stats.chi2.sf(c2, len(rows))))
    out["brant"] = dict(splitB=[f[1:].tolist() for f in fits], splitSe=[np.sqrt(np.diag(c))[1:].tolist() for c in covs],
                        omnibus=wald(list(range(k))), perVar=[wald([v]) for v in range(k)])
    return out


def sim_ologit(n, J, betas, cuts, xs, dist="logistic"):
    eta = sum(b * x for b, x in zip(betas, xs))
    lat = eta + rng.logistic(size=n)
    return np.digitize(lat, cuts) + 1


cases = []
# O1: J=4, 숫자 1개 + 이진 1개 (비례오즈 성립)
n = 200
x1 = np.round(rng.normal(0, 1, n), 2); x2 = rng.integers(0, 2, n).astype(float)
y = sim_ologit(n, 4, [0.8, 0.6], [-1.0, 0.0, 1.3], [x1, x2])
cases.append(analyse("4범주, 숫자·이진 독립변수 (비례오즈 성립)", y, [x1, x2], ["x1", "x2"], [1, 2, 3, 4]))
# O2: 5점 척도, 숫자 2개 + 3범주 범주형(더미)
n = 300
a = np.round(rng.normal(50, 10, n), 1); b = rng.integers(1, 6, n).astype(float); g = rng.integers(1, 4, n)
y = sim_ologit(n, 5, [0.05, 0.3, 0.5, 1.0], [-1.0 + 2.5 * 0, 2.5, 4.0, 5.2][:4], [np.round(a - 50, 1) * 1.0, b, (g == 2) * 1.0, (g == 3) * 1.0])
y = np.clip(y, 1, 5)
cats5 = sorted(set(y.tolist()))
cases.append(dict(analyse("5점 척도, 숫자 2개 + 3범주 요인(더미)", y, [a, b, (g == 2) * 1.0, (g == 3) * 1.0], ["a", "b", "g = 2", "g = 3"], cats5), g=g.tolist()))
# O3: 비례오즈 위반 (효과가 문턱값마다 다름)
n = 400
x = np.round(rng.normal(0, 1, n), 2); x2 = np.round(rng.normal(0, 1, n), 2)
u = rng.uniform(size=n)
lat = 0.6 * x2 + rng.logistic(size=n) * 0.9 + np.where(x > 0, 1.6 * x, 0.1 * x)
y = np.digitize(lat, [-0.8, 0.9]) + 1
cases.append(analyse("3범주, 비례오즈 위반 (한 변수의 효과가 비선형)", y, [x, x2], ["x", "x2"], sorted(set(y.tolist()))))
# O4: 7범주, 불균형
n = 220
x = np.round(rng.normal(0, 1, n), 2)
lat = 1.1 * x + rng.logistic(size=n)
y = np.digitize(lat, [-2.2, -1.2, -0.3, 0.6, 1.8, 3.0]) + 1
cases.append(analyse("7범주, 범주 분포 불균형, 독립변수 1개", y, [x], ["x"], sorted(set(y.tolist()))))
# O5: 정수 점수 독립변수 3개 (5점 척도 설문), 90% 신뢰구간
n = 250
s1 = rng.integers(1, 6, n).astype(float); s2 = rng.integers(1, 6, n).astype(float); s3 = rng.integers(1, 8, n).astype(float)
lat = 0.5 * s1 + 0.3 * s2 - 0.2 * s3 + rng.logistic(size=n)
y = np.digitize(lat, [-0.2, 0.9, 1.7, 2.6]) + 1
cases.append(analyse("5범주, 정수 점수 독립변수 3개, 90% 신뢰구간", y, [s1, s2, s3], ["s1", "s2", "s3"], sorted(set(y.tolist())), conf=90))

path = os.path.join(os.path.dirname(__file__), "reference", "ologit.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"statsmodels": statsmodels.__version__, "scipy": scipy.__version__, "cases": cases}, fh, ensure_ascii=False, indent=0, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("wrote", len(cases), "ologit cases (statsmodels OrderedModel = numpy 로그가능도·헤세 일치 확인 통과);",
      [(c["n"], c["J"], round(c["nk"], 3), round(c["brant"]["omnibus"]["p"], 4)) for c in cases])
