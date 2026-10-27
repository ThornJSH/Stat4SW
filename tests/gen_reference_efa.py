"""탐색적 요인분석 기준값을 tests/reference/efa.json 에 저장한다.

실행: py -3.12 tests/gen_reference_efa.py   (statsmodels, factor_analyzer 가 개발용으로 설치돼 있어야 한다)
기준값은 세 갈래로 얻고 서로 일치하는지 이 스크립트가 확인한다.
 1) numpy 로 정의대로 직접 계산 (이 도구의 JS 구현과 다른 경로)
    - 상관행렬 고유분해(PCA), 주축분해(SMC 초기 공통성, numpy eigh 반복, 수렴 1e-10),
      베리맥스는 R stats::varimax 와 같은 알고리즘을 numpy SVD 로 구현 (JS 는 극분해를 Jacobi 고유분해로 구함)
    - KMO: 상관행렬의 역행렬에서 편상관, Bartlett: −(n−1−(2p+5)/6) ln|R|
 2) statsmodels: PCA 고유값, Factor(method='pa') 주축분해 적재값·공통성, (statsmodels 의 varimax 는 Kaiser 정규화를 하지 않아 회전은 대조하지 않음)
 3) factor_analyzer: calculate_kmo, calculate_bartlett_sphericity, Rotator(varimax, Kaiser 정규화, promax 패턴·요인 상관). (FactorAnalyzer.fit 은 설치된 scikit-learn 과 맞지 않아 쓰지 않는다)
열은 제곱합 큰 순으로 정렬하고 각 열의 합이 양수가 되도록 부호를 맞춘 형태(canon)로 저장한다 (도구와 같은 규약).
저장하는 값: R 의 고유값, KMO와 변수별 MSA, Bartlett(χ², df, p), 요인 수, 회전 전 적재값, 초기·추출 공통성, 제곱합 적재값, 회전 뒤 적재값
"""
import json
import os
import warnings

import numpy as np
import scipy
import statsmodels
from factor_analyzer import calculate_bartlett_sphericity, calculate_kmo
from factor_analyzer.rotator import Rotator
from scipy import stats
from statsmodels.multivariate.factor import Factor
from statsmodels.multivariate.pca import PCA

warnings.filterwarnings("ignore")
rng = np.random.default_rng(20261027)


def canon(L):
    L = np.asarray(L, float)
    order = np.argsort(-(L ** 2).sum(0), kind="stable")
    L = L[:, order]
    sign = np.where(L.sum(0) < 0, -1.0, 1.0)
    return L * sign


def varimax(L, eps=1e-13, maxit=5000):
    """R 의 stats::varimax(normalize = TRUE) 와 같은 알고리즘, 수렴을 더 엄격하게."""
    sc = np.sqrt((L ** 2).sum(1))
    sc[sc == 0] = 1
    x = L / sc[:, None]
    p, nc = x.shape
    TT = np.eye(nc)
    d = 0.0
    for _ in range(maxit):
        z = x @ TT
        B = x.T @ (z ** 3 - z @ np.diag((z ** 2).sum(0)) / p)
        u, s, vt = np.linalg.svd(B)
        TT = u @ vt
        dpast, d = d, s.sum()
        if d < dpast * (1 + eps):
            break
    return (x @ TT) * sc[:, None]


def promax(L, kappa=4):
    """Kaiser 정규화 프로맥스 (factor_analyzer 와 같은 방식; R 의 stats::promax 는 정규화하지 않는다). numpy lstsq 로. 패턴 행렬과 요인 상관을 정준화해 돌려준다."""
    rl = np.sqrt((L ** 2).sum(1))
    x = varimax(L / rl[:, None])
    Q = x * np.abs(x) ** (kappa - 1)
    U = np.linalg.lstsq(x, Q, rcond=None)[0]
    d = np.diag(np.linalg.inv(U.T @ U))
    U = U * np.sqrt(d)
    P = (x @ U) * rl[:, None]
    phi = np.linalg.inv(U.T @ U)
    order = np.argsort(-(P ** 2).sum(0), kind="stable")
    P = P[:, order]
    sign = np.where(P.sum(0) < 0, -1.0, 1.0)
    phi = phi[np.ix_(order, order)] * np.outer(sign, sign)
    return P * sign, phi


def paf(R, m, tol=1e-10, maxit=5000):
    h2 = 1 - 1 / np.diag(np.linalg.inv(R))
    init = h2.copy()
    for it in range(1, maxit + 1):
        Rr = R.copy()
        np.fill_diagonal(Rr, h2)
        w, v = np.linalg.eigh(Rr)
        idx = np.argsort(-w)[:m]
        L = v[:, idx] * np.sqrt(np.maximum(w[idx], 0))
        nh = (L ** 2).sum(1)
        done = np.abs(nh - h2).max() < tol
        h2 = nh
        if done:
            break
    return L, h2, init, it


def kmo_bartlett(R, n):
    p = R.shape[0]
    Ri = np.linalg.inv(R)
    Q = -Ri / np.sqrt(np.outer(np.diag(Ri), np.diag(Ri)))
    off = ~np.eye(p, dtype=bool)
    r2 = np.where(off, R ** 2, 0)
    q2 = np.where(off, Q ** 2, 0)
    kmo = r2.sum() / (r2.sum() + q2.sum())
    msa = r2.sum(1) / (r2.sum(1) + q2.sum(1))
    sign, logdet = np.linalg.slogdet(R)
    chi2 = -(n - 1 - (2 * p + 5) / 6) * logdet
    df = p * (p - 1) / 2
    return float(kmo), msa.tolist(), dict(chi2=float(chi2), df=int(df), p=float(stats.chi2.sf(chi2, df)))


def analyse(name, X, method, rotation, nfactors=None):
    X = np.asarray(X, float)
    n, p = X.shape
    R = np.corrcoef(X, rowvar=False)
    w, v = np.linalg.eigh(R)
    idx = np.argsort(-w)
    eig, vec = w[idx], v[:, idx]
    kaiser = int((eig > 1).sum())
    m = nfactors if nfactors else kaiser
    mauto = nfactors is None
    kmo, msa, bart = kmo_bartlett(R, n)
    if method == "pca":
        L = vec[:, :m] * np.sqrt(eig[:m])
        h2_init = np.ones(p)
        iters = None
    else:
        L, h2, h2_init, iters = paf(R, m)
    L = canon(L)
    h2 = (L ** 2).sum(1)
    out = dict(name=name, n=n, p=p, method=method, rotation=rotation, nfactors=nfactors, m=m, mauto=bool(mauto), kaiser=kaiser,
               data=X.tolist(), eig=eig.tolist(), kmo=kmo, msa=msa, bartlett=bart, L=L.tolist(), h2=h2.tolist(), h2_init=h2_init.tolist(),
               ss_ext=(L ** 2).sum(0).tolist(), paf_iter=iters, rotated=None, ss_rot=None)
    if rotation == "varimax" and m > 1:
        Rot = canon(varimax(L))
        out["rotated"] = Rot.tolist()
        out["ss_rot"] = (Rot ** 2).sum(0).tolist()

    # ---- 프로맥스, 요인 점수 ----
    out["promax"] = None
    if m > 1:
        P, phi = promax(L)
        S = P @ phi
        out["promax"] = dict(P=P.tolist(), phi=phi.tolist(), S=S.tolist())
        # 공통성은 회전해도 같다: diag(P Φ P')
        np.testing.assert_allclose(np.einsum("ij,jk,ik->i", P, phi, P), h2, rtol=1e-8, err_msg=name + " 프로맥스 공통성")
        rot = Rotator(method="promax", max_iter=5000, tol=1e-12).fit_transform(L)
        np.testing.assert_allclose(canon(rot), P, rtol=1e-4, atol=1e-5, err_msg=name + " factor_analyzer promax 패턴")
        phi_fa = np.asarray(Rotator(method="promax", max_iter=5000, tol=1e-12).fit(L).phi_)
        o = np.argsort(-(rot ** 2).sum(0), kind="stable")
        sg = np.where(rot[:, o].sum(0) < 0, -1.0, 1.0)
        np.testing.assert_allclose(phi_fa[np.ix_(o, o)] * np.outer(sg, sg), phi, rtol=1e-4, atol=1e-5, err_msg=name + " factor_analyzer promax 요인 상관")
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    Rinv = np.linalg.inv(R)
    shown = np.array(out["rotated"]) if out["rotated"] is not None else L
    sc = Z @ (Rinv @ shown)
    out["scores"] = dict(first=sc[:6].tolist(), var=sc.var(0, ddof=1).tolist(), mean=sc.mean(0).tolist())
    if m > 1:
        sc2 = Z @ (Rinv @ np.array(out["promax"]["S"]))
        out["scores_promax"] = dict(first=sc2[:6].tolist(), var=sc2.var(0, ddof=1).tolist())

    # ---- 교차 확인 ----
    chk = lambda a, b, what, tol=1e-7, atol=1e-9: np.testing.assert_allclose(np.asarray(a, float), np.asarray(b, float), rtol=tol, atol=atol, err_msg=f"{name} {what}")
    fk, ft = calculate_kmo(X)
    chk(ft, kmo, "factor_analyzer KMO"); chk(fk, msa, "factor_analyzer MSA")
    bc, bp = calculate_bartlett_sphericity(X)
    chk(bc, bart["chi2"], "factor_analyzer Bartlett chi2"); chk(bp, bart["p"], "factor_analyzer Bartlett p", atol=1e-12)
    pc = PCA(X, standardize=True, demean=True, method="eig")
    chk(pc.eigenvals / pc.eigenvals.sum(), eig / p, "statsmodels PCA 고유값 비율")
    if method == "pca":
        # statsmodels PCA 의 고유벡터 방향 (단위 벡터)에 √고유값을 곱해 비교 (고유값 비율은 위에서 따로 확인)
        chk(canon(pc.loadings.values[:, :m] * np.sqrt(eig[:m])) if hasattr(pc.loadings, "values") else canon(pc.loadings[:, :m] * np.sqrt(eig[:m])), L, "statsmodels PCA 적재값", 1e-6)
    else:
        sm = Factor(X, n_factor=m, method="pa").fit(maxiter=5000, tol=1e-10)
        chk(canon(sm.loadings), L, "statsmodels pa 적재값", 1e-5, 1e-6)
        chk(sm.communality, h2, "statsmodels pa 공통성", 1e-5, 1e-6)
    if out["rotated"] is not None:
        # factor_analyzer 의 Rotator: Kaiser 정규화 베리맥스 (수렴 허용 1e-5 라 1e-4 로 비교; statsmodels 의 varimax 는 정규화를 하지 않아 대조하지 않는다)
        rot = Rotator(method="varimax", max_iter=5000, tol=1e-12).fit_transform(L)
        chk(canon(rot), out["rotated"], "factor_analyzer Rotator varimax", 1e-5, 1e-6)
    return out


def latent(n, loads, phi=0.0, noise=0.5, likert=False, seed=None):
    """loads: p × k 적재값 (0 이면 해당 요인과 무관). 요인 상관 phi. 반환 n × p (소수 3자리 또는 1~5 정수)"""
    loads = np.asarray(loads, float)
    k = loads.shape[1]
    C = np.full((k, k), phi) + np.eye(k) * (1 - phi)
    F = rng.multivariate_normal(np.zeros(k), C, n)
    E = rng.normal(0, 1, (n, loads.shape[0])) * noise
    X = F @ loads.T + E
    if likert:
        X = np.clip(np.round(2.5 + X * 0.9), 1, 5)
    return np.round(X, 3)


cases = []
# A: 직교 3요인, 9변수, 300명, PAF + 베리맥스 + 고유값 1 이상
L9 = np.zeros((9, 3))
for f in range(3):
    L9[3 * f:3 * f + 3, f] = [0.8, 0.7, 0.6]
cases.append(analyse("300명 9변수 직교 3요인, 주축분해·베리맥스", latent(300, L9, 0.0), "paf", "varimax"))
# B: 상관 2요인(r=.4) + 교차적재, 8변수, 200명, PAF 요인 수 2 지정
L8 = np.array([[0.8, 0], [0.7, 0], [0.7, 0.3], [0.6, 0], [0, 0.8], [0, 0.7], [0.3, 0.6], [0, 0.6]])
XB = latent(200, L8, 0.4)
cases.append(analyse("200명 8변수 상관 2요인(교차적재 포함), 주축분해·베리맥스, 요인 수 2", XB, "paf", "varimax", 2))
# C: 5점 척도 정수 자료 12변수 3요인, PCA + 베리맥스 + 고유값 1 이상
L12 = np.zeros((12, 3))
for f in range(3):
    L12[4 * f:4 * f + 4, f] = [0.8, 0.75, 0.7, 0.65]
cases.append(analyse("150명 12변수 5점 척도(정수), 주성분·베리맥스", latent(150, L12, 0.2, likert=True), "pca", "varimax"))
# D: 120명 6변수 2요인, PAF, 회전 없음
L6 = np.array([[0.8, 0], [0.75, 0], [0.7, 0.1], [0, 0.8], [0, 0.7], [0.1, 0.65]])
cases.append(analyse("120명 6변수 2요인, 주축분해·회전 없음", latent(120, L6, 0.0), "paf", "none", 2))
# E: 500명 20변수 4요인, PAF + 베리맥스 (큰 문제)
L20 = np.zeros((20, 4))
for f in range(4):
    L20[5 * f:5 * f + 5, f] = [0.75, 0.7, 0.65, 0.6, 0.55]
cases.append(analyse("500명 20변수 4요인, 주축분해·베리맥스", latent(500, L20, 0.1), "paf", "varimax"))
# F: 80명 6변수 1요인, PAF (요인 1개, 회전하지 않음)
cases.append(analyse("80명 6변수 1요인, 주축분해", latent(80, np.array([[0.8], [0.75], [0.7], [0.65], [0.6], [0.55]]), 0.0), "paf", "varimax"))
# G: B 와 같은 자료, PCA + 베리맥스 + 요인 수 2
cases.append(analyse("200명 8변수 상관 2요인, 주성분·베리맥스, 성분 수 2", XB, "pca", "varimax", 2))

def parallel(n, p, sims=4000, seed=7):
    """병렬분석 기준값: numpy 난수로 sims 번 모의한 상관행렬 고유값의 평균과 95백분위 (JS 와 난수가 달라 허용오차 비교용)"""
    g = np.random.default_rng(seed)
    E = np.empty((sims, p))
    for s in range(sims):
        E[s] = np.sort(np.linalg.eigvalsh(np.corrcoef(g.standard_normal((n, p)), rowvar=False)))[::-1]
    return dict(n=n, p=p, sims=sims, mean=E.mean(0).tolist(), p95=np.quantile(E, 0.95, axis=0).tolist())


par = [parallel(c["n"], c["p"]) for c in {(c["n"], c["p"]): c for c in cases}.values()]

path = os.path.join(os.path.dirname(__file__), "reference", "efa.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"statsmodels": statsmodels.__version__, "scipy": scipy.__version__, "cases": cases, "parallel": par}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "efa cases (numpy = statsmodels = factor_analyzer 일치 확인 통과);",
      [(c["n"], c["p"], c["m"], round(c["kmo"], 3), c["paf_iter"], round(max(c["h2"]), 3)) for c in cases])
