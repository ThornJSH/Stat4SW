"""평가자 일치도(Cohen·Fleiss kappa)와 급내상관(ICC) 기준값을 tests/reference/agreement.json 에 저장한다.

실행: py -3.12 tests/gen_reference_agreement.py   (statsmodels, pingouin, scikit-learn 이 개발용으로 설치돼 있어야 한다)
기준값은 세 갈래로 얻고 서로 일치하는지 이 스크립트가 확인한다.
 1) numpy 로 정의대로 직접 계산 (kappa: 교차표, 평균제곱 ICC, Fleiss: 사례별 일치도)
 2) scikit-learn cohen_kappa_score (가중치 없음·선형·이차), statsmodels cohens_kappa (kappa, 표준오차(비영가설·영가설), 신뢰구간, z, 유의확률, kappa 최대값),
    statsmodels fleiss_kappa
 3) pingouin intraclass_corr (ICC 6유형의 ICC, F, 자유도, 유의확률, 95% 신뢰구간)
Fleiss kappa 의 표준오차·z(영가설)는 Fleiss(1971) 식이며 비교할 도구가 없어 numpy 식만 저장한다 (범주별 kappa 와 전체 kappa 의 가중평균 항등식으로 내부 확인).
"""
import json
import os
import warnings

import numpy as np
import pandas as pd
import pingouin as pg
import scipy
import statsmodels
from scipy import stats
from sklearn.metrics import cohen_kappa_score
from statsmodels.stats.inter_rater import aggregate_raters, cohens_kappa, fleiss_kappa

warnings.filterwarnings("ignore")
pg.options["round"] = None      # pingouin 은 기본으로 결과를 반올림한다 (신뢰구간 열은 2자리)
pg.options["round.column.CI95"] = None
rng = np.random.default_rng(20261101)
chk = lambda x, y, what, tol=1e-9, atol=1e-12: np.testing.assert_allclose(np.asarray(x, float), np.asarray(y, float), rtol=tol, atol=atol, err_msg=what)


def table_of(a, b, cats):
    t = np.zeros((len(cats), len(cats)))
    for x, y in zip(a, b):
        t[cats.index(x), cats.index(y)] += 1
    return t


def cohen(a, b, cats, conf=95):
    """2명 평가자의 Cohen kappa(가중치 없음, 선형, 이차)를 정의대로 계산하고 sklearn·statsmodels 와 대조."""
    r = len(cats)
    T = table_of(a, b, cats)
    N = T.sum()
    P = T / N
    row, col = P.sum(1), P.sum(0)
    out = {}
    z = stats.norm.ppf(1 - (1 - conf / 100) / 2)
    for name, wt in (("none", None), ("linear", "linear"), ("quadratic", "quadratic")):
        if wt is None:
            d = 1 - np.eye(r)
        elif wt == "linear":
            d = np.abs(np.subtract.outer(np.arange(r), np.arange(r))) / (r - 1)
        else:
            d = (np.subtract.outer(np.arange(r), np.arange(r)) / (r - 1)) ** 2
        obs_d = (d * P).sum()
        exp_d = (d * np.outer(row, col)).sum()
        kappa = 1 - obs_d / exp_d
        # 표준오차 (Fleiss, Cohen, Everitt 1969; 가중 확장)
        w = 1 - d
        w_row = (col[None, :] * w).sum(1)       # i 행의 가중 기대
        w_col = (row[:, None] * w).sum(0)
        agree_wexp = (w * np.outer(row, col)).sum()
        fac = 1.0 / ((1 - agree_wexp) ** 2 * N)
        if wt is None:
            ta = (np.diag(P) * (1 - (row + col) * (1 - kappa)) ** 2).sum()
            tb = P * (col[:, None] + row) ** 2
            np.fill_diagonal(tb, 0)
            tb = (1 - kappa) ** 2 * tb.sum()
            tc = (kappa - agree_wexp * (1 - kappa)) ** 2
            var = (ta + tb - tc) * fac
            tc0 = (row * col * (row + col)).sum()
            var0 = (agree_wexp + agree_wexp ** 2 - tc0) * fac
        else:
            ta = (P * (w - (w_col[None, :] + w_row[:, None]) * (1 - kappa)) ** 2).sum()
            var = (ta - (kappa - agree_wexp * (1 - kappa)) ** 2) * fac
            var0 = ((np.outer(row, col) * (w - (w_col[None, :] + w_row[:, None])) ** 2).sum() - agree_wexp ** 2) * fac
        se, se0 = np.sqrt(var), np.sqrt(var0)
        z0 = kappa / se0
        po = 1 - obs_d
        pe = 1 - exp_d
        kmax = (np.minimum(row, col).sum() - np.diag(np.outer(row, col)).sum()) / (1 - np.diag(np.outer(row, col)).sum())
        res = dict(n=int(N), po=float(po), pe=float(pe), kappa=float(kappa), se=float(se), se0=float(se0), z=float(z0),
                   p=float(2 * stats.norm.sf(abs(z0))), lo=float(kappa - z * se), hi=float(kappa + z * se), kmax=float(kmax))
        # 교차 확인
        sk = cohen_kappa_score(a, b, labels=cats, weights=wt)
        chk(sk, kappa, f"sklearn kappa {name}", 1e-10)
        sm = cohens_kappa(T, wt=wt if wt else None)
        chk(sm.kappa, kappa, f"statsmodels kappa {name}", 1e-10)
        chk(sm.std_kappa, se, f"statsmodels 표준오차 {name}", 1e-9)
        chk(sm.std_kappa0, se0, f"statsmodels 영가설 표준오차 {name}", 1e-9)
        chk(sm.z_value, z0, f"statsmodels z {name}", 1e-8)
        chk(sm.pvalue_two_sided, res["p"], f"statsmodels p {name}", 1e-7, 1e-14)
        if wt is None:
            chk(sm.kappa_max, kmax, "statsmodels kappa_max", 1e-10)
        out[name] = res
    return dict(table=T.tolist(), cats=cats, **out)


def fleiss(R, cats, conf=95):
    """R: 사례 × 평가자 (범주 값). Fleiss kappa, 범주별 kappa, 표준오차(영가설), z, p, 평가자 쌍별 Cohen kappa"""
    R = np.asarray(R)
    N, n = R.shape
    X = np.zeros((N, len(cats)))
    for i in range(N):
        for v in R[i]:
            X[i, cats.index(v)] += 1
    p = X.sum(0) / (N * n)
    q = 1 - p
    Pi = ((X ** 2).sum(1) - n) / (n * (n - 1))
    Pbar, Pe = Pi.mean(), (p ** 2).sum()
    kappa = (Pbar - Pe) / (1 - Pe)
    chk(kappa, fleiss_kappa(aggregate_raters(R, n_cat=len(cats))[0] if False else X), "statsmodels fleiss_kappa", 1e-10)
    kj = 1 - (X * (n - X)).sum(0) / (N * n * (n - 1) * p * q)
    chk((p * q * kj).sum() / (p * q).sum(), kappa, "Fleiss 항등식(범주별 가중평균)", 1e-9)
    sp = (p * q).sum()
    se_j = np.sqrt(2.0 / (N * n * (n - 1)))
    se = np.sqrt(2.0 / (N * n * (n - 1))) * np.sqrt(sp ** 2 - (p * q * (q - p)).sum()) / sp
    z = kappa / se
    pair = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            pair[i, j] = 1.0 if i == j else cohen_kappa_score(R[:, i], R[:, j], labels=cats)
    iu = np.triu_indices(n, 1)
    return dict(n=int(N), raters=int(n), cats=cats, counts=X.sum(0).tolist(), p=p.tolist(), pbar=float(Pbar), pe=float(Pe), kappa=float(kappa),
                se=float(se), z=float(z), pval=float(2 * stats.norm.sf(abs(z))), kappa_j=kj.tolist(), se_j=float(se_j),
                z_j=(kj / se_j).tolist(), p_j=(2 * stats.norm.sf(np.abs(kj / se_j))).tolist(), pair=pair.tolist(), light=float(pair[iu].mean()))


def icc(Y, conf=95):
    """Y: 대상 × 평가자. ICC 6유형 (McGraw & Wong), 신뢰구간 식은 pingouin 과 같다 (신뢰수준만 일반화)."""
    Y = np.asarray(Y, float)
    n, k = Y.shape
    gm = Y.mean()
    ssb = k * ((Y.mean(1) - gm) ** 2).sum()
    ssj = n * ((Y.mean(0) - gm) ** 2).sum()
    sst = ((Y - gm) ** 2).sum()
    sse = sst - ssb - ssj
    dfb, dfj, dfe = n - 1, k - 1, (n - 1) * (k - 1)
    msb, msj, mse = ssb / dfb, ssj / dfj, sse / dfe
    msw = (ssj + sse) / (dfj + dfe)
    alpha = 1 - conf / 100
    icc1 = (msb - msw) / (msb + (k - 1) * msw)
    icc2 = (msb - mse) / (msb + (k - 1) * mse + k * (msj - mse) / n)
    icc3 = (msb - mse) / (msb + (k - 1) * mse)
    icc1k = (msb - msw) / msb
    icc2k = (msb - mse) / (msb + (msj - mse) / n)
    icc3k = (msb - mse) / msb
    F1, F2 = msb / msw, msb / mse
    df1kd, df2kd = n * (k - 1), (n - 1) * (k - 1)
    ppf = lambda q, a, b: stats.f.ppf(q, a, b)
    f1l = F1 / ppf(1 - alpha / 2, dfb, df1kd); f1u = F1 * ppf(1 - alpha / 2, df1kd, dfb)
    l1, u1 = (f1l - 1) / (f1l + (k - 1)), (f1u - 1) / (f1u + (k - 1))
    f3l = F2 / ppf(1 - alpha / 2, dfb, df2kd); f3u = F2 * ppf(1 - alpha / 2, df2kd, dfb)
    l3, u3 = (f3l - 1) / (f3l + (k - 1)), (f3u - 1) / (f3u + (k - 1))
    fj = msj / mse
    vn = df2kd * (k * icc2 * fj + n * (1 + (k - 1) * icc2) - k * icc2) ** 2
    vd = dfb * k ** 2 * icc2 ** 2 * fj ** 2 + (n * (1 + (k - 1) * icc2) - k * icc2) ** 2
    v = vn / vd
    f2u, f2l = ppf(1 - alpha / 2, n - 1, v), ppf(1 - alpha / 2, v, n - 1)
    l2 = n * (msb - f2u * mse) / (f2u * (k * msj + (k * n - k - n) * mse) + n * msb)
    u2 = n * (f2l * msb - mse) / (k * msj + (k * n - k - n) * mse + n * f2l * msb)
    rows = [
        dict(key="icc1", icc=icc1, F=F1, df1=dfb, df2=df1kd, lo=l1, hi=u1),
        dict(key="icc2", icc=icc2, F=F2, df1=dfb, df2=df2kd, lo=l2, hi=u2),
        dict(key="icc3", icc=icc3, F=F2, df1=dfb, df2=df2kd, lo=l3, hi=u3),
        dict(key="icc1k", icc=icc1k, F=F1, df1=dfb, df2=df1kd, lo=1 - 1 / f1l, hi=1 - 1 / f1u),
        dict(key="icc2k", icc=icc2k, F=F2, df1=dfb, df2=df2kd, lo=l2 * k / (1 + l2 * (k - 1)), hi=u2 * k / (1 + u2 * (k - 1))),
        dict(key="icc3k", icc=icc3k, F=F2, df1=dfb, df2=df2kd, lo=1 - 1 / f3l, hi=1 - 1 / f3u),
    ]
    for r in rows:
        r["p"] = float(stats.f.sf(r["F"], r["df1"], r["df2"]))
    if abs(conf - 95) < 1e-12:
        long = pd.DataFrame([dict(s=i, r=j, y=Y[i, j]) for i in range(n) for j in range(k)])
        pa = pg.intraclass_corr(data=long, targets="s", raters="r", ratings="y")
        for r, (_, prow) in zip(rows, pa.iterrows()):
            chk(r["icc"], prow["ICC"], "pingouin ICC " + r["key"], 1e-9)
            chk(r["F"], prow["F"], "pingouin F " + r["key"], 1e-9)
            chk(r["p"], prow["pval"], "pingouin p " + r["key"], 1e-7, 1e-14)
            chk(r["lo"], prow["CI95"][0], "pingouin 하한 " + r["key"], 1e-7, 1e-9)
            chk(r["hi"], prow["CI95"][1], "pingouin 상한 " + r["key"], 1e-7, 1e-9)
    return dict(n=int(n), k=int(k), conf=conf, Y=Y.tolist(), ssb=float(ssb), ssj=float(ssj), sse=float(sse), sst=float(sst), dfb=dfb, dfj=dfj, dfe=dfe,
                msb=float(msb), msj=float(msj), mse=float(mse), msw=float(msw), means=Y.mean(0).tolist(), sds=Y.std(0, ddof=1).tolist(), rows=rows)


def noisy(truth, cats, agree, size=None):
    """진짜 범주에서 agree 확률로 같은 범주, 아니면 임의 범주 (선택적으로 편향)"""
    out = []
    for t in truth:
        out.append(t if rng.random() < agree else cats[rng.integers(0, len(cats))])
    return out


kappa_cases, fleiss_cases, icc_cases = [], [], []
# K1: 이진 (유병률이 한쪽으로 쏠림)
cats = [0, 1]
truth = list(rng.choice(cats, 100, p=[0.85, 0.15]))
a, b = noisy(truth, cats, 0.88), noisy(truth, cats, 0.85)
kappa_cases.append(dict(name="이진 평가 100건, 한쪽 범주가 많음(kappa 역설 가능)", a=a, b=b, **cohen(a, b, cats)))
# K2: 5점 서열
cats = [1, 2, 3, 4, 5]
truth = list(rng.integers(1, 6, 120))
a = [int(min(5, max(1, t + rng.integers(-1, 2)))) for t in truth]
b = [int(min(5, max(1, t + rng.integers(-1, 2) + (1 if rng.random() < 0.2 else 0)))) for t in truth]
kappa_cases.append(dict(name="5점 서열 120건 (가중 kappa)", a=a, b=b, **cohen(a, b, cats)))
# K3: 명목 3범주, 한 평가자는 한 범주를 쓰지 않음 (주변합 0)
cats = [1, 2, 3]
truth = list(rng.choice(cats, 50, p=[0.5, 0.3, 0.2]))
a = noisy(truth, cats, 0.7)
b = [t if t != 3 else 2 for t in noisy(truth, cats, 0.8)]
kappa_cases.append(dict(name="명목 3범주 50건, 평가자 B는 3을 쓰지 않음", a=a, b=b, **cohen(a, b, cats)))
# F1: 3명 3범주 60건
cats = [1, 2, 3]
truth = list(rng.choice(cats, 60, p=[0.4, 0.4, 0.2]))
R = np.array([noisy(truth, cats, 0.75) for _ in range(3)]).T
fleiss_cases.append(dict(name="평가자 3명 × 3범주 60건", R=R.tolist(), **fleiss(R, cats)))
# F2: 4명 4범주 80건
cats = [1, 2, 3, 4]
truth = list(rng.choice(cats, 80))
R = np.array([noisy(truth, cats, 0.6) for _ in range(4)]).T
fleiss_cases.append(dict(name="평가자 4명 × 4범주 80건", R=R.tolist(), **fleiss(R, cats)))
# F3: 2명 (Fleiss kappa 는 Scott 의 π 와 같다)
cats = [0, 1]
truth = list(rng.choice(cats, 70, p=[0.6, 0.4]))
R = np.array([noisy(truth, cats, 0.8) for _ in range(2)]).T
fleiss_cases.append(dict(name="평가자 2명 이진 70건 (Fleiss = Scott π)", R=R.tolist(), **fleiss(R, cats)))

# ICC
def sim_icc(n, k, subj_sd, rater_bias, noise, dec=2):
    s = rng.normal(0, subj_sd, n)
    return np.round(10 + s[:, None] + np.array(rater_bias)[None, :] + rng.normal(0, noise, (n, k)), dec)

icc_cases.append(dict(name="4명 평가자 × 30대상, 평가자 편향 작음", **icc(sim_icc(30, 4, 3, [0, 0.2, -0.1, 0.1], 1.0))))
icc_cases.append(dict(name="2회 측정 × 40대상 (검사-재검사), 일치 높음", **icc(sim_icc(40, 2, 4, [0, 0.1], 0.6))))
icc_cases.append(dict(name="3명 평가자 × 12대상, 평가자 편향 큼", **icc(sim_icc(12, 3, 2, [0, 1.5, 3.0], 0.8))))
icc_cases.append(dict(name="5명 평가자 × 20대상, 일치 낮음", **icc(sim_icc(20, 5, 0.4, [0, 0, 0, 0, 0], 2.0))))
icc_cases.append(dict(name="3회 측정 5점 척도(정수) × 50대상", **icc(np.clip(np.round(sim_icc(50, 3, 1.2, [0, 0.2, 0.4], 0.7) - 7), 1, 5), conf=95)))
icc_cases.append(dict(name="3명 평가자 × 15대상, 90% 신뢰구간", **icc(sim_icc(15, 3, 2, [0, 0.5, 1.0], 1.0), conf=90)))

path = os.path.join(os.path.dirname(__file__), "reference", "agreement.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"statsmodels": statsmodels.__version__, "pingouin": pg.__version__, "scipy": scipy.__version__, "kappa": kappa_cases, "fleiss": fleiss_cases, "icc": icc_cases}, fh, ensure_ascii=False, indent=0, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("wrote kappa", len(kappa_cases), "fleiss", len(fleiss_cases), "icc", len(icc_cases), "(numpy = sklearn = statsmodels = pingouin 일치 확인 통과);",
      [round(c["none"]["kappa"], 3) for c in kappa_cases], [round(c["kappa"], 3) for c in fleiss_cases], [round(c["rows"][1]["icc"], 3) for c in icc_cases])
