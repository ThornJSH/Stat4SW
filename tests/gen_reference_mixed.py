"""혼합 설계 분산분석(집단 1개 × 반복측정 1개) 기준값을 tests/reference/mixed.json 에 저장한다.

실행: py -3.12 tests/gen_reference_mixed.py   (statsmodels, pingouin 이 개발용으로 설치돼 있어야 한다)
기준값은 세 갈래로 얻고 서로 일치하는지 이 스크립트가 확인한다.
 1) numpy 로 정의대로 직접 계산: 사람별 평균(√k 배)으로 집단 효과, 직교(Helmert) 대비 점수에서
    시간 = 절편(비가중 평균, 유형 III), 시간×집단 = 집단 효과, 오차 = 집단 안 잔차; 구형성은 집단 안에서 합친 공분산(자유도 N − g)
 2) statsmodels: 대비 점수·사람별 평균에 대한 OLS(합 0 대비) + anova_lm(typ=3) — 집단 크기가 달라도(불균형) 유형 III 제곱합을 독립적으로 확인
 3) pingouin: mixed_anova (균형 자료의 제곱합·F·p·부분 η²·ε·Mauchly), 집단별 rm_anova·pairwise_tests, scipy f_oneway·levene
Huynh–Feldt ε 는 Lecoutre 보정식 ((N−g+1)·p·ε − 2)/(p·(N−g − p·ε)) (R car) — 다른 도구로 대조하지 못해 numpy 식만 쓴다.
저장하는 값: 집단·시점 평균·SD·N, 개체-간(절편·집단·오차), 개체-내(시간·시간×집단·오차), Mauchly, ε 3종, 교정 p, 시점별 Levene,
집단별 시간 효과(rm ANOVA, 쌍별 비교), 시점별 집단 효과(일원 ANOVA)
"""
import json
import os
import warnings

import numpy as np
import pandas as pd
import pingouin as pg
import scipy
import statsmodels
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.anova import anova_lm

warnings.filterwarnings("ignore")
rng = np.random.default_rng(20261031)


def helmert(k):
    C = np.zeros((k - 1, k))
    for a in range(1, k):
        C[a - 1, :a] = 1 / np.sqrt(a * (a + 1))
        C[a - 1, a] = -a / np.sqrt(a * (a + 1))
    return C


def mauchly(Sig, dfe, k):
    p = k - 1
    tr = np.trace(Sig)
    W = float(np.linalg.det(Sig) / (tr / p) ** p)
    d = 1 - (2 * p * p + p + 2) / (6 * p * dfe)
    chi2 = float(-dfe * d * np.log(W))
    df = k * (k - 1) // 2 - 1
    w2 = (p + 2) * (p - 1) * (p - 2) * (2 * p ** 3 + 6 * p ** 2 + 3 * k + 2) / (288 * (dfe * p * d) ** 2)
    q1, q2 = stats.chi2.sf(chi2, df), stats.chi2.sf(chi2, df + 4)
    return dict(W=W, chi2=chi2, df=df, p=float(q1 + w2 * (q2 - q1)))


def analyse(name, Y, G, g, conf=95):
    Y = np.asarray(Y, float)
    G = np.asarray(G) - 1
    N, k = Y.shape
    p = k - 1
    nh = np.array([(G == h).sum() for h in range(g)])
    Hm = helmert(k)
    Z = Y @ Hm.T
    zbar = np.array([Z[G == h].mean(0) for h in range(g)])
    zw = Z.mean(0)
    inv = (1.0 / nh).sum()
    # 개체-간 (사람별 평균 × √k)
    w = np.sqrt(k) * Y.mean(1)
    wh = np.array([w[G == h].mean() for h in range(g)])
    ss_group = float(sum(nh[h] * (wh[h] - w.mean()) ** 2 for h in range(g)))
    ss_errb = float(((w - wh[G]) ** 2).sum())
    ss_int = float(g ** 2 * wh.mean() ** 2 / inv)
    df_g, df_eb = g - 1, N - g
    # 개체-내
    ss_time = float((g ** 2 * zbar.mean(0) ** 2 / inv).sum())
    ss_tg = float(sum(nh[h] * ((zbar[h] - zw) ** 2).sum() for h in range(g)))
    res = Z - zbar[G]
    ss_errw = float((res ** 2).sum())
    df_t, df_tg, df_ew = p, (g - 1) * p, (N - g) * p
    Sig = res.T @ res / (N - g)
    tr = np.trace(Sig)
    gg = float(tr ** 2 / (p * np.trace(Sig @ Sig)))
    dfe = N - g
    hf = float(min(1.0, ((dfe + 1) * p * gg - 2) / (p * (dfe - p * gg)))) if p > 1 else 1.0
    eps = dict(gg=gg, hf=max(hf, gg) if p > 1 else 1.0, lb=1 / p)
    mau = mauchly(Sig, dfe, k) if (k >= 3 and N - g >= p) else None

    def eff(ss, df, ss_err, df_err):
        F = (ss / df) / (ss_err / df_err)
        out = dict(ss=ss, df=df, ms=ss / df, F=float(F), p=float(stats.f.sf(F, df, df_err)), peta=ss / (ss + ss_err),
                   corr={key: dict(eps=e, df1=e * df, df2=e * df_err, p=float(stats.f.sf(F, e * df, e * df_err))) for key, e in eps.items()})
        return out

    out = dict(name=name, g=g, k=k, N=N, conf=conf, Y=Y.tolist(), G=(G + 1).tolist(), nh=nh.tolist(),
               between=dict(intercept=dict(ss=ss_int, df=1, F=ss_int / (ss_errb / df_eb), p=float(stats.f.sf(ss_int / (ss_errb / df_eb), 1, df_eb))),
                            group=dict(ss=ss_group, df=df_g, ms=ss_group / df_g, F=(ss_group / df_g) / (ss_errb / df_eb), p=float(stats.f.sf((ss_group / df_g) / (ss_errb / df_eb), df_g, df_eb)), peta=ss_group / (ss_group + ss_errb)),
                            err=dict(ss=ss_errb, df=df_eb, ms=ss_errb / df_eb)),
               time=eff(ss_time, df_t, ss_errw, df_ew), inter=eff(ss_tg, df_tg, ss_errw, df_ew), errw=dict(ss=ss_errw, df=df_ew, ms=ss_errw / df_ew),
               eps=eps, mauchly=mau, sigma=Sig.tolist(),
               cells=[[dict(n=int(nh[h]), mean=float(Y[G == h][:, j].mean()), sd=float(Y[G == h][:, j].std(ddof=1)) if nh[h] > 1 else None) for j in range(k)] for h in range(g)],
               col_mean=Y.mean(0).tolist(), col_sd=Y.std(0, ddof=1).tolist())

    # ---- 시점별 Levene, 시점별 집단 효과(일원 ANOVA) ----
    out["levene"] = []
    out["group_at_time"] = []
    for j in range(k):
        grp = [Y[G == h][:, j] for h in range(g)]
        W, pl = stats.levene(*grp, center="mean")
        out["levene"].append(dict(W=float(W), df1=g - 1, df2=N - g, p=float(pl)))
        Fo, po = stats.f_oneway(*grp)
        out["group_at_time"].append(dict(F=float(Fo), df1=g - 1, df2=N - g, p=float(po)))

    # ---- 집단별 시간 효과(반복측정 일원 ANOVA)와 쌍별 비교 ----
    out["time_in_group"] = []
    m_pairs = k * (k - 1) / 2
    for h in range(g):
        Yh = Y[G == h]
        n = len(Yh)
        gm = Yh.mean()
        ss_t = n * ((Yh.mean(0) - gm) ** 2).sum()
        ss_s = k * ((Yh.mean(1) - gm) ** 2).sum()
        ss_e = ((Yh - gm) ** 2).sum() - ss_t - ss_s
        dft, dfe_h = k - 1, (n - 1) * (k - 1)
        F = (ss_t / dft) / (ss_e / dfe_h)
        rec = dict(n=n, ss=float(ss_t), ss_err=float(ss_e), df1=dft, df2=dfe_h, F=float(F), p=float(stats.f.sf(F, dft, dfe_h)))
        if k >= 3:
            long = pd.DataFrame({"id": np.repeat(np.arange(n), k), "t": np.tile(np.arange(k), n), "y": Yh.ravel()})
            a = pg.rm_anova(long, dv="y", within="t", subject="id", detailed=True, correction=True).iloc[0]
            np.testing.assert_allclose(a["F"], F, rtol=1e-7, err_msg=name + " 집단별 rm F (pingouin)")
            np.testing.assert_allclose(a["p_unc"], rec["p"], rtol=1e-6, atol=1e-12, err_msg=name + " 집단별 rm p")
            rec["p_gg"] = float(a["p_GG_corr"])
        pairs = []
        tcrit = stats.t.ppf(1 - (1 - conf / 100) / (2 * m_pairs), n - 1)
        for i in range(k):
            for j in range(i + 1, k):
                d = Yh[:, i] - Yh[:, j]
                se = d.std(ddof=1) / np.sqrt(n)
                t = d.mean() / se
                pu = float(2 * stats.t.sf(abs(t), n - 1))
                pairs.append(dict(i=i, j=j, diff=float(d.mean()), se=float(se), t=float(t), df=n - 1, p=pu, p_bonf=min(1.0, pu * m_pairs),
                                  lo=float(d.mean() - tcrit * se), hi=float(d.mean() + tcrit * se)))
        rec["pairs"] = pairs
        out["time_in_group"].append(rec)

    # ---- 교차 확인: statsmodels (대비 점수와 사람별 평균에 유형 III) ----
    gdf = pd.Categorical(G + 1)
    chk = lambda x, y, what, tol=1e-8, atol=1e-10: np.testing.assert_allclose(np.asarray(x, float), np.asarray(y, float), rtol=tol, atol=atol, err_msg=f"{name} {what}")
    sm_time = sm_tg = sm_err = 0.0
    for c in range(p):
        d = pd.DataFrame(dict(z=Z[:, c], grp=gdf))
        t3 = anova_lm(smf.ols("z ~ C(grp, Sum)", data=d).fit(), typ=3)
        sm_time += t3.loc["Intercept", "sum_sq"]
        sm_tg += t3.loc["C(grp, Sum)", "sum_sq"]
        sm_err += t3.loc["Residual", "sum_sq"]
    chk(sm_time, ss_time, "statsmodels 시간 SS"); chk(sm_tg, ss_tg, "statsmodels 시간×집단 SS"); chk(sm_err, ss_errw, "statsmodels 개체-내 오차 SS")
    d = pd.DataFrame(dict(z=w, grp=gdf))
    t3 = anova_lm(smf.ols("z ~ C(grp, Sum)", data=d).fit(), typ=3)
    chk(t3.loc["Intercept", "sum_sq"], ss_int, "statsmodels 절편 SS"); chk(t3.loc["C(grp, Sum)", "sum_sq"], ss_group, "statsmodels 집단 SS"); chk(t3.loc["Residual", "sum_sq"], ss_errb, "statsmodels 개체-간 오차 SS")

    # ---- 교차 확인: pingouin mixed_anova (균형 자료) ----
    if len(set(nh.tolist())) == 1:
        long = pd.DataFrame([dict(id=i, grp=int(G[i]) + 1, t=j, y=Y[i, j]) for i in range(N) for j in range(k)])
        pa = pg.mixed_anova(data=long, dv="y", within="t", between="grp", subject="id", correction=True, effsize="np2").set_index("Source")
        chk(pa.loc["grp", "SS"], ss_group, "pingouin 집단 SS"); chk(pa.loc["grp", "F"], out["between"]["group"]["F"], "pingouin 집단 F")
        chk(pa.loc["t", "SS"], ss_time, "pingouin 시간 SS"); chk(pa.loc["t", "F"], out["time"]["F"], "pingouin 시간 F")
        chk(pa.loc["Interaction", "SS"], ss_tg, "pingouin 상호작용 SS"); chk(pa.loc["Interaction", "F"], out["inter"]["F"], "pingouin 상호작용 F")
        chk(pa.loc["grp", "np2"], out["between"]["group"]["peta"], "pingouin 집단 np2"); chk(pa.loc["t", "np2"], out["time"]["peta"], "pingouin 시간 np2")
        chk(pa.loc["Interaction", "np2"], out["inter"]["peta"], "pingouin 상호작용 np2")
        chk(pa.loc["t", "p_unc"], out["time"]["p"], "pingouin 시간 p", 1e-6, 1e-12); chk(pa.loc["Interaction", "p_unc"], out["inter"]["p"], "pingouin 상호작용 p", 1e-6, 1e-12)
        if k >= 3:
            chk(pa.loc["t", "eps"], gg, "pingouin ε(GG, 집단 안 합친 공분산)")
            chk(pa.loc["t", "p_GG_corr"], out["time"]["corr"]["gg"]["p"], "pingouin 시간 p(GG)", 1e-6, 1e-12)
            chk(pa.loc["Interaction", "p_GG_corr"], out["inter"]["corr"]["gg"]["p"], "pingouin 상호작용 p(GG)", 1e-6, 1e-12)
            chk(pa.loc["t", "W_spher"], mau["W"], "pingouin Mauchly W", 1e-7)
            chk(pa.loc["t", "p_spher"], mau["p"], "pingouin Mauchly p", 1e-6, 1e-12)
    return out


def sim(name, ns, k, shifts, slopes, sub_sd, noise, g):
    Y, G = [], []
    for h in range(g):
        for _ in range(ns[h]):
            s = rng.normal(0, sub_sd)
            Y.append([round(10 + shifts[h] + slopes[h] * j + s + rng.normal(0, 1) * noise[j], 2) for j in range(k)])
            G.append(h + 1)
    return name, Y, G, g


cases = []
cases.append(analyse(*sim("2집단 × 3시점 균형 (집단당 12명)", [12, 12], 3, [0, 1], [0.5, 1.3], 2, [1.5, 1.5, 1.5], 2)))
cases.append(analyse(*sim("3집단 × 4시점 불균형 (8, 12, 10명), 오차가 점점 커짐", [8, 12, 10], 4, [0, 1, 2], [0.2, 0.8, 1.5], 2, [0.6, 1.0, 2.0, 3.5], 3)))
cases.append(analyse(*sim("2집단 × 2시점 균형 (집단당 15명)", [15, 15], 2, [0, 0.5], [0.0, 1.6], 2, [1.2, 1.2], 2)))
cases.append(analyse(*sim("3집단 × 3시점 균형 (집단당 10명), 상호작용 강함", [10, 10, 10], 3, [0, 0, 0], [-1.0, 0.0, 1.5], 2.5, [1.3, 1.3, 1.3], 3)))
cases.append(analyse(*sim("2집단 × 4시점 불균형 (7, 14명)", [7, 14], 4, [0, 0.4], [0.4, 0.7], 1.5, [1.0, 1.4, 1.2, 1.6], 2)))

path = os.path.join(os.path.dirname(__file__), "reference", "mixed.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"statsmodels": statsmodels.__version__, "pingouin": pg.__version__, "scipy": scipy.__version__, "cases": cases}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "mixed cases (numpy = statsmodels = pingouin 일치 확인 통과);",
      [(c["N"], c["g"], c["k"], round(c["inter"]["F"], 3), None if c["mauchly"] is None else round(c["mauchly"]["p"], 4), round(c["eps"]["gg"], 3)) for c in cases])
