"""일원 반복측정 분산분석 기준값을 tests/reference/rmanova.json 에 저장한다.

실행: py -3.12 tests/gen_reference_rmanova.py   (pingouin, statsmodels 가 개발용으로 설치돼 있어야 한다)
기준값은 세 갈래로 얻고 서로 일치하는지 이 스크립트가 확인한다.
 1) numpy 로 정의대로 직접 계산 (이 도구의 JS 구현과 다른 경로: 행렬 대각화 없이 직교 대비의 공분산 행렬과 행렬식 사용)
 2) pingouin (rm_anova, sphericity, epsilon, pairwise_tests) — F, 유의확률, 구형성(Mauchly), Greenhouse–Geisser·Huynh–Feldt·하한 ε, 부분·일반화 η², Bonferroni 사후비교
 3) statsmodels AnovaRM — F, 자유도, 유의확률(구형성 가정)
저장하는 값: 평균·SD, SS(시간·사람·오차), F, p, Mauchly W·χ²·df·p, ε 세 가지와 교정 p, 부분 η², 일반화 η², 쌍별 t 와 Bonferroni 보정 p·신뢰구간
"""
import json
import os
import warnings

import numpy as np
import pandas as pd
import pingouin as pg
import scipy
from scipy import stats
from statsmodels.stats.anova import AnovaRM

warnings.filterwarnings("ignore")
rng = np.random.default_rng(20261038)


def helmert(k):
    C = np.zeros((k - 1, k))
    for j in range(1, k):
        C[j - 1, :j] = 1 / np.sqrt(j * (j + 1))
        C[j - 1, j] = -j / np.sqrt(j * (j + 1))
    return C


def analyse(name, Y, conf=95):
    Y = np.asarray(Y, float)
    n, k = Y.shape
    M = Y.mean()
    ss_time = n * ((Y.mean(0) - M) ** 2).sum()
    ss_subj = k * ((Y.mean(1) - M) ** 2).sum()
    ss_total = ((Y - M) ** 2).sum()
    ss_err = ss_total - ss_time - ss_subj
    df_t, df_e = k - 1, (n - 1) * (k - 1)
    F = (ss_time / df_t) / (ss_err / df_e)
    p = float(stats.f.sf(F, df_t, df_e))
    S = np.cov(Y, rowvar=False)
    C = helmert(k)
    Sc = C @ S @ C.T
    p1 = k - 1
    W = float(np.linalg.det(Sc) / (np.trace(Sc) / p1) ** p1)
    d = 1 - (2 * p1 ** 2 + p1 + 2) / (6 * p1 * (n - 1))
    chi2 = float(-(n - 1) * d * np.log(W))
    dfm = k * (k - 1) // 2 - 1
    # 2차 근사 (R mauchly.test, pingouin 과 같음; k = 3 이면 보정항이 0 이라 1차 근사와 같다)
    # 2차 항의 3k (k = 측정 횟수)는 R 의 mauchly.test 를 따른 것이다 (pingouin 도 같은 값을 낸다; Anderson 의 원식은 3(k − 1)).
    w2 = (p1 + 2) * (p1 - 1) * (p1 - 2) * (2 * p1 ** 3 + 6 * p1 ** 2 + 3 * k + 2) / (288 * ((n - 1) * p1 * d) ** 2)
    pr1, pr2 = float(stats.chi2.sf(chi2, dfm)), float(stats.chi2.sf(chi2, dfm + 4))
    p_m = pr1 + w2 * (pr2 - pr1)
    gg = float(np.trace(Sc) ** 2 / (p1 * np.trace(Sc @ Sc)))
    hf = float(min(1.0, (n * p1 * gg - 2) / (p1 * (n - 1 - p1 * gg))))
    lb = 1 / p1
    corr = {}
    for key, e in [("gg", gg), ("hf", hf), ("lb", lb)]:
        corr[key] = dict(eps=e, df1=e * df_t, df2=e * df_e, p=float(stats.f.sf(F, e * df_t, e * df_e)))
    peta2 = ss_time / (ss_time + ss_err)
    geta2 = ss_time / ss_total
    m = k * (k - 1) / 2
    alpha = 1 - conf / 100
    tcrit = float(stats.t.ppf(1 - alpha / (2 * m), n - 1))
    pairs = []
    for i in range(k):
        for j in range(i + 1, k):
            diff = Y[:, i] - Y[:, j]
            sd = diff.std(ddof=1)
            se = sd / np.sqrt(n)
            t = diff.mean() / se
            pu = float(2 * stats.t.sf(abs(t), n - 1))
            pairs.append(dict(i=i, j=j, diff=float(diff.mean()), se=float(se), t=float(t), df=n - 1, p=pu, p_bonf=min(1.0, pu * m),
                              lo=float(diff.mean() - tcrit * se), hi=float(diff.mean() + tcrit * se)))

    # ---- 교차 확인: pingouin, statsmodels ----
    long = pd.DataFrame({"id": np.repeat(np.arange(n), k), "t": np.tile(np.arange(k), n), "y": Y.ravel()})
    a = pg.rm_anova(long, dv="y", within="t", subject="id", detailed=True, correction=True, effsize="np2").iloc[0]
    a2 = pg.rm_anova(long, dv="y", within="t", subject="id", detailed=True, correction=True, effsize="ng2").iloc[0]
    sp = pg.sphericity(long, dv="y", within="t", subject="id")
    check = lambda x, y, what, tol=1e-7: np.testing.assert_allclose(float(x), float(y), rtol=tol, atol=1e-10, err_msg=name + " " + what)
    check(a["SS"], ss_time, "SS"); check(a["F"], F, "F"); check(a["p_unc"], p, "p"); check(a["np2"], peta2, "np2"); check(a2["ng2"], geta2, "ng2")
    check(a["eps"], gg, "GG eps"); check(a["p_GG_corr"], corr["gg"]["p"], "p_GG")
    check(sp.W, W, "Mauchly W"); check(sp.chi2, chi2, "Mauchly chi2"); check(sp.pval, p_m, "Mauchly p")
    check(pg.epsilon(long, dv="y", within="t", subject="id", correction="hf"), hf, "HF eps")
    check(pg.epsilon(long, dv="y", within="t", subject="id", correction="lb"), lb, "LB eps")
    sm = AnovaRM(long, "y", "id", within=["t"]).fit().anova_table.iloc[0]
    check(sm["F Value"], F, "statsmodels F"); check(sm["Pr > F"], p, "statsmodels p")
    pw = pg.pairwise_tests(long, dv="y", within="t", subject="id", padjust="bonf")
    for row, pr in zip(pw.itertuples(), pairs):
        check(row.T, -pr["t"], "pairwise t (pingouin 은 A−B 가 아니라 부호가 다를 수 있어 절댓값으로)") if False else check(abs(row.T), abs(pr["t"]), "pairwise |t|")
        check(row.p_corr, pr["p_bonf"], "pairwise p_bonf")

    return dict(name=name, n=n, k=k, conf=conf, data=Y.tolist(), means=Y.mean(0).tolist(), sds=Y.std(0, ddof=1).tolist(),
                ss_time=float(ss_time), ss_subj=float(ss_subj), ss_err=float(ss_err), ss_total=float(ss_total), df_time=df_t, df_err=df_e,
                ms_time=float(ss_time / df_t), ms_err=float(ss_err / df_e), F=float(F), p=p,
                mauchly=dict(W=W, chi2=chi2, df=dfm, p=p_m, p_first_order=pr1, w2=float(w2)), corr=corr, peta2=float(peta2), geta2=float(geta2), pairs=pairs, tcrit=tcrit)


cases = []
# A: 구형성이 대체로 성립 (20명 × 4회)
Y = rng.normal(0, 1, (20, 1)) * 1.2 + np.array([0, 0.5, 1.0, 1.5]) + rng.normal(0, 1, (20, 4))
cases.append(analyse("20명 × 4회, 구형성 대체로 성립", np.round(Y, 3)))
# B: 구형성 위반 (조건마다 오차 크기가 크게 다름, 25명 × 4회)
Y = rng.normal(0, 1, (25, 1)) + np.array([0, 0.3, 0.9, 1.6]) + rng.normal(0, 1, (25, 4)) * np.array([0.5, 0.5, 2.5, 4.0])
cases.append(analyse("25명 × 4회, 구형성 위반", np.round(Y, 3)))
# C: 12명 × 5회
Y = rng.normal(0, 1.5, (12, 1)) + np.array([0, 0.4, 0.8, 1.0, 1.4]) + rng.normal(0, 1, (12, 5))
cases.append(analyse("12명 × 5회", np.round(Y, 3)))
# D: 15명 × 3회, 시간에 따른 오차 증가
Y = rng.normal(0, 1, (15, 1)) + np.array([0, 0.6, 1.1]) + rng.normal(0, 1, (15, 3)) * np.array([0.6, 1.0, 2.0])
cases.append(analyse("15명 × 3회, 오차가 점점 커짐", np.round(Y, 2)))
# E: 8명 × 3회 (소표본), 신뢰수준 99
Y = rng.normal(50, 6, (8, 1)) + np.array([0, 2.0, 5.0]) + rng.normal(0, 3, (8, 3))
cases.append(analyse("8명 × 3회, 99% 구간", np.round(Y, 1), conf=99))

# F: 구형성이 깨졌지만 ε(G-G) > .75 인 경우 (Huynh–Feldt 행을 읽는 경로) — 조건을 만족하는 자료를 찾는다
for attempt in range(200):
    Y = np.round(rng.normal(0, 1, (60, 1)) + np.array([0, 0.3, 0.6, 0.9]) + rng.normal(0, 1, (60, 4)) * np.array([1.0, 1.0, 1.6, 1.0]), 3)
    S_ = np.cov(Y, rowvar=False)
    C_ = helmert(4)
    Sc_ = C_ @ S_ @ C_.T
    eps_ = np.trace(Sc_) ** 2 / (3 * np.trace(Sc_ @ Sc_))
    p1_, n_ = 3, 60
    W_ = np.linalg.det(Sc_) / (np.trace(Sc_) / 3) ** 3
    d_ = 1 - (2 * 9 + 3 + 2) / (6 * 3 * (n_ - 1))
    pm_ = stats.chi2.sf(-(n_ - 1) * d_ * np.log(W_), 5)
    if eps_ > 0.80 and pm_ < 0.04:
        break
else:
    raise SystemExit("조건을 만족하는 자료를 못 찾음")
cases.append(analyse("60명 × 4회, 구형성 위반이지만 ε(G-G) > .75", Y))

path = os.path.join(os.path.dirname(__file__), "reference", "rmanova.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"pingouin": pg.__version__, "scipy": scipy.__version__, "cases": cases}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "rmanova cases (numpy = pingouin = statsmodels 일치 확인 통과);",
      [(c["n"], c["k"], round(c["F"], 3), round(c["mauchly"]["p"], 4), round(c["corr"]["gg"]["eps"], 3), round(c["corr"]["hf"]["eps"], 3)) for c in cases])
