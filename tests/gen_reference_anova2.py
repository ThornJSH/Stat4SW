"""이원 분산분석 기준값을 tests/reference/anova2.json 에 저장한다.

실행: py -3.12 tests/gen_reference_anova2.py   (statsmodels, pingouin 이 개발용으로 설치돼 있어야 한다)
기준값은 세 갈래로 얻고 서로 일치하는지 이 스크립트가 확인한다.
 1) numpy 로 정의대로 직접 계산 (효과 코딩 설계 행렬의 lstsq 로 축소·전체 모형 잔차제곱합; 유형 II·III, 균형 자료는 닫힌 식 제곱합도)
 2) statsmodels: OLS(합 0 대비) + anova_lm(typ=2, 3), 추정 주변평균은 t_test(선형결합), 단순 주효과는 셀 평균 대비의 f_test
 3) pingouin: anova(between=[A, B], ss_type=2, 3) 의 제곱합·F·p·부분 η²
사후검정은 scipy.stats.studentized_range (Tukey) 와 t 분포(Bonferroni), Levene 은 scipy.stats.levene(center='mean').
저장하는 값: 셀 N·평균·SD, 주변 평균, 유형 III·II 의 항별 제곱합·자유도·F·p·부분 η², 절편(유형 III), 오차, R², 추정 주변평균, 사후검정, 단순 주효과, Levene
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
from patsy import dmatrix
from scipy import stats
from statsmodels.stats.anova import anova_lm

warnings.filterwarnings("ignore")
rng = np.random.default_rng(20261030)


def effect_cols(idx, k):
    return [np.where(idx == l, 1.0, np.where(idx == k - 1, -1.0, 0.0)) for l in range(k - 1)]


def sse_of(cols, y):
    n = len(y)
    X = np.column_stack([np.ones(n)] + cols) if cols else np.ones((n, 1))
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    return float(r @ r)


def analyse(name, y, ia, ib, a, b, conf=95):
    y = np.asarray(y, float)
    ia = np.asarray(ia) - 1
    ib = np.asarray(ib) - 1
    n = len(y)
    df = pd.DataFrame(dict(y=y, A=pd.Categorical(ia + 1), B=pd.Categorical(ib + 1)))
    sst = float(((y - y.mean()) ** 2).sum())
    A, B = effect_cols(ia, a), effect_cols(ib, b)
    AB = [ca * cb for ca in A for cb in B]
    sse_full = sse_of(A + B + AB, y)
    dfe = n - a * b
    mse = sse_full / dfe
    dfs = dict(A=a - 1, B=b - 1, AB=(a - 1) * (b - 1))
    ss3 = dict(A=sse_of(B + AB, y) - sse_full, B=sse_of(A + AB, y) - sse_full, AB=sse_of(A + B, y) - sse_full)
    sse_ab = sse_of(A + B, y)
    ss2 = dict(A=sse_of(B, y) - sse_ab, B=sse_of(A, y) - sse_ab, AB=sse_ab - sse_full)

    def pack(ss):
        out = {}
        for k in ("A", "B", "AB"):
            F = (ss[k] / dfs[k]) / mse
            out[k] = dict(ss=ss[k], df=dfs[k], ms=ss[k] / dfs[k], F=F, p=float(stats.f.sf(F, dfs[k], dfe)), peta=ss[k] / (ss[k] + sse_full))
        return out

    # 셀 통계
    cells = [[None] * b for _ in range(a)]
    for i in range(a):
        for j in range(b):
            v = y[(ia == i) & (ib == j)]
            cells[i][j] = dict(n=int(len(v)), mean=float(v.mean()), sd=float(v.std(ddof=1)) if len(v) > 1 else None)
    row_n = [int((ia == i).sum()) for i in range(a)]
    col_n = [int((ib == j).sum()) for j in range(b)]
    row_mean = [float(y[ia == i].mean()) for i in range(a)]
    col_mean = [float(y[ib == j].mean()) for j in range(b)]
    cm = np.array([[cells[i][j]["mean"] for j in range(b)] for i in range(a)])
    cn = np.array([[cells[i][j]["n"] for j in range(b)] for i in range(a)])
    # 균형 자료의 닫힌 식 제곱합 (유형 II = III = I)
    balanced = len(set(cn.ravel())) == 1
    if balanced:
        m = cn[0, 0]
        gm = y.mean()
        np.testing.assert_allclose(ss3["A"], m * b * ((cm.mean(1) - gm) ** 2).sum(), rtol=1e-8, err_msg=name + " 균형 SS_A")
        np.testing.assert_allclose(ss3["B"], m * a * ((cm.mean(0) - gm) ** 2).sum(), rtol=1e-8, err_msg=name + " 균형 SS_B")
        np.testing.assert_allclose(ss3["AB"], m * ((cm - cm.mean(1, keepdims=True) - cm.mean(0, keepdims=True) + gm) ** 2).sum(), rtol=1e-6, atol=1e-9, err_msg=name + " 균형 SS_AB")
        np.testing.assert_allclose(ss3["A"], ss2["A"], rtol=1e-8)
    # 절편 (유형 III): μ̄² (ab)² / Σ(1/n_ij)
    mu = cm.mean()
    ss_int = float(mu ** 2 * (a * b) ** 2 / (1.0 / cn).sum())

    # ---- statsmodels ----
    mod = smf.ols("y ~ C(A, Sum) * C(B, Sum)", data=df).fit()
    t3 = anova_lm(mod, typ=3)
    t2 = anova_lm(mod, typ=2)
    chk = lambda x, z, what, tol=1e-8, atol=1e-10: np.testing.assert_allclose(np.asarray(x, float), np.asarray(z, float), rtol=tol, atol=atol, err_msg=f"{name} {what}")
    for key, row in (("A", "C(A, Sum)"), ("B", "C(B, Sum)"), ("AB", "C(A, Sum):C(B, Sum)")):
        chk(t3.loc[row, "sum_sq"], ss3[key], f"statsmodels III SS {key}")
        chk(t3.loc[row, "PR(>F)"], pack(ss3)[key]["p"], f"statsmodels III p {key}", 1e-6, 1e-12)
        chk(t2.loc[row, "sum_sq"], ss2[key], f"statsmodels II SS {key}")
        chk(t2.loc[row, "PR(>F)"], pack(ss2)[key]["p"], f"statsmodels II p {key}", 1e-6, 1e-12)
    chk(t3.loc["Intercept", "sum_sq"], ss_int, "statsmodels III 절편")
    chk(t3.loc["Residual", "sum_sq"], sse_full, "statsmodels 오차")
    chk(mod.rsquared, 1 - sse_full / sst, "R2")
    chk(mod.rsquared_adj, 1 - (1 - (1 - sse_full / sst)) * (n - 1) / dfe, "adj R2")
    # ---- pingouin ----
    for typ, ss in ((3, ss3), (2, ss2)):
        pa = pg.anova(data=df, dv="y", between=["A", "B"], ss_type=typ, detailed=True, effsize="np2").set_index("Source")
        for key, row in (("A", "A"), ("B", "B"), ("AB", "A * B")):
            chk(pa.loc[row, "SS"], ss[key], f"pingouin type{typ} SS {key}")
            chk(pa.loc[row, "np2"], ss[key] / (ss[key] + sse_full), f"pingouin type{typ} np2 {key}")
            chk(pa.loc[row, "p_unc"], pack(ss)[key]["p"], f"pingouin type{typ} p {key}", 1e-6, 1e-12)

    # ---- 추정 주변평균 (statsmodels t_test 로 선형결합의 표준오차) ----
    grid = pd.DataFrame([dict(A=pd.Categorical([i + 1], categories=list(range(1, a + 1)))[0], B=pd.Categorical([j + 1], categories=list(range(1, b + 1)))[0]) for i in range(a) for j in range(b)])
    Xc = np.asarray(dmatrix("C(A, Sum) * C(B, Sum)", grid, return_type="dataframe").reindex(columns=mod.params.index)).reshape(a, b, -1)
    tc = stats.t.ppf(1 - (1 - conf / 100) / 2, dfe)

    def emm_of(L):
        tt = mod.t_test(L.reshape(1, -1))
        est, se = float(np.squeeze(tt.effect)), float(np.squeeze(tt.sd))
        return dict(mean=est, se=se, lo=est - tc * se, hi=est + tc * se)

    emm = dict(rows=[emm_of(Xc[i].mean(0)) for i in range(a)], cols=[emm_of(Xc[:, j].mean(0)) for j in range(b)],
               cells=[[emm_of(Xc[i, j]) for j in range(b)] for i in range(a)])
    # 닫힌 식과 대조
    for i in range(a):
        chk(emm["rows"][i]["mean"], cm[i].mean(), "EMM 평균(행)")
        chk(emm["rows"][i]["se"], np.sqrt(mse * (1.0 / cn[i]).sum()) / b, "EMM SE(행)")
    for j in range(b):
        chk(emm["cols"][j]["se"], np.sqrt(mse * (1.0 / cn[:, j]).sum()) / a, "EMM SE(열)")

    # ---- 사후검정 (관측 주변평균 + 이원 모형의 MSE) ----
    def post(means, ns, k):
        rows = []
        qc = stats.studentized_range.ppf(1 - (1 - conf / 100), k, dfe)
        tcb = stats.t.ppf(1 - (1 - conf / 100) / (2 * k * (k - 1) / 2), dfe)
        for i in range(k - 1):
            for j in range(i + 1, k):
                diff = means[i] - means[j]
                se = float(np.sqrt(mse * (1 / ns[i] + 1 / ns[j])))
                t = diff / se
                rows.append(dict(i=i, j=j, diff=float(diff), se=se,
                                 p_tukey=float(stats.studentized_range.sf(abs(t) * np.sqrt(2), k, dfe)),
                                 p_bonf=float(min(1.0, k * (k - 1) / 2 * 2 * stats.t.sf(abs(t), dfe))),
                                 lo_t=float(diff - qc * se / np.sqrt(2)), hi_t=float(diff + qc * se / np.sqrt(2)),
                                 lo_b=float(diff - tcb * se), hi_b=float(diff + tcb * se)))
        return rows

    posthoc = dict(A=post(row_mean, row_n, a) if a >= 3 else None, B=post(col_mean, col_n, b) if b >= 3 else None)

    # ---- 단순 주효과: 닫힌 식 제곱합, F·p 는 statsmodels f_test 의 셀 평균 대비 ----
    def one(means, ns):
        g = (np.asarray(means) * np.asarray(ns)).sum() / np.sum(ns)
        return float((np.asarray(ns) * (np.asarray(means) - g) ** 2).sum())

    def contrast_test(rows):
        R = np.array([rows[k] - rows[0] for k in range(1, len(rows))])
        ft = mod.f_test(R)
        return float(np.squeeze(ft.fvalue)), float(np.squeeze(ft.pvalue)), int(ft.df_num)

    simple = dict(aAtB=[], bAtA=[])
    for j in range(b):
        ssv = one(cm[:, j], cn[:, j])
        F, p, dn = contrast_test([Xc[i, j] for i in range(a)])
        chk(F, (ssv / (a - 1)) / mse, f"단순 주효과 A@B{j} F")
        simple["aAtB"].append(dict(ss=ssv, df=dn, F=F, p=p))
    for i in range(a):
        ssv = one(cm[i], cn[i])
        F, p, dn = contrast_test([Xc[i, j] for j in range(b)])
        chk(F, (ssv / (b - 1)) / mse, f"단순 주효과 B@A{i} F")
        simple["bAtA"].append(dict(ss=ssv, df=dn, F=F, p=p))

    # ---- Levene (평균 기준) ----
    groups = [y[(ia == i) & (ib == j)] for i in range(a) for j in range(b)]
    W, pl = stats.levene(*groups, center="mean")

    return dict(name=name, a=a, b=b, n=n, conf=conf, balanced=bool(balanced), A=(ia + 1).tolist(), B=(ib + 1).tolist(), y=y.tolist(),
                cells=cells, row_n=row_n, col_n=col_n, row_mean=row_mean, col_mean=col_mean, sst=sst, sumsq=float((y ** 2).sum()),
                sse=sse_full, dfe=dfe, mse=mse, r2=1 - sse_full / sst, adj=1 - (1 - (1 - sse_full / sst)) * (n - 1) / dfe,
                III=dict(terms=pack(ss3), intercept=dict(ss=ss_int, F=ss_int / mse, p=float(stats.f.sf(ss_int / mse, 1, dfe))),
                         model=dict(ss=sst - sse_full, df=a * b - 1, F=((sst - sse_full) / (a * b - 1)) / mse, p=float(stats.f.sf(((sst - sse_full) / (a * b - 1)) / mse, a * b - 1, dfe)), peta=(sst - sse_full) / sst)),
                II=dict(terms=pack(ss2)), emm=emm, posthoc=posthoc, simple=simple, levene=dict(W=float(W), df1=a * b - 1, df2=n - a * b, p=float(pl)))


def make(name, sizes, mu, sd, a, b, dec=2):
    ia, ib, y = [], [], []
    for i in range(a):
        for j in range(b):
            for _ in range(sizes[i][j]):
                ia.append(i + 1); ib.append(j + 1); y.append(round(float(rng.normal(mu[i][j], sd)), dec))
    return name, y, ia, ib, a, b


cases = []
# 1: 균형 2×3, n=10, 상호작용 있음
cases.append(analyse(*make("균형 2×3 (셀당 10건), 상호작용 있음", [[10] * 3] * 2, [[50, 55, 60], [52, 51, 49]], 5, 2, 3)))
# 2: 불균형 3×3
sz = rng.integers(4, 13, (3, 3)).tolist()
cases.append(analyse(*make("불균형 3×3 (셀당 4~12건)", sz, [[10, 12, 14], [11, 15, 13], [9, 10, 16]], 3, 3, 3)))
# 3: 균형 2×2, n=15, 상호작용 없음
cases.append(analyse(*make("균형 2×2 (셀당 15건), 상호작용 없음", [[15] * 2] * 2, [[20, 24], [26, 30]], 4, 2, 2)))
# 4: 불균형 2×2, 작은 셀
cases.append(analyse(*make("불균형 2×2 (셀당 3~9건)", [[3, 9], [7, 4]], [[5.0, 6.5], [7.0, 6.0]], 1.5, 2, 2)))
# 5: 불균형 4×3
sz5 = rng.integers(5, 10, (4, 3)).tolist()
cases.append(analyse(*make("불균형 4×3 (셀당 5~9건)", sz5, [[30, 32, 35], [31, 36, 33], [38, 34, 30], [29, 33, 40]], 4, 4, 3, dec=1)))

path = os.path.join(os.path.dirname(__file__), "reference", "anova2.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"statsmodels": statsmodels.__version__, "pingouin": pg.__version__, "scipy": scipy.__version__, "cases": cases}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "anova2 cases (numpy = statsmodels = pingouin 일치 확인 통과);",
      [(c["n"], c["balanced"], round(c["III"]["terms"]["AB"]["F"], 3), round(c["III"]["terms"]["A"]["ss"] - c["II"]["terms"]["A"]["ss"], 3)) for c in cases])
