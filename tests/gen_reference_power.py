"""검정력·표본 크기 기준값을 tests/reference/power.json 에 저장한다.

실행: py -3.12 tests/gen_reference_power.py   (statsmodels, pingouin 이 개발용으로 설치돼 있어야 한다)
기준값은 세 갈래로 얻고 서로 일치하는지 이 스크립트가 확인한다.
 1) scipy 분포로 직접: nct(독립·대응·일표본 t), ncf(분산분석·회귀), norm(상관 Fisher z 보정판, 두 비율 Cohen h)
 2) statsmodels.stats.power: TTestIndPower, TTestPower, FTestAnovaPower, NormalIndPower, FTestPower — 검정력과 solve_power(효과크기)
 3) pingouin: power_ttest, power_ttest2n 대신 power_ttest(paired/one-sample)·power_anova·power_corr
표본 크기는 "목표 검정력 이상이 되는 가장 작은 정수"를 statsmodels 의 검정력 함수로 1씩 늘려 찾아 저장한다.
비중심 F 누적확률은 scipy ncf.cdf 와 직접 대조하는 격자도 저장한다.
"""
import json
import os
import warnings

import numpy as np
import pingouin as pg
import scipy
import statsmodels
from scipy import stats
from statsmodels.stats.power import FTestAnovaPower, FTestPower, NormalIndPower, TTestIndPower, TTestPower

warnings.filterwarnings("ignore")
chk = lambda x, y, what, tol=1e-7, atol=1e-10: np.testing.assert_allclose(np.asarray(x, float), np.asarray(y, float), rtol=tol, atol=atol, err_msg=what)
alt = lambda tails: "two-sided" if tails == 2 else "larger"


def p_t2(d, n1, ratio, alpha, tails):
    n2 = int(np.ceil(ratio * n1))
    df = n1 + n2 - 2
    ncp = d * np.sqrt(n1 * n2 / (n1 + n2))
    if tails == 2:
        tc = stats.t.isf(alpha / 2, df)
        return float(stats.nct.sf(tc, df, ncp) + stats.nct.cdf(-tc, df, ncp))
    return float(stats.nct.sf(stats.t.isf(alpha, df), df, ncp))


def p_t1(d, n, alpha, tails):
    df, ncp = n - 1, d * np.sqrt(n)
    if tails == 2:
        tc = stats.t.isf(alpha / 2, df)
        return float(stats.nct.sf(tc, df, ncp) + stats.nct.cdf(-tc, df, ncp))
    return float(stats.nct.sf(stats.t.isf(alpha, df), df, ncp))


def p_anova(f, n, k, alpha):
    N = k * n
    df1, df2 = k - 1, N - k
    return float(stats.ncf.sf(stats.f.isf(alpha, df1, df2), df1, df2, f * f * N))


def p_corr(r, n, alpha, tails):
    """Fisher z 근사 (임계 상관 r_c 와 편향 보정 r/(2(n−1)) 포함 — pingouin power_corr 와 같은 식)"""
    dof = n - 2
    ttt = stats.t.isf(alpha / 2 if tails == 2 else alpha, dof)
    rc = np.sqrt(ttt ** 2 / (ttt ** 2 + dof))
    zr = np.arctanh(r) + r / (2 * (n - 1))
    zrc = np.arctanh(rc)
    pw = stats.norm.cdf((zr - zrc) * np.sqrt(n - 3))
    if tails == 2:
        pw += stats.norm.cdf((-zr - zrc) * np.sqrt(n - 3))
    return float(pw)


def p_prop(h, n1, ratio, alpha, tails):
    n2 = int(np.ceil(ratio * n1))
    z = stats.norm.isf(alpha / 2 if tails == 2 else alpha)
    d = abs(h) * np.sqrt(n1 * n2 / (n1 + n2))
    return float(stats.norm.sf(z - d) + (stats.norm.cdf(-z - d) if tails == 2 else 0.0))


def p_reg(f2, N, u, alpha):
    df2 = N - u - 1
    return float(stats.ncf.sf(stats.f.isf(alpha, u, df2), u, df2, f2 * N))


def smallest(fn, lo, target):
    n = lo
    while fn(n) < target:
        n += 1
        if n > 200000:
            raise RuntimeError("너무 큼")
    return n


rows = []
# --- 독립표본 t ---
for d in (0.2, 0.5, 0.8, 1.2):
    for ratio in (1, 2):
        for alpha, tails in ((0.05, 2), (0.01, 2), (0.05, 1)):
            for n1 in (10, 30, 64, 200):
                pw = p_t2(d, n1, ratio, alpha, tails)
                sm = TTestIndPower().power(effect_size=d, nobs1=n1, alpha=alpha, ratio=ratio, alternative=alt(tails))
                chk(pw, sm, f"t2 power d={d} n={n1} r={ratio}", 1e-6, 1e-9)
                rows.append(dict(test="ttest2", es=d, n=n1, ratio=ratio, alpha=alpha, tails=tails, power=pw))
            tgt = 0.8
            nreq = smallest(lambda n: TTestIndPower().power(effect_size=d, nobs1=n, alpha=alpha, ratio=ratio, alternative=alt(tails)), 2, tgt)
            rows.append(dict(test="ttest2", kind="nreq", es=d, ratio=ratio, alpha=alpha, tails=tails, target=tgt, nreq=nreq))
            if ratio == 1:
                eff = TTestIndPower().solve_power(effect_size=None, nobs1=40, alpha=alpha, power=0.8, ratio=1, alternative=alt(tails))
                rows.append(dict(test="ttest2", kind="effect", n=40, ratio=1, alpha=alpha, tails=tails, target=0.8, effect=float(eff)))
if True:
    for d, n, a in ((0.5, 34, 0.05), (0.3, 90, 0.05), (0.8, 15, 0.01)):
        chk(p_t2(d, n, 1, a, 2), pg.power_ttest2n(nx=n, ny=n, d=d, alpha=a), "pingouin power_ttest2n", 1e-5, 1e-8)
# --- 대응표본·일표본 t ---
for d in (0.2, 0.5, 0.8):
    for alpha, tails in ((0.05, 2), (0.05, 1), (0.01, 2)):
        for n in (8, 20, 34, 100):
            pw = p_t1(d, n, alpha, tails)
            chk(pw, TTestPower().power(effect_size=d, nobs=n, alpha=alpha, alternative=alt(tails)), f"t1 power d={d} n={n}", 1e-6, 1e-9)
            rows.append(dict(test="ttest1", es=d, n=n, ratio=1, alpha=alpha, tails=tails, power=pw))
        nreq = smallest(lambda n: TTestPower().power(effect_size=d, nobs=n, alpha=alpha, alternative=alt(tails)), 2, 0.8)
        rows.append(dict(test="ttest1", kind="nreq", es=d, ratio=1, alpha=alpha, tails=tails, target=0.8, nreq=nreq))
    eff = TTestPower().solve_power(effect_size=None, nobs=25, alpha=0.05, power=0.8, alternative="two-sided")
    rows.append(dict(test="ttest1", kind="effect", n=25, ratio=1, alpha=0.05, tails=2, target=0.8, effect=float(eff)))
for d, n in ((0.5, 34), (0.3, 80)):
    chk(p_t1(d, n, 0.05, 2), pg.power_ttest(d=d, n=n, alpha=0.05, contrast="paired"), "pingouin power_ttest paired", 1e-5, 1e-8)
# --- 일원 분산분석 ---
for f in (0.1, 0.25, 0.4):
    for k in (3, 4, 5):
        for alpha in (0.05, 0.01):
            for n in (10, 30, 53, 100):
                pw = p_anova(f, n, k, alpha)
                chk(pw, FTestAnovaPower().power(effect_size=f, nobs=k * n, alpha=alpha, k_groups=k), f"anova power f={f} k={k} n={n}", 1e-6, 1e-9)
                rows.append(dict(test="anova", es=f, n=n, k=k, ratio=1, alpha=alpha, tails=2, power=pw))
            nreq = smallest(lambda n: FTestAnovaPower().power(effect_size=f, nobs=k * n, alpha=alpha, k_groups=k), 2, 0.8)
            rows.append(dict(test="anova", kind="nreq", es=f, k=k, ratio=1, alpha=alpha, tails=2, target=0.8, nreq=nreq))
    eff = FTestAnovaPower().solve_power(effect_size=None, nobs=3 * 40, alpha=0.05, power=0.8, k_groups=3)
    rows.append(dict(test="anova", kind="effect", n=40, k=3, ratio=1, alpha=0.05, tails=2, target=0.8, effect=float(eff)))
for f, k, n in ((0.25, 3, 30), (0.4, 4, 15)):
    eta2 = f * f / (1 + f * f)
    chk(p_anova(f, n, k, 0.05), pg.power_anova(eta_squared=eta2, k=k, n=n, alpha=0.05), "pingouin power_anova", 1e-5, 1e-8)
# --- 상관 (Fisher z) ---
for r in (0.1, 0.3, 0.5, 0.7):
    for alpha, tails in ((0.05, 2), (0.05, 1), (0.01, 2)):
        for n in (15, 40, 85, 300):
            pw = p_corr(r, n, alpha, tails)
            rows.append(dict(test="corr", es=r, n=n, ratio=1, alpha=alpha, tails=tails, power=pw))
            pgp = pg.power_corr(r=r, n=n, alpha=alpha, alternative="two-sided" if tails == 2 else "greater")
            chk(pw, pgp, f"pingouin power_corr r={r} n={n} tails={tails}", 1e-5, 1e-8)
        nreq = smallest(lambda n: p_corr(r, n, alpha, tails), 4, 0.8)
        rows.append(dict(test="corr", kind="nreq", es=r, ratio=1, alpha=alpha, tails=tails, target=0.8, nreq=nreq))
    # 효과크기: 이분법으로 직접
    lo, hi = 1e-6, 0.9999
    for _ in range(200):
        mid = (lo + hi) / 2
        if p_corr(mid, 60, 0.05, 2) >= 0.8:
            hi = mid
        else:
            lo = mid
rows.append(dict(test="corr", kind="effect", n=60, ratio=1, alpha=0.05, tails=2, target=0.8, effect=float(hi)))
# --- 두 비율 (Cohen h) ---
for h in (0.2, 0.5, 0.8):
    for ratio in (1, 2):
        for alpha, tails in ((0.05, 2), (0.05, 1)):
            for n1 in (20, 50, 100, 300):
                pw = p_prop(h, n1, ratio, alpha, tails)
                chk(pw, NormalIndPower().power(effect_size=h, nobs1=n1, alpha=alpha, ratio=ratio, alternative=alt(tails)), f"prop power h={h} n={n1}", 1e-6, 1e-9)
                rows.append(dict(test="prop", es=h, n=n1, ratio=ratio, alpha=alpha, tails=tails, power=pw))
            nreq = smallest(lambda n: NormalIndPower().power(effect_size=h, nobs1=n, alpha=alpha, ratio=ratio, alternative=alt(tails)), 2, 0.8)
            rows.append(dict(test="prop", kind="nreq", es=h, ratio=ratio, alpha=alpha, tails=tails, target=0.8, nreq=nreq))
# 두 비율에서 h
for p1, p2 in ((0.5, 0.3), (0.2, 0.1), (0.7, 0.55)):
    h = abs(2 * np.arcsin(np.sqrt(p1)) - 2 * np.arcsin(np.sqrt(p2)))
    rows.append(dict(test="prop", kind="h", p1=p1, p2=p2, h=float(h)))
# --- 선형 회귀 전체 R² ---
for f2 in (0.02, 0.15, 0.35):
    for u in (1, 3, 6):
        for alpha in (0.05, 0.01):
            for N in (30, 77, 150, 400):
                if N < u + 3:
                    continue
                pw = p_reg(f2, N, u, alpha)
                rows.append(dict(test="reg", es=f2, n=N, u=u, ratio=1, alpha=alpha, tails=2, power=pw))
                # statsmodels 의 FTestPower 는 df_num 과 df_denom 의 뜻이 서로 바뀌어 있다 (내부 ftest_power(effect_size, df2, df1) 의 인자 순서)
                chk(pw, FTestPower().power(effect_size=np.sqrt(f2), df_num=N - u - 1, df_denom=u, alpha=alpha), f"reg power f2={f2} u={u} N={N}", 1e-6, 1e-9)
            nreq = smallest(lambda n: p_reg(f2, n, u, alpha), u + 2, 0.8)
            rows.append(dict(test="reg", kind="nreq", es=f2, u=u, ratio=1, alpha=alpha, tails=2, target=0.8, nreq=nreq))
# --- 비중심 F 누적확률 격자 ---
grid = []
for d1, d2, lam in ((2, 57, 4.0), (3, 120, 12.5), (5, 40, 0.5), (1, 20, 30.0), (4, 400, 200.0), (6, 1000, 600.0), (2, 8, 1e-3), (10, 100, 80.0)):
    for q in (0.1, 0.5, 0.9, 0.99):
        x = float(stats.ncf.ppf(q, d1, d2, lam))
        grid.append(dict(x=x, d1=d1, d2=d2, lam=lam, cdf=float(stats.ncf.cdf(x, d1, d2, lam))))

# scipy nct 는 비중심도가 매우 큰 경우(검정력이 거의 1)에 NaN 을 돌려주므로 그런 행은 기준값에서 뺀다
n_before = len(rows)
rows = [r for r in rows if not any(isinstance(v, float) and v != v for v in r.values())]
print("NaN 으로 뺀 행:", n_before - len(rows))
path = os.path.join(os.path.dirname(__file__), "reference", "power.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"statsmodels": statsmodels.__version__, "pingouin": pg.__version__, "scipy": scipy.__version__, "rows": rows, "ncf": grid}, fh, ensure_ascii=False, indent=0)
print("wrote", len(rows), "power rows,", len(grid), "ncf points (scipy = statsmodels = pingouin 일치 확인 통과)")
