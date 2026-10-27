"""일원 분산분석과 사후검정의 기준값을 tests/reference/anova.json 에 저장한다.

실행: py -3.12 tests/gen_reference_anova.py

기준값의 출처
- F, p: scipy.stats.f_oneway
- Tukey HSD: scipy.stats.tukey_hsd (평균차이, 유의확률, 신뢰구간)
- LSD, Šidák: 같은 합동분산 t 통계량에 scipy t 분포로 계산한다 (LSD 는 보정 없음, Šidák 은 1−(1−p)^m 과 α' = 1−(1−α)^(1/m)).
- Dunnett(통제집단 대비 양측): scipy.stats.dunnett 은 준-몬테카를로(qmvt)라 유의확률이 약 1e-3 오차이므로 보조 대조로만 쓰고,
  기준값은 다변량 t 확률을 합동 표준편차 u 와 통제집단 평균 z0 에 대한 이중 적분(scipy.integrate.quad + 300점 가우스-르장드르)으로
  정밀하게 계산한 값과 brentq 로 푼 임계값이다 (구현과 다른 적분 방식).
- Welch, Brown-Forsythe, Games-Howell, Bonferroni, Scheffé: scipy에 없어서 이 파일에서 공식대로 구현하고
  분포 함수만 scipy(f, t, studentized_range)를 쓴다. k = 2 일 때 Welch·Brown-Forsythe 통계량이
  scipy의 Welch t검정(ttest_ind, equal_var=False)의 t² 및 자유도와 같아야 한다는 항등식으로 교차 확인한다.
평균차이는 집단 i − 집단 j (i < j), 사후검정은 모든 쌍을 한 번씩.
"""
import itertools
import json
import os

import numpy as np
import scipy
from scipy import stats

rng = np.random.default_rng(20261023)
f = lambda v: float(v)


def welch_anova(groups):
    k = len(groups)
    n = np.array([len(g) for g in groups], float)
    m = np.array([np.mean(g) for g in groups])
    v = np.array([np.var(g, ddof=1) for g in groups])
    w = n / v
    W = w.sum()
    wm = (w * m).sum() / W
    lam = (((1 - w / W) ** 2) / (n - 1)).sum()
    F = ((w * (m - wm) ** 2).sum() / (k - 1)) / (1 + 2 * (k - 2) / (k * k - 1) * lam)
    df2 = (k * k - 1) / (3 * lam)
    return F, k - 1, df2, stats.f.sf(F, k - 1, df2)


def brown_forsythe(groups):
    k = len(groups)
    n = np.array([len(g) for g in groups], float)
    N = n.sum()
    m = np.array([np.mean(g) for g in groups])
    v = np.array([np.var(g, ddof=1) for g in groups])
    grand = (n * m).sum() / N
    den = ((1 - n / N) * v).sum()
    F = (n * (m - grand) ** 2).sum() / den
    c = (1 - n / N) * v / den
    df2 = 1 / ((c ** 2) / (n - 1)).sum()
    return F, k - 1, df2, stats.f.sf(F, k - 1, df2)


def posthoc(groups, conf):
    k = len(groups)
    n = np.array([len(g) for g in groups], float)
    N = n.sum()
    m = np.array([np.mean(g) for g in groups])
    v = np.array([np.var(g, ddof=1) for g in groups])
    dfw = N - k
    msw = sum((len(g) - 1) * np.var(g, ddof=1) for g in groups) / dfw
    alpha = 1 - conf / 100
    pairs = k * (k - 1) / 2

    tk = stats.tukey_hsd(*groups)
    tk_ci = tk.confidence_interval(conf / 100)

    out = {"tukey": [], "gameshowell": [], "bonferroni": [], "scheffe": [], "lsd": [], "sidak": []}
    for i, j in itertools.combinations(range(k), 2):
        diff = m[i] - m[j]
        se = np.sqrt(msw * (1 / n[i] + 1 / n[j]))
        # Tukey (scipy)
        out["tukey"].append({"i": i, "j": j, "diff": f(tk.statistic[i, j]), "se": f(se), "p": f(tk.pvalue[i, j]),
                             "lo": f(tk_ci.low[i, j]), "hi": f(tk_ci.high[i, j])})
        # Games-Howell
        se_gh = np.sqrt(v[i] / n[i] + v[j] / n[j])
        df_gh = (v[i] / n[i] + v[j] / n[j]) ** 2 / ((v[i] / n[i]) ** 2 / (n[i] - 1) + (v[j] / n[j]) ** 2 / (n[j] - 1))
        q = abs(diff) * np.sqrt(2) / se_gh
        hw = stats.studentized_range.isf(alpha, k, df_gh) * se_gh / np.sqrt(2)
        out["gameshowell"].append({"i": i, "j": j, "diff": f(diff), "se": f(se_gh),
                                   "p": f(stats.studentized_range.sf(q, k, df_gh)),
                                   "lo": f(diff - hw), "hi": f(diff + hw)})
        # Bonferroni
        t = diff / se
        tc = stats.t.isf(alpha / (2 * pairs), dfw)
        out["bonferroni"].append({"i": i, "j": j, "diff": f(diff), "se": f(se),
                                  "p": f(min(1.0, pairs * 2 * stats.t.sf(abs(t), dfw))),
                                  "lo": f(diff - tc * se), "hi": f(diff + tc * se)})
        # LSD (보정 없음)와 Šidák
        tc_lsd = stats.t.isf(alpha / 2, dfw)
        out["lsd"].append({"i": i, "j": j, "diff": f(diff), "se": f(se), "p": f(2 * stats.t.sf(abs(t), dfw)),
                          "lo": f(diff - tc_lsd * se), "hi": f(diff + tc_lsd * se)})
        tc_sid = stats.t.isf((1 - (1 - alpha) ** (1 / pairs)) / 2, dfw)
        out["sidak"].append({"i": i, "j": j, "diff": f(diff), "se": f(se), "p": f(1 - (1 - 2 * stats.t.sf(abs(t), dfw)) ** pairs),
                            "lo": f(diff - tc_sid * se), "hi": f(diff + tc_sid * se)})
        # Scheffé
        crit = np.sqrt((k - 1) * stats.f.isf(alpha, k - 1, dfw))
        out["scheffe"].append({"i": i, "j": j, "diff": f(diff), "se": f(se),
                               "p": f(stats.f.sf(t * t / (k - 1), k - 1, dfw)),
                               "lo": f(diff - crit * se), "hi": f(diff + crit * se)})
    return out


from scipy import integrate
from scipy.stats import norm, chi2


def dunnett_prob(q, a, b, df):
    """P(모든 |T_i| < q): T_i = (a_i Z_i - b_i z0) / u, u = sqrt(chi2_df / df). 정밀 이중 적분."""
    zs, zw = np.polynomial.legendre.leggauss(300)
    z = 9 * zs
    w = 9 * zw * norm.pdf(z)

    def inner(u):
        hi = (q * u + b[:, None] * z[None, :]) / a[:, None]
        lo = (-q * u + b[:, None] * z[None, :]) / a[:, None]
        return (w * np.prod(norm.cdf(hi) - norm.cdf(lo), axis=0)).sum()

    fu = lambda u: chi2.pdf(df * u * u, df) * 2 * df * u
    u0 = np.sqrt(chi2.ppf(1e-14, df) / df)
    u1 = np.sqrt(chi2.isf(1e-14, df) / df)
    val, _ = integrate.quad(lambda u: inner(u) * fu(u), u0, u1, epsabs=1e-13, epsrel=1e-12, limit=400)
    return val


def dunnett(groups, ctrl, conf):
    from scipy.optimize import brentq
    k = len(groups)
    n = np.array([len(g) for g in groups], float)
    m = np.array([np.mean(g) for g in groups])
    dfw = n.sum() - k
    msw = sum((len(g) - 1) * np.var(g, ddof=1) for g in groups) / dfw
    others = [i for i in range(k) if i != ctrl]
    a = np.array([np.sqrt(n[ctrl] / (n[i] + n[ctrl])) for i in others])
    b = np.array([np.sqrt(n[i] / (n[i] + n[ctrl])) for i in others])
    alpha = 1 - conf / 100
    qd = brentq(lambda q: dunnett_prob(q, a, b, dfw) - (1 - alpha), stats.t.isf(alpha / 2, dfw) * 0.999, stats.t.isf(alpha / (2 * len(others)), dfw) * 1.001, xtol=1e-13, rtol=1e-13)
    rows = []
    for i in others:
        diff = m[i] - m[ctrl]
        se = np.sqrt(msw * (1 / n[i] + 1 / n[ctrl]))
        p = 1 - dunnett_prob(abs(diff / se), a, b, dfw)
        rows.append({"i": i, "j": ctrl, "diff": f(diff), "se": f(se), "p": f(p), "lo": f(diff - qd * se), "hi": f(diff + qd * se)})
    # scipy 의 준-몬테카를로 값은 보조 대조 (오차 약 1e-3)
    others_samples = [groups[i] for i in others]
    sc = stats.dunnett(*others_samples, control=groups[ctrl], rng=np.random.default_rng(1))
    sci = sc.confidence_interval(conf / 100)
    for r, p, lo, hi in zip(rows, sc.pvalue, sci.low, sci.high):
        r["scipy_p"] = f(p)
        r["scipy_lo"] = f(lo)
        r["scipy_hi"] = f(hi)
    return {"q": f(qd), "rows": rows}


sets = {
    "정규 3집단 (크기·분산 다름)": [rng.normal(10, 1, 12), rng.normal(12, 3, 25), rng.normal(11, 2, 18)],
    "정규 4집단 (균형)": [rng.normal(50 + 2 * i, 5, 15) for i in range(4)],
    "로그정규 3집단": [rng.lognormal(0, s, n) for s, n in [(0.4, 14), (0.8, 20), (1.2, 9)]],
    "정수 점수 5집단": [rng.integers(1, 8, n).astype(float) + i for i, n in enumerate([10, 14, 9, 17, 12])],
    "2집단": [rng.normal(5, 1, 20), rng.normal(6, 2.5, 13)],
    "작은 표본 3집단 (3~4건, 오차 자유도 7)": [rng.normal(10, 2, 3), rng.normal(13, 2, 4), rng.normal(11, 2, 3)],
    "불균형 7집단 (크기 4~30)": [rng.normal(20 + i, 4, n) for i, n in enumerate([4, 30, 9, 15, 6, 22, 11])],
}

cases = []
for name, groups in sets.items():
    groups = [np.asarray(g, float) for g in groups]
    F, p = stats.f_oneway(*groups)
    k = len(groups)
    allv = np.concatenate(groups)
    grand = allv.mean()
    ssb = sum(len(g) * (g.mean() - grand) ** 2 for g in groups)
    ssw = sum(((g - g.mean()) ** 2).sum() for g in groups)
    wl = welch_anova(groups)
    bf = brown_forsythe(groups)
    if k == 2:  # 항등식 교차 확인
        t = stats.ttest_ind(groups[0], groups[1], equal_var=False)
        for label, got in (("welch", wl), ("brown-forsythe", bf)):
            assert abs(got[0] - t.statistic ** 2) < 1e-9 * got[0], (name, label, got[0], t.statistic ** 2)
            assert abs(got[2] - t.df) < 1e-9 * t.df, (name, label, got[2], t.df)
            assert abs(got[3] - t.pvalue) < 1e-9, (name, label, got[3], t.pvalue)
    msw = ssw / (len(allv) - k)
    cases.append({
        "name": name, "groups": [[f(v) for v in g] for g in groups],
        "F": f(F), "p": f(p), "ssb": f(ssb), "ssw": f(ssw), "dfb": k - 1, "dfw": len(allv) - k,
        "eta2": f(ssb / (ssb + ssw)), "omega2": f((ssb - (k - 1) * msw) / (ssb + ssw + msw)),
        "welch": {"F": f(wl[0]), "df1": wl[1], "df2": f(wl[2]), "p": f(wl[3])},
        "bf": {"F": f(bf[0]), "df1": bf[1], "df2": f(bf[2]), "p": f(bf[3])},
        "posthoc": {str(conf): posthoc(groups, conf) for conf in (95, 99)},
        "dunnett": {str(conf): {str(ctrl): dunnett(groups, ctrl, conf) for ctrl in sorted({0, k // 2, k - 1})} for conf in (95, 99)},
    })

path = os.path.join(os.path.dirname(__file__), "reference", "anova.json")
os.makedirs(os.path.dirname(path), exist_ok=True)
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "cases": cases}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "ANOVA cases to", path)
