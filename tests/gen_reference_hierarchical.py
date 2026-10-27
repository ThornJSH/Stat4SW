"""위계적(블록) 회귀 기준값을 tests/reference/hierarchical.json 에 저장한다.

실행: py -3.12 tests/gen_reference_hierarchical.py   (statsmodels 가 개발용으로 설치돼 있어야 한다)
기준값은 두 갈래로 얻고 서로 일치하는지 이 스크립트가 확인한다.
 1) numpy 로 정의대로 직접 계산 (lstsq, (X'X)^-1, 다른 변수로의 회귀로 VIF, 잔차 제곱합 차이로 변화 F)
 2) statsmodels (OLS 의 계수·표준오차·R²·수정 R², anova_lm 의 중첩모형 F 검정, durbin_watson)
저장하는 값: 단계별 R, R², 수정 R², 추정값의 표준오차, 변화량(ΔR², ΔF, 자유도, 유의확률), 모형 F·p, 제곱합,
계수(B, 표준오차, t, p, 신뢰구간, 베타, 공차, VIF), 더빈-왓슨
범주형은 첫 범주를 기준으로 한 더미로 바꿔 계산한다 (도구와 같은 규칙).
"""
import json
import os
import warnings

import numpy as np
import scipy
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.stattools import durbin_watson

warnings.filterwarnings("ignore")
rng = np.random.default_rng(20261026)


def dummies(vals):
    cats = sorted(set(vals))
    return [(f"= {c}", [1.0 if v == c else 0.0 for v in vals]) for c in cats[1:]]


def fit_numpy(cols, y, conf):
    n = len(y)
    X = np.column_stack([np.ones(n)] + cols)
    p = X.shape[1]
    B, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ B
    sse = float(resid @ resid)
    sst = float(((y - y.mean()) ** 2).sum())
    ssr = sst - sse
    dfm, dfr = p - 1, n - p
    mse = sse / dfr
    inv = np.linalg.inv(X.T @ X)
    se = np.sqrt(mse * np.diag(inv))
    t = B / se
    pv = 2 * stats.t.sf(np.abs(t), dfr)
    tc = stats.t.ppf(1 - (1 - conf / 100) / 2, dfr)
    sdy = y.std(ddof=1)
    beta, tol, vif = [None], [None], [None]
    for j, c in enumerate(cols):
        beta.append(float(B[j + 1] * c.std(ddof=1) / sdy))
        others = [np.ones(n)] + [cc for q, cc in enumerate(cols) if q != j]
        Z = np.column_stack(others)
        g, *_ = np.linalg.lstsq(Z, c, rcond=None)
        r = c - Z @ g
        r2j = 1 - float(r @ r) / float(((c - c.mean()) ** 2).sum())
        vif.append(1 / (1 - r2j))
        tol.append(1 - r2j)
    F = (ssr / dfm) / mse
    return dict(n=n, p=p, B=B.tolist(), se=se.tolist(), t=t.tolist(), pval=pv.tolist(),
                lo=(B - tc * se).tolist(), hi=(B + tc * se).tolist(), beta=beta, tol=tol, vif=vif,
                sse=sse, ssr=ssr, sst=sst, dfm=dfm, dfr=dfr, mse=mse, s=float(np.sqrt(mse)),
                r2=ssr / sst, r=float(np.sqrt(ssr / sst)), adj=1 - (1 - ssr / sst) * (n - 1) / dfr,
                F=float(F), pF=float(stats.f.sf(F, dfm, dfr)),
                dw=float(((resid[1:] - resid[:-1]) ** 2).sum() / sse)), X


def analyse(name, y, blocks, conf=95):
    """blocks: [[(변수이름, 종류, 값 목록)]], 종류는 'num' 또는 'factor'."""
    y = np.asarray(y, float)
    cols, names, steps, models = [], ["(상수)"], [], []
    prev = None
    for bi, block in enumerate(blocks):
        for vname, kind, vals in block:
            if kind == "factor":
                for suffix, col in dummies(list(vals)):
                    cols.append(np.asarray(col)); names.append(f"{vname} {suffix}")
            else:
                cols.append(np.asarray(vals, float)); names.append(vname)
        f, X = fit_numpy(cols, y, conf)
        f["names"] = list(names)
        # 변화량: 잔차 제곱합 차이로 직접
        if prev is None:
            f.update(dR2=f["r2"], dF=f["F"], df1=f["dfm"], df2=f["dfr"], pDF=f["pF"])
        else:
            m = f["p"] - prev["p"]
            dF = ((prev["sse"] - f["sse"]) / m) / (f["sse"] / f["dfr"])
            f.update(dR2=f["r2"] - prev["r2"], dF=float(dF), df1=m, df2=f["dfr"], pDF=float(stats.f.sf(dF, m, f["dfr"])))
        # statsmodels 교차 확인
        ols = sm.OLS(y, X).fit()
        chk = lambda a, b, what, tol=1e-7: np.testing.assert_allclose(np.asarray(a, float), np.asarray(b, float), rtol=tol, atol=1e-10, err_msg=f"{name} 단계{bi + 1} {what}")
        chk(ols.params, f["B"], "B"); chk(ols.bse, f["se"], "se"); chk(ols.pvalues, f["pval"], "p")
        chk(ols.rsquared, f["r2"], "R2"); chk(ols.rsquared_adj, f["adj"], "adjR2"); chk(ols.fvalue, f["F"], "F")
        chk(durbin_watson(ols.resid), f["dw"], "DW")
        if prev is not None:
            a = anova_lm(models[-1], ols).iloc[1]
            chk(a["F"], f["dF"], "변화 F"); chk(a["Pr(>F)"], f["pDF"], "변화 F 의 p"); chk(a["df_diff"], f["df1"], "df1")
        models.append(ols)
        steps.append(f)
        prev = f
    return dict(name=name, conf=conf, n=len(y), y=y.tolist(),
                blocks=[[dict(name=v, kind=k, values=list(map(float, vals))) for v, k, vals in b] for b in blocks],
                steps=steps)


cases = []

# A: 3단계, 모두 숫자 (60명): 통제변수 → 핵심 변수 → 추가 변수 2개
n = 60
age = np.round(rng.normal(40, 8, n), 1)
ten = np.round(0.5 * (age - 25) + rng.normal(0, 3, n), 1)
x1 = np.round(rng.normal(0, 1, n) + 0.2 * (age - 40) / 8, 2)
x2 = np.round(rng.normal(0, 1, n), 2)
x3 = np.round(rng.normal(0, 1, n) + 0.4 * x1, 2)
y = np.round(50 + 0.3 * age + 0.5 * ten + 4 * x1 + 1.5 * x2 + rng.normal(0, 5, n), 1)
cases.append(analyse("60명, 숫자 변수 3단계", y, [[("age", "num", age), ("ten", "num", ten)], [("x1", "num", x1)], [("x2", "num", x2), ("x3", "num", x3)]]))

# B: 2단계, 1단계에 3범주 + 숫자, 2단계에 숫자 2개 (45명)
n = 45
grp = rng.integers(1, 4, n)
age = np.round(rng.normal(35, 9, n), 1)
z1 = np.round(rng.normal(0, 1, n), 2)
z2 = np.round(rng.normal(0, 1, n) + 0.3 * z1, 2)
y = np.round(20 + 2.0 * (grp == 2) + 4.5 * (grp == 3) + 0.2 * age + 2.5 * z1 + 0.6 * z2 + rng.normal(0, 4, n), 1)
cases.append(analyse("45명, 1단계에 3범주 + 숫자, 2단계 숫자 2개", y, [[("grp", "factor", grp), ("age", "num", age)], [("z1", "num", z1), ("z2", "num", z2)]]))

# C: 4단계, 단계마다 변수 1개 (30명), 마지막 단계는 설명력이 거의 늘지 않음
n = 30
v = [np.round(rng.normal(0, 1, n), 2) for _ in range(4)]
y = np.round(10 + 2.0 * v[0] + 1.2 * v[1] + 0.6 * v[2] + rng.normal(0, 2, n), 2)
cases.append(analyse("30명, 4단계(단계마다 1개), 마지막 단계는 변화가 거의 없음", y, [[(f"v{i + 1}", "num", v[i])] for i in range(4)]))

# D: 2단계, 리커트(1~5) 점수 변수와 99% 신뢰구간 (80명)
n = 80
l1 = rng.integers(1, 6, n).astype(float)
l2 = np.clip(np.round(l1 + rng.normal(0, 1.2, n)), 1, 5)
l3 = rng.integers(1, 6, n).astype(float)
l4 = np.clip(np.round(0.5 * l3 + rng.normal(2, 1, n)), 1, 5)
y = np.round(3 + 0.4 * l1 + 0.5 * l3 + 0.3 * l4 + rng.normal(0, 1, n), 2)
cases.append(analyse("80명, 리커트 점수 변수 2단계, 99% 신뢰구간", y, [[("l1", "num", l1), ("l2", "num", l2)], [("l3", "num", l3), ("l4", "num", l4)]], conf=99))

path = os.path.join(os.path.dirname(__file__), "reference", "hierarchical.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"statsmodels": sm.__version__, "scipy": scipy.__version__, "cases": cases}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "hierarchical cases (numpy = statsmodels 일치 확인 통과);",
      [[(round(s["r2"], 3), round(s["dF"], 3), round(s["pDF"], 4)) for s in c["steps"]] for c in cases])
