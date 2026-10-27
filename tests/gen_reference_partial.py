"""편상관 기준값을 tests/reference/partial.json 에 저장한다.

실행: py -3.12 tests/gen_reference_partial.py   (pingouin 이 개발용으로 설치돼 있어야 한다)
기준값은 세 갈래로 얻고 서로 일치하는지 이 스크립트가 확인한다.
 1) numpy 로 정의대로: 통제변수에 회귀한 잔차끼리의 상관 (lstsq)
 2) numpy 로 다른 경로: 변수+통제변수 상관행렬의 역행렬에서 −P_ij/√(P_ii P_jj) (통제변수만 있는 부분행렬을 쓰는 방식이 아니라 전체 변수 중 두 변수를 뺀 나머지를 통제하는 식)
 3) pingouin partial_corr (피어슨·스피어만): r, 신뢰구간, 유의확률, n
스피어만은 모든 변수를 순위(동순위 평균)로 바꾼 뒤 같은 계산이다 (pingouin 과 같은 방식).
저장하는 값: 통제 전 상관 행렬, 편상관 r·p·95% 신뢰구간·자유도 (변수 쌍마다, 피어슨·스피어만)
"""
import json
import os
import warnings
from itertools import combinations

import numpy as np
import pandas as pd
import pingouin as pg
import scipy
from scipy import stats

warnings.filterwarnings("ignore")
pg.options["round"] = None
pg.options["round.column.CI95"] = None
rng = np.random.default_rng(20261103)
chk = lambda x, y, what, tol=1e-9, atol=1e-12: np.testing.assert_allclose(np.asarray(x, float), np.asarray(y, float), rtol=tol, atol=atol, err_msg=what)


def rank(a):
    return stats.rankdata(a)


def analyse(name, V, C, vnames, cnames, conf=95):
    V = np.asarray(V, float)
    C = np.asarray(C, float)
    n, k = V.shape
    c = C.shape[1]
    out = dict(name=name, n=n, k=k, c=c, vnames=vnames, cnames=cnames, V=V.T.tolist(), C=C.T.tolist(), conf=conf)
    for method in ("pearson", "spearman"):
        Vm = np.apply_along_axis(rank, 0, V) if method == "spearman" else V
        Cm = np.apply_along_axis(rank, 0, C) if method == "spearman" else C
        zero = np.corrcoef(Vm, rowvar=False)
        pr = np.full((k, k), np.nan)
        pv = np.full((k, k), np.nan)
        lo = np.full((k, k), np.nan)
        hi = np.full((k, k), np.nan)
        df = n - 2 - c
        zc = stats.norm.ppf(1 - (1 - conf / 100) / 2)
        Zc = np.column_stack([np.ones(n), Cm])
        # 방법 1: 잔차의 상관
        res = Vm - Zc @ np.linalg.lstsq(Zc, Vm, rcond=None)[0]
        pr1 = np.corrcoef(res, rowvar=False)
        for i, j in combinations(range(k), 2):
            # 방법 2: 두 변수 + 통제변수의 역상관행렬
            M = np.column_stack([Vm[:, i], Vm[:, j], Cm])
            P = np.linalg.inv(np.corrcoef(M, rowvar=False))
            r2 = -P[0, 1] / np.sqrt(P[0, 0] * P[1, 1])
            chk(r2, pr1[i, j], f"{name} {method} 두 방법 ({i},{j})", 1e-8)
            r = float(pr1[i, j])
            t = r * np.sqrt(df / (1 - r * r))
            p = float(2 * stats.t.sf(abs(t), df))
            z = np.arctanh(r)
            se = 1 / np.sqrt(n - c - 3)
            pr[i, j] = pr[j, i] = r
            pv[i, j] = pv[j, i] = p
            lo[i, j] = lo[j, i] = np.tanh(z - zc * se)
            hi[i, j] = hi[j, i] = np.tanh(z + zc * se)
            # pingouin
            data = pd.DataFrame(np.column_stack([V[:, i], V[:, j], C]), columns=["x", "y"] + [f"z{q}" for q in range(c)])
            pgr = pg.partial_corr(data=data, x="x", y="y", covar=[f"z{q}" for q in range(c)], method=method).iloc[0]
            chk(pgr["r"], r, f"{name} {method} pingouin r ({i},{j})", 1e-8)
            chk(pgr["p_val"], p, f"{name} {method} pingouin p ({i},{j})", 1e-6, 1e-14)
            if abs(conf - 95) < 1e-12:
                chk(pgr["CI95"][0], lo[i, j], f"{name} {method} pingouin 하한", 1e-5, 2e-6)
                chk(pgr["CI95"][1], hi[i, j], f"{name} {method} pingouin 상한", 1e-5, 2e-6)
        np.fill_diagonal(pr, 1.0)
        out[method] = dict(zero=zero.tolist(), r=pr.tolist(), p=pv.tolist(), lo=lo.tolist(), hi=hi.tolist(), df=df)
    return out


def sim(n, k, c, load_c, extra=0.0):
    C = rng.normal(0, 1, (n, c))
    V = np.column_stack([C @ np.array(load_c[j]) + rng.normal(0, 1, n) + (extra * rng.normal(0, 1, n) if j == 0 else 0) for j in range(k)])
    return np.round(V, 2), np.round(C, 2)


cases = []
V, C = sim(60, 3, 1, [[0.8], [0.7], [-0.4]])
cases.append(analyse("변수 3개, 통제변수 1개", V, C, ["x", "y", "w"], ["z"]))
V, C = sim(80, 4, 2, [[0.7, 0.2], [0.5, 0.5], [0.1, 0.8], [-0.3, 0.3]])
cases.append(analyse("변수 4개, 통제변수 2개", V, C, ["a", "b", "c", "d"], ["z1", "z2"]))
V, C = sim(50, 3, 3, [[0.6, 0.3, 0.1], [0.2, 0.6, 0.4], [0.5, 0.5, 0.5]])
cases.append(analyse("변수 3개, 통제변수 3개", V, C, ["x", "y", "w"], ["z1", "z2", "z3"]))
# 스피어만, 동순위가 많은 정수 점수
Cz = rng.integers(1, 6, 70).astype(float)[:, None]
V = np.column_stack([np.clip(np.round(Cz[:, 0] * 0.6 + rng.normal(0, 1.2, 70) + 1), 1, 5), np.clip(np.round(Cz[:, 0] * 0.5 + rng.normal(0, 1.1, 70) + 1), 1, 5), np.clip(np.round(rng.normal(3, 1.2, 70)), 1, 5)])
cases.append(analyse("5점 척도 정수 점수 3개, 통제변수 1개 (동순위 많음)", V, Cz, ["q1", "q2", "q3"], ["age"]))
# 억제: 통제하면 상관의 부호가 바뀐다
z = rng.normal(0, 1, 100)
x = np.round(z + rng.normal(0, 0.5, 100), 2)
y = np.round(-0.6 * x + 1.4 * z + rng.normal(0, 0.6, 100), 2)
cases.append(analyse("통제 후 상관의 부호가 바뀌는 억제 사례", np.column_stack([x, y]), np.round(z, 2)[:, None], ["x", "y"], ["z"]))
# 작은 표본, 90% 신뢰구간 (pingouin 대조는 95%만)
V, C = sim(12, 2, 2, [[0.6, 0.2], [0.3, 0.6]])
cases.append(analyse("작은 표본 12명, 변수 2개, 통제변수 2개, 90% 구간", V, C, ["x", "y"], ["z1", "z2"], conf=90))

def scrub(o):
    """NaN 은 JSON 에 쓸 수 없어 null 로 바꾼다 (대각선의 유의확률·신뢰구간)."""
    if isinstance(o, dict):
        return {k: scrub(v) for k, v in o.items()}
    if isinstance(o, list):
        return [scrub(v) for v in o]
    if isinstance(o, float) and o != o:
        return None
    return o


cases = scrub(cases)
path = os.path.join(os.path.dirname(__file__), "reference", "partial.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"pingouin": pg.__version__, "scipy": scipy.__version__, "cases": cases}, fh, ensure_ascii=False, indent=0, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("wrote", len(cases), "partial cases (numpy 두 경로 = pingouin 일치 확인 통과);",
      [(c["n"], c["c"], round(c["pearson"]["zero"][0][1], 3), round(c["pearson"]["r"][0][1], 3), round(c["spearman"]["r"][0][1], 3)) for c in cases])
