"""scipy로 확률분포 기준값을 만들어 tests/reference/distributions.json 에 저장한다.

실행: py -3.12 tests/gen_reference.py
기준값 파일은 저장소에 커밋하므로, scipy가 없는 환경에서도 Node 테스트가 돈다.
"""
import json
import math
import os

import scipy
from scipy import stats

cases = []


skipped = []


def add(fn, args, lower, expected):
    if expected is None or not math.isfinite(expected):
        return
    cases.append({"fn": fn, "args": args, "lower": lower, "expected": float(expected)})


def addq(fn, dist, p, extra, lower):
    """역함수 기준값. scipy 자신의 역산(cdf/sf로 되돌린 값)이 p와 1e-10 안에서 맞지 않으면
    scipy의 역함수가 부정확한 사례이므로 기준값에서 제외한다. (예: f.isf(1e-12, 1, 1))"""
    x = dist.ppf(p, *extra) if lower else dist.isf(p, *extra)
    if not math.isfinite(x) or x == 0:
        return
    back = dist.cdf(x, *extra) if lower else dist.sf(x, *extra)
    if not abs(back - p) <= 1e-10 * p:
        skipped.append((fn, [p, *extra], lower))
        return
    add(fn, [p, *extra], lower, x)


# p 함수는 아래꼬리(cdf)와 위꼬리(sf)를 모두 기록한다.
# q 함수는 ppf(아래꼬리 확률)와 isf(위꼬리 확률)를 모두 기록한다.
P_GRID = [1e-300, 1e-100, 1e-20, 1e-12, 1e-8, 1e-5, 0.001, 0.01, 0.025, 0.05, 0.1, 0.3, 0.5,
          0.7, 0.9, 0.95, 0.975, 0.99, 0.999, 0.99999]

# 정규
for x in [-38, -30, -10, -5, -3, -1.96, -1, -0.1, 0, 0.5, 1, 1.96, 2, 3, 5, 8, 10, 20, 37]:
    add("pnorm", [x], True, stats.norm.cdf(x))
    add("pnorm", [x], False, stats.norm.sf(x))
for p in P_GRID:
    addq("qnorm", stats.norm, p, [], True)
    addq("qnorm", stats.norm, p, [], False)

# t
T_DF = [0.5, 1, 2, 3, 5, 10, 30, 100, 1000, 1e5, 1e7]
for df in T_DF:
    for t in [-50, -10, -5, -2, -1, -0.1, -1e-4, 0, 1e-4, 0.1, 0.5, 1, 2, 3, 5, 10, 50]:
        add("pt", [t, df], True, stats.t.cdf(t, df))
        add("pt", [t, df], False, stats.t.sf(t, df))
    for p in P_GRID:
        addq("qt", stats.t, p, [df], True)
        addq("qt", stats.t, p, [df], False)

# χ²
C_DF = [0.1, 0.5, 1, 2, 3, 5, 10, 50, 100, 1000, 1e5]
for df in C_DF:
    for m in [1e-5, 0.01, 0.1, 0.5, 1, 2, 5, 20]:
        x = df * m
        add("pchisq", [x, df], True, stats.chi2.cdf(x, df))
        add("pchisq", [x, df], False, stats.chi2.sf(x, df))
    for p in P_GRID:
        addq("qchisq", stats.chi2, p, [df], True)
        addq("qchisq", stats.chi2, p, [df], False)

# F
F_DF = [(1, 1), (1, 10), (2, 5), (3, 30), (5, 100), (10, 10), (30, 5), (100, 1000), (1, 1e5), (1000, 1000), (0.5, 3)]
for d1, d2 in F_DF:
    for f in [0.001, 0.01, 0.1, 0.5, 1, 2, 5, 20, 100, 1000]:
        add("pf", [f, d1, d2], True, stats.f.cdf(f, d1, d2))
        add("pf", [f, d1, d2], False, stats.f.sf(f, d1, d2))
    for p in P_GRID:
        addq("qf", stats.f, p, [d1, d2], True)
        addq("qf", stats.f, p, [d1, d2], False)

out = os.path.join(os.path.dirname(__file__), "reference", "distributions.json")
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "cases": cases}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "cases to", out)
print("scipy 역산 불일치로 제외한 역함수 사례:", len(skipped))
for item in skipped:
    print("  ", item)
