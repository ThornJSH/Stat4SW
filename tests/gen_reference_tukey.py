"""scipy.stats.studentized_range 로 Studentized range 분포의 기준값을 tests/reference/tukey.json 에 저장한다.

실행: py -3.12 tests/gen_reference_tukey.py
scipy 자신의 역함수(ppf/isf)가 cdf/sf로 되돌렸을 때 p와 맞지 않는 사례는 기준값에서 제외한다.
"""
import json
import os

import scipy
from scipy import stats

K = [2, 3, 4, 5, 8, 10, 20]
DF = [2, 3, 5, 10, 20, 60, 120, 1000]
Q = [0.1, 0.5, 1, 2, 3, 3.5, 4, 5, 7, 10]
P = [0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99]

cdf, ppf, skipped = [], [], []
for k in K:
    for df in DF:
        d = stats.studentized_range(k, df)
        for q in Q:
            cdf.append({"q": q, "k": k, "df": df, "lower": float(d.cdf(q)), "upper": float(d.sf(q))})
        for p in P:
            x = float(d.ppf(p))
            back = float(d.cdf(x))
            if abs(back - p) > 1e-9:
                skipped.append((k, df, p))
                continue
            ppf.append({"p": p, "k": k, "df": df, "q": x})

path = os.path.join(os.path.dirname(__file__), "reference", "tukey.json")
os.makedirs(os.path.dirname(path), exist_ok=True)
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "cdf": cdf, "ppf": ppf}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cdf), "cdf and", len(ppf), "ppf cases; skipped (scipy 역산 불일치):", skipped)
