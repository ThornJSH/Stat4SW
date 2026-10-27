"""McNemar 검정 기준값을 tests/reference/mcnemar.json 에 저장한다.

실행: py -3.12 tests/gen_reference_mcnemar.py
- 정확 유의확률(양측): scipy.stats.binomtest(min(b, c), b + c, 0.5).pvalue
- 카이제곱(연속성 수정) (|b − c| − 1)² / (b + c), 자유도 1, p = chi2.sf
- 양성(큰 값) 비율의 변화 (b − c)/n 와 Wald 신뢰구간 ± z·√(b + c − (b − c)²/n)/n
- 오즈비 b/c 의 정확 신뢰구간: 불일치 쌍 중 b 의 비율 π 에 Clopper–Pearson (scipy.stats.beta.ppf) 후 π/(1 − π)
표 모양: a = (전 작은 값, 후 작은 값), b = (전 작은 값, 후 큰 값), c = (전 큰 값, 후 작은 값), d = (전 큰 값, 후 큰 값)
"""
import json
import os

import numpy as np
import scipy
from scipy import stats

rng = np.random.default_rng(20261035)


def case(name, a, b, c, d, lo_code=0, hi_code=1, conf=95):
    n = a + b + c + d
    nd = b + c
    before, after = [], []
    for cnt, (x, y) in [(a, (lo_code, lo_code)), (b, (lo_code, hi_code)), (c, (hi_code, lo_code)), (d, (hi_code, hi_code))]:
        before += [x] * cnt
        after += [y] * cnt
    order = rng.permutation(n)
    before = [before[i] for i in order]
    after = [after[i] for i in order]
    exact = float(stats.binomtest(min(b, c), nd, 0.5).pvalue) if nd > 0 else 1.0
    chi2 = (abs(b - c) - 1) ** 2 / nd if nd > 0 else None
    chi_p = float(stats.chi2.sf(chi2, 1)) if nd > 0 else None
    diff = (b - c) / n
    z = float(stats.norm.ppf(1 - (1 - conf / 100) / 2))
    se = float(np.sqrt(nd - (b - c) ** 2 / n) / n)
    pi_lo = float(stats.beta.ppf((1 - conf / 100) / 2, b, c + 1)) if b > 0 else 0.0
    pi_hi = float(stats.beta.ppf(1 - (1 - conf / 100) / 2, b + 1, c)) if c > 0 else 1.0
    odds = lambda p: None if p >= 1 else p / (1 - p)
    return dict(name=name, a=a, b=b, c=c, d=d, n=n, nd=nd, lo_code=lo_code, hi_code=hi_code, before=before, after=after,
                exact=exact, chi2=chi2, chi_p=chi_p,
                p_before=(c + d) / n, p_after=(b + d) / n, diff=diff, diff_lo=diff - z * se, diff_hi=diff + z * se,
                odds=(b / c if c > 0 else None), odds_lo=odds(pi_lo), odds_hi=odds(pi_hi))


cases = [
    case("보통 크기 (불일치 21쌍)", 40, 6, 15, 39),
    case("불일치가 매우 적음 (7쌍)", 30, 1, 6, 23),
    case("한쪽 방향만 변함 (c = 0)", 25, 9, 0, 16),
    case("불일치가 많음 (100쌍, 점근 근사 가능)", 120, 70, 30, 80),
    case("값이 1·2로 코딩됨, b = c", 18, 8, 8, 26, lo_code=1, hi_code=2),
]
path = os.path.join(os.path.dirname(__file__), "reference", "mcnemar.json")
with open(path, "w", encoding="utf-8") as fh:
    json.dump({"scipy": scipy.__version__, "cases": cases}, fh, ensure_ascii=False, indent=0)
print("wrote", len(cases), "mcnemar cases;", [(c["exact"].__round__(5), None if c["chi_p"] is None else round(c["chi_p"], 5)) for c in cases])
