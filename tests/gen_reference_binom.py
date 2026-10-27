"""이항 검정 기준값을 tests/reference/binom.json 에 저장한다 (개발용; py -3.12, scipy·statsmodels 필요).

실행: py -3.12 tests/gen_reference_binom.py

- 유의확률(양측·단측): scipy.stats.binomtest (양측은 R binom.test 와 같은 "관측값보다 확률이 같거나 작은 결과의 합")
- 신뢰구간: scipy 의 proportion_ci (exact = Clopper–Pearson, wilson) 와 statsmodels proportion_confint(beta, wilson)가
  서로 일치하는 것을 이 스크립트에서 확인한다.
"""
import json
import os

import numpy as np
from scipy import stats
from statsmodels.stats.proportion import proportion_confint

cases_in = [
    (27, 40, 0.5, 95), (0, 20, 0.3, 95), (20, 20, 0.9, 95), (3, 10, 0.5, 95), (5, 10, 0.5, 95),
    (18, 200, 0.15, 95), (45, 48, 0.8, 99), (3000, 5000, 0.58, 95), (7, 50, 0.1, 90), (1, 3, 0.5, 95),
    (12, 15, 0.5, 99), (2, 40, 0.01, 95), (49, 50, 0.99, 95), (6, 11, 0.25, 90), (30, 30, 0.5, 95),
    (24999, 50000, 0.5, 95), (130, 400, 0.3, 95), (1, 1000, 0.002, 99), (8, 8, 0.7, 95), (0, 5, 0.5, 90),
]

cases = []
for k, n, p, conf in cases_in:
    bt = stats.binomtest(k, n, p)
    ex = bt.proportion_ci(confidence_level=conf / 100, method='exact')
    wi = bt.proportion_ci(confidence_level=conf / 100, method='wilson')
    alpha = 1 - conf / 100
    sb = proportion_confint(k, n, alpha=alpha, method='beta')
    sw = proportion_confint(k, n, alpha=alpha, method='wilson')
    assert abs(ex.low - sb[0]) < 1e-9 and abs(ex.high - sb[1]) < 1e-9, (k, n, ex, sb)
    assert abs(wi.low - sw[0]) < 1e-9 and abs(wi.high - sw[1]) < 1e-9, (k, n, wi, sw)
    cases.append({
        'k': k, 'n': n, 'p0': p, 'conf': conf,
        'two': float(bt.pvalue),
        'greater': float(stats.binomtest(k, n, p, alternative='greater').pvalue),
        'less': float(stats.binomtest(k, n, p, alternative='less').pvalue),
        'exact': [float(ex.low), float(ex.high)], 'wilson': [float(wi.low), float(wi.high)],
    })

path = os.path.join(os.path.dirname(__file__), 'reference', 'binom.json')
with open(path, 'w', encoding='utf-8') as fh:
    json.dump({'cases': cases}, fh, ensure_ascii=False)
print('wrote', len(cases), 'cases')
