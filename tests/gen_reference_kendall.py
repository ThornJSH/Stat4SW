"""켄달 타우-b 기준값을 tests/reference/kendall.json 에 저장한다 (개발용; py -3.12, scipy 필요).

실행: py -3.12 tests/gen_reference_kendall.py

- 타우-b, 유의확률: scipy.stats.kendalltau (variant='b')
  · method='auto' (기본: 동순위가 없고 n ≤ 33 이면 정확, 아니면 점근)
  · method='asymptotic' (동순위 보정 분산의 정규 근사)
  · method='exact' (동순위가 없는 경우만)
- pingouin 의 corr(method='kendall') 도 같은 타우-b 를 낸다는 것을 확인한다.
"""
import json
import os

import numpy as np
import pingouin as pg
from scipy import stats

rng = np.random.default_rng(20261027)


def case(name, x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    auto = stats.kendalltau(x, y, variant='b')
    asym = stats.kendalltau(x, y, variant='b', method='asymptotic')
    out = {'name': name, 'n': n, 'x': x.tolist(), 'y': y.tolist(),
           'tau': float(auto.statistic), 'p_auto': float(auto.pvalue),
           'p_asym': float(asym.pvalue)}
    has_ties = len(np.unique(x)) < n or len(np.unique(y)) < n
    out['ties'] = bool(has_ties)
    if not has_ties:
        out['p_exact'] = float(stats.kendalltau(x, y, variant='b', method='exact').pvalue)
    pgr = pg.corr(x, y, method='kendall')
    assert abs(float(pgr['r'].iloc[0]) - out['tau']) < 1e-9, (name, pgr['r'].iloc[0], out['tau'])
    return out


cases = []
# 동순위 없음, 작은 표본 (정확 검정)
for n in (5, 8, 12, 20, 33):
    z = rng.normal(0, 1, n)
    cases.append(case('연속 n=%d (정확)' % n, z, 0.5 * z + rng.normal(0, 1, n)))
# 동순위 없음, 34 이상 (점근)
for n in (34, 60, 300):
    z = rng.normal(0, 1, n)
    cases.append(case('연속 n=%d (점근)' % n, z, -0.3 * z + rng.normal(0, 1, n)))
# 동순위가 많은 정수 점수
for n, lv in ((15, 4), (40, 5), (120, 5), (500, 7), (2000, 5)):
    z = rng.normal(0, 1, n)
    x = np.clip(np.round(2.5 + z + rng.normal(0, 0.8, n)), 1, lv)
    y = np.clip(np.round(2.5 + 0.7 * z + rng.normal(0, 1.0, n)), 1, lv)
    cases.append(case('정수 1~%d 점수 n=%d (동순위 많음)' % (lv, n), x, y))
# 완전 일치, 완전 반대, 한쪽만 동순위
x = np.arange(1, 11, dtype=float)
cases.append(case('완전 일치 n=10', x, x * 2 + 1))
cases.append(case('완전 반대 n=10', x, -x))
cases.append(case('한쪽만 동순위 n=25', np.round(rng.normal(0, 1, 25), 2), rng.integers(1, 4, 25)))
cases.append(case('거의 독립 n=3000', rng.normal(0, 1, 3000), rng.normal(0, 1, 3000)))

path = os.path.join(os.path.dirname(__file__), 'reference', 'kendall.json')
with open(path, 'w', encoding='utf-8') as fh:
    json.dump({'cases': cases}, fh, ensure_ascii=False)
print('wrote', len(cases), 'cases')
for c in cases:
    print(c['name'], 'tau', round(c['tau'], 4), 'p_auto', c['p_auto'], 'p_asym', c['p_asym'], c.get('p_exact'))
