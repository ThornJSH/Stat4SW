"""Cochran Q 검정 기준값을 tests/reference/cochran.json 에 저장한다 (개발용; py -3.12, scipy·statsmodels 필요).

실행: py -3.12 tests/gen_reference_cochran.py

- Q, 유의확률: statsmodels.stats.contingency_tables.cochrans_q
- 같은 자료에 scipy.stats.friedmanchisquare 를 적용하면 (0/1 자료에서) Q 와 같다는 항등식으로 교차 확인한다.
- 쌍별 비교: statsmodels mcnemar(exact=True) 의 정확 유의확률
- 효과크기 Kendall's W = Q / (N (k - 1)) 는 식으로 계산한다.
"""
import itertools
import json
import os

import numpy as np
from scipy import stats
from statsmodels.stats.contingency_tables import cochrans_q, mcnemar

rng = np.random.default_rng(20261026)


def make(n, k, probs, rho=0.0):
    """사람마다 성향(잠재 효과)을 두어 반응이 상관되도록 만든 n × k 이진 자료"""
    z = rng.normal(0, 1, (n, 1))
    lat = np.array(probs)[None, :] + rho * z + rng.normal(0, 1, (n, k))
    return (lat > 0).astype(int)


sets = {
    '3회 측정 (N=30, 성공 비율이 점점 오른다)': make(30, 3, [-0.4, 0.0, 0.5], 0.8),
    '4회 측정 (N=50, 한 시점만 다르다)': make(50, 4, [0.0, 0.0, 0.0, 1.0], 1.0),
    '5회 측정 (N=80, 차이가 없다)': make(80, 5, [0.1, 0.1, 0.1, 0.1, 0.1], 1.2),
    '3회 측정 (N=12, 모두 같은 응답인 사례가 많다)': make(12, 3, [0.0, 0.2, 1.0], 2.5),
    '6회 측정 (N=200)': make(200, 6, [-0.3, -0.1, 0.0, 0.2, 0.3, 0.6], 0.7),
}

cases = []
for name, x in sets.items():
    n, k = x.shape
    q = cochrans_q(x)
    fr = stats.friedmanchisquare(*[x[:, j].astype(float) for j in range(k)])
    assert abs(q.statistic - fr.statistic) < 1e-9 * max(1, q.statistic), (name, q.statistic, fr.statistic)
    assert abs(q.pvalue - fr.pvalue) < 1e-9, (name, q.pvalue, fr.pvalue)
    pairs = []
    for i, j in itertools.combinations(range(k), 2):
        a = int(((x[:, i] == 0) & (x[:, j] == 0)).sum())
        b = int(((x[:, i] == 0) & (x[:, j] == 1)).sum())
        c = int(((x[:, i] == 1) & (x[:, j] == 0)).sum())
        d = int(((x[:, i] == 1) & (x[:, j] == 1)).sum())
        m = mcnemar([[a, b], [c, d]], exact=True)
        pairs.append({'i': i, 'j': j, 'ps_fail': c, 'pf_succ': b, 'p': float(m.pvalue)})
    informative = int(((x.sum(axis=1) > 0) & (x.sum(axis=1) < k)).sum())
    cases.append({
        'name': name, 'n': n, 'k': k, 'x': x.tolist(),
        'q': float(q.statistic), 'df': int(q.df), 'p': float(q.pvalue),
        'w': float(q.statistic / (n * (k - 1))), 'informative': informative,
        'colsum': x.sum(axis=0).tolist(), 'pairs': pairs,
    })

path = os.path.join(os.path.dirname(__file__), 'reference', 'cochran.json')
with open(path, 'w', encoding='utf-8') as fh:
    json.dump({'cases': cases}, fh, ensure_ascii=False)
print('wrote', len(cases), 'cases')
for c in cases:
    print(c['name'], 'Q', round(c['q'], 4), 'p', round(c['p'], 5), 'informative', c['informative'], c['colsum'])
