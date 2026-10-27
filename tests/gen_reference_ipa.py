"""IPA(중요도–수행도 분석) 기준값을 tests/reference/ipa.json 에 저장한다 (개발용; py -3.12, numpy·scipy 필요).

실행: py -3.12 tests/gen_reference_ipa.py

- 문항별 평균·표준편차: numpy (ddof=1), 결측은 문항마다 유효한 값만 쓴다.
- 기준선: 문항 평균들의 평균 / 중앙값 + 이동량 / 직접 지정 값.
- 사분면: 중요도 ≥ 기준선이면 "높음", 수행도 ≥ 기준선이면 "높음" (같으면 높음). 집중 개선 = 중요도 높음·수행도 낮음.
- 격차 검정: scipy.stats.ttest_rel (짝이 모두 있는 사례), 효과크기 d_z = 차이의 평균 / 차이의 표준편차.
"""
import json
import os

import numpy as np
from scipy import stats

rng = np.random.default_rng(20261027)


def likert(mu, n, sd=1.0, levels=5):
    return np.clip(np.round(rng.normal(mu, sd, n)), 1, levels)


def make(n, imp_mu, perf_mu, miss=0.0, levels=5):
    k = len(imp_mu)
    z = rng.normal(0, 0.6, n)
    imp = np.column_stack([np.clip(np.round(imp_mu[j] + z + rng.normal(0, 0.8, n)), 1, levels) for j in range(k)])
    perf = np.column_stack([np.clip(np.round(perf_mu[j] + 0.5 * z + rng.normal(0, 1.0, n)), 1, levels) for j in range(k)])
    imp = imp.astype(object)
    perf = perf.astype(object)
    if miss > 0:
        for arr in (imp, perf):
            m = rng.uniform(size=arr.shape) < miss
            arr[m] = None
    return imp, perf


def col_stats(col):
    v = np.array([x for x in col if x is not None], float)
    return {'n': int(len(v)), 'mean': float(v.mean()), 'sd': float(v.std(ddof=1))}


def analyse(imp, perf, scenario, listwise=False):
    n, k = imp.shape
    use = np.ones(n, bool)
    if listwise:
        for arr in (imp, perf):
            use &= np.array([all(x is not None for x in row) for row in arr])
    imp_u, perf_u = imp[use], perf[use]
    items = []
    for j in range(k):
        si, sp = col_stats(imp_u[:, j]), col_stats(perf_u[:, j])
        both = [(a, b) for a, b in zip(imp_u[:, j], perf_u[:, j]) if a is not None and b is not None]
        a = np.array([p[0] for p in both], float)
        b = np.array([p[1] for p in both], float)
        d = b - a
        item = {'imp': si, 'perf': sp, 'pairs': len(both)}
        if len(both) >= 3 and d.std(ddof=1) > 0:
            tt = stats.ttest_rel(b, a)
            item.update({'mdiff': float(d.mean()), 'sddiff': float(d.std(ddof=1)), 't': float(tt.statistic),
                         'df': len(both) - 1, 'p': float(tt.pvalue), 'dz': float(d.mean() / d.std(ddof=1))})
        else:
            item.update({'mdiff': float(d.mean()) if len(both) else None, 'sddiff': None, 't': None, 'df': None, 'p': None, 'dz': None})
        items.append(item)
    im = np.array([i['imp']['mean'] for i in items])
    pm = np.array([i['perf']['mean'] for i in items])
    base = scenario['baseline']
    if base == 'mean':
        bi, bp = im.mean(), pm.mean()
    elif base == 'median':
        bi, bp = float(np.median(im)), float(np.median(pm))
    else:
        bi, bp = scenario['impLine'], scenario['perfLine']
    si = scenario.get('impShift', 0.0) if base != 'custom' else 0.0
    sp = scenario.get('perfShift', 0.0) if base != 'custom' else 0.0
    li, lp = float(bi) + si, float(bp) + sp
    names = ['집중 개선', '현상 유지', '낮은 우선순위', '과잉 투자']
    for it, m_i, m_p in zip(items, im, pm):
        hi, hp = m_i >= li - 1e-12, m_p >= lp - 1e-12
        it['quadrant'] = names[0] if (hi and not hp) else names[1] if (hi and hp) else names[3] if (not hi and hp) else names[2]
        it['edge'] = bool(abs(m_i - li) <= 0.1 + 1e-12 or abs(m_p - lp) <= 0.1 + 1e-12)
        it['gap'] = float(m_p - m_i)
    return {'scenario': scenario, 'listwise': listwise, 'n': int(use.sum()), 'base_imp': float(bi), 'base_perf': float(bp),
            'line_imp': li, 'line_perf': lp, 'items': items}


cases = []
for name, n, imp_mu, perf_mu, miss in [
    ('5문항, 결측 없음 (N=100)', 100, [4.2, 3.4, 3.9, 2.8, 3.6], [3.1, 3.8, 4.1, 3.0, 2.7], 0.0),
    ('8문항, 결측 있음 (N=150, 약 6%)', 150, [4.4, 3.2, 3.8, 4.1, 2.6, 3.5, 4.0, 3.0], [3.0, 3.9, 4.2, 3.5, 3.1, 2.8, 4.0, 3.6], 0.06),
    ('4문항, 작은 표본 (N=25)', 25, [4.0, 3.0, 3.6, 4.4], [3.4, 3.6, 2.8, 3.9], 0.0),
]:
    imp, perf = make(n, imp_mu, perf_mu, miss)
    scenarios = [
        {'baseline': 'mean', 'impShift': 0.0, 'perfShift': 0.0},
        {'baseline': 'median', 'impShift': -0.1, 'perfShift': 0.2},
        {'baseline': 'custom', 'impLine': 3.0, 'perfLine': 3.5},
        {'baseline': 'mean', 'impShift': 0.15, 'perfShift': -0.05},
    ]
    res = [analyse(imp, perf, s) for s in scenarios]
    res.append(analyse(imp, perf, scenarios[0], listwise=True))
    cases.append({'name': name, 'n': n, 'k': imp.shape[1], 'imp': imp.tolist(), 'perf': perf.tolist(), 'results': res})

path = os.path.join(os.path.dirname(__file__), 'reference', 'ipa.json')
with open(path, 'w', encoding='utf-8') as fh:
    json.dump({'cases': cases}, fh, ensure_ascii=False)
print('wrote', len(cases), 'cases')
for c in cases:
    r = c['results'][0]
    print(c['name'], 'lines', round(r['line_imp'], 3), round(r['line_perf'], 3), [i['quadrant'][:2] + ('*' if i['edge'] else '') for i in r['items']])
