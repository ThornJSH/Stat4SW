"""Kaplan-Meier + log-rank 기준값 (개발용; py -3.12 에 statsmodels, lifelines, scipy, numpy 필요).
실행: py -3.12 tests/gen_reference_km.py  → tests/reference/km.json

- 생존함수·Greenwood 표준오차·위험집단 수: statsmodels SurvfuncRight (surv_times 는 사건 시점만)
- log-log 신뢰구간, 백분위 생존시간(중앙값 등)과 그 신뢰구간, 제한 평균: lifelines KaplanMeierFitter
- 집단 비교: statsmodels survdiff (log-rank, gb = Breslow, tw = Tarone-Ware), log-rank 는 lifelines 와도 대조
"""
import json
import math
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.duration.survfunc import SurvfuncRight, survdiff
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from lifelines.utils import qth_survival_times, restricted_mean_survival_time

rng = np.random.default_rng(20261026)


def scrub(x):
    if isinstance(x, np.ndarray):
        x = x.tolist()
    if isinstance(x, dict):
        return {k: scrub(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [scrub(v) for v in x]
    if isinstance(x, (np.floating, float)):
        x = float(x)
        return None if (math.isnan(x) or math.isinf(x)) else x
    if isinstance(x, (np.integer,)):
        return int(x)
    return x


def one_group(time, event, conf):
    time = np.asarray(time, float)
    event = np.asarray(event, int)
    sf = SurvfuncRight(time, event)
    z = norm.ppf(1 - (1 - conf / 100) / 2)
    tt = np.asarray(sf.surv_times, float)
    S = np.asarray(sf.surv_prob, float)
    se = np.asarray(sf.surv_prob_se, float)
    nr = np.asarray(sf.n_risk, float)
    dd = np.asarray(sf.n_events, float)
    kmf = KaplanMeierFitter(alpha=1 - conf / 100)
    kmf.fit(time, event)
    ci = kmf.confidence_interval_
    # lifelines 의 구간을 사건 시점에 맞춘다
    lo = [float(ci.iloc[:, 0].loc[t]) for t in tt]
    hi = [float(ci.iloc[:, 1].loc[t]) for t in tt]
    tau = float(time.max())
    out = {
        'n': int(len(time)), 'events': int(event.sum()), 'tau': tau,
        'tt': tt, 'S': S, 'se': se, 'nr': nr, 'dd': dd, 'lo': lo, 'hi': hi,
        'rmst': float(restricted_mean_survival_time(kmf, t=tau)),
        'quantiles': {}
    }
    for p in (0.25, 0.5, 0.75):
        q = float(qth_survival_times(1 - p, kmf.survival_function_))  # lifelines 의 q 는 생존율 수준(S ≤ q)
        qci = qth_survival_times(1 - p, ci)
        sm_q = sf.quantile(p)
        tie = bool(np.any(np.isclose(S, 1 - p, rtol=0, atol=1e-12)))
        try:
            sm_ci = sf.quantile_ci(p, method='cloglog', alpha=1 - conf / 100)
        except Exception:
            sm_ci = (float('nan'), float('nan'))
        out['quantiles'][str(p)] = {
            't': q, 'lo': float(qci.iloc[0, 0]) if hasattr(qci, 'iloc') else float(qci[0]), 'hi': float(qci.iloc[0, 1]) if hasattr(qci, 'iloc') else float(qci[1]),
            'sm_t': float(sm_q), 'tie': tie, 'sm_lo': float(sm_ci[0]), 'sm_hi': float(sm_ci[1])
        }
    # 임의 시점의 생존율 (statsmodels 단계함수)
    qs = sorted(set([float(np.round(t, 6)) for t in np.linspace(0, tau * 1.1, 9)] + [float(tt[0]) if len(tt) else 0.0]))
    pts = []
    for t in qs:
        k = int(np.searchsorted(tt, t, side='right')) - 1
        s = 1.0 if k < 0 else float(S[k])
        pts.append({'t': t, 's': s, 'nrisk': int((time >= t).sum()), 'k': k})
    out['at'] = pts
    return out


def compare(time, event, group, labels):
    time = np.asarray(time, float)
    event = np.asarray(event, int)
    group = np.asarray(group)
    res = {}
    for key, w in (('logrank', None), ('breslow', 'gb'), ('tarone', 'tw')):
        chi, p = survdiff(time, event, group, weight_type=w)
        res[key] = {'chi2': float(chi), 'p': float(p), 'df': int(len(np.unique(group)) - 1)}
    df = pd.DataFrame({'t': time, 'e': event, 'g': group})
    ll = multivariate_logrank_test(df['t'], df['g'], df['e'])
    res['ll_logrank_chi2'] = float(ll.test_statistic)
    res['ll_logrank_p'] = float(ll.p_value)
    pairs = []
    ug = list(np.unique(group))
    for i in range(len(ug)):
        for j in range(i + 1, len(ug)):
            m = (group == ug[i]) | (group == ug[j])
            chi, p = survdiff(time[m], event[m], group[m])
            pairs.append({'i': int(i), 'j': int(j), 'chi2': float(chi), 'p': float(p)})
    res['pairs'] = pairs
    return res


cases = []


def add(name, time, event, group=None, conf=95, note=''):
    time = np.asarray(time, float)
    event = np.asarray(event, int)
    c = {'name': name, 'note': note, 'conf': conf, 'time': time, 'event': event}
    if group is None:
        c['group'] = None
        c['total'] = one_group(time, event, conf)
    else:
        group = np.asarray(group, int)
        c['group'] = group
        ug = sorted(set(group.tolist()))
        c['groups'] = [one_group(time[group == g], event[group == g], conf) for g in ug]
        c['total'] = one_group(time, event, conf)
        c['compare'] = compare(time, event, group, ug)
    cases.append(c)


# 1. 중도절단 없음, 집단 없음
t = np.round(rng.exponential(20, 40), 2)
add('km-nocens', t, np.ones(40, int), None, 95, '중도절단 없음, 마지막 관측이 사건이라 생존율이 0 이 된다')

# 2. 60% 중도절단, 집단 없음, 99% 신뢰구간
t = np.round(rng.exponential(30, 80), 2)
c = np.round(rng.uniform(5, 40, 80), 2)
obs = np.minimum(t, c)
ev = (t <= c).astype(int)
add('km-heavy-cens', obs, ev, None, 99, '중도절단이 많은 자료')

# 3. 정수 시간(동점 많음), 두 집단 불균형
n1, n2 = 60, 40
t1 = np.ceil(rng.exponential(12, n1)); t2 = np.ceil(rng.exponential(7, n2))
c1 = np.ceil(rng.uniform(1, 25, n1)); c2 = np.ceil(rng.uniform(1, 25, n2))
time = np.concatenate([np.minimum(t1, c1), np.minimum(t2, c2)])
event = np.concatenate([(t1 <= c1), (t2 <= c2)]).astype(int)
group = np.array([0] * n1 + [1] * n2)
add('km-ties-2g', time, event, group, 95, '정수 시간으로 동점이 많고 두 집단 크기가 다르다')

# 4. 세 집단 불균형
sizes = [70, 50, 30]
scales = [20, 12, 8]
tl, el, gl = [], [], []
for g, (nn, sc) in enumerate(zip(sizes, scales)):
    t = np.round(rng.exponential(sc, nn), 1)
    c = np.round(rng.uniform(2, 40, nn), 1)
    tl.append(np.minimum(t, c)); el.append((t <= c).astype(int)); gl.append(np.full(nn, g))
add('km-3g', np.concatenate(tl), np.concatenate(el), np.concatenate(gl), 90, '세 집단, 90% 신뢰구간')

# 5. 사건이 없는 집단이 있는 경우 (집단 2 는 전부 중도절단, 늦은 시점까지 위험집단에 있음)
t = np.round(rng.exponential(10, 60), 1)
time = np.concatenate([t[:30], t[30:60], np.round(rng.uniform(30, 50, 20), 1)])
event = np.concatenate([np.ones(30, int), np.ones(30, int), np.zeros(20, int)])
group = np.array([0] * 30 + [1] * 30 + [2] * 20)
add('km-noevent-group', time, event, group, 95, '한 집단은 사건이 없고 모두 중도절단')

# 6. 마지막 관측이 중도절단
t = np.round(rng.exponential(15, 50), 2)
ev = np.ones(50, int)
i = int(np.argmax(t))
ev[i] = 0
add('km-last-cens', t, ev, None, 95, '가장 큰 시간이 중도절단이라 생존율이 0 으로 내려가지 않는다')

# 7. 시간 0 이 있는 자료와 두 집단
t = np.round(rng.exponential(10, 50), 1)
t[:3] = 0.0
e = (rng.uniform(size=50) < 0.7).astype(int)
e[:3] = [1, 0, 1]
add('km-zero-time', t, e, np.array([0, 1] * 25), 95, '시간 0 이 있는 자료, 두 집단')

# 8. 세 번째 집단이 첫 사건 전에 모두 중도절단 (위험집단에 한 번도 들지 않음) → 비교는 두 집단과 같다
t01 = np.round(rng.exponential(10, 40) + 5, 1)
e01 = (rng.uniform(size=40) < 0.8).astype(int)
t2c = np.round(rng.uniform(0.5, 1.5, 8), 1)
time = np.concatenate([t01, t2c])
event = np.concatenate([e01, np.zeros(8, int)])
group = np.array([0] * 20 + [1] * 20 + [2] * 8)
# survdiff 는 위험집단에 들지 않는 집단 때문에 특이행렬이 되므로, 기준은 그 집단을 뺀 두 집단 비교로 만든다
m = group < 2
case = {'name': 'km-early-cens-group', 'note': '세 번째 집단이 첫 사건 전에 모두 중도절단되어 비교에서 빠진다 (기준은 앞의 두 집단)', 'conf': 95,
        'time': time, 'event': event, 'group': group}
case['groups'] = [one_group(time[group == g], event[group == g], 95) for g in (0, 1, 2)]
case['total'] = one_group(time, event, 95)
case['compare'] = compare(time[m], event[m], group[m], [0, 1])
case['compare_df'] = 1
cases.append(case)

with open('tests/reference/km.json', 'w', encoding='utf-8') as f:
    json.dump(scrub({'cases': cases}), f, ensure_ascii=False)
print('cases', len(cases))
for c in cases:
    tot = c['total']
    q = tot['quantiles']['0.5']
    print(c['name'], 'n', tot['n'], 'events', tot['events'], 'median', q['t'], 'lifelines CI', q['lo'], q['hi'], 'sm', q['sm_t'], q['sm_lo'], q['sm_hi'], 'tie', q['tie'])
    if c.get('compare'):
        print('   ', {k: (round(v['chi2'], 4), round(v['p'], 5)) for k, v in c['compare'].items() if k in ('logrank', 'breslow', 'tarone')}, 'll', round(c['compare']['ll_logrank_chi2'], 4))
