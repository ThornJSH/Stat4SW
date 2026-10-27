const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const S = require('../src/Stats_Survival.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');
const { qnorm, pchisq } = require('../src/Distributions.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'km.json'), 'utf8'));
const near = (a, e, tol, abs = 1e-12) => (a === null || e === null ? a === e : Math.abs(a - e) <= Math.max(abs, tol * Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const def = (name, scale, vl) => ({ name, label: '', scale: scale || 'scale', valueLabels: vl || {}, missing: [] });
const zOf = conf => qnorm(1 - (1 - conf / 100) / 2);

function dsOfCase(c, groupLabels) {
  const rows = [['t', 'e', 'g']];
  c.time.forEach((t, i) => rows.push([t, c.event[i], c.group ? c.group[i] : 1]));
  return Data.buildDataset(rows, [def('t'), def('e', 'nominal'), def('g', 'nominal', groupLabels || {})]);
}

function checkGroup(km, e, m) {
  assert.equal(km.n, e.n, m + ' n');
  assert.equal(km.events, e.events, m + ' events');
  assert.equal(km.tt.length, e.tt.length, m + ' 사건 시점 수');
  e.tt.forEach((t, k) => {
    assert.equal(km.tt[k], t, `${m} tt[${k}]`);
    assert.equal(km.nr[k], e.nr[k], `${m} nr[${k}]`);
    assert.equal(km.dd[k], e.dd[k], `${m} dd[${k}]`);
    assert.ok(near(km.S[k], e.S[k], 1e-12, 1e-14), `${m} S[${k}]`);
    if (e.S[k] > 0) {
      assert.ok(near(km.se[k], e.se[k], 1e-9, 1e-12), `${m} se[${k}]: ${km.se[k]} vs ${e.se[k]}`);
      assert.ok(near(km.lo[k], e.lo[k], 1e-8, 1e-10), `${m} lo[${k}]: ${km.lo[k]} vs ${e.lo[k]}`);
      assert.ok(near(km.hi[k], e.hi[k], 1e-8, 1e-10), `${m} hi[${k}]: ${km.hi[k]} vs ${e.hi[k]}`);
    }
  });
  assert.ok(near(S.svRestrictedMean(km, e.tau), e.rmst, 1e-10, 1e-12), `${m} 제한 평균 ${S.svRestrictedMean(km, e.tau)} vs ${e.rmst}`);
  Object.keys(e.quantiles).forEach(p => {
    const q = S.svQuantile(km, 1 - Number(p)), eq = e.quantiles[p];
    if (eq.tie) {   // 생존율이 정확히 그 수준이 되는 경우: 부동소수점 잡음 때문에 lifelines 는 다음 시점을 고를 수 있어, 수학적으로 맞는 이른 시점인지 본다
      assert.ok(q.t <= eq.t && Math.abs(km.S[km.tt.indexOf(q.t)] - (1 - Number(p))) < 1e-12, `${m} ${p} 동률 시점`);
    } else assert.ok(near(q.t, eq.t, 1e-12), `${m} ${p} 시점 ${q.t} vs ${eq.t}`);
    assert.ok(near(q.lo, eq.lo, 1e-12), `${m} ${p} 하한 ${q.lo} vs ${eq.lo}`);
    assert.ok(near(q.hi, eq.hi, 1e-12), `${m} ${p} 상한 ${q.hi} vs ${eq.hi}`);
    if (!eq.tie) assert.ok(near(q.t, eq.sm_t, 1e-12), `${m} ${p} statsmodels 시점`);
    if (eq.sm_lo !== null && eq.sm_hi !== null && q.t !== null) {
      assert.ok(near(q.lo, eq.sm_lo, 1e-12) && (q.hi === eq.sm_hi || (q.hi === null && eq.sm_hi === null)), `${m} ${p} statsmodels 신뢰구간 ${q.lo},${q.hi} vs ${eq.sm_lo},${eq.sm_hi}`);
    }
  });
  e.at.forEach(a => {
    const r = S.svSurvivalAt(km, a.t);
    if (r.s !== null) assert.ok(near(r.s, a.s, 1e-12, 1e-14), `${m} S(${a.t})`);
  });
}

test('Kaplan–Meier: statsmodels(생존율·Greenwood·위험집단)·lifelines(log-log 신뢰구간·백분위·제한 평균)와 같다 — 8가지 자료의 모든 집단', () => {
  assert.equal(ref.cases.length, 8);
  ref.cases.forEach(c => {
    const z = zOf(c.conf);
    const tot = S.svKm(c.time, c.event, z);
    checkGroup(tot, c.total, c.name + ' 전체');
    if (c.group) {
      const ug = Array.from(new Set(c.group)).sort((a, b) => a - b);
      ug.forEach((g, q) => {
        const idx = c.group.map((x, i) => (x === g ? i : -1)).filter(i => i >= 0);
        checkGroup(S.svKm(idx.map(i => c.time[i]), idx.map(i => c.event[i]), z), c.groups[q], `${c.name} 집단${g}`);
      });
    }
  });
});

test('Kaplan–Meier: 흔히 쓰는 교과서 예(6-MP 백혈병 21명 중 치료군)의 값', () => {
  const time = [6, 6, 6, 7, 10, 13, 16, 22, 23, 6, 9, 10, 11, 17, 19, 20, 25, 32, 32, 34, 35];
  const ev = [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
  const km = S.svKm(time, ev, 1.96);
  [0.8571, 0.8067, 0.7529, 0.6902, 0.6275, 0.5378, 0.4482].forEach((v, k) => assert.ok(Math.abs(km.S[k] - v) < 5e-5));
  assert.equal(S.svQuantile(km, 0.5).t, 23);
});

test('집단 비교: statsmodels survdiff(log-rank, Breslow=gb, Tarone–Ware=tw)와 lifelines log-rank, 쌍별 검정이 같다', () => {
  ref.cases.filter(c => c.compare).forEach(c => {
    const ug = Array.from(new Set(c.group)).sort((a, b) => a - b);
    const grp = c.group, all = ug.map((g, q) => q);
    const r = S.svCompare(c.time, c.event, grp.map(g => ug.indexOf(g)), all, pchisq);
    const e = c.compare;
    ['logrank', 'breslow', 'tarone'].forEach((key, k) => {
      assert.ok(near(r[k].chi2, e[key].chi2, 1e-9, 1e-10), `${c.name} ${key} χ² ${r[k].chi2} vs ${e[key].chi2}`);
      assert.ok(near(r[k].p, e[key].p, 1e-7, 1e-14), `${c.name} ${key} p`);
      assert.equal(r[k].df, c.compare_df || e[key].df, `${c.name} ${key} df`);
    });
    if (!c.compare_df) assert.ok(near(r[0].chi2, e.ll_logrank_chi2, 1e-9, 1e-10));
    if (!c.compare_df) e.pairs.forEach(pr => {
      const rr = S.svCompare(c.time, c.event, grp.map(g => ug.indexOf(g)), [pr.i, pr.j], pchisq)[0];
      assert.ok(near(rr.chi2, pr.chi2, 1e-9, 1e-10) && near(rr.p, pr.p, 1e-7, 1e-14), `${c.name} 쌍 ${pr.i}-${pr.j}`);
    });
  });
});

test('runKm: 집단별 처리 요약·시점별 생존율·백분위·검정·쌍별(Bonferroni)·그림이 기준값과 같다 (세 집단, 90% 신뢰구간)', () => {
  const c = ref.cases.find(x => x.name === 'km-3g');
  const ds = dsOfCase(c, { 0: '대조', 1: '약 A', 2: '약 B' });
  const res = S.runKm(ds, { analysis: 'km', time: ['t'], event: ['e'], group: ['g'], options: { eventValue: '1', conf: '90', times: '5, 10, 20', table: false, tests: true, pairs: true, plot: true } });
  const sum = table(res, '처리 요약');
  assert.deepEqual(sum.rows.map(r => r.label), ['대조', '약 A', '약 B', '전체']);
  c.groups.forEach((g, q) => {
    assert.equal(sum.rows[q].n, g.n); assert.equal(sum.rows[q].events, g.events);
    assert.ok(near(sum.rows[q].mean, g.rmst, 1e-10));
  });
  const at = table(res, '시점별 생존율');
  assert.equal(at.rows.length, 9);
  const z = zOf(90);
  c.groups.forEach((g, q) => [5, 10, 20].forEach((t, k) => {
    const row = at.rows[q * 3 + k], km = S.svKm(c.time.filter((_, i) => c.group[i] === q), c.event.filter((_, i) => c.group[i] === q), z);
    const e = g.at.find(a => a.t === t);
    if (e) assert.ok(near(row.s, e.s, 1e-10));
    assert.equal(row.nr, km.sorted.filter(x => x >= t).length);
  }));
  const qt = table(res, '생존시간 백분위수');
  assert.equal(qt.rows.length, 9);
  c.groups.forEach((g, q) => {
    assert.ok(near(qt.rows[q * 3 + 1].t, g.quantiles['0.5'].t, 1e-12));
    assert.ok(near(qt.rows[q * 3].t, g.quantiles['0.25'].t, 1e-12));
  });
  assert.equal(qt.rows[1].pct.startsWith('중앙값'), true);
  const cmp = table(res, '집단 비교 검정');
  assert.ok(near(cmp.rows[0].chi2, c.compare.logrank.chi2, 1e-9) && cmp.rows[0].df === 2);
  assert.ok(near(cmp.rows[1].chi2, c.compare.breslow.chi2, 1e-9) && near(cmp.rows[2].chi2, c.compare.tarone.chi2, 1e-9));
  const pw = table(res, '집단 쌍별');
  assert.equal(pw.rows.length, 3);
  pw.rows.forEach((r, k) => assert.ok(near(r.pAdj, Math.min(1, c.compare.pairs[k].p * 3), 1e-7, 1e-14)));
  const chart = res.blocks.find(b => b.type === 'chart');
  assert.equal(chart.kind, 'line');
  assert.equal(chart.series.length, 3);
  chart.series.forEach((s, q) => {
    assert.deepEqual(s.points[0], [0, 1]);
    for (let k = 1; k < s.points.length; k++) assert.ok(s.points[k][1] <= s.points[k - 1][1] + 1e-15 && s.points[k][0] >= s.points[k - 1][0]);
    assert.ok(near(s.points[s.points.length - 1][1], c.groups[q].S[c.groups[q].S.length - 1], 1e-12));
  });
  assert.ok(R.renderResultHtml(res).includes('Kaplan–Meier') && Ex.resultToGrid(res).rows.length > 20);
});

test('runKm: 집단 없이 하나의 곡선, 사건 시점 표 옵션, 시점 자동 선택, 사건값 지정(문자)', () => {
  const c = ref.cases.find(x => x.name === 'km-heavy-cens');
  const ds = dsOfCase(c);
  const res = S.runKm(ds, { analysis: 'km', time: ['t'], event: ['e'], group: [], options: { eventValue: '1', conf: '99', table: true } });
  assert.ok(!table(res, '집단 비교 검정'));
  assert.equal(table(res, '처리 요약').rows.length, 1);
  const full = table(res, '사건 시점별 생존표');
  assert.equal(full.rows.length, c.total.tt.length);
  assert.ok(near(full.rows[3].se, c.total.se[3], 1e-9) && near(full.rows[3].lo, c.total.lo[3], 1e-8));
  const at = table(res, '시점별 생존율');
  assert.ok(at.rows.length >= 3 && at.rows.length <= 8);
  assert.ok(at.notes.some(n => /자동으로/.test(n)));
  // 문자 사건값
  const rows = [['t', 'ev']].concat(c.time.map((t, i) => [t, c.event[i] ? '사망' : '생존']));
  const ds2 = Data.buildDataset(rows, [def('t'), def('ev', 'nominal')]);
  const r2 = S.runKm(ds2, { analysis: 'km', time: ['t'], event: ['ev'], group: [], options: { eventValue: '사망', conf: '99' } });
  assert.equal(table(r2, '처리 요약').rows[0].events, c.total.events);
});

test('runKm: 오류와 경고 — 음수 시간, 사건 없음, 집단 1개, 집단 11개, 시점 입력 오류, 결측 제외, 사건 변수 값 3가지', () => {
  const mk = (rows, defs) => Data.buildDataset([['t', 'e', 'g']].concat(rows), defs || [def('t'), def('e', 'nominal'), def('g', 'nominal')]);
  const req = (extra, opts) => Object.assign({ analysis: 'km', time: ['t'], event: ['e'], group: [], options: Object.assign({ eventValue: '1' }, opts || {}) }, extra || {});
  assert.throws(() => S.runKm(mk([[-1, 1, 1], [2, 1, 1], [3, 0, 1]]), req()), /음수가 1건/);
  assert.throws(() => S.runKm(mk([[1, 0, 1], [2, 0, 1], [3, 0, 1]]), req()), /사건이 한 건도 없습니다/);
  assert.throws(() => S.runKm(mk([[1, 1, 1], [2, 1, 1], [3, 0, 1]]), req({ group: ['g'] })), /집단이 1개/);
  const many = []; for (let i = 0; i < 44; i++) many.push([i + 1, i % 2, (i % 11) + 1]);
  assert.throws(() => S.runKm(mk(many), req({ group: ['g'] })), /11개.*최대 10개/);
  assert.throws(() => S.runKm(mk([[1, 1, 1], [2, 1, 1], [3, 0, 1]]), req({}, { times: '5, abc' })), /시점/);
  assert.throws(() => S.runKm(mk([[1, 1, 1], [2, 1, 1], [3, 0, 1]]), req({}, { conf: '40' })), /신뢰수준/);
  assert.throws(() => S.runKm(mk([[1, 1, 1]]), req()), /2건 이상/);
  const miss = S.runKm(mk([[1, 1, 1], [2, 1, 1], [3, 0, 1], [4, 1, null], [null, 1, 1], [5, 1, 2], [6, 0, 2], [7, 1, 2]]), req({ group: ['g'] }));
  assert.ok(miss.warnings.some(w => /2건을 제외/.test(w)));
  const three = S.runKm(mk([[1, 1, 1], [2, 2, 1], [3, 0, 1], [4, 1, 1], [5, 2, 1], [6, 0, 1], [7, 1, 1]]), req());
  assert.ok(three.warnings.some(w => /값이 3가지/.test(w)));
  assert.ok(three.warnings.some(w => /사건이 3건뿐/.test(w)));
});

test('runKm: 사건 없는 집단·마지막 관측이 중도절단인 경우 경고와 "도달하지 않음" 표시', () => {
  const c = ref.cases.find(x => x.name === 'km-noevent-group');
  const res = S.runKm(dsOfCase(c), { analysis: 'km', time: ['t'], event: ['e'], group: ['g'], options: { eventValue: '1' } });
  assert.ok(res.warnings.some(w => /사건이 한 건도 없는 집단/.test(w)));
  const qt = table(res, '생존시간 백분위수');
  assert.equal(qt.rows[7].t, '도달하지 않음');
  const html = R.renderResultHtml(res);
  assert.ok(html.includes('도달하지 않음'));
});

test('runKm: 첫 사건 전에 모두 중도절단된 집단은 비교에서 빠지고 자유도가 줄어든다', () => {
  const c = ref.cases.find(x => x.name === 'km-early-cens-group');
  const res = S.runKm(dsOfCase(c), { analysis: 'km', time: ['t'], event: ['e'], group: ['g'], options: { eventValue: '1' } });
  const cmp = table(res, '집단 비교 검정');
  assert.equal(cmp.rows[0].df, 1);
  assert.ok(near(cmp.rows[0].chi2, c.compare.logrank.chi2, 1e-9));
  assert.ok(cmp.notes.some(n => /자유도를 1로/.test(n)));
});
