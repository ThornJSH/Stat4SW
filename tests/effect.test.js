const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const dist = require('../src/Distributions.js');
const E = require('../src/Stats_Effect.js');
const C = require('../src/Stats_Compare.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'effect.json'), 'utf8'));
const near = (a, e, tol) => (a === null || e === null ? a === e : Math.abs(a - e) <= tol * Math.max(1e-10, Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const stats = v => { const n = v.length, m = v.reduce((s, x) => s + x, 0) / n; return { n, mean: m, variance: v.reduce((s, x) => s + (x - m) * (x - m), 0) / (n - 1), sd: Math.sqrt(v.reduce((s, x) => s + (x - m) * (x - m), 0) / (n - 1)) }; };

test('pnct: scipy.stats.nct.cdf 와 436개 조합에서 일치 (자유도 1~50000, 음수 t·비중심 모수 포함)', () => {
  assert.ok(ref.grid.length >= 400);
  let worst = 0;
  ref.grid.forEach(g => {
    const v = dist.pnct(g.t, g.df, g.ncp);
    const err = Math.abs(v - g.p);
    worst = Math.max(worst, err);
    assert.ok(err <= 1e-9 + 1e-8 * Math.abs(g.p), `pnct(${g.t}, ${g.df}, ${g.ncp}) = ${v} vs ${g.p}`);
  });
  assert.ok(worst < 1e-8);
  // ncp = 0 이면 중심 t 와 같고, 대칭성 P(T ≤ t; δ) = 1 − P(T ≤ −t; −δ)
  assert.ok(Math.abs(dist.pnct(1.7, 12, 0) - dist.pt(1.7, 12)) < 1e-12);
  assert.ok(Math.abs(dist.pnct(-1.7, 12, 0.9) - (1 - dist.pnct(1.7, 12, -0.9))) < 1e-12);
  assert.equal(dist.pnct(Infinity, 5, 1), 1);
  assert.equal(dist.pnct(-Infinity, 5, 1), 0);
  assert.ok(Number.isNaN(dist.pnct(1, 0, 1)));
});

test('효과크기: 9가지 자료의 d, Hedges g, 정확 신뢰구간, d_av 가 scipy 기준값과 같다', () => {
  assert.equal(ref.cases.length, 9);
  ref.cases.forEach(c => {
    let e;
    if (c.kind === 'independent') e = E.effectIndependent(stats(c.a), stats(c.b), c.conf);
    else if (c.kind === 'paired') {
      const diff = stats(c.x.map((v, i) => v - c.y[i]));
      e = E.effectPaired(diff, stats(c.x).sd, stats(c.y).sd, c.conf);
      assert.ok(near(e.dav, c.d_av, 1e-10), c.name + ' d_av');
      assert.ok(near(e.std, c.sd_diff, 1e-10));
    } else e = E.effectOneSample(stats(c.values), c.mu, c.conf);
    const w = `${c.kind} ${c.name}`;
    assert.ok(near(e.d, c.d, 1e-10) && near(e.J, c.J, 1e-12) && near(e.g, c.g, 1e-10), w + ' 점추정');
    assert.ok(near(e.dLo, c.d_lo, 1e-7) && near(e.dHi, c.d_hi, 1e-7), `${w} d 구간 ${e.dLo}..${e.dHi} vs ${c.d_lo}..${c.d_hi}`);
    assert.ok(near(e.gLo, c.g_lo, 1e-7) && near(e.gHi, c.g_hi, 1e-7), w + ' g 구간');
    assert.ok(near(e.t, c.t, 1e-10) && e.df === c.df, w + ' t, df');
    assert.ok(e.dLo < e.d && e.d < e.dHi && Math.abs(e.g) < Math.abs(e.d), w + ' 순서');
  });
  assert.equal(E.effectSizeWord(0.1), '거의 없음');
  assert.equal(E.effectSizeWord(-0.2), '작음');
  assert.equal(E.effectSizeWord(0.5), '중간');
  assert.equal(E.effectSizeWord(-1.2), '큼');
  assert.equal(E.hedgesJ(10), 1 - 3 / 39);
});

test('효과크기: 계산할 수 없으면 null (사례 2개 미만, 분산 0)', () => {
  assert.equal(E.effectIndependent({ n: 1, mean: 1, variance: 0 }, { n: 5, mean: 2, variance: 1 }, 95), null);
  assert.equal(E.effectIndependent({ n: 5, mean: 1, variance: 0 }, { n: 5, mean: 2, variance: 0 }, 95), null);
  assert.equal(E.effectPaired({ n: 5, mean: 1, sd: 0 }, 1, 1, 95), null);
  assert.equal(E.effectOneSample({ n: 1, mean: 3, sd: NaN }, 0, 95), null);
  // t 가 매우 크면 정규 근사 표시
  assert.equal(E.effectCI(9, 100, 99, 95).approx, true);
  assert.equal(E.effectCI(0.5, 40, 78, 95).approx, false);
});

function dsOf(rows, defs) { return Data.buildDataset(rows, defs || []); }

test('독립표본 t검정: 효과크기 표(결과 아래)에 d, g 와 구간이 기준값과 같고, 옵션을 끄면 표가 없다', () => {
  ref.cases.filter(c => c.kind === 'independent').forEach(c => {
    const rows = [['값', '집단']];
    c.a.forEach(v => rows.push([v, 1]));
    c.b.forEach(v => rows.push([v, 2]));
    const ds = dsOf(rows, [{ name: '집단', label: '', scale: 'nominal', valueLabels: { '1': '참여', '2': '미참여' }, missing: [] }]);
    const res = C.runTTestIndependent(ds, { dependent: ['값'], factor: ['집단'], options: { effect: true, conf: String(c.conf) } });
    const t = table(res, '효과크기');
    assert.ok(t, c.name);
    assert.deepEqual(t.rows.map(r => r.kind), ["Cohen's d", "Hedges' g"]);
    assert.ok(near(t.rows[0].est, c.d, 1e-9) && near(t.rows[0].lo, c.d_lo, 1e-7) && near(t.rows[0].hi, c.d_hi, 1e-7), c.name);
    assert.ok(near(t.rows[1].est, c.g, 1e-9) && near(t.rows[1].lo, c.g_lo, 1e-7), c.name + ' g');
    assert.ok(near(t.rows[0].std, c.sp, 1e-10));
    assert.ok(t.columns[4].label.includes(c.conf + '%'));
    assert.ok(res.blocks.findIndex(b => b === t) === res.blocks.findIndex(b => b.caption.startsWith('독립표본 검정')) + 1);   // 검정 표 바로 아래
    assert.equal(table(C.runTTestIndependent(ds, { dependent: ['값'], factor: ['집단'], options: { conf: String(c.conf) } }), '효과크기'), undefined);
  });
  const rows = [['값', '집단'], [1, 1], [1, 1], [2, 2], [2, 2]];                       // 집단 안 변동이 없음
  const w = C.runTTestIndependent(dsOf(rows), { dependent: ['값'], factor: ['집단'], options: { effect: true } });
  assert.ok(!table(w, '효과크기') && w.warnings.some(x => x.includes('효과크기를 계산할 수 없습니다')));
});

test('대응표본 t검정: d_z, g_z, d_av 세 줄과 구간이 기준값과 같다 (d_av 는 구간 없음)', () => {
  ref.cases.filter(c => c.kind === 'paired').forEach(c => {
    const rows = [['전', '후']];
    c.x.forEach((v, i) => rows.push([v, c.y[i]]));
    const res = C.runTTestPaired(dsOf(rows), { first: ['전'], second: ['후'], options: { effect: true, conf: String(c.conf) } });
    const t = table(res, '효과크기');
    assert.deepEqual(t.rows.map(r => r.kind), ["Cohen's d_z", "Hedges' g_z", "Cohen's d_av"]);
    assert.ok(near(t.rows[0].est, c.d, 1e-9) && near(t.rows[0].lo, c.d_lo, 1e-7) && near(t.rows[0].hi, c.d_hi, 1e-7), c.name);
    assert.ok(near(t.rows[1].est, c.g, 1e-9));
    assert.ok(near(t.rows[2].est, c.d_av, 1e-9) && t.rows[2].lo === null && t.rows[2].hi === null);
    assert.ok(t.notes.some(n => n.includes('d_av')));
  });
});

test('일표본 t검정: 변수마다 d, g 두 줄, 렌더·내보내기에 포함', () => {
  const cs = ref.cases.filter(c => c.kind === 'one');
  const rows = [['가', '나']];
  const n = Math.max(cs[0].values.length, cs[1].values.length);
  for (let i = 0; i < n; i++) rows.push([i < cs[0].values.length ? cs[0].values[i] : '', i < cs[1].values.length ? cs[1].values[i] : '']);
  const res = C.runTTestOne(dsOf(rows), { vars: ['가'], options: { testValue: String(cs[0].mu), effect: true } });
  const t = table(res, '효과크기');
  assert.ok(near(t.rows[0].est, cs[0].d, 1e-9) && near(t.rows[0].lo, cs[0].d_lo, 1e-7) && near(t.rows[0].hi, cs[0].d_hi, 1e-7));
  assert.ok(near(t.rows[1].est, cs[0].g, 1e-9));
  const html = R.renderResultHtml(res);
  assert.ok(html.includes('Hedges') && html.includes('크기(참고)') && html.includes('효과크기'));
  assert.ok(Ex.resultToGrid(res).rows.length > 8);
  // 여러 변수: 변수마다 두 줄
  const two = C.runTTestOne(dsOf(rows), { vars: ['가', '나'], options: { testValue: '0', effect: true } });
  assert.equal(table(two, '효과크기').rows.length, 4);
});
