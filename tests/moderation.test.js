const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const S = require('../src/Sample.js');
const M = require('../src/Stats_Moderation.js');
const Ch = require('../src/Charts.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'moderation.json'), 'utf8'));
const near = (a, e, tol) => (a === null || e === null ? a === e : Math.abs(a - e) <= tol * Math.max(1e-12, Math.abs(e)) + 1e-300);
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));

/* ---------- numpy 독립 계산과 대조 ---------- */

test('moderationFit: numpy 기준값 대조 (계수, SE, t, p, 신뢰구간, R², ΔR², 단순기울기, J-N, 그림 예측값)', () => {
  assert.ok(ref.cases.length >= 7);
  for (const c of ref.cases) {
    const i = c.input, e = c.expected;
    const inp = { y: i.y, x: i.x, xName: 'X', covs: i.covs.map((v, k) => ({ name: 'C' + k, values: v })), center: i.center, conf: i.conf };
    inp.w = i.wtype === 'continuous'
      ? { type: 'continuous', name: 'W', values: i.w }
      : { type: 'categorical', name: 'W', codes: i.w, labels: Array.from({ length: i.ncat }, (_, k) => 'cat' + k), refIdx: i.ref };
    const r = M.moderationFit(inp);
    assert.ok(!r.error, c.name);
    const where = c.name;
    e.coef.forEach((x, k) => {
      const y = r.fit.coef[k];
      for (const key of ['B', 'se', 't', 'p', 'lo', 'hi']) assert.ok(near(y[key], x[key], 1e-9), `${where} coef[${k}].${key} ${y[key]} vs ${x[key]}`);
    });
    for (const [a, b, name] of [[r.fit.r2, e.r2, 'R²'], [r.fit.adjR2, e.adj, '수정 R²'], [r.fit.F, e.F, 'F'], [r.fit.pF, e.pF, 'pF'],
      [r.deltaR2, e.dr2, 'ΔR²'], [r.F, e.Fchg, 'F 변화량'], [r.pF, e.pchg, 'p 변화량']]) {
      assert.ok(near(a, b, 1e-8), `${where} ${name} ${a} vs ${b}`);
    }
    assert.equal(r.df1, e.q, where);
    assert.equal(r.df2, e.dfr, where);
    e.slopes.forEach((x, k) => {
      const y = r.slopes[k];
      for (const key of ['b', 'se', 't', 'p', 'lo', 'hi']) assert.ok(near(y[key], x[key], 1e-8), `${where} slope[${k}].${key} ${y[key]} vs ${x[key]}`);
    });
    if (e.jn) {
      assert.equal(r.jn.roots.length, e.jn.roots.length, where);
      e.jn.roots.forEach((v, k) => assert.ok(near(r.jn.roots[k], v, 1e-8), `${where} J-N 경계 ${r.jn.roots[k]} vs ${v}`));
      assert.deepEqual(r.jn.regions.map(g => g.sig), e.jn.regions.map(g => g.sig), where);
      assert.deepEqual(r.jn.regions.map(g => g.sign), e.jn.regions.map(g => g.sign), where);
    } else {
      assert.equal(r.jn, null, where);
    }
    e.plot.forEach((pl, k) => pl.forEach((p, j) => {
      assert.ok(near(r.plot.series[k].points[j][0], p[0], 1e-12) && near(r.plot.series[k].points[j][1], p[1], 1e-8), `${where} 그림 ${k}/${j}`);
    }));
  }
});

test('johnsonNeyman: 경계는 t = 임계값이 되는 지점이고 구간 판정이 일관된다', () => {
  // b1 = 0.3, b3 = 0.8, 분산은 작게: 효과 0.3 + 0.8w 가 w = -0.375 근처에서 0
  const tcrit = 1.98;
  const jn = M.johnsonNeyman(0.3, 0.8, 0.01, 0.0002, 0.004, tcrit, -2, 2);
  assert.equal(jn.roots.length, 2);
  jn.roots.forEach(w => {
    const t = Math.abs(0.3 + 0.8 * w) / Math.sqrt(0.01 + 2 * w * 0.0002 + w * w * 0.004);
    assert.ok(Math.abs(t - tcrit) < 1e-9, `경계에서 |t| = ${t}`);
  });
  assert.deepEqual(jn.regions.map(g => g.sig), [true, false, true]);
  assert.equal(jn.regions[0].sign, -1);                     // 첫 구간(w가 작음)은 효과가 음수, 마지막은 양수
  assert.equal(jn.regions[2].sign, 1);
  // 이차항 계수가 0인 경우 (A = 0): 경계 하나
  const lin = M.johnsonNeyman(1, 0.5, 0.04, 0, 0.01 * 0 + (0.5 * 0.5) / (1.98 * 1.98), 1.98, -3, 3);
  assert.ok(lin.roots.length <= 1);
  // 범위 밖 경계는 버린다
  const out = M.johnsonNeyman(0.3, 0.8, 0.01, 0.0002, 0.004, tcrit, 1, 2);
  assert.equal(out.roots.length, 0);
  assert.equal(out.regions.length, 1);
  assert.equal(out.regions[0].sig, true);
});

/* ---------- 결과 객체 ---------- */

const defs = [
  { name: 'y', label: '성과', scale: 'scale', valueLabels: {}, missing: [999] },
  { name: 'x', label: '경력', scale: 'scale', valueLabels: {}, missing: [] },
  { name: 'w', label: '자율성', scale: 'scale', valueLabels: {}, missing: [] },
  { name: 'g', label: '집단', scale: 'nominal', valueLabels: { '1': '가', '2': '나', '3': '다' }, missing: [] },
  { name: 'c', label: '공변량', scale: 'scale', valueLabels: {}, missing: [] },
  { name: 't', label: '글자', scale: 'nominal', valueLabels: {}, missing: [] }
];
function makeData(n) {
  let seed = 7;
  const u = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return (seed + 0.5) / 4294967296; };
  const norm = () => Math.sqrt(-2 * Math.log(u())) * Math.cos(2 * Math.PI * u());
  const rows = [['y', 'x', 'w', 'g', 'c', 't']];
  for (let i = 0; i < n; i++) {
    const x = 10 + 3 * norm(), w = 50 + 8 * norm(), g = 1 + (i % 3), c = norm();
    const y = 5 + 0.4 * x + 0.1 * w + 0.06 * (x - 10) * (w - 50) + [0, 0.5, 1][g - 1] * (x - 10) + 0.5 * c + norm();
    rows.push([i === 3 ? 999 : y, x, w, g, c, i % 2 ? 'a' : 'b']);
  }
  return Data.buildDataset(rows, defs);
}
const ds = makeData(90);

test('runModeration: 표 구성(모형 요약, 계수, 상호작용 검정, 단순기울기, J-N)과 그림, 결측 경고', () => {
  const res = M.runModeration(ds, { dependent: ['y'], independent: ['x'], moderator: ['w'], covariates: ['c'], options: {} });
  assert.equal(res.title, '조절효과 분석');
  assert.deepEqual(table(res, '모형 요약').columns.map(c => c.key), ['r', 'r2', 'adj', 's', 'F', 'df1', 'df2', 'p']);
  const coef = table(res, '계수');
  assert.deepEqual(coef.rows.map(r => r.name), ['(상수)', '경력', '자율성', '경력 × 자율성', '공변량']);
  assert.deepEqual(coef.rows.map(r => !!r._bold), [false, false, false, true, false]);              // 상호작용 행 굵게
  assert.ok(table(res, '상호작용 검정').rows[0].p < 0.001);
  const slopes = table(res, '조건부 효과');
  assert.deepEqual(slopes.rows.map(r => r.label), ['평균 − 1SD', '평균', '평균 + 1SD']);
  assert.ok(slopes.rows[0].b < slopes.rows[1].b && slopes.rows[1].b < slopes.rows[2].b);
  assert.ok(table(res, 'Johnson-Neyman'));
  assert.ok(res.blocks.some(b => b.type === 'chart' && b.kind === 'line' && b.series.length === 3));
  assert.ok(res.warnings.some(w => w.includes('목록별 결측 제거로 1건')));
  assert.ok(R.renderResultHtml(res).includes('<svg'));
  assert.ok(Ex.resultToGrid(res).rows.length > 15);
});

test('runModeration: 범주형 조절변수(척도 명목) — 범주별 단순기울기, 기준 범주 옵션, 전체 상호작용 검정', () => {
  const res = M.runModeration(ds, { dependent: ['y'], independent: ['x'], moderator: ['g'], options: {} });
  const coef = table(res, '계수');
  assert.deepEqual(coef.rows.map(r => r.name), ['(상수)', '경력', '집단 = 나', '집단 = 다', '경력 × 집단 = 나', '경력 × 집단 = 다']);
  const slopes = table(res, '조건부 효과');
  assert.deepEqual(slopes.rows.map(r => r.label), ['가 (기준)', '나', '다']);
  assert.ok(slopes.rows[0].b < slopes.rows[1].b && slopes.rows[1].b < slopes.rows[2].b);
  assert.equal(table(res, '상호작용 검정').rows[0].df1, 2);
  assert.ok(!table(res, 'Johnson-Neyman'));
  assert.ok(res.blocks.some(b => b.kind === 'line' && b.series.length === 3 && b.series[0].label === '가'));
  const last = M.runModeration(ds, { dependent: ['y'], independent: ['x'], moderator: ['g'], options: { lastRef: true } });
  assert.equal(table(last, '조건부 효과').rows[2].label, '다 (기준)');
});

test('runModeration: 옵션 — 중심화 끄기, 그림·J-N 끄기, 신뢰수준', () => {
  const on = M.runModeration(ds, { dependent: ['y'], independent: ['x'], moderator: ['w'], options: {} });
  const off = M.runModeration(ds, { dependent: ['y'], independent: ['x'], moderator: ['w'], options: { center: false, plot: false, jn: false, conf: '90' } });
  // 상호작용 계수와 ΔR² 는 중심화와 무관
  const b3 = r => table(r, '계수').rows[3].B;
  assert.ok(Math.abs(b3(on) - b3(off)) < 1e-9);
  assert.ok(Math.abs(table(on, '상호작용 검정').rows[0].dr2 - table(off, '상호작용 검정').rows[0].dr2) < 1e-12);
  // 중심화하지 않으면 X의 계수(W=0에서의 효과)가 달라지고 VIF 경고가 나온다
  assert.ok(Math.abs(table(on, '계수').rows[1].B - table(off, '계수').rows[1].B) > 0.1);
  assert.ok(off.warnings.some(w => w.includes('VIF') && w.includes('평균 중심화')));
  assert.ok(!off.blocks.some(b => b.type === 'chart') && !table(off, 'Johnson-Neyman'));
  assert.ok(table(off, '계수').columns[5].label.includes('90%'));
  // 단순기울기 값은 중심화 여부와 무관
  table(on, '조건부 효과').rows.forEach((r, k) => assert.ok(Math.abs(r.b - table(off, '조건부 효과').rows[k].b) < 1e-9));
});

test('runModeration: 입력 검증', () => {
  const run = req => M.runModeration(ds, Object.assign({ options: {} }, req));
  assert.throws(() => run({ dependent: [], independent: ['x'], moderator: ['w'] }), /변수를 선택/);
  assert.throws(() => run({ dependent: ['y'], independent: ['x', 'w'], moderator: ['g'] }), /독립변수는 하나만/);
  assert.throws(() => run({ dependent: ['y'], independent: ['x'], moderator: ['x'] }), /동시에 쓸 수 없습니다/);
  assert.throws(() => run({ dependent: ['y'], independent: ['x'], moderator: ['w'], covariates: ['w'] }), /동시에 쓸 수 없습니다/);
  assert.throws(() => run({ dependent: ['t'], independent: ['x'], moderator: ['w'] }), /숫자가 아닌 값/);
  assert.throws(() => run({ dependent: ['y'], independent: ['x'], moderator: ['w'], options: { conf: '100' } }), /신뢰수준/);
  const tiny = Data.buildDataset([['y', 'x', 'w'], [1, 2, 3], [2, 3, 5], [3, 5, 4], [5, 4, 9]], []);
  assert.throws(() => M.runModeration(tiny, { dependent: ['y'], independent: ['x'], moderator: ['w'], options: {} }), /사례 수/);
  const flat = Data.buildDataset([['y', 'x', 'w'], [1, 2, 3], [2, 2, 5], [3, 2, 4], [5, 2, 9], [4, 2, 1], [8, 2, 2]], []);
  assert.throws(() => M.runModeration(flat, { dependent: ['y'], independent: ['x'], moderator: ['w'], options: {} }), /모두 같아/);
});

test('선 그래프: 그림 SVG 와 비어 있을 때 안내', () => {
  const block = Ch.lineChartBlock('상호작용 그림', 'X', 'Y', [{ label: '낮음', points: [[0, 1], [10, 2]] }, { label: '높음', points: [[0, 1.5], [10, 6]] }]);
  assert.equal(block.kind, 'line');
  const html = Ch.renderChartHtml(block);
  assert.ok(html.includes('<svg') && html.includes('<polyline') && html.includes('낮음') && html.includes('높음'));
  assert.ok(Ch.renderChartHtml(Ch.lineChartBlock('빈', 'X', 'Y', [])).includes('그릴 데이터가 없습니다'));
  assert.equal((Ch.renderChartHtml(block).match(/<polyline/g) || []).length, 2);
});

/* ---------- 예제 데이터 ---------- */

test('예제 데이터: 근무연수 → 직무몰입이 자율성과 집단에 따라 달라지고, 성별에는 달라지지 않는다', () => {
  const sample = S.sampleDataset_();
  const parsed = Data.parseVarDefs([['변수명', '라벨', '척도', '값 레이블', '결측값 코드']].concat(sample.varRows));
  const full = Data.buildDataset([sample.header].concat(sample.rows), parsed.defs);
  const run = w => M.runModeration(full, { dependent: ['직무몰입'], independent: ['근무연수'], moderator: [w], options: {} });
  const bySig = w => table(run(w), '상호작용 검정').rows[0].p;
  assert.ok(bySig('직무자율성') < 0.001);
  assert.ok(bySig('집단') < 0.05);
  assert.ok(bySig('성별') > 0.05);
  const res = run('직무자율성');
  const jn = table(res, 'Johnson-Neyman');
  assert.ok(jn.rows.some(r => r.result.startsWith('유의함')));
  assert.deepEqual(table(res, '조건부 효과').rows.map(r => r.p < 0.05), [true, true, true]);
  // 직무자율성·직무몰입은 별도 난수열 (그 뒤에 열이 더 붙어도 값이 변하지 않는다)
  assert.deepEqual(sample.header.slice(20, 22), ['직무자율성', '직무몰입']);
});
