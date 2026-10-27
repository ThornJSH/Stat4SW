const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Ch = require('../src/Charts.js');
const R = require('../src/Render.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'charts.json'), 'utf8'));

function close(a, e, msg, tol = 1e-12) {
  assert.ok(Math.abs(a - e) <= tol * Math.max(1, Math.abs(e)), `${msg}: ${a} vs ${e}`);
}
const count = (s, re) => (s.match(re) || []).length;

/* ---------- 눈금 ---------- */

test('niceTicks_: 간격은 1·2·5×10^k, 구간을 덮고 하나 이상', () => {
  for (const [lo, hi] of [[0, 7], [0, 103], [-3.2, 4.8], [0.001, 0.0093], [1e5, 1e5 + 12], [5, 5]]) {
    const t = Ch.niceTicks_(lo, hi, 5);
    assert.ok(t.ticks.length >= 2, `${lo}~${hi}`);
    assert.ok(t.min <= lo + 1e-12 && t.max >= hi - 1e-12, `${lo}~${hi} 덮기`);
    const mant = t.step / Math.pow(10, Math.floor(Math.log10(t.step)));
    assert.ok([1, 2, 5].some(m => Math.abs(mant - m) < 1e-9), `${lo}~${hi} 간격 ${t.step}`);
  }
  assert.deepEqual(Ch.niceTicks_(0, 10, 5).ticks, [0, 2, 4, 6, 8, 10]);
});

/* ---------- 상자그림 통계 ---------- */

test('boxStats: 독립 공식(fivenum)과 정의대로 계산한 기준값 대조', () => {
  for (const c of ref.cases) {
    const r = Ch.boxStats(c.values), e = c.expected;
    assert.equal(r.n, e.n, c.name);
    close(r.q1, e.q1, `${c.name} q1`); close(r.median, e.median, `${c.name} median`); close(r.q3, e.q3, `${c.name} q3`);
    close(r.lo, e.lo, `${c.name} lo`); close(r.hi, e.hi, `${c.name} hi`);
    assert.deepEqual(r.outliers, e.outliers, `${c.name} 이상치`);
    assert.deepEqual(r.extremes, e.extremes, `${c.name} 극단값`);
  }
});

test('boxStats: 손으로 확인한 예 (1, 2, 4, 8, 100)', () => {
  const r = Ch.boxStats([100, 1, 8, 2, 4]);
  assert.equal(r.median, 4);
  assert.equal(r.q1, 2);              // 아래 절반 (중앙값 포함) (1, 2, 4) 의 중앙값
  assert.equal(r.q3, 8);              // 위 절반 (중앙값 포함) (4, 8, 100) 의 중앙값
  assert.deepEqual(r.extremes, [100]);
  assert.equal(r.hi, 8);
});

test('boxStats: 빈 배열, 값이 하나', () => {
  assert.equal(Ch.boxStats([]).median, null);
  const one = Ch.boxStats([5]);
  assert.deepEqual([one.q1, one.median, one.q3, one.lo, one.hi], [5, 5, 5, 5, 5]);
});

/* ---------- 히스토그램 구간 ---------- */

test('histogramBins: 사례 수 보존, 범위 덮기, 간격, 직접 센 값과 일치', () => {
  const sets = [
    Array.from({ length: 50 }, (_, i) => i * 1.7 + (i % 3)),
    Array.from({ length: 200 }, (_, i) => Math.sin(i) * 40 + 100),
    [1, 1, 2, 2, 2, 3, 9], [3.5, 3.5, 3.5], [0.001, 0.004, 0.0023, 0.0099],
    Array.from({ length: 1000 }, (_, i) => (i * 7919) % 1000 / 10)
  ];
  for (const vals of sets) {
    const bins = Ch.histogramBins(vals);
    assert.equal(bins.reduce((s, b) => s + b.count, 0), vals.length);
    assert.ok(bins[0].lo <= Math.min(...vals) + 1e-9 && bins[bins.length - 1].hi >= Math.max(...vals) - 1e-9);
    if (vals.length && Math.max(...vals) > Math.min(...vals)) {
      const step = bins[0].hi - bins[0].lo;
      const mant = step / Math.pow(10, Math.floor(Math.log10(step) + 1e-12));
      assert.ok([1, 2, 2.5, 5, 10].some(m => Math.abs(mant - m) < 1e-6), `간격 ${step}`);
      bins.forEach((b, i) => {
        let n = 0;
        for (const v of vals) if (v >= b.lo - 1e-12 && (v < b.hi - 1e-12 || (i === bins.length - 1 && v <= b.hi + 1e-12))) n++;
        assert.equal(b.count, n, `구간 ${i} [${b.lo}, ${b.hi})`);
      });
      for (let i = 1; i < bins.length; i++) assert.ok(Math.abs(bins[i].lo - bins[i - 1].hi) < 1e-9);   // 이어져 있다
    }
  }
  assert.equal(Ch.histogramBins([3.5, 3.5, 3.5]).length, 1);                 // 값이 모두 같으면 구간 하나
  assert.equal(Ch.histogramBins([3.5, 3.5, 3.5])[0].count, 3);
});

/* ---------- 차트 블록 ---------- */

test('scatterBlock: 회귀선 기울기·절편, 표시 점 수 제한', () => {
  const x = [1, 2, 3, 4, 5], y = x.map(v => 2 * v + 1);
  const b = Ch.scatterBlock('t', x, y, 'x', 'y', 1);
  close(b.line.slope, 2, 'slope'); close(b.line.intercept, 1, 'intercept');
  assert.equal(b.points.length, 5);
  const big = Array.from({ length: 10000 }, (_, i) => i);
  const bb = Ch.scatterBlock('t', big, big.map(v => v * 0.5), 'x', 'y', 1);
  assert.ok(bb.points.length <= 3000 && bb.n === 10000);
  assert.equal(Ch.scatterBlock('t', [2, 2, 2], [1, 2, 3], 'x', 'y', null).line, null);   // x 분산 0
});

test('histogramBlock, boxplotBlock: 요약 값과 구조', () => {
  const h = Ch.histogramBlock('h', [1, 2, 3, 4, 5], 'v');
  assert.equal(h.n, 5); close(h.mean, 3, 'mean'); close(h.sd, Math.sqrt(2.5), 'sd');
  const b = Ch.boxplotBlock('b', [{ label: 'A', values: [1, 2, 3, 4, 50] }, { label: 'B', values: [] }], 'y');
  assert.equal(b.groups.length, 2);
  assert.deepEqual(b.groups[0].extremes, [50]);
  assert.equal(b.groups[1].n, 0);
});

/* ---------- SVG 렌더링 ---------- */

test('renderChartHtml: 막대그래프 (막대 수, 값 레이블, 이스케이프, 회전)', () => {
  const html = Ch.renderChartHtml(Ch.barChartBlock('<b>제목</b>', ['남', '여', '<img src=x>'], [35, 25, 4], 'X & Y'));
  assert.ok(html.startsWith('<div class="chart"><svg') && html.endsWith('</svg></div>'));
  assert.equal(count(html, /<rect x="[\d.]+" y="[\d.]+" width="[\d.]+" height="[\d.]+" fill="#4e79a7"/g), 3);
  assert.ok(html.includes('>35<') && html.includes('X &amp; Y'));
  assert.ok(!html.includes('<b>') && !html.includes('<img') && html.includes('&lt;b&gt;제목'));
  const many = Ch.renderChartHtml(Ch.barChartBlock('c', Array.from({ length: 12 }, (_, i) => '범주' + i), Array.from({ length: 12 }, (_, i) => i + 1), ''));
  assert.ok(many.includes('rotate(-35'));
  assert.ok(!html.includes('NaN') && !many.includes('NaN'));
});

test('renderChartHtml: 히스토그램 (막대 수, 정규곡선, 요약 문구)', () => {
  const vals = Array.from({ length: 80 }, (_, i) => 50 + Math.sin(i * 1.3) * 12 + (i % 7));
  const block = Ch.histogramBlock('분포', vals, '점수');
  const html = Ch.renderChartHtml(block);
  assert.equal(count(html, /fill="#4e79a7" stroke="#ffffff"/g), block.bins.length);
  assert.ok(html.includes('<polyline') && html.includes('N = 80') && html.includes('평균 =') && html.includes('표준편차 ='));
  assert.ok(!html.includes('NaN') && !html.includes('Infinity'));
  const flat = Ch.renderChartHtml(Ch.histogramBlock('상수', [4, 4, 4], 'v'));    // sd = 0 → 곡선 없음
  assert.ok(!flat.includes('<polyline') && !flat.includes('NaN'));
});

test('renderChartHtml: 산점도 (점 수, 회귀선, r)', () => {
  const x = Array.from({ length: 30 }, (_, i) => i), y = x.map(v => v * 0.8 + (v % 4));
  const html = Ch.renderChartHtml(Ch.scatterBlock('산점도', x, y, 'x', 'y', 0.987));
  assert.equal(count(html, /<circle /g), 30);
  assert.ok(html.includes('stroke="#e15759"') && html.includes('r = 0.987') && html.includes('N = 30'));
  assert.ok(!html.includes('NaN'));
});

test('renderChartHtml: 상자그림 (상자 수, 이상치와 극단값 표시, 빈 집단)', () => {
  const block = Ch.boxplotBlock('상자그림', [
    { label: 'A', values: [1, 2, 3, 4, 5, 6, 7, 15, 60] },
    { label: 'B', values: [4, 5, 6, 7, 8] },
    { label: 'C', values: [] }
  ], '점수');
  const html = Ch.renderChartHtml(block);
  assert.equal(count(html, /fill="#cfe0f1"/g), 2);                 // 빈 집단은 상자가 없다
  assert.ok(count(html, /<circle /g) >= 1);                         // 15 가 이상치 (울타리 13 초과, 극단 경계 19 이하)
  assert.ok(html.includes('>*<'));                                  // 60 이 극단값
  assert.ok(html.includes('N=9') && html.includes('N=0') && !html.includes('NaN'));
});

test('renderChartHtml: 데이터가 없으면 안내 문구', () => {
  assert.ok(Ch.renderChartHtml(Ch.barChartBlock('빈 막대', [], [], '')).includes('그릴 데이터가 없습니다'));
  assert.ok(Ch.renderChartHtml(Ch.boxplotBlock('빈 상자', [{ label: 'A', values: [] }], 'y')).includes('그릴 데이터가 없습니다'));
});

test('Render.renderResultHtml: 차트 블록을 포함해 렌더링', () => {
  const html = R.renderResultHtml({
    analysis: 'x', title: 't', warnings: [],
    blocks: [Ch.barChartBlock('막대', ['a', 'b'], [3, 4], '')]
  });
  assert.ok(html.includes('<div class="chart"><svg'));
});

/* ---------- 분석의 차트 옵션 연결 ---------- */

const Data = require('../src/Data.js');
const Desc = require('../src/Stats_Descriptive.js');
const Rel = require('../src/Stats_Relation.js');
const Cmp = require('../src/Stats_Compare.js');
const An = require('../src/Stats_Anova.js');
const charts = res => res.blocks.filter(b => b.type === 'chart');

function chartData() {
  const defs = Data.parseVarDefs([['h'], ['g', '집단', '명목', '1=대조, 2=처치', '']]).defs;
  const rows = [['x', 'y', 'z', 'g']];
  for (let i = 0; i < 30; i++) rows.push([i % 7 + i * 0.1, (i * 5) % 11 + i * 0.3, (i * 3) % 13, 1 + (i % 2)]);
  rows.push(['', 4, 5, 1], [3, '', 5, 2]);
  return Data.buildDataset(rows, defs);
}

test('빈도분석: 막대그래프 옵션 (값 레이블, 빈도, 값이 너무 많으면 경고)', () => {
  const d = chartData();
  const without = Desc.runFrequencies(d, { vars: ['g'], options: {} });
  assert.equal(charts(without).length, 0);
  const res = Desc.runFrequencies(d, { vars: ['g'], options: { barChart: true } });
  const c = charts(res);
  assert.equal(c.length, 1);
  assert.equal(c[0].kind, 'bar');
  assert.deepEqual(c[0].labels, ['대조', '처치']);
  assert.deepEqual(c[0].values, [16, 16]);          // 반복 30건은 15/15, 추가한 2건이 각각 한 집단에 더해짐
  const many = Desc.runFrequencies(Data.buildDataset([['v'], ...Array.from({ length: 40 }, (_, i) => [i])]), { vars: ['v'], options: { barChart: true } });
  assert.equal(charts(many).length, 0);
  assert.ok(many.warnings.some(w => w.includes('막대그래프를 그리지 않았습니다')));
});

test('기술통계: 히스토그램 옵션 (변수마다, 유효값만)', () => {
  const d = chartData();
  const res = Desc.runDescriptives(d, { vars: ['x', 'y'], options: { histogram: true } });
  const c = charts(res);
  assert.equal(c.length, 2);
  assert.deepEqual(c.map(b => b.kind), ['histogram', 'histogram']);
  assert.equal(c[0].n, 31);                          // x: 32건 중 결측 1건 제외
  assert.equal(c[1].n, 31);                          // y: 32건 중 결측 1건 제외
  assert.equal(c[0].bins.reduce((s, b) => s + b.count, 0), 31);
  assert.equal(charts(Desc.runDescriptives(d, { vars: ['x'], options: {} })).length, 0);
  assert.ok(R.renderResultHtml(res).includes('<svg'));
});

test('상관분석: 산점도 옵션 (쌍별, 회귀선, r, 쌍이 많으면 10쌍까지)', () => {
  const d = chartData();
  const res = Rel.runCorrelation(d, { vars: ['x', 'y', 'z'], options: { pearson: true, scatter: true } });
  const c = charts(res);
  assert.equal(c.length, 3);
  assert.ok(c[0].caption.includes('x × y'));
  assert.ok(c[0].line !== null && typeof c[0].r === 'number');
  assert.equal(c[0].n, 30);                          // x와 y가 모두 유효한 사례
  const rows = [Array.from({ length: 6 }, (_, i) => 'v' + i)];
  for (let t = 0; t < 12; t++) rows.push(Array.from({ length: 6 }, (_, i) => (t * (i + 3)) % 7 + t * 0.1));
  const wide = Rel.runCorrelation(Data.buildDataset(rows), { vars: ['v0', 'v1', 'v2', 'v3', 'v4', 'v5'], options: { pearson: true, scatter: true } });
  assert.equal(charts(wide).length, 10);             // 15쌍 중 10쌍
  assert.ok(wide.warnings.some(w => w.includes('처음 10쌍')));
  assert.equal(charts(Rel.runCorrelation(d, { vars: ['x', 'y'], options: { pearson: true } })).length, 0);
});

test('독립표본 t검정, ANOVA: 상자그림 옵션 (집단마다 상자, 이상치)', () => {
  const d = chartData();
  const t = Cmp.runTTestIndependent(d, { dependent: ['x', 'y'], factor: ['g'], options: { boxplot: true } });
  const ct = charts(t);
  assert.equal(ct.length, 2);
  assert.deepEqual(ct[0].groups.map(g => g.label), ['대조', '처치']);
  assert.ok(ct[0].caption.includes('집단변수: 집단'));
  const a = An.runOneWayAnova(d, { dependent: ['x'], factor: ['g'], options: { boxplot: true, tukey: true } });
  assert.equal(charts(a).length, 1);
  assert.equal(a.blocks[a.blocks.length - 1].type, 'chart');      // 차트는 표들 뒤에
  assert.equal(charts(An.runOneWayAnova(d, { dependent: ['x'], factor: ['g'], options: {} })).length, 0);
  assert.equal(charts(Cmp.runTTestIndependent(d, { dependent: ['x'], factor: ['g'], options: {} })).length, 0);
  const html = R.renderResultHtml(a);
  assert.ok(html.includes('<svg') && html.includes('N='));
});
