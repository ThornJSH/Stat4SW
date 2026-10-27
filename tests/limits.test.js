const test = require('node:test');
const assert = require('node:assert/strict');
const Data = require('../src/Data.js');
const Ch = require('../src/Charts.js');
const Desc = require('../src/Stats_Descriptive.js');
const An = require('../src/Stats_Anova.js');
const Cmp = require('../src/Stats_Compare.js');
const Rel = require('../src/Stats_Relation.js');
const Reg = require('../src/Stats_Regression.js');
const R = require('../src/Render.js');

/**
 * 사례 수가 많을 때의 안정성. 배열이 클 때 터지는 코드(Math.min.apply 등)와
 * 큰 n에서의 수치 정확도를 확인한다. 한계 값은 Config.js 의 MAX_ROWS 와 함께 정한 근거다.
 */

let seed = 987654;
const rnd = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };
const gauss = () => (rnd() + rnd() + rnd() + rnd() - 2) * 1.7;

function bigDataset(n) {
  const rows = [['g3', 'g2', 'x1', 'x2', 'x3', 'y']];
  for (let i = 0; i < n; i++) {
    const g3 = 1 + Math.floor(rnd() * 3), g2 = 1 + Math.floor(rnd() * 2);
    const x1 = gauss() * 10 + 50, x2 = gauss() * 5 + 20, x3 = gauss() * 3;
    rows.push([g3, g2, x1, x2, x3, 0.3 * x1 - 0.2 * x2 + g3 * 2 + gauss() * 5]);
  }
  return Data.buildDataset(rows, []);
}

const N = 200000;     // Math.min.apply 가 실패하던 크기(약 12만)보다 크게
const ds = bigDataset(N);

test('큰 배열: 히스토그램 구간, 상자그림, 집단 요약이 호출 스택 오류 없이 계산된다', () => {
  const x = ds.data.x1;
  assert.doesNotThrow(() => Ch.histogramBins(x));
  assert.equal(Ch.histogramBins(x).reduce((s, b) => s + b.count, 0), N);
  assert.doesNotThrow(() => Ch.boxStats(x));
  const block = Ch.histogramBlock('h', x, 'x1');
  assert.equal(block.n, N);
  assert.ok(R.renderResultHtml({ analysis: 'x', title: 't', warnings: [], blocks: [block] }).includes('<svg'));
});

test('큰 n: 기술통계(+히스토그램), ANOVA(+Tukey·상자그림), t검정, 상관, 회귀가 끝까지 계산된다', () => {
  assert.doesNotThrow(() => Desc.runDescriptives(ds, { vars: ['x1', 'x2'], options: { skewness: true, kurtosis: true, histogram: true } }));
  const an = An.runOneWayAnova(ds, { dependent: ['y'], factor: ['g3'], options: { levene: true, effect: true, welch: true, boxplot: true } });
  assert.ok(an.blocks.some(b => b.type === 'chart'));
  assert.doesNotThrow(() => Cmp.runTTestIndependent(ds, { dependent: ['y'], factor: ['g2'], options: { boxplot: true } }));
  assert.doesNotThrow(() => Rel.runCorrelation(ds, { vars: ['x1', 'x2', 'x3'], options: { pearson: true, spearman: true, scatter: true } }));
  assert.doesNotThrow(() => Rel.runCrosstab(ds, { rowvar: ['g3'], colvar: ['g2'], options: { expected: true } }));
  const reg = Reg.runRegression(ds, { dependent: ['y'], covariates: ['x1', 'x2', 'x3'], factors: ['g3'], options: { collinearity: true } });
  const coef = reg.blocks.find(b => b.caption === '계수').rows;
  assert.ok(Math.abs(coef.find(r => r.name === 'x1').B - 0.3) < 0.01);                // 생성한 모형의 계수 0.3 을 되찾는다
  assert.ok(Math.abs(coef.find(r => r.name === 'x2').B + 0.2) < 0.01);
});

test('큰 n: 평균·분산·왜도·첨도가 정확 산술과 같은 자릿수로 맞는다 (200000건, 평균 1e6)', () => {
  // 평균이 큰 자료에서도 보정된 중심 적률이 안정적이어야 한다. 정확 산술(BigInt 분수)과 비교한다.
  const m = 200000, vals = [];
  for (let i = 0; i < m; i++) vals.push(1e6 + gauss());
  const got = Desc.describeColumn(vals);
  // 정확한 값: 입력 double 을 2^-30 단위의 정수로 바꿔 BigInt 로 계산 (이 입력은 2^-30 배수가 아니므로 먼저 반올림한 값을 입력으로 쓴다)
  const SCALE = 2n ** 30n;
  const q = vals.map(v => BigInt(Math.round(v * 2 ** 30)));
  const exactVals = q.map(b => Number(b) / 2 ** 30);                                 // 위 정수가 나타내는 값 (double 로 정확히 표현됨)
  const got2 = Desc.describeColumn(exactVals);
  const n = BigInt(m);
  const sum = q.reduce((a, b) => a + b, 0n);
  const meanNum = sum, meanDen = n;                                                    // 평균 = sum / n (단위 SCALE)
  // Σ(x - mean)^2 = Σ x^2 - sum^2 / n  → (n Σ x^2 - sum^2) / n  (단위 SCALE^2)
  const sumSq = q.reduce((a, b) => a + b * b, 0n);
  const m2num = n * sumSq - sum * sum;                                                 // 단위 SCALE^2, 분모 n
  const varianceExact = Number(m2num * 10n ** 30n / (n * (n - 1n) * SCALE * SCALE)) / 1e30;
  const meanExact = Number(meanNum * 10n ** 30n / (meanDen * SCALE)) / 1e30;
  assert.ok(Math.abs(got2.mean - meanExact) / meanExact < 1e-14, 'mean ' + got2.mean + ' vs ' + meanExact);
  assert.ok(Math.abs(got2.variance - varianceExact) / varianceExact < 1e-9, 'variance ' + got2.variance + ' vs ' + varianceExact);
  // 균등분포 4개의 합: 왜도 0, 초과첨도의 이론값은 -6/(5·4) = -0.3 (표본 20만 건이면 소수 둘째 자리까지 맞는다)
  assert.ok(Math.abs(got.skewness) < 0.02, 'skew ' + got.skewness);
  assert.ok(Math.abs(got.kurtosis + 0.3) < 0.03, 'kurt ' + got.kurtosis);
});

test('큰 문제: 요인분석 변수 50개 × 사례 5,000건(주축분해, 베리맥스, 요인 6개)이 Node 에서 10초 안에 끝난다', () => {
  const p = 50, n = 5000, k = 6;
  const rows = [Array.from({ length: p }, (_, j) => 'v' + j)];
  for (let i = 0; i < n; i++) {
    const f = Array.from({ length: k }, gauss);
    rows.push(Array.from({ length: p }, (_, j) => Math.round((f[j % k] * (0.5 + (j % 5) * 0.1) + gauss() * 0.8) * 1000) / 1000));
  }
  const big = Data.buildDataset(rows, []);
  const Efa = require('../src/Stats_Efa.js');
  const t0 = Date.now();
  const res = Efa.runEfa(big, { vars: rows[0], options: { method: 'paf', rotation: 'varimax' } });
  const sec = (Date.now() - t0) / 1000;
  const comm = res.blocks.find(b => b.caption === '공통성').rows;
  assert.equal(comm.length, p);
  assert.ok(res.blocks.some(b => b.caption === '회전된 요인 행렬'));
  assert.ok(sec < 10, `요인분석 ${sec}초`);
  const t1 = Date.now();
  Efa.runEfa(big, { vars: rows[0], options: { method: 'pca', rotation: 'varimax', nfactors: '8' } });
  assert.ok((Date.now() - t1) / 1000 < 10);
});

test('큰 문제: 요인분석 프로맥스 + 병렬분석 + 점수 계산(변수 50개 × 사례 5,000건)이 Node 에서 10초 안에 끝난다', () => {
  const p = 50, n = 5000, k = 6;
  const rows = [Array.from({ length: p }, (_, j) => 'v' + j)];
  for (let i = 0; i < n; i++) {
    const f = Array.from({ length: k }, gauss);
    rows.push(Array.from({ length: p }, (_, j) => Math.round((f[j % k] * (0.5 + (j % 5) * 0.1) + gauss() * 0.8) * 1000) / 1000));
  }
  const big = Data.buildDataset(rows, []);
  const t0 = Date.now();
  const res = require('../src/Stats_Efa.js').runEfa(big, { vars: rows[0], options: { method: 'paf', rotation: 'promax', parallel: true, save: true } });
  const sec = (Date.now() - t0) / 1000;
  assert.ok(res.blocks.some(b => b.caption === '패턴 행렬') && res.blocks.some(b => b.caption === '병렬분석'));
  assert.equal(res.saveColumns.length, 6);
  assert.equal(res.saveColumns[0].values.length, n);
  assert.ok(sec < 10, `요인분석(프로맥스+병렬분석+점수) ${sec}초`);
});

test('큰 문제: 이원 분산분석 10×10 수준 × 사례 100,000건(유형 III, 사후검정·단순 주효과·주변평균)이 Node 에서 10초 안에 끝난다', () => {
  const n = 100000;
  const rows = [['y', 'A', 'B']];
  for (let i = 0; i < n; i++) {
    const a = 1 + Math.floor(rnd() * 10), b = 1 + Math.floor(rnd() * 10);
    rows.push([a * 0.3 + b * 0.2 + (a === 3 && b === 4 ? 2 : 0) + gauss() * 5, a, b]);
  }
  const big = Data.buildDataset(rows, []);
  const t0 = Date.now();
  const res = require('../src/Stats_Anova2.js').runAnova2(big, { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: { emmeans: true } });
  const sec = (Date.now() - t0) / 1000;
  const tt = res.blocks.find(b => b.caption.startsWith('개체-간 효과 검정')).rows;
  assert.ok(tt.find(r => r.source === 'A').F > 100);
  assert.ok(res.blocks.some(b => b.caption.startsWith('사후검정: A')));
  assert.ok(sec < 10, `이원 분산분석 ${sec}초`);
});

test('큰 문제: 혼합 설계 분산분석 10집단 × 12회 측정 × 5만 명(단순 효과·쌍별 비교·Levene 포함)이 Node 에서 10초 안에 끝난다', () => {
  const n = 50000, k = 12, g = 10;
  const names = Array.from({ length: k }, (_, j) => 't' + j);
  const rows = [names.concat(['g'])];
  for (let i = 0; i < n; i++) {
    const s = gauss(), gg = 1 + Math.floor(rnd() * g);
    rows.push(names.map((_, j) => s + j * 0.01 * gg + gauss()).concat([gg]));
  }
  const big = Data.buildDataset(rows, [{ name: 'g', label: '', scale: 'nominal', valueLabels: {}, missing: [] }]);
  const t0 = Date.now();
  const res = require('../src/Stats_RepAnova.js').runRepeatedAnova(big, { vars: names, factor: ['g'], options: { posthoc: true, plot: true, simple: true, levene: true } });
  const sec = (Date.now() - t0) / 1000;
  assert.ok(res.blocks.some(b => b.caption === '쌍별 비교 (집단별, Bonferroni 보정)'));
  assert.ok(sec < 10, `혼합 설계 ${sec}초`);
});

test('큰 문제: 서열 로지스틱 10범주 × 독립변수 10개 × 5만 건(Brant·분류표·저장 포함)이 Node 에서 10초 안에 끝난다', () => {
  const n = 50000, k = 10, J = 10;
  const names = Array.from({ length: k }, (_, i) => 'x' + i);
  const rows = [['y'].concat(names)];
  const lg = () => { const u = rnd() || 1e-9; return Math.log(u / (1 - u)); };
  for (let i = 0; i < n; i++) {
    const x = names.map(() => (rnd() + rnd() + rnd() - 1.5) * 2);
    let lat = lg();
    x.forEach((v, j) => { lat += v * (0.3 + 0.05 * j); });
    let y = 1;
    [-3, -2, -1, -0.3, 0.4, 1.1, 1.9, 2.8, 3.8].forEach(c => { if (lat > c) y++; });
    rows.push([y].concat(x));
  }
  const big = Data.buildDataset(rows, [{ name: 'y', label: '', scale: 'ordinal', valueLabels: {}, missing: [] }]);
  const t0 = Date.now();
  const res = require('../src/Stats_Ologit.js').runOrdinalLogit(big, { dependent: ['y'], covariates: names, options: { parallel: true, classify: true, save: true } });
  const sec = (Date.now() - t0) / 1000;
  assert.ok(res.blocks.some(b => b.caption.startsWith('비례오즈 가정 점검')) && res.saveColumns.length === 1 + J);
  assert.ok(sec < 10, `서열 로지스틱 ${sec}초`);
});

test('큰 문제: Kaplan–Meier 10집단 × 5만 건(고유 시간이 거의 5만 개, 전체 표·쌍별 45쌍·그림 포함)이 Node 에서 10초 안에 끝난다', () => {
  const n = 50000, G = 10;
  const rows = [['t', 'e', 'g']];
  for (let i = 0; i < n; i++) {
    const g = i % G;
    const tt = -Math.log(rnd() || 1e-9) * (10 + g);
    const cc = rnd() * 40;
    rows.push([Math.min(tt, cc), tt <= cc ? 1 : 0, g + 1]);
  }
  const def = (name, scale) => ({ name, label: '', scale, valueLabels: {}, missing: [] });
  const big = Data.buildDataset(rows, [def('t', 'scale'), def('e', 'nominal'), def('g', 'nominal')]);
  const t0 = Date.now();
  const res = require('../src/Stats_Survival.js').runKm(big, { time: ['t'], event: ['e'], group: ['g'], options: { eventValue: '1', tests: true, pairs: true, plot: true, table: true } });
  const sec = (Date.now() - t0) / 1000;
  assert.ok(res.blocks.some(b => b.caption.startsWith('집단 쌍별')) && res.blocks.filter(b => b.caption && b.caption.startsWith('사건 시점별 생존표')).every(b => b.rows.length === 500));
  assert.ok(res.blocks.find(b => b.type === 'chart').series.every(s => s.points.length <= 2 * 240 + 2));
  assert.ok(sec < 10, `Kaplan–Meier ${sec}초`);
});
