const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const Data = require('../src/Data.js');
const S = require('../src/Sample.js');
const A = require('../src/Stats_Assumptions.js');
const C = require('../src/Stats_Compare.js');
const V = require('../src/Stats_Anova.js');
const NP = require('../src/Stats_Nonparam.js');
const Rel = require('../src/Stats_Relation.js');
const Reg = require('../src/Stats_Regression.js');
const Desc = require('../src/Stats_Descriptive.js');
const K = require('../src/Compute.js');

const HEADER = ['변수명', '라벨', '척도', '값 레이블', '결측값 코드'];
const sample = S.sampleDataset_();
const parsed = Data.parseVarDefs([HEADER].concat(sample.varRows));
const ds = Data.buildDataset([sample.header].concat(sample.rows), parsed.defs);
const col = name => ds.data[name].filter(v => v !== null);
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const texts = res => res.blocks.filter(b => b.type === 'text').map(b => b.text);

function pearson(a, b) {
  const x = [], y = [];
  ds.data[a].forEach((v, i) => { const w = ds.data[b][i]; if (v !== null && w !== null) { x.push(v); y.push(w); } });
  return Rel.pearsonR(x, y).r;
}
function groups(dep, factor) {
  const g = {};
  ds.data[factor].forEach((k, i) => { const v = ds.data[dep][i]; if (v !== null && k !== null) (g[k] = g[k] || []).push(v); });
  return Object.keys(g).sort().map(k => g[k]);
}

/* ---------- 구조 ---------- */

test('예제 데이터: 60명, 27개 변수, 이름이 겹치지 않고 변수정의가 경고 없이 읽힌다', () => {
  assert.equal(ds.n, 60);
  assert.equal(sample.header.length, 27);
  assert.equal(new Set(sample.header).size, 27);
  assert.ok(sample.header.length <= 40);                   // 기본 26열을 넘으면 예제 채우기가 열을 늘린다
  assert.deepEqual(parsed.warnings, []);
  assert.deepEqual(ds.warnings, []);
  assert.deepEqual(sample.varRows.map(r => r[0]), sample.header);
  sample.rows.forEach(r => assert.equal(r.length, 27));
  // 항상 같은 값 (고정 난수)
  assert.deepEqual(S.sampleDataset_().rows, sample.rows);
});

test('예제 데이터: 앞의 다섯 열은 처음 예제와 값이 똑같다 (확인사항의 기준 결과가 의존)', () => {
  let seed = 20260101;
  const rnd = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };
  for (let i = 0; i < 60; i++) {
    const g = 1 + Math.floor(rnd() * 3);
    const s = Math.round(60 + g * 4 + rnd() * 30);
    const legacy = [1 + Math.floor(rnd() * 2), g, 20 + Math.floor(rnd() * 30), rnd() < 0.05 ? 999 : s, rnd() < 0.08 ? 99 : 1 + Math.floor(rnd() * 5)];
    assert.deepEqual(sample.rows[i].slice(0, 5), legacy, `행 ${i + 2}`);
  }
  const sw = A.shapiroWilk(col('점수'));                    // 사용자가 시트에서 확인한 값
  assert.equal(sw.n, 57);
  assert.ok(Math.abs(sw.w - 0.97179) < 1e-5 && Math.abs(sw.p - 0.20308) < 1e-4);
});

test('예제 데이터: 결측값 코드와 빈 칸이 결측으로 처리되고, 척도·글자 변수가 의도대로 인식된다', () => {
  const by = Object.fromEntries(ds.vars.map(v => [v.name, v]));
  assert.ok(sample.rows.some(r => r[3] === 999) && sample.rows.some(r => r[4] === 99));
  assert.ok(!col('점수').includes(999) && !col('만족도').includes(99) && !col('월소득').includes(9999));
  assert.ok(col('사전점수').length < 60 && col('사후점수').length < 60);            // 빈 칸
  assert.equal(by['지역'].scale, 'nominal');
  assert.ok(col('지역').every(v => typeof v === 'string'));
  assert.deepEqual([by['학력'].scale, by['문항1'].scale, by['합격'].scale, by['월소득'].scale], ['ordinal', 'ordinal', 'nominal', 'scale']);
  assert.equal(by['문항3'].valueLabels['5'], '매우 그렇다');
  [1, 2, 3, 4].forEach(k => assert.ok(col('학력').includes(k), `학력 ${k}`));
  ['서울', '경기', '부산', '기타'].forEach(r => assert.ok(col('지역').filter(v => v === r).length >= 5, r));
  const pass = col('합격').filter(v => v === 1).length;
  assert.ok(pass >= 20 && pass <= 40, `합격 ${pass}`);
});

/* ---------- 변수마다 의도한 특성 ---------- */

test('예제 데이터: 월소득은 치우쳐서 정규성이 기각되고 비모수·스피어만 추천이 나온다', () => {
  assert.ok(A.shapiroWilk(col('월소득')).p < 0.001);
  const res = C.runTTestIndependent(ds, { dependent: ['월소득'], factor: ['성별'], options: { assumptions: true } });
  assert.ok(texts(res).some(t => t.includes('Mann-Whitney U 검정')), texts(res).join('\n'));
  const cor = Rel.runCorrelation(ds, { vars: ['연령', '월소득'], options: { pearson: true, assumptions: true } });
  assert.ok(texts(cor).some(t => t.includes('스피어만')), texts(cor).join('\n'));
});

test('예제 데이터: 불안점수는 성별 분산이 달라 "등분산을 가정하지 않음" 행이, 점수는 "가정함" 행이 굵게', () => {
  assert.ok(A.leveneTest(groups('불안점수', '성별'), false).p < 0.01);
  const bold = res => table(res, '독립표본 검정').rows.filter(r => r._bold).map(r => r.assume);
  const anx = C.runTTestIndependent(ds, { dependent: ['불안점수'], factor: ['성별'], options: {} });
  assert.deepEqual(bold(anx), ['등분산을 가정하지 않음']);
  const sc = C.runTTestIndependent(ds, { dependent: ['점수'], factor: ['성별'], options: {} });
  assert.deepEqual(bold(sc), ['등분산을 가정함']);
});

test('예제 데이터: 스트레스는 집단 분산이 달라 분산분석 Levene이 기각되고 Welch를 안내', () => {
  assert.ok(A.leveneTest(groups('스트레스', '집단'), false).p < 0.01);
  const res = V.runOneWayAnova(ds, { dependent: ['스트레스'], factor: ['집단'], options: { assumptions: true, levene: true, welch: true } });
  assert.ok(table(res, '등분산 검정').rows[0].p < 0.01);
  assert.ok(texts(res).some(t => t.includes('Welch')), texts(res).join('\n'));
  const ok = V.runOneWayAnova(ds, { dependent: ['점수'], factor: ['집단'], options: { levene: true } });
  assert.ok(table(ok, '등분산 검정').rows[0].p > 0.05);
});

test('예제 데이터: 사전·사후점수는 집단마다 향상 폭이 다르고, 대응표본 검정이 유의하다', () => {
  const diffByGroup = [0, 0, 0], n = [0, 0, 0];
  ds.data['집단'].forEach((g, i) => {
    const a = ds.data['사전점수'][i], b = ds.data['사후점수'][i];
    if (a !== null && b !== null) { diffByGroup[g - 1] += b - a; n[g - 1]++; }
  });
  const m = diffByGroup.map((s, i) => s / n[i]);
  assert.ok(m[0] < m[1] && m[1] < m[2], m.join(', '));
  assert.ok(m[2] - m[0] > 4);
  const res = C.runTTestPaired(ds, { first: ['사후점수'], second: ['사전점수'], options: { assumptions: true } });
  assert.ok(table(res, '대응표본 검정').rows[0].p2 < 0.001);
  assert.ok(pearson('사전점수', '사후점수') > 0.7);
  const w = NP.runWilcoxonPaired(ds, { first: ['사후점수'], second: ['사전점수'], options: { effect: true } });
  assert.ok(table(w, 'Wilcoxon').rows[0].p < 0.001);
});

test('예제 데이터: 연령과 근무연수는 거의 같이 움직이고, 성과점수는 회귀로 잘 설명된다', () => {
  assert.ok(pearson('연령', '근무연수') > 0.9);
  const res = Reg.runRegression(ds, { dependent: ['성과점수'], covariates: ['근무연수', '학력'], factors: ['집단'], options: { collinearity: true } });
  const summary = table(res, '모형 요약');
  assert.ok(summary.rows[0].r2 > 0.6, `R² ${summary.rows[0].r2}`);
  const coll = Reg.runRegression(ds, { dependent: ['성과점수'], covariates: ['연령', '근무연수'], options: { collinearity: true } });
  const coeffs = table(coll, '계수');
  assert.ok(coeffs.rows.some(r => r.vif > 5), '연령과 근무연수를 함께 넣으면 VIF가 커야 한다');
});

test('예제 데이터: 문항3·5는 역문항 — 뒤집기 전에는 음의 상관, 뒤집은 뒤에는 양의 상관이고 α가 높다', () => {
  assert.ok(pearson('문항1', '문항3') < -0.3 && pearson('문항1', '문항5') < -0.3);
  assert.ok(pearson('문항1', '문항2') > 0.3 && pearson('문항1', '문항4') > 0.3);
  const plan = K.planCompute_(ds, { kind: 'reverse', vars: ['문항3', '문항5'], min: 1, max: 5, suffix: '_r' });
  assert.deepEqual(plan.columns.map(c => c.name), ['문항3_r', '문항5_r']);
  // 새 열을 데이터셋에 더해 상관과 크론바흐 α를 계산한다
  const items = ['문항1', '문항2', '문항4'].map(n => ds.data[n]).concat(plan.columns.map(c => c.values));
  const valid = ds.data['문항1'].map((_, i) => items.every(c => c[i] !== null));
  const cols = items.map(c => c.filter((_, i) => valid[i]));
  const k = cols.length, n = cols[0].length;
  const variance = c => { const m = c.reduce((s, v) => s + v, 0) / c.length; return c.reduce((s, v) => s + (v - m) * (v - m), 0) / (c.length - 1); };
  const totals = Array.from({ length: n }, (_, i) => cols.reduce((s, c) => s + c[i], 0));
  const alpha = k / (k - 1) * (1 - cols.reduce((s, c) => s + variance(c), 0) / variance(totals));
  assert.ok(alpha > 0.7, `α ${alpha}`);
  assert.ok(n >= 40, `완전한 사례 ${n}`);
  const mean = K.planCompute_(ds, { kind: 'mean', vars: ['문항1', '문항2', '문항4'], name: '만족척도', minValid: 2 });
  assert.ok(mean.columns[0].valid >= 55);
});

/* ---------- 모든 분석이 새 변수로 돌아간다 ---------- */

test('예제 데이터: 지금까지의 모든 분석이 새 변수로 오류 없이 실행된다', () => {
  const run = (fn, req) => { const res = fn(ds, req); assert.ok(res.blocks.length > 0, res.title); return res; };
  run(Desc.runFrequencies, { vars: ['지역', '학력', '합격'], options: { barChart: true } });
  run(Desc.runDescriptives, { vars: ['월소득', '근무연수', '성과점수'], options: { skewness: true, kurtosis: true, histogram: true } });
  run(A.runNormality, { vars: ['월소득', '성과점수', '불안점수'] });
  run(A.runLevene, { dependent: ['불안점수', '스트레스'], factor: ['성별'], options: {} });
  run(C.runTTestOne, { vars: ['성과점수'], options: { testValue: '70', assumptions: true } });
  run(C.runTTestPaired, { first: ['사후점수'], second: ['사전점수'], options: { assumptions: true } });
  run(V.runOneWayAnova, { dependent: ['성과점수', '스트레스'], factor: ['지역'], options: { assumptions: true, levene: true, effect: true, welch: true, bonferroni: true, boxplot: true } });
  run(NP.runMannWhitney, { dependent: ['월소득', '불안점수'], factor: ['성별'], options: { effect: true, boxplot: true } });
  run(NP.runWilcoxonOne, { vars: ['월소득'], options: { testValue: '300', effect: true } });
  const kw = run(NP.runKruskalWallis, { dependent: ['월소득'], factor: ['지역'], options: { dunn: true, effect: true } });
  assert.equal(table(kw, '다중비교: Dunn').rows.length, 6);                 // 4개 지역 → 6쌍
  run(Rel.runCorrelation, { vars: ['연령', '근무연수', '월소득', '성과점수', '문항1'], options: { pearson: true, spearman: true, ci: true, assumptions: true } });
  const ct = run(Rel.runCrosstab, { rowvar: ['합격'], colvar: ['집단'], options: { expected: true, rowPct: true, effect: true } });
  assert.ok(ct.blocks.length > 1);
  run(Rel.runCrosstab, { rowvar: ['지역'], colvar: ['합격'], options: { rowPct: true } });
  run(Reg.runRegression, { dependent: ['성과점수'], covariates: ['근무연수', '문항1'], factors: ['지역', '집단'], options: { ci: true, collinearity: true, dw: true } });
});

/* ---------- CSV 파일 ---------- */

test('예제 데이터 CSV 파일(doc/예제데이터)이 Sample.js 와 같다 (바뀌면 node tests/gen_sample_csv.js 로 다시 만든다)', () => {
  const { buildSampleCsv } = require('./gen_sample_csv.js');
  const csv = buildSampleCsv();
  const dir = path.join(__dirname, '..', 'doc', '예제데이터');
  assert.equal(fs.readFileSync(path.join(dir, '예제_데이터.csv'), 'utf8'), csv.data);
  assert.equal(fs.readFileSync(path.join(dir, '예제_변수정의.csv'), 'utf8'), csv.vars);
  assert.equal(csv.data.split('\r\n').filter(Boolean).length, 61);          // 머리글 + 60명
  assert.equal(csv.vars.split('\r\n').filter(Boolean).length, 28);          // 머리글 + 27개 변수
  assert.ok(csv.vars.includes('"1=남, 2=여"'));                              // 쉼표가 든 값 레이블은 따옴표로 묶는다
});

/* ---------- 예제 데이터 채우기 (시트 쓰기) ---------- */

test('시뮬레이션: 예제 데이터 채우기가 데이터 탭과 변수정의 탭에 모든 열과 정의를 쓴다', () => {
  const written = { data: null, vars: null, calls: [] };
  const sheetMock = key => ({
    getMaxRows: () => 1000, setFrozenRows() {}, insertColumnsAfter() { written.calls.push(key + ':insertColumns'); }, clearContents() { written.calls.push(key + ':clear'); },
    getRange: (r, c, nr, nc) => {
      const rng = {
        clearContent() { written.calls.push(key + ':clearContent'); return rng; },
        setNumberFormat(f) { written.calls.push(key + ':format:' + f + ':' + c); return rng; },
        setValues(v) { written[key] = { r, c, nr, nc, v }; if (v.length !== nr || v.some(a => a.length !== nc)) throw new Error('setValues 크기 불일치'); return rng; },
        setFontWeight() { return rng; }, setBackground() { return rng; }, setDataValidation() { return rng; }
      };
      return rng;
    }
  });
  const sheets = { '데이터': sheetMock('data'), '변수정의': sheetMock('vars') };
  const env = {
    console: { log() {}, error() {} }, Utilities: {}, Session: {}, PropertiesService: {},
    SpreadsheetApp: {
      getUi: () => ({ alert: () => 'YES', ButtonSet: { YES_NO: 'YES_NO' }, Button: { YES: 'YES' } }),
      getActiveSpreadsheet: () => ({ getSheetByName: n => sheets[n] }),
      newDataValidation: () => { const b = { requireValueInRange: () => b, requireValueInList: () => b, setAllowInvalid: () => b, setHelpText: () => b, build: () => ({}) }; return b; }
    }
  };
  sheets['데이터'].getMaxColumns = () => 26;
  const ctx = vm.createContext(env);
  const SRC = path.join(__dirname, '..', 'src');
  for (const f of fs.readdirSync(SRC).filter(x => x.endsWith('.js')).sort()) vm.runInContext(fs.readFileSync(path.join(SRC, f), 'utf8'), ctx, { filename: f });
  vm.runInContext('fillSampleData()', ctx);
  assert.equal(written.data.nr, 61);
  assert.equal(written.data.nc, 27);
  assert.deepEqual(JSON.parse(JSON.stringify(written.data.v[0])), sample.header);
  assert.deepEqual(JSON.parse(JSON.stringify(written.data.v.slice(1))), sample.rows);
  assert.equal(written.vars.nr, 27);
  assert.equal(written.vars.nc, 5);
  assert.deepEqual(JSON.parse(JSON.stringify(written.vars.v)), sample.varRows);
  assert.ok(written.calls.includes('vars:format:@:4'));
});
