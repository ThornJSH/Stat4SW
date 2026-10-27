const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const Data = require('../src/Data.js');
const TD = require('../src/TestData.js');
const Gen = require('./gen_test_datasets.js');

const SRC = path.join(__dirname, '..', 'src');
const RUN = {
  ipa: require('../src/Stats_Ipa.js').runIpa,
  cochran: require('../src/Stats_Cochran.js').runCochran,
  binom: require('../src/Stats_Binomial.js').runBinomial,
  anova: require('../src/Stats_Anova.js').runOneWayAnova,
  mannWhitney: require('../src/Stats_Nonparam.js').runMannWhitney,
  wilcoxonOne: require('../src/Stats_Nonparam.js').runWilcoxonOne,
  kruskalWallis: require('../src/Stats_Nonparam.js').runKruskalWallis,
  reliability: require('../src/Stats_Reliability.js').runReliability,
  moderation: require('../src/Stats_Moderation.js').runModeration,
  mediation: require('../src/Stats_Mediation.js').runMediation,
  logistic: require('../src/Stats_Logistic.js').runLogistic,
  distance: require('../src/Stats_Distance.js').runDistance,
  psm: require('../src/Stats_Psm.js').runPsm,
  kmeans: require('../src/Stats_Kmeans.js').runKmeans,
  cart: require('../src/Stats_Cart.js').runCart,
  mcnemar: require('../src/Stats_Mcnemar.js').runMcnemar,
  ttestInd: require('../src/Stats_Compare.js').runTTestIndependent,
  ttestPaired: require('../src/Stats_Compare.js').runTTestPaired,
  ttest1: require('../src/Stats_Compare.js').runTTestOne,
  friedman: require('../src/Stats_Nonparam.js').runFriedman,
  rmanova: require('../src/Stats_RepAnova.js').runRepeatedAnova,
  hierarchical: require('../src/Stats_Hierarchical.js').runHierarchical,
  efa: require('../src/Stats_Efa.js').runEfa,
  anova2: require('../src/Stats_Anova2.js').runAnova2,
  kappa: require('../src/Stats_Agreement.js').runKappa,
  icc: require('../src/Stats_Agreement.js').runIcc,
  ologit: require('../src/Stats_Ologit.js').runOrdinalLogit,
  correlation: require('../src/Stats_Relation.js').runCorrelation,
  km: require('../src/Stats_Survival.js').runKm
};
const HEADER = ['변수명', '라벨', '척도', '값 레이블', '결측값 코드'];
const blank = rows => rows.map(r => r.map(v => (v === null || v === undefined ? '' : v)));

function datasetOf(entry) {
  const parsed = Data.parseVarDefs([HEADER].concat(entry.varRows));
  return Data.buildDataset([entry.header].concat(blank(entry.rows)), parsed.defs);
}
function pick(result, spec) {
  const t = result.blocks.find(b => b.type === 'table' && b.caption.startsWith(spec.table));
  assert.ok(t, `표 "${spec.table}" 없음`);
  const row = t.rows[spec.row];
  assert.ok(row, `"${spec.table}"의 ${spec.row}번 행 없음`);
  return row[spec.key];
}

/* ---------- 레지스트리의 규칙 ---------- */

test('시험 데이터: client.html 의 모든 분석에 시험 데이터 항목이 있거나 기본 예제로 시험한다고 적혀 있다 (새 기법을 만들면 항목을 추가해야 한다)', () => {
  const html = fs.readFileSync(path.join(SRC, 'client.html'), 'utf8');
  const start = html.indexOf('var ANALYSES = {');
  assert.ok(start >= 0);
  const keys = [...html.slice(start).matchAll(/^    ([A-Za-z0-9]+): \{\n      label:/gm)].map(m => m[1]);
  assert.ok(keys.length >= 15, '분석 키를 찾지 못함');
  const withData = new Set(TD.TEST_DATASETS_.map(d => d.analysis));
  const missing = keys.filter(k => !withData.has(k) && TD.TEST_SAMPLE_COVERED_.indexOf(k) < 0 && TD.TEST_NO_DATA_.indexOf(k) < 0);
  assert.deepEqual(missing, [], '시험 데이터 항목이 없는 분석: ' + missing.join(', ') + ' → tests/gen_test_datasets.js 에 항목을 추가하고 node tests/gen_test_datasets.js 를 실행하세요');
  TD.TEST_SAMPLE_COVERED_.forEach(k => assert.ok(keys.indexOf(k) >= 0, '기본 예제 목록에 없는 분석 키: ' + k));
  TD.TEST_NO_DATA_.forEach(k => assert.ok(keys.indexOf(k) >= 0, '데이터 필요 없음 목록에 없는 분석 키: ' + k));
  withData.forEach(k => assert.ok(keys.indexOf(k) >= 0 && TD.TEST_ANALYSIS_LABELS_[k], '알 수 없는 분석: ' + k));
});

test('시험 데이터: src/TestData.js 가 생성기와 같다 (바뀌면 node tests/gen_test_datasets.js 로 다시 만든다)', () => {
  assert.equal(fs.readFileSync(path.join(SRC, 'TestData.js'), 'utf8'), Gen.build());
});

test('시험 데이터: 항목 구조 (id 유일, 머리글과 변수정의 일치, 요청의 변수가 자료에 있음, 안내와 기대값이 있음)', () => {
  const ids = new Set();
  TD.TEST_DATASETS_.forEach(d => {
    assert.ok(!ids.has(d.id), '중복 id ' + d.id); ids.add(d.id);
    assert.deepEqual(d.varRows.map(r => r[0]), d.header, d.id);
    d.rows.forEach(r => assert.equal(r.length, d.header.length, d.id));
    const req = d.request;
    ['dependent', 'independent', 'moderator', 'mediators', 'covariates', 'factors', 'factor', 'vars', 'block1', 'block2', 'block3', 'block4', 'factor1', 'factor2', 'controls', 'time', 'event', 'group', 'importance', 'performance'].forEach(k => {
      (req[k] || []).forEach(v => assert.ok(d.header.indexOf(v) >= 0, `${d.id}: 요청의 변수 ${v} 가 자료에 없음`));
    });
    assert.ok(d.guide.length >= 2 && d.expected.length >= 2 && d.title && d.description, d.id);
    assert.equal(req.analysis, d.analysis, d.id);
  });
});

/* ---------- 기대값(파이썬 기준값)과 이 도구의 계산이 같다 ---------- */

for (const entry of TD.TEST_DATASETS_) {
  test(`시험 데이터 ${entry.id}: 파이썬 기준값과 이 도구의 결과가 일치한다 (${entry.title})`, () => {
    const ds = datasetOf(entry);
    const res = RUN[entry.analysis](ds, JSON.parse(JSON.stringify(entry.request)));
    let checked = 0;
    entry.expected.forEach(spec => {
      if (spec.info) { assert.ok(pick(res, spec) !== undefined, spec.label); return; }
      let got = pick(res, spec);
      let want = spec.value;
      if (spec.abs) { got = Math.abs(got); want = Math.abs(want); }
      assert.equal(typeof got, 'number', `${entry.id} ${spec.label}: ${got}`);
      assert.ok(Math.abs(got - want) <= spec.tol * Math.max(1e-9, Math.abs(want)), `${entry.id} ${spec.label}: 도구 ${got} vs 기준값 ${want}`);
      checked++;
    });
    assert.ok(checked >= 2, entry.id);
  });
}
