/**
 * Apps Script 환경 흉내: src/*.js 를 module 없이 한 전역 범위(vm 컨텍스트)에 불러온다.
 * - 파일 사이의 전역 함수 이름이 겹치거나(서로 덮어쓰기) 로드 순서에 의존하면 여기서 실패한다
 * - Api.js 의 runAnalysis 가 모든 분석에 대해 끝까지 동작하는지, 결과가 JSON 으로 오갈 수 있는지 확인한다
 * 실제 SpreadsheetApp, PropertiesService 는 mock 이므로 Apps Script 자체의 동작 검증은 아니다.
 */
const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const Safety = require('../src/Safety.js');

const SRC = path.join(__dirname, '..', 'src');

/* ---------- 예제 데이터 (Template.js 의 예제와 같은 생성 방식) ---------- */

function sampleSheets() {
  let seed = 20260101;
  const rnd = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };
  const data = [['성별', '집단', '연령', '점수', '만족도']];
  for (let i = 0; i < 60; i++) {
    const g = 1 + Math.floor(rnd() * 3);
    const s = Math.round(60 + g * 4 + rnd() * 30);
    data.push([1 + Math.floor(rnd() * 2), g, 20 + Math.floor(rnd() * 30), rnd() < 0.05 ? 999 : s, rnd() < 0.08 ? 99 : 1 + Math.floor(rnd() * 5)]);
  }
  const vars = [['변수명', '라벨', '척도', '값 레이블', '결측값 코드'],
    ['성별', '성별', '명목', '1=남, 2=여', ''], ['집단', '실험 집단', '명목', '1=A, 2=B, 3=C', ''],
    ['연령', '연령(세)', '연속', '', ''], ['점수', '시험 점수', '연속', '', 999],
    ['만족도', '만족도(5점)', '서열', '1=매우 낮음, 5=매우 높음', 99]];
  return { data, vars };
}

/** 시트 mock: 읽기용 탭 + 결과 탭 쓰기 기록 */
function makeEnv(overrides) {
  const sheets = sampleSheets();
  const log = { errors: [], resultCalls: [], props: {}, images: [] };
  const resultSheet = {
    getLastRow: () => 0, getMaxRows: () => 5000, getMaxColumns: () => 30,
    insertRowsAfter() {}, insertColumnsAfter() {}, setColumnWidths() {}, activate() {},
    insertImage: (blob, col, row) => {
      const rec = { blobType: blob.contentType, name: blob.name, col, row };
      log.images.push(rec);
      const img = { setWidth: (w) => { rec.width = w; return img; }, setHeight: (h) => { rec.height = h; return img; } };
      return img;
    },
    getRange: (r, c, nr, nc) => {
      const rec = (name) => (arr) => { if (arr.length !== (nr || 1) || arr[0].length !== (nc || 1)) throw new Error(name + ' 크기 불일치'); log.resultCalls.push(name); return rng; };
      const rng = {
        setNumberFormats: rec('nf'), setValues: rec('v'), setFontWeights: rec('fw'), setFontSizes: rec('fs'), setFontColors: rec('fc'),
        setFontStyles: rec('fst'), setBackgrounds: rec('bg'), setHorizontalAlignments: rec('ha'),
        setBorder() { log.resultCalls.push('border'); return rng; }, activate() { return rng; }
      };
      return rng;
    }
  };
  const tab = rows => ({ getLastRow: () => rows.length, getDataRange: () => ({ getValues: () => JSON.parse(JSON.stringify(rows)) }) });
  const env = {
    console: { log() {}, error: (m) => log.errors.push(String(m)) },
    SpreadsheetApp: {
      BorderStyle: { SOLID: 'SOLID', SOLID_MEDIUM: 'SOLID_MEDIUM' },
      getActiveSpreadsheet: () => ({
        getSheetByName: (n) => (n === '데이터' ? tab(sheets.data) : n === '변수정의' ? tab(sheets.vars) : null),
        insertSheet: () => resultSheet,
        getSheets: () => [{ getName: () => '데이터' }, { getName: () => '변수정의' }]
      })
    },
    Utilities: {
      formatDate: () => '2026-10-21 12:00:00',
      base64Decode: (s) => Buffer.from(s, 'base64'),
      newBlob: (data, contentType, name) => ({ data, contentType, name })
    },
    Buffer,
    Session: { getScriptTimeZone: () => 'Asia/Seoul' },
    PropertiesService: { getUserProperties: () => ({ getProperty: (k) => (k in log.props ? log.props[k] : null), setProperty: (k, v) => { log.props[k] = v; } }) }
  };
  Object.assign(env, overrides || {});
  return { env, log };
}

/** src/*.js 를 한 전역 범위에 불러온다 (Apps Script처럼). reverse: 반대 순서 */
function loadGas(reverse, overrides) {
  const { env, log } = makeEnv(overrides);
  const ctx = vm.createContext(env);
  let files = fs.readdirSync(SRC).filter(f => f.endsWith('.js')).sort();
  if (reverse) files = files.reverse();
  for (const f of files) vm.runInContext(fs.readFileSync(path.join(SRC, f), 'utf8'), ctx, { filename: f });
  return { ctx, log };
}
const call = (ctx, expr) => vm.runInContext(expr, ctx);

// vm 컨텍스트에서는 전역 변수 접근이 느려서(실제 Apps Script와 무관) Studentized range 를 쓰는 사후검정(Tukey, Games-Howell)은
// 반복 실행하는 요청에서 뺀다. 그 전역 연결은 아래 별도 테스트가 Tukey 로 한 번 확인하고, 계산 자체는 anova.test.js 가 확인한다.
const REQUESTS = {
  frequencies: { analysis: 'frequencies', vars: ['성별', '만족도'], options: { barChart: true } },
  descriptives: { analysis: 'descriptives', vars: ['연령', '점수'], options: { skewness: true, kurtosis: true, histogram: true } },
  normality: { analysis: 'normality', vars: ['점수', '연령'], options: {} },
  levene: { analysis: 'levene', dependent: ['점수'], factor: ['집단'], options: {} },
  ttest1: { analysis: 'ttest1', vars: ['점수'], options: { testValue: '80', assumptions: true, conf: '95' } },
  ttestInd: { analysis: 'ttestInd', dependent: ['점수'], factor: ['성별'], options: { assumptions: true, boxplot: true, conf: '95' } },
  ttestPaired: { analysis: 'ttestPaired', first: ['연령'], second: ['점수'], options: { assumptions: true, conf: '95' } },
  anova: { analysis: 'anova', dependent: ['점수'], factor: ['집단'], options: { assumptions: true, levene: true, effect: true, welch: true, bonferroni: true, scheffe: true, lsd: true, sidak: true, dunnett: true, control: '2', boxplot: true, conf: '95' } },
  mannWhitney: { analysis: 'mannWhitney', dependent: ['점수'], factor: ['성별'], options: { effect: true, boxplot: true } },
  wilcoxonPaired: { analysis: 'wilcoxonPaired', first: ['연령'], second: ['점수'], options: { effect: true } },
  wilcoxonOne: { analysis: 'wilcoxonOne', vars: ['점수'], options: { testValue: '80', effect: true } },
  kruskalWallis: { analysis: 'kruskalWallis', dependent: ['점수'], factor: ['집단'], options: { dunn: true, effect: true, boxplot: true } },
  correlation: { analysis: 'correlation', vars: ['연령', '점수', '만족도'], options: { assumptions: true, pearson: true, spearman: true, kendall: true, ci: true, scatter: true, conf: '95' } },
  crosstab: { analysis: 'crosstab', rowvar: ['성별'], colvar: ['집단'], options: { expected: true, rowPct: true, colPct: true, totalPct: true, residual: true, effect: true } },
  reliability: { analysis: 'reliability', vars: ['연령', '점수', '만족도'], options: { itemTotal: true, itemDesc: true, interitem: true, scale: true, ci: true, conf: '95' } },
  moderation: { analysis: 'moderation', dependent: ['점수'], independent: ['연령'], moderator: ['만족도'], covariates: [], options: { center: true, jn: true, plot: true, conf: '95' } },
  distance: { analysis: 'distance', vars: ['연령', '점수', '만족도'], options: { target: 'case', measure: 'mahalanobis', standardize: 'z', outlier: true } },
  cart: { analysis: 'cart', dependent: ['성별'], covariates: ['연령', '점수', '만족도'], factors: ['집단'], options: { maxDepth: '3', prune: 'cv1se', folds: '5', seed: '7', save: false } },
  rmanova: { analysis: 'rmanova', vars: ['연령', '점수', '만족도'], options: { posthoc: true, plot: true } },
  friedman: { analysis: 'friedman', vars: ['연령', '점수', '만족도'], options: { posthoc: true, effect: true, boxplot: true } },
  kmeans: { analysis: 'kmeans', vars: ['연령', '점수', '만족도'], options: { k: '3', standardize: 'z', elbow: true, starts: '3', seed: '7' } },
  psm: { analysis: 'psm', treatment: ['성별'], covariates: ['연령', '점수'], factors: ['집단'], outcome: ['만족도'], options: { method: 'ps', caliper: '0.5', ratio: '1', replace: false, conf: '95' } },
  mediation: { analysis: 'mediation', dependent: ['점수'], independent: ['연령'], mediators: ['만족도'], covariates: [], options: { boot: '500', seed: '7', sobel: true, diagram: true, conf: '95' } },
  logistic: { analysis: 'logistic', dependent: ['성별'], covariates: ['연령', '점수'], factors: ['집단'], options: { hl: true, roc: true, cutoff: '0.5', conf: '95' } },
  anova2: { analysis: 'anova2', dependent: ['점수'], factor1: ['성별'], factor2: ['집단'], options: { ss: 'III', levene: true, plot: true, simple: true, posthoc: 'tukey', emmeans: true, conf: '95' } },
  cochran: { analysis: 'cochran', vars: ['성별', '집단', '만족도'], options: { success: '1', posthoc: true, plot: true } },
  binom: { analysis: 'binom', vars: ['성별', '집단'], options: { success: '1', p0: '0.5', tails: 'two', conf: '95' } },
  km: { analysis: 'km', time: ['연령'], event: ['성별'], group: ['집단'], options: { eventValue: '1', times: '', tests: true, pairs: true, plot: true, table: true, conf: '95' } },
  power: { analysis: 'power', options: { test: 'ttest2', solve: 'n', effect: '0.5', alpha: '0.05', power: '0.8', n: '', k: '3', u: '3', p1: '', p2: '', tails: 'two', ratio: '1', dropout: '10' } },
  ologit: { analysis: 'ologit', dependent: ['집단'], covariates: ['연령', '점수'], factors: ['성별'], options: { parallel: true, classify: true, save: false, conf: '95' } },
  kappa: { analysis: 'kappa', vars: ['성별', '집단'], options: { weights: 'none', conf: '95' } },
  icc: { analysis: 'icc', vars: ['연령', '점수', '만족도'], options: { conf: '95' } },
  efa: { analysis: 'efa', vars: ['연령', '점수', '만족도'], options: { method: 'paf', rotation: 'varimax', nfactors: '1', scree: true, parallel: true, sort: true, save: false, seed: '2026', suppress: '0.3' } },
  hierarchical: { analysis: 'hierarchical', dependent: ['점수'], block1: ['연령', '성별'], block2: ['만족도'], options: { ci: true, collinearity: true, dw: true } },
  regression: { analysis: 'regression', dependent: ['점수'], covariates: ['연령'], factors: ['성별', '집단'], options: { ci: true, collinearity: true, dw: true, conf: '95' } }
};

for (const reverse of [false, true]) {
  const order = reverse ? '역순' : '정순';

  test(`GAS 시뮬레이션(${order} 로드): 모든 분석이 끝까지 실행되고 결과가 JSON 으로 오갈 수 있다`, () => {
    const { ctx, log } = loadGas(reverse);
    ctx.__reqs = REQUESTS;
    for (const key of Object.keys(REQUESTS)) {
      const res = call(ctx, `runAnalysis(__reqs['${key}'])`);
      assert.equal(typeof res.title, 'string', key);
      assert.ok(res.html.includes('<table'), `${key}: 표가 없음`);
      assert.ok(Array.isArray(res.result.blocks) && res.result.blocks.length > 0, key);
      const roundTrip = JSON.parse(JSON.stringify(res.result));              // NaN, Infinity 가 남아 있지 않아야 한다
      assert.deepEqual(roundTrip, JSON.parse(JSON.stringify(roundTrip)));
      assert.ok(!JSON.stringify(res.result).includes('NaN'), key);
    }
    assert.deepEqual(JSON.parse(JSON.stringify(log.errors)), []);
  });

  test(`GAS 시뮬레이션(${order} 로드): 변수 목록, 결과 내보내기, 옵션 저장·복원`, () => {
    const { ctx, log } = loadGas(reverse);
    const info = call(ctx, 'getDatasetInfo()');
    assert.equal(info.n, 60);
    assert.deepEqual(JSON.parse(JSON.stringify(info.vars.map(v => v.name))), ['성별', '집단', '연령', '점수', '만족도']);   // vm 컨텍스트의 배열은 다른 realm 이라 JSON 으로 정규화한다
    ctx.__req = REQUESTS.anova;
    const res = call(ctx, 'runAnalysis(__req)');
    ctx.__res = JSON.parse(JSON.stringify(res.result));                      // 클라이언트 왕복
    const out = call(ctx, 'exportResults([__res])');
    assert.equal(out.sheet, '결과');
    assert.ok(out.lastRow > 20 && log.resultCalls.includes('v') && log.resultCalls.includes('border'));
    ctx.__opts = { anova: { tukey: true, conf: '99' }, bad: 5, 'x y': { a: 1 } };
    assert.equal(call(ctx, 'saveOptions(__opts)'), true);
    assert.deepEqual(JSON.parse(JSON.stringify(call(ctx, 'getSavedOptions()'))), { anova: { tukey: true, conf: '99' } });
  });
}

test('GAS 시뮬레이션: ANOVA 사후검정(Tukey)이 ptukey/qtukey 전역 함수를 거쳐 끝까지 동작한다', () => {
  const { ctx, log } = loadGas(true);
  ctx.__r = { analysis: 'anova', dependent: ['점수'], factor: ['집단'], options: { tukey: true, conf: '95' } };
  const res = call(ctx, 'runAnalysis(__r)');
  assert.ok(res.html.includes('다중비교: Tukey HSD'));
  const tukey = res.result.blocks.find(b => b.caption === '다중비교: Tukey HSD');
  assert.equal(tukey.rows.length, 3);
  assert.ok(tukey.rows.every(r => typeof r.p === 'number' && r.p > 0 && r.p < 1));
  assert.equal(tukey.rows.length === 3 && Math.abs(tukey.rows[0].p - 0.160) < 0.0005, true);   // scipy tukey_hsd 와 같은 값
  assert.deepEqual(JSON.parse(JSON.stringify(log.errors)), []);
});

test('GAS 시뮬레이션: 차트 그림(PNG)이 insertImage 로 들어가고, 그림이 없는 차트는 안내 줄만 남는다', () => {
  const PNG = 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwAEhQGAhKmMIQAAAABJRU5ErkJggg==';
  for (const reverse of [false, true]) {
    const { ctx, log } = loadGas(reverse);
    ctx.__r = { analysis: 'descriptives', vars: ['연령', '점수'], options: { histogram: true } };          // 차트 2개
    const res = call(ctx, 'runAnalysis(__r)');
    const result = JSON.parse(JSON.stringify(res.result));
    assert.equal(result.blocks.filter(b => b.type === 'chart').length, 2);
    ctx.__res = [result];
    ctx.__imgs = [[{ png: PNG, width: 480, height: 300 }, null]];                                         // 둘째는 변환 실패(null)
    const out = call(ctx, 'exportResults(__res, __imgs)');
    assert.equal(out.images, 1);
    assert.equal(log.images.length, 1);
    assert.deepEqual(JSON.parse(JSON.stringify(log.images[0])), { blobType: 'image/png', name: 'chart-1-1.png', col: 1, row: log.images[0].row, width: 480, height: 300 });
    assert.ok(log.images[0].row > 3);
    assert.deepEqual(JSON.parse(JSON.stringify(log.errors)), []);
  }
});

test('GAS 시뮬레이션: 정순과 역순 로드가 같은 HTML 을 만든다', () => {
  const a = loadGas(false), b = loadGas(true);
  for (const key of ['descriptives', 'regression', 'crosstab', 'correlation', 'ttestInd']) {
    a.ctx.__r = REQUESTS[key]; b.ctx.__r = REQUESTS[key];
    assert.equal(call(a.ctx, 'runAnalysis(__r)').html, call(b.ctx, 'runAnalysis(__r)').html, key);
  }
});

test('GAS 시뮬레이션: 전역 함수 이름이 파일 사이에 겹치지 않는다', () => {
  // 각 파일이 최상위에 선언하는 함수 이름을 모아 중복을 찾는다 (겹치면 GAS에서 나중 파일이 앞 파일을 덮어쓴다)
  const seen = {};
  for (const f of fs.readdirSync(SRC).filter(x => x.endsWith('.js'))) {
    const code = fs.readFileSync(path.join(SRC, f), 'utf8');
    for (const m of code.matchAll(/^function\s+([A-Za-z0-9_]+)\s*\(/gm)) {
      (seen[m[1]] = seen[m[1]] || []).push(f);
    }
    for (const m of code.matchAll(/^var\s+([A-Za-z0-9_]+)\s*=/gm)) {
      (seen[m[1]] = seen[m[1]] || []).push(f + '(var)');
    }
  }
  const dup = Object.entries(seen).filter(([, files]) => files.length > 1);
  assert.deepEqual(dup, [], '중복 선언: ' + JSON.stringify(dup));
});

test('오류 안내: 사용자 입력 오류는 그대로, 프로그램 오류는 안내 문구로 바뀌고 로그에 남는다', () => {
  const { ctx, log } = loadGas(false);
  assert.throws(() => call(ctx, "runAnalysis({ analysis: 'nope' })"), /알 수 없는 분석입니다/);
  assert.throws(() => call(ctx, "runAnalysis({ analysis: 'frequencies', vars: [] })"), /변수를 선택/);
  assert.deepEqual(log.errors, []);                               // 입력 오류는 로그에 남기지 않는다
  // 프로그램 오류: 시트가 이상한 값을 돌려주는 경우 (TypeError)
  const broken = loadGas(false, {
    SpreadsheetApp: { getActiveSpreadsheet: () => ({ getSheetByName: () => ({ getDataRange: () => ({ getValues: () => null }) }) }) }
  });
  let err;
  try { call(broken.ctx, "runAnalysis({ analysis: 'frequencies', vars: ['성별'] })"); } catch (e) { err = e; }
  assert.ok(err && err.message.startsWith('예기치 않은 오류가 발생했습니다 ('), err && err.message);
  assert.ok(err.message.includes('TypeError'));
  assert.equal(broken.log.errors.length, 1);                      // 원인이 실행 로그에 남는다
  assert.ok(broken.log.errors[0].includes('runAnalysis'));
});

/* ---------- Safety.js 순수 함수 ---------- */

test('sanitizeSavedOptions_: 허용하는 모양만 남긴다', () => {
  const S = Safety.sanitizeSavedOptions_;
  assert.deepEqual(S({ anova: { tukey: true, conf: '99', n: 5 } }), { anova: { tukey: true, conf: '99', n: 5 } });
  assert.deepEqual(S('{"a":{"b":true}}'), { a: { b: true } });                       // JSON 문자열도 받는다
  assert.deepEqual(S('깨진 JSON'), {});
  assert.deepEqual(S(null), {});
  assert.deepEqual(S([1, 2]), {});
  assert.deepEqual(S({ a: [1], b: 'x', 'bad key': { x: true }, c: { 'bad key': true, ok: false } }), { c: { ok: false } });
  assert.deepEqual(S({ a: { x: 'ㅇ'.repeat(30), y: NaN, z: Infinity, w: { n: 1 }, v: null, u: undefined, t: true } }), { a: { t: true } });
  const many = {};
  for (let i = 0; i < 50; i++) many['k' + i] = { o: true };
  assert.equal(Object.keys(S(many)).length, 30);                                    // 분석 개수 상한
});

test('asUserError_: 입력 오류(Error)는 그대로, 내장 오류는 변환', () => {
  const logs = [];
  const logger = m => logs.push(m);
  const input = new Error('변수를 선택하세요.');
  assert.equal(Safety.asUserError_(input, 'x', logger), input);
  assert.equal(logs.length, 0);
  for (const E of [TypeError, ReferenceError, RangeError, SyntaxError]) {
    const out = Safety.asUserError_(new E('원인 문구'), '어딘가', logger);
    assert.ok(out instanceof Error && !(out instanceof E));
    assert.ok(out.message.includes(E.name) && out.message.includes('원인 문구') && out.message.includes('계속되면'));
  }
  assert.equal(logs.length, 4);
  assert.ok(logs[0].includes('[어딘가]'));
});

test('largeDataWarning_: 2만 건을 넘으면 경고', () => {
  assert.equal(Safety.largeDataWarning_(20000), null);
  assert.ok(Safety.largeDataWarning_(20001).includes('20001건'));
});

/* ---------- 배포 준비: 빈 데이터 안내, 템플릿 점검 ---------- */

function emptyTemplateApp(extraTabs) {
  const tab = rows => ({ getLastRow: () => rows.length, getDataRange: () => ({ getValues: () => JSON.parse(JSON.stringify(rows)) }) });
  const sheets = { '데이터': tab([['']]), '변수정의': tab([['변수명', '라벨', '척도', '값 레이블', '결측값 코드'], ['', '', '', '', '']]) };
  return {
    getActiveSpreadsheet: () => ({
      getSheetByName: (n) => sheets[n] || null,
      getSheets: () => Object.keys(sheets).concat(extraTabs || []).map(n => ({ getName: () => n }))
    })
  };
}

test('배포 준비: 데이터 탭이 비어 있으면 어떻게 하면 되는지 알려 준다 (변수명이 아닌 안내 문구가 변수로 읽히지 않는다)', () => {
  const { ctx } = loadGas(false, { SpreadsheetApp: emptyTemplateApp() });
  assert.throws(() => call(ctx, 'getDatasetInfo()'), /1행에 변수명을 입력하고 2행부터 데이터를 넣으세요/);
  assert.throws(() => call(ctx, "runAnalysis({ analysis: 'frequencies', vars: ['x'] })"), /1행에 변수명을 입력하고/);
});

test('배포 준비: 템플릿 점검이 시트를 읽어 요약하고, 예제가 남은 시트와 비운 시트를 구분한다', () => {
  const dev = loadGas(false);                                        // 예제 데이터가 들어 있는 시트
  const s1 = JSON.parse(JSON.stringify(call(dev.ctx, 'collectTemplateSummary_()')));
  assert.equal(s1.dataRows, 60);
  assert.equal(s1.dataHeaderCells, 5);
  assert.equal(s1.varRows, 5);
  assert.equal(s1.hasResult, false);
  dev.ctx.__s = s1;
  const r1 = JSON.parse(JSON.stringify(call(dev.ctx, 'evaluateTemplateReadiness_(__s)')));
  assert.equal(r1.ok, false);
  const empty = loadGas(false, { SpreadsheetApp: emptyTemplateApp() });
  const s2 = JSON.parse(JSON.stringify(call(empty.ctx, 'collectTemplateSummary_()')));
  assert.deepEqual([s2.dataRows, s2.dataHeaderCells, s2.varRows, s2.hasData, s2.hasVars, s2.otherTabs], [0, 0, 0, true, true, []]);
  empty.ctx.__s = s2;
  assert.deepEqual(JSON.parse(JSON.stringify(call(empty.ctx, 'evaluateTemplateReadiness_(__s)'))), { ok: true, problems: [] });
  const withExtra = loadGas(false, { SpreadsheetApp: emptyTemplateApp(['시트1']) });
  assert.deepEqual(JSON.parse(JSON.stringify(call(withExtra.ctx, 'collectTemplateSummary_()').otherTabs)), ['시트1']);
});

/* ---------- 행 수 제한, 입력 시 경고 ---------- */

function bigSheetApp(rowCount, reads) {
  const header = [['x', 'y']];
  const rows = header.concat(Array.from({ length: rowCount }, (_, i) => [i, i % 7]));
  const data = { getLastRow: () => rows.length, getDataRange: () => ({ getValues: () => { reads.count++; return rows; } }) };
  const vars = { getLastRow: () => 1, getDataRange: () => ({ getValues: () => [['변수명', '라벨', '척도', '값 레이블', '결측값 코드']] }) };
  return { getActiveSpreadsheet: () => ({ getSheetByName: (n) => (n === '데이터' ? data : n === '변수정의' ? vars : null) }) };
}

test('행 수 제한: 한도 이하는 분석하고, 한도를 넘으면 거부하며 줄이는 방법을 알려 준다', () => {
  const reads = { count: 0 };
  const okEnv = loadGas(false, { SpreadsheetApp: bigSheetApp(50000, reads) });          // 정확히 한도 (50000건)
  assert.equal(call(okEnv.ctx, 'getDatasetInfo()').n, 50000);
  const reads2 = { count: 0 };
  const over = loadGas(false, { SpreadsheetApp: bigSheetApp(50001, reads2) });
  assert.throws(() => call(over.ctx, 'getDatasetInfo()'), /사례가 50001건 있습니다.*최대 50000건까지.*줄여 주세요/);
  assert.throws(() => call(over.ctx, "runAnalysis({ analysis: 'descriptives', vars: ['x'] })"), /최대 50000건/);
  // 한도의 4배를 넘으면 시트를 읽기 전에 거부한다 (읽기만으로 오래 걸리지 않도록)
  const reads3 = { count: 0 };
  const huge = loadGas(false, { SpreadsheetApp: bigSheetApp(200001, reads3) });
  assert.throws(() => call(huge.ctx, 'getDatasetInfo()'), /최대 50000건/);
  assert.equal(reads3.count, 0);
});

test('입력 시 경고: 데이터 탭에서 한도를 넘는 행에 입력하면 토스트, 그 밖에는 조용', () => {
  const { ctx } = loadGas(false);
  const toasts = [];
  ctx.__ev = (sheetName, lastRow) => ({
    range: { getSheet: () => ({ getName: () => sheetName }), getLastRow: () => lastRow },
    source: { toast: (msg, title, sec) => toasts.push({ msg, title, sec }) }
  });
  call(ctx, "onEdit(__ev('데이터', 100))");
  call(ctx, "onEdit(__ev('데이터', 50001))");                 // 머리글 1행 + 사례 50000건: 한도 안
  call(ctx, "onEdit(__ev('결과', 90000))");                   // 다른 탭
  assert.equal(toasts.length, 0);
  call(ctx, "onEdit(__ev('데이터', 50002))");                 // 사례 50001건째
  call(ctx, "onEdit(__ev('데이터', 120000))");                // 한꺼번에 붙여 넣은 경우
  assert.equal(toasts.length, 2);
  assert.ok(toasts[0].msg.includes('50002행까지') && toasts[0].msg.includes('최대 50000건') && toasts[0].msg.includes('분석이 거부'));
  assert.ok(toasts[0].title.includes('50000건') && toasts[0].sec > 0);
  assert.doesNotThrow(() => call(ctx, 'onEdit(undefined)'));       // 이벤트가 없어도 편집에 영향이 없다
  assert.doesNotThrow(() => call(ctx, "onEdit({ range: null })"));
});

test('결과 탭 지우기: 값·서식·그림을 지우고 다른 탭은 건드리지 않으며, 확인창에서 아니오를 누르면 아무것도 하지 않는다', () => {
  function resultApp(state) {
    const sheet = {
      getLastRow: () => state.rows,
      getImages: () => state.images.map(() => ({ remove: () => { state.removed++; } })),
      clear: () => { state.cleared++; state.rows = 0; }
    };
    return {
      BorderStyle: {},
      getUi: () => ({
        ButtonSet: { OK: 'OK', YES_NO: 'YES_NO' }, Button: { YES: 'YES', NO: 'NO' },
        alert: (title, msg, buttons) => { state.alerts.push({ title, msg, buttons }); return buttons === 'YES_NO' ? state.answer : 'OK'; }
      }),
      getActiveSpreadsheet: () => ({ getSheetByName: (n) => { state.asked.push(n); return n === '결과' && state.hasResult ? sheet : null; } })
    };
  }
  const base = () => ({ rows: 40, images: [1, 2], removed: 0, cleared: 0, alerts: [], asked: [], answer: 'YES', hasResult: true });

  let s = base();
  let env = loadGas(false, { SpreadsheetApp: resultApp(s) });
  call(env.ctx, 'clearResultSheet()');
  assert.equal(s.removed, 2);
  assert.equal(s.cleared, 1);
  assert.deepEqual(JSON.parse(JSON.stringify(s.asked)), ['결과']);                   // 결과 탭만 찾는다
  assert.ok(s.alerts[0].msg.includes('40행, 그림 2개') && s.alerts[0].msg.includes('다른 탭은 그대로'));
  assert.ok(s.alerts[1].msg.includes('비웠습니다'));

  s = base(); s.answer = 'NO';
  env = loadGas(false, { SpreadsheetApp: resultApp(s) });
  call(env.ctx, 'clearResultSheet()');
  assert.deepEqual([s.removed, s.cleared, s.alerts.length], [0, 0, 1]);                // 확인창만 뜨고 지우지 않는다

  s = base(); s.hasResult = false;
  env = loadGas(false, { SpreadsheetApp: resultApp(s) });
  call(env.ctx, 'clearResultSheet()');
  assert.ok(s.alerts[0].msg.includes('탭이 없습니다'));
  assert.equal(s.cleared, 0);

  s = base(); s.rows = 0; s.images = [];
  env = loadGas(false, { SpreadsheetApp: resultApp(s) });
  call(env.ctx, 'clearResultSheet()');
  assert.ok(s.alerts[0].msg.includes('이미 비어 있습니다'));
  assert.equal(s.cleared, 0);
});

test('결과 탭 지우기: 안내에 지울 탭을 모두 적고, 예를 누르면 거리·유사도 행렬 탭(_2 포함)도 삭제하며 사용자의 다른 탭은 지우지 않는다', () => {
  function app(state) {
    const mkTab = (name, rows, cols) => ({ getName: () => name, getLastRow: () => rows, getLastColumn: () => cols });
    const result = { getName: () => '결과', getLastRow: () => state.rows, getImages: () => [], clear: () => { state.cleared++; state.rows = 0; } };
    const sheets = [mkTab('데이터', 61, 5), result, mkTab('거리행렬', 41, 41), mkTab('유사도행렬_2', 11, 11), mkTab('거리행렬표_내메모', 3, 3), mkTab('내 거리행렬', 3, 3)];
    state.sheets = sheets;
    return {
      BorderStyle: {},
      getUi: () => ({
        ButtonSet: { OK: 'OK', YES_NO: 'YES_NO' }, Button: { YES: 'YES', NO: 'NO' },
        alert: (title, msg, buttons) => { state.alerts.push({ title, msg, buttons }); return buttons === 'YES_NO' ? state.answer : 'OK'; }
      }),
      getActiveSpreadsheet: () => ({
        getSheetByName: (n) => sheets.find((s) => s.getName() === n) || null,
        getSheets: () => sheets.slice(),
        deleteSheet: (s) => { state.deleted.push(s.getName()); sheets.splice(sheets.indexOf(s), 1); }
      })
    };
  }
  const base = () => ({ rows: 40, cleared: 0, alerts: [], deleted: [], answer: 'YES' });
  let s = base();
  let env = loadGas(false, { SpreadsheetApp: app(s) });
  call(env.ctx, 'clearResultSheet()');
  const ask = s.alerts[0].msg;
  assert.ok(ask.includes('"결과" 탭의 내용 전체 (40행') && ask.includes('"거리행렬" 탭 (거리·유사도 행렬, 41행 × 41열) — 탭 자체를 삭제') && ask.includes('"유사도행렬_2" 탭'));
  assert.ok(ask.includes('되돌릴 수 없습니다') && ask.includes('다른 탭은 그대로'));
  assert.ok(!ask.includes('내메모') && !ask.includes('내 거리행렬'));                              // 이름이 비슷한 사용자 탭은 목록에 없다
  assert.deepEqual(JSON.parse(JSON.stringify(s.deleted)), ['거리행렬', '유사도행렬_2']);
  assert.equal(s.cleared, 1);
  assert.ok(s.alerts[1].msg.includes('행렬 탭 2개를 삭제했습니다: 거리행렬, 유사도행렬_2'));
  assert.deepEqual(s.sheets.map((x) => x.getName()), ['데이터', '결과', '거리행렬표_내메모', '내 거리행렬']);

  s = base(); s.answer = 'NO';                                                                      // 아니오: 아무것도 지우지 않는다
  env = loadGas(false, { SpreadsheetApp: app(s) });
  call(env.ctx, 'clearResultSheet()');
  assert.deepEqual([s.cleared, s.deleted.length, s.alerts.length], [0, 0, 1]);

  s = base(); s.rows = 0;                                                                          // 결과 탭은 비었지만 행렬 탭이 있으면 그것만 지운다
  env = loadGas(false, { SpreadsheetApp: app(s) });
  call(env.ctx, 'clearResultSheet()');
  assert.ok(!s.alerts[0].msg.includes('"결과" 탭의 내용'));
  assert.deepEqual([s.cleared, JSON.parse(JSON.stringify(s.deleted))], [0, ['거리행렬', '유사도행렬_2']]);
});

test('메뉴 "Stat4SW 통계분석": 자주 쓰는 항목은 위에, 고급 분석 켜기/끄기와 관리 기능은 하위 메뉴 "설정·관리"에 있고 모든 항목이 실제 함수를 가리킨다', () => {
  const tree = [];
  const mkMenu = (name) => {
    const m = { name, items: [], addItem(label, fn) { m.items.push({ label, fn }); return m; }, addSeparator() { m.items.push({ sep: true }); return m; }, addSubMenu(sub) { m.items.push({ sub }); return m; }, addToUi() { tree.push(m); return m; } };
    return m;
  };
  const env = loadGas(false, { SpreadsheetApp: { getUi: () => ({ createMenu: mkMenu }), getActiveSpreadsheet: () => ({}) } });
  call(env.ctx, 'onOpen()');
  const menu = tree[0];
  assert.equal(menu.name, 'Stat4SW 통계분석');
  const top = menu.items.filter((i) => i.label).map((i) => i.label);
  assert.deepEqual(top, ['분석 창 열기', '변수 설정 편집', '변수 계산 (역문항, 합계·평균, z점수)', '결과 탭 지우기 (거리·유사도 행렬 탭 포함)', '새 작업 시작 (데이터 모두 지우기)', '예제 데이터 채우기', '분석 선택 도우미 (넓게 보기)', '사용 안내', '정보·저작권']);
  const sub = menu.items.find((i) => i.sub).sub;
  assert.equal(sub.name, '설정·관리');
  assert.deepEqual(sub.items.map((i) => i.label), ['고급 분석 켜기/끄기', '시험 데이터 (기법별 점검용)']);
  const all = menu.items.concat(sub.items).filter((i) => i.fn);
  all.forEach((i) => assert.equal(call(env.ctx, 'typeof ' + i.fn), 'function', i.label + ' → ' + i.fn));
  assert.equal(menu.items[0].label, '분석 창 열기');
});

test('rowLimitMessage_: 한도 이하는 null, 초과는 건수·한도·줄이는 방법을 담은 안내', () => {
  assert.equal(Safety.rowLimitMessage_(50000, 50000, '데이터'), null);
  assert.equal(Safety.rowLimitMessage_(0, 50000, '데이터'), null);
  const m = Safety.rowLimitMessage_(50001, 50000, '데이터');
  assert.ok(m.includes('"데이터" 탭에 사례가 50001건') && m.includes('최대 50000건') && m.includes('나누거나 표본을 추출'));
});

test('GAS 시뮬레이션: IPA(변수 6개가 필요해 별도 시트로)가 전역 함수 연결과 그림 만들기까지 동작한다', () => {
  const data = [['중요1', '중요2', '중요3', '수행1', '수행2', '수행3']];
  for (let i = 0; i < 40; i++) data.push([1 + (i % 5), 1 + ((i * 3) % 5), 1 + ((i * 7) % 5), 1 + ((i + 2) % 5), 1 + ((i * 2 + 1) % 5), 1 + ((i * 5 + 3) % 5)]);
  const vars = [['변수명', '라벨', '척도', '값 레이블', '결측값 코드']].concat(data[0].map(n => [n, n, '서열', '', '']));
  const tab = rows => ({ getLastRow: () => rows.length, getDataRange: () => ({ getValues: () => JSON.parse(JSON.stringify(rows)) }) });
  const SpreadsheetApp = {
    BorderStyle: { SOLID: 'SOLID', SOLID_MEDIUM: 'SOLID_MEDIUM' },
    getActiveSpreadsheet: () => ({ getSheetByName: n => (n === '데이터' ? tab(data) : n === '변수정의' ? tab(vars) : null), getSheets: () => [] })
  };
  const { ctx, log } = loadGas(false, { SpreadsheetApp });
  ctx.__r = { analysis: 'ipa', importance: ['중요1', '중요2', '중요3'], performance: ['수행1', '수행2', '수행3'], options: { baseline: 'median', impShift: '-0.1', perfShift: '0.2', gap: true, plot: true } };
  const res = call(ctx, 'runAnalysis(__r)');
  assert.ok(res.html.includes('<svg') && res.html.includes('집중 개선'));
  const chart = res.result.blocks.find(b => b.type === 'chart');
  assert.equal(chart.kind, 'ipa');
  assert.ok(!JSON.stringify(res.result).includes('NaN'));
  assert.deepEqual(JSON.parse(JSON.stringify(log.errors)), []);
});

test('GAS 시뮬레이션: 고급 분석 켜기/끄기 — 기본은 꺼짐, 켜면 분석 창·도우미에 알리고, 관리 기능(시험 데이터, 템플릿 시트)은 켠 뒤에만 열린다', () => {
  const calls = { alerts: [], dialogs: [] };
  const base = sampleSheets();
  const tab = rows => ({ getLastRow: () => rows.length, getDataRange: () => ({ getValues: () => JSON.parse(JSON.stringify(rows)) }) });
  const ui = { alert: (t, m) => { calls.alerts.push({ t, m }); return 'OK'; }, ButtonSet: { OK: 'OK', YES_NO: 'YES_NO' }, Button: { YES: 'YES' }, showModalDialog: (h, t) => { calls.dialogs.push(t); }, showModelessDialog: (h, t) => { calls.dialogs.push(t); } };
  const SpreadsheetApp = {
    BorderStyle: { SOLID: 'SOLID', SOLID_MEDIUM: 'SOLID_MEDIUM' },
    getUi: () => ui,
    getActiveSpreadsheet: () => ({ getSheetByName: n => (n === '데이터' ? tab(base.data) : n === '변수정의' ? tab(base.vars) : null), getSheets: () => [] })
  };
  const out = () => ({ setWidth() { return out(); }, setHeight() { return out(); } });
  const HtmlService = { createHtmlOutputFromFile: () => out(), createTemplateFromFile: () => ({ evaluate: () => out() }) };
  const { ctx } = loadGas(false, { SpreadsheetApp, HtmlService });
  const info0 = call(ctx, 'getDatasetInfo()');
  assert.equal(info0.advanced, false);
  assert.equal(info0.advancedKeys.length, 19);
  assert.equal(call(ctx, 'getAdvisor()').advanced, false);
  call(ctx, 'menuOpenTestData()');
  assert.deepEqual(calls.dialogs, []);                                  // 꺼져 있으면 열리지 않는다
  assert.ok(calls.alerts[0].m.includes('고급 분석을 켠 뒤') && calls.alerts[0].m.includes('설정·관리 > 고급 분석 켜기/끄기'));
  call(ctx, 'toggleAdvancedMode()');                                    // 켠다
  assert.ok(calls.alerts[1].t.includes('켜짐') && calls.alerts[1].m.includes('19개'));
  assert.equal(call(ctx, 'getDatasetInfo()').advanced, true);
  assert.equal(call(ctx, 'getAdvisor()').advanced, true);
  call(ctx, 'menuOpenTestData()');
  assert.deepEqual(calls.dialogs, ['시험 데이터']);                      // 켠 뒤에는 열린다
  call(ctx, 'toggleAdvancedMode()');                                    // 끈다
  assert.ok(calls.alerts[2].t.includes('꺼짐') && calls.alerts[2].m.includes('18개'));
  assert.equal(call(ctx, 'getDatasetInfo()').advanced, false);
  call(ctx, 'showAbout()');
  assert.equal(calls.dialogs[calls.dialogs.length - 1], 'Stat4SW 정보·저작권');
});

test('GAS 시뮬레이션: 개인 사본 안내는 사용자마다 처음 분석 창을 열 때 한 번만 보이고, 안내를 못 보여 줘도 분석 창은 열린다', () => {
  const calls = { alerts: [], dialogs: [] };
  const ui = { alert: (t, m) => { calls.alerts.push({ t, m }); return 'OK'; }, ButtonSet: { OK: 'OK' }, showModelessDialog: (h, t) => { calls.dialogs.push(t); } };
  const SpreadsheetApp = { getUi: () => ui, getActiveSpreadsheet: () => ({ getSheetByName: () => null, getSheets: () => [] }) };
  const out = () => ({ setWidth() { return out(); }, setHeight() { return out(); } });
  const HtmlService = { createTemplateFromFile: () => ({ evaluate: () => out() }) };
  const { ctx, log } = loadGas(false, { SpreadsheetApp, HtmlService });
  call(ctx, 'openAnalysis()');
  assert.equal(calls.alerts.length, 1);
  assert.equal(calls.alerts[0].t, '개인 사본 안내');
  assert.ok(calls.alerts[0].m.includes('개인 사본') && calls.alerts[0].m.includes('공유하면 데이터도 함께 보입니다') && calls.alerts[0].m.includes('익명화'));
  assert.equal(log.props['stat4sw.privacyNoticed'], '1');
  assert.deepEqual(calls.dialogs, ['Stat4SW']);
  call(ctx, 'openAnalysis()');                                          // 두 번째부터는 안내 없이 열린다
  assert.equal(calls.alerts.length, 1);
  assert.equal(calls.dialogs.length, 2);
  const bad = loadGas(false, { SpreadsheetApp: { getUi: () => ({ showModelessDialog: (h, t) => { calls.dialogs.push('bad:' + t); } }), getActiveSpreadsheet: () => ({}) }, HtmlService });
  call(bad.ctx, 'openAnalysis()');                                      // alert 가 없는 환경에서도 분석 창은 열린다
  assert.ok(calls.dialogs.includes('bad:Stat4SW'));
});

test('GAS 시뮬레이션: 데이터 탭이 없어 getDatasetInfo 가 실패해도 고급 모드 정보는 따로 받을 수 있고 분석 창 HTML 에 끼워진다', () => {
  const SpreadsheetApp = {
    BorderStyle: { SOLID: 'SOLID', SOLID_MEDIUM: 'SOLID_MEDIUM' },
    getActiveSpreadsheet: () => ({ getSheetByName: () => null, getSheets: () => [] })
  };
  const { ctx, log } = loadGas(false, { SpreadsheetApp });
  assert.throws(() => call(ctx, 'getDatasetInfo()'), /탭이 없습니다/);
  const m = JSON.parse(JSON.stringify(call(ctx, 'getModeInfo()')));
  assert.equal(m.advanced, false);
  assert.equal(m.advancedKeys.length, 19);
  const boot = JSON.parse(call(ctx, 'modeBootJson_()'));
  assert.deepEqual(boot, m);
  log.props['stat4sw.advanced'] = '1';                                  // 사용자 속성에 켜짐이 저장돼 있으면
  assert.equal(JSON.parse(call(ctx, 'modeBootJson_()')).advanced, true);
});

test('GAS 시뮬레이션: 진단 정보(getDiagnostics)에 버전·모드·탭 크기가 있고 데이터 값, 변수 이름, 시트 제목은 없다', () => {
  const { ctx } = loadGas(false);
  const text = call(ctx, 'getDiagnostics()');
  assert.ok(text.startsWith('Stat4SW 진단 정보'));
  assert.ok(text.includes('버전: ' + call(ctx, 'APP_VERSION')) && text.includes('고급 분석: 꺼짐'));
  assert.ok(text.includes('데이터 탭: 있음 (사례 60행') && text.includes('변수정의 탭: 있음 (정의 5개)') && text.includes('결과 탭: 없음'));
  ['성별', '집단', '연령', '점수', '만족도'].forEach(name => assert.ok(!text.includes(name), '변수 이름이 들어 있음: ' + name));
  const bad = loadGas(false, { SpreadsheetApp: { getActiveSpreadsheet: () => { throw new Error('boom'); } } });
  assert.ok(call(bad.ctx, 'getDiagnostics()').includes('시트 정보를 읽지 못했습니다.'));   // 시트를 못 읽어도 진단 정보는 만들어진다
});
