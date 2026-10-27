const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const A = require('../src/Advisor.js');
const S = require('../src/Sample.js');

const SRC = path.join(__dirname, '..', 'src');
const html = fs.readFileSync(path.join(SRC, 'client.html'), 'utf8');
const start = html.indexOf('  var ANALYSES = {');
const ANALYSIS_KEYS = [...html.slice(start).matchAll(/^    ([A-Za-z0-9]+): \{\n      label:/gm)].map(m => m[1]);
// 역할 키 (분석 열기 때 변수를 넣을 칸)
const ROLE_KEYS = {};
ANALYSIS_KEYS.forEach((k, i) => {
  const from = html.indexOf('\n    ' + k + ': {\n'), to = i + 1 < ANALYSIS_KEYS.length ? html.indexOf('\n    ' + ANALYSIS_KEYS[i + 1] + ': {\n') : html.length;
  const seg = html.slice(from, to), rs = seg.indexOf('roles: ['), re = seg.indexOf('options: [');
  ROLE_KEYS[k] = [...seg.slice(rs, re).matchAll(/key: '([A-Za-z0-9]+)'/g)].map(m => m[1]);
});

function walk(node, trail, out) {
  if (node.q) {
    assert.ok(node.choices.length >= 2, '질문 "' + node.q + '"의 선택지가 2개 미만');
    assert.equal(new Set(node.choices.map(c => c.label)).size, node.choices.length, '선택지 이름 중복: ' + node.q);
    node.choices.forEach((c, i) => walk(c.next, trail.concat([node.q + ' → ' + c.label]), out));
  } else {
    assert.ok(node.rec, '잎에 추천이 없음: ' + trail.join(' / '));
    out.push({ node: node, trail: trail });
  }
  return out;
}
const leaves = walk(A.ADVISOR_TREE_, [], []);

test('질문 나무: 모든 잎의 추천·대안은 실제 분석 키이고, 미지원 분석은 쓸 수 있는 대체 분석을 가리킨다', () => {
  assert.equal(ANALYSIS_KEYS.length >= 19, true);
  leaves.forEach(({ node, trail }) => {
    const where = trail.join(' / ');
    if (node.rec.unsupported) {
      assert.ok(node.rec.name && node.rec.why && ANALYSIS_KEYS.includes(node.rec.fallback), '미지원 추천의 대체 분석: ' + where);
    } else {
      assert.ok(ANALYSIS_KEYS.includes(node.rec.key), '추천 키: ' + node.rec.key + ' (' + where + ')');
      assert.ok(node.rec.why.length >= 10, '이유 문구: ' + where);
    }
    node.alt.forEach(a => assert.ok(ANALYSIS_KEYS.includes(a.key) && a.why, '대안 키: ' + a.key + ' (' + where + ')'));
    assert.ok(trail.length <= 4, '질문이 4단계를 넘음: ' + where);
  });
});

test('질문 나무: 이 도구의 모든 분석이 추천 또는 대안으로 한 번 이상 나온다 (새 분석을 만들면 규칙을 추가해야 한다)', () => {
  const used = new Set();
  leaves.forEach(({ node }) => { if (node.rec.key) used.add(node.rec.key); node.alt.forEach(a => used.add(a.key)); });
  const missing = ANALYSIS_KEYS.filter(k => !used.has(k));
  assert.deepEqual(missing, [], '선택 도우미 규칙에 없는 분석: ' + missing.join(', ') + ' → src/Advisor.js 의 advBuildTree_ 에 추가하세요');
});

test('질문 나무: 같은 질문 경로가 같은 추천으로 끝나고, 비모수 짝이 있는 모수 검정에는 parametric 표시와 대안이 있다', () => {
  const para = leaves.filter(l => l.node.rec.parametric);
  assert.ok(para.length >= 6);
  para.forEach(l => assert.ok(l.node.alt.length >= 1 || l.node.notes.some(n => n.includes('스피어만')), '비모수 대안이 없음: ' + l.trail.join(' / ')));
  const keys = para.map(l => l.node.rec.key);
  ['ttestInd', 'anova', 'ttestPaired', 'ttest1', 'correlation'].forEach(k => assert.ok(keys.includes(k), k));
  // 서열 점수는 비모수가 기본 추천
  const ord = leaves.filter(l => l.trail.some(t => t.includes('서열·리커트 (5점')));
  assert.ok(ord.some(l => l.node.rec.key === 'mannWhitney') && ord.some(l => l.node.rec.key === 'kruskalWallis') && ord.some(l => l.node.rec.key === 'wilcoxonPaired'));
  // 미지원: McNemar, 서열 로지스틱
  const unsupported = leaves.filter(l => l.node.rec.unsupported).map(l => l.node.rec.name);
  assert.ok(!unsupported.some(n => n.includes('이항')));                                       // 이항 검정은 이제 지원한다
  assert.ok(!unsupported.some(n => n.includes('순서형 로지스틱')));                          // 서열 로지스틱은 이제 지원한다
  assert.ok(!unsupported.some(n => n.includes('McNemar')));                                  // McNemar 는 이제 지원한다
  const mc = leaves.find(l => l.node.rec.key === 'mcnemar');
  assert.ok(mc && mc.trail.some(t => t.includes('사전·사후')) && mc.node.alt[0].key === 'crosstab');
});

test('getAdvisor: JSON 으로 오갈 수 있고(순환 없음) 표본 크기 안내를 담는다', () => {
  const r = A.getAdvisor();
  const back = JSON.parse(JSON.stringify(r));
  assert.deepEqual(back.tree.q, A.ADVISOR_TREE_.q);
  assert.equal(back.smallN, 30);
  assert.ok(back.smallNote.includes('가정 점검'));
  assert.equal(back.tree.choices.length, 14);
});

/* ---------- 변수 기반 추천 ---------- */

const sc = (name, scale, groups) => ({ name, scale, groups: groups || 30 });
test('adviseFromItems_: 변수 1개, 2개의 조합별 추천과 역할 배치', () => {
  const f = A.adviseFromItems_;
  assert.equal(f([]).none, true);
  assert.equal(f([sc('성별', 'nominal', 2)]).key, 'frequencies');
  assert.equal(f([sc('점수', 'scale')]).key, 'descriptives');
  const t = f([sc('점수', 'scale'), sc('성별', 'nominal', 2)]);
  assert.equal(t.key, 'ttestInd');
  assert.deepEqual(t.roles, { dependent: ['점수'], factor: ['성별'] });
  assert.equal(t.alt[0].key, 'mannWhitney');
  assert.equal(f([sc('성별', 'nominal', 2), sc('점수', 'scale')]).key, 'ttestInd');                 // 순서와 무관
  assert.equal(f([sc('점수', 'scale'), sc('집단', 'nominal', 3)]).key, 'anova');
  assert.equal(f([sc('만족도', 'ordinal', 5), sc('성별', 'nominal', 2)]).key, 'mannWhitney');
  assert.equal(f([sc('만족도', 'ordinal', 5), sc('집단', 'nominal', 3)]).key, 'kruskalWallis');
  const ct = f([sc('성별', 'nominal', 2), sc('집단', 'nominal', 3)]);
  assert.deepEqual([ct.key, ct.roles], ['crosstab', { rowvar: ['성별'], colvar: ['집단'] }]);
  assert.deepEqual(ct.alt, []);                                                                // 집단이 3개라 McNemar 대안이 없다
  const ct2 = f([sc('사전취업', 'nominal', 2), sc('사후취업', 'nominal', 2)]);
  assert.equal(ct2.alt[0].key, 'mcnemar');
  assert.deepEqual(ct2.alt[0].roles, { first: ['사전취업'], second: ['사후취업'] });
  const cor = f([sc('사전', 'scale'), sc('사후', 'scale')]);
  assert.equal(cor.key, 'correlation');
  assert.deepEqual(cor.alt[0].roles, { first: ['사전'], second: ['사후'] });                              // 사전·사후일 수 있다
  assert.ok(f([sc('a', 'ordinal', 5), sc('b', 'scale')]).why.includes('스피어만'));
  assert.equal(f([sc('점수', 'scale'), sc('x', 'scale'), sc('y', 'scale')]).alt[0].key, 'reliability');
  // 추천할 수 없는 경우는 안내만
  assert.ok(f([sc('점수', 'scale'), sc('ID', 'nominal', 40)]).message.includes('집단 비교로 보기 어렵습니다'));
  assert.ok(f([sc('점수', 'scale'), sc('성별', 'nominal', 2), sc('집단', 'nominal', 3)]).none);
  assert.ok(f([sc('점수', 'scale'), sc('성별', 'nominal', 1)]).message);                                  // 값이 하나뿐인 집단변수
});

test('adviseFromItems_: 추천·대안의 키와 역할 키가 모두 실제 분석 정의에 있다', () => {
  const combos = [
    [sc('a', 'scale')], [sc('a', 'nominal', 2)], [sc('a', 'scale'), sc('b', 'nominal', 2)], [sc('a', 'ordinal', 5), sc('b', 'nominal', 2)],
    [sc('a', 'scale'), sc('b', 'nominal', 4)], [sc('a', 'ordinal', 5), sc('b', 'nominal', 4)], [sc('a', 'nominal', 2), sc('b', 'nominal', 3)],
    [sc('a', 'scale'), sc('b', 'scale')], [sc('a', 'scale'), sc('b', 'scale'), sc('c', 'scale')]
  ];
  combos.forEach(items => {
    const r = A.adviseFromItems_(items);
    assert.ok(!r.none, JSON.stringify(items));
    [r].concat(r.alt).forEach(x => {
      assert.ok(ANALYSIS_KEYS.includes(x.key), x.key);
      Object.keys(x.roles).forEach(rk => assert.ok(ROLE_KEYS[x.key].includes(rk), `${x.key} 에 역할 "${rk}" 없음 (${ROLE_KEYS[x.key]})`));
    });
  });
});

/* ---------- 시트 시뮬레이션: adviseFromVariables ---------- */

test('GAS 시뮬레이션: adviseFromVariables 가 데이터에서 척도와 집단 수를 읽어 추천한다', () => {
  const sample = S.sampleDataset_();
  const data = [sample.header].concat(sample.rows), vars = [['변수명', '라벨', '척도', '값 레이블', '결측값 코드']].concat(sample.varRows);
  const mk = rows => ({ getLastRow: () => rows.length, getDataRange: () => ({ getValues: () => rows.map(r => r.slice()) }) });
  const ss = { getSheetByName: n => (n === '데이터' ? mk(data) : n === '변수정의' ? mk(vars) : null) };
  const ctx = vm.createContext({ console: { log() {}, error() {} }, SpreadsheetApp: { getActiveSpreadsheet: () => ss }, Utilities: {}, Session: {}, PropertiesService: {} });
  for (const f of fs.readdirSync(SRC).filter(x => x.endsWith('.js')).sort()) vm.runInContext(fs.readFileSync(path.join(SRC, f), 'utf8'), ctx, { filename: f });
  const run = names => { ctx.__n = names; return JSON.parse(JSON.stringify(vm.runInContext('adviseFromVariables(__n)', ctx))); };
  const t = run(['점수', '성별']);
  assert.equal(t.key, 'ttestInd');
  assert.deepEqual(t.roles, { dependent: ['점수'], factor: ['성별'] });
  assert.equal(t.n, 60);
  assert.equal(run(['점수', '집단']).key, 'anova');
  assert.equal(run(['성별', '집단']).key, 'crosstab');
  assert.equal(run(['만족도', '성별']).key, 'mannWhitney');
  assert.equal(run(['연령', '점수']).key, 'correlation');
  assert.equal(run(['성별']).key, 'frequencies');
  assert.throws(() => run(['없는변수']), /찾을 수 없습니다/);
  const tree = JSON.parse(JSON.stringify(vm.runInContext('getAdvisor()', ctx)));
  assert.equal(tree.tree.choices.length, 14);
});

/* ---------- 넓게 보는 도우미 대화상자와 분석 창 연결 ---------- */

test('분석 이름표(ADVISOR_LABELS_)는 client.html 의 분석 이름과 같다 (새 분석을 만들면 함께 추가해야 한다)', () => {
  const s = html.indexOf('  var ANALYSIS_GROUPS = '), e = html.indexOf('\n  };\n', html.indexOf('  var ANALYSES = {')) + 5;
  const ctx = vm.createContext({});
  vm.runInContext(html.slice(s, e) + '\nthis.out = ANALYSES;', ctx);
  const want = {};
  Object.keys(ctx.out).forEach(k => { want[k] = ctx.out[k].label; });
  assert.deepEqual(JSON.parse(JSON.stringify(A.ADVISOR_LABELS_)), want);
});

test('advisor.html: 스크립트가 문법 오류 없이 컴파일되고, 부르는 서버 함수가 모두 존재한다', () => {
  const page = fs.readFileSync(path.join(SRC, 'advisor.html'), 'utf8');
  const code = page.match(/<script>([\s\S]*)<\/script>/)[1];
  assert.doesNotThrow(() => new vm.Script(code));
  const calls = [...code.matchAll(/\.(getAdvisor|getDatasetInfo|adviseFromVariables|openAnalysisWith)\(/g)].map(m => m[1]);
  assert.deepEqual([...new Set(calls)].sort(), ['adviseFromVariables', 'getAdvisor', 'getDatasetInfo', 'openAnalysisWith']);
  const ctx = vm.createContext({ console: { log() {}, error() {} }, SpreadsheetApp: {}, Utilities: {}, Session: {}, PropertiesService: {} });
  for (const f of fs.readdirSync(SRC).filter(x => x.endsWith('.js')).sort()) vm.runInContext(fs.readFileSync(path.join(SRC, f), 'utf8'), ctx, { filename: f });
  ['getAdvisor', 'getDatasetInfo', 'adviseFromVariables', 'openAnalysisWith', 'takePendingAdvice', 'openAdvisor'].forEach(fn => assert.equal(vm.runInContext('typeof ' + fn, ctx), 'function', fn));
  // 분석 창도 대기 요청을 받는다
  assert.ok(html.includes('takePendingAdvice') && html.includes('openAdvised(p.key, p.roles)'));
});

test('openAnalysisWith·takePendingAdvice: 요청을 저장하고 분석 창을 열며, 한 번만 돌려주고, 오래된 것·잘못된 것은 버린다', () => {
  const store = {};
  let opened = 0;
  const mkCtx = () => {
    const ctx = vm.createContext({
      console: { log() {}, error() {} }, Utilities: {}, Session: {},
      PropertiesService: { getUserProperties: () => ({ getProperty: k => (k in store ? store[k] : null), setProperty: (k, v) => { store[k] = v; }, deleteProperty: k => { delete store[k]; } }) },
      HtmlService: { createTemplateFromFile: () => ({ evaluate: () => ({ setWidth() { return this; }, setHeight() { return this; } }) }) },
      SpreadsheetApp: { getUi: () => ({ showModelessDialog: () => { opened++; } }) }
    });
    for (const f of fs.readdirSync(SRC).filter(x => x.endsWith('.js')).sort()) vm.runInContext(fs.readFileSync(path.join(SRC, f), 'utf8'), ctx, { filename: f });
    return ctx;
  };
  const ctx = mkCtx();
  ctx.__roles = { dependent: ['점수'], factor: ['성별'], bad: 'x', 'a b': ['q'] };
  vm.runInContext("openAnalysisWith('ttestInd', __roles)", ctx);
  assert.equal(opened, 1);
  const saved = JSON.parse(store['stat4sw.pendingAdvice']);
  assert.deepEqual(saved.roles, { dependent: ['점수'], factor: ['성별'] });                      // 형식이 맞지 않은 역할은 버린다
  const got = JSON.parse(JSON.stringify(vm.runInContext('takePendingAdvice()', ctx)));
  assert.deepEqual(got, { key: 'ttestInd', roles: { dependent: ['점수'], factor: ['성별'] } });
  assert.equal(vm.runInContext('takePendingAdvice()', ctx), null);                              // 한 번만
  assert.throws(() => vm.runInContext("openAnalysisWith('../x')", ctx), /알 수 없는 분석/);
  assert.throws(() => vm.runInContext('openAnalysisWith(5)', ctx), /알 수 없는 분석/);
  store['stat4sw.pendingAdvice'] = JSON.stringify({ key: 'anova', roles: {}, t: Date.now() - 10 * 60 * 1000 });
  assert.equal(vm.runInContext('takePendingAdvice()', ctx), null);                              // 2분이 지나면 버린다
  assert.equal('stat4sw.pendingAdvice' in store, false);
  store['stat4sw.pendingAdvice'] = '{깨진 json';
  assert.equal(vm.runInContext('takePendingAdvice()', ctx), null);
});

test('반복측정 경로: 사전·사후에서 측정 횟수를 묻고, 3회 이상은 Friedman(서열), Cochran Q(예/아니오), 반복측정 ANOVA(숫자, 대안 Friedman)로 안내한다', () => {
  const rep = leaves.filter(l => l.trail.some(t => t.includes('3번 이상')));
  assert.equal(rep.length, 3);
  assert.deepEqual(rep.map(l => l.node.rec.key).sort(), ['cochran', 'friedman', 'rmanova']);
  assert.equal(rep.find(l => l.node.rec.key === 'cochran').node.alt[0].key, 'friedman');   // 0/1 자료의 Friedman 은 Cochran Q 와 같다
  const num = rep.find(l => l.node.rec.key === 'rmanova');
  assert.equal(num.node.rec.parametric, true);
  assert.equal(num.node.alt[0].key, 'friedman');                                           // 숫자는 반복측정 ANOVA(대안 Friedman)
  const two = leaves.filter(l => l.trail.some(t => t.includes('2번 (사전·사후)')));
  assert.deepEqual(two.map(l => l.node.rec.key).sort(), ['mcnemar', 'ttestPaired', 'wilcoxonPaired']);
  const items = [sc('a', 'scale'), sc('b', 'scale'), sc('c', 'scale')];
  const adv = A.adviseFromItems_(items);
  assert.deepEqual(adv.alt.map(x => x.key), ['reliability', 'rmanova', 'friedman']);
  assert.deepEqual(adv.alt[1].roles, { vars: ['a', 'b', 'c'] });
});
