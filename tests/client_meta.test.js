const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

/** client.html 의 분석 정의(ANALYSES)와 묶음(ANALYSIS_GROUPS)만 꺼내 평가한다 (화면 코드는 실행하지 않는다) */
function loadDefs() {
  const html = fs.readFileSync(path.join(__dirname, '..', 'src', 'client.html'), 'utf8');
  const start = html.indexOf('  var ANALYSIS_GROUPS = ');
  const end = html.indexOf('\n  };\n', html.indexOf('  var ANALYSES = {')) + 5;
  assert.ok(start >= 0 && end > start);
  const ctx = vm.createContext({});
  vm.runInContext(html.slice(start, end) + '\nthis.out = { ANALYSIS_GROUPS: ANALYSIS_GROUPS, ANALYSES: ANALYSES };', ctx);
  return JSON.parse(JSON.stringify(ctx.out));
}
const { ANALYSIS_GROUPS, ANALYSES } = loadDefs();

test('분석 선택 화면: 모든 분석이 목록의 묶음에 속하고, 언제 쓰는지 한 줄 설명이 있다', () => {
  assert.equal(new Set(ANALYSIS_GROUPS).size, ANALYSIS_GROUPS.length);
  Object.keys(ANALYSES).forEach(k => {
    const a = ANALYSES[k];
    assert.ok(ANALYSIS_GROUPS.includes(a.group), k + ': group 이 ANALYSIS_GROUPS 에 없음 (' + a.group + ')');
    assert.ok(typeof a.desc === 'string' && a.desc.length >= 15 && a.desc.length <= 110, k + ': desc 는 15~110자 한 줄 (' + (a.desc || '').length + ')');
  });
  ANALYSIS_GROUPS.forEach(g => assert.ok(Object.values(ANALYSES).some(a => a.group === g), '빈 묶음: ' + g));
});

test('분석 선택 화면: 세부 옵션은 키가 겹치지 않고, 기본 옵션이 하나 이상 남거나 옵션이 없는 분석만 있다', () => {
  Object.keys(ANALYSES).forEach(k => {
    const keys = ANALYSES[k].options.map(o => o.key);
    assert.equal(new Set(keys).size, keys.length, k + ': 옵션 키 중복');
    const basic = ANALYSES[k].options.filter(o => !o.adv);
    // 기본 옵션이 전혀 없으면 "옵션" 제목 아래에 접힌 세부 옵션만 남아 어색하다 (신뢰수준만 있는 분석은 예외로 두지 않는다)
    if (ANALYSES[k].options.length) assert.ok(basic.length >= 1, k + ': 기본 옵션이 없음');
  });
  // 저장·확인창이 걸린 옵션은 숨기지 않는다
  assert.ok(!ANALYSES.logistic.options.find(o => o.key === 'save').adv);
});

test('분석 선택 화면: 선택 목록(select) 옵션은 기본값이 선택지에 있고, 거리 측도·표준화 선택지는 서버가 아는 키다', () => {
  const D = require('../src/Stats_Distance.js');
  Object.keys(ANALYSES).forEach(k => ANALYSES[k].options.filter(o => o.type === 'select').forEach(o => {
    assert.ok(o.choices.length >= 2 && o.choices.some(c => c[0] === o.def), k + '.' + o.key + ': 기본값이 선택지에 없음');
    assert.ok(o.def.length <= 20 && o.choices.every(c => c[0].length <= 20), k + '.' + o.key + ': 저장할 수 없는 긴 값');   // 옵션 기억은 20자 이하 문자열만 저장한다
  }));
  const pick = key => ANALYSES.distance.options.find(o => o.key === key).choices.map(c => c[0]);
  assert.deepEqual(pick('measure').slice().sort(), Object.keys(D.DIST_MEASURES_).sort());
  assert.deepEqual(pick('standardize').slice().sort(), Object.keys(D.DIST_STD_).sort());
});
