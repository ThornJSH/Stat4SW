const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const M = require('../src/Mode.js');

const ROOT = path.join(__dirname, '..');
const read = f => fs.readFileSync(path.join(ROOT, f), 'utf8');

/** client.html 의 분석 정의(ANALYSES)만 꺼낸다 */
function analysisKeys() {
  const html = read('src/client.html');
  const start = html.indexOf('  var ANALYSIS_GROUPS = ');
  const end = html.indexOf('\n  };\n', html.indexOf('  var ANALYSES = {')) + 5;
  const ctx = vm.createContext({});
  vm.runInContext(html.slice(start, end) + '\nthis.out = Object.keys(ANALYSES);', ctx);
  return JSON.parse(JSON.stringify(ctx.out));
}

test('일반/고급 분류: client.html 의 모든 분석이 기본 또는 고급 중 한 곳에만 들어 있고, 목록에 없는 키는 없다 (새 분석을 만들면 분류해야 한다)', () => {
  const keys = analysisKeys();
  const basic = M.BASIC_ANALYSES_, adv = M.ADVANCED_ANALYSES_;
  assert.equal(new Set(basic).size, basic.length);
  assert.equal(new Set(adv).size, adv.length);
  keys.forEach(k => assert.equal(basic.includes(k) + adv.includes(k), 1, `${k}: 기본·고급 어느 한 곳에만 있어야 함 → src/Mode.js`));
  basic.concat(adv).forEach(k => assert.ok(keys.includes(k), `${k}: client.html 에 없는 분석 키`));
  assert.equal(basic.length + adv.length, keys.length);
  assert.equal(basic.length, 18);
  assert.equal(adv.length, 19);
});

test('정보·저작권 화면과 LICENSE: 저작권자 표기, MIT 문구, 상표 고지, 직접 구현 안내가 있고 두 곳의 표기가 같다', () => {
  const page = read('src/about.html'), lic = read('LICENSE'), doc = read('doc/저작권·라이선스.md');
  const holder = 'thornjsh.github.io with Claude';
  [page, lic, doc].forEach(t => assert.ok(t.includes(holder), '저작권자 표기'));
  [page, lic].forEach(t => assert.ok(t.includes('Copyright (c) 2026 ' + holder) && t.includes('Permission is hereby granted, free of charge') && t.includes('THE SOFTWARE IS PROVIDED "AS IS"')));
  assert.ok(lic.startsWith('MIT License'));
  assert.ok(page.includes('MIT 라이선스') && page.includes('Google LLC의 상표') && page.includes('외부 라이브러리 없이 직접 구현') && page.includes('<?= APP_VERSION ?>'));
  assert.ok(page.includes('계산 방법의 출처') && page.includes('Royston') && page.includes('Lenth'));
  assert.ok(read('src/Menu.js').includes("addItem('정보·저작권', 'showAbout')") && read('src/Menu.js').includes("createTemplateFromFile('about')"));
});

test('상용 통계 제품 이름 금지: 코드(src, tests)와 사용자에게 보이는 화면에는 특정 제품 이름이 없다 (문서의 사실 비교는 허용)', () => {
  const word = 'SP' + 'SS';
  const files = [];
  const walk = d => fs.readdirSync(path.join(ROOT, d), { withFileTypes: true }).forEach(e => {
    if (e.isDirectory()) { if (!['reference', 'ui', 'node_modules'].includes(e.name)) walk(path.join(d, e.name)); }
    else if (/\.(js|html|py)$/.test(e.name)) files.push(path.join(d, e.name));
  });
  walk('src'); walk('tests');
  const hits = files.filter(f => read(f).includes(word));
  assert.deepEqual(hits, [], '특정 상용 제품 이름이 들어 있음: ' + hits.join(', '));
});

test('고급 분석 표시: 분석 목록은 모드에 따라 걸러지고, 도우미 카드는 잠긴 분석을 열지 않는다 (client.html, advisor.html 코드 확인)', () => {
  const c = read('src/client.html'), a = read('src/advisor.html');
  assert.ok(c.includes('function isLocked(key)') && c.includes('function buildAnalysisList()') && c.includes('mode.advanced = !!info.advanced'));
  assert.ok(c.includes('STAT4SW_MODE_BOOT') && c.includes('.getModeInfo()') && c.includes('function applyMode(m)'));   // 데이터를 못 읽어도 모드가 적용되어야 한다
  assert.ok(read('src/index.html').includes('STAT4SW_MODE_BOOT = <?!= modeBootJson_() ?>'));
  assert.ok(c.includes('if (!ANALYSES[key] || isLocked(key)) return;'));
  assert.ok(a.includes('advisor.advanced === false') && a.includes('advancedKeys'));
  assert.ok(read('src/index.html').includes('id="modeHint"'));
});

/** client.html 의 분석 목록 만들기(buildAnalysisList)를 가짜 화면 요소로 실제 실행한다 */
function runBuildList(advanced) {
  const html = read('src/client.html');
  const defs = html.slice(html.indexOf('  var ANALYSIS_GROUPS = '), html.indexOf('\n  };\n', html.indexOf('  var ANALYSES = {')) + 5);
  const isLocked = html.slice(html.indexOf('  function isLocked(key)'), html.indexOf('\n', html.indexOf('  function isLocked(key)')));
  const start = html.indexOf('  function buildAnalysisList() {');
  const end = html.indexOf('\n  buildAnalysisList();\n', start);
  assert.ok(start >= 0 && end > start);
  const els = {};
  const el = id => els[id] || (els[id] = {
    innerHTML: '', value: '', textContent: '',
    querySelector(sel) { const m = /<option value="([^"]+)"/.exec(this.innerHTML); return m ? { value: m[1] } : null; }
  });
  const ctx = vm.createContext({ $: el, esc: s => String(s) });
  const advKeys = M.ADVANCED_ANALYSES_;
  vm.runInContext(defs + '\nvar state = { analysis: "frequencies" };\nvar mode = { advanced: ' + advanced + ', keys: ' + JSON.stringify(advKeys) + ', hint: "H" };\n' + isLocked + '\n' + html.slice(start, end) +
    '\nthis.run = function (a) { state.analysis = a; buildAnalysisList(); return state.analysis; };\nthis.state = state;', ctx);
  return { ctx, els };
}

test('고급 분석 켜기/끄기: 분석 창의 분석 목록이 모드에 따라 실제로 걸러진다 (꺼짐 18개, 켜짐 37개, 고급 분석을 보던 중 꺼지면 첫 기본 분석으로)', () => {
  const count = els => (els.analysis.innerHTML.match(/<option value="/g) || []).length;
  const off = runBuildList(false);
  const sel = off.ctx.run('frequencies');
  assert.equal(count(off.els), 18);
  M.ADVANCED_ANALYSES_.forEach(k => assert.ok(!off.els.analysis.innerHTML.includes('value="' + k + '"'), k + ' 가 꺼진 목록에 보임'));
  M.BASIC_ANALYSES_.forEach(k => assert.ok(off.els.analysis.innerHTML.includes('value="' + k + '"'), k + ' 가 기본 목록에 없음'));
  assert.equal(sel, 'frequencies');
  ['연구 설계', '생존 분석', '다변량 탐색'].forEach(g => assert.ok(!off.els.analysis.innerHTML.includes('label="' + g), g + ' 묶음이 숨겨지지 않음'));   // 기본 분석이 하나도 없는 묶음은 통째로 숨긴다
  assert.ok(off.els.modeHint.textContent.includes('고급 분석 19개는 숨겨져 있습니다'));
  assert.equal(off.ctx.run('anova2'), 'frequencies');                    // 고급 분석을 보던 중 꺼지면 목록의 첫 분석으로
  assert.equal(off.els.analysis.value, 'frequencies');
  const on = runBuildList(true);
  on.ctx.run('anova2');
  assert.equal(count(on.els), 37);
  assert.equal(on.els.analysis.value, 'anova2');
  assert.ok(on.els.modeHint.textContent.includes('고급 분석이 켜져 있습니다'));
});

test('정보 화면: 새 버전 안내(업데이트 링크 자리), 진단 정보 복사 단추, 서버 함수가 연결돼 있다', () => {
  const page = read('src/about.html'), cfg = read('src/Config.js');
  assert.ok(cfg.includes("var RELEASES_URL_ = '';") || /var RELEASES_URL_ = 'https:\/\//.test(cfg));
  assert.ok(page.includes('<? if (RELEASES_URL_) { ?>') && page.includes('target="_blank" rel="noopener"'));
  assert.ok(page.includes('id="diag"') && page.includes('.getDiagnostics()') && page.includes('execCommand'));
  assert.ok(read('src/Safety.js').includes('function getDiagnostics()'));
  assert.ok(/var APP_VERSION = '\d+\.\d+\.\d+';/.test(cfg));
});
