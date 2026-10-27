const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

const read = f => fs.readFileSync(path.join(__dirname, '..', f), 'utf8');

/** client.html 의 "화면 도우미" 구간(순수 함수)만 꺼내 실행한다 */
function loadHelpers() {
  const html = read('src/client.html');
  const a = html.indexOf('/* ---------- 화면 도우미');
  const b = html.indexOf('/* ---------- /화면 도우미 ---------- */');
  assert.ok(a >= 0 && b > a, '화면 도우미 구간을 찾지 못함');
  const ctx = vm.createContext({});
  vm.runInContext(html.slice(a, b) + '\nthis.out = { friendlyError: friendlyError, matchVar: matchVar, START: START_GUIDE_HTML, NEWWORK: NEW_WORK_GUIDE_HTML };', ctx);
  return ctx.out;
}
const H = loadHelpers();

test('friendlyError: 구글이 영어로 돌려주는 오류를 한글 안내로 바꾸고, 모르는 오류는 그대로 둔다', () => {
  const cases = [
    ['Exceeded maximum execution time', /계산 시간이 길어 중단.*약 6분.*줄이거나/],
    ['Authorization is required to perform that action.', /권한 승인이 필요/],
    ['You do not have permission to call SpreadsheetApp.openById', /권한 승인이 필요/],
    ['Service Spreadsheets timed out while accessing document with id abc', /시간 초과/],
    ['Service invoked too many times for one day: premium', /너무 많이 호출/],
    ['Failed to fetch', /네트워크 연결/],
    ['Script function not found: runAnalysis', /버전이 서로 맞지 않습니다/]
  ];
  cases.forEach(([raw, re]) => {
    assert.match(H.friendlyError(new Error(raw)), re, raw);
    assert.match(H.friendlyError({ message: raw }), re);
    assert.match(H.friendlyError(raw), re);
  });
  assert.equal(H.friendlyError(new Error('집단변수는 하나만 선택하세요.')), '집단변수는 하나만 선택하세요.');   // 사용자가 고칠 수 있는 한글 오류는 그대로
  assert.equal(H.friendlyError(null), '');
  assert.equal(H.friendlyError(undefined), '');
  assert.equal(H.friendlyError({}), '[object Object]');
});

test('matchVar: 변수 이름이나 라벨에 검색어가 들어 있으면 보이고, 대소문자·앞뒤 공백을 무시하며 빈 검색어는 모두 보인다', () => {
  const v = { name: 'Q12_만족', label: '서비스 만족도(5점)' };
  assert.ok(H.matchVar(v, ''));
  assert.ok(H.matchVar(v, '   '));
  assert.ok(H.matchVar(v, undefined));
  assert.ok(H.matchVar(v, 'q12'));
  assert.ok(H.matchVar(v, ' Q12_ '));
  assert.ok(H.matchVar(v, '만족'));
  assert.ok(H.matchVar(v, '서비스'));          // 라벨로도 찾는다
  assert.ok(!H.matchVar(v, '불안'));
  assert.ok(H.matchVar({ name: 'age' }, 'AG'));  // 라벨이 없어도 된다
  assert.ok(!H.matchVar({ name: 'age' }, 'x'));
});

test('시작 안내 문구: 데이터가 없을 때는 3단계와 예제 데이터, 탭이 없을 때는 새 작업 시작을 알려 준다', () => {
  assert.ok(H.START.includes('처음 시작하기') && H.START.includes('데이터</b> 탭 1행') && H.START.includes('예제 데이터 채우기') && H.START.includes('새로고침'));
  assert.ok(H.NEWWORK.includes('새 작업 시작 (데이터 모두 지우기)') && H.NEWWORK.includes('만들어집니다') && H.NEWWORK.includes('내용이 있으면 지워지므로'));
});

test('화면 연결: 변수 검색, 시작 안내, 오류 한글화, 결과 자동 저장, Enter 실행, 접근성 속성이 화면 코드에 연결돼 있다', () => {
  const c = read('src/client.html'), i = read('src/index.html'), s = read('src/style.html');
  assert.ok(c.includes('matchVar(v, q)') && c.includes("$('varSearch').addEventListener('input', renderVars)") && c.includes('검색 결과가 없습니다'));
  assert.ok(c.includes('info.n === 0 ? START_GUIDE_HTML') && c.includes('NEW_WORK_GUIDE_HTML') && c.includes('friendlyError(err)'));
  assert.ok(c.includes("if ($('autoExport').checked) exportIds([String(resultCount)])") && c.includes('saved.ui.autoExport'));
  assert.ok(c.includes("e.key === 'Enter'") && c.includes('run();'));
  ['id="varSearch"', 'id="autoExport"', 'aria-label="분석 선택"', 'role="status" aria-live="polite"', 'aria-label="변수 목록"'].forEach(x => assert.ok(i.includes(x), x));
  assert.ok(s.includes(':focus-visible'));
});

test('옵션 기억: 결과 자동 저장 설정(ui)은 저장·복원 대상이고 모르는 모양은 버려진다', () => {
  const S = require('../src/Safety.js');
  assert.deepEqual(S.sanitizeSavedOptions_({ ui: { autoExport: true }, anova: { tukey: true } }), { ui: { autoExport: true }, anova: { tukey: true } });
  assert.deepEqual(S.sanitizeSavedOptions_({ ui: { autoExport: { x: 1 } } }), { ui: {} });
});
