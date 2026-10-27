const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const { evaluateTemplateReadiness_: check } = require('../src/Release.js');

const clean = () => ({
  hasData: true, hasVars: true, hasResult: false,
  dataHeaderCells: 0, dataRows: 0,
  varHeader: ['변수명', '라벨', '척도', '값 레이블', '결측값 코드'], varRows: 0,
  otherTabs: []
});

test('evaluateTemplateReadiness_: 비워 둔 표준 템플릿은 문제 없음', () => {
  assert.deepEqual(check(clean()), { ok: true, problems: [] });
});

test('evaluateTemplateReadiness_: 남아 있는 데이터, 변수정의, 결과 탭, 다른 탭을 모두 알려 준다', () => {
  const s = Object.assign(clean(), { dataHeaderCells: 5, dataRows: 60, varRows: 5, hasResult: true, otherTabs: ['시트1', '메모'] });
  const r = check(s);
  assert.equal(r.ok, false);
  assert.equal(r.problems.length, 5);
  const text = r.problems.join('\n');
  assert.ok(text.includes('데이터가 60행') && text.includes('5칸') && text.includes('변수 정의가 5개'));
  assert.ok(text.includes('"결과" 탭이 있습니다') && text.includes('시트1, 메모'));
});

test('evaluateTemplateReadiness_: 필수 탭 누락과 머리글 불일치', () => {
  assert.deepEqual(check(Object.assign(clean(), { hasData: false })).problems.length, 1);
  const noVars = check(Object.assign(clean(), { hasVars: false, varHeader: [] }));
  assert.equal(noVars.problems.length, 1);                        // 변수정의가 없으면 머리글 검사는 건너뛴다
  assert.ok(noVars.problems[0].includes('"변수정의" 탭이 없습니다'));
  const bad = check(Object.assign(clean(), { varHeader: ['이름', '라벨', '척도', '값 레이블', '결측값 코드'] }));
  assert.equal(bad.ok, false);
  assert.ok(bad.problems[0].includes('표준 머리글'));
  const spaced = check(Object.assign(clean(), { varHeader: [' 변수명 ', '라벨', '척도', '값 레이블', '결측값 코드'] }));
  assert.equal(spaced.ok, true);                                  // 앞뒤 공백은 허용
});

test('사용 안내 문서(help.html): 템플릿 표현식은 버전 하나뿐이고 핵심 안내가 들어 있다', () => {
  const html = fs.readFileSync(path.join(__dirname, '..', 'src', 'help.html'), 'utf8');
  assert.deepEqual(html.match(/<\?[^>]*\?>/g), ['<?= APP_VERSION ?>']);
  for (const phrase of ['데이터', '변수정의', '분석 창 열기', '시트로 내보내기', '익명화', '권한 승인', '예제 데이터 채우기']) {
    assert.ok(html.includes(phrase), phrase);
  }
});
