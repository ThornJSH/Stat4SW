const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const G = require('./gen_recheck.js');

const DOC = path.join(__dirname, '..', 'doc');
const src = fs.readFileSync(path.join(DOC, '확인사항.md'), 'utf8');

test('재확인 기록부: 확인사항.md 요약표의 모든 항목이 기록부에 있고 순서와 이름이 같다', () => {
  const tableIds = [...src.matchAll(/^\| ([PMUA]-\d+) \|/gm)].map(m => m[1]);
  const items = G.parse(src);
  assert.deepEqual(items.map(x => x.id).sort(), tableIds.slice().sort(), '요약표와 절차의 항목이 다름');
  assert.ok(items.length >= 60);
  items.forEach(x => assert.ok(x.title && x.version && x.status && x.steps.length >= 1, `${x.id}: 파싱 실패 (대상 버전·상태·수행 방법이 있어야 함)`));
});

test('재확인 기록부: doc/재확인기록부.md 에 모든 항목이 들어 있다 (확인사항.md 를 고쳤다면 node tests/gen_recheck.js 로 다시 만든다)', () => {
  const doc = fs.readFileSync(path.join(DOC, '재확인기록부.md'), 'utf8');
  G.parse(src).forEach(x => assert.ok(doc.includes('### ' + x.id + '. ' + x.title), x.id + ' 가 기록부에 없음 → node tests/gen_recheck.js'));
});
