/**
 * 재확인 기록부 생성기: doc/확인사항.md(점검 절차서)의 모든 항목을 처음부터 한 문서에서 다시 확인할 수 있게
 * doc/재확인기록부.md 로 만든다. 항목마다 수행 방법의 줄이 체크칸(☐)이 되고 "재확인 결과" 칸이 붙는다.
 *
 * 실행: node tests/gen_recheck.js   (확인사항.md 를 고친 뒤 다시 실행. 기록부 안의 재확인 결과 칸은 다시 생성하면 지워지므로,
 *       재확인을 시작하면 기록부를 복사해 쓰거나 결과는 확인사항.md의 확인 이력에 옮긴다)
 * tests/recheck.test.js 가 두 문서의 항목이 같은지 확인한다 (항목이 빠지지 않게; 내용이 달라져도 실패하지 않는다, 기록부는 사용자가 복사해 쓰는 용도).
 */
const fs = require('node:fs');
const path = require('node:path');

const SRC = path.join(__dirname, '..', 'doc', '확인사항.md');
const OUT = path.join(__dirname, '..', 'doc', '재확인기록부.md');

/** 확인사항.md → [{ section, id, title, version, status, data, prep, steps: [줄...], history }] */
function parse(text) {
  const lines = text.replace(/\r\n/g, '\n').split('\n');
  const start = lines.findIndex(l => l.startsWith('## 2. 점검 절차'));
  const end = lines.findIndex((l, i) => i > start && l.startsWith('## 3.'));
  const items = [];
  let section = '', cur = null, inSteps = false;
  for (let i = start + 1; i < end; i++) {
    const l = lines[i];
    const sec = l.match(/^### (\S+) (.+)$/);
    if (sec) { section = sec[1] + ' ' + sec[2]; cur = null; continue; }
    const it = l.match(/^#### ([PMUA]-\d+)\. (.+)$/);
    if (it) { cur = { section, id: it[1], title: it[2], version: '', status: '', data: '', prep: '', steps: [], history: '' }; items.push(cur); inSteps = false; continue; }
    if (!cur) continue;
    let m;
    if ((m = l.match(/^- \*\*대상 버전\*\*: (.*?) · \*\*상태\*\*: (.*?) · \*\*시험 데이터\*\*: (.*)$/))) { cur.version = m[1]; cur.status = m[2]; cur.data = m[3]; inSteps = false; }
    else if ((m = l.match(/^- \*\*준비\*\*: (.*)$/))) { cur.prep = m[1]; inSteps = false; }
    else if (/^- \*\*수행 방법과 기대 결과\*\*:/.test(l)) inSteps = true;
    else if ((m = l.match(/^- \*\*확인 이력\*\*: (.*)$/))) { cur.history = m[1]; inSteps = false; }
    else if (inSteps && /^ {2}- /.test(l)) cur.steps.push(l.replace(/^ {2}- /, ''));
    else if (inSteps && /^ {4}/.test(l) && cur.steps.length) cur.steps[cur.steps.length - 1] += ' ' + l.trim();
  }
  return items;
}

function build() {
  const items = parse(fs.readFileSync(SRC, 'utf8'));
  const out = [];
  out.push('# 재확인 기록부', '');
  out.push('이 문서는 [확인사항.md](확인사항.md)(점검 절차서)의 **모든 항목을 처음부터 다시 확인**할 때 쓰는 별도 기록부입니다. 확인사항.md에서 자동으로 만들어지며(`node tests/gen_recheck.js`), 항목마다 수행 방법의 각 줄이 체크칸(☐)이고 맨 아래에 **재확인 결과** 칸이 있습니다.', '');
  out.push('**쓰는 법**: ① 이 파일을 **복사**해서(예: `재확인기록부_2026-11-01.md`) 복사본에 표시하세요 (원본은 다시 생성하면 덮어써집니다). ② 항목마다 시험 데이터와 준비를 따르고 ☐를 ☑로 바꿉니다. ③ 기대 결과와 다르면 결과 칸에 항목 ID와 증상을 적어 알려 주세요. ④ 끝나면 확인사항.md의 확인 이력에 날짜를 남깁니다.', '');
  out.push('상태는 확인사항.md의 값입니다 (✅ = 사용자가 시트에서 확인했다고 알려 준 날짜). 이 문서에는 **' + items.length + '개 항목**이 있습니다.', '');
  out.push('## 항목 목록', '', '| ID | 기능 | 대상 버전 | 상태(최근 확인) | 재확인 날짜 | 결과 |', '| --- | --- | --- | --- | --- | --- |');
  items.forEach(x => out.push('| ' + [x.id, x.title.replace(/\|/g, '/'), x.version.replace(/\|/g, '/'), x.status, '', '☐ 정상 ☐ 문제'].join(' | ') + ' |'));
  out.push('');
  let sec = '';
  items.forEach(x => {
    if (x.section !== sec) { sec = x.section; out.push('## ' + sec, ''); }
    out.push('### ' + x.id + '. ' + x.title, '');
    out.push('- 대상 버전: ' + x.version + ' · 최근 상태: ' + x.status + ' · 시험 데이터: ' + x.data);
    if (x.prep) out.push('- 준비: ' + x.prep);
    out.push('');
    x.steps.forEach(s => out.push('- ☐ ' + s));
    out.push('', '**재확인 결과**: ☐ 모두 기대 결과와 같음  ☐ 다른 점 있음  · 날짜: ____  · 메모: ________________________________', '');
    if (x.history) out.push('> 지금까지의 확인 이력: ' + x.history, '');
  });
  return out.join('\n');
}

if (require.main === module) {
  fs.writeFileSync(OUT, build());
  console.log('wrote', OUT, parse(fs.readFileSync(SRC, 'utf8')).length, '항목');
}
module.exports = { parse, build };
