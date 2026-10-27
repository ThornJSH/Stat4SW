/**
 * 예제 데이터를 CSV 파일로 저장한다 (시트의 "예제 데이터 채우기"를 쓰지 않고 파일로 가져오고 싶을 때).
 * 실행: node tests/gen_sample_csv.js
 * 결과: doc/예제데이터/예제_데이터.csv, doc/예제데이터/예제_변수정의.csv (UTF-8 BOM, 쉼표 구분)
 * Sample.js 가 바뀌면 다시 실행한다. tests/sample.test.js 가 파일이 최신인지 확인한다.
 */
const fs = require('node:fs');
const path = require('node:path');
const { sampleDataset_ } = require('../src/Sample.js');

const VAR_HEADER = ['변수명', '라벨', '척도', '값 레이블', '결측값 코드'];

function csvCell(v) {
  const s = v === null || v === undefined ? '' : String(v);
  return /[",\r\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
}
function toCsv(rows) {
  return '﻿' + rows.map(r => r.map(csvCell).join(',')).join('\r\n') + '\r\n';
}

function buildSampleCsv() {
  const s = sampleDataset_();
  return {
    data: toCsv([s.header].concat(s.rows)),
    vars: toCsv([VAR_HEADER].concat(s.varRows))
  };
}

if (require.main === module) {
  const dir = path.join(__dirname, '..', 'doc', '예제데이터');
  fs.mkdirSync(dir, { recursive: true });
  const csv = buildSampleCsv();
  fs.writeFileSync(path.join(dir, '예제_데이터.csv'), csv.data);
  fs.writeFileSync(path.join(dir, '예제_변수정의.csv'), csv.vars);
  console.log('wrote', dir);
}

module.exports = { buildSampleCsv, toCsv };
