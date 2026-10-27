const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const D = require('../src/Distributions.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'distributions.json'), 'utf8'));

// p 함수는 상대 1e-10, 역함수는 상대 1e-9 (꼬리 확률이 아주 작아도 상대 오차로 비교)
// t, F 는 자유도가 아주 크면 불완전 베타의 지수부 상쇄로 오차가 df에 비례해 커진다 (≈ 3e-17 × df).
// df ≤ 1e6 에서는 기본 허용오차 안에 들어오고, 그 이상은 5e-17 × df 까지 허용한다.
const TOL = { p: 1e-10, q: 1e-9 };
function tolFor(fn, args) {
  const base = fn.startsWith('p') ? TOL.p : TOL.q;
  const dfArgs = { pt: [1], qt: [1], pf: [1, 2], qf: [1, 2] }[fn];
  if (!dfArgs) return base;
  const dfMax = Math.max(...dfArgs.map(i => args[i]));
  return Math.max(base, 5e-17 * dfMax);
}

function relErr(actual, expected) {
  if (actual === expected) return 0;
  return Math.abs(actual - expected) / Math.max(Math.abs(expected), 1e-300);
}

for (const fn of ['pnorm', 'qnorm', 'pt', 'qt', 'pchisq', 'qchisq', 'pf', 'qf']) {
  test(`scipy 기준값 대조: ${fn}`, () => {
    const cases = ref.cases.filter(c => c.fn === fn);
    assert.ok(cases.length > 0);
    const bad = [];
    let worst = 0;
    for (const c of cases) {
      const actual = D[fn](...c.args, c.lower);
      const e = relErr(actual, c.expected);
      worst = Math.max(worst, e);
      if (!(e <= tolFor(fn, c.args))) bad.push(`${fn}(${c.args.join(', ')}, lower=${c.lower}) = ${actual}, scipy ${c.expected}, 상대오차 ${e.toExponential(2)}`);
    }
    assert.equal(bad.length, 0, `${bad.length}/${cases.length}건 실패 (최대 상대오차 ${worst.toExponential(2)})\n` + bad.slice(0, 12).join('\n'));
  });
}

test('잘 알려진 값', () => {
  const close = (a, b, tol = 1e-12) => assert.ok(relErr(a, b) < tol, `${a} vs ${b}`);
  close(D.qnorm(0.975), 1.959963984540054);
  close(D.qt(0.975, 10), 2.228138851986274);
  close(D.qchisq(0.95, 1), 3.841458820694124);
  close(D.qf(0.95, 3, 10), 3.7082648190468435);
  close(D.pValueT(2.228138851986274, 10), 0.05, 1e-10);
  close(D.pValueChisq(3.841458820694124, 1), 0.05, 1e-10);
});

test('경계와 잘못된 입력', () => {
  assert.equal(D.qnorm(0), -Infinity);
  assert.equal(D.qnorm(1), Infinity);
  assert.equal(D.qnorm(0.5), 0);
  assert.equal(D.qt(0.5, 7), 0);
  assert.equal(D.qt(0, 5), -Infinity);
  assert.equal(D.qchisq(0, 3), 0);
  assert.equal(D.qf(1, 3, 4), Infinity);
  assert.equal(D.pt(Infinity, 5), 1);
  assert.equal(D.pt(-Infinity, 5), 0);
  assert.equal(D.pchisq(0, 3), 0);
  assert.equal(D.pchisq(-1, 3, false), 1);
  assert.equal(D.pf(0, 3, 4), 0);
  assert.equal(D.pf(Infinity, 3, 4, false), 0);
  assert.ok(Number.isNaN(D.pt(1, 0)));
  assert.ok(Number.isNaN(D.pchisq(1, -1)));
  assert.ok(Number.isNaN(D.pf(1, 0, 3)));
  assert.ok(Number.isNaN(D.qnorm(1.5)));
  assert.ok(Number.isNaN(D.qt(NaN, 3)));
});

test('p함수와 q함수는 서로 역함수', () => {
  for (const p of [1e-9, 0.01, 0.3, 0.5, 0.8, 0.999]) {
    assert.ok(relErr(D.pnorm(D.qnorm(p)), p) < 1e-12);
    assert.ok(relErr(D.pt(D.qt(p, 7), 7), p) < 1e-12);
    assert.ok(relErr(D.pchisq(D.qchisq(p, 4), 4), p) < 1e-12);
    assert.ok(relErr(D.pf(D.qf(p, 3, 9), 3, 9), p) < 1e-12);
    assert.ok(relErr(D.pf(D.qf(p, 3, 9, false), 3, 9, false), p) < 1e-12);
  }
});

test('scipy 역함수가 부정확해 기준값에서 뺀 영역도 p함수로 되돌려 확인 (F 위꼬리)', () => {
  // 기준값(scipy sf)과 1e-10 안에서 일치하는 pf를 믿고, qf(p) → pf 가 p로 돌아오는지 본다.
  const pairs = [[1, 1], [1, 10], [2, 5], [3, 30], [5, 100], [10, 10], [30, 5], [100, 1000], [1, 1e5], [1000, 1000], [0.5, 3]];
  for (const [d1, d2] of pairs) {
    for (const p of [1e-12, 1e-8, 1e-5, 0.01]) {
      const x = D.qf(p, d1, d2, false);
      assert.ok(relErr(D.pf(x, d1, d2, false), p) < 1e-11, `qf(${p}, ${d1}, ${d2}, upper) = ${x}`);
    }
  }
  // F(1,1) = t² 이므로 위꼬리는 (2/π)·atan(1/√x) 로 닫힌 형태가 있다.
  const x = D.qf(1e-12, 1, 1, false);
  assert.ok(relErr((2 / Math.PI) * Math.atan(1 / Math.sqrt(x)), 1e-12) < 1e-10);
});

test('lower와 upper 꼬리의 합은 1, 대칭성', () => {
  for (const [t, df] of [[0.3, 4], [2.5, 12], [-1.7, 30]]) {
    assert.ok(Math.abs(D.pt(t, df) + D.pt(t, df, false) - 1) < 1e-15);
    assert.ok(Math.abs(D.pt(t, df) - D.pt(-t, df, false)) < 1e-15);
  }
  assert.ok(Math.abs(D.pnorm(1.3) + D.pnorm(1.3, false) - 1) < 1e-15);
});

/* ---------- Studentized range (Tukey) ---------- */

const tukeyRef = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'tukey.json'), 'utf8'));

test('ptukey: scipy studentized_range 기준값 대조 (아래·위꼬리 절대오차 1e-10)', () => {
  let worst = 0;
  for (const c of tukeyRef.cdf) {
    const lo = Math.abs(D.ptukey(c.q, c.k, c.df) - c.lower);
    const up = Math.abs(D.ptukey(c.q, c.k, c.df, false) - c.upper);
    worst = Math.max(worst, lo, up);
    assert.ok(lo < 1e-10 && up < 1e-10, `ptukey(${c.q}, ${c.k}, ${c.df}): 오차 ${Math.max(lo, up).toExponential(2)}`);
  }
  assert.ok(worst < 1e-10);
});

test('qtukey: scipy 기준값 대조 (상대오차 1e-9, 3건 중 1건 — 한 건에 약 50ms가 들어 전체를 돌리면 테스트가 느려진다)', () => {
  for (const c of tukeyRef.ppf.filter((_, i) => i % 3 === 0)) {
    const x = D.qtukey(c.p, c.k, c.df);
    assert.ok(relErr(x, c.q) < 1e-9, `qtukey(${c.p}, ${c.k}, ${c.df}) = ${x}, scipy ${c.q}`);
  }
});

test('qtukey: 위꼬리 입력과 ptukey 왕복', () => {
  for (const [k, df] of [[3, 12], [6, 40], [10, 5]]) {
    for (const p of [0.01, 0.05, 0.2]) {
      const q = D.qtukey(p, k, df, false);
      assert.ok(Math.abs(D.ptukey(q, k, df, false) - p) < 1e-11, `k=${k} df=${df} p=${p}`);
      assert.ok(relErr(q, D.qtukey(1 - p, k, df)) < 1e-9);
    }
  }
});

test('ptukey: 항등식과 경계', () => {
  // k = 2 이면 Q = √2·|T|:  P(Q ≤ q) = 1 − 2·P(T > q/√2)
  for (const df of [3, 10, 50]) {
    for (const q of [0.5, 2, 4]) {
      assert.ok(Math.abs(D.ptukey(q, 2, df) - (1 - 2 * D.pt(q / Math.SQRT2, df, false))) < 1e-12);
    }
  }
  // 자유도 무한대: k = 2 이면 Q = √2·|Z|
  assert.ok(Math.abs(D.ptukey(3, 2, Infinity) - (2 * D.pnorm(3 / Math.SQRT2) - 1)) < 1e-12);
  assert.equal(D.ptukey(0, 3, 10), 0);
  assert.equal(D.ptukey(-1, 3, 10), 0);
  assert.equal(D.ptukey(Infinity, 3, 10), 1);
  assert.equal(D.qtukey(0, 3, 10), 0);
  assert.equal(D.qtukey(1, 3, 10), Infinity);
  assert.ok(Number.isNaN(D.ptukey(1, 1, 10)));   // k < 2
  assert.ok(Number.isNaN(D.ptukey(1, 3, 0)));
  assert.ok(Number.isNaN(D.qtukey(1.5, 3, 10)));
  // 집단 수가 많을수록 같은 q 에서 누적확률은 작아진다
  assert.ok(D.ptukey(3, 3, 20) > D.ptukey(3, 6, 20));
});
