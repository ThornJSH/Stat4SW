const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const An = require('../src/Stats_Anova.js');
const R = require('../src/Render.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'anova.json'), 'utf8'));

function relErr(a, e) {
  if (a === e) return 0;
  return Math.abs(a - e) / Math.max(Math.abs(e), 1e-300);
}
/** 상대오차 tol, 또는 0 근처 값은 절대오차 absTol */
function close(a, e, tol, msg, absTol = 1e-10) {
  assert.ok(relErr(a, e) < tol || Math.abs(a - e) < absTol, `${msg}: ${a} vs ${e}`);
}

const ds = (values, defs) => Data.buildDataset(values, defs || []);
const tbl = (res, caption) => res.blocks.find(b => b.caption && b.caption.startsWith(caption));

/* ---------- 계산: scipy 대조 ---------- */

test('oneWayAnova: scipy f_oneway 대조 (F, p, 제곱합, 효과크기)', () => {
  for (const c of ref.cases) {
    const r = An.oneWayAnova(c.groups);
    close(r.F, c.F, 1e-9, `${c.name} F`);
    close(r.p, c.p, 1e-8, `${c.name} p`);
    close(r.ssb, c.ssb, 1e-10, `${c.name} ssb`);
    close(r.ssw, c.ssw, 1e-10, `${c.name} ssw`);
    assert.equal(r.dfb, c.dfb);
    assert.equal(r.dfw, c.dfw);
    close(r.eta2, c.eta2, 1e-10, `${c.name} eta2`);
    close(r.omega2, c.omega2, 1e-9, `${c.name} omega2`);
  }
});

test('welchAnova, brownForsythe: 독립 구현 기준값 대조 (k=2는 Welch t² 항등식으로 교차 확인됨)', () => {
  for (const c of ref.cases) {
    const w = An.welchAnova(c.groups), b = An.brownForsythe(c.groups);
    close(w.F, c.welch.F, 1e-9, `${c.name} Welch F`);
    close(w.df2, c.welch.df2, 1e-9, `${c.name} Welch df2`);
    close(w.p, c.welch.p, 1e-8, `${c.name} Welch p`);
    assert.equal(w.df1, c.welch.df1);
    close(b.F, c.bf.F, 1e-9, `${c.name} BF F`);
    close(b.df2, c.bf.df2, 1e-9, `${c.name} BF df2`);
    close(b.p, c.bf.p, 1e-8, `${c.name} BF p`);
  }
});

test('postHoc: Tukey(scipy tukey_hsd), Games-Howell, Bonferroni, Scheffé, LSD, Šidák 기준값 대조', () => {
  for (const c of ref.cases) {
    for (const conf of [95, 99]) {
      for (const method of ['tukey', 'gameshowell', 'bonferroni', 'scheffe', 'lsd', 'sidak']) {
        const got = An.postHoc(method, c.groups, conf);
        const exp = c.posthoc[String(conf)][method];
        assert.equal(got.length, exp.length, `${c.name} ${method}`);
        got.forEach((g, idx) => {
          const e = exp[idx];
          const m = `${c.name} ${method} conf=${conf} (${e.i},${e.j})`;
          assert.equal(g.i, e.i); assert.equal(g.j, e.j);
          close(g.diff, e.diff, 1e-10, m + ' diff');
          close(g.se, e.se, 1e-10, m + ' se');
          // 유의확률은 절대오차 1e-9 (studentized range 누적확률의 정확도), 신뢰구간은 상대 1e-8
          assert.ok(Math.abs(g.p - e.p) < 1e-9, `${m} p: ${g.p} vs ${e.p}`);
          close(g.lo, e.lo, 1e-8, m + ' lo');
          close(g.hi, e.hi, 1e-8, m + ' hi');
        });
      }
    }
  }
});

test('계산할 수 없는 경우는 null (집단 내 변동 없음, 사례 부족)', () => {
  const flat = An.oneWayAnova([[3, 3, 3], [5, 5, 5]]);
  assert.equal(flat.F, null);
  assert.equal(flat.p, null);
  assert.equal(An.welchAnova([[3, 3, 3], [5, 5, 5]]).F, null);
  assert.equal(An.welchAnova([[1], [2, 3, 5]]).F, null);
  assert.equal(An.brownForsythe([[1], [2, 3, 5]]).F, null);
  const one = An.postHoc('tukey', [[3, 3, 3], [5, 5, 5]], 95);
  assert.equal(one[0].p, null);
  const gh = An.postHoc('gameshowell', [[1], [2, 3, 5]], 95);
  assert.equal(gh[0].p, null);
  assert.equal(An.oneWayAnova([[1], [2]]).F, null);   // 오차 자유도 0
});

test('postHoc: 알 수 없는 방법은 오류', () => {
  assert.throws(() => An.postHoc('nope', [[1, 2, 4], [3, 5, 8]], 95), /알 수 없는/);
});

/* ---------- 실행 함수 ---------- */

const defs = Data.parseVarDefs([['h'], ['g', '집단', '명목', '1=A, 2=B, 3=C', '']]).defs;

function sample() {
  return ds([
    ['y', 'z', 'g', 'one'],
    [10, 5, 1, 1], [12, 6, 1, 1], [11, 7, 1, 1], [13, 5, 1, 1],
    [15, 9, 2, 1], [18, 8, 2, 1], [14, 10, 2, 1], [20, 9, 2, 1], [17, 11, 2, 1],
    [22, 12, 3, 1], [25, 14, 3, 1], [21, 13, 3, 1],
    ['', 3, 3, 1], [9, '', 1, 1], [7, 7, '', 1]
  ], defs);
}

const ALL = { levene: true, effect: true, welch: true, tukey: true, gameshowell: true, bonferroni: true, scheffe: true };

test('runOneWayAnova: 기본 표들(기술통계, 분산분석)과 값 레이블, 결측 제외', () => {
  const res = An.runOneWayAnova(sample(), { dependent: ['y'], factor: ['g'], options: { conf: 95 } });
  const d = tbl(res, '기술통계');
  assert.deepEqual(d.rows.map(r => r.group), ['A', 'B', 'C', '전체']);
  assert.deepEqual(d.rows.map(r => r.n), [5, 5, 3, 13]);     // A: 행 14(y=9, g=1) 포함, C: 행 13은 y 결측
  assert.equal(d.rows[3]._kind, 'total');
  const a = tbl(res, '분산분석');
  assert.deepEqual(a.rows.map(r => r.source), ['집단-간', '집단-내', '합계']);
  assert.deepEqual(a.rows.map(r => r.df), [2, 10, 12]);
  assert.ok(Math.abs(a.rows[0].ss + a.rows[1].ss - a.rows[2].ss) < 1e-9);
  assert.ok(a.rows[0].p < 0.05);
  assert.ok(!res.blocks.some(b => b.caption.startsWith('등분산')));           // 옵션을 켜지 않으면 없음
  assert.ok(!res.blocks.some(b => b.caption.startsWith('다중비교')));
});

test('runOneWayAnova: 옵션을 모두 켜면 모든 표가 만들어진다', () => {
  const res = An.runOneWayAnova(sample(), { dependent: ['y', 'z'], factor: ['g'], options: ALL });
  assert.deepEqual(res.blocks.map(b => b.caption), [
    '기술통계 (집단변수: 집단)', '등분산 검정 (Levene, 평균 기준)', '분산분석', '효과크기', '평균 동일성의 강건 검정',
    '다중비교: Tukey HSD', '다중비교: Games-Howell', '다중비교: Bonferroni', '다중비교: Scheffé'
  ]);
  const tukey = tbl(res, '다중비교: Tukey');
  assert.equal(tukey.rows.length, 6);                    // 종속변수 2개 × 쌍 3개
  assert.equal(tukey.rows[0].name, 'y');
  assert.equal(tukey.rows[1].name, '');
  assert.equal(tukey.rows[3].name, 'z');
  assert.deepEqual([tukey.rows[0].gi, tukey.rows[0].gj], ['A', 'B']);
  assert.ok(tukey.rows.some(r => r.sig === '*'));
  assert.equal(tbl(res, '평균 동일성').rows.length, 4);   // 종속변수 2 × (Welch, Brown-Forsythe)
  const html = R.renderResultHtml(res);
  assert.ok(html.includes('다중비교: Games-Howell') && html.includes('Brown-Forsythe') && html.includes('η²'));
});

test('runOneWayAnova: 유의수준·신뢰수준이 열 제목과 별표 기준에 반영', () => {
  const r99 = An.runOneWayAnova(sample(), { dependent: ['y'], factor: ['g'], options: { tukey: true, conf: 99 } });
  const t = tbl(r99, '다중비교');
  assert.ok(t.columns.some(c => c.label === '차이의 99% 신뢰구간 하한'));
  assert.ok(t.notes.some(n => n.includes('0.01 미만')));
  assert.ok(tbl(r99, '기술통계').columns.some(c => c.label === '평균의 99% 신뢰구간 하한'));
});

test('runOneWayAnova: 입력 오류', () => {
  const d = sample();
  assert.throws(() => An.runOneWayAnova(d, { dependent: ['y'], factor: [] }), /변수를 선택/);
  assert.throws(() => An.runOneWayAnova(d, { dependent: ['y'], factor: ['y'] }), /같은 변수/);
  assert.throws(() => An.runOneWayAnova(d, { dependent: ['y'], factor: ['g', 'z'] }), /하나만/);
  assert.throws(() => An.runOneWayAnova(d, { dependent: ['y'], factor: ['one'] }), /1개입니다/);
  assert.throws(() => An.runOneWayAnova(d, { dependent: ['y'], factor: ['g'], options: { conf: 100 } }), /신뢰수준/);
  const rows = [['y', 'g']];
  for (let i = 0; i < 25; i++) rows.push([i, i], [i + 1, i]);
  assert.throws(() => An.runOneWayAnova(ds(rows), { dependent: ['y'], factor: ['g'] }), /25개/);
});

test('runOneWayAnova: 검정할 수 없으면 경고', () => {
  const d = ds([['y', 'g'], [1, 'a'], [1, 'a'], [1, 'a'], [3, 'b'], [3, 'b'], [3, 'b']]);
  const res = An.runOneWayAnova(d, { dependent: ['y'], factor: ['g'], options: { tukey: true, welch: true } });
  assert.ok(res.warnings.length >= 2);
  const a = tbl(res, '분산분석');
  assert.equal(a.rows[0].F, null);
  assert.ok(R.renderResultHtml(res).includes('분산분석'));
});

/* ---------- 사후검정 추가: LSD, Šidák, Dunnett ---------- */

test('postHoc: Dunnett(통제집단 대비) — 이중 적분 정밀값(임계값, 유의확률, 신뢰구간)과 scipy.stats.dunnett(준-몬테카를로)와 일치', () => {
  assert.equal(ref.cases.length, 7);
  for (const c of ref.cases) {
    for (const conf of [95, 99]) {
      for (const ctrl of Object.keys(c.dunnett[String(conf)])) {
        const exp = c.dunnett[String(conf)][ctrl].rows;
        const got = An.postHoc('dunnett', c.groups, conf, Number(ctrl));
        assert.equal(got.length, c.groups.length - 1);
        got.forEach((g, idx) => {
          const e = exp[idx], m = `${c.name} dunnett conf=${conf} 통제=${ctrl} (${e.i} vs ${e.j})`;
          assert.equal(g.i, e.i); assert.equal(g.j, Number(ctrl));
          close(g.diff, e.diff, 1e-10, m + ' diff');
          close(g.se, e.se, 1e-10, m + ' se');
          assert.ok(Math.abs(g.p - e.p) < 1e-9, `${m} p: ${g.p} vs ${e.p}`);
          close(g.lo, e.lo, 1e-8, m + ' lo');
          close(g.hi, e.hi, 1e-8, m + ' hi');
          // scipy 의 값은 준-몬테카를로라 오차가 크다 (유의확률 약 1e-4, 신뢰구간 약 1e-3)
          assert.ok(Math.abs(g.p - e.scipy_p) < 5e-4, `${m} scipy p`);
          if (conf === 95) assert.ok(Math.abs(g.lo - e.scipy_lo) < 5e-3 * Math.max(1, Math.abs(e.lo)), `${m} scipy lo`);   // scipy 의 99% 신뢰구간은 수렴 경고가 나서 95% 만 본다
        });
      }
    }
  }
});

test('postHoc: Dunnett — 집단이 2개면 합동분산 t검정과 같고, Tukey·Bonferroni보다 통제집단 비교의 유의확률이 작거나 같다', () => {
  const g = [[10, 12, 11, 13, 9], [15, 18, 14, 20, 17, 16]];
  const d = An.postHoc('dunnett', g, 95)[0], l = An.postHoc('lsd', g, 95)[0];
  close(d.p, l.p, 1e-12, 'p');
  close(d.lo, -l.hi, 1e-12, 'lo');   // Dunnett 의 차이는 (집단 2 − 집단 1)이라 LSD 쌍(집단 1 − 집단 2)과 부호가 반대
  const c = ref.cases[1];
  const dn = An.postHoc('dunnett', c.groups, 95, 0), bo = An.postHoc('bonferroni', c.groups, 95);
  dn.forEach(r => {
    const same = bo.find(b => b.i === 0 && b.j === r.i);
    assert.ok(r.p <= same.p + 1e-12 && r.p >= An.postHoc('lsd', c.groups, 95).find(b => b.i === 0 && b.j === r.i).p - 1e-12);
  });
});

test('postHoc: LSD·Šidák은 보정 없음/Šidák 보정이고, 계산할 수 없으면 null', () => {
  const flat = [[3, 3, 3], [5, 5, 5], [6, 6, 6]];
  ['lsd', 'sidak', 'dunnett'].forEach(m => assert.equal(An.postHoc(m, flat, 95)[0].p, null));
  assert.throws(() => An.postHoc('dunnett', [[1, 2, 4], [3, 5, 8]], 95, 5), /통제집단 번호/);
  const c = ref.cases[0];
  const lsd = An.postHoc('lsd', c.groups, 95), sid = An.postHoc('sidak', c.groups, 95), bon = An.postHoc('bonferroni', c.groups, 95);
  lsd.forEach((r, i) => { assert.ok(r.p <= sid[i].p && sid[i].p <= bon[i].p + 1e-12); assert.ok(r.hi - r.lo <= sid[i].hi - sid[i].lo); });
});

test('runOneWayAnova: LSD·Šidák·Dunnett 옵션 — 표 이름, 통제집단(비움=첫 집단, 값·레이블로 지정), 오류', () => {
  const base = { dependent: ['y', 'z'], factor: ['g'], options: { lsd: true, sidak: true, dunnett: true, conf: 95 } };
  const res = An.runOneWayAnova(sample(), base);
  assert.deepEqual(res.blocks.map(b => b.caption).filter(c => c.startsWith('다중비교')), ['다중비교: LSD (최소유의차)', '다중비교: Šidák', '다중비교: Dunnett (통제집단 대비)']);
  const dn = tbl(res, '다중비교: Dunnett');
  assert.equal(dn.rows.length, 4);                       // 종속변수 2개 × (집단 3개 − 통제 1개)
  assert.deepEqual([dn.rows[0].gi, dn.rows[0].gj, dn.rows[1].gi, dn.rows[1].gj], ['B', 'A', 'C', 'A']);
  assert.ok(dn.notes.some(n => /통제집단은 "A"\(J\)이다 \(옵션을 비워 첫 집단을 썼다\)/.test(n)));
  assert.ok(tbl(res, '다중비교: LSD').notes.some(n => /모든 쌍을 한 번씩/.test(n)));
  const byLabel = An.runOneWayAnova(sample(), Object.assign({}, base, { options: { dunnett: true, control: 'C' } }));
  const d2 = tbl(byLabel, '다중비교: Dunnett');
  assert.deepEqual([d2.rows[0].gi, d2.rows[0].gj, d2.rows[1].gi], ['A', 'C', 'B']);
  assert.ok(d2.notes.some(n => /통제집단은 "C"\(J\)이다\./.test(n)));
  const byValue = An.runOneWayAnova(sample(), Object.assign({}, base, { options: { dunnett: true, control: '2' } }));
  assert.equal(tbl(byValue, '다중비교: Dunnett').rows[0].gj, 'B');
  assert.throws(() => An.runOneWayAnova(sample(), Object.assign({}, base, { options: { dunnett: true, control: '없음' } })), /통제집단 "없음".*A, B, C/);
  // control 은 dunnett 를 켤 때만 해석한다
  assert.doesNotThrow(() => An.runOneWayAnova(sample(), Object.assign({}, base, { options: { tukey: true, control: '없음' } })));
});

test('큰 문제: Dunnett 20집단(크기가 모두 다름) × 5만 건이 Node 에서 10초 안에 끝난다', () => {
  let seed = 7;
  const rnd = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };
  const groups = [];
  for (let k = 0; k < 20; k++) {
    const n = 1500 + k * 53, arr = [];
    for (let i = 0; i < n; i++) arr.push(k * 0.02 + (rnd() + rnd() + rnd() - 1.5));
    groups.push(arr);
  }
  const t0 = Date.now();
  const rows = An.postHoc('dunnett', groups, 95, 0);
  const sec = (Date.now() - t0) / 1000;
  assert.equal(rows.length, 19);
  assert.ok(rows.every(r => r.p >= 0 && r.p <= 1 && r.lo < r.hi));
  assert.ok(sec < 10, `Dunnett ${sec}초`);
});
