/**
 * 시험 데이터 레지스트리 생성기: 기준값 JSON(tests/reference/*.json, 파이썬이 독립 계산한 값)의 자료와 기대 결과를
 * 시트의 메뉴 "Stat4SW 통계분석 > 준비·연습 > 시험 데이터"로 넣을 수 있는 src/TestData.js 로 바꾼다.
 *
 * 실행: node tests/gen_test_datasets.js   (새 기법을 만들 때 기준값 JSON을 만든 뒤 이 파일에 항목을 추가하고 다시 실행)
 *
 * 새 기법을 추가할 때 지켜야 할 것 (tests/testdata.test.js 가 강제한다):
 *  1. client.html 의 분석 키마다 여기에 시험 데이터 항목이 있거나, 기본 예제로 충분하면 SAMPLE_COVERED 에 키를 적는다.
 *  2. 항목의 expected 는 파이썬 기준값에서 가져온 숫자이고, 테스트가 이 도구의 계산 결과와 자동으로 비교한다.
 *  3. guide 는 분석 창에서 무엇을 넣고 어떤 옵션을 켜는지 사용자가 따라 할 수 있게 적는다.
 */
const fs = require('node:fs');
const path = require('node:path');

const REF = name => JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', name), 'utf8'));
const OUT = path.join(__dirname, '..', 'src', 'TestData.js');

/** 기본 예제(Sample.js)로 시험하는 분석 (이미 시트에서 확인됨). 새 분석은 여기에 넣지 말고 항목을 만든다. */
/** 시트의 데이터가 필요 없는 분석 (client.html 의 noData). 시험은 점검 절차서의 기대값으로 한다 (테스트는 tests/power.test.js). */
const NO_DATA = ['power'];
const SAMPLE_COVERED = ['frequencies', 'descriptives', 'normality', 'levene', 'ttest1', 'ttestInd', 'ttestPaired', 'anova', 'correlation', 'crosstab', 'regression', 'wilcoxonPaired'];

const ANALYSIS_LABELS = {
  mannWhitney: 'Mann-Whitney U 검정', wilcoxonOne: 'Wilcoxon 부호순위 검정 (일표본)', kruskalWallis: 'Kruskal-Wallis 검정',
  reliability: '신뢰도 분석 (크론바흐 α)', moderation: '조절효과 분석', mediation: '매개효과 분석', logistic: '로지스틱 회귀', distance: '거리·유사도 (마할라노비스 이상치)', psm: '성향점수 매칭', kmeans: 'k-평균 군집분석', cart: '의사결정나무 (CART)', mcnemar: 'McNemar 검정', ttestInd: '독립표본 t검정', ttestPaired: '대응표본 t검정', ttest1: '일표본 t검정', friedman: 'Friedman 검정', rmanova: '반복측정 분산분석 (집단변수는 선택)', kappa: '평가자 일치도 (kappa)', icc: '급내상관 (ICC)', ologit: '서열 로지스틱 회귀 (비례오즈)', correlation: '상관분석 (편상관, 켄달 타우-b)', km: 'Kaplan–Meier 생존분석 (log-rank)', anova: '일원 분산분석 (LSD·Šidák·Dunnett 사후검정)', binom: '이항 검정 (성공 비율 검정)', cochran: 'Cochran Q 검정 (예/아니오, 3회 이상 반복측정)', ipa: 'IPA 중요도–수행도 분석', hierarchical: '위계적 회귀 (단계별 변수 추가)', efa: '탐색적 요인분석 (EFA)', anova2: '이원 분산분석 (요인 2개)'
};

const entries = [];
const fmt = v => (v === null || v === undefined ? '' : Number(v.toPrecision(6)).toString());

/** 기대값 한 개: table(표 제목의 앞부분), row(0부터), key(열 키), value, abs(부호 무시), info(참고용, 자동 비교하지 않음) */
function ex(label, table, row, key, value, extra) {
  return Object.assign({ label: label, table: table, row: row, key: key, value: value, tol: 1e-6 }, extra || {});
}
function rowsOf(columns) {
  const n = columns[0].length, out = [];
  for (let i = 0; i < n; i++) out.push(columns.map(c => (c[i] === null || c[i] === undefined ? null : c[i])));
  return out;
}

/* ---------- Mann-Whitney U ---------- */
const np = REF('nonparam.json');
[[0, '작은 표본, 동순위 없음(정확 유의확률도 나옴)'], [4, '동순위가 많은 정수 자료(점근 유의확률만 나옴)']].forEach(([idx, desc], k) => {
  const c = np.mw[idx];
  const n1 = c.a.length, n2 = c.b.length;
  const group = c.a.map(() => 1).concat(c.b.map(() => 2));
  const u = Math.min(c.u1, n1 * n2 - c.u1);
  const expected = [ex('Mann-Whitney U', 'Mann-Whitney U 검정 통계량', 0, 'u', u), ex('점근 유의확률(양측)', 'Mann-Whitney U 검정 통계량', 0, 'p', c.pAsym)];
  if (c.pExact !== null) expected.push(ex('정확 유의확률(양측)', 'Mann-Whitney U 검정 통계량', 0, 'pExact', c.pExact));
  entries.push({
    id: 'mw-' + (k + 1), analysis: 'mannWhitney', title: 'Mann-Whitney U: ' + c.name, description: desc,
    header: ['집단', '값'], rows: group.map((g, i) => [g, c.a.concat(c.b)[i]]),
    varRows: [['집단', '집단', '명목', '1=집단1, 2=집단2', ''], ['값', '측정값', '연속', '', '']],
    request: { analysis: 'mannWhitney', dependent: ['값'], factor: ['집단'], options: { effect: true } },
    guide: ['분석: Mann-Whitney U 검정', '종속변수: 값', '집단변수: 집단 (집단 2개)', '옵션: 기본 (효과크기 r 켬)'],
    expected: expected
  });
});

/* ---------- Wilcoxon 부호순위 (일표본, 검정값 0) ---------- */
{
  const c = np.wilcoxon[1];
  const expected = [ex('|Z| (부호는 도구 방식: 양의 순위합 기준)', 'Wilcoxon 부호순위 검정 통계량', 0, 'z', Math.abs(c.z), { abs: true }), ex('점근 유의확률(양측)', 'Wilcoxon 부호순위 검정 통계량', 0, 'p', c.pAsym)];
  if (c.pExact !== null) expected.push(ex('정확 유의확률(양측)', 'Wilcoxon 부호순위 검정 통계량', 0, 'pExact', c.pExact));
  entries.push({
    id: 'wilcoxon-1', analysis: 'wilcoxonOne', title: 'Wilcoxon 부호순위(일표본): ' + c.name, description: '차이 값 30개를 검정값 0과 비교한다. 정확 유의확률이 나온다.',
    header: ['차이'], rows: c.diffs.map(d => [d]),
    varRows: [['차이', '차이 값', '연속', '', '']],
    request: { analysis: 'wilcoxonOne', vars: ['차이'], options: { testValue: 0, effect: true } },
    guide: ['분석: Wilcoxon 부호순위 검정 (일표본)', '변수: 차이', '옵션: 검정값 0 (기본), 효과크기 켬'],
    expected: expected
  });
}

/* ---------- Kruskal-Wallis + Dunn ---------- */
{
  const c = np.kw[0];
  const group = [], values = [];
  c.groups.forEach((g, gi) => g.forEach(v => { group.push(gi + 1); values.push(v); }));
  const expected = [ex('H', 'Kruskal-Wallis 검정 통계량', 0, 'h', c.H), ex('유의확률(점근)', 'Kruskal-Wallis 검정 통계량', 0, 'p', c.p)];
  c.dunn.forEach((d, i) => { expected.push(ex('Dunn 보정 유의확률 (집단 ' + (d.i + 1) + ' vs ' + (d.j + 1) + ')', '다중비교: Dunn', i, 'pAdj', d.pAdj)); });
  entries.push({
    id: 'kw-1', analysis: 'kruskalWallis', title: 'Kruskal-Wallis: ' + c.name, description: '3집단 비교와 Dunn 사후검정(Bonferroni 보정).',
    header: ['집단', '값'], rows: group.map((g, i) => [g, values[i]]),
    varRows: [['집단', '집단', '명목', '1=집단1, 2=집단2, 3=집단3', ''], ['값', '측정값', '연속', '', '']],
    request: { analysis: 'kruskalWallis', dependent: ['값'], factor: ['집단'], options: { dunn: true, effect: true } },
    guide: ['분석: Kruskal-Wallis 검정', '종속변수: 값', '집단변수: 집단', '옵션: Dunn 사후검정 켬(기본), 효과크기 켬'],
    expected: expected
  });
}

/* ---------- 신뢰도 (크론바흐 α) ---------- */
const rel = REF('reliability.json');
[[0, '리커트 5문항 (결측 없음)'], [2, '역문항 두 개가 섞여 α가 음수로 나오고 역문항 처리 안내가 보인다'], [1, '결측이 섞인 5문항 (목록별 제거로 사례가 줄어든다)']].forEach(([idx, desc], k) => {
  const c = rel.cases[idx], e = c.expected, kItems = c.columns.length;
  const names = c.columns.map((_, i) => '문항' + (i + 1));
  const expected = [ex('크론바흐 α', '신뢰도 통계량', 0, 'alpha', e.alpha), ex('표준화된 문항에 기초한 α', '신뢰도 통계량', 0, 'alphaStd', e.alphaStd)];
  if (e.ci) { expected.push(ex('α의 95% 신뢰구간 하한', '신뢰도 통계량', 0, 'lo', e.ci.lo), ex('α의 95% 신뢰구간 상한', '신뢰도 통계량', 0, 'hi', e.ci.hi)); }
  expected.push(ex('문항1의 수정된 문항-전체 상관', '문항-전체 통계량', 0, 'corrected', e.items[0].corrected));
  if (e.items[0].alphaIfDeleted !== null) expected.push(ex('문항1을 뺀 α', '문항-전체 통계량', 0, 'alphaDel', e.items[0].alphaIfDeleted));
  entries.push({
    id: 'reliability-' + (k + 1), analysis: 'reliability', title: '신뢰도: ' + c.name, description: desc + ' (사용 사례 ' + e.n + '건)',
    header: names, rows: rowsOf(c.columns),
    varRows: names.map(n => [n, n, '서열', '', '']),
    request: { analysis: 'reliability', vars: names, options: { itemTotal: true, itemDesc: true, ci: true, conf: 95 } },
    guide: ['분석: 신뢰도 분석 (크론바흐 α)', '변수: ' + names.join(', ') + ' (모두 선택)', '옵션: 문항-전체 통계량, 문항 통계량 켬, α의 신뢰구간 켬(95%)'],
    expected: expected
  });
});

/* ---------- 조절효과 ---------- */
const mod = REF('moderation.json');
[[0, '연속형 조절변수, 평균 중심화, 상호작용이 있는 자료'], [3, '범주형 조절변수(3범주, 기준 첫 범주)'], [5, 'Johnson-Neyman 경계가 관측 범위 안에 있다 (효과의 부호가 W에 따라 바뀜)']].forEach(([idx, desc], k) => {
  const c = mod.cases[idx], i = c.input, e = c.expected;
  const cat = i.wtype === 'categorical';
  const covNames = i.covs.map((_, j) => 'C' + (j + 1));
  const header = ['Y', 'X', 'W'].concat(covNames);
  const cols = [i.y, i.x, i.w].concat(i.covs);
  const varRows = [['Y', '종속변수', '연속', '', ''], ['X', '독립변수', '연속', '', ''], ['W', cat ? '조절변수(범주)' : '조절변수', cat ? '명목' : '연속', cat ? Array.from({ length: i.ncat }, (_, q) => q + '=범주' + q).join(', ') : '', '']]
    .concat(covNames.map(n => [n, '공변량 ' + n, '연속', '', '']));
  const expected = [ex('R² 변화량(ΔR²)', '상호작용 검정', 0, 'dr2', e.dr2), ex('F 변화량', '상호작용 검정', 0, 'F', e.Fchg), ex('상호작용 유의확률', '상호작용 검정', 0, 'p', e.pchg)];
  e.slopes.forEach((s, q) => expected.push(ex('조건부 효과 ' + (q + 1) + ' (' + (cat ? '범주 ' + q : ['평균 −1SD', '평균', '평균 +1SD'][q]) + ')', '조건부 효과', q, 'b', s.b)));
  if (e.jn) expected.push(ex('Johnson-Neyman 경계값(참고)', 'Johnson-Neyman', 0, 'range', e.jn.roots.map(v => Number(v.toPrecision(5))).join(', '), { info: true }));
  const options = { center: i.center, conf: i.conf };
  if (cat && i.ref > 0) options.lastRef = true;
  entries.push({
    id: 'moderation-' + (k + 1), analysis: 'moderation', title: '조절효과: ' + c.name, description: desc + ' (사례 ' + e.n + '건)',
    header: header, rows: rowsOf(cols), varRows: varRows,
    request: { analysis: 'moderation', dependent: ['Y'], independent: ['X'], moderator: ['W'], covariates: covNames, options: options },
    guide: ['분석: 조절효과 분석', '종속변수: Y, 독립변수: X, 조절변수: W' + (covNames.length ? ', 공변량: ' + covNames.join(', ') : ''),
      '옵션: 평균 중심화 ' + (i.center ? '켬' : '끔') + ', 신뢰수준 ' + i.conf + (options.lastRef ? ', "범주형 W의 기준 범주를 마지막으로" 켬' : '')],
    expected: expected
  });
});

/* ---------- 매개효과 ---------- */
const med = REF('mediation.json');
[[0, '매개변수 1개, 부분 매개'], [1, '매개변수 2개 + 공변량 1개']].forEach(([idx, desc], k) => {
  const c = med.cases[idx], i = c.input, e = c.expected, m = i.ms.length;
  const mNames = i.ms.map((_, q) => 'M' + (q + 1)), covNames = i.covs.map((_, q) => 'C' + (q + 1));
  const header = ['X'].concat(mNames, ['Y'], covNames);
  const cols = [i.x].concat(i.ms, [i.y], i.covs);
  const T = '효과 분해';
  const expected = [ex('총효과 c', T, 0, 'B', e.cTotal.B), ex('직접효과 c\'', T, 1, 'B', e.cDirect.B)];
  if (m > 1) expected.push(ex('총 간접효과', T, 2, 'B', e.totalIndirect), ex('총 간접효과 부트스트랩 하한', T, 2, 'lo', e.bootTotal.lo), ex('총 간접효과 부트스트랩 상한', T, 2, 'hi', e.bootTotal.hi));
  e.effects.forEach((ef, q) => {
    const row = m > 1 ? 3 + q : 2;
    expected.push(ex('간접효과 ' + (m > 1 ? 'M' + (q + 1) + ' ' : '') + '(a·b)', T, row, 'B', ef.ab), ex('간접효과 부트스트랩 하한', T, row, 'lo', e.bootInd[q].lo), ex('간접효과 부트스트랩 상한', T, row, 'hi', e.bootInd[q].hi));
  });
  e.a.forEach((a, q) => expected.push(ex('a 경로 ' + (m > 1 ? (q + 1) : ''), '경로 계수', q, 'B', a.B)));
  e.b.forEach((b, q) => expected.push(ex('b 경로 ' + (m > 1 ? (q + 1) : ''), '경로 계수', m + q, 'B', b.B)));
  entries.push({
    id: 'mediation-' + (k + 1), analysis: 'mediation', title: '매개효과: ' + c.name, description: desc + ' (사례 ' + e.n + '건, 부트스트랩 ' + i.boot + '회, 시드 ' + i.seed + ' — 시드를 같게 하면 신뢰구간까지 똑같이 나온다)',
    header: header, rows: rowsOf(cols),
    varRows: [['X', '독립변수', '연속', '', '']].concat(mNames.map(n => [n, '매개변수 ' + n, '연속', '', '']), [['Y', '종속변수', '연속', '', '']], covNames.map(n => [n, '공변량 ' + n, '연속', '', ''])),
    request: { analysis: 'mediation', dependent: ['Y'], independent: ['X'], mediators: mNames, covariates: covNames, options: { boot: i.boot, seed: i.seed, conf: i.conf } },
    guide: ['분석: 매개효과 분석', '종속변수: Y, 독립변수: X, 매개변수: ' + mNames.join(', ') + (covNames.length ? ', 공변량: ' + covNames.join(', ') : ''),
      '옵션: 부트스트랩 횟수 ' + i.boot + ', 시드 ' + i.seed + ', 신뢰수준 ' + i.conf],
    expected: expected
  });
});

/* ---------- 로지스틱 회귀 ---------- */
const lg = REF('logistic.json');
/** 기준값 JSON 의 로지스틱 사례 하나를 시험 데이터 항목으로 만든다. recode: 분석 창에서 바로 하는 재부호화(종속변수 원자료와 옵션) */
function logisticEntry(idx, id, desc) {
  const c = lg.cases[idx], i = c.input, e = c.expected;
  const raw = i.raw, depName = i.eventAtLeast !== undefined ? '점수' : i.eventValues ? '등급' : 'y';
  const depVar = i.eventAtLeast !== undefined ? ['점수', '측정 점수', '연속', '', ''] : i.eventValues ? ['등급', '등급(1~4)', '서열', '1=하, 2=중하, 3=중상, 4=상', ''] : ['y', '결과(0/1)', '명목', '0=아니오, 1=예', ''];
  const recode = [];
  if (i.eventValues) recode.push('"사건(1)으로 볼 값"에서 "값을 직접 고르기"를 누르고 ' + i.eventValues.join(', ') + '에 체크 (시트에 0/1 열을 만들 필요 없음)');
  if (i.eventAtLeast !== undefined) recode.push('"사건(1)으로 볼 값"에서 "이 값 이상이면 사건" 아래 칸에 ' + i.eventAtLeast + ' 입력');
  if (i.refCategory !== undefined) recode.push('"범주형 독립변수의 기준 범주"에서 F를 ' + i.refCategory + '으로 선택');
  const options = { cutoff: i.cutoff, conf: i.conf, lastRef: !!i.factor && i.ref > 0 && i.refCategory === undefined };
  if (i.eventValues) options.eventValues = i.eventValues;
  if (i.eventAtLeast !== undefined) options.eventAtLeast = String(i.eventAtLeast);
  if (i.refCategory !== undefined) options.refCategories = { F: i.refCategory };
  const covNames = i.covs.map((_, q) => 'x' + q), hasF = !!i.factor;
  const header = [depName].concat(covNames, hasF ? ['F'] : []);
  const cols = [raw || i.y].concat(i.covs, hasF ? [i.factor] : []);
  const T = '방정식의 변수';
  const names = ['상수'].concat(covNames, hasF ? Array.from({ length: i.ncat }, (_, q) => q).filter(q => q !== i.ref).map(q => 'F = ' + q) : []);
  const expected = [];
  e.B.forEach((b, q) => expected.push(ex('B (' + names[q] + ')', T, q, 'B', b)));
  e.expB.forEach((v, q) => { if (q > 0) expected.push(ex('Exp(B) (' + names[q] + ')', T, q, 'or', v)); });
  expected.push(ex('모형 검정 카이제곱', '모형 검정', 0, 'chi2', e.lr), ex('Nagelkerke R²', '모형 요약', 0, 'nk', e.nk),
    ex('AUC', 'ROC 곡선 아래 면적', 0, 'auc', e.roc.auc), ex('AUC 표준오차(DeLong)', 'ROC 곡선 아래 면적', 0, 'se', e.roc.se),
    ex('Hosmer-Lemeshow 카이제곱', 'Hosmer-Lemeshow 검정', 0, 'chi2', e.hl.chi2), ex('Hosmer-Lemeshow 유의확률', 'Hosmer-Lemeshow 검정', 0, 'p', e.hl.p));
  entries.push({
    id: id, analysis: 'logistic', title: '로지스틱 회귀: ' + c.name, description: desc + ' (사례 ' + e.n + '건, 사건 ' + e.n1 + '건)',
    header: header, rows: rowsOf(cols),
    varRows: [depVar].concat(covNames.map(n => [n, '독립변수 ' + n, '연속', '', '']), hasF ? [['F', '범주 F', '명목', '', '']] : []),
    request: { analysis: 'logistic', dependent: [depName], covariates: covNames, factors: hasF ? ['F'] : [], options: options },
    guide: ['분석: 로지스틱 회귀 (이분형)', '종속변수: ' + depName + ', 독립변수(숫자): ' + covNames.join(', ') + (hasF ? ', 범주형 독립변수: F' : '')]
      .concat(recode, ['옵션: 분류 절단값 ' + i.cutoff + ', 신뢰수준 ' + i.conf + (options.lastRef ? ', "범주형의 기준 범주를 마지막으로" 켬' : '')]),
    expected: expected
  });
}
logisticEntry(0, 'logistic-1', '연속 독립변수 2개');
logisticEntry(1, 'logistic-2', '공변량 1개 + 범주형 3범주(기준 첫 범주)');
logisticEntry(2, 'logistic-3', '사건이 드문 경우(약 10%), 절단값 0.2');
logisticEntry(6, 'logistic-4', '종속변수가 1~4등급인 자료에서 3·4등급을 사건으로 고르기 — 시트에 0/1 열 없이 분석');
logisticEntry(7, 'logistic-5', '연속 점수를 "70 이상이면 사건"으로 나누기 — 시트에 0/1 열 없이 분석');
logisticEntry(8, 'logistic-6', '범주형 3범주에서 가운데 범주(1)를 기준 범주로 지정');

/* ---------- 거리·유사도, 마할라노비스 이상치 ---------- */
const dst = REF('distance.json');
{
  const pair = nm => dst.pairs.find(p => p.name === nm);
  const cellEx = (p, label, tableName, cells) => cells.map(([i, j]) => ex(label + ' (' + (i + 1) + '번째 사례 — ' + (j + 1) + '번째 사례)', tableName, i, 'c' + j, p.expected[i][j]));
  // 1) 단위가 다른 변수는 z점수로 표준화한 뒤 유클리드 거리
  {
    const p = pair('사례 간, z점수, euclid'), names = ['연령', '소득', '점수', '만족도'];
    entries.push({
      id: 'distance-1', analysis: 'distance', title: '거리·유사도: 유클리드 거리, z점수 표준화', description: '단위가 다른 4개 변수(연령, 소득, 점수, 만족도)의 12사례. 표준화하면 소득이 거리를 지배하지 않는다. 표준화 없이 돌리면 단위 경고가 나온다.',
      header: names, rows: p.data,
      varRows: names.map((n, q) => [n, n, '연속', '', '']),
      request: { analysis: 'distance', vars: names, options: { target: 'case', measure: 'euclid', standardize: 'z' } },
      guide: ['분석: 거리·유사도 (마할라노비스 이상치)', '변수: ' + names.join(', '), '옵션: 계산 대상 "사례 간", 측도 "유클리드 거리", 표준화 "z점수"', '"표준화 하지 않음"으로 바꿔 돌리면 "변수의 표준편차가 크게 다릅니다" 경고가 나오는지도 본다.'],
      expected: cellEx(p, '거리', '유클리드 거리 행렬', [[0, 1], [0, 11], [3, 7], [5, 9]])
    });
  }
  // 2) 마할라노비스 이상치 점검 (사례 8번이 이상치: 각 변수는 정상 범위, 조합이 이상)
  {
    const m = dst.maha, names = ['x1', 'x2', 'x3'], T = '마할라노비스 이상치 점검';
    entries.push({
      id: 'distance-2', analysis: 'distance', title: '마할라노비스 이상치 점검: 40사례 3변수', description: '상관이 있는 세 변수 40사례에 한 사례(데이터 탭 ' + (m.top + 2) + '행)만 조합이 이상하다. 변수 하나씩 보면 정상 범위지만 D²로는 이상치 후보가 된다.',
      header: names, rows: m.data,
      varRows: names.map(n => [n, n, '연속', '', '']),
      request: { analysis: 'distance', vars: names, options: { target: 'case', measure: 'euclid', standardize: 'none', outlier: true } },
      guide: ['분석: 거리·유사도 (마할라노비스 이상치)', '변수: ' + names.join(', '), '옵션: "마할라노비스 이상치 점검 (D²)" 켬 (측도는 기본 그대로 두어도 된다)', '결과의 "마할라노비스 이상치 점검"과 "이상치 후보" 표를 기대 결과와 비교한다.'],
      expected: [
        ex('임계값 (χ², p = .001)', T, 2, 'value', m.critical), ex('최대 D²', T, 3, 'value', m.d2[m.top]),
        ex('이상치 후보 수', T, 4, 'value', m.flagged.length), ex('후보 사례의 D²', '이상치 후보', 0, 'd2', m.d2[m.top]),
        ex('후보 사례의 유의확률', '이상치 후보', 0, 'p', m.p[m.top]), ex('참고: 이상치 후보 사례의 데이터 탭 행', '이상치 후보', 0, 'label', '행 ' + (m.top + 2), { info: true })
      ]
    });
  }
  // 3) 이분(0/1) 자료의 Jaccard 유사도
  {
    const p = dst.pairs.find(q => q.binary && q.measure === 'jaccard'), names = p.data[0].map((_, q) => 'v' + (q + 1));
    entries.push({
      id: 'distance-3', analysis: 'distance', title: '거리·유사도: 이분 자료의 Jaccard 유사도', description: '0/1로 코딩한 8개 항목 10사례. 두 사례가 함께 1인 항목의 비율(둘 중 하나라도 1인 항목 대비)을 유사도로 본다.',
      header: names, rows: p.data,
      varRows: names.map(n => [n, n, '명목', '0=아니오, 1=예', '']),
      request: { analysis: 'distance', vars: names, options: { target: 'case', measure: 'jaccard', standardize: 'none' } },
      guide: ['분석: 거리·유사도 (마할라노비스 이상치)', '변수: ' + names.join(', '), '옵션: 계산 대상 "사례 간", 측도 "Jaccard (0/1 자료)"', '값이 클수록 비슷하다 (대각선은 1).'],
      expected: cellEx(p, '유사도', 'Jaccard 계수 (유사도) 행렬', [[0, 1], [2, 5], [4, 9]])
    });
  }
}

/* ---------- 성향점수 매칭 ---------- */
const psmRef = REF('psm.json');
[[0, 'psm-1', '성향점수 최근접 1:1, 캘리퍼 0.2, 비복원(기본 설정)'], [1, 'psm-2', '캘리퍼 안 마할라노비스 매칭 1:1']].forEach(([idx, id, desc]) => {
  const c = psmRef.cases[idx], D = psmRef.data, n = D.y.length;
  const names = ['t', 'x1', 'x2', 'g', 'y'];
  const methodLabel = c.options.method === 'maha' ? '캘리퍼 안 마할라노비스' : '성향점수 최근접';
  const expected = [
    ex('짝지은 처치 수', '매칭 요약', 1, 'n', c.matched.length), ex('버린 처치 수', '매칭 요약', 2, 'n', c.n_treated - c.matched.length),
    ex('사용한 대조 수', '매칭 요약', 4, 'n', c.n_used_controls),
    ex('x1의 SMD (매칭 전)', '균형 진단', 0, 'smd0', c.smd[0].smd0), ex('x1의 SMD (매칭 후)', '균형 진단', 0, 'smd1', c.smd[0].smd1),
    ex('x2의 SMD (매칭 후)', '균형 진단', 1, 'smd1', c.smd[1].smd1), ex('성향점수(logit)의 SMD (매칭 후)', '균형 진단', 4, 'smd1', c.smd[4].smd1),
    ex('단순 비교의 평균 차이 (매칭 전)', '처치 효과 추정', 0, 'diff', c.effect.naive_diff), ex('ATT (매칭 후 평균 차이)', '처치 효과 추정', 1, 'diff', c.effect.att),
    ex('ATT 표준오차', '처치 효과 추정', 1, 'se', c.effect.se), ex('ATT 유의확률', '처치 효과 추정', 1, 'p', c.effect.p)
  ];
  entries.push({
    id: id, analysis: 'psm', title: '성향점수 매칭: ' + c.name, description: desc + ' — 150사례(처치 ' + c.n_treated + ', 대조 ' + c.n_control + '), 처치는 공변량에 따라 선택되고 진짜 처치 효과는 2이다. 매칭 전 단순 비교는 편향돼 있다.',
    header: names, rows: rowsOf([D.treat, D.x1, D.x2, D.g, D.y]),
    varRows: [['t', '처치(프로그램 참여)', '명목', '0=미참여, 1=참여', ''], ['x1', '공변량 1', '연속', '', ''], ['x2', '공변량 2', '연속', '', ''], ['g', '범주 공변량', '명목', '0=A, 1=B, 2=C', ''], ['y', '결과 변수', '연속', '', '']],
    request: { analysis: 'psm', treatment: ['t'], covariates: ['x1', 'x2'], factors: ['g'], outcome: ['y'],
      options: { method: c.options.method, caliper: String(c.options.caliper), ratio: String(c.options.ratio), replace: c.options.replace, conf: '95' } },
    guide: ['분석: 성향점수 매칭 (처치 효과 추정)', '처치 변수: t, 공변량(숫자): x1, x2, 범주형 공변량: g, 결과 변수: y',
      '옵션: 매칭 방식 "' + methodLabel + '", 캘리퍼 ' + c.options.caliper + ', 처치 1명당 대조 수 ' + c.options.ratio + (c.options.replace ? ', 복원 켬' : ''),
      '"균형 진단"에서 매칭 뒤 |SMD|가 0.1 아래로 줄었는지, "처치 효과 추정"에서 매칭 전 단순 비교(편향)와 매칭 후 ATT가 기대 결과와 같은지 본다.'],
    expected: expected
  });
});

/* ---------- k-평균 군집분석 ---------- */
const kmRef = REF('kmeans.json');
[[0, 'kmeans-1', '단위가 다른 세 변수를 z점수로 표준화해 3개 군집으로 나눈다 (군집 수 판단 보조 포함)', true], [2, 'kmeans-2', '표준화하지 않고 군집하면 단위가 큰 소득이 군집을 좌우한다 (경고가 나온다)', false]].forEach(([idx, id, desc, elbow]) => {
  const c = kmRef.cases[idx], names = ['연령', '소득', '점수'];
  const std = c.standardize === 'none' ? '하지 않음' : 'z점수';
  const expected = [
    ex('군집 1의 사례 수', '군집별 사례 수', 0, 'n', c.sizes[0]), ex('군집 2의 사례 수', '군집별 사례 수', 1, 'n', c.sizes[1]), ex('군집 3의 사례 수', '군집별 사례 수', 2, 'n', c.sizes[2]),
    ex('군집 내 제곱합(WSS)', '군집 결과 요약', 5, 'value', c.wss), ex('설명 비율', '군집 결과 요약', 6, 'value', 1 - c.wss / c.tss),
    ex('군집 1 중심: 연령 (표준화한 값 또는 원래 값)', '최종 군집 중심', 0, 'v0', c.centers[0][0]), ex('군집 3 중심: 점수', '최종 군집 중심', 2, 'v2', c.centers[2][2]),
    ex('변수 소득의 F', '변수별 군집 간 차이', 1, 'F', c.f[1].F), ex('군집 1과 2 중심 사이 거리', '군집 중심 사이의 거리', 0, 'c1', c.center_dist[0][1])
  ];
  if (elbow) {
    const e3 = c.elbow.find(e => e.k === 3), e2 = c.elbow.find(e => e.k === 2);
    expected.push(ex('k=2의 군집 내 제곱합', '군집 수 판단 보조', 0, 'wss', e2.wss), ex('k=3의 평균 실루엣', '군집 수 판단 보조', 1, 'sil', e3.silhouette));
  }
  entries.push({
    id: id, analysis: 'kmeans', title: 'k-평균 군집분석: ' + c.name, description: desc + ' — 90사례, 세 변수(연령, 소득, 점수)에 실제로 3개의 집단이 섞여 있다.',
    header: names, rows: kmRef.data,
    varRows: [['연령', '연령', '연속', '', ''], ['소득', '월소득(만원)', '연속', '', ''], ['점수', '점수', '연속', '', '']],
    request: { analysis: 'kmeans', vars: names, options: { k: String(c.K), seed: String(c.seed), starts: String(c.starts), standardize: c.standardize, elbow: elbow } },
    guide: ['분석: k-평균 군집분석', '변수: ' + names.join(', '), '옵션: 군집 수 ' + c.K + ', 표준화 "' + std + '"' + (elbow ? ', 군집 수 판단 보조 켬' : ', 군집 수 판단 보조 끔') + ', 세부 옵션의 시드 ' + c.seed + ', 시작 횟수 ' + c.starts,
      '시드가 같으면 항상 같은 결과이다. 군집 번호는 큰 군집부터 매겨진다. 기대 결과와 같은지 본다.'],
    expected: expected
  });
});

/* ---------- 의사결정나무 (CART) ---------- */
const cartRef = REF('cart.json');
[[0, 'A', 'cart-1', '이진 분류 나무: 결과(0/1)를 x1, x2, 잡음 변수, 범주 c로 예측. 가지치기로 작은 나무가 남고 잡음 변수는 쓰이지 않는다', ['x1', 'x2', 'noise'], '결과', '0=아니오, 1=예'],
  [3, 'B', 'cart-2', '3분류 나무: 범주 c(4개)의 부분집합 분할과 숫자 변수 분할을 함께 탐색한다', ['x1', 'x2'], '유형', '0=가형, 1=나형, 2=다형'],
  [4, 'C', 'cart-3', '회귀 나무: 연속 결과를 예측한다 (나무가 비교적 크게 남고 RMSE와 R²를 본다)', ['x1', 'x2'], '점수', '']].forEach(([ci, key, id, desc, numNames, yLabel, yVl]) => {
  const c = cartRef.cases[ci], d = cartRef.data[key], n = d.y.length;
  const cls = c.mode === 'class';
  const header = ['y'].concat(numNames, ['c']);
  const cols = [d.y].concat(numNames.map(nm => d[nm]), [d.c]);
  const sel = c.selected, T = '가지치기 곡선', trainErr0 = c.train_err[0] / n;
  const expected = [
    ex('선택한 나무의 잎 수', '나무 요약', 3, 'value', c.leaves[sel]), ex('키운 나무의 잎 수', '나무 요약', 2, 'value', c.leaves[0]),
    ex('키운 나무의 훈련 ' + (cls ? '오분류율' : 'MSE'), T, 0, 'train', trainErr0),
    ex('선택한 나무의 CV ' + (cls ? '오분류율' : 'MSE'), T, sel, 'cv', c.xerror[sel]), ex('선택한 나무의 CV 표준오차', T, sel, 'se', c.se[sel]),
    ex('가장 큰 복잡도 α (루트만 남는 지점)', T, c.alphas.length, 'alpha', c.alphas[c.alphas.length - 1])
  ];
  if (cls) expected.push(ex('재대입 전체 정확도(%)', '분류표', C_ROWS(c), 'pct', c.accuracy * 100));
  else expected.push(ex('R² (재대입)', '모형 성능', 0, 'value', c.r2), ex('RMSE (재대입)', '모형 성능', 1, 'value', c.rmse), ex('RMSE (교차검증)', '모형 성능', 2, 'value', Math.sqrt(c.xerror[sel])));
  const top = Object.keys(c.importance).sort((a, b) => c.importance[b] - c.importance[a])[0];
  expected.push(ex('가장 중요한 변수의 불순도 감소 합', '변수 중요도', 0, 'gain', c.importance[top]), ex('참고: 가장 중요한 변수', '변수 중요도', 0, 'name', top === 'c' ? '구분' : top, { info: true }));
  const covNames = numNames;
  const opts = { maxDepth: String(c.params.maxDepth), prune: c.rule, minSplit: String(c.params.minSplit), minLeaf: String(c.params.minLeaf), folds: String(c.nfold), seed: String(c.seed) };
  const varRows = [[ 'y', yLabel, cls ? '명목' : '연속', yVl, '' ]].concat(covNames.map(nm => [nm, '예측변수 ' + nm, '연속', '', '']), [['c', '구분', '명목', c.name.indexOf('4개') >= 0 ? '0=a, 1=b, 2=c, 3=d' : '0=A, 1=B, 2=C', '']]);
  entries.push({
    id: id, analysis: 'cart', title: '의사결정나무: ' + c.name, description: desc + ' — ' + n + '사례.',
    header: header, rows: rowsOf(cols), varRows: varRows,
    request: { analysis: 'cart', dependent: ['y'], covariates: covNames, factors: ['c'], options: opts },
    guide: ['분석: 의사결정나무 (CART)', '종속변수: y, 예측변수(숫자): ' + covNames.join(', ') + ', 범주형 예측변수: c',
      '옵션: 최대 깊이 ' + c.params.maxDepth + ', 가지치기 "' + (c.rule === 'cv1se' ? '교차검증 1 표준오차 규칙' : c.rule === 'cv' ? '교차검증 최소 오차' : '하지 않음') + '"',
      '세부 옵션: 분할 최소 사례 수 ' + c.params.minSplit + ', 잎 최소 사례 수 ' + c.params.minLeaf + ', 폴드 수 ' + c.nfold + ', 시드 ' + c.seed,
      '결과의 가지치기 곡선(어느 크기의 나무가 선택됐는지), 나무 구조 표와 나무 그림, 변수 중요도, ' + (cls ? '분류표' : '모형 성능') + '이 기대 결과와 같은지 본다. 시드가 같으면 항상 같은 결과이다.'],
    expected: expected
  });
});
function C_ROWS(c) { return c.C; }

/* ---------- McNemar 검정 ---------- */
const mcRef = REF('mcnemar.json');
[[0, 'mcnemar-1', '프로그램 전후 취업 여부(0=미취업, 1=취업) 변화. 값이 바뀐 쌍이 21쌍이라 정확 유의확률과 점근 유의확률이 비슷하다'],
  [1, 'mcnemar-2', '값이 바뀐 쌍이 7쌍뿐인 작은 표본 — 정확 유의확률을 보고 검정력이 낮다는 경고가 나온다']].forEach(([idx, id, desc]) => {
  const c = mcRef.cases[idx];
  const T = 'McNemar 검정 통계량', V = '변화 요약';
  entries.push({
    id: id, analysis: 'mcnemar', title: 'McNemar 검정: ' + c.name, description: desc + ' — ' + c.n + '명.',
    header: ['사전취업', '사후취업'], rows: c.before.map((v, i) => [v, c.after[i]]),
    varRows: [['사전취업', '프로그램 전 취업 여부', '명목', '0=미취업, 1=취업', ''], ['사후취업', '프로그램 후 취업 여부', '명목', '0=미취업, 1=취업', '']],
    request: { analysis: 'mcnemar', first: ['사전취업'], second: ['사후취업'], options: {} },
    guide: ['분석: McNemar 검정 (비모수, 전·후 예/아니오)', '변수1(전): 사전취업, 변수2(후): 사후취업', '옵션 없음',
      '교차표의 칸 수(미취업→취업 ' + c.b + '명, 취업→미취업 ' + c.c + '명)와 정확 유의확률, 카이제곱, 비율 차이, 오즈비가 기대 결과와 같은지 본다.'],
    expected: [
      ex('미취업→취업 (b)', '교차표', 0, 'c1', c.b), ex('취업→미취업 (c)', '교차표', 1, 'c0', c.c), ex('값이 바뀐 쌍', T, 0, 'nd', c.nd),
      ex('정확 유의확률(양측)', T, 0, 'exactP', c.exact), ex('카이제곱(연속성 수정)', T, 0, 'chi2', c.chi2), ex('점근 유의확률', T, 0, 'chiP', c.chi_p),
      ex('비율 차이 (후 − 전)', V, 2, 'value', c.diff), ex('비율 차이 하한', V, 2, 'lo', c.diff_lo), ex('오즈비', V, 3, 'value', c.odds), ex('오즈비 하한', V, 3, 'lo', c.odds_lo)
    ]
  });
});

/* ---------- 효과크기 (Cohen d, Hedges g): t검정 3종 ---------- */
const efRef = REF('effect.json');
{
  const E = '효과크기';
  const c1 = efRef.cases.find(c => c.kind === 'independent' && c.name === '보통 효과');
  entries.push({
    id: 'effect-1', analysis: 'ttestInd', title: '독립표본 t검정의 효과크기 (Cohen d, Hedges g)', description: '참여군 ' + c1.a.length + '명과 미참여군 ' + c1.b.length + '명의 점수. t검정 결과 아래의 효과크기 표(d, g와 정확 신뢰구간)를 scipy 값과 비교한다.',
    header: ['점수', '집단'], rows: c1.a.map(v => [v, 1]).concat(c1.b.map(v => [v, 2])),
    varRows: [['점수', '점수', '연속', '', ''], ['집단', '집단', '명목', '1=참여, 2=미참여', '']],
    request: { analysis: 'ttestInd', dependent: ['점수'], factor: ['집단'], options: { effect: true, conf: '95' } },
    guide: ['분석: 독립표본 t검정', '종속변수: 점수, 집단변수: 집단', '옵션: "효과크기" 켬 (기본값), 신뢰수준 95 (세부 옵션)', '결과의 "효과크기" 표에서 Cohen d, Hedges g와 95% 신뢰구간이 기대 결과와 같은지 본다. 방향은 집단1(참여) − 집단2(미참여)이다.'],
    expected: [ex('Cohen d', E, 0, 'est', c1.d), ex('d 신뢰구간 하한', E, 0, 'lo', c1.d_lo), ex('d 신뢰구간 상한', E, 0, 'hi', c1.d_hi), ex('Hedges g', E, 1, 'est', c1.g), ex('g 신뢰구간 하한', E, 1, 'lo', c1.g_lo), ex('표준화 기준값(합동 표준편차)', E, 0, 'std', c1.sp)]
  });
  const c2 = efRef.cases.find(c => c.kind === 'paired' && c.name === '사전·사후 향상');
  entries.push({
    id: 'effect-2', analysis: 'ttestPaired', title: '대응표본 t검정의 효과크기 (d_z, g_z, d_av)', description: c2.x.length + '명의 사전·사후 점수. 효과크기 표의 d_z(대응표본의 Cohen d), Hedges g_z, d_av를 scipy 값과 비교한다.',
    header: ['사전', '사후'], rows: c2.x.map((v, i) => [v, c2.y[i]]),
    varRows: [['사전', '사전 점수', '연속', '', ''], ['사후', '사후 점수', '연속', '', '']],
    request: { analysis: 'ttestPaired', first: ['사전'], second: ['사후'], options: { effect: true, conf: '95' } },
    guide: ['분석: 대응표본 t검정', '변수1: 사전, 변수2: 사후', '옵션: "효과크기" 켬 (기본값)', '효과크기 표의 세 줄(d_z, g_z, d_av)이 기대 결과와 같은지 본다. d_av는 신뢰구간이 없다. 방향은 변수1 − 변수2 (사전 − 사후)라 향상이면 음수이다.'],
    expected: [ex("Cohen's d_z", E, 0, 'est', c2.d), ex('d_z 신뢰구간 하한', E, 0, 'lo', c2.d_lo), ex('d_z 신뢰구간 상한', E, 0, 'hi', c2.d_hi), ex("Hedges' g_z", E, 1, 'est', c2.g), ex('d_av', E, 2, 'est', c2.d_av), ex('표준화 기준값(차이의 표준편차)', E, 0, 'std', c2.sd_diff)]
  });
  const c3 = efRef.cases.find(c => c.kind === 'one' && c.name === '기준값보다 높음');
  entries.push({
    id: 'effect-3', analysis: 'ttest1', title: '일표본 t검정의 효과크기 (Cohen d)', description: c3.values.length + '명의 점수를 기준값 ' + c3.mu + '점과 비교한다. d와 정확 신뢰구간을 scipy 값과 비교한다.',
    header: ['점수'], rows: c3.values.map(v => [v]),
    varRows: [['점수', '점수', '연속', '', '']],
    request: { analysis: 'ttest1', vars: ['점수'], options: { testValue: String(c3.mu), effect: true, assumptions: true, conf: '95' } },
    guide: ['분석: 일표본 t검정', '변수: 점수', '옵션: 검정값 ' + c3.mu + ', "효과크기" 켬 (기본값)', '효과크기 표의 d와 g, 95% 신뢰구간이 기대 결과와 같은지 본다 (d = (평균 − 검정값) / 표본 표준편차).'],
    expected: [ex('Cohen d', E, 0, 'est', c3.d), ex('d 신뢰구간 하한', E, 0, 'lo', c3.d_lo), ex('d 신뢰구간 상한', E, 0, 'hi', c3.d_hi), ex('Hedges g', E, 1, 'est', c3.g)]
  });
}

/* ---------- Friedman 검정 (3회 이상 반복측정) ---------- */
const frRef = REF('friedman.json');
[[1, 'friedman-1', '5점 척도 만족도를 사전·중간·사후 3번 측정한 30명 (동순위가 많아 동순위 보정이 걸린다)', ['사전', '중간', '사후']],
  [0, 'friedman-2', '점수를 4번 측정한 20명 (동순위 없음, 시간이 갈수록 오르는 경향)', ['1회', '2회', '3회', '4회']]].forEach(([idx, id, desc, names]) => {
  const c = frRef.cases[idx], T = 'Friedman 검정 통계량', P = '사후비교';
  const expected = [ex('Friedman χ²', T, 0, 'chi2', c.chi2), ex('점근 유의확률', T, 0, 'p', c.p), ex("Kendall's W", T, 0, 'w', c.w)];
  names.forEach((nm, j) => expected.push(ex('평균 순위 (' + nm + ')', '순위', j, 'meanRank', c.mean_ranks[j])));
  [0, c.pairs.length - 1].forEach(q => expected.push(ex('사후비교 ' + names[c.pairs[q].i] + ' vs ' + names[c.pairs[q].j] + ' Z', P, q, 'z', c.pairs[q].z), ex('사후비교 ' + names[c.pairs[q].i] + ' vs ' + names[c.pairs[q].j] + ' 보정 유의확률', P, q, 'pAdj', c.pairs[q].p_adj)));
  entries.push({
    id: id, analysis: 'friedman', title: 'Friedman 검정: ' + c.name, description: desc + ' — ' + c.n + '명.',
    header: names, rows: c.data,
    varRows: names.map(nm => [nm, nm + ' 측정', idx === 1 ? '서열' : '연속', '', '']),
    request: { analysis: 'friedman', vars: names, options: { posthoc: true, effect: true } },
    guide: ['분석: Friedman 검정 (비모수, 3회 이상 반복측정)', '변수: ' + names.join(', ') + ' (이 순서대로 보내세요)', '옵션: 사후비교 Dunn–Bonferroni, 효과크기 Kendall W (모두 기본값)',
      '평균 순위, χ²와 유의확률, W, 사후비교의 Z와 보정 유의확률이 기대 결과와 같은지 본다.' + (idx === 1 ? ' 동순위가 많은 자료라 χ²에 동순위 보정이 적용된다.' : '')],
    expected: expected
  });
});

/* ---------- 일원 반복측정 분산분석 ---------- */
const rmRef = REF('rmanova.json');
[[1, 'rmanova-1', '프로그램 전·사후 1·3·6개월 점수를 측정한 25명. 시간이 갈수록 오차가 커져 **구형성이 깨진다** (Greenhouse–Geisser 행을 읽는 경우)', ['사전', '1개월', '3개월', '6개월']],
  [0, 'rmanova-2', '같은 형식의 20명, 구형성이 대체로 성립해 "구형성 가정" 행을 읽는 경우', ['사전', '1개월', '3개월', '6개월']]].forEach(([idx, id, desc, names]) => {
  const c = rmRef.cases[idx], M = 'Mauchly', T = '개체-내 효과 검정', P = '쌍별 비교';
  const last = c.pairs.length - 1;
  entries.push({
    id: id, analysis: 'rmanova', title: '일원 반복측정 분산분석: ' + c.name, description: desc + ' — ' + c.n + '명.',
    header: names, rows: c.data,
    varRows: names.map(nm => [nm, nm + ' 점수', '연속', '', '']),
    request: { analysis: 'rmanova', vars: names, options: { posthoc: true, plot: true, conf: '95' } },
    guide: ['분석: 일원 반복측정 분산분석 (3회 이상 측정)', '변수: ' + names.join(', ') + ' (이 순서대로 보내세요)', '옵션: 쌍별 비교와 평균 선 그림 켬 (기본값)',
      'Mauchly 구형성 검정(W, 유의확률, ε), 개체-내 효과 검정의 네 행(교정에 따라 자유도·유의확률이 다르고 읽을 행이 굵게 표시), 부분 η², 쌍별 비교가 기대 결과와 같은지 본다.'],
    expected: [
      ex("Mauchly's W", M, 0, 'W', c.mauchly.W), ex('Mauchly 유의확률', M, 0, 'p', c.mauchly.p), ex('ε Greenhouse–Geisser', M, 0, 'gg', c.corr.gg.eps), ex('ε Huynh–Feldt', M, 0, 'hf', c.corr.hf.eps),
      ex('F (구형성 가정)', T, 0, 'F', c.F), ex('유의확률 (구형성 가정)', T, 0, 'p', c.p), ex('유의확률 (Greenhouse–Geisser)', T, 1, 'p', c.corr.gg.p), ex('유의확률 (Huynh–Feldt)', T, 2, 'p', c.corr.hf.p),
      ex('부분 η²', T, 0, 'peta', c.peta2), ex('일반화 η²', T, 0, 'geta', c.geta2),
      ex('쌍별 비교 ' + names[c.pairs[0].i] + '−' + names[c.pairs[0].j] + ' 평균차이', P, 0, 'diff', c.pairs[0].diff), ex('쌍별 비교 ' + names[c.pairs[0].i] + '−' + names[c.pairs[0].j] + ' 보정 유의확률', P, 0, 'pAdj', c.pairs[0].p_bonf),
      ex('쌍별 비교 ' + names[c.pairs[last].i] + '−' + names[c.pairs[last].j] + ' 보정 유의확률', P, last, 'pAdj', c.pairs[last].p_bonf)
    ]
  });
});

/* ---------- 위계적 회귀 ---------- */
const hrRef = REF('hierarchical.json');
[[0, 'hierarchical-1', '통제변수 2개를 1단계에, 핵심 변수를 2단계에, 나머지 2개를 3단계에 넣는다. 3단계에서는 R² 증가가 유의하지 않다', 'scale'],
  [1, 'hierarchical-2', '1단계에 3범주 변수(명목, 더미 자동 변환)와 숫자 변수를 넣고 2단계에 숫자 2개를 추가한다', 'scale'],
  [3, 'hierarchical-3', '리커트(1~5) 점수 변수를 서열 척도로 두고 숫자 점수로 투입한다. 신뢰수준 99%', 'ordinal']].forEach(([idx, id, desc, numScale]) => {
  const c = hrRef.cases[idx];
  const vars = [];
  c.blocks.forEach(b => b.forEach(v => vars.push(v)));
  const header = ['y'].concat(vars.map(v => v.name));
  const nStep = c.steps.length, S = '모형 요약', Cf = '계수';
  const expected = [];
  c.steps.forEach((s, k) => {
    expected.push(ex('R² (' + (k + 1) + '단계)', S, k, 'r2', s.r2), ex('R² 변화량 (' + (k + 1) + '단계)', S, k, 'dR2', s.dR2),
      ex('F 변화량 (' + (k + 1) + '단계)', S, k, 'dF', s.dF), ex('유의확률 F 변화량 (' + (k + 1) + '단계)', S, k, 'pF', s.pDF));
  });
  let off = 0;
  c.steps.slice(0, nStep - 1).forEach(s => { off += s.names.length; });
  const fin = c.steps[nStep - 1];
  fin.names.forEach((nm, j) => {
    expected.push(ex('마지막 단계 B (' + nm + ')', Cf, off + j, 'B', fin.B[j]));
    if (j > 0) expected.push(ex('마지막 단계 베타 (' + nm + ')', Cf, off + j, 'beta', fin.beta[j]));
  });
  const request = { analysis: 'hierarchical', dependent: ['y'], options: { ci: true, conf: String(c.conf) } };
  c.blocks.forEach((b, bi) => { request['block' + (bi + 1)] = b.map(v => v.name); });
  entries.push({
    id: id, analysis: 'hierarchical', title: '위계적 회귀: ' + c.name, description: desc + ' — ' + c.n + '명.',
    header: header, rows: rowsOf([c.y].concat(vars.map(v => v.values))),
    varRows: [['y', '결과 점수', '연속', '', '']].concat(vars.map(v => [v.name, '독립변수 ' + v.name, v.kind === 'factor' ? '명목' : (numScale === 'ordinal' ? '서열' : '연속'), '', ''])),
    request: request,
    guide: ['분석: 위계적 회귀 (단계별 변수 추가)', '종속변수: y'].concat(c.blocks.map((b, bi) => (bi + 1) + '단계 독립변수: ' + b.map(v => v.name).join(', ')),
      ['옵션: B의 신뢰구간 켬' + (c.conf !== 95 ? ', 신뢰수준 ' + c.conf : ''),
        '모형 요약의 단계별 R², R 제곱 변화량, F 변화량, 유의확률 F 변화량과 계수 표의 마지막 단계 B·베타가 기대 결과와 같은지 본다.']),
    expected: expected
  });
});

/* ---------- 탐색적 요인분석 ---------- */
const efaRef = REF('efa.json');
const EfaMod = require('../src/Stats_Efa.js');
[[0, 'efa-1', '문항 9개가 서로 직교인 3요인 구조. 주축분해와 베리맥스로 요인 수를 고유값 1 이상으로 정하면 3개가 나오고 각 요인에 문항 3개씩 묶인다', 'scale'],
  [1, 'efa-2', '상관이 있는 2요인(r=.4)에 교차적재 문항이 있는 8문항. 요인 수를 2로 직접 지정한다', 'scale'],
  [2, 'efa-3', '5점 척도(정수) 12문항, 서열 척도. 주성분과 베리맥스(성분 3개)', 'ordinal']].forEach(([idx, id, desc, scale]) => {
  const c = efaRef.cases[idx];
  const names = c.data[0].map((_, j) => 'q' + (j + 1));
  const word = c.method === 'pca' ? '성분' : '요인';
  const tol = c.method === 'paf' ? { tol: 3e-5 } : {};
  const expected = [ex('KMO', 'KMO와 Bartlett', 0, 'kmo', c.kmo), ex('Bartlett 카이제곱', 'KMO와 Bartlett', 0, 'chi2', c.bartlett.chi2), ex('Bartlett 유의확률', 'KMO와 Bartlett', 0, 'p', c.bartlett.p, { tol: 1e-4 })];
  [0, c.p - 1].forEach(k => expected.push(ex('변수별 MSA (q' + (k + 1) + ')', '공통성', k, 'msa', c.msa[k]), ex('추출 공통성 (q' + (k + 1) + ')', '공통성', k, 'ext', c.h2[k], tol)));
  for (let k = 0; k < c.m; k++) expected.push(ex('초기 고유값 합계 (' + word + ' ' + (k + 1) + ')', '설명된 총분산', k, 'iTot', c.eig[k]), ex('추출 제곱합 적재값 합계 (' + word + ' ' + (k + 1) + ')', '설명된 총분산', k, 'eTot', c.ss_ext[k], tol));
  const fin = c.rotated || c.L;
  if (c.rotated) for (let k = 0; k < c.m; k++) expected.push(ex('회전 제곱합 적재값 합계 (' + word + ' ' + (k + 1) + ')', '설명된 총분산', k, 'rTot', c.ss_rot[k], tol));
  const order = EfaMod.efaSortOrder(fin);
  [0, Math.floor(c.p / 2), c.p - 1].forEach(r => {
    const k = order[r], row = fin[k];
    let best = 0;
    row.forEach((v, f) => { if (Math.abs(v) > Math.abs(row[best])) best = f; });
    expected.push(ex('회전된 적재값 (q' + (k + 1) + ', ' + word + ' ' + (best + 1) + ')', c.rotated ? '회전된 ' : (word + ' 행렬'), r, 'f' + best, row[best], { tol: c.method === 'paf' ? 1e-4 : 1e-5 }));
  });
  const options = { method: c.method, rotation: c.rotation, nfactors: c.nfactors ? String(c.nfactors) : '', scree: true, sort: true, suppress: '0.3' };
  entries.push({
    id: id, analysis: 'efa', title: '탐색적 요인분석: ' + c.name, description: desc + ' — ' + c.n + '명, ' + c.p + '문항.',
    header: names, rows: c.data,
    varRows: names.map((nm, j) => [nm, '문항 ' + (j + 1), scale === 'ordinal' ? '서열' : '연속', '', '']),
    request: { analysis: 'efa', vars: names, options: options },
    guide: ['분석: 탐색적 요인분석 (EFA)', '문항: ' + names.join(', '),
      '옵션: 추출 방법 ' + (c.method === 'pca' ? '주성분' : '주축분해') + ', 회전 ' + (c.rotation === 'none' ? '회전 안 함' : '베리맥스') + ', 요인 수 ' + (c.nfactors ? c.nfactors + '(직접 입력)' : '비움(고유값 1 이상 → ' + c.m + '개)') + ', 스크리·정렬 켬',
      'KMO·Bartlett, 공통성(MSA), 설명된 총분산, ' + (c.rotated ? '회전 전·후 ' : '') + word + ' 행렬이 기대 결과와 같은지 본다. 숨긴 적재값(0.3 미만)은 빈칸이다.'],
    expected: expected
  });
});

/* ---------- 탐색적 요인분석 2단계: 프로맥스·병렬분석·요인 점수 ---------- */
{
  const c = efaRef.cases[1];                                  // 상관 2요인 + 교차적재 (프로맥스에 알맞다)
  const names = c.data[0].map((_, j) => 'q' + (j + 1));
  const pm = c.promax;
  const T = '패턴 행렬';
  const expected = [ex('요인 상관 (요인 1, 요인 2)', '요인 상관 행렬', 0, 'f1', pm.phi[0][1], { tol: 1e-4 })];
  const order = EfaMod.efaSortOrder(pm.P);
  order.forEach((k, r) => {
    const row = pm.P[k];
    let best = 0;
    row.forEach((v, f) => { if (Math.abs(v) > Math.abs(row[best])) best = f; });
    if (r === 0 || r === 3 || r === 7) expected.push(ex('패턴 행렬 (q' + (k + 1) + ', 요인 ' + (best + 1) + ')', T, r, 'f' + best, row[best], { tol: 1e-3 }));
  });
  const kSort = order;                                       // 구조 행렬도 패턴 행렬의 정렬 순서를 따른다
  [0, 5].forEach(r => {
    const k = kSort[r], row = pm.S[k];
    let best = 0;
    row.forEach((v, f) => { if (Math.abs(v) > Math.abs(row[best])) best = f; });
    expected.push(ex('구조 행렬 (q' + (k + 1) + ', 요인 ' + (best + 1) + ')', '구조 행렬', r, 'f' + best, row[best], { tol: 1e-3 }));
  });
  entries.push({
    id: 'efa-4', analysis: 'efa', title: '탐색적 요인분석: 프로맥스(사각 회전), 병렬분석, 요인 점수 저장',
    description: '요인 사이 상관이 .4인 8문항(교차적재 포함). 프로맥스로 회전해 패턴 행렬·구조 행렬·요인 상관을 보고 병렬분석이 권하는 요인 수를 확인한다. 요인 점수를 데이터 탭에 저장하는 옵션도 시험한다 — ' + c.n + '명, ' + c.p + '문항.',
    header: names, rows: c.data,
    varRows: names.map((nm, j) => [nm, '문항 ' + (j + 1), '연속', '', '']),
    request: { analysis: 'efa', vars: names, options: { method: 'paf', rotation: 'promax', nfactors: '2', scree: true, parallel: true, sort: true, save: false, seed: '2026', suppress: '0.3' } },
    guide: ['분석: 탐색적 요인분석 (EFA)', '문항: ' + names.join(', '),
      '옵션: 추출 방법 주축분해, 회전 프로맥스, 요인 수 2(직접 입력), 스크리 그림·병렬분석·크기순 정렬 켬',
      '요인 상관 행렬의 요인 1–요인 2 상관, 패턴 행렬과 구조 행렬(회전 전 요인 행렬은 그대로), 병렬분석 표(권하는 요인 수 2개)가 기대 결과와 같은지 본다. 프로맥스에서는 회전 뒤 제곱합 열이 없다.',
      '추가 시험: "요인 점수를 데이터 탭에 저장"을 켜고 다시 실행하면 데이터 탭 맨 오른쪽에 요인점수1, 요인점수2 열이 생기고 변수정의에도 등록된다.'],
    expected: expected
  });
}

/* ---------- 이원 분산분석 ---------- */
const a2Ref = REF('anova2.json');
[[0, 'anova2-1', '셀마다 10건인 균형 2×3 설계, 상호작용이 있다. 유형 III와 II의 값이 같다', 'III', '1=낮음, 2=높음', '1=가, 2=나, 3=다'],
  [1, 'anova2-2', '셀마다 사례 수가 다른 불균형 3×3 설계(4~12건). 유형 III 기본, 추정 주변평균과 사후검정(Tukey)도 본다', 'III', '1=낮음, 2=중간, 3=높음', '1=가, 2=나, 3=다'],
  [3, 'anova2-3', '작은 셀(3~9건)의 불균형 2×2 설계. 제곱합 유형을 II로 두면 주효과 값이 III와 다르다', 'II', '1=남, 2=여', '1=통제, 2=프로그램']].forEach(([idx, id, desc, ssType, vlA, vlB]) => {
  const c = a2Ref.cases[idx];
  const tm = c[ssType].terms, T = '개체-간 효과 검정', iA = ssType === 'III' ? 2 : 1;
  const expected = [
    ex('요인 A F', T, iA, 'F', tm.A.F), ex('요인 A 유의확률', T, iA, 'p', tm.A.p, { tol: 1e-5 }), ex('요인 A 부분 η²', T, iA, 'peta', tm.A.peta),
    ex('요인 B F', T, iA + 1, 'F', tm.B.F), ex('요인 B 유의확률', T, iA + 1, 'p', tm.B.p, { tol: 1e-5 }),
    ex('상호작용 F', T, iA + 2, 'F', tm.AB.F), ex('상호작용 유의확률', T, iA + 2, 'p', tm.AB.p, { tol: 1e-5 }), ex('상호작용 부분 η²', T, iA + 2, 'peta', tm.AB.peta),
    ex('오차 제곱합', T, iA + 3, 'ss', c.sse), ex('Levene 통계량', 'Levene', 0, 'W', c.levene.W), ex('Levene 유의확률', 'Levene', 0, 'p', c.levene.p, { tol: 1e-5 }),
    ex('첫 셀 평균', '기술통계', 0, 'mean', c.cells[0][0].mean), ex('단순 주효과 (요인 B 첫 수준에서 A) F', '단순 주효과', 0, 'F', c.simple.aAtB[0].F),
    ex('단순 주효과 (요인 A 첫 수준에서 B) 유의확률', '단순 주효과', c.b, 'p', c.simple.bAtA[0].p, { tol: 1e-5 })
  ];
  if (ssType === 'III') expected.push(ex('절편 제곱합', T, 1, 'ss', c.III.intercept.ss), ex('수정된 모형 F', T, 0, 'F', c.III.model.F));
  if (c.posthoc.A) expected.push(ex('요인 A 사후검정 첫 쌍 평균차이', '사후검정: 요인 A', 0, 'diff', c.posthoc.A[0].diff), ex('요인 A 사후검정 첫 쌍 Tukey 유의확률', '사후검정: 요인 A', 0, 'p', c.posthoc.A[0].p_tukey, { tol: 1e-4 }));
  if (id === 'anova2-2') expected.push(ex('추정 주변평균 (요인 A 첫 수준)', '추정 주변평균', 0, 'mean', c.emm.rows[0].mean), ex('추정 주변평균 표준오차 (요인 A 첫 수준)', '추정 주변평균', 0, 'se', c.emm.rows[0].se));
  entries.push({
    id: id, analysis: 'anova2', title: '이원 분산분석: ' + c.name, description: desc + ' — ' + c.n + '건.',
    header: ['y', 'A', 'B'], rows: c.y.map((v, i) => [v, c.A[i], c.B[i]]),
    varRows: [['y', '결과 점수', '연속', '', ''], ['A', '요인 A', '명목', vlA, ''], ['B', '요인 B', '명목', vlB, '']],
    request: { analysis: 'anova2', dependent: ['y'], factor1: ['A'], factor2: ['B'], options: { ss: ssType, levene: true, plot: true, simple: true, posthoc: 'tukey', emmeans: id === 'anova2-2', conf: '95' } },
    guide: ['분석: 이원 분산분석 (요인 2개)', '종속변수: y, 요인 1: A, 요인 2: B', '옵션: 제곱합 유형 ' + ssType + ', Levene·상호작용 그림·단순 주효과 켬, 사후검정 Tukey' + (id === 'anova2-2' ? ', 추정 주변평균 켬' : ''),
      '개체-간 효과 검정의 F·유의확률·부분 η², 오차 제곱합, Levene, 단순 주효과가 기대 결과와 같은지 본다. 상호작용 그림의 x축에 요인 1의 수준 이름이 나오는지도 본다.'],
    expected: expected
  });
});

/* ---------- 혼합 설계 분산분석 (집단 × 반복측정) ---------- */
const mxRef = REF('mixed.json');
[[0, 'mixed-1', '집단 2개 × 3회 측정, 집단당 12명으로 균형이고 구형성이 성립한다. 시간×집단 상호작용이 유의하다', ['사전', '중간', '사후'], '1=통제, 2=프로그램'],
  [1, 'mixed-2', '집단 3개 × 4회 측정, 집단 크기가 8·12·10명으로 다르고 측정이 늘수록 오차가 커져 **구형성이 깨진다** (Greenhouse–Geisser 행)', ['사전', '1개월', '3개월', '6개월'], '1=통제, 2=프로그램A, 3=프로그램B'],
  [2, 'mixed-3', '집단 2개 × 2회 측정(사전·사후), 집단당 15명. 측정이 2회라 구형성 검정이 없고 구형성 가정 행만 나온다', ['사전', '사후'], '1=통제, 2=프로그램']].forEach(([idx, id, desc, names, vl]) => {
  const c = mxRef.cases[idx];
  const T = '개체-내 효과 검정 (측정)', X = '개체-내 효과 검정 (측정 × 집단)', B = '개체-간 효과 검정';
  const expected = [
    ex('집단 F (개체-간)', B, 1, 'F', c.between.group.F), ex('집단 유의확률 (개체-간)', B, 1, 'p', c.between.group.p, { tol: 1e-5 }), ex('집단 부분 η²', B, 1, 'peta', c.between.group.peta),
    ex('개체-간 오차 제곱합', B, 2, 'ss', c.between.err.ss),
    ex('측정 F', T, 0, 'F', c.time.F), ex('측정 유의확률 (구형성 가정)', T, 0, 'p', c.time.p, { tol: 1e-5 }), ex('측정 부분 η²', T, 0, 'peta', c.time.peta),
    ex('측정×집단 F', X, 0, 'F', c.inter.F), ex('측정×집단 유의확률 (구형성 가정)', X, 0, 'p', c.inter.p, { tol: 1e-5 }), ex('측정×집단 부분 η²', X, 0, 'peta', c.inter.peta),
    ex('개체-내 오차 제곱합', X, 0, 'ssE', c.errw.ss),
    ex('첫 셀(집단 1, 첫 측정) 평균', '기술통계', 0, 'mean', c.cells[0][0].mean), ex('시점별 Levene (첫 측정)', 'Levene', 0, 'W', c.levene[0].W),
    ex('집단별 측정 효과 (집단 1) F', '측정의 단순 효과', 0, 'F', c.time_in_group[0].F), ex('시점별 집단 효과 (첫 측정) F', '집단의 단순 효과', 0, 'F', c.group_at_time[0].F)
  ];
  if (c.k >= 3) {
    expected.push(ex("Mauchly's W", 'Mauchly', 0, 'W', c.mauchly.W), ex('Mauchly 유의확률', 'Mauchly', 0, 'p', c.mauchly.p, { tol: 1e-4 }),
      ex('ε Greenhouse–Geisser', 'Mauchly', 0, 'gg', c.eps.gg), ex('ε Huynh–Feldt (Lecoutre)', 'Mauchly', 0, 'hf', c.eps.hf),
      ex('측정 유의확률 (Greenhouse–Geisser)', T, 1, 'p', c.time.corr.gg.p, { tol: 1e-5 }), ex('측정×집단 유의확률 (Greenhouse–Geisser)', X, 1, 'p', c.inter.corr.gg.p, { tol: 1e-5 }));
  }
  const pr = c.time_in_group[0].pairs[0];
  expected.push(ex('쌍별 비교 (집단 1, 첫 쌍) 평균차이', '쌍별 비교', 0, 'diff', pr.diff), ex('쌍별 비교 (집단 1, 첫 쌍) 보정 유의확률', '쌍별 비교', 0, 'pAdj', pr.p_bonf, { tol: 1e-5 }));
  const header = names.concat(['집단']);
  entries.push({
    id: id, analysis: 'rmanova', title: '혼합 설계 분산분석: ' + c.name, description: desc + ' — ' + c.N + '명.',
    header: header, rows: c.Y.map((r, i) => r.concat([c.G[i]])),
    varRows: names.map(nm => [nm, nm + ' 점수', '연속', '', '']).concat([['집단', '집단', '명목', vl, '']]),
    request: { analysis: 'rmanova', vars: names, factor: ['집단'], options: { posthoc: true, plot: true, simple: true, levene: true, conf: '95' } },
    guide: ['분석: 반복측정 분산분석 (집단변수는 선택)', '측정 변수: ' + names.join(', ') + ' (이 순서대로 보내세요)', '집단변수: 집단 (넣으면 혼합 설계)', '옵션: 쌍별 비교, 상호작용 그림, 단순 효과, 등분산 검정 켬 (기본값)',
      '개체-간 효과(집단), 개체-내 효과 두 표(측정, 측정 × 집단; 읽을 행이 굵게), ' + (c.k >= 3 ? 'Mauchly·ε, ' : '') + '단순 효과, 집단별 쌍별 비교가 기대 결과와 같은지 본다.'],
    expected: expected
  });
});

/* ---------- 평가자 일치도 (kappa), 급내상관 (ICC) ---------- */
const agRef = REF('agreement.json');
{
  const k1 = agRef.kappa[0], k2 = agRef.kappa[1], f1 = agRef.fleiss[0];
  const mk = (id, title, desc, header, rows, varRows, request, guide, expected) => entries.push({ id: id, analysis: 'kappa', title: title, description: desc, header: header, rows: rows, varRows: varRows, request: request, guide: guide, expected: expected });
  mk('kappa-1', '평가자 일치도: ' + k1.name, '두 평가자의 예(1)/아니오(0) 평가 100건. 한쪽 범주(아니오)가 많아 kappa(.733)가 관찰 일치율(91%)보다 많이 낮다 — ' + k1.none.n + '건.',
    ['평가A', '평가B'], k1.a.map((v, i) => [v, k1.b[i]]), [['평가A', '평가자 A', '명목', '0=아니오, 1=예', ''], ['평가B', '평가자 B', '명목', '0=아니오, 1=예', '']],
    { analysis: 'kappa', vars: ['평가A', '평가B'], options: { weights: 'none', conf: '95' } },
    ['분석: 평가자 일치도 (kappa)', '평가자: 평가A, 평가B', '옵션: 가중치 없음 (기본)', '교차표의 대각선과 kappa, 표준오차, 신뢰구간, z, 관찰·우연 일치율이 기대 결과와 같은지 본다.'],
    [ex('관찰 일치율', 'Cohen kappa', 0, 'po', k1.none.po), ex('우연 일치율', 'Cohen kappa', 0, 'pe', k1.none.pe), ex('Cohen kappa', 'Cohen kappa', 0, 'kappa', k1.none.kappa), ex('표준오차', 'Cohen kappa', 0, 'se', k1.none.se),
      ex('신뢰구간 하한', 'Cohen kappa', 0, 'lo', k1.none.lo), ex('신뢰구간 상한', 'Cohen kappa', 0, 'hi', k1.none.hi), ex('z (kappa = 0)', 'Cohen kappa', 0, 'z', k1.none.z), ex('kappa 최대값', 'Cohen kappa', 0, 'kmax', k1.none.kmax),
      ex('교차표 (아니오, 아니오)', '평가 일치 교차표', 0, 'c0', k1.table[0][0])]);
  mk('kappa-2', '평가자 일치도: ' + k2.name, '두 평가자의 5점 척도 평가 120건. 옵션을 "이차 가중"으로 두면 가중 kappa가 나온다 — ' + k2.none.n + '건.',
    ['평가A', '평가B'], k2.a.map((v, i) => [v, k2.b[i]]), [['평가A', '평가자 A', '서열', '', ''], ['평가B', '평가자 B', '서열', '', '']],
    { analysis: 'kappa', vars: ['평가A', '평가B'], options: { weights: 'quadratic', conf: '95' } },
    ['분석: 평가자 일치도 (kappa)', '평가자: 평가A, 평가B', '옵션: 가중치 이차 가중', '이차 가중 kappa, 표준오차, z가 기대 결과와 같은지 본다. 가중치 없음·선형으로 바꾸면 값이 달라지는지도 본다 (가중치 없음 ' + k2.none.kappa.toFixed(3) + ', 선형 ' + k2.linear.kappa.toFixed(3) + ').'],
    [ex('이차 가중 kappa', '이차 가중 kappa', 0, 'kappa', k2.quadratic.kappa), ex('표준오차', '이차 가중 kappa', 0, 'se', k2.quadratic.se), ex('신뢰구간 하한', '이차 가중 kappa', 0, 'lo', k2.quadratic.lo),
      ex('z (kappa = 0)', '이차 가중 kappa', 0, 'z', k2.quadratic.z), ex('관찰 일치도(가중)', '이차 가중 kappa', 0, 'po', k2.quadratic.po)]);
  const names3 = ['평가A', '평가B', '평가C'];
  mk('kappa-3', '평가자 일치도: ' + f1.name, '세 평가자의 3범주 평가 60건. 평가자가 3명이라 자동으로 Fleiss kappa와 평가자 쌍별 kappa 행렬이 나온다 — ' + f1.n + '건.',
    names3, f1.R, names3.map(nm => [nm, '평가자 ' + nm.slice(2), '명목', '', '']),
    { analysis: 'kappa', vars: names3, options: { weights: 'none', conf: '95' } },
    ['분석: 평가자 일치도 (kappa)', '평가자: 평가A, 평가B, 평가C', '옵션: 기본 (3명 이상은 Fleiss kappa로 자동 전환)', '전체 Fleiss kappa, 범주별 kappa, 평가자 쌍별 Cohen kappa가 기대 결과와 같은지 본다.'],
    [ex('Fleiss kappa (전체)', 'Fleiss kappa', 3, 'kappa', f1.kappa), ex('표준오차 (영가설)', 'Fleiss kappa', 3, 'se', f1.se), ex('z', 'Fleiss kappa', 3, 'z', f1.z),
      ex('범주 1 kappa', 'Fleiss kappa', 0, 'kappa', f1.kappa_j[0]), ex('범주 3 kappa', 'Fleiss kappa', 2, 'kappa', f1.kappa_j[2]), ex('평가자 A–B 쌍 kappa', '평가자 쌍별', 0, 'c1', f1.pair[0][1]), ex('평가자 B–C 쌍 kappa', '평가자 쌍별', 1, 'c2', f1.pair[1][2])]);
  [[agRef.icc[0], 'icc-1', '4명 평가자가 30명의 대상을 연속 점수로 평가. 평가자 편향이 작아 절대 일치도와 일관성의 ICC가 비슷하다'],
    [agRef.icc[2], 'icc-2', '3명 평가자 × 12명 대상, 평가자마다 점수 수준이 크게 달라(편향 큼) 절대 일치도 ICC(2,1)가 일관성 ICC(3,1)보다 훨씬 낮다'],
    [agRef.icc[1], 'icc-3', '같은 사람 40명을 2번 측정한 검사-재검사 자료, 일치가 높다']].forEach(([c, id, desc]) => {
    const names = c.Y[0].map((_, j) => '평가' + (j + 1));
    const T = '급내상관계수';
    const e = [];
    [0, 1, 2, 3, 4, 5].forEach(q => { const lab = ['ICC(1,1)', 'ICC(2,1)', 'ICC(3,1)', 'ICC(1,k)', 'ICC(2,k)', 'ICC(3,k)'][q]; e.push(ex(lab, T, q, 'icc', c.rows[q].icc)); });
    e.push(ex('ICC(2,1) 신뢰구간 하한', T, 1, 'lo', c.rows[1].lo, { tol: 1e-5 }), ex('ICC(2,1) 신뢰구간 상한', T, 1, 'hi', c.rows[1].hi, { tol: 1e-5 }), ex('ICC(3,1) 신뢰구간 하한', T, 2, 'lo', c.rows[2].lo, { tol: 1e-5 }),
      ex('ICC(2,1) F', T, 1, 'F', c.rows[1].F), ex('ICC(2,1) 유의확률', T, 1, 'p', c.rows[1].p, { tol: 1e-5 }), ex('대상 사이 평균제곱', '이원 분산분석표', 0, 'ms', c.msb), ex('평가자 사이 평균제곱', '이원 분산분석표', 1, 'ms', c.msj), ex('잔차 평균제곱', '이원 분산분석표', 2, 'ms', c.mse));
    entries.push({
      id: id, analysis: 'icc', title: '급내상관 (ICC): ' + c.name, description: desc + ' — 대상 ' + c.n + '명, 평가 ' + c.k + '개.',
      header: names, rows: c.Y, varRows: names.map(nm => [nm, nm + ' 점수', '연속', '', '']),
      request: { analysis: 'icc', vars: names, options: { conf: '95' } },
      guide: ['분석: 급내상관 (ICC)', '평가자·측정: ' + names.join(', '), '옵션: 기본 (신뢰수준 95%)', 'ICC 6유형의 값(굵게 표시된 ICC(2,1), ICC(3,1) 포함), 신뢰구간, F와 유의확률, 이원 분산분석표의 평균제곱이 기대 결과와 같은지 본다.'],
      expected: e
    });
  });
}

/* ---------- 서열(비례오즈) 로지스틱 회귀 ---------- */
const olRef = REF('ologit.json');
[[0, 'ologit-1', '4범주 서열 결과, 숫자 독립변수 1개와 0/1 독립변수 1개. 비례오즈 가정이 성립해 Brant 검정이 유의하지 않다'],
  [2, 'ologit-2', '3범주 서열 결과이고 한 독립변수(x)의 효과가 문턱값마다 크게 달라 **비례오즈 가정이 깨진다** (Brant 전체 검정이 유의, 경고가 나온다)'],
  [3, 'ologit-3', '7범주 서열 결과, 범주 분포가 한쪽으로 치우친 자료에서 독립변수 1개의 효과를 본다']].forEach(([idx, id, desc]) => {
  const c = olRef.cases[idx];
  const names = c.names;
  const T = '모수 추정값', M = '모형 적합 정보', Br = '비례오즈 가정 점검';
  const expected = [ex('모형 카이제곱', M, 1, 'chi2', c.chi2), ex('모형 유의확률', M, 1, 'p', c.p, { tol: 1e-4 }), ex('−2 로그 가능도 (최종)', M, 1, 'm2ll', -2 * c.ll),
    ex('Cox와 Snell', '유사 R 제곱', 0, 'cs', c.cs), ex('Nagelkerke', '유사 R 제곱', 0, 'nk', c.nk), ex('McFadden', '유사 R 제곱', 0, 'mf', c.mf)];
  c.theta.forEach((th, j) => expected.push(ex('문턱값 ' + (j + 1), T, j, 'B', th, { tol: 1e-5 })));
  c.beta.forEach((b, i) => {
    const row = c.J - 1 + i;
    expected.push(ex('B (' + names[i] + ')', T, row, 'B', b, { tol: 1e-5 }), ex('표준오차 (' + names[i] + ')', T, row, 'se', c.se[i], { tol: 1e-4 }), ex('Exp(B) (' + names[i] + ')', T, row, 'or', Math.exp(b), { tol: 1e-5 }));
  });
  expected.push(ex('Brant 전체 카이제곱', Br, c.k, 'chi2', c.brant.omnibus.chi2, { tol: 1e-5 }), ex('Brant 전체 유의확률', Br, c.k, 'p', c.brant.omnibus.p, { tol: 1e-4 }));
  entries.push({
    id: id, analysis: 'ologit', title: '서열 로지스틱 회귀: ' + c.name, description: desc + ' — ' + c.n + '건, 범주 ' + c.J + '개.',
    header: ['y'].concat(names), rows: c.y.map((v, i) => [v].concat(c.cols.map(col => col[i]))),
    varRows: [['y', '서열 결과', '서열', c.cats.map((v, i) => v + '=' + (i + 1) + '단계').join(', '), '']].concat(names.map(nm => [nm, '독립변수 ' + nm, '연속', '', ''])),
    request: { analysis: 'ologit', dependent: ['y'], covariates: names, options: { parallel: true, classify: false, save: false, conf: '95' } },
    guide: ['분석: 서열 로지스틱 회귀 (비례오즈)', '종속변수: y, 독립변수(숫자): ' + names.join(', '), '옵션: 비례오즈 가정 점검 켬 (기본), 분류표·저장 끔',
      '모형 적합 정보, 유사 R², 모수 추정값(문턱값 ' + (c.J - 1) + '개와 독립변수 B, 표준오차, Exp(B)), 비례오즈 가정 점검(Brant 전체)이 기대 결과와 같은지 본다.'],
    expected: expected
  });
});

/* ---------- 편상관 (상관분석의 통제변수 칸) ---------- */
const ptRef = REF('partial.json');
[[0, 'partcorr-1', '변수 3개(x, y, w)의 상관을 통제변수 z 하나의 영향을 뺀 피어슨 편상관으로 본다 (통제하면 x–y 상관의 부호가 바뀐다)', 'pearson', false],
  [3, 'partcorr-2', '5점 척도 정수 점수 3개를 통제변수(연령) 하나의 영향을 뺀 스피어만 순위 편상관으로 본다. 동순위가 많고 신뢰구간 옵션을 켠다', 'spearman', true]].forEach(([idx, id, desc, method, ci]) => {
  const c = ptRef.cases[idx];
  const names = c.vnames.concat(c.cnames), cols = c.V.concat(c.C), e = c[method];
  const M = method === 'pearson' ? '피어슨 편상관 행렬' : '스피어만 순위 편상관 행렬', Cmp = (method === 'pearson' ? '피어슨' : '스피어만 순위') + ' 편상관: 통제 전';
  const expected = [ex('편상관 (' + c.vnames[0] + ', ' + c.vnames[1] + ')', M, 0, 'c1', e.r[0][1]), ex('유의확률 (' + c.vnames[0] + ', ' + c.vnames[1] + ')', M, 1, 'c1', e.p[0][1], { tol: 1e-4 }),
    ex('자유도', M, 2, 'c1', e.df), ex('편상관 (' + c.vnames[0] + ', ' + c.vnames[2] + ')', M, 0, 'c2', e.r[0][2]), ex('편상관 (' + c.vnames[1] + ', ' + c.vnames[2] + ')', M, 3, 'c2', e.r[1][2]),
    ex('통제 전 상관 (' + c.vnames[0] + ', ' + c.vnames[1] + ')', Cmp, 0, 'r0', e.zero[0][1]), ex('통제 후 편상관 (비교표)', Cmp, 0, 'r', e.r[0][1])];
  if (ci) expected.push(ex('신뢰구간 하한 (' + c.vnames[0] + ', ' + c.vnames[1] + ')', Cmp, 0, 'lo', e.lo[0][1], { tol: 1e-4 }), ex('신뢰구간 상한 (' + c.vnames[0] + ', ' + c.vnames[1] + ')', Cmp, 0, 'hi', e.hi[0][1], { tol: 1e-4 }));
  entries.push({
    id: id, analysis: 'correlation', title: '편상관: ' + c.name, description: desc + ' — ' + c.n + '건.',
    header: names, rows: cols[0].map((_, i) => cols.map(col => col[i])), varRows: names.map(nm => [nm, nm + ' 점수', '연속', '', '']),
    request: { analysis: 'correlation', vars: c.vnames, controls: c.cnames, options: { pearson: method === 'pearson', spearman: method === 'spearman', assumptions: false, ci: ci, scatter: false, conf: '95' } },
    guide: ['분석: 상관분석', '변수: ' + c.vnames.join(', ') + ', 통제변수: ' + c.cnames.join(', '), '옵션: ' + (method === 'pearson' ? '피어슨' : '스피어만 순위상관') + (ci ? ', 피어슨 신뢰구간 옵션(편상관에서는 해당 상관의 Fisher z 구간)' : ''),
      '편상관 행렬(계수, 유의확률, 자유도)과 "통제 전 → 통제 후 비교" 표의 값이 기대 결과와 같은지 본다. 통제변수를 비우고 같은 변수로 실행하면 통제 전 상관과 같은 값이 나오는지도 본다.'],
    expected: expected
  });
});

/* ---------- Kaplan–Meier + log-rank ---------- */
const svRef = REF('km.json');
/** 기준값(statsmodels 생존함수)에서 시점 t 의 생존율·표준오차·신뢰구간 (사건 시점의 단계함수) */
function svAt(g, t) {
  let k = -1;
  g.tt.forEach((x, i) => { if (x <= t) k = i; });
  return k < 0 ? null : { s: g.S[k], se: g.se[k], lo: g.lo[k], hi: g.hi[k] };
}
[['km-ties-2g', 'km-1', '정수 시간(동점이 많음)과 두 집단(크기 60, 40). log-rank 세 가지 검정과 중앙값·신뢰구간을 본다', { 1: '집단 A', 2: '집단 B' }, [4, 8, 12], false],
  ['km-3g', 'km-2', '세 집단(크기 70, 50, 30), 90% 신뢰구간. 전체 검정과 Bonferroni 쌍별 검정, 시점별 생존율을 본다', { 1: '대조', 2: '약 A', 3: '약 B' }, [5, 10, 20], false],
  ['km-heavy-cens', 'km-3', '집단 없이 하나의 곡선. 중도절단이 47.5%로 많고(80명 중 사건 42건) 신뢰수준 99%. 사건 시점별 생존표 옵션을 켠다', null, [10, 20, 30], true]].forEach(([cname, id, desc, glabels, times, fullTable]) => {
  const c = svRef.cases.find(x => x.name === cname);
  const grouped = !!c.group;
  const groups = grouped ? c.groups : [c.total];
  const S = '처리 요약', At = '시점별 생존율', Q = '생존시간 백분위수', Cmp = '집단 비교 검정', Pw = '집단 쌍별';
  const expected = [];
  groups.forEach((g, q) => {
    const nm = grouped ? glabels[q + 1] : '전체';
    expected.push(ex('사건 수 (' + nm + ')', S, q, 'events', g.events, { tol: 0 }), ex('제한 평균 생존시간 (' + nm + ')', S, q, 'mean', g.rmst, { tol: 1e-6 }));
    const md = g.quantiles['0.5'];
    if (md.t !== null && !md.tie) {
      expected.push(ex('중앙값 (' + nm + ')', Q, q * 3 + 1, 't', md.t, { tol: 1e-9 }));
      if (md.lo !== null) expected.push(ex('중앙값 하한 (' + nm + ')', Q, q * 3 + 1, 'lo', md.lo, { tol: 1e-9 }));
      if (md.hi !== null) expected.push(ex('중앙값 상한 (' + nm + ')', Q, q * 3 + 1, 'hi', md.hi, { tol: 1e-9 }));
    }
    times.forEach((t, k) => {
      const a = svAt(g, t);
      if (a && a.s > 0 && t <= g.tau) expected.push(ex('생존율 ' + nm + ' t=' + t, At, q * times.length + k, 's', a.s, { tol: 1e-6 }), ex('Greenwood 표준오차 ' + nm + ' t=' + t, At, q * times.length + k, 'se', a.se, { tol: 1e-5 }),
        ex('신뢰구간 하한 ' + nm + ' t=' + t, At, q * times.length + k, 'lo', a.lo, { tol: 1e-5 }), ex('신뢰구간 상한 ' + nm + ' t=' + t, At, q * times.length + k, 'hi', a.hi, { tol: 1e-5 }));
    });
  });
  if (grouped) {
    expected.push(ex('log-rank 카이제곱', Cmp, 0, 'chi2', c.compare.logrank.chi2, { tol: 1e-6 }), ex('log-rank 유의확률', Cmp, 0, 'p', c.compare.logrank.p, { tol: 1e-4 }),
      ex('Breslow 카이제곱', Cmp, 1, 'chi2', c.compare.breslow.chi2, { tol: 1e-6 }), ex('Tarone–Ware 카이제곱', Cmp, 2, 'chi2', c.compare.tarone.chi2, { tol: 1e-6 }), ex('Tarone–Ware 유의확률', Cmp, 2, 'p', c.compare.tarone.p, { tol: 1e-4 }));
    if (groups.length >= 3) c.compare.pairs.forEach((pr, k) => expected.push(ex('쌍별 log-rank 카이제곱 (' + glabels[pr.i + 1] + ' vs ' + glabels[pr.j + 1] + ')', Pw, k, 'chi2', pr.chi2, { tol: 1e-6 }),
      ex('보정 유의확률 (' + glabels[pr.i + 1] + ' vs ' + glabels[pr.j + 1] + ')', Pw, k, 'pAdj', Math.min(1, pr.p * c.compare.pairs.length), { tol: 1e-4 })));
  }
  if (fullTable) expected.push(ex('사건 시점별 생존표 첫 행 생존율', '사건 시점별 생존표', 0, 's', c.total.S[0], { tol: 1e-9 }), ex('사건 시점별 생존표 마지막 행 생존율', '사건 시점별 생존표', c.total.tt.length - 1, 's', c.total.S[c.total.S.length - 1], { tol: 1e-9 }));
  const header = grouped ? ['집단', '시간', '사건'] : ['시간', '사건'];
  const rows = c.time.map((t, i) => (grouped ? [c.group[i] + 1, t, c.event[i]] : [t, c.event[i]]));
  const varRows = [['시간', '사건 또는 관찰 종료까지의 시간', '연속', '', ''], ['사건', '사건 여부', '명목', '1=사건 발생, 0=중도절단', '']];
  if (grouped) varRows.unshift(['집단', '집단', '명목', Object.keys(glabels).map(k => k + '=' + glabels[k]).join(', '), '']);
  entries.push({
    id: id, analysis: 'km', title: 'Kaplan–Meier: ' + c.note, description: desc + ' — ' + c.time.length + '건.',
    header: header, rows: rows, varRows: varRows,
    request: { analysis: 'km', time: ['시간'], event: ['사건'], group: grouped ? ['집단'] : [], options: { eventValue: '1', times: times.join(', '), tests: true, pairs: true, plot: true, table: fullTable, conf: String(c.conf) } },
    guide: ['분석: Kaplan–Meier 생존분석 (log-rank)', '시간 변수: 시간, 사건 변수: 사건' + (grouped ? ', 집단변수: 집단' : ' (집단변수는 비움)'), '옵션: 사건으로 볼 값 1, 생존율을 볼 시점 ' + times.join(', ') + (fullTable ? ', 사건 시점별 생존표 켬' : '') + ', 신뢰수준 ' + c.conf + '%',
      '처리 요약(사건 수, 제한 평균), 시점별 생존율(Greenwood 표준오차와 신뢰구간), 생존시간 백분위수의 중앙값' + (grouped ? ', 집단 비교 검정' + (groups.length >= 3 ? '과 쌍별 비교' : '') : '') + '이 기대 결과와 같은지 본다.'],
    expected: expected
  });
});

/* ---------- 일원 분산분석: LSD, Šidák, Dunnett 사후검정 ---------- */
const anRef = REF('anova.json');
[[0, 'anova-1', '크기와 분산이 다른 세 집단(12, 25, 18건). LSD·Šidák·Dunnett을 켜고 Dunnett의 통제집단은 비워 첫 집단으로 둔다', '', '95'],
  [5, 'anova-2', '크기가 작고 다른 세 집단(3, 4, 3건, 오차 자유도 7). 신뢰수준 99%이고 Dunnett의 통제집단을 2번 집단으로 지정한다', '2', '99'],
  [6, 'anova-3', '크기가 4~30건으로 다른 7집단. Dunnett의 통제집단을 4번 집단으로 지정하면 나머지 6집단과 비교한다', '4', '95']].forEach(([idx, id, desc, control, conf]) => {
  const c = anRef.cases[idx], k = c.groups.length;
  const header = ['집단', '값'];
  const rows = [];
  c.groups.forEach((g, q) => g.forEach(v => rows.push([q + 1, v])));
  const dn = c.dunnett[conf][control === '' ? '0' : String(Number(control) - 1)].rows;
  const L = c.posthoc[conf].lsd, S = c.posthoc[conf].sidak;
  const expected = [ex('LSD 평균차이 (집단1−집단2)', '다중비교: LSD', 0, 'diff', L[0].diff), ex('LSD 유의확률 (집단1−집단2)', '다중비교: LSD', 0, 'p', L[0].p, { tol: 1e-5 }),
    ex('LSD 신뢰구간 하한 (집단1−집단2)', '다중비교: LSD', 0, 'lo', L[0].lo), ex('Šidák 유의확률 (집단1−집단2)', '다중비교: Šidák', 0, 'p', S[0].p, { tol: 1e-5 }),
    ex('Šidák 신뢰구간 상한 (집단1−집단2)', '다중비교: Šidák', 0, 'hi', S[0].hi)];
  dn.forEach((r, q) => expected.push(ex('Dunnett 평균차이 (집단' + (r.i + 1) + '−통제)', '다중비교: Dunnett', q, 'diff', r.diff), ex('Dunnett 유의확률 (집단' + (r.i + 1) + ' vs 통제)', '다중비교: Dunnett', q, 'p', r.p, { tol: 1e-4 }),
    ex('Dunnett 신뢰구간 하한 (집단' + (r.i + 1) + ')', '다중비교: Dunnett', q, 'lo', r.lo, { tol: 1e-5 }), ex('Dunnett 신뢰구간 상한 (집단' + (r.i + 1) + ')', '다중비교: Dunnett', q, 'hi', r.hi, { tol: 1e-5 })));
  const ctrlNo = control === '' ? '1' : control;
  entries.push({
    id: id, analysis: 'anova', title: '일원 분산분석 사후검정: ' + c.name, description: desc + ' — ' + rows.length + '건, ' + k + '집단.',
    header: header, rows: rows,
    varRows: [['집단', '집단', '명목', c.groups.map((_, q) => (q + 1) + '=집단' + (q + 1)).join(', '), ''], ['값', '측정값', '연속', '', '']],
    request: { analysis: 'anova', dependent: ['값'], factor: ['집단'], options: { levene: false, effect: false, assumptions: false, lsd: true, sidak: true, dunnett: true, control: control, conf: conf } },
    guide: ['분석: 일원 분산분석 (ANOVA)', '종속변수: 값, 집단변수: 집단', '세부 옵션: 사후검정 LSD, Šidák, Dunnett를 켠다' + (control === '' ? ' (Dunnett의 통제집단은 비움 = 첫 집단)' : ', Dunnett의 통제집단에 ' + control + ' 입력') + ', 신뢰수준 ' + conf,
      '다중비교 표 3개의 값이 기대 결과와 같은지 본다. Dunnett 표는 통제집단(' + ctrlNo + '번)을 제외한 집단만 (I) 집단으로 나오고 (J) 집단이 항상 통제집단이어야 한다.'],
    expected: expected
  });
});

/* ---------- 이항 검정 ---------- */
const bnRef = REF('binom.json');
[[0, 'binom-1', '40명 중 27명이 성공(합격). 기준 비율 0.5와 다른지 양측으로 본다', 'two'],
  [5, 'binom-2', '200명 중 18명이 성공(9%). 기준 비율 0.15보다 **작은지** 단측으로 본다', 'less'],
  [6, 'binom-3', '48명 중 45명이 성공. 기준 비율 0.8과 다른지 양측으로, 신뢰수준 99%로 본다', 'two']].forEach(([idx, id, desc, alt]) => {
  const c = bnRef.cases[idx];
  const rows = [];
  for (let i = 0; i < c.n; i++) rows.push([i < c.k ? 1 : 0]);
  const T = '이항 검정', Ci = '성공 비율의';
  entries.push({
    id: id, analysis: 'binom', title: '이항 검정: ' + c.k + '/' + c.n + ', 검정 비율 ' + c.p0, description: desc + ' — ' + c.n + '건.',
    header: ['합격'], rows: rows, varRows: [['합격', '합격 여부', '명목', '1=합격, 0=불합격', '']],
    request: { analysis: 'binom', vars: ['합격'], options: { success: '1', p0: String(c.p0), tails: alt, conf: String(c.conf) } },
    guide: ['분석: 이항 검정 (성공 비율 검정)', '변수: 합격', '옵션: 성공으로 볼 값 1, 검정 비율 ' + c.p0 + ', 검정 방향 ' + (alt === 'two' ? '양측' : '단측: 성공 비율이 검정 비율보다 작다') + ', 신뢰수준 ' + c.conf,
      '빈도(성공 ' + c.k + ', 실패 ' + (c.n - c.k) + '), 이항 검정의 정확 유의확률, 정확(Clopper–Pearson)·Wilson 신뢰구간이 기대 결과와 같은지 본다.'],
    expected: [ex('성공 수', '성공·실패', 0, 'x', c.k, { tol: 0 }), ex('실패 수', '성공·실패', 0, 'y', c.n - c.k, { tol: 0 }),
      ex('관측 비율', T, 0, 'ph', c.k / c.n, { tol: 1e-9 }), ex('정확 유의확률 (' + (alt === 'two' ? '양측' : '단측') + ')', T, 0, 'p', c[alt], { tol: 1e-6 }),
      ex('Clopper–Pearson 하한', Ci, 0, 'lo', c.exact[0], { tol: 1e-6 }), ex('Clopper–Pearson 상한', Ci, 0, 'hi', c.exact[1], { tol: 1e-6 }),
      ex('Wilson 하한', Ci, 0, 'wlo', c.wilson[0], { tol: 1e-6 }), ex('Wilson 상한', Ci, 0, 'whi', c.wilson[1], { tol: 1e-6 })]
  });
});

/* ---------- Cochran Q ---------- */
const cqRef = REF('cochran.json');
[[0, 'cochran-1', '30명을 3번 측정한 합격(1)/불합격(0). 성공 비율이 측정마다 오르고 Q 검정이 유의하다'],
  [1, 'cochran-2', '50명을 4번 측정. 마지막 측정만 성공 비율이 높다. 보정한 쌍별 비교에서 1-4회, 2-4회 쌍이 유의하다'],
  [4, 'cochran-3', '200명을 6번 측정(쌍 15개). Q가 매우 크고 쌍별 보정 유의확률이 몇 쌍만 유의하다']].forEach(([idx, id, desc]) => {
  const c = cqRef.cases[idx];
  const names = c.x[0].map((_, j) => (j + 1) + '회');
  const m = c.k * (c.k - 1) / 2, last = c.pairs[c.pairs.length - 1], first = c.pairs[0];
  const T = 'Cochran Q 검정', P = '쌍별 비교';
  entries.push({
    id: id, analysis: 'cochran', title: 'Cochran Q 검정: ' + c.name, description: desc + ' — ' + c.n + '명, ' + c.k + '회.',
    header: names, rows: c.x, varRows: names.map(nm => [nm, nm + ' 측정', '명목', '1=예, 0=아니오', '']),
    request: { analysis: 'cochran', vars: names, options: { success: '1', posthoc: true, plot: true } },
    guide: ['분석: Cochran Q 검정 (예/아니오, 3회 이상 반복측정)', '측정 변수: ' + names.join(', ') + ' (순서대로)', '옵션: 성공으로 볼 값 1, 쌍별 비교·그림 켬 (기본)',
      '측정별 성공 빈도, Cochran Q·자유도·유의확률·Kendall의 W, 쌍별 McNemar(정확, Bonferroni ×' + m + ') 표가 기대 결과와 같은지 본다.'],
    expected: [ex('Cochran Q', T, 0, 'q', c.q, { tol: 1e-6 }), ex('자유도', T, 0, 'df', c.df, { tol: 0 }), ex('유의확률', T, 0, 'p', c.p, { tol: 1e-4 }), ex("Kendall의 W", T, 0, 'w', c.w, { tol: 1e-6 }),
      ex('첫 측정의 성공 수', '측정별 성공 빈도', 0, 'x', c.colsum[0], { tol: 0 }), ex('마지막 측정의 성공 수', '측정별 성공 빈도', c.k - 1, 'x', c.colsum[c.k - 1], { tol: 0 }),
      ex('첫 쌍 정확 유의확률', P, 0, 'p', first.p, { tol: 1e-5 }), ex('첫 쌍 보정 유의확률', P, 0, 'pAdj', Math.min(1, first.p * m), { tol: 1e-5 }),
      ex('마지막 쌍 정확 유의확률', P, m - 1, 'p', last.p, { tol: 1e-5 }), ex('마지막 쌍 보정 유의확률', P, m - 1, 'pAdj', Math.min(1, last.p * m), { tol: 1e-5 })]
  });
});

/* ---------- 켄달 타우-b (상관분석의 옵션) ---------- */
const kdRef = REF('kendall.json');
[[3, 'kendall-1', '동순위가 없는 연속 점수 20쌍. 정확 분포로 유의확률을 구한다 (N이 33 이하이고 동순위 없음)'],
  [10, 'kendall-2', '1~5점 정수 점수 120쌍으로 동순위가 매우 많다. 동순위 보정 분산의 정규 근사로 유의확률을 구한다']].forEach(([idx, id, desc]) => {
  const c = kdRef.cases[idx];
  const M = '켄달 타우-b 상관 행렬';
  entries.push({
    id: id, analysis: 'correlation', title: '켄달 타우-b: ' + c.name, description: desc + ' — ' + c.n + '건.',
    header: ['x', 'y'], rows: c.x.map((v, i) => [v, c.y[i]]), varRows: [['x', 'x 점수', '연속', '', ''], ['y', 'y 점수', '연속', '', '']],
    request: { analysis: 'correlation', vars: ['x', 'y'], options: { pearson: false, spearman: false, kendall: true, assumptions: false, ci: false, scatter: false, conf: '95' } },
    guide: ['분석: 상관분석', '변수: x, y (통제변수는 비움)', '옵션: 피어슨을 끄고 켄달 타우-b 순위상관을 켠다 (가정 점검·산점도 끔)',
      '켄달 타우-b 상관 행렬의 상관계수와 유의확률(양측)이 기대 결과와 같은지 본다. ' + (c.ties ? '동순위가 많아 각주에 정규 근사로 구했다는 설명이 나온다.' : '동순위가 없고 N이 33 이하라 각주에 정확 분포로 구했다는 설명이 나온다.')],
    expected: [ex('켄달 타우-b', M, 0, 'c1', c.tau, { tol: 1e-6 }), ex('유의확률 (양측)', M, 1, 'c1', c.p_auto, { tol: 1e-5 }), ex('N', M, 2, 'c1', c.n, { tol: 0 })]
  });
});

/* ---------- IPA (중요도–수행도 분석) ---------- */
const ipRef = REF('ipa.json');
[[0, 0, 'ipa-1', '5문항, 100명의 5점 척도 응답. 기준선은 기본(문항 평균들의 평균)이고 집중 개선·현상 유지·낮은 우선순위·과잉 투자 문항이 모두 나온다'],
  [1, 1, 'ipa-2', '8문항, 150명, 약 6% 결측. 기준선 방식을 중앙값으로 하고 중요도 −0.1, 수행도 +0.2만큼 이동한다 (기본 기준선과 사분면이 달라지는 문항이 있다)'],
  [2, 2, 'ipa-3', '4문항, 25명. 기준선을 직접 지정한다 (중요도 3.0, 수행도 3.5). 이동량은 무시된다']].forEach(([ci, si, id, desc]) => {
  const c = ipRef.cases[ci], e = c.results[si], k = c.k;
  const inames = Array.from({ length: k }, (_, j) => '중요' + (j + 1)), pnames = Array.from({ length: k }, (_, j) => '수행' + (j + 1));
  const It = '문항별', Ln = '기준선', Qd = '사분면별', Gt = '격차 검정';
  const sc = e.scenario;
  const options = { baseline: sc.baseline, gap: true, plot: true, impShift: String(sc.impShift || 0), perfShift: String(sc.perfShift || 0), impLine: sc.baseline === 'custom' ? String(sc.impLine) : '', perfLine: sc.baseline === 'custom' ? String(sc.perfLine) : '', listwise: false };
  const expected = [ex('중요도 최종 기준선', Ln, 0, 'line', e.line_imp, { tol: 1e-6 }), ex('수행도 최종 기준선', Ln, 1, 'line', e.line_perf, { tol: 1e-6 })];
  e.items.forEach((x, j) => expected.push(ex('문항 ' + (j + 1) + ' 중요도 평균', It, j, 'mi', x.imp.mean, { tol: 1e-6 }), ex('문항 ' + (j + 1) + ' 수행도 평균', It, j, 'mp', x.perf.mean, { tol: 1e-6 })));
  ['집중 개선', '현상 유지', '낮은 우선순위', '과잉 투자'].forEach((q, qi) => expected.push(ex(q + ' 문항 수', Qd, qi, 'n', e.items.filter(x => x.quadrant === q).length, { tol: 0 })));
  const x0 = e.items[0];
  if (x0.p !== null) expected.push(ex('문항 1 격차 검정 t', Gt, 0, 't', x0.t, { tol: 1e-6 }), ex('문항 1 격차 검정 유의확률', Gt, 0, 'p', x0.p, { tol: 1e-4 }), ex('문항 1 효과크기 d_z', Gt, 0, 'dz', x0.dz, { tol: 1e-6 }));
  entries.push({
    id: id, analysis: 'ipa', title: 'IPA 중요도–수행도: ' + c.name, description: desc + ' — ' + c.n + '명, ' + k + '문항.',
    header: inames.concat(pnames), rows: c.imp.map((r, i) => r.concat(c.perf[i])), varRows: inames.map((nm, j) => [nm, '문항' + (j + 1), '서열', '', '']).concat(pnames.map((nm, j) => [nm, '문항' + (j + 1) + '(수행도)', '서열', '', ''])),
    request: { analysis: 'ipa', importance: inames, performance: pnames, options: options },
    guide: ['분석: IPA 중요도–수행도 분석', '중요도 문항: ' + inames.join(', ') + ' (순서대로), 수행도 문항: ' + pnames.join(', ') + ' (같은 순서)',
      '옵션: 기준선 방식 ' + (sc.baseline === 'mean' ? '문항 평균들의 평균 (기본)' : sc.baseline === 'median' ? '문항 평균들의 중앙값, 세부 옵션의 중요도 이동 ' + sc.impShift + ', 수행도 이동 ' + sc.perfShift : '직접 지정, 세부 옵션의 중요도 기준선 값 ' + sc.impLine + ', 수행도 기준선 값 ' + sc.perfLine),
      '기준선 표의 최종 값, 문항별 평균, 사분면별 문항 수, 격차 검정이 기대 결과와 같은지 보고 그림에 십자 기준선과 문항 이름이 보이는지 확인한다.'],
    expected: expected
  });
});

/* ---------- 파일 쓰기 ---------- */
function build() {
  const body = JSON.stringify(entries);
  return `/**
 * 시험 데이터 레지스트리 (자동 생성 파일 — 직접 고치지 말 것).
 * 만드는 법: node tests/gen_test_datasets.js  (기준값 JSON과 이 생성기의 항목에서 만든다)
 * 시트에서는 메뉴 "Stat4SW 통계분석 > 설정·관리 > 시험 데이터"(고급 분석을 켠 뒤)로 넣는다 (TestTool.js).
 */
var TEST_ANALYSIS_LABELS_ = ${JSON.stringify(ANALYSIS_LABELS)};
var TEST_SAMPLE_COVERED_ = ${JSON.stringify(SAMPLE_COVERED)};
var TEST_NO_DATA_ = ${JSON.stringify(NO_DATA)};
var TEST_DATASETS_ = ${body};

if (typeof module !== 'undefined') {
  module.exports = { TEST_DATASETS_: TEST_DATASETS_, TEST_ANALYSIS_LABELS_: TEST_ANALYSIS_LABELS_, TEST_SAMPLE_COVERED_: TEST_SAMPLE_COVERED_, TEST_NO_DATA_: TEST_NO_DATA_ };
}
`;
}

if (require.main === module) {
  fs.writeFileSync(OUT, build());
  console.log('wrote', OUT, entries.length, 'datasets,', Math.round(fs.statSync(OUT).size / 1024) + ' KB');
}

module.exports = { entries, build, SAMPLE_COVERED, NO_DATA, ANALYSIS_LABELS };
