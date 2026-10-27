/**
 * 매개효과 분석 (Hayes PROCESS 모델 4, 병렬 다중 매개): 순수 함수 (dataset, request) → 결과 객체.
 *
 *   M_i = a0_i + a_i·X (+ 공변량)                      (a 경로)
 *   Y   = c0' + c'·X + Σ b_i·M_i (+ 공변량)             (b 경로, 직접효과 c')
 *   Y   = c0 + c·X (+ 공변량)                          (총효과 c = c' + Σ a_i·b_i)
 *
 * 특정 간접효과 = a_i·b_i, 총 간접효과 = Σ a_i·b_i. 유의성은 **부트스트랩 퍼센타일 신뢰구간**(기본 5,000회, 시드로 재현 가능)이 기본이고
 * Sobel 검정은 참고로 보여 준다. 회귀 계산은 Stats_Regression.js의 olsFit(QR)을 쓰고, 부트스트랩은 표본마다 교차곱 합과
 * 촐레스키 풀이로 계수만 빠르게 구한다. 통계적 매개 검정은 인과관계를 증명하지 않는다.
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function medDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    var desc = require('./Stats_Descriptive.js');
    var reg = require('./Stats_Regression.js');
    return {
      getVars: data.getVars, pnorm: dist.pnorm, pt: dist.pt, centralSums: desc.centralSums, olsFit: reg.olsFit,
      pathDiagramBlock: require('./Charts.js').pathDiagramBlock
    };
  }
  return { getVars: getVars, pnorm: pnorm, pt: pt, centralSums: centralSums, olsFit: olsFit, pathDiagramBlock: pathDiagramBlock };
}

var MED_BOOT_DEFAULT_ = 5000;
var MED_BOOT_MIN_ = 500;
var MED_BOOT_MAX_ = 20000;
var MED_BUDGET_ = 2.5e7;          // 부트스트랩 횟수 × 사례 수의 상한 (큰 표본에서는 횟수를 줄인다)
var MED_BOOT_FLOOR_ = 1000;       // 예산 때문에 줄여도 이 횟수 아래로는 내리지 않는다
var MED_MAX_MEDIATORS_ = 8;
function medName_(v) { return v.label ? v.label : v.name; }

/** 32비트 난수 생성기 mulberry32 (시드가 같으면 항상 같은 수열) */
function mulberry32_(seed) {
  var a = seed | 0;
  return function () {
    a = (a + 0x6D2B79F5) | 0;
    var t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** 정렬된 배열의 분위수 (선형 보간, numpy 기본과 같음). q 는 0~1 */
function quantileSorted_(sorted, q) {
  var pos = (sorted.length - 1) * q, lo = Math.floor(pos), hi = Math.ceil(pos);
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo);
}

/** 대칭 양의 정부호 행렬 A(p×p)로 A·x = r 을 촐레스키로 푼다. rhs: r 벡터들의 배열. 실패(특이)하면 null. */
function choleskySolve_(A, rhs) {
  var p = A.length, L = [], i, j, k;
  var maxDiag = 0;
  for (i = 0; i < p; i++) if (A[i][i] > maxDiag) maxDiag = A[i][i];
  for (i = 0; i < p; i++) {
    L.push(new Array(p).fill(0));
    for (j = 0; j <= i; j++) {
      var s = A[i][j];
      for (k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      if (i === j) {
        if (!(s > 1e-12 * maxDiag)) return null;
        L[i][i] = Math.sqrt(s);
      } else {
        L[i][j] = s / L[j][j];
      }
    }
  }
  return rhs.map(function (r) {
    var z = new Array(p), x = new Array(p);
    for (i = 0; i < p; i++) {
      var s2 = r[i];
      for (k = 0; k < i; k++) s2 -= L[i][k] * z[k];
      z[i] = s2 / L[i][i];
    }
    for (i = p - 1; i >= 0; i--) {
      var s3 = z[i];
      for (k = i + 1; k < p; k++) s3 -= L[k][i] * x[k];
      x[i] = s3 / L[i][i];
    }
    return x;
  });
}

/**
 * 부트스트랩: 사례를 복원추출해 a_i·b_i 를 B번 계산한다.
 * x, ms(배열의 배열), covs(배열의 배열), y: 같은 길이의 숫자 배열. 반환: { ind: [[B회 a_i·b_i] × m], total: [B회 합], invalid }
 * 변수를 전체 평균으로 중심화한 뒤 표본마다 교차곱 합으로 정규방정식을 푼다 (평균을 옮겨도 기울기는 같다).
 */
function mediationBootstrap(x, ms, covs, y, B, seed) {
  var n = y.length, m = ms.length, c = covs.length, i, j, t, b;
  var vars = [x].concat(ms, covs, [y]);                      // 0: X, 1..m: M, m+1..m+c: 공변량, 마지막: Y
  var K = vars.length;
  var cols = vars.map(function (v) {
    var s = 0;
    for (t = 0; t < n; t++) s += v[t];
    var mean = s / n, out = new Float64Array(n);
    for (t = 0; t < n; t++) out[t] = v[t] - mean;
    return out;
  });
  var rnd = mulberry32_(seed);
  var idx = new Int32Array(n);
  var sums = new Float64Array(K), prod = new Float64Array(K * K);
  var ind = [], total = [], invalid = 0;
  for (i = 0; i < m; i++) ind.push([]);
  var yi = K - 1;
  var xcIdx = [0].concat(covs.map(function (cc, k) { return m + 1 + k; }));          // a 경로의 설명변수: X, 공변량
  var yIdx = [0].concat(ms.map(function (mm, k) { return 1 + k; }), covs.map(function (cc, k) { return m + 1 + k; }));   // b 경로: X, M, 공변량
  for (b = 0; b < B; b++) {
    for (t = 0; t < n; t++) idx[t] = Math.floor(rnd() * n);
    for (i = 0; i < K; i++) sums[i] = 0;
    for (i = 0; i < K * K; i++) prod[i] = 0;
    for (t = 0; t < n; t++) {
      var r = idx[t];
      for (i = 0; i < K; i++) {
        var vi = cols[i][r];
        sums[i] += vi;
        for (j = i; j < K; j++) prod[i * K + j] += vi * cols[j][r];
      }
    }
    var S = function (p, q) {
      var lo = p <= q ? p : q, hi = p <= q ? q : p;
      return prod[lo * K + hi] - sums[lo] * sums[hi] / n;
    };
    var A1 = xcIdx.map(function (p) { return xcIdx.map(function (q) { return S(p, q); }); });
    var sol1 = choleskySolve_(A1, ms.map(function (mm, k) { return xcIdx.map(function (p) { return S(p, 1 + k); }); }));
    var A2 = yIdx.map(function (p) { return yIdx.map(function (q) { return S(p, q); }); });
    var sol2 = choleskySolve_(A2, [yIdx.map(function (p) { return S(p, yi); })]);
    if (!sol1 || !sol2) { invalid++; continue; }
    var tot = 0;
    for (i = 0; i < m; i++) {
      var ab = sol1[i][0] * sol2[0][1 + i];                // a_i · b_i
      ind[i].push(ab);
      tot += ab;
    }
    total.push(tot);
  }
  return { ind: ind, total: total, invalid: invalid };
}

function bootCi_(values, conf) {
  var s = values.slice().sort(function (a, b) { return a - b; });
  var alpha = 1 - conf / 100, n = s.length, mean = 0, ss = 0, i;
  for (i = 0; i < n; i++) mean += s[i];
  mean /= n;
  for (i = 0; i < n; i++) ss += (s[i] - mean) * (s[i] - mean);
  return { lo: quantileSorted_(s, alpha / 2), hi: quantileSorted_(s, 1 - alpha / 2), se: n > 1 ? Math.sqrt(ss / (n - 1)) : null, n: n };
}

/**
 * 핵심 계산. input = { x: [..], xName, ms: [{ name, values }], covs: [{ name, values }], y: [..], yName, conf, boot, seed }
 * (모든 배열은 결측이 없는 같은 길이). 반환: { aModels, yModel, tModel, paths, effects, sobel, boot } 또는 { error }
 */
function mediationFit(input) {
  var d = medDeps_();
  var n = input.y.length, m = input.ms.length, i;
  var conf = input.conf;
  var covCols = input.covs.map(function (c) { return { name: c.name, values: c.values }; });
  var xCol = { name: input.xName, values: input.x };

  var aModels = input.ms.map(function (mm) { return d.olsFit([xCol].concat(covCols), mm.values, conf); });
  var yModel = d.olsFit([xCol].concat(input.ms.map(function (mm) { return { name: mm.name, values: mm.values }; }), covCols), input.y, conf);
  var tModel = d.olsFit([xCol].concat(covCols), input.y, conf);
  var all = aModels.concat([yModel, tModel]);
  for (i = 0; i < all.length; i++) {
    if (all[i].tooFew || all[i].constant || all[i].collinear) return { error: all[i] };
  }

  var sx = d.centralSums(input.x), sy = d.centralSums(input.y);
  var sdX = Math.sqrt(sx.c2 / (n - 1)), sdY = Math.sqrt(sy.c2 / (n - 1));
  var a = aModels.map(function (f) { return f.coef[1]; });
  var b = [];
  for (i = 0; i < m; i++) b.push(yModel.coef[2 + i]);
  var cDirect = yModel.coef[1], cTotal = tModel.coef[1];
  var effects = [], sobel = [], totalInd = 0;
  for (i = 0; i < m; i++) {
    var ab = a[i].B * b[i].B;
    totalInd += ab;
    var se = a[i].se !== null && b[i].se !== null ? Math.sqrt(a[i].B * a[i].B * b[i].se * b[i].se + b[i].B * b[i].B * a[i].se * a[i].se) : null;
    var z = se !== null && se > 0 ? ab / se : null;
    effects.push({ name: input.ms[i].name, ab: ab, std: sdY > 0 ? ab * sdX / sdY : null });
    sobel.push({ name: input.ms[i].name, ab: ab, se: se, z: z, p: z !== null ? 2 * d.pnorm(Math.abs(z), false) : null });
  }
  var res = {
    n: n, m: m, aModels: aModels, yModel: yModel, tModel: tModel, a: a, b: b, cDirect: cDirect, cTotal: cTotal,
    effects: effects, totalIndirect: totalInd, totalStd: sdY > 0 ? totalInd * sdX / sdY : null, sobel: sobel,
    sdX: sdX, sdY: sdY, conf: conf, boot: null
  };
  if (input.boot > 0) {
    var bs = mediationBootstrap(input.x, input.ms.map(function (mm) { return mm.values; }), input.covs.map(function (c) { return c.values; }), input.y, input.boot, input.seed);
    res.boot = { B: input.boot, seed: input.seed, invalid: bs.invalid, ind: bs.ind.map(function (v) { return bootCi_(v, conf); }), total: bootCi_(bs.total, conf) };
  }
  return res;
}

function medSig_(p) { return p === null || p === undefined ? '' : p < 0.001 ? '***' : p < 0.01 ? '**' : p < 0.05 ? '*' : ''; }

function runMediation(dataset, request) {
  var d = medDeps_();
  var dep = d.getVars(dataset, request.dependent, { numeric: true });
  var ind = d.getVars(dataset, request.independent, { numeric: true });
  var med = d.getVars(dataset, request.mediators, { numeric: true });
  var covs = request.covariates && request.covariates.length ? d.getVars(dataset, request.covariates, { numeric: true }) : [];
  if (dep.length !== 1) throw new Error('종속변수는 하나만 선택하세요.');
  if (ind.length !== 1) throw new Error('독립변수는 하나만 선택하세요.');
  if (med.length < 1) throw new Error('매개변수를 하나 이상 선택하세요.');
  if (med.length > MED_MAX_MEDIATORS_) throw new Error('매개변수는 ' + MED_MAX_MEDIATORS_ + '개까지 선택할 수 있습니다.');
  var seen = {};
  seen[dep[0].name] = '종속변수';
  seen[ind[0].name] = '독립변수';
  med.forEach(function (v) {
    if (seen[v.name]) throw new Error('변수 "' + v.name + '"를 ' + seen[v.name] + '와(과) 매개변수에 동시에 쓸 수 없습니다.');
    seen[v.name] = '매개변수';
  });
  covs.forEach(function (v) {
    if (seen[v.name]) throw new Error('변수 "' + v.name + '"를 ' + seen[v.name] + '와(과) 공변량에 동시에 쓸 수 없습니다.');
    seen[v.name] = '공변량';
  });
  var opt = request.options || {};
  var conf = opt.conf === undefined || opt.conf === null || opt.conf === '' ? 95 : Number(opt.conf);
  if (!isFinite(conf) || conf <= 0 || conf >= 100) throw new Error('"신뢰수준(%)"은(는) 0 초과 100 미만이어야 합니다.');
  var boot = opt.boot === undefined || opt.boot === null || opt.boot === '' ? MED_BOOT_DEFAULT_ : Number(opt.boot);
  if (!isFinite(boot) || Math.floor(boot) !== boot || boot < MED_BOOT_MIN_ || boot > MED_BOOT_MAX_) {
    throw new Error('"부트스트랩 횟수"는 ' + MED_BOOT_MIN_ + ' 이상 ' + MED_BOOT_MAX_ + ' 이하의 정수여야 합니다.');
  }
  var seed = opt.seed === undefined || opt.seed === null || opt.seed === '' ? 12345 : Number(opt.seed);
  if (!isFinite(seed) || Math.floor(seed) !== seed) throw new Error('"시드"는 정수여야 합니다.');

  var all = dep.concat(ind, med, covs);
  var rowsOk = [];
  for (var i = 0; i < dataset.n; i++) {
    var ok = true;
    for (var a = 0; a < all.length; a++) if (dataset.data[all[a].name][i] === null) { ok = false; break; }
    if (ok) rowsOk.push(i);
  }
  var n = rowsOk.length;
  var pick = function (name) { return rowsOk.map(function (r) { return dataset.data[name][r]; }); };
  var warnings = [];
  if (dataset.n > n) warnings.push('목록별 결측 제거로 ' + (dataset.n - n) + '건을 제외했습니다 (사용한 사례 ' + n + '건).');
  var bootUsed = boot;
  if (bootUsed * n > MED_BUDGET_) {
    bootUsed = Math.max(MED_BOOT_FLOOR_, Math.floor(MED_BUDGET_ / n));
    if (bootUsed < boot) warnings.push('사례가 ' + n + '건으로 많아 계산 시간을 줄이려고 부트스트랩 횟수를 ' + boot + '회에서 ' + bootUsed + '회로 줄였습니다.');
  }
  if (n < 50) warnings.push('사례가 ' + n + '건으로 적어 부트스트랩 신뢰구간이 불안정할 수 있습니다.');

  var input = {
    x: pick(ind[0].name), xName: medName_(ind[0]), y: pick(dep[0].name), yName: medName_(dep[0]),
    ms: med.map(function (v) { return { name: medName_(v), values: pick(v.name) }; }),
    covs: covs.map(function (v) { return { name: medName_(v), values: pick(v.name) }; }),
    conf: conf, boot: bootUsed, seed: seed
  };
  var r = mediationFit(input);
  if (r.error) {
    var e = r.error;
    if (e.tooFew) throw new Error('사례 수(' + e.n + ')가 계수의 수(' + e.p + ', 상수 포함)보다 많아야 합니다.');
    if (e.constant) throw new Error('변수 "' + e.constant + '"의 값이 모두 같아 분석에 쓸 수 없습니다.');
    throw new Error('변수 "' + e.collinear + '"이(가) 다른 변수들의 선형결합이라 계산할 수 없습니다 (완전 공선성). 같은 정보를 담은 변수를 함께 쓰지 않았는지 확인하세요.');
  }
  if (r.boot && r.boot.invalid > 0) warnings.push('부트스트랩 표본 ' + r.boot.invalid + '개에서 회귀를 풀 수 없어(공선성) 제외했습니다.');
  r.yModel.coef.forEach(function (c) {
    if (c.vif !== null && c.vif > 10) warnings.push('"' + c.name + '"의 VIF가 ' + c.vif.toFixed(1) + '로 큽니다. 독립변수와 매개변수가 거의 같은 정보를 담고 있으면 b 경로 추정이 불안정합니다.');
  });

  var m = med.length;
  var F3 = { type: 'fixed', digits: 3 };
  var blocks = [];

  // 모형 요약
  var modelRows = [];
  r.aModels.forEach(function (f, k) {
    modelRows.push({ model: input.ms[k].name + ' ~ ' + input.xName + (covs.length ? ' + 공변량' : '') + ' (a 경로)', r: f.r, r2: f.r2, F: f.F, df1: f.dfm, df2: f.dfr, p: f.pF });
  });
  modelRows.push({ model: input.yName + ' ~ ' + input.xName + ' + ' + input.ms.map(function (mm) { return mm.name; }).join(' + ') + (covs.length ? ' + 공변량' : '') + ' (b, c\' 경로)', r: r.yModel.r, r2: r.yModel.r2, F: r.yModel.F, df1: r.yModel.dfm, df2: r.yModel.dfr, p: r.yModel.pF });
  modelRows.push({ model: input.yName + ' ~ ' + input.xName + (covs.length ? ' + 공변량' : '') + ' (총효과 c)', r: r.tModel.r, r2: r.tModel.r2, F: r.tModel.F, df1: r.tModel.dfm, df2: r.tModel.dfr, p: r.tModel.pF });
  blocks.push({
    type: 'table', caption: '모형 요약',
    columns: [
      { key: 'model', label: '회귀모형' }, { key: 'r', label: 'R', format: F3 }, { key: 'r2', label: 'R 제곱', format: F3 },
      { key: 'F', label: 'F', format: F3 }, { key: 'df1', label: '자유도1', format: { type: 'int' } }, { key: 'df2', label: '자유도2', format: { type: 'int' } },
      { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
    ],
    rows: modelRows,
    notes: ['독립변수 ' + input.xName + ', 매개변수 ' + input.ms.map(function (mm) { return mm.name; }).join(', ') + ', 종속변수 ' + input.yName + (covs.length ? ', 공변량 ' + input.covs.map(function (c) { return c.name; }).join(', ') : '') + '. 유효 사례 ' + n + '건.']
  });

  // 경로 계수
  var pathRows = [];
  function pathRow(label, c) { return { path: label, B: c.B, se: c.se, t: c.t, p: c.p, lo: c.lo, hi: c.hi }; }
  r.a.forEach(function (c, k) { pathRows.push(pathRow('a' + (m > 1 ? k + 1 : '') + ': ' + input.xName + ' → ' + input.ms[k].name, c)); });
  r.b.forEach(function (c, k) { pathRows.push(pathRow('b' + (m > 1 ? k + 1 : '') + ': ' + input.ms[k].name + ' → ' + input.yName, c)); });
  pathRows.push(pathRow('c\' (직접효과): ' + input.xName + ' → ' + input.yName, r.cDirect));
  pathRows.push(pathRow('c (총효과): ' + input.xName + ' → ' + input.yName, r.cTotal));
  blocks.push({
    type: 'table', caption: '경로 계수',
    columns: [
      { key: 'path', label: '경로' }, { key: 'B', label: 'B', format: F3 }, { key: 'se', label: '표준오차', format: F3 },
      { key: 't', label: 't', format: F3 }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } },
      { key: 'lo', label: conf + '% 신뢰구간 하한', format: F3 }, { key: 'hi', label: conf + '% 신뢰구간 상한', format: F3 }
    ],
    rows: pathRows,
    notes: ['a: 독립변수가 매개변수에 미치는 효과, b: 독립변수를 통제한 매개변수의 효과, c\': 매개변수를 통제한 독립변수의 직접효과, c: 독립변수의 총효과 (c = c\' + 간접효과의 합). B는 비표준화 계수이다.']
  });

  // 효과 분해
  var effRows = [];
  var bootNote;
  effRows.push({ _bold: false, effect: '총효과 (c)', B: r.cTotal.B, se: r.cTotal.se, lo: r.cTotal.lo, hi: r.cTotal.hi, std: null, ratio: null, method: 't 분포', sig: medSig_(r.cTotal.p) });
  effRows.push({ _bold: false, effect: '직접효과 (c\')', B: r.cDirect.B, se: r.cDirect.se, lo: r.cDirect.lo, hi: r.cDirect.hi, std: null, ratio: null, method: 't 분포', sig: medSig_(r.cDirect.p) });
  var totalOk = r.cTotal.B !== 0 && isFinite(r.cTotal.B);
  var excludes0 = function (ci) { return ci && ci.lo !== null && (ci.lo > 0 || ci.hi < 0); };
  if (r.boot) {
    if (m > 1) {
      effRows.push({ _bold: true, effect: '총 간접효과 (합)', B: r.totalIndirect, se: r.boot.total.se, lo: r.boot.total.lo, hi: r.boot.total.hi, std: r.totalStd, ratio: totalOk ? r.totalIndirect / r.cTotal.B : null, method: '부트스트랩', sig: excludes0(r.boot.total) ? '*' : '' });
    }
    r.effects.forEach(function (ef, k) {
      effRows.push({ _bold: true, effect: (m > 1 ? '간접효과: ' : '간접효과 (a·b): ') + input.xName + ' → ' + ef.name + ' → ' + input.yName, B: ef.ab, se: r.boot.ind[k].se, lo: r.boot.ind[k].lo, hi: r.boot.ind[k].hi, std: ef.std, ratio: totalOk ? ef.ab / r.cTotal.B : null, method: '부트스트랩', sig: excludes0(r.boot.ind[k]) ? '*' : '' });
    });
    bootNote = '간접효과의 유의성은 부트스트랩 퍼센타일 신뢰구간으로 판단한다 (복원추출 ' + r.boot.B + '회, 시드 ' + r.boot.seed + '). 신뢰구간이 0을 포함하지 않으면 * 로 표시했다 (굵게 표시한 행이 간접효과).';
  }
  blocks.push({
    type: 'table', caption: '효과 분해 (총효과, 직접효과, 간접효과)',
    columns: [
      { key: 'effect', label: '효과' }, { key: 'B', label: '효과 크기', format: F3 }, { key: 'se', label: '표준오차', format: F3 },
      { key: 'lo', label: conf + '% 신뢰구간 하한', format: F3 }, { key: 'hi', label: conf + '% 신뢰구간 상한', format: F3 },
      { key: 'std', label: '완전 표준화 효과', format: F3 }, { key: 'ratio', label: '총효과 대비 비율', format: F3 },
      { key: 'method', label: '방법' }, { key: 'sig', label: '', align: 'left' }
    ],
    rows: effRows,
    notes: [bootNote,
      '총효과와 직접효과의 표준오차·신뢰구간은 t 분포 기준이고 * 는 유의확률 .05 미만(** .01, *** .001)이다. 간접효과의 표준오차는 부트스트랩 표본의 표준편차이다.',
      '완전 표준화 효과 = 효과 × (독립변수의 표준편차 / 종속변수의 표준편차). 총효과 대비 비율(간접효과 / 총효과)은 총효과가 0에 가깝거나 직접효과와 부호가 다르면(억제 효과) 해석하기 어렵다.',
      '통계적 매개 검정은 시간 순서와 이론적 근거가 있을 때 의미가 있으며 인과관계를 증명하지 않는다. 직접효과가 유의하지 않다고 "완전 매개"로 단정하지 않는다.'
    ].filter(Boolean)
  });

  if (opt.sobel !== false) {
    blocks.push({
      type: 'table', caption: 'Sobel 검정 (참고)',
      columns: [
        { key: 'name', label: '매개변수' }, { key: 'ab', label: '간접효과 (a·b)', format: F3 }, { key: 'se', label: 'Sobel 표준오차', format: F3 },
        { key: 'z', label: 'Z', format: F3 }, { key: 'p', label: '유의확률(양측)', format: { type: 'pvalue' } }
      ],
      rows: r.sobel,
      notes: ['Sobel 표준오차 = √(a²·SE(b)² + b²·SE(a)²), Z = a·b / SE, 정규분포 근사이다. 간접효과의 분포는 대체로 비대칭이라 표본이 작으면 부정확하므로 부트스트랩 결과를 우선한다.']
    });
  }

  if (opt.diagram !== false) {
    function stars(c) { return Number(c.B.toFixed(3)) + medSig_(c.p); }
    blocks.push(d.pathDiagramBlock('경로도: ' + input.xName + ' → ' + input.yName, {
      x: input.xName, y: input.yName,
      mediators: input.ms.map(function (mm, k) { return { label: mm.name, a: stars(r.a[k]), b: stars(r.b[k]) }; }),
      direct: stars(r.cDirect), total: stars(r.cTotal)
    }));
  }
  return { analysis: 'mediation', title: '매개효과 분석', blocks: blocks, warnings: warnings };
}

if (typeof module !== 'undefined') {
  module.exports = {
    mulberry32_: mulberry32_,
    quantileSorted_: quantileSorted_,
    choleskySolve_: choleskySolve_,
    mediationBootstrap: mediationBootstrap,
    mediationFit: mediationFit,
    runMediation: runMediation
  };
}
