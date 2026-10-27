/**
 * 생존분석: Kaplan–Meier 생존함수(집단별)와 집단 비교 검정(log-rank, Breslow(Gehan), Tarone–Ware).
 * 순수 함수: (dataset, request) → 결과 객체. 시간 변수(숫자, 0 이상), 사건 변수(사건값과 같으면 사건, 나머지는 중도절단), 집단변수(선택).
 *
 * - 생존함수: S(t) = Π (1 − d_j/n_j), 같은 시점이면 사건이 중도절단보다 먼저 일어난 것으로 본다 (위험집단 n_j는 그 시점 이상인 사례).
 * - 표준오차: Greenwood, SE = S·√Σ d/(n(n−d)). 신뢰구간은 log-log(cloglog) 변환: S^exp(±z·SE_ll), SE_ll = √Σ d/(n(n−d)) / |ln S|.
 *   (lifelines 와 R survfit conf.type="log-log" 와 같다)
 * - 백분위 생존시간: S(t) ≤ 1−p 가 되는 가장 이른 시점. 신뢰구간은 시점별 신뢰구간이 1−p 를 지나는 시점으로 구한다
 *   (하한 = 신뢰구간 하한 곡선이 처음 1−p 이하가 되는 시점, 상한 = 상한 곡선이 처음 1−p 이하가 되는 시점; 도달하지 않으면 "도달하지 않음").
 * - 평균 생존시간: 곡선 아래 넓이를 가장 큰 관측 시간까지 적분한 제한 평균 (가장 큰 시간이 중도절단이면 그 시점까지의 값).
 * - 비교 검정: 사건 시점마다 집단별 기대 사건 수와 초기하 공분산으로 U = Σ w(d_k − E_k), V 를 만들고 χ² = U'V⁻¹U (자유도 집단 수 − 1).
 *   가중치 w: log-rank 1, Breslow(일반화 Wilcoxon, Gehan) n_j, Tarone–Ware √n_j. 3집단 이상이면 두 집단씩의 log-rank 를 Bonferroni 로 보정한다.
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function svDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    return { getVars: data.getVars, compareValues_: data.compareValues_, qnorm: dist.qnorm, pchisq: dist.pchisq, lineChartBlock: require('./Charts.js').lineChartBlock };
  }
  return { getVars: getVars, compareValues_: compareValues_, qnorm: qnorm, pchisq: pchisq, lineChartBlock: lineChartBlock };
}

var SV_MAX_GROUPS_ = 10;
var SV_TABLE_CAP_ = 500;
var SV_CURVE_CAP_ = 240;
var SV_MAX_TIMES_ = 20;
var SV_F3_ = { type: 'fixed', digits: 3 };
var SV_T_ = { type: 'trim', digits: 3 };
function svName_(v) { return v.label ? v.label : v.name; }

/* ---------- Kaplan–Meier ---------- */

/**
 * time, event(1 사건 / 0 중도절단) 배열에서 생존함수. z: 신뢰구간의 정규 분위수.
 * 반환: { n, events, censored, tmax, tt, nr, dd, S, se, lo, hi, sorted }
 *   tt: 사건이 있는 시점들, nr: 그 시점의 위험집단 수, dd: 사건 수, S: 직후 생존율, se/lo/hi: Greenwood 표준오차·log-log 신뢰구간 (S = 0 이면 null, lo=hi=0)
 */
function svKm(time, event, z) {
  var n = time.length, i;
  var idx = new Array(n);
  for (i = 0; i < n; i++) idx[i] = i;
  idx.sort(function (a, b) { return time[a] - time[b]; });
  var sorted = new Array(n), events = 0;
  for (i = 0; i < n; i++) { sorted[i] = time[idx[i]]; if (event[idx[i]]) events++; }
  var res = { n: n, events: events, censored: n - events, tmax: n ? sorted[n - 1] : 0, tt: [], nr: [], dd: [], S: [], se: [], lo: [], hi: [], sorted: sorted };
  var s = 1, g = 0, pos = 0;
  while (pos < n) {
    var t = sorted[pos], d = 0, j = pos;
    while (j < n && sorted[j] === t) { if (event[idx[j]]) d++; j++; }
    var atRisk = n - pos;
    if (d > 0) {
      s *= 1 - d / atRisk;
      g += atRisk > d ? d / (atRisk * (atRisk - d)) : Infinity;
      res.tt.push(t); res.nr.push(atRisk); res.dd.push(d); res.S.push(s);
      if (s > 0 && isFinite(g)) {
        var ll = Math.log(-Math.log(s)), sell = Math.sqrt(g) / Math.abs(Math.log(s));
        res.se.push(s * Math.sqrt(g));
        res.lo.push(Math.exp(-Math.exp(ll + z * sell)));
        res.hi.push(Math.exp(-Math.exp(ll - z * sell)));
      } else { res.se.push(null); res.lo.push(0); res.hi.push(0); }
    }
    pos = j;
  }
  return res;
}

/** 시점 t 에서 마지막으로 tt ≤ t 인 단계의 번호 (없으면 −1) */
function svStepAt_(km, t) {
  var lo = 0, hi = km.tt.length - 1, ans = -1;
  while (lo <= hi) {
    var mid = (lo + hi) >> 1;
    if (km.tt[mid] <= t) { ans = mid; lo = mid + 1; } else hi = mid - 1;
  }
  return ans;
}

/** 시간이 t 이상인 사례 수 (그 시점의 위험집단) */
function svAtRisk_(km, t) {
  var lo = 0, hi = km.sorted.length;
  while (lo < hi) {
    var mid = (lo + hi) >> 1;
    if (km.sorted[mid] < t) lo = mid + 1; else hi = mid;
  }
  return km.sorted.length - lo;
}

/** 시점 t 의 생존율: { s, se, lo, hi } (마지막 관측 뒤인데 생존율이 0 이 아니면 s = null) */
function svSurvivalAt(km, t) {
  var k = svStepAt_(km, t);
  if (k < 0) return { s: 1, se: 0, lo: null, hi: null };
  var s = km.S[k];
  if (t > km.tmax && s > 0) return { s: null, se: null, lo: null, hi: null };
  return { s: s, se: km.se[k], lo: s > 0 ? km.lo[k] : null, hi: s > 0 ? km.hi[k] : null };
}

/** 백분위 생존시간: S(t) ≤ s0 인 가장 이른 시점과 그 신뢰구간. 도달하지 않으면 null. */
function svQuantile(km, s0) {
  var eps = 1e-12, i, out = { t: null, lo: null, hi: null };
  for (i = 0; i < km.tt.length; i++) if (km.S[i] <= s0 + eps) { out.t = km.tt[i]; break; }
  for (i = 0; i < km.tt.length; i++) if (km.lo[i] <= s0 + eps) { out.lo = km.tt[i]; break; }
  for (i = 0; i < km.tt.length; i++) if (km.hi[i] <= s0 + eps) { out.hi = km.tt[i]; break; }
  return out;
}

/** 제한 평균 생존시간: 0 부터 가장 큰 관측 시간 tau 까지 S(t) 아래 넓이 */
function svRestrictedMean(km, tau) {
  var area = 0, prevT = 0, prevS = 1, i;
  for (i = 0; i < km.tt.length; i++) {
    if (km.tt[i] > tau) break;
    area += prevS * (km.tt[i] - prevT);
    prevT = km.tt[i]; prevS = km.S[i];
  }
  return area + prevS * (tau - prevT);
}

/* ---------- 집단 비교 ---------- */

/** 작은 연립방정식 풀이 (부분 피벗). 특이하면 null */
function svSolve_(A, b) {
  var n = b.length, M = A.map(function (r, i) { return r.slice().concat([b[i]]); }), i, j, k;
  for (i = 0; i < n; i++) {
    var p = i;
    for (j = i + 1; j < n; j++) if (Math.abs(M[j][i]) > Math.abs(M[p][i])) p = j;
    if (Math.abs(M[p][i]) < 1e-12) return null;
    var tmp = M[i]; M[i] = M[p]; M[p] = tmp;
    for (j = i + 1; j < n; j++) {
      var f = M[j][i] / M[i][i];
      for (k = i; k <= n; k++) M[j][k] -= f * M[i][k];
    }
  }
  var x = new Array(n);
  for (i = n - 1; i >= 0; i--) {
    var s = M[i][n];
    for (j = i + 1; j < n; j++) s -= M[i][j] * x[j];
    x[i] = s / M[i][i];
  }
  return x;
}

var SV_WEIGHTS_ = [
  { key: 'logrank', label: 'Log-rank (Mantel–Cox)', w: function () { return 1; } },
  { key: 'breslow', label: 'Breslow (일반화 Wilcoxon, Gehan)', w: function (n) { return n; } },
  { key: 'tarone', label: 'Tarone–Ware', w: function (n) { return Math.sqrt(n); } }
];

/**
 * 집단 비교 검정. time, event, grp(0..K−1 번호), active(비교할 집단 번호 배열, 2개 이상).
 * 반환: [{ key, label, chi2, df, p }] (가중치 3종), 계산 못 하면 chi2 = null
 */
function svCompare(time, event, grp, active, pchisq) {
  var K = active.length, pos = {}, i, a, b;
  active.forEach(function (g, q) { pos[g] = q; });
  var idx = [];
  for (i = 0; i < time.length; i++) if (pos[grp[i]] !== undefined) idx.push(i);
  idx.sort(function (x, y) { return time[x] - time[y]; });
  var N = idx.length, nk = new Array(K).fill(0);
  idx.forEach(function (r) { nk[pos[grp[r]]]++; });
  var W = SV_WEIGHTS_.length, U = [], V = [];
  for (var w = 0; w < W; w++) {
    U.push(new Array(K).fill(0));
    V.push(Array.apply(null, Array(K)).map(function () { return new Array(K).fill(0); }));
  }
  var p = 0, alive = N;
  while (p < N) {
    var t = time[idx[p]], dk = new Array(K).fill(0), d = 0, q = p, cnt = new Array(K).fill(0);
    while (q < N && time[idx[q]] === t) {
      var gq = pos[grp[idx[q]]];
      cnt[gq]++;
      if (event[idx[q]]) { dk[gq]++; d++; }
      q++;
    }
    if (d > 0) {
      for (var wi = 0; wi < W; wi++) {
        var wt = SV_WEIGHTS_[wi].w(alive);
        for (a = 0; a < K; a++) U[wi][a] += wt * (dk[a] - d * nk[a] / alive);
        if (alive > 1) {
          var c = wt * wt * d * (alive - d) / (alive - 1);
          for (a = 0; a < K; a++) {
            for (b = 0; b < K; b++) {
              V[wi][a][b] += c * (nk[a] / alive) * ((a === b ? 1 : 0) - nk[b] / alive);
            }
          }
        }
      }
    }
    for (a = 0; a < K; a++) nk[a] -= cnt[a];
    alive -= q - p;
    p = q;
  }
  // 분산이 0 인 집단(사건 시점에 위험집단에 한 번도 없었던 집단)은 빼고, 마지막 하나를 기준으로 K−1 개를 쓴다
  var use = [];
  for (a = 0; a < K; a++) if (V[0][a][a] > 1e-12) use.push(a);
  var dfUse = use.length - 1;
  return SV_WEIGHTS_.map(function (sw, wi) {
    var out = { key: sw.key, label: sw.label, chi2: null, df: dfUse, p: null };
    if (dfUse < 1) return out;
    var sub = use.slice(0, dfUse);
    var A = sub.map(function (r) { return sub.map(function (c2) { return V[wi][r][c2]; }); });
    var x = svSolve_(A, sub.map(function (r) { return U[wi][r]; }));
    if (!x) return out;
    var chi2 = 0;
    sub.forEach(function (r, q2) { chi2 += U[wi][r] * x[q2]; });
    out.chi2 = Math.max(0, chi2);
    out.p = pchisq(out.chi2, dfUse, false);
    return out;
  });
}

/* ---------- 입력 처리 ---------- */

function svIsEvent_(v, ev) {
  var evNum = Number(ev);
  if (typeof v === 'number' && ev !== '' && isFinite(evNum)) return v === evNum;
  return String(v).trim() === String(ev).trim();
}

/** 표시용 "보기 좋은" 시점: 가장 큰 시간의 1/5 근처의 1, 2, 2.5, 5 × 10^k 간격 */
function svAutoTimes_(tmax) {
  if (!(tmax > 0)) return [];
  var raw = tmax / 5, mag = Math.pow(10, Math.floor(Math.log(raw) / Math.LN10)), step = mag;
  [1, 2, 2.5, 5, 10].some(function (m) { step = m * mag; return step >= raw - 1e-12; });
  var out = [];
  for (var k = 1; k * step <= tmax + 1e-9 && out.length < 8; k++) out.push(Number((k * step).toPrecision(12)));
  return out;
}

function svParseTimes_(text) {
  var parts = String(text).split(/[\s,;]+/).filter(function (s) { return s !== ''; });
  var out = [];
  parts.forEach(function (s) {
    var x = Number(s);
    if (!isFinite(x) || x < 0) throw new Error('"생존율을 볼 시점"에는 0 이상의 숫자를 쉼표로 구분해 쓰세요 (예: 12, 24, 36). "' + s + '"을(를) 읽을 수 없습니다.');
    if (out.indexOf(x) < 0) out.push(x);
  });
  if (out.length > SV_MAX_TIMES_) throw new Error('생존율을 볼 시점은 ' + SV_MAX_TIMES_ + '개까지 쓸 수 있습니다.');
  return out.sort(function (a, b) { return a - b; });
}

function svPct_(a, b) { return b > 0 ? 100 * a / b : null; }

function svCurveSeries_(label, km) {
  var pts = [[0, 1]];
  var T = km.tt.length, pickIdx = [], i;
  if (T <= SV_CURVE_CAP_) { for (i = 0; i < T; i++) pickIdx.push(i); }
  else {
    for (i = 0; i < SV_CURVE_CAP_; i++) pickIdx.push(Math.round(i * (T - 1) / (SV_CURVE_CAP_ - 1)));
  }
  var prevS = 1;
  pickIdx.forEach(function (k) {
    pts.push([km.tt[k], prevS]);
    pts.push([km.tt[k], km.S[k]]);
    prevS = km.S[k];
  });
  if (km.tmax > pts[pts.length - 1][0]) pts.push([km.tmax, prevS]);
  return { label: label, points: pts };
}

/* ---------- 실행 ---------- */

function runKm(dataset, request) {
  var d = svDeps_();
  var opt = request.options || {};
  var timeVars = d.getVars(dataset, request.time, { numeric: true });
  var eventVars = d.getVars(dataset, request.event);
  if (timeVars.length !== 1) throw new Error('시간 변수는 하나만 선택하세요.');
  if (eventVars.length !== 1) throw new Error('사건 변수는 하나만 선택하세요.');
  var groupVars = request.group && request.group.length ? d.getVars(dataset, request.group) : [];
  if (groupVars.length > 1) throw new Error('집단변수는 하나만 선택하세요.');
  var tv = timeVars[0], ev = eventVars[0], gv = groupVars[0] || null;
  if (tv.name === ev.name || (gv && (gv.name === tv.name || gv.name === ev.name))) throw new Error('시간 변수, 사건 변수, 집단변수는 서로 다른 변수여야 합니다.');

  var evText = opt.eventValue === undefined || opt.eventValue === null || String(opt.eventValue).trim() === '' ? '1' : String(opt.eventValue).trim();
  var conf = opt.conf === undefined || opt.conf === null || opt.conf === '' ? 95 : Number(opt.conf);
  if (!isFinite(conf) || conf <= 50 || conf >= 100) throw new Error('"신뢰수준(%)"은(는) 50보다 크고 100보다 작아야 합니다.');
  var z = d.qnorm(1 - (1 - conf / 100) / 2);

  // 사례 모으기 (시간·사건·집단이 모두 있는 사례)
  var time = [], event = [], gkey = [], seenEvVals = {}, evOrder = [], neg = 0, i;
  var tcol = dataset.data[tv.name], ecol = dataset.data[ev.name], gcol = gv ? dataset.data[gv.name] : null;
  for (i = 0; i < dataset.n; i++) {
    if (tcol[i] === null || ecol[i] === null || (gcol && gcol[i] === null)) continue;
    if (tcol[i] < 0) { neg++; continue; }
    time.push(tcol[i]);
    var isE = svIsEvent_(ecol[i], evText);
    event.push(isE ? 1 : 0);
    var ek = String(ecol[i]);
    if (!(ek in seenEvVals)) { seenEvVals[ek] = 0; evOrder.push(ecol[i]); }
    seenEvVals[ek]++;
    gkey.push(gcol ? gcol[i] : 0);
  }
  if (neg > 0) throw new Error('시간 변수 "' + svName_(tv) + '"에 음수가 ' + neg + '건 있습니다. 생존 시간은 0 이상이어야 합니다.');
  var used = time.length;
  if (used < 2) throw new Error('시간·사건' + (gv ? '·집단' : '') + ' 변수가 모두 있는 사례가 ' + used + '건입니다 (2건 이상 필요).');
  var totalEvents = event.reduce(function (a, b) { return a + b; }, 0);
  if (totalEvents === 0) throw new Error('사건이 한 건도 없습니다. "사건으로 볼 값"(' + evText + ')이 사건 변수 "' + svName_(ev) + '"에 있는지 확인하세요. 변수에 있는 값: ' + evOrder.slice(0, 10).map(String).join(', ') + '.');

  // 집단 나누기
  var groups = [], grp = new Array(used).fill(0);
  if (gv) {
    var tally = {}, order = [];
    for (i = 0; i < used; i++) {
      var key = typeof gkey[i] + ':' + gkey[i];
      if (!(key in tally)) { tally[key] = { value: gkey[i], rows: [] }; order.push(key); }
      tally[key].rows.push(i);
    }
    var gl = order.map(function (k) { return tally[k]; }).sort(d.compareValues_);
    if (gl.length < 2) throw new Error('집단변수 "' + svName_(gv) + '"의 집단이 ' + gl.length + '개입니다. 집단이 2개 이상이어야 합니다 (집단을 나누지 않으려면 집단변수 칸을 비우세요).');
    if (gl.length > SV_MAX_GROUPS_) throw new Error('집단변수 "' + svName_(gv) + '"의 집단이 ' + gl.length + '개입니다 (최대 ' + SV_MAX_GROUPS_ + '개). 범주가 적은 변수를 선택하세요.');
    gl.forEach(function (g, q) {
      var lab = gv.valueLabels[String(g.value)];
      groups.push({ label: lab !== undefined ? lab : String(g.value), rows: g.rows });
      g.rows.forEach(function (r) { grp[r] = q; });
    });
  } else {
    var all = []; for (i = 0; i < used; i++) all.push(i);
    groups.push({ label: '전체', rows: all });
  }
  groups.forEach(function (g) {
    g.time = g.rows.map(function (r) { return time[r]; });
    g.event = g.rows.map(function (r) { return event[r]; });
    g.km = svKm(g.time, g.event, z);
  });

  var warnings = [];
  if (dataset.n > used) warnings.push('시간·사건' + (gv ? '·집단' : '') + ' 변수 중 하나라도 비어 있는 ' + (dataset.n - used) + '건을 제외했습니다. 분석에 쓴 사례는 ' + used + '건입니다.');
  var noEvent = groups.filter(function (g) { return g.km.events === 0; });
  if (noEvent.length) warnings.push('사건이 한 건도 없는 집단이 있습니다 (' + noEvent.map(function (g) { return g.label; }).join(', ') + '). 그 집단은 생존율이 계속 100%이고 중앙값을 구할 수 없습니다.');
  if (totalEvents < 10) warnings.push('사건이 ' + totalEvents + '건뿐이라 생존율의 신뢰구간이 매우 넓고 검정의 근사가 부정확할 수 있습니다 (사건 10건 이상을 권장하는 기준이 흔합니다).');
  var small = groups.filter(function (g) { return g.km.events > 0 && g.km.events < 5; });
  if (gv && small.length) warnings.push('사건이 5건 미만인 집단이 있습니다 (' + small.map(function (g) { return g.label; }).join(', ') + '). 그 집단의 추정은 불안정합니다.');
  var distinctEv = evOrder.length;
  var evNote = '사건 변수 "' + svName_(ev) + '": 값이 ' + evText + '이면 사건, 그 밖의 값은 모두 중도절단으로 처리했다 (변수의 값: ' + evOrder.slice(0, 8).map(function (v) { var l = ev.valueLabels[String(v)]; return String(v) + (l !== undefined ? '(' + l + ')' : '') + ' ' + seenEvVals[String(v)] + '건'; }).join(', ') + (distinctEv > 8 ? ', …' : '') + ').';
  if (distinctEv > 2) warnings.push('사건 변수에 값이 ' + distinctEv + '가지 있습니다. ' + evText + '인 사례만 사건이고 나머지는 모두 중도절단으로 처리했습니다. 의도와 다르면 "사건으로 볼 값"을 확인하세요.');

  var blocks = [];

  // 처리 요약
  var sumRows = groups.map(function (g) {
    var tau = g.km.tmax;
    return { label: g.label, n: g.km.n, events: g.km.events, censored: g.km.censored, pct: svPct_(g.km.censored, g.km.n), mean: svRestrictedMean(g.km, tau), tau: tau };
  });
  if (groups.length > 1) sumRows.push({ _kind: 'total', label: '전체', n: used, events: totalEvents, censored: used - totalEvents, pct: svPct_(used - totalEvents, used), mean: null, tau: null });
  blocks.push({
    type: 'table', caption: '처리 요약',
    columns: [
      { key: 'label', label: gv ? svName_(gv) : '' }, { key: 'n', label: 'N', format: { type: 'int' } },
      { key: 'events', label: '사건 수', format: { type: 'int' } }, { key: 'censored', label: '중도절단 수', format: { type: 'int' } },
      { key: 'pct', label: '중도절단 %', format: { type: 'fixed', digits: 1 } },
      { key: 'mean', label: '평균 생존시간 (제한)', format: SV_F3_ }, { key: 'tau', label: '제한 시점', format: SV_T_ }
    ],
    rows: sumRows,
    notes: [evNote, '평균 생존시간은 곡선 아래 넓이를 그 집단에서 가장 큰 관측 시간(제한 시점)까지 적분한 값이다. 마지막 관측이 중도절단이면 실제 평균은 더 클 수 있어 집단 사이 비교에는 같은 제한 시점을 쓰는 것이 좋다. 중앙값을 우선 보세요.']
  });

  // 시점별 생존율
  var tList = opt.times === undefined || opt.times === null || String(opt.times).trim() === '' ? null : svParseTimes_(opt.times);
  var tmaxAll = Math.max.apply(null, groups.map(function (g) { return g.km.tmax; }));
  var times = tList || svAutoTimes_(tmaxAll);
  if (times.length) {
    var atRows = [], beyond = false;
    groups.forEach(function (g) {
      times.forEach(function (t, q) {
        var r = svSurvivalAt(g.km, t);
        if (r.s === null) beyond = true;
        atRows.push({ label: q === 0 ? g.label : '', t: t, nr: svAtRisk_(g.km, t), s: r.s, se: r.se, lo: r.lo, hi: r.hi });
      });
    });
    blocks.push({
      type: 'table', caption: '시점별 생존율 (Kaplan–Meier)',
      columns: [
        { key: 'label', label: gv ? svName_(gv) : '' }, { key: 't', label: '시점', format: SV_T_ }, { key: 'nr', label: '위험집단 수', format: { type: 'int' } },
        { key: 's', label: '생존율', format: SV_F3_ }, { key: 'se', label: '표준오차 (Greenwood)', format: SV_F3_ },
        { key: 'lo', label: conf + '% 신뢰구간 하한', format: SV_F3_ }, { key: 'hi', label: conf + '% 신뢰구간 상한', format: SV_F3_ }
      ],
      rows: atRows,
      notes: ['생존율은 그 시점까지 사건이 없을 확률의 추정값이다. 신뢰구간은 log-log 변환(생존율이 0과 1 사이에 있도록)을 쓴다. 위험집단 수는 그 시점에 시간이 그 시점 이상인 사례 수이다.']
        .concat(tList ? [] : ['시점을 따로 지정하지 않아 가장 큰 관측 시간을 5등분하는 눈금 근처의 시점을 자동으로 골랐다.'])
        .concat(beyond ? ['빈 칸은 그 집단의 마지막 관측 시점 뒤라 추정하지 않은 값이다. 첫 사건 전의 생존율은 1이고 신뢰구간이 없다.'] : ['첫 사건 전의 생존율은 1이고 신뢰구간이 없어 신뢰구간 칸이 비어 있다.'])
    });
  }

  // 백분위 생존시간
  var qRows = [], reachedAll = true;
  groups.forEach(function (g) {
    [[0.5, '중앙값 (50%)'], [0.75, '75% (S ≤ 25%)'], [0.25, '25% (S ≤ 75%)']].forEach(function (pq, q) {
      var r = svQuantile(g.km, 1 - pq[0]);
      var nr = function (v) { return v === null ? '도달하지 않음' : v; };
      if (r.t === null) reachedAll = false;
      qRows.push({ label: q === 0 ? g.label : '', pct: pq[1], t: nr(r.t), lo: r.t === null ? '' : nr(r.lo), hi: r.t === null ? '' : nr(r.hi) });
    });
  });
  // 보기 좋게 25%, 50%, 75% 순서로
  var qOrdered = [];
  for (var gi = 0; gi < groups.length; gi++) {
    var base = gi * 3;
    var trio = [qRows[base + 2], qRows[base], qRows[base + 1]];
    trio[0].label = groups[gi].label; trio[1].label = ''; trio[2].label = '';
    qOrdered = qOrdered.concat(trio);
  }
  blocks.push({
    type: 'table', caption: '생존시간 백분위수',
    columns: [
      { key: 'label', label: gv ? svName_(gv) : '' }, { key: 'pct', label: '백분위수 (생존율이 처음 이 이하가 되는 시점)' }, { key: 't', label: '추정', format: SV_T_ },
      { key: 'lo', label: conf + '% 신뢰구간 하한', format: SV_T_ }, { key: 'hi', label: conf + '% 신뢰구간 상한', format: SV_T_ }
    ],
    rows: qOrdered,
    notes: ['중앙값은 생존율이 처음 50% 이하로 내려가는 시점이다 (사례의 절반이 사건을 겪는 시점). 곡선이 그 수준까지 내려가지 않으면 "도달하지 않음"으로 표시한다. 신뢰구간은 시점별 신뢰구간 곡선이 그 수준을 지나는 시점으로 구한다 (구간의 상한이 비어 있으면 상한 곡선도 도달하지 않은 것이다).']
      .concat(reachedAll ? [] : ['중도절단이 많거나 추적 기간이 짧으면 중앙값에 도달하지 못합니다. 이럴 때는 특정 시점의 생존율로 비교하세요.'])
  });

  // 집단 비교
  if (gv && opt.tests !== false) {
    var allIdx = groups.map(function (g, q) { return q; });
    var cmp = svCompare(time, event, grp, allIdx, d.pchisq);
    var tRows = cmp.map(function (r) { return { label: r.label, chi2: r.chi2, df: r.df, p: r.p }; });
    var tNotes = ['세 검정 모두 영가설은 "집단 사이에 생존곡선(위험)의 차이가 없다"이다. 대표 검정은 log-rank이고 전 기간의 차이에 고르게 민감하다. Breslow와 Tarone–Ware는 초기 사건에 더 큰 가중을 둔다 (곡선이 교차하거나 초기에 차이가 큰 자료에서 결과가 달라질 수 있다).',
      'log-rank는 위험비가 시간에 따라 일정(비례위험)할 때 검정력이 가장 높다. 곡선이 교차하면 검정력이 크게 떨어질 수 있으니 그림을 함께 보세요. 검정은 중도절단이 집단과 무관(비정보적)하다고 가정한다.'];
    if (cmp[0].df !== groups.length - 1) tNotes.push('사건 시점에 한 번도 위험집단에 들지 않은 집단이 있어 자유도를 ' + cmp[0].df + '로 줄였다.');
    blocks.push({
      type: 'table', caption: '집단 비교 검정 (전체 ' + groups.length + '개 집단)',
      columns: [{ key: 'label', label: '검정' }, { key: 'chi2', label: '카이제곱', format: SV_F3_ }, { key: 'df', label: '자유도', format: { type: 'int' } }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } }],
      rows: tRows.map(function (r, q) { return Object.assign({ _bold: q === 0 }, r); }),
      notes: tNotes
    });
    if (groups.length >= 3 && opt.pairs !== false) {
      var pRows = [], m = groups.length * (groups.length - 1) / 2;
      for (var a = 0; a < groups.length; a++) {
        for (var b = a + 1; b < groups.length; b++) {
          var rr = svCompare(time, event, grp, [a, b], d.pchisq)[0];
          pRows.push({ gi: groups[a].label, gj: groups[b].label, chi2: rr.chi2, df: rr.df, p: rr.p, pAdj: rr.p === null ? null : Math.min(1, rr.p * m) });
        }
      }
      blocks.push({
        type: 'table', caption: '집단 쌍별 log-rank 검정 (Bonferroni 보정)',
        columns: [
          { key: 'gi', label: '집단 A' }, { key: 'gj', label: '집단 B' }, { key: 'chi2', label: '카이제곱', format: SV_F3_ },
          { key: 'df', label: '자유도', format: { type: 'int' } }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } }, { key: 'pAdj', label: '보정 유의확률 (×' + m + ')', format: { type: 'pvalue' } }
        ],
        rows: pRows,
        notes: ['두 집단씩 따로 log-rank 검정을 하고 유의확률에 비교 횟수(' + m + ')를 곱했다 (1을 넘으면 1). 전체 검정이 유의할 때 어느 쌍이 다른지 보는 용도이다.']
      });
    }
  }

  // 전체 사건 시점 표
  if (opt.table === true) {
    groups.forEach(function (g) {
      var km = g.km, rows = [], cap = Math.min(km.tt.length, SV_TABLE_CAP_);
      for (var k = 0; k < cap; k++) rows.push({ t: km.tt[k], nr: km.nr[k], dd: km.dd[k], s: km.S[k], se: km.se[k], lo: km.S[k] > 0 ? km.lo[k] : null, hi: km.S[k] > 0 ? km.hi[k] : null });
      blocks.push({
        type: 'table', caption: '사건 시점별 생존표' + (gv ? ' — ' + g.label : ''),
        columns: [
          { key: 't', label: '시점', format: SV_T_ }, { key: 'nr', label: '위험집단 수', format: { type: 'int' } }, { key: 'dd', label: '사건 수', format: { type: 'int' } },
          { key: 's', label: '생존율', format: SV_F3_ }, { key: 'se', label: '표준오차', format: SV_F3_ },
          { key: 'lo', label: conf + '% 신뢰구간 하한', format: SV_F3_ }, { key: 'hi', label: conf + '% 신뢰구간 상한', format: SV_F3_ }
        ],
        rows: rows,
        notes: ['사건이 일어난 시점만 나열했다 (중도절단만 있는 시점은 위험집단 수의 감소로 반영된다).']
          .concat(km.tt.length > cap ? ['사건 시점이 ' + km.tt.length + '개라 처음 ' + cap + '개만 표시했다.'] : [])
      });
    });
  }

  // 생존곡선
  if (opt.plot !== false) {
    var series = groups.map(function (g) { return svCurveSeries_(g.label, g.km); });
    blocks.push(d.lineChartBlock('Kaplan–Meier 생존곡선', '시간 (' + svName_(tv) + ')', '생존율', series));
    var thinned = groups.some(function (g) { return g.km.tt.length > SV_CURVE_CAP_; });
    if (thinned) blocks.push({ type: 'text', text: '사건 시점이 많아 그림은 시점을 줄여 그렸다 (표의 값은 줄이지 않았다).' });
  }

  return { analysis: 'km', title: 'Kaplan–Meier 생존분석' + (gv ? ' (집단: ' + svName_(gv) + ')' : ''), blocks: blocks, warnings: warnings };
}

if (typeof module !== 'undefined') {
  module.exports = { svKm: svKm, svSurvivalAt: svSurvivalAt, svQuantile: svQuantile, svRestrictedMean: svRestrictedMean, svCompare: svCompare, runKm: runKm };
}
