/**
 * 의사결정나무(CART). 순수 함수: (dataset, request) → 결과 객체.
 *
 * - 종속변수가 연속(scale)이면 회귀 나무(SSE), 명목·서열이면 분류 나무(지니 불순도). 이진 분할만 한다.
 * - 불순도: 분류 1 − Σp², 회귀 SSE/n. 노드의 위험 R(t): 분류는 오분류 건수, 회귀는 SSE. 분할 이득 Δ = n·i(t) − nL·i(L) − nR·i(R).
 * - 분할 후보: 숫자 변수는 노드 안 서로 다른 값의 중간값(≤ 왼쪽). 범주형은 (회귀, 2분류) 응답 통계(평균, 둘째 범주 비율)로 범주를 정렬해 앞에서부터 자른 집합,
 *   3분류 이상은 첫 범주를 포함하는 모든 부분집합(범주 12개까지). 양쪽이 최소 잎 크기 이상이어야 한다.
 * - 최고 이득은 변수 순서(숫자 변수들, 범주형들), 후보 순서대로 보며 a − b > 1e-10·max(1, |b|)일 때만 바꾼다 (동점은 앞선 후보).
 * - 중지: n < minSplit, 깊이 ≥ maxDepth, 불순도 ≤ 1e-12, 유효한 분할 없음 또는 이득 ≤ 1e-12. 노드 번호는 전위 순회.
 * - 비용복잡도 가지치기: 약한 고리 g(t) = (R(t) − R(T_t)) / (잎 수 − 1)가 최소인 노드들을 한꺼번에 접어 부분나무열 S_0(전체) … S_m(루트)을 만든다.
 * - 교차검증으로 부분나무를 고른다 (시드로 재현): 각 폴드에서 나무를 키워 β_k = √(α_k·α_{k+1})에 해당하는 부분나무로 검증 사례를 예측한다.
 *   선택: cv1se(기본) = 최소 CV 오차 + 1 표준오차 이하 중 가장 작은 나무, cv = 최소 CV 오차, none = 가지치기 안 함.
 * - 변수 중요도는 최종 나무에서 그 변수로 분할한 이득의 합 (대리 분할은 쓰지 않는다).
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function cartDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var charts = require('./Charts.js');
    return {
      getVars: data.getVars, compareValues_: data.compareValues_, mulberry32_: require('./Stats_Mediation.js').mulberry32_,
      lineChartBlock: charts.lineChartBlock, treeChartBlock: charts.treeChartBlock
    };
  }
  return { getVars: getVars, compareValues_: compareValues_, mulberry32_: mulberry32_, lineChartBlock: lineChartBlock, treeChartBlock: treeChartBlock };
}

var CART_MAX_N_ = 5000;
var CART_MAX_CLASSES_ = 10;
var CART_MAX_CATS_ = 12;
var CART_MAX_DEPTH_ = 8;
var CART_DIAGRAM_MAX_LEAVES_ = 40;
function cartName_(v) { return v.label ? v.label : v.name; }
function cartFmt_(x) { return String(Number(Number(x).toPrecision(6))); }
function cartBetter_(a, b) { return a - b > 1e-10 * Math.max(1, Math.abs(b)); }

/** 노드 통계 (idx: 사례 인덱스 배열). 분류는 counts, 회귀는 mean */
function cartStats_(T, idx) {
  var n = idx.length, i;
  if (T.mode === 'class') {
    var counts = new Array(T.C).fill(0);
    for (i = 0; i < n; i++) counts[T.y[idx[i]]]++;
    var mx = 0, pred = 0, sq = 0;
    for (i = 0; i < T.C; i++) { sq += counts[i] * counts[i]; if (counts[i] > mx) { mx = counts[i]; pred = i; } }
    return { n: n, counts: counts, imp: 1 - sq / (n * n), pred: pred, err: n - mx, nimp: n - sq / n };
  }
  var s = 0;
  for (i = 0; i < n; i++) s += T.y[idx[i]];
  var m = s / n, sse = 0;
  for (i = 0; i < n; i++) sse += (T.y[idx[i]] - m) * (T.y[idx[i]] - m);
  return { n: n, mean: m, imp: sse / n, pred: m, err: sse, nimp: sse };
}

/** 최고 분할을 찾는다 (없으면 null). parentNimp: 부모의 n·i */
function cartBestSplit_(T, idx, parentNimp) {
  var best = null, n = idx.length, minLeaf = T.params.minLeaf, j, t, c, k;
  var consider = function (gain, cand) {
    if (best === null || cartBetter_(gain, best.gain)) { cand.gain = gain; best = cand; }
  };
  // 회귀는 노드 평균으로 중심화해 접두합의 소거 오차를 줄인다
  var center = 0;
  if (T.mode === 'reg') { for (t = 0; t < n; t++) center += T.y[idx[t]]; center /= n; }
  for (j = 0; j < T.P.length; j++) {
    var col = T.P[j].values;
    if (T.P[j].type === 'num') {
      var order = idx.slice().sort(function (a, b) { return col[a] - col[b] || a - b; });
      var lc = null, rc = null, ls = 0, lq = 0, S = 0, Q = 0;
      if (T.mode === 'class') { lc = new Array(T.C).fill(0); rc = new Array(T.C).fill(0); for (t = 0; t < n; t++) rc[T.y[order[t]]]++; }
      else { for (t = 0; t < n; t++) { var yc = T.y[order[t]] - center; S += yc; Q += yc * yc; } }
      for (t = 1; t < n; t++) {
        var prev = order[t - 1];
        if (T.mode === 'class') { lc[T.y[prev]]++; rc[T.y[prev]]--; }
        else { var yp = T.y[prev] - center; ls += yp; lq += yp * yp; }
        if (!(col[order[t - 1]] < col[order[t]])) continue;
        if (t < minLeaf || n - t < minLeaf) continue;
        var nl = t, nr = n - t, imp;
        if (T.mode === 'class') {
          var sl = 0, sr = 0;
          for (c = 0; c < T.C; c++) { sl += lc[c] * lc[c]; sr += rc[c] * rc[c]; }
          imp = (nl - sl / nl) + (nr - sr / nr);
        } else {
          imp = (lq - ls * ls / nl) + ((Q - lq) - (S - ls) * (S - ls) / nr);
        }
        consider(parentNimp - imp, { j: j, kind: 'num', arg: (col[order[t - 1]] + col[order[t]]) / 2 });
      }
    } else {
      // 노드에 있는 범주별 집계
      var agg = {}, cats = [];
      for (t = 0; t < n; t++) {
        var code = col[idx[t]];
        if (!agg[code]) { agg[code] = { n: 0, s: 0, q: 0, counts: T.mode === 'class' ? new Array(T.C).fill(0) : null }; cats.push(code); }
        var a = agg[code];
        a.n++;
        if (T.mode === 'class') a.counts[T.y[idx[t]]]++;
        else { var yy = T.y[idx[t]] - center; a.s += yy; a.q += yy * yy; }
      }
      cats.sort(function (a2, b2) { return a2 - b2; });
      var m = cats.length;
      if (m < 2) continue;
      var subsets = [];
      if (T.mode === 'reg' || T.C === 2) {
        var stat = {};
        cats.forEach(function (cc) { stat[cc] = T.mode === 'reg' ? agg[cc].s / agg[cc].n + center : agg[cc].counts[1] / agg[cc].n; });
        var ord = cats.slice().sort(function (a2, b2) { return stat[a2] - stat[b2] || a2 - b2; });
        for (t = 1; t < m; t++) subsets.push(ord.slice(0, t));
      } else {
        var total = Math.pow(2, m - 1) - 1;
        for (var mask = 0; mask < total; mask++) {
          var sub = [cats[0]];
          for (k = 0; k < m - 1; k++) if ((mask >> k) & 1) sub.push(cats[k + 1]);
          subsets.push(sub);
        }
      }
      var allN = n, allS = 0, allQ = 0, allC = T.mode === 'class' ? new Array(T.C).fill(0) : null;
      cats.forEach(function (cc) {
        allS += agg[cc].s; allQ += agg[cc].q;
        if (allC) for (var z = 0; z < T.C; z++) allC[z] += agg[cc].counts[z];
      });
      subsets.forEach(function (sub) {
        var inSet = {};
        sub.forEach(function (cc) { inSet[cc] = true; });
        var nl2 = 0, s2 = 0, q2 = 0, lc2 = allC ? new Array(T.C).fill(0) : null;
        sub.forEach(function (cc) {
          nl2 += agg[cc].n; s2 += agg[cc].s; q2 += agg[cc].q;
          if (lc2) for (var z = 0; z < T.C; z++) lc2[z] += agg[cc].counts[z];
        });
        var nr2 = allN - nl2;
        if (nl2 < minLeaf || nr2 < minLeaf) return;
        var imp2;
        if (allC) {
          var a1 = 0, a2 = 0;
          for (var z2 = 0; z2 < T.C; z2++) { a1 += lc2[z2] * lc2[z2]; a2 += (allC[z2] - lc2[z2]) * (allC[z2] - lc2[z2]); }
          imp2 = (nl2 - a1 / nl2) + (nr2 - a2 / nr2);
        } else {
          imp2 = (q2 - s2 * s2 / nl2) + ((allQ - q2) - (allS - s2) * (allS - s2) / nr2);
        }
        consider(parentNimp - imp2, { j: j, kind: 'cat', arg: sub.slice().sort(function (a3, b3) { return a3 - b3; }) });
      });
    }
  }
  return best;
}

function cartGoesLeft_(T, node, i) {
  var s = node.split, v = T.P[s.j].values[i];
  if (s.kind === 'num') return v <= s.arg;
  if (s.argSet[v]) return true;
  if (s.known && !s.known[v]) return s.majorLeft;
  return false;
}

/** 나무 키우기. idx: 학습 사례 인덱스. 반환: { nodes (id-1 로 접근), P, y, mode, C, params } */
function cartBuild(P, y, mode, C, params, idx) {
  var T = { P: P, y: y, mode: mode, C: C, params: params, nodes: [] };
  var grow = function (ids, depth) {
    var st = cartStats_(T, ids);
    var node = { id: T.nodes.length + 1, depth: depth, left: 0, right: 0, split: null, n: st.n, imp: st.imp, pred: st.pred, err: st.err, counts: st.counts || null, mean: st.mean };
    T.nodes.push(node);
    if (st.n < params.minSplit || depth >= params.maxDepth || st.imp <= 1e-12) return node;
    var sp = cartBestSplit_(T, ids, st.nimp);
    if (sp === null || sp.gain <= 1e-12) return node;
    var argSet = {};
    if (sp.kind === 'cat') sp.arg.forEach(function (c) { argSet[c] = true; });
    node.split = { j: sp.j, kind: sp.kind, arg: sp.arg, argSet: argSet, gain: sp.gain, known: null, majorLeft: false };
    var left = [], right = [], known = {};
    ids.forEach(function (i) {
      var isLeft = cartGoesLeft_(T, node, i);
      (isLeft ? left : right).push(i);
      if (sp.kind === 'cat') known[T.P[sp.j].values[i]] = true;
    });
    node.split.majorLeft = left.length >= right.length;
    if (sp.kind === 'cat') node.split.known = known;
    node.left = grow(left, depth + 1).id;
    node.right = grow(right, depth + 1).id;
    return node;
  };
  grow(idx, 0);
  return T;
}

function cartLeaf_(T, i, collapsed) {
  var node = T.nodes[0];
  while (node.split !== null && !collapsed[node.id]) node = T.nodes[(cartGoesLeft_(T, node, i) ? node.left : node.right) - 1];
  return node;
}

function cartSubtree_(T, node, collapsed) {
  if (node.split === null || collapsed[node.id]) return { r: node.err, leaves: 1 };
  var a = cartSubtree_(T, T.nodes[node.left - 1], collapsed), b = cartSubtree_(T, T.nodes[node.right - 1], collapsed);
  return { r: a.r + b.r, leaves: a.leaves + b.leaves };
}

function cartReach_(T, collapsed) {
  var seen = {};
  var rec = function (node) {
    seen[node.id] = true;
    if (node.split !== null && !collapsed[node.id]) { rec(T.nodes[node.left - 1]); rec(T.nodes[node.right - 1]); }
  };
  rec(T.nodes[0]);
  return seen;
}

/** 부분나무열. 반환: { alphas (S_{k-1}→S_k 의 α, 길이 m), sets (접힌 노드 집합, 길이 m+1), leaves (길이 m+1) } */
function cartPruneSequence(T) {
  var collapsed = {}, alphas = [], sets = [{}], leaves = [cartSubtree_(T, T.nodes[0], collapsed).leaves], i;
  for (;;) {
    var reach = cartReach_(T, collapsed), internal = [];
    for (i = 1; i <= T.nodes.length; i++) if (reach[i] && T.nodes[i - 1].split !== null && !collapsed[i]) internal.push(i);
    if (!internal.length) break;
    var g = {}, gmin = Infinity;
    internal.forEach(function (id) {
      var node = T.nodes[id - 1], st = cartSubtree_(T, node, collapsed);
      g[id] = Math.max(0, (node.err - st.r) / (st.leaves - 1));
      if (g[id] < gmin) gmin = g[id];
    });
    internal.forEach(function (id) {
      if (g[id] <= gmin + 1e-12 * Math.max(1, Math.abs(gmin)) && cartReach_(T, collapsed)[id]) collapsed[id] = true;
    });
    alphas.push(gmin);
    var copy = {};
    Object.keys(collapsed).forEach(function (k) { copy[k] = true; });
    sets.push(copy);
    leaves.push(cartSubtree_(T, T.nodes[0], collapsed).leaves);
  }
  return { alphas: alphas, sets: sets, leaves: leaves };
}

/** 대표 복잡도 β_0 = 0, β_k = √(α_k·α_{k+1}), β_m = ∞ */
function cartBetas_(alphas) {
  var m = alphas.length, out = [0], k;
  for (k = 1; k < m; k++) out.push(alphas[k - 1] > 0 ? Math.sqrt(alphas[k - 1] * alphas[k]) : 0);
  out.push(Infinity);
  return out.slice(0, m + 1);
}

/** 교차검증. 반환: { xerror, se, betas, losses(부분나무 × 사례) } */
function cartCrossValidate(P, y, mode, C, params, alphas, seed, nfold) {
  var d = cartDeps_(), n = y.length, i, k, f;
  var rnd = d.mulberry32_(seed), order = [];
  for (i = 0; i < n; i++) order.push(i);
  for (i = n - 1; i > 0; i--) { var j = Math.floor(rnd() * (i + 1)); var tmp = order[i]; order[i] = order[j]; order[j] = tmp; }
  var fold = new Array(n);
  order.forEach(function (idx, pos) { fold[idx] = pos % nfold; });
  var bt = cartBetas_(alphas), losses = bt.map(function () { return new Array(n).fill(0); });
  for (f = 0; f < nfold; f++) {
    var train = [], test = [];
    for (i = 0; i < n; i++) (fold[i] === f ? test : train).push(i);
    var T = cartBuild(P, y, mode, C, params, train), seq = cartPruneSequence(T);
    for (k = 0; k < bt.length; k++) {
      var jj = 0;
      seq.alphas.forEach(function (a) { if (a <= bt[k]) jj++; });
      test.forEach(function (ti) {
        var leaf = cartLeaf_(T, ti, seq.sets[jj]);
        losses[k][ti] = mode === 'class' ? (leaf.pred === y[ti] ? 0 : 1) : (leaf.pred - y[ti]) * (leaf.pred - y[ti]);
      });
    }
  }
  var xerror = [], se = [];
  losses.forEach(function (ls) {
    var s = 0, m, ss = 0;
    for (i = 0; i < n; i++) s += ls[i];
    m = s / n;
    for (i = 0; i < n; i++) ss += (ls[i] - m) * (ls[i] - m);
    xerror.push(m); se.push(Math.sqrt(ss / (n - 1)) / Math.sqrt(n));
  });
  return { xerror: xerror, se: se, betas: bt, losses: losses };
}

/** 부분나무 선택. rule: 'cv' | 'cv1se' */
function cartChoose(xerror, se, rule) {
  var m = xerror.length - 1, kmin = 0, k;
  for (k = 1; k <= m; k++) if (xerror[k] <= xerror[kmin]) kmin = k;           // 동점은 작은 나무(큰 k)
  if (rule === 'cv') return kmin;
  var thr = xerror[kmin] + se[kmin], pick = kmin;
  for (k = 0; k <= m; k++) if (xerror[k] <= thr + 1e-12) pick = k;
  return pick;
}

/** 보이는 나무를 전위 순회로 다시 번호 매긴 노드 목록 (접힌 노드는 잎) */
function cartExport(T, collapsed) {
  var out = [];
  var rec = function (node, parent, side) {
    var nid = out.length + 1;
    var rn = { id: nid, old: node.id, parent: parent, side: side, depth: node.depth, n: node.n, pred: node.pred, imp: node.imp, err: node.err, counts: node.counts, mean: node.mean, left: 0, right: 0, split: null, known: null };
    out.push(rn);
    if (node.split !== null && !collapsed[node.id]) {
      rn.split = { j: node.split.j, kind: node.split.kind, arg: node.split.arg, gain: node.split.gain };
      rn.known = node.split.known;
      rn.left = rec(T.nodes[node.left - 1], nid, 'L');
      rn.right = rec(T.nodes[node.right - 1], nid, 'R');
    }
    return nid;
  };
  rec(T.nodes[0], 0, '');
  return out;
}

function runCart(dataset, request) {
  var d = cartDeps_();
  var dep = d.getVars(dataset, request.dependent);
  if (dep.length !== 1) throw new Error('종속변수는 하나만 선택하세요.');
  var covs = request.covariates && request.covariates.length ? d.getVars(dataset, request.covariates, { numeric: true }) : [];
  var facs = request.factors && request.factors.length ? d.getVars(dataset, request.factors) : [];
  if (covs.length + facs.length === 0) throw new Error('예측변수(숫자 또는 범주형)를 하나 이상 선택하세요.');
  var seen = {};
  seen[dep[0].name] = true;
  covs.concat(facs).forEach(function (v) {
    if (seen[v.name]) throw new Error('변수 "' + v.name + '"를 종속변수와 예측변수에 동시에 쓸 수 없습니다.');
    seen[v.name] = true;
  });
  var opt = request.options || {};
  var num = function (v, def) { return v === undefined || v === null || v === '' ? def : Number(v); };
  var isInt = function (v, lo, hi) { return isFinite(v) && v === Math.floor(v) && v >= lo && v <= hi; };
  var maxDepth = num(opt.maxDepth, 4), minSplit = num(opt.minSplit, 20), minLeaf = num(opt.minLeaf, 7), nfold = num(opt.folds, 10), seed = num(opt.seed, 12345);
  if (!isInt(maxDepth, 1, CART_MAX_DEPTH_)) throw new Error('최대 깊이는 1~' + CART_MAX_DEPTH_ + '의 정수여야 합니다.');
  if (!isInt(minSplit, 2, 100000)) throw new Error('분할 최소 사례 수는 2 이상의 정수여야 합니다.');
  if (!isInt(minLeaf, 1, 100000)) throw new Error('잎 최소 사례 수는 1 이상의 정수여야 합니다.');
  if (!isInt(nfold, 2, 20)) throw new Error('교차검증 폴드 수는 2~20의 정수여야 합니다.');
  if (!isInt(seed, -2147483648, 2147483647)) throw new Error('시드는 정수여야 합니다.');
  var prune = opt.prune === 'cv' || opt.prune === 'none' ? opt.prune : 'cv1se';
  var mode = opt.mode === 'class' || opt.mode === 'reg' ? opt.mode : (dep[0].scale === 'scale' ? 'reg' : 'class');
  if (mode === 'reg') d.getVars(dataset, [dep[0].name], { numeric: true });

  var all = dep.concat(covs, facs), rowsOk = [], i, a;
  for (i = 0; i < dataset.n; i++) {
    var ok = true;
    for (a = 0; a < all.length; a++) if (dataset.data[all[a].name][i] === null) { ok = false; break; }
    if (ok) rowsOk.push(i);
  }
  var n = rowsOk.length, warnings = [];
  if (dataset.n > n) warnings.push('목록별 결측 제거로 ' + (dataset.n - n) + '건을 제외했습니다 (사용한 사례 ' + n + '건).');
  if (n < 10) throw new Error('사용할 수 있는 사례가 10건 이상 있어야 합니다 (현재 ' + n + '건).');
  if (n > CART_MAX_N_) throw new Error('사례가 ' + n + '건입니다. 의사결정나무는 교차검증 때문에 최대 ' + CART_MAX_N_ + '건까지 계산합니다. 사례를 줄이거나 나눠 주세요.');
  var pick = function (name) { return rowsOk.map(function (r) { return dataset.data[name][r]; }); };

  // 종속변수
  var yRaw = pick(dep[0].name), y, C = 0, levels = [];
  if (mode === 'class') {
    var seenL = {};
    yRaw.forEach(function (x) { var key = typeof x + ':' + x; if (!seenL[key]) { seenL[key] = true; levels.push({ value: x, key: key, count: 0 }); } });
    levels.sort(d.compareValues_);
    if (levels.length < 2) throw new Error('종속변수 "' + cartName_(dep[0]) + '"의 값이 한 가지뿐입니다. 분류 나무에는 두 가지 이상이 필요합니다.');
    if (levels.length > CART_MAX_CLASSES_) throw new Error('종속변수 "' + cartName_(dep[0]) + '"의 값이 ' + levels.length + '가지입니다 (분류 나무는 최대 ' + CART_MAX_CLASSES_ + '가지). 숫자 변수이면 척도를 연속으로 바꾸거나 값을 묶어 보세요.');
    var codeOf = {};
    levels.forEach(function (l, k) {
      codeOf[l.key] = k;
      var lab = dep[0].valueLabels[String(l.value)];
      l.label = lab !== undefined ? lab : String(l.value);
    });
    C = levels.length;
    y = yRaw.map(function (x) { var l = levels[codeOf[typeof x + ':' + x]]; l.count++; return codeOf[typeof x + ':' + x]; });
  } else {
    y = yRaw.slice();
  }

  // 예측변수: 숫자(순서가 있는 임계값 분할), 범주형(범주 집합 분할)
  var P = [], catInfo = [];
  covs.forEach(function (v) { P.push({ name: cartName_(v), type: 'num', values: pick(v.name), cats: null }); });
  facs.forEach(function (v) {
    var raw = pick(v.name), cs = [], seenC = {};
    raw.forEach(function (x) { var key = typeof x + ':' + x; if (!seenC[key]) { seenC[key] = true; cs.push({ value: x, key: key }); } });
    cs.sort(d.compareValues_);
    if (cs.length > CART_MAX_CATS_) throw new Error('범주형 예측변수 "' + cartName_(v) + '"의 범주가 ' + cs.length + '개입니다 (최대 ' + CART_MAX_CATS_ + '개). 범주를 묶거나 다른 변수를 선택하세요.');
    var idxOf = {};
    cs.forEach(function (c, k) {
      idxOf[c.key] = k;
      var lab = v.valueLabels[String(c.value)];
      c.label = lab !== undefined ? lab : String(c.value);
    });
    P.push({ name: cartName_(v), type: 'cat', values: raw.map(function (x) { return idxOf[typeof x + ':' + x]; }), cats: cs.map(function (c) { return c.label; }) });
  });
  var params = { maxDepth: maxDepth, minSplit: minSplit, minLeaf: minLeaf };
  var idxAll = [];
  for (i = 0; i < n; i++) idxAll.push(i);
  var T = cartBuild(P, y, mode, C, params, idxAll);
  var seq = cartPruneSequence(T);
  var m = seq.alphas.length, cv = null, k = 0;
  if (prune !== 'none' && m > 0) {
    cv = cartCrossValidate(P, y, mode, C, params, seq.alphas, seed, Math.min(nfold, n));
    k = cartChoose(cv.xerror, cv.se, prune);
  }
  var collapsed = seq.sets[k];
  var nodes = cartExport(T, collapsed);
  var leafCount = nodes.filter(function (nd) { return nd.split === null; }).length;
  if (leafCount === 1) warnings.push(m === 0 ? '유효한 분할이 없어 나무가 루트 하나뿐입니다 (최소 사례 수와 깊이 설정을 확인하세요).' : '교차검증으로 가지치기한 결과 루트 하나만 남았습니다. 예측변수가 종속변수를 설명하지 못하거나 가지치기 방식을 "최소 CV 오차"로 바꿔 볼 수 있습니다.');

  // 규칙 문구
  var catLabels = function (nd, side) {
    var P0 = P[nd.split.j], inSet = {};
    nd.split.arg.forEach(function (c) { inSet[c] = true; });
    var codes = Object.keys(nd.known || {}).map(Number).sort(function (a2, b2) { return a2 - b2; });
    var sel = codes.filter(function (c) { return side === 'L' ? inSet[c] : !inSet[c]; });
    return '{' + sel.map(function (c) { return P0.cats[c]; }).join(', ') + '}';
  };
  var edgeText = function (nd) {
    if (nd.parent === 0) return '';
    var p = nodes[nd.parent - 1], s = p.split;
    if (s.kind === 'num') return nd.side === 'L' ? P[s.j].name + ' ≤ ' + cartFmt_(s.arg) : P[s.j].name + ' > ' + cartFmt_(s.arg);
    return P[s.j].name + ' ∈ ' + catLabels(p, nd.side);
  };
  var predText = function (nd) { return mode === 'class' ? levels[nd.pred].label : cartFmt_(nd.mean); };
  var F3 = { type: 'fixed', digits: 3 };
  var blocks = [];
  blocks.push({
    type: 'table', caption: '처리 요약',
    columns: [{ key: 'label', label: '' }, { key: 'n', label: 'N', format: { type: 'int' } }],
    rows: [{ label: '분석에 포함', n: n }, { label: '결측으로 제외', n: dataset.n - n }, { _kind: 'total', label: '합계', n: dataset.n }],
    notes: ['종속변수 "' + cartName_(dep[0]) + '": ' + (mode === 'class' ? '분류 나무 (지니 불순도, ' + C + '개 범주)' : '회귀 나무 (제곱오차)') + '. 예측변수 ' + P.length + '개 (숫자 ' + covs.length + ', 범주형 ' + facs.length + '). 목록별 결측 제거.']
  });
  if (mode === 'class') {
    blocks.push({
      type: 'table', caption: '종속변수 분포 (' + cartName_(dep[0]) + ')',
      columns: [{ key: 'orig', label: '값' }, { key: 'n', label: 'N', format: { type: 'int' } }, { key: 'pct', label: '%', format: { type: 'fixed', digits: 1 } }],
      rows: levels.map(function (l) { return { orig: l.label, n: l.count, pct: l.count / n * 100 }; })
    });
  }
  var fullLeaves = seq.leaves[0];
  var selRow = [
    { label: '최대 깊이 / 분할 최소 사례 / 잎 최소 사례', value: maxDepth + ' / ' + minSplit + ' / ' + minLeaf },
    { label: '가지치기', value: prune === 'none' ? '하지 않음' : (prune === 'cv' ? '최소 CV 오차' : '1 표준오차 규칙') + ' (' + Math.min(nfold, n) + '겹 교차검증, 시드 ' + seed + ')' },
    { label: '키운 나무의 잎 수', value: fullLeaves },
    { _kind: 'total', label: '선택한 나무의 잎 수', value: leafCount }
  ];
  blocks.push({ type: 'table', caption: '나무 요약', columns: [{ key: 'label', label: '' }, { key: 'value', label: '값' }], rows: selRow });

  // 가지치기 곡선
  if (cv) {
    var curve = [], ptsCv = [], ptsTr = [];
    for (var q = 0; q <= m; q++) {
      var trErr = cartSubtree_(T, T.nodes[0], seq.sets[q]).r / n;
      var cvErr = mode === 'class' ? cv.xerror[q] : cv.xerror[q];
      curve.push({ leaves: seq.leaves[q], alpha: q === 0 ? null : seq.alphas[q - 1], train: trErr, cv: cvErr, se: cv.se[q], sel: q === k ? '◀ 선택' : '', _bold: q === k });
      ptsCv.push([seq.leaves[q], cvErr]); ptsTr.push([seq.leaves[q], trErr]);
    }
    blocks.push({
      type: 'table', caption: '가지치기 곡선 (비용복잡도 + 교차검증)',
      columns: [{ key: 'leaves', label: '잎 수', format: { type: 'int' } }, { key: 'alpha', label: '복잡도 α', format: { type: 'fixed', digits: 4 } },
        { key: 'train', label: mode === 'class' ? '훈련 오분류율' : '훈련 MSE', format: { type: 'fixed', digits: 4 } },
        { key: 'cv', label: mode === 'class' ? 'CV 오분류율' : 'CV MSE', format: { type: 'fixed', digits: 4 } }, { key: 'se', label: 'CV 표준오차', format: { type: 'fixed', digits: 4 } }, { key: 'sel', label: '' }],
      rows: curve,
      notes: ['복잡도 α가 클수록 가지를 많이 칩니다 (첫 줄은 키운 그대로). CV 오차는 ' + Math.min(nfold, n) + '겹 교차검증으로 추정한 새 자료에서의 오차이고, 훈련 오차는 항상 더 낙관적이다. ' +
        (prune === 'cv1se' ? '1 표준오차 규칙: 최소 CV 오차에서 표준오차 하나 이내인 가장 작은 나무를 골랐다.' : '최소 CV 오차인 나무(동점이면 더 작은 나무)를 골랐다.')]
    });
    blocks.push(d.lineChartBlock('잎 수와 오차 (가지치기 곡선)', '잎 수', mode === 'class' ? '오분류율' : 'MSE', [{ label: 'CV 오차', points: ptsCv }, { label: '훈련 오차', points: ptsTr }]));
  }

  // 나무 구조 표
  var treeCols = [{ key: 'id', label: '노드', format: { type: 'int' } }, { key: 'rule', label: '들어오는 조건' }, { key: 'n', label: 'N', format: { type: 'int' } }, { key: 'pred', label: mode === 'class' ? '예측' : '평균' }];
  if (mode === 'class') levels.forEach(function (l, c) { treeCols.push({ key: 'c' + c, label: l.label + ' %', format: { type: 'fixed', digits: 1 } }); });
  treeCols.push({ key: 'imp', label: mode === 'class' ? '지니 불순도' : '분산', format: F3 }, { key: 'split', label: '이 노드의 분할' });
  var splitText = function (nd) {
    if (nd.split === null) return '잎';
    var s = nd.split;
    return s.kind === 'num' ? P[s.j].name + ' ≤ ' + cartFmt_(s.arg) : P[s.j].name + ' ∈ ' + catLabels(nd, 'L');
  };
  blocks.push({
    type: 'table', caption: '나무 구조 (전위 순서)', columns: treeCols,
    rows: nodes.map(function (nd) {
      var r = { id: nd.id, rule: nd.parent === 0 ? '(전체)' : edgeText(nd), n: nd.n, pred: predText(nd), imp: nd.imp, split: splitText(nd), _bold: nd.split === null };
      if (mode === 'class') levels.forEach(function (l, c) { r['c' + c] = nd.counts[c] / nd.n * 100; });
      return r;
    }),
    notes: ['굵은 줄이 잎(마지막 노드)이다. "들어오는 조건"은 부모 노드의 분할 중 이 노드로 오는 쪽이다. 같은 변수에서 숫자 변수는 값이 작거나 같으면(≤) 왼쪽, 범주형은 목록에 들면(∈) 왼쪽으로 간다.']
  });

  // 나무 그림
  if (leafCount <= CART_DIAGRAM_MAX_LEAVES_) {
    blocks.push(d.treeChartBlock('나무 그림', nodes.map(function (nd) {
      var lines = [];
      if (nd.split !== null) lines.push(splitText(nd) + '?');
      lines.push('N=' + nd.n);
      lines.push((mode === 'class' ? '예측: ' : '평균 ') + predText(nd));
      if (mode === 'class') lines.push(levels.map(function (l, c) { return nd.counts[c] ? l.label + ' ' + nd.counts[c] : null; }).filter(Boolean).join(' / '));
      return { id: nd.id, left: nd.left, right: nd.right, edge: edgeText(nd), lines: lines, cls: mode === 'class' ? nd.pred : -1 };
    })));
  } else {
    warnings.push('잎이 ' + leafCount + '개로 많아 나무 그림은 생략했습니다 (표는 그대로). 최대 깊이를 줄이거나 가지치기를 켜 보세요.');
  }

  // 변수 중요도
  var impo = {}, order = [], total = 0;
  nodes.forEach(function (nd) {
    if (nd.split === null) return;
    var nm = P[nd.split.j].name;
    if (!(nm in impo)) { impo[nm] = 0; order.push(nm); }
    impo[nm] += nd.split.gain; total += nd.split.gain;
  });
  if (order.length) {
    order.sort(function (a2, b2) { return impo[b2] - impo[a2]; });
    blocks.push({
      type: 'table', caption: '변수 중요도',
      columns: [{ key: 'name', label: '변수' }, { key: 'gain', label: '불순도 감소 합', format: F3 }, { key: 'pct', label: '비율 %', format: { type: 'fixed', digits: 1 } }],
      rows: order.map(function (nm) { return { name: nm, gain: impo[nm], pct: impo[nm] / total * 100 }; }),
      notes: ['선택한 나무에서 그 변수로 분할할 때 줄어든 불순도(' + (mode === 'class' ? 'n × 지니' : '제곱합') + ')의 합이다. 대리 분할은 쓰지 않아 쓰이지 않은 변수는 0이다.']
    });
  }

  // 성능: 재대입(훈련) + 교차검증
  var predAll = rowsOk.map(function (r, ii) { return cartLeaf_(T, ii, collapsed); });
  if (mode === 'class') {
    var conf = [], hit = 0;
    for (var a1 = 0; a1 < C; a1++) conf.push(new Array(C).fill(0));
    predAll.forEach(function (lf, ii) { conf[y[ii]][lf.pred]++; if (y[ii] === lf.pred) hit++; });
    var crows = levels.map(function (l, a2) {
      var r = { obs: l.label }, tot = 0;
      for (var b2 = 0; b2 < C; b2++) { r['p' + b2] = conf[a2][b2]; tot += conf[a2][b2]; }
      r.pct = tot ? conf[a2][a2] / tot * 100 : null;
      return r;
    });
    var ccols = [{ key: 'obs', label: '관측' }];
    levels.forEach(function (l, b2) { ccols.push({ key: 'p' + b2, label: '예측: ' + l.label, format: { type: 'int' } }); });
    ccols.push({ key: 'pct', label: '정확도(%)', format: { type: 'fixed', digits: 1 } });
    crows.push({ _kind: 'total', obs: '전체 정확도', pct: hit / n * 100 });
    var cnotes = ['선택한 나무로 훈련 자료를 다시 예측한 결과(재대입)라 새 자료에서의 정확도보다 낙관적이다.'];
    if (cv) cnotes.push('교차검증 정확도(새 자료 추정): ' + ((1 - cv.xerror[k]) * 100).toFixed(1) + '% (표준오차 ' + (cv.se[k] * 100).toFixed(1) + '%p).');
    if (C === 2) {
      var tp = conf[1][1], fn = conf[1][0], fp = conf[0][1], tn = conf[0][0];
      cnotes.push('"' + levels[1].label + '"을(를) 양성으로 보면 민감도 ' + (tp + fn ? (tp / (tp + fn) * 100).toFixed(1) : '.') + '%, 특이도 ' + (tn + fp ? (tn / (tn + fp) * 100).toFixed(1) : '.') + '%.');
    }
    blocks.push({ type: 'table', caption: '분류표 (재대입)', columns: ccols, rows: crows, notes: cnotes });
  } else {
    var sse = 0, sst = 0, mean = 0;
    y.forEach(function (v) { mean += v / n; });
    predAll.forEach(function (lf, ii) { sse += (y[ii] - lf.pred) * (y[ii] - lf.pred); sst += (y[ii] - mean) * (y[ii] - mean); });
    var rrows = [{ label: 'R² (재대입)', value: sst > 0 ? 1 - sse / sst : null }, { label: 'RMSE (재대입)', value: Math.sqrt(sse / n) }];
    if (cv) {
      rrows.push({ label: 'RMSE (교차검증)', value: Math.sqrt(cv.xerror[k]) });
      rrows.push({ label: 'R² (교차검증)', value: sst > 0 ? 1 - cv.xerror[k] * n / sst : null });
    }
    blocks.push({
      type: 'table', caption: '모형 성능', columns: [{ key: 'label', label: '' }, { key: 'value', label: '값', format: F3 }], rows: rrows,
      notes: ['재대입 값은 훈련 자료를 다시 예측한 것이라 새 자료에서의 성능보다 낙관적이다. 교차검증 값이 새 자료에서의 추정이다.']
    });
  }
  warnings.push('의사결정나무는 자료가 조금만 바뀌어도 모양이 크게 달라질 수 있고(불안정성) 깊을수록 과적합하기 쉽습니다. 교차검증 오차를 함께 보고 해석하세요.');

  var result = { analysis: 'cart', title: '의사결정나무', blocks: blocks, warnings: warnings };
  if (opt.save) {
    var pv = new Array(dataset.n).fill(null), nd2 = new Array(dataset.n).fill(null), vl = null;
    rowsOk.forEach(function (r, ii) {
      var lf = predAll[ii], exp = null;
      for (var z = 0; z < nodes.length; z++) if (nodes[z].old === lf.id) { exp = nodes[z]; break; }
      pv[r] = mode === 'class' ? levels[lf.pred].value : lf.mean;
      nd2[r] = exp.id;
    });
    if (mode === 'class') vl = levels.map(function (l) { return String(l.value).replace(/[,;=\n]+/g, ' ') + '=' + String(l.label).replace(/[,;=\n]+/g, ' '); }).join(', ');
    result.saveColumns = [
      { name: '나무예측', label: '의사결정나무 예측 (' + cartName_(dep[0]) + ')', scale: mode === 'class' ? '명목' : '연속', values: pv, fixed: mode === 'class' ? undefined : 4, valueLabels: vl || undefined },
      { name: '나무노드', label: '의사결정나무 소속 잎 노드 번호', scale: '명목', values: nd2 }
    ];
  }
  return result;
}

if (typeof module !== 'undefined') {
  module.exports = { cartBuild: cartBuild, cartPruneSequence: cartPruneSequence, cartCrossValidate: cartCrossValidate, cartChoose: cartChoose, cartExport: cartExport, cartLeaf_: cartLeaf_, runCart: runCart };
}
