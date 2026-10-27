/**
 * 빈도분석, 기술통계. 순수 함수: (dataset, request) → 결과 객체 (doc/설계-결과객체.md)
 * 계산식: 표본 표준편차(n-1), 왜도 G1, 첨도 G2(초과첨도)와 그 표준오차.
 */

/**
 * Node 테스트에서는 require, GAS에서는 같은 전역 범위의 함수를 쓴다 (호출 시점에 찾는다).
 * 이름은 파일마다 고유해야 한다: GAS는 모든 파일이 한 전역 범위라 같은 이름은 서로 덮어쓴다.
 */
function descDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var ch = require('./Charts.js');
    return {
      getVars: data.getVars, selectComplete: data.selectComplete, compareValues_: data.compareValues_,
      barChartBlock: ch.barChartBlock, histogramBlock: ch.histogramBlock, CHART_MAX_BARS_: ch.CHART_MAX_BARS_
    };
  }
  return {
    getVars: getVars, selectComplete: selectComplete, compareValues_: compareValues_,
    barChartBlock: barChartBlock, histogramBlock: histogramBlock, CHART_MAX_BARS_: CHART_MAX_BARS_
  };
}

var FREQ_MAX_VALUES_ = 500;   // 빈도표 한 개의 최대 줄 수 (고유한 값의 개수)

function displayName_(v) {
  return v.label ? v.label : v.name;
}

/**
 * 평균과 중심 적률의 합 c2 = Σ(x-m)², c3 = Σ(x-m)³, c4 = Σ(x-m)⁴ (n ≥ 1).
 * 평균을 한 번 보정한 뒤 편차의 거듭제곱합으로 구하고, 평균이 double로 반올림된 오차(e = Σd/n)는
 * 적률 항등식으로 보정한다. 평균이 1e6 같은 큰 값이어도 정확 산술 대비 1e-12 수준을 유지한다.
 */
function centralSums(values) {
  var n = values.length, i, sum = 0;
  for (i = 0; i < n; i++) sum += values[i];
  var mean = sum / n;
  var corr = 0;
  for (i = 0; i < n; i++) corr += values[i] - mean;
  mean += corr / n;

  var p1 = 0, p2 = 0, p3 = 0, p4 = 0;
  for (i = 0; i < n; i++) {
    var d = values[i] - mean, d2 = d * d;
    p1 += d; p2 += d2; p3 += d2 * d; p4 += d2 * d2;
  }
  var e = p1 / n;
  return {
    n: n,
    sum: sum,
    mean: mean,
    c2: p2 - n * e * e,
    c3: p3 - 3 * e * p2 + 2 * n * e * e * e,
    c4: p4 - 4 * e * p3 + 6 * e * e * p2 - 3 * n * e * e * e * e
  };
}

/** 숫자 배열(결측 제외)의 기술통계. 계산할 수 없는 값은 null. */
function describeColumn(values) {
  var n = values.length;
  var r = {
    n: n, min: null, max: null, range: null, sum: null, mean: null,
    variance: null, sd: null, seMean: null,
    skewness: null, seSkewness: null, kurtosis: null, seKurtosis: null
  };
  if (n === 0) return r;

  var min = values[0], max = values[0], i;
  for (i = 1; i < n; i++) {
    if (values[i] < min) min = values[i];
    if (values[i] > max) max = values[i];
  }
  var cs = centralSums(values);
  var c2 = cs.c2, c3 = cs.c3, c4 = cs.c4;
  r.min = min; r.max = max; r.range = max - min; r.sum = cs.sum; r.mean = cs.mean;

  if (n > 1) {
    r.variance = c2 / (n - 1);
    r.sd = Math.sqrt(r.variance);
    r.seMean = r.sd / Math.sqrt(n);
  }
  if (n > 2) {
    r.seSkewness = Math.sqrt(6 * n * (n - 1) / ((n - 2) * (n + 1) * (n + 3)));
  }
  if (n > 3) {
    r.seKurtosis = 2 * r.seSkewness * Math.sqrt((n * n - 1) / ((n - 3) * (n + 5)));
  }
  if (r.sd !== null && r.sd > 0) {
    if (n > 2) r.skewness = n / ((n - 1) * (n - 2)) * c3 / (r.sd * r.sd * r.sd);
    if (n > 3) {
      r.kurtosis = n * (n + 1) / ((n - 1) * (n - 2) * (n - 3)) * c4 / (r.variance * r.variance) -
        3 * (n - 1) * (n - 1) / ((n - 2) * (n - 3));
    }
  }
  return r;
}

/* ---------- 빈도분석 ---------- */

function frequencyTable_(v, column, N) {
  var tally = {}, keys = [], missing = 0;
  column.forEach(function (x) {
    if (x === null) { missing++; return; }
    var k = typeof x + ':' + x;
    if (!(k in tally)) { tally[k] = { value: x, count: 0 }; keys.push(k); }
    tally[k].count++;
  });
  var entries = keys.map(function (k) { return tally[k]; }).sort(descDeps_().compareValues_);
  var validN = N - missing;

  var rows = [], cum = 0;
  entries.forEach(function (e, i) {
    cum += e.count;
    var label = v.valueLabels[String(e.value)];
    rows.push({
      group: i === 0 ? '유효' : '',
      value: label !== undefined ? label : String(e.value),
      freq: e.count,
      pct: e.count / N * 100,
      validPct: e.count / validN * 100,
      cumPct: cum / validN * 100
    });
  });
  if (missing > 0) {
    if (validN > 0) rows.push({ _kind: 'subtotal', group: '', value: '합계', freq: validN, pct: validN / N * 100, validPct: 100 });
    rows.push({ group: '결측', value: '결측', freq: missing, pct: missing / N * 100 });
  }
  rows.push({ _kind: 'total', group: '합계', value: '', freq: N, pct: 100 });

  return {
    unique: entries.length,
    chart: {
      labels: entries.map(function (e) {
        var label = v.valueLabels[String(e.value)];
        return label !== undefined ? label : String(e.value);
      }),
      counts: entries.map(function (e) { return e.count; })
    },
    block: {
      type: 'table',
      caption: displayName_(v),
      columns: [
        { key: 'group', label: '' },
        { key: 'value', label: '값' },
        { key: 'freq', label: '빈도', format: { type: 'int' } },
        { key: 'pct', label: '퍼센트', format: { type: 'fixed', digits: 1 } },
        { key: 'validPct', label: '유효 퍼센트', format: { type: 'fixed', digits: 1 } },
        { key: 'cumPct', label: '누적 퍼센트', format: { type: 'fixed', digits: 1 } }
      ],
      rows: rows,
      notes: []
    }
  };
}

function runFrequencies(dataset, request) {
  var dd = descDeps_();
  var opt = request.options || {};
  var vars = dd.getVars(dataset, request.vars);
  var N = dataset.n;
  if (N === 0) throw new Error('데이터가 없습니다.');
  var warnings = [];

  var summaryRows = vars.map(function (v) {
    var missing = dataset.data[v.name].filter(function (x) { return x === null; }).length;
    return { name: displayName_(v), valid: N - missing, missing: missing };
  });
  var blocks = [{
    type: 'table',
    caption: '통계량',
    columns: [
      { key: 'name', label: '변수' },
      { key: 'valid', label: '유효', format: { type: 'int' } },
      { key: 'missing', label: '결측', format: { type: 'int' } }
    ],
    rows: summaryRows,
    notes: []
  }];

  vars.forEach(function (v) {
    var t = frequencyTable_(v, dataset.data[v.name], N);
    if (t.unique > FREQ_MAX_VALUES_) {
      // 연속형 변수로 빈도표를 만들면 사례 수만큼 줄이 생겨 화면과 전송이 무거워진다
      warnings.push('"' + displayName_(v) + '"의 고유한 값이 ' + t.unique + '개라서 빈도표를 만들지 않았습니다 (최대 ' + FREQ_MAX_VALUES_ + '개). 연속형 변수는 기술통계(히스토그램)를 쓰세요.');
      return;
    }
    if (t.unique > 100) {
      warnings.push('"' + displayName_(v) + '"의 고유한 값이 ' + t.unique + '개입니다. 연속형 변수라면 기술통계가 더 적합할 수 있습니다.');
    }
    blocks.push(t.block);
    if (opt.barChart) {
      if (t.unique > dd.CHART_MAX_BARS_) {
        warnings.push('"' + displayName_(v) + '"의 값이 ' + t.unique + '개라서 막대그래프를 그리지 않았습니다 (최대 ' + dd.CHART_MAX_BARS_ + '개).');
      } else if (t.unique > 0) {
        blocks.push(dd.barChartBlock('막대그래프: ' + displayName_(v), t.chart.labels, t.chart.counts, displayName_(v)));
      }
    }
  });
  return { analysis: 'frequencies', title: '빈도분석', blocks: blocks, warnings: warnings };
}

/* ---------- 기술통계 ---------- */

function runDescriptives(dataset, request) {
  var d = descDeps_();
  var vars = d.getVars(dataset, request.vars, { numeric: true });
  var names = vars.map(function (v) { return v.name; });
  var opt = request.options || {};
  var warnings = [];

  var complete = d.selectComplete(dataset, names);
  var columns = {};
  if (opt.listwise) {
    columns = complete.columns;
    if (complete.dropped > 0) {
      warnings.push('목록별 결측 제거로 ' + complete.dropped + '건을 제외했습니다.');
    }
  } else {
    names.forEach(function (nm) {
      columns[nm] = dataset.data[nm].filter(function (x) { return x !== null; });
    });
  }

  var cols = [
    { key: 'name', label: '변수' },
    { key: 'n', label: 'N', format: { type: 'int' } }
  ];
  if (opt.range) cols.push({ key: 'range', label: '범위', format: { type: 'trim', digits: 3 } });
  cols.push(
    { key: 'min', label: '최솟값', format: { type: 'trim', digits: 3 } },
    { key: 'max', label: '최댓값', format: { type: 'trim', digits: 3 } }
  );
  if (opt.sum) cols.push({ key: 'sum', label: '합계', format: { type: 'trim', digits: 3 } });
  cols.push({ key: 'mean', label: '평균', format: { type: 'fixed', digits: 3 } });
  if (opt.seMean) cols.push({ key: 'seMean', label: '평균의 표준오차', format: { type: 'fixed', digits: 3 } });
  cols.push({ key: 'sd', label: '표준편차', format: { type: 'fixed', digits: 3 } });
  if (opt.variance) cols.push({ key: 'variance', label: '분산', format: { type: 'fixed', digits: 3 } });
  if (opt.skewness) {
    cols.push({ key: 'skewness', label: '왜도', format: { type: 'fixed', digits: 3 } });
    cols.push({ key: 'seSkewness', label: '왜도 표준오차', format: { type: 'fixed', digits: 3 } });
  }
  if (opt.kurtosis) {
    cols.push({ key: 'kurtosis', label: '첨도', format: { type: 'fixed', digits: 3 } });
    cols.push({ key: 'seKurtosis', label: '첨도 표준오차', format: { type: 'fixed', digits: 3 } });
  }

  var rows = vars.map(function (v) {
    var stats = describeColumn(columns[v.name]);
    if (stats.n === 0) warnings.push('"' + displayName_(v) + '"에 유효한 값이 없습니다.');
    stats.name = displayName_(v);
    return stats;
  });
  rows.push({ _kind: 'total', name: '유효 N(목록별)', n: complete.n });

  var result = {
    analysis: 'descriptives',
    title: '기술통계',
    blocks: [{
      type: 'table',
      caption: '기술통계량',
      columns: cols,
      rows: rows,
      notes: [opt.listwise
        ? '목록별 결측 제거: 선택한 모든 변수가 유효한 사례만 사용했습니다.'
        : '변수별로 유효한 사례를 사용했습니다. 마지막 행은 모든 변수가 유효한 사례 수입니다.']
    }],
    warnings: warnings
  };
  if (opt.histogram) {
    vars.forEach(function (v) {
      if (columns[v.name].length > 0) {
        result.blocks.push(d.histogramBlock('히스토그램: ' + displayName_(v), columns[v.name], displayName_(v)));
      }
    });
  }
  return result;
}

if (typeof module !== 'undefined') {
  module.exports = {
    centralSums: centralSums,
    describeColumn: describeColumn,
    runFrequencies: runFrequencies,
    runDescriptives: runDescriptives
  };
}
