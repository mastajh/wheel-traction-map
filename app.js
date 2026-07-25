/* Bekker–Wong 트랙션 맵 — 농업 보조 4륜차
 * 압력-침하: p = (kc/b + kφ)·zⁿ  (Bekker)
 * 전단응력:  τ = (c + σ·tanφ)·(1 − e^(−j/K)),  j = s·x  (Janosi–Hanamoto)
 * 주행 모드: 순견인력 T = H(총견인력) − Rc(압실저항), 바퀴당
 * 버티기 모드(능동 서스펜션으로 차체 수평, 정적하중):
 *   수직하중 W를 지반법선 N = W·cosθ, 경사방향 요구전단 F = W·sinθ로 분해
 *   안전율 SF = (c·A + N·tanφ) / F  —  SF < 1 이면 실속(미끄러짐)
 */
"use strict";

/* ---------- 감마 함수 (Lanczos 근사) : ∫₀¹(1−u²)ⁿdu 계산용 ---------- */
function gamma(z) {
  const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  if (z < 0.5) return Math.PI / (Math.sin(Math.PI * z) * gamma(1 - z));
  z -= 1;
  let x = c[0];
  for (let i = 1; i < 9; i++) x += c[i] / (z + i);
  const t = z + 7.5;
  return Math.sqrt(2 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;
}
const In = n => Math.sqrt(Math.PI) * gamma(n + 1) / (2 * gamma(n + 1.5));

/* ---------- 토양 파라미터 (SI: N, m, Pa, rad) ----------
 * 문헌: J.Y. Wong, "Terramechanics and Off-Road Vehicle Engineering" 2nd ed. (2009)
 *       압력-침하·전단 파라미터 표 (광물 지반, LLL/Michigan/WES 실측 시리즈)
 *       검색 인용: Li et al. (2015, doi:10.1155/2015/293125), Guo et al. (2025, MDPI Appl. Sci. 15:6566) */
const SOILS = {
  loam:  { name: "양토",     kc: 52.53e3, kphi: 1127.97e3, n: 0.9, c: 4.83e3, phi: 20 * Math.PI / 180, K: 0.025 },  // Wong — Sandy loam (Michigan)
  clay:  { name: "점토",     kc: 13.19e3, kphi: 692.15e3,  n: 0.5, c: 4.14e3, phi: 13 * Math.PI / 180, K: 0.010 },  // Wong — Clayey soil
  sand:  { name: "마른 모래", kc: 0.99e3,  kphi: 1528.43e3, n: 1.1, c: 1.04e3, phi: 28 * Math.PI / 180, K: 0.025 },  // Wong — Dry sand (LLL)
  paddy: { name: "습한 논흙", kc: 3.0e3,   kphi: 45.0e3,    n: 0.7, c: 1.50e3, phi: 20 * Math.PI / 180, K: 0.020 },  // 추정치 (문헌값 아님)
};

/* ---------- 트레드 프리셋 ----------
 * m  : 마찰각 보정계수 — 유효 마찰각 δ = m·φ
 *      흙-금속/고무 마찰비 δ/φ ≈ 0.5~0.9 (Duncan & Mokwa 2001; Fine 2011 권장표)
 *      매끈한 바퀴 0.70 (강·고무 smooth), 일반 트레드 0.85 (rough),
 *      깊은 러그 1.00 (러그 끝에서 흙-흙 전단, δ→φ)
 * lug: 러그 수동저항 부가견인 계수 — 그루서가 흙을 파고 옹벽처럼 미는 추력
 *      H → H·(1+lug)  (Wong 2009 그루서 추력 모델의 1차 근사; 계수는 가정치,
 *      연약지반에서 유효, 단단한 지면에서는 과대평가될 수 있음) */
const TREADS = {
  slick:  { name: "매끈 강성륜", m: 0.70, lug: 0.00 },
  rib:    { name: "일반 트레드", m: 0.85, lug: 0.10 },
  aglug:  { name: "농업용 러그", m: 1.00, lug: 0.25 },
  paddle: { name: "패들 러그",   m: 1.00, lug: 0.40 },
};

/* ---------- 접지패치 기하 (경량) : 침하·접지길이·면적 ---------- */
function contactPatch(W, D, b, soil) {
  const { kc, kphi, n } = soil;
  const k = kc / b + kphi;
  const z0 = Math.pow(W / (b * k * Math.sqrt(D) * In(n)), 1 / (n + 0.5));
  if (!isFinite(z0) || z0 >= D * 0.5) return null;   // 매몰 → 무효
  const L = Math.sqrt(D * z0 - z0 * z0);
  return { z0, L, A: b * L, k };
}

/* ---------- 바퀴 1개 해석 ----------
 * W: 바퀴하중(N, 지반법선 성분), D: 직경(m), b: 폭(m), s: 슬립률
 * 접지패치: 진입점(x=0) → 최저점(x=L), 침하 프로파일 z(x) = z0 − (L−x)²/D */
function wheelTraction(W, D, b, soil, s, tread = TREADS.aglug) {
  const cp = contactPatch(W, D, b, soil);
  if (!cp) return null;
  const { z0, L, A, k } = cp;
  const { n, c, phi, K } = soil;
  const tanp = Math.tan(phi * tread.m);
  const N = 120, dx = L / N;
  let H = 0;
  for (let i = 0; i <= N; i++) {
    const x = i * dx;
    let z = z0 - (L - x) * (L - x) / D;
    if (z < 0) z = 0;
    const tau = (c + k * Math.pow(z, n) * tanp) * (1 - Math.exp(-s * x / K));
    H += tau * (i === 0 || i === N ? 0.5 : 1);
  }
  H *= b * dx * (1 + tread.lug);
  const Rc = b * k * Math.pow(z0, n + 1) / (n + 1);
  return { z0, L, A, H, Rc, T: H - Rc };
}

/* ---------- 버티기: 최대 경사각 θ_max (이분탐색, SF=1 해) ---------- */
function thetaMax(Wv, D, b, soil, tread = TREADS.aglug) {
  const lo0 = 0.5 * Math.PI / 180, hi0 = 85 * Math.PI / 180;
  if (!contactPatch(Wv * Math.cos(lo0), D, b, soil)) return null;  // 평지에서도 매몰
  const f = th => {
    const N = Wv * Math.cos(th);
    const cp = contactPatch(N, D, b, soil);
    if (!cp) return -1e9;
    return (soil.c * cp.A + N * Math.tan(soil.phi * tread.m)) * (1 + tread.lug) - Wv * Math.sin(th);
  };
  if (f(hi0) > 0) return 85;
  if (f(lo0) <= 0) return 0;
  let lo = lo0, hi = hi0;
  for (let i = 0; i < 16; i++) {
    const mid = (lo + hi) / 2;
    if (f(mid) > 0) lo = mid; else hi = mid;
  }
  return (lo + hi) / 2 * 180 / Math.PI;
}

/* ---------- 격자 스윕 ---------- */
function computeGrid(p) {
  const soil = SOILS[p.soil];
  const tread = TREADS[p.tread];
  const W = p.mass * 9.81 / 4;                    // 바퀴당 수직하중 (차체 수평, 균등분배)
  const th = p.slope * Math.PI / 180;
  const Nn = W * Math.cos(th);                    // 지반법선 하중 (버티기 모드)
  const Fv = W * Math.sin(th);                    // 경사방향 요구전단 (버티기 모드)
  const ND = 46, NB = 41;
  const Ds = [], bs = [];
  for (let j = 0; j < ND; j++) Ds.push(p.dmin + (p.dmax - p.dmin) * j / (ND - 1));
  for (let i = 0; i < NB; i++) bs.push(p.bmin + (p.bmax - p.bmin) * i / (NB - 1));
  const T = [], Z = [], TM = [];
  for (let i = 0; i < NB; i++) {
    T.push([]); Z.push([]); TM.push([]);
    for (let j = 0; j < ND; j++) {
      if (p.mode === "drive") {
        const r = wheelTraction(W, Ds[j], bs[i], soil, p.slip, tread);
        T[i].push(r ? r.T : null);
        Z[i].push(r ? r.z0 / Ds[j] : null);
        TM[i].push(null);
      } else {
        const cp = contactPatch(Nn, Ds[j], bs[i], soil);
        if (!cp || Fv < 1e-9) { T[i].push(null); Z[i].push(null); TM[i].push(null); continue; }
        const Hhold = (soil.c * cp.A + Nn * Math.tan(soil.phi * tread.m)) * (1 + tread.lug);
        T[i].push(Hhold / Fv);                                  // 안전율 SF
        Z[i].push(cp.z0 / Ds[j]);
        TM[i].push(thetaMax(W, Ds[j], bs[i], soil, tread));
      }
    }
  }
  return { Ds, bs, T, Z, TM, W, Nn, Fv };
}

/* ---------- 등값선 추출: 각 D 열에서 b축 선형보간 ---------- */
function extractLine(Ds, bs, grid, level) {
  const outD = [], outB = [];
  for (let j = 0; j < Ds.length; j++) {
    let prev = null, prevB = null;
    for (let i = 0; i < bs.length; i++) {
      const v = grid[i][j];
      if (v === null) continue;
      if (prev !== null && (prev - level) * (v - level) < 0) {
        outD.push(Ds[j]);
        outB.push(prevB + (bs[i] - prevB) * (level - prev) / (v - prev));
        break;
      }
      prev = v; prevB = bs[i];
    }
  }
  return { D: outD, b: outB };
}

/* ---------- 파워로 피팅 b = α·D^β (log-log 최소제곱) ---------- */
function fitPower(Dl, bl) {
  const n = Dl.length;
  if (n < 3) return null;
  let sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (let i = 0; i < n; i++) {
    const x = Math.log(Dl[i]), y = Math.log(bl[i]);
    sx += x; sy += y; sxx += x * x; sxy += x * y;
  }
  const beta = (n * sxy - sx * sy) / (n * sxx - sx * sx);
  const alpha = Math.exp((sy - beta * sx) / n);
  let maxErr = 0;
  for (let i = 0; i < n; i++)
    maxErr = Math.max(maxErr, Math.abs(alpha * Math.pow(Dl[i], beta) - bl[i]) / bl[i]);
  return { alpha, beta, maxErr };
}

/* ---------- 상태 ---------- */
const state = { mode: "drive", soil: "loam", tread: "aglug", mass: 1000, slip: 0.20, tgt: 0.15,
                slope: 12, sfTgt: 1.5,
                dmin: 0.30, dmax: 1.20, bmin: 0.10, bmax: 0.50 };

const $ = id => document.getElementById(id);
const fmt = (v, d = 2) => v.toFixed(d);

const DIVERGE = [[0, "#c2261f"], [0.25, "#e3695f"], [0.5, "#f5f5f7"], [0.75, "#8fd6a0"], [1, "#2f9e44"]];
const SFSCALE = [[0, "#c2261f"], [0.4, "#e3695f"], [0.5, "#f5f5f7"], [0.7, "#8fd6a0"], [1, "#2f9e44"]];
const PLOTLT = { paper_bgcolor: "#f5f5f7", plot_bgcolor: "#f5f5f7",
                 font: { family: "-apple-system, BlinkMacSystemFont, 'Apple SD Gothic Neo', 'Noto Sans KR', sans-serif",
                         color: "#86868b", size: 13 } };
const LC = { red: "#ff3b30", blue: "#0066cc", teal: "#30b0c7" };   // Apple 시스템 팔레트
const CFG = { displayModeBar: false, responsive: true };

/* 모드별 UI 텍스트 */
const MODE_TXT = {
  drive: {
    sub: "농업 보조 4륜차 · 휠 직경 × 폭 → 바퀴당 순견인력 표면",
    lg0: "자력주행 경계 — 순견인력 0",
    lg1: "목표 견인 달성선",
    eq0t: "자력주행 경계 — 이 선 아래면 견인력 &lt; 저항 (못 움직임)",
    d0: "이 선 위쪽이어야 바퀴가 스스로 굴릴 수 있습니다.",
    d1: "이 선 위쪽이면 목표 견인력을 확보합니다.",
    zAxis: "순견인력 (N)", cbar: "N/바퀴",
    t2d: "설계 평면 — D × b 가이드 맵 (숫자: 순견인력 N/바퀴)",
  },
  hold: {
    sub: "능동 서스펜션 차체 수평 · 정적 버티기 — 안전율 SF 표면",
    lg0: "실속 경계 — SF = 1",
    lg1: "목표 안전율선",
    eq0t: "실속 경계 — 이 선 아래는 미끄러짐 (SF &lt; 1)",
    d0: "이 선 위쪽이어야 경사면에서 바퀴가 버팁니다.",
    d1: "설계 권장 마진을 확보하는 선입니다.",
    zAxis: "안전율 SF (−)", cbar: "SF",
    t2d: "설계 평면 — D × b 가이드 맵 (숫자: 안전율 SF)",
  },
};

/* ---------- 메인 계산 + 렌더 ---------- */
function update() {
  if (state.dmin >= state.dmax) state.dmax = state.dmin + 0.05;
  if (state.bmin >= state.bmax) state.bmax = state.bmin + 0.05;
  const p = { ...state };
  const g = computeGrid(p);
  const hold = p.mode === "hold";
  const lv0 = hold ? 1.0 : 0;                 // 경계선 레벨
  const lv1 = hold ? p.sfTgt : p.tgt * g.W;   // 목표선 레벨

  const l0 = extractLine(g.Ds, g.bs, g.T, lv0);
  const l1 = extractLine(g.Ds, g.bs, g.T, lv1);
  const l2 = extractLine(g.Ds, g.bs, g.Z, 0.25);
  const f0 = fitPower(l0.D, l0.b);
  const f1 = fitPower(l1.D, l1.b);
  const f2 = fitPower(l2.D, l2.b);

  const valid = g.T.flat().filter(v => v !== null);
  const stats = valid.length
    ? { n: valid.length, min: Math.min(...valid), max: Math.max(...valid) }
    : { n: 0, min: NaN, max: NaN };

  draw3d(g, l0, l1, l2, hold);
  draw2d(g, l0, l1, l2, f0, f1, f2, hold);
  const txt = MODE_TXT[p.mode];
  showEq("eq0", "eq0d", f0, txt.d0, stats, lv0);
  showEq("eq1", "eq1d", f1, txt.d1, stats, lv1);
  showEq("eq2", "eq2d", f2, "이 선 아래는 침하가 깊어(직경의 25%↑) 주행성이 급락합니다.", null, null);
  showRec(g, l1, hold, stats, lv1);
}

/* 등값선 위의 표면 높이 값 */
function lineZ(g, line) {
  return line.D.map((D, k) => {
    const j = g.Ds.indexOf(D);
    return g.T[nearestIdx(g.bs, line.b[k])][j];
  });
}
function nearestIdx(arr, v) {
  let bi = 0, bd = 1e9;
  arr.forEach((a, i) => { const d = Math.abs(a - v); if (d < bd) { bd = d; bi = i; } });
  return bi;
}

function draw3d(g, l0, l1, l2, hold) {
  const txt = MODE_TXT[state.mode];
  let colorscale, cmin, cmax, hover;
  if (hold) {
    colorscale = SFSCALE; cmin = 0; cmax = 2;
    hover = "D %{x:.2f} m · b %{y:.2f} m<br>SF %{z:.2f} · θ<sub>max</sub> %{customdata:.1f}°<extra></extra>";
  } else {
    const absMax = Math.max(1, ...g.T.flat().filter(v => v !== null).map(Math.abs));
    colorscale = DIVERGE; cmin = -absMax; cmax = absMax;
    hover = "D %{x:.2f} m · b %{y:.2f} m<br>순견인력 %{z:.0f} N<extra></extra>";
  }
  const surf = {
    type: "surface", x: g.Ds, y: g.bs, z: g.T, customdata: hold ? g.TM : undefined,
    colorscale, cmin, cmax,
    colorbar: { title: txt.cbar, tickfont: { color: "#86868b" }, titlefont: { color: "#86868b" },
                thickness: 12, len: 0.6, outlinewidth: 0 },
    contours: { z: { show: true, color: "rgba(0,0,0,0.16)", usecolormap: false } },
    lighting: { ambient: 0.65, diffuse: 0.8, specular: 0.15, roughness: 0.9 },
    hovertemplate: hover,
  };
  const mk = (line, color, name) => ({
    type: "scatter3d", mode: "lines", x: line.D, y: line.b, z: lineZ(g, line),
    line: { color, width: 7 }, name, hoverinfo: "name", showlegend: true,
  });
  const layout = { ...PLOTLT, margin: { l: 0, r: 0, t: 10, b: 0 },
    legend: { x: 0.02, y: 0.95, bgcolor: "#ffffffd9", font: { size: 12, color: "#1d1d1f" },
              bordercolor: "#e0e0e0", borderwidth: 1 },
    scene: {
      xaxis: { title: "직경 D (m)", color: "#6e6e73", gridcolor: "#e3e3e6", backgroundcolor: "#f5f5f7" },
      yaxis: { title: "폭 b (m)", color: "#6e6e73", gridcolor: "#e3e3e6", backgroundcolor: "#f5f5f7" },
      zaxis: { title: txt.zAxis, color: "#6e6e73", gridcolor: "#e3e3e6", backgroundcolor: "#f5f5f7" },
      camera: { eye: { x: -1.55, y: 1.55, z: 0.9 } },
    } };
  Plotly.react("plot3d", [surf,
    mk(l0, LC.red, hold ? "실속 경계 SF=1" : "자력주행 경계"),
    mk(l1, LC.blue, hold ? "목표 안전율선" : "목표 견인선"),
    mk(l2, LC.teal, "침하 한계선")], layout, CFG);
}

function draw2d(g, l0, l1, l2, f0, f1, f2, hold) {
  const txt = MODE_TXT[state.mode];
  let colorscale, zmin, zmax, hover, zfmt;
  if (hold) {
    colorscale = SFSCALE; zmin = 0; zmax = 2; zfmt = ".2f";
    hover = "D %{x:.2f} · b %{y:.2f}<br>SF %{z:.2f} · θ<sub>max</sub> %{customdata:.1f}°<extra></extra>";
  } else {
    const absMax = Math.max(1, ...g.T.flat().filter(v => v !== null).map(Math.abs));
    colorscale = DIVERGE; zmin = -absMax; zmax = absMax; zfmt = ".0f";
    hover = "D %{x:.2f} · b %{y:.2f}<br>%{z:.0f} N<extra></extra>";
  }
  const heat = {
    type: "contour", x: g.Ds, y: g.bs, z: g.T, customdata: hold ? g.TM : undefined,
    colorscale, zmin, zmax,
    contours: { coloring: "heatmap", showlabels: true,
                labelfont: { color: "#6e6e73", size: 10.5 }, labelformat: zfmt },
    line: { color: "rgba(0,0,0,0.12)", width: 0.5 },
    colorbar: { title: hold ? "SF" : "N", tickfont: { color: "#86868b" }, titlefont: { color: "#86868b" },
                thickness: 10, len: 0.8, outlinewidth: 0 },
    hovertemplate: hover,
  };
  const solid = (line, color) => ({
    type: "scatter", mode: "lines", x: line.D, y: line.b,
    line: { color, width: 2.5 }, hoverinfo: "skip", showlegend: false,
  });
  const fitLine = (f, color) => {
    if (!f) return { type: "scatter", x: [], y: [], showlegend: false };
    const xs = [], ys = [];
    for (let j = 0; j < g.Ds.length; j++) {
      const y = f.alpha * Math.pow(g.Ds[j], f.beta);
      if (y >= g.bs[0] && y <= g.bs[g.bs.length - 1]) { xs.push(g.Ds[j]); ys.push(y); }
    }
    return { type: "scatter", mode: "lines", x: xs, y: ys,
             line: { color, width: 1.5, dash: "dash" }, hoverinfo: "skip", showlegend: false };
  };
  const layout = { ...PLOTLT, margin: { l: 55, r: 10, t: 24, b: 45 },
    title: { text: txt.t2d, font: { size: 13, color: "#86868b" } },
    showlegend: false,
    xaxis: { title: "휠 직경 D (m)", color: "#6e6e73", gridcolor: "#e8e8ea", zeroline: false },
    yaxis: { title: "휠 폭 b (m)", color: "#6e6e73", gridcolor: "#e8e8ea", zeroline: false } };
  Plotly.react("plot2d", [heat,
    solid(l0, LC.red), fitLine(f0, "#ff3b3077"),
    solid(l1, LC.blue), fitLine(f1, "#0066cc77"),
    solid(l2, LC.teal), fitLine(f2, "#30b0c777")], layout, CFG);
}

function showEq(fId, dId, f, desc, stats, level) {
  if (!f) {
    let note = "";
    if (stats && stats.n) {
      if (stats.min >= level) note = " — 전 구간 조건 만족";
      else if (stats.max <= level) note = " — 전 구간 미달";
    }
    $(fId).innerHTML = `<span class="na">현재 범위 내 선 없음${note}</span>`;
    $(dId).textContent = "";
    return;
  }
  $(fId).textContent = `b = ${f.alpha.toFixed(3)} · D^(${f.beta.toFixed(2)})`;
  const ratio = Math.pow(2, f.beta);
  $(dId).innerHTML = `${desc}<br>직경 2배 → 필요 폭 ≈ <b>${ratio >= 1 ? fmt(ratio, 1) + "배" : "1/" + fmt(1 / ratio, 1) + "배"}</b>` +
    ` <span class="err">· 피팅오차 ${(f.maxErr * 100).toFixed(1)}%</span>`;
}

function showRec(g, l1, hold, stats, lv1) {
  const tbl = $("recTable");
  let D, b, headRow = null;
  if (!l1.D.length) {
    /* 선이 범위 밖: 전 구간 만족이면 최소 유효 사양을 제시, 아니면 불가 안내 */
    if (stats.n && stats.min >= lv1) {
      let found = false;
      /* 1차: 목표 + 침하한계(z0/D≤0.25) 동시 만족, 없으면 2차: 목표만 */
      for (const useSink of [true, false]) {
        for (let j = 0; j < g.Ds.length && !found; j++)
          for (let i = 0; i < g.bs.length && !found; i++)
            if (g.T[i][j] !== null && g.T[i][j] >= lv1 &&
                (!useSink || (g.Z[i][j] !== null && g.Z[i][j] <= 0.25))) {
              D = g.Ds[j]; b = g.bs[i]; found = true;
            }
        if (found) break;
      }
      headRow = ["달성 상태", "범위 전체 만족 (최소 유효 사양)"];
    } else {
      const msg = hold
        ? `현재 범위에서 목표 안전율 ${state.sfTgt.toFixed(1)} 달성 불가 — 범위를 넓혀보세요`
        : `현재 범위에서 목표 견인 달성 불가 — 범위를 넓혀보세요`;
      tbl.innerHTML = `<tr><td colspan="2" class="na">${msg}</td></tr>`;
      return;
    }
  } else {
    D = l1.D[0]; b = l1.b[0];
  }
  let rows;
  const tread = TREADS[state.tread];
  if (hold) {
    const th = state.slope * Math.PI / 180;
    const Nn = g.W * Math.cos(th), Fv = g.W * Math.sin(th);
    const r = contactPatch(Nn, D, b, SOILS[state.soil]);
    const Hhold = (SOILS[state.soil].c * r.A + Nn * Math.tan(SOILS[state.soil].phi * tread.m)) * (1 + tread.lug);
    const tm = thetaMax(g.W, D, b, SOILS[state.soil], tread);
    rows = [
      ["휠 직경 D", `≥ ${fmt(D)} m`],
      ["휠 폭 b", `≥ ${fmt(b)} m`],
      ["직경/폭 비", fmt(D / b, 1) + " : 1"],
      [`침하량 (θ=${state.slope}°)`, fmt(r.z0 * 100, 1) + " cm"],
      ["접지 면적 (바퀴당)", fmt(r.A * 1e4, 0) + " cm²"],
      ["안전율 SF", fmt(Hhold / Fv, 2)],
      ["최대 버티기 각도 θ_max", tm === null ? "매몰 한계 (사양 상향 필요)" : fmt(tm, 1) + "°"],
      ["바퀴당 버티기 여유력", fmt(Hhold - Fv, 0) + " N"],
    ];
  } else {
    const r = wheelTraction(g.W, D, b, SOILS[state.soil], state.slip, tread);
    rows = [
      ["휠 직경 D", `≥ ${fmt(D)} m`],
      ["휠 폭 b", `≥ ${fmt(b)} m`],
      ["직경/폭 비", fmt(D / b, 1) + " : 1"],
      ["예상 침하량", fmt(r.z0 * 100, 1) + " cm"],
      ["접지 면적 (바퀴당)", fmt(r.A * 1e4, 0) + " cm²"],
      ["바퀴당 순견인력", fmt(r.T, 0) + " N"],
      ["차량 총 견인력 (4륜)", fmt(4 * r.T, 0) + " N"],
      ["차량 견인계수", fmt(4 * r.T / (state.mass * 9.81), 2)],
    ];
  }
  if (headRow) rows.unshift(headRow);
  tbl.innerHTML = rows.map(([k, v]) => `<tr><td>${k}</td><td>${v}</td></tr>`).join("");
}

/* ---------- 모드 전환 ---------- */
function setMode(mode) {
  state.mode = mode;
  document.querySelectorAll(".mode-btn").forEach(b =>
    b.classList.toggle("on", b.dataset.mode === mode));
  $("driveCtl").style.display = mode === "drive" ? "" : "none";
  $("holdCtl").style.display = mode === "hold" ? "" : "none";
  const t = MODE_TXT[mode];
  $("subTitle").textContent = t.sub;
  $("lg0").textContent = t.lg0;
  $("lg1").textContent = t.lg1;
  $("eq0t").innerHTML = t.eq0t;
  $("eq1t").innerHTML = mode === "drive"
    ? `목표 견인 달성선 — 바퀴하중의 ${Math.round(state.tgt * 100)}%`
    : `목표 안전율선 — SF = ${state.sfTgt.toFixed(1)}`;
  update();
}
document.querySelectorAll(".mode-btn").forEach(btn =>
  btn.addEventListener("click", () => setMode(btn.dataset.mode)));

/* ---------- 컨트롤 이벤트 ---------- */
function bind(id, outId, fmtFn, apply) {
  const el = $(id);
  el.addEventListener("input", () => {
    $(outId).textContent = fmtFn(+el.value);
    apply(+el.value);
    debounce();
  });
}
let timer = null;
function debounce() { clearTimeout(timer); timer = setTimeout(update, 160); }

document.querySelectorAll(".soil-btn:not(.mode-btn):not(.tread-btn)").forEach(btn => {
  btn.addEventListener("click", () => {
    document.querySelectorAll(".soil-btn:not(.mode-btn):not(.tread-btn)").forEach(b => b.classList.remove("on"));
    btn.classList.add("on");
    state.soil = btn.dataset.soil;
    $("soilBadge").textContent = SOILS[state.soil].name;
    update();
  });
});

document.querySelectorAll(".tread-btn").forEach(btn => {
  btn.addEventListener("click", () => {
    document.querySelectorAll(".tread-btn").forEach(b => b.classList.remove("on"));
    btn.classList.add("on");
    state.tread = btn.dataset.tread;
    update();
  });
});

bind("mass", "oMass", v => `${v} kg`, v => state.mass = v);
bind("slip", "oSlip", v => `${v}%`, v => state.slip = v / 100);
bind("tgt", "oTgt", v => `${v}%`, v => {
  state.tgt = v / 100;
  if (state.mode === "drive")
    $("eq1t").innerHTML = `목표 견인 달성선 — 바퀴하중의 ${v}%`;
});
bind("slope", "oSlope", v => `${v}°`, v => state.slope = v);
bind("sfTgt", "oSf", v => (v / 10).toFixed(1), v => {
  state.sfTgt = v / 10;
  if (state.mode === "hold")
    $("eq1t").innerHTML = `목표 안전율선 — SF = ${state.sfTgt.toFixed(1)}`;
});
bind("dmin", "oDmin", v => fmt(v / 100) + " m", v => state.dmin = v / 100);
bind("dmax", "oDmax", v => fmt(v / 100) + " m", v => state.dmax = v / 100);
bind("bmin", "oBmin", v => fmt(v / 100) + " m", v => state.bmin = v / 100);
bind("bmax", "oBmax", v => fmt(v / 100) + " m", v => state.bmax = v / 100);

if (window.Plotly) update();
else {
  document.getElementById("wrap").style.display = "none";
  document.getElementById("fallback").style.display = "grid";
}
