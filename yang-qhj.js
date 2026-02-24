"use strict";
// ═══════════════════════════════════════════════════════════════
//  YANG QHJ DYNAMICS MODULE
//  ─────────────────────────
//  Based on: Yang, C.D. (2005)
//  "Solving Quantum Trajectories in Coulomb Potential
//   by Quantum Hamilton–Jacobi Theory"
//  Int. J. Quantum Chem. 106, 1620–1639
//
//  Implements the total potential V̄(ρ,θ) from equations.tex §4:
//
//    V̄ = −2Z/ρ  +  (4 + cot²θ)/(4ρ²)
//         − d²ln R_{nl}/dρ²  − (1/ρ²) d²ln Θ_{lm}/dθ²
//
//  Forces use the ANALYTICAL formulas from equations.tex §5:
//
//    f̄^r = −2/ρ²+(4+cot²θ)/(2ρ³)+d³lnR/dρ³ −(2/ρ³)d²lnΘ/dθ²
//    f̄^θ = (1/ρ²)d³lnΘ/dθ³ + cosθ/(2ρ²sin³θ)
//    f̄^φ = 0
//
//  All internal coordinates are in atomic units (a₀, Hartree).
// ═══════════════════════════════════════════════════════════════

// ─── Associated Legendre function P_l^m(x) ───────────────────
// Condon–Shortley phase convention.  m must be ≥ 0.
// For m < 0, use P_l^{-m} ∝ P_l^m (sign irrelevant for |P|).
function assocLegendreP(l, m, x) {
  m = Math.abs(m);
  if (m > l) return 0;

  // P_m^m via double-factorial recurrence
  let pmm = 1.0;
  if (m > 0) {
    const somx2 = Math.sqrt(Math.max(0, (1 - x) * (1 + x)));
    let fact = 1.0;
    for (let i = 1; i <= m; i++) {
      pmm *= -fact * somx2;
      fact += 2.0;
    }
  }
  if (l === m) return pmm;

  // P_{m+1}^m
  let pmmp1 = x * (2 * m + 1) * pmm;
  if (l === m + 1) return pmmp1;

  // Upward recurrence to P_l^m
  let pll = 0;
  for (let ll = m + 2; ll <= l; ll++) {
    pll = ((2 * ll - 1) * x * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
    pmm = pmmp1;
    pmmp1 = pll;
  }
  return pll;
}

// ─── Generalised Laguerre L_k^α(x) ──────────────────────────
function laguerreQHJ(k, alpha, x) {
  if (k === 0) return 1;
  if (k === 1) return 1 + alpha - x;
  let Lm2 = 1, Lm1 = 1 + alpha - x;
  for (let j = 2; j <= k; j++) {
    const Lj = ((2 * j - 1 + alpha - x) * Lm1 - (j - 1 + alpha) * Lm2) / j;
    Lm2 = Lm1;
    Lm1 = Lj;
  }
  return Lm1;
}

// ─── d²ln R_{nl}(ρ) / dρ²  (analytical) ─────────────────────
// R_{nl}(ρ) ∝ ρ^l · exp(−Zρ/n) · L_{n−l−1}^{2l+1}(2Zρ/n)
// Result: −l/ρ² + (2Z/n)² · (L″L − L′²) / L²
function d2lnR_yang(n, l, rho, Z) {
  if (n < 1 || l >= n || rho <= 1e-10) return 0;
  const k = n - l - 1;
  const alpha = 2 * l + 1;
  const x = 2 * Z * rho / n;
  const c = 2 * Z / n;

  const L   = laguerreQHJ(k, alpha, x);
  const Lp  = (k >= 1) ? -laguerreQHJ(k - 1, alpha + 1, x) : 0;
  const Lpp = (k >= 2) ?  laguerreQHJ(k - 2, alpha + 2, x) : 0;

  // Near a radial node: use relative threshold scaled by polynomial magnitude
  const scale = Math.max(Math.abs(Lp), Math.abs(Lpp), 1);
  if (Math.abs(L) < 1e-12 * scale) return -1e6;  // large negative → barrier

  const ratio = (Lpp * L - Lp * Lp) / (L * L);
  if (!Number.isFinite(ratio)) return -1e6;

  return -l / (rho * rho) + c * c * ratio;
}

// ─── d²ln Θ_{lm}(θ) / dθ²  (semi-analytical) ───────────────
// Uses the associated Legendre differential equation to avoid
// numerical differentiation.  Result:
//   −cotθ · u − l(l+1) + m²/sin²θ − u²
// where u = d ln|Θ|/dθ is computed from P_l^m and P_{l−1}^m.
function d2lnTheta_yang(l, m, theta) {
  if (l === 0) return 0;

  const absm = Math.abs(m);
  const ct = Math.cos(theta);
  const st = Math.sin(theta);

  // Near poles: for m≠0 the Legendre polynomial vanishes → genuine singularity.
  // For m=0 the function is finite at the poles (u→0, d²lnΘ/dθ² → finite),
  // so we let the computation proceed normally.
  if (Math.abs(st) < 1e-8 && absm > 0) return -1e6;

  const P = assocLegendreP(l, absm, ct);

  // Near angular node: use relative threshold
  let Pm1 = 0;
  if (l - 1 >= absm) Pm1 = assocLegendreP(l - 1, absm, ct);
  const angScale = Math.max(Math.abs(Pm1), 1);
  if (Math.abs(P) < 1e-12 * angScale) return -1e6;

  // u = dlnΘ/dθ = (l·cosθ·P - (l+|m|)·P_{l-1}) / (sinθ·P)
  // Written in factored form to avoid catastrophic cancellation near sinθ→0 for m=0.
  let u;
  if (l - 1 >= absm) {
    u = (l * ct * P - (l + absm) * Pm1) / (st * P);
  } else {
    // P_{l-1}^{|m|} = 0 when l-1 < |m|  (e.g. l=|m|, so l-1 < |m|)
    u = l * ct / st;
  }

  // Cap u to prevent overflow
  const U_CAP = 1e4;
  u = Math.max(-U_CAP, Math.min(U_CAP, u));

  return -ct / st * u - l * (l + 1) + absm * absm / (st * st) - u * u;
}

// ─── d³ln R_{nl}(ρ) / dρ³  (analytical) ─────────────────────
// Used in the analytical radial force f̄^r (equations.tex §5).
// d³lnR/dρ³ = 2l/ρ³ + c³ · (L'''/L − 3(L'/L)(L''/L) + 2(L'/L)³)
// where c = 2Z/n, L''' = −L_{k-3}^{α+3}  for k≥3, else 0
function d3lnR_yang(n, l, rho, Z) {
  if (n < 1 || l >= n || rho <= 1e-10) return 0;
  const k = n - l - 1;
  const alpha = 2 * l + 1;
  const x = 2 * Z * rho / n;
  const c = 2 * Z / n;

  const L    = laguerreQHJ(k, alpha, x);
  const Lp   = (k >= 1) ? -laguerreQHJ(k - 1, alpha + 1, x) : 0;
  const Lpp  = (k >= 2) ?  laguerreQHJ(k - 2, alpha + 2, x) : 0;
  const Lppp = (k >= 3) ? -laguerreQHJ(k - 3, alpha + 3, x) : 0;

  const scale = Math.max(Math.abs(Lp), Math.abs(Lpp), 1);
  if (Math.abs(L) < 1e-12 * scale) return 0;  // at node: force singular, return 0 (capped separately)

  const fL  = Lp / L;
  const gL  = Lpp / L;
  const hL  = Lppp / L;
  // d³lnL/dx³ = L'''/L − 3(L'/L)(L''/L) + 2(L'/L)³
  const fppp = hL - 3 * fL * gL + 2 * fL * fL * fL;

  const res = 2 * l / (rho * rho * rho) + c * c * c * fppp;
  return Number.isFinite(res) ? res : 0;
}

// ─── d³ln Θ_{lm}(θ) / dθ³  (analytical) ────────────────────
// equations.tex §5: f̄^θ involves d³lnΘ/dθ³.
// d³lnΘ/dθ³ = u/sin²θ + 2m²cosθ/sin³θ − (cotθ + 2u)·d²lnΘ/dθ²
// where u = dlnΘ/dθ  (same u computed in d2lnTheta_yang).
function d3lnTheta_yang(l, m, theta) {
  if (l === 0) return 0;
  const absm = Math.abs(m);
  const ct   = Math.cos(theta);
  const st   = Math.sin(theta);

  if (Math.abs(st) < 1e-8 && absm > 0) return 1e6;  // singular at poles for m≠0
  if (Math.abs(st) < 1e-8 && absm === 0) return 0;  // finite limit for m=0

  const P = assocLegendreP(l, absm, ct);
  let Pm1 = 0;
  if (l - 1 >= absm) Pm1 = assocLegendreP(l - 1, absm, ct);
  const angScale = Math.max(Math.abs(Pm1), 1);
  if (Math.abs(P) < 1e-12 * angScale) return 1e6;  // near angular node: large repulsion

  // u = dlnΘ/dθ  (stable factored form)
  let u = (l - 1 >= absm)
    ? (l * ct * P - (l + absm) * Pm1) / (st * P)
    : l * ct / st;
  const U_CAP = 1e4;
  u = Math.max(-U_CAP, Math.min(U_CAP, u));

  // f'' = d²lnΘ/dθ²
  const fpp = -ct / st * u - l * (l + 1) + absm * absm / (st * st) - u * u;

  // f''' = u/sin²θ + 2m²cosθ/sin³θ − (cotθ + 2u)·f''
  const res = u / (st * st)
    + 2 * absm * absm * ct / (st * st * st)
    - (ct / st + 2 * u) * fpp;

  return Number.isFinite(res) ? Math.max(-1e6, Math.min(1e6, res)) : 0;
}

// ─── Yang's total potential V̄(ρ,θ)  [dimensionless, Eq. 40] ─
// Units: ℏ²/(2m_e a₀²) = 0.5 Hartree
function yangPotentialDimless(rho, theta, n, l, m, Z) {
  // Regularise
  rho   = Math.max(rho, 1e-6);
  theta = Math.max(1e-6, Math.min(Math.PI - 1e-6, theta));

  const rho2 = rho * rho;
  const cotTh = Math.cos(theta) / Math.sin(theta);

  // (1) Coulomb
  const vCoulomb = -2 * Z / rho;

  // (2) Universal quantum term  (4 + cot²θ)/(4ρ²)
  const vUniversal = (4 + cotTh * cotTh) / (4 * rho2);

  // (3) Radial correction  −d²ln R/dρ²
  const radCorr = d2lnR_yang(n, l, rho, Z);

  // (4) Angular correction  −(1/ρ²) d²ln Θ/dθ²
  const angCorr = d2lnTheta_yang(l, m, theta);

  let V = vCoulomb + vUniversal - radCorr - angCorr / rho2;

  // Cap to prevent numerical overflow
  if (!Number.isFinite(V)) V = 2000;
  V = Math.max(-2000, Math.min(2000, V));

  return V;
}

// ─── Potential at Cartesian coordinates  [Hartree] ───────────
// Conversion: V_Hartree = 0.5 × V̄_dimless
function yangPotentialXYZ(x, y, z, n, l, m, Z) {
  const r = Math.sqrt(x * x + y * y + z * z);
  if (r < 1e-6) return 0.5 * yangPotentialDimless(1e-6, Math.PI / 2, n, l, m, Z);
  const theta = Math.acos(Math.max(-1, Math.min(1, z / r)));
  return 0.5 * yangPotentialDimless(r, theta, n, l, m, Z);
}

// ─── Force on electron in Cartesian coords [Hartree/a₀] ─────
// (x,y,z) = electron position relative to nucleus, in a₀.
//
// Uses the ANALYTICAL force formulas from equations.tex §5 (Yang 2005):
//
//   f̄^r = −2/ρ² + (4+cot²θ)/(2ρ³) + d³lnR/dρ³ − (2/ρ³)·d²lnΘ/dθ²
//   f̄^θ = (1/ρ²)·d³lnΘ/dθ³ + cosθ/(2ρ²sin³θ)
//   f̄^φ = 0   (no azimuthal force, conserved Lz)
//
// Physical forces [Hartree/a₀]:
//   F_r = 0.5 · f̄^r
//   F_θ = 0.5 · f̄^θ / ρ   (metric factor from -∂V/∂θ → force)
//
// Converted to Cartesian via:
//   F⃗ = F_r · r̂ + (F_θ/sinθ) · θ̂_proj
// where θ̂_proj = (xz, yz, -rxy²) / (r · rxy)  and  F_θ includes 1/r.
//
// Falls back to numerical gradient when near the z-axis (rxy < 1e-8·r)
// to avoid the coordinate singularity of spherical coordinates.

const YANG_FORCE_CAP_BASE = 200;   // base max force magnitude (Hartree/a₀)

function yangForceCartesian(x, y, z, n, l, m, Z) {
  const forceCap = YANG_FORCE_CAP_BASE * Math.max(1, n * n);

  const r2  = x * x + y * y + z * z;
  const r   = Math.sqrt(r2);
  if (r < 1e-10) return { fx: 0, fy: 0, fz: 0 };

  const rho  = r;  // in atomic units, ρ = r/a₀ = r
  const rho2 = rho * rho;
  const rho3 = rho2 * rho;

  const cosTheta = Math.max(-1, Math.min(1, z / r));
  const theta    = Math.acos(cosTheta);
  const sinTheta = Math.sin(theta);

  // rxy = projected radius in xy-plane
  const rxy = Math.sqrt(x * x + y * y);

  // ── Analytical forces in spherical coordinates ──────────────
  //
  // Check for z-axis proximity: fall back to numerical gradient
  // to avoid the 1/sinθ coordinate singularity.
  const USE_NUMERICAL = rxy < 1e-8 * r;

  let fx, fy, fz;

  if (!USE_NUMERICAL) {
    const cotTh  = cosTheta / sinTheta;
    const cot2   = cotTh * cotTh;

    // d³ of lnR and d² and d³ of lnΘ (d²lnR is only needed for the potential, not the force)
    const d3R    = d3lnR_yang(n, l, rho, Z);
    const d2Th   = d2lnTheta_yang(l, m, theta);
    const d3Th   = d3lnTheta_yang(l, m, theta);

    // f̄^r = −2/ρ² + (4+cot²θ)/(2ρ³) + d³lnR/dρ³ − (2/ρ³)·d²lnΘ/dθ²
    const fbar_r = -2 / rho2
      + (4 + cot2) / (2 * rho3)
      + d3R
      - 2 * d2Th / rho3;

    // f̄^θ = (1/ρ²)·d³lnΘ/dθ³ + cosθ/(2ρ²·sin³θ)
    const fbar_th = d3Th / rho2
      + cosTheta / (2 * rho2 * sinTheta * sinTheta * sinTheta);

    // Physical forces [Hartree/a₀]:  F_r = 0.5·f̄^r
    //   F_θ_eff = 0.5·f̄^θ / r   (the r factor comes from converting
    //             the generalised θ-force to Cartesian: F_x += F_θ·∂θ/∂x)
    const Fr      = 0.5 * fbar_r;
    const Fth_eff = 0.5 * fbar_th / r;   // [Hartree/a₀²] times the geometric factors below

    // Conversion sph → Cartesian:
    //   r̂ components: (x/r, y/r, z/r)
    //   (∂θ/∂x, ∂θ/∂y, ∂θ/∂z) = (xz, yz, −rxy²) / (r² · rxy)
    const geom = 1.0 / (r2 * rxy);   // 1/(r²·rxy)
    fx = Fr * (x / r) + Fth_eff * (x * z) * geom;
    fy = Fr * (y / r) + Fth_eff * (y * z) * geom;
    fz = Fr * (z / r) - Fth_eff * (rxy * rxy) * geom;

    // NaN guard
    if (!Number.isFinite(fx)) fx = 0;
    if (!Number.isFinite(fy)) fy = 0;
    if (!Number.isFinite(fz)) fz = 0;

  } else {
    // ── Numerical gradient fallback (near z-axis) ────────────
    let h = Math.max(r * 1e-4, 1e-6);

    let Vxp = yangPotentialXYZ(x + h, y, z, n, l, m, Z);
    let Vxm = yangPotentialXYZ(x - h, y, z, n, l, m, Z);
    let Vyp = yangPotentialXYZ(x, y + h, z, n, l, m, Z);
    let Vym = yangPotentialXYZ(x, y - h, z, n, l, m, Z);
    let Vzp = yangPotentialXYZ(x, y, z + h, n, l, m, Z);
    let Vzm = yangPotentialXYZ(x, y, z - h, n, l, m, Z);

    const maxDV = Math.max(Math.abs(Vxp - Vxm), Math.abs(Vyp - Vym), Math.abs(Vzp - Vzm));
    if (maxDV > 10) {
      h *= 0.01;
      Vxp = yangPotentialXYZ(x + h, y, z, n, l, m, Z);
      Vxm = yangPotentialXYZ(x - h, y, z, n, l, m, Z);
      Vyp = yangPotentialXYZ(x, y + h, z, n, l, m, Z);
      Vym = yangPotentialXYZ(x, y - h, z, n, l, m, Z);
      Vzp = yangPotentialXYZ(x, y, z + h, n, l, m, Z);
      Vzm = yangPotentialXYZ(x, y, z - h, n, l, m, Z);
    }

    fx = -(Vxp - Vxm) / (2 * h);
    fy = -(Vyp - Vym) / (2 * h);
    fz = -(Vzp - Vzm) / (2 * h);

    if (!Number.isFinite(fx)) fx = 0;
    if (!Number.isFinite(fy)) fy = 0;
    if (!Number.isFinite(fz)) fz = 0;
  }

  // Cap force magnitude
  const fMag = Math.sqrt(fx * fx + fy * fy + fz * fz);
  if (fMag > forceCap) {
    const s = forceCap / fMag;
    fx *= s;  fy *= s;  fz *= s;
  }

  return { fx, fy, fz };
}

// ─── Equilibrium θ from angular potential ─────────────────────
// The θ-dependent terms of V̄ scale as 1/ρ², so θ_eq is independent
// of ρ.  We minimise A(θ) = cot²θ/4 − d²lnΘ_{lm}/dθ² using a
// grid search + golden-section refinement.
//
// Known analytical results (Yang 2005):
//   l=0  or |m|=l  →  θ_eq = π/2
//   l=1, m=0       →  θ_eq = arccos √(2/3) ≈ 35.26°
//   l=2, m=1       →  θ_eq = arccos √(2/5) ≈ 50.77°
function yangEquilibriumTheta(l, m) {
  if (l === 0) return Math.PI / 2;
  const absm = Math.abs(m);
  if (absm === l) return Math.PI / 2;   // |P_l^l| ∝ sinˡθ → max at π/2

  // Angular potential (ρ-independent part)
  function angPot(theta) {
    const st = Math.sin(theta);
    if (Math.abs(st) < 1e-10) return 1e10;
    const ct = Math.cos(theta);
    return ct * ct / (4 * st * st) - d2lnTheta_yang(l, absm, theta);
  }

  // Grid search (find the well with lowest minimum)
  const N = 2000;
  let bestTheta = Math.PI / 4;
  let bestVal = angPot(bestTheta);
  for (let i = 1; i < N; i++) {
    const theta = (i / N) * Math.PI;
    const val = angPot(theta);
    if (val < bestVal) { bestVal = val; bestTheta = theta; }
  }

  // Golden-section refinement
  let lo = Math.max(0.01, bestTheta - 2 * Math.PI / N);
  let hi = Math.min(Math.PI - 0.01, bestTheta + 2 * Math.PI / N);
  for (let iter = 0; iter < 60; iter++) {
    const m1 = lo + 0.381966 * (hi - lo);
    const m2 = lo + 0.618034 * (hi - lo);
    if (angPot(m1) < angPot(m2)) hi = m2; else lo = m1;
    if (hi - lo < 1e-12) break;
  }
  return (lo + hi) / 2;
}

// ─── Equilibrium radius from Yang's potential ────────────────
// Searches for dV̄/dρ = 0 at θ = θ_eq using Newton's method.
function yangEquilibriumRadius(n, l, m, Z) {
  const theta = yangEquilibriumTheta(l, m);
  let rho = n * n / Z;  // initial guess

  for (let iter = 0; iter < 60; iter++) {
    rho = Math.max(rho, 0.05);
    const h = Math.max(rho * 1e-5, 1e-7);

    const Vp = yangPotentialDimless(rho + h, theta, n, l, m, Z);
    const Vm = yangPotentialDimless(rho - h, theta, n, l, m, Z);
    const F = -(Vp - Vm) / (2 * h);   // radial force (dimless)

    if (Math.abs(F) < 1e-9) break;

    // dF/dρ for Newton step
    const V2p = yangPotentialDimless(rho + 2 * h, theta, n, l, m, Z);
    const V2m = yangPotentialDimless(rho - 2 * h, theta, n, l, m, Z);
    const V0  = yangPotentialDimless(rho, theta, n, l, m, Z);
    const dF  = -((V2p - 2 * V0 + V2m) / (h * h)) / 2; // approximate

    if (Math.abs(dF) < 1e-18) break;
    const dr = -F / dF;
    rho += Math.max(-rho * 0.4, Math.min(rho * 0.4, dr));
  }

  return Math.max(rho, 0.3);
}

// ─── Export ──────────────────────────────────────────────────
window.assocLegendreP        = assocLegendreP;
window.laguerreQHJ           = laguerreQHJ;
window.d2lnR_yang            = d2lnR_yang;
window.d3lnR_yang            = d3lnR_yang;
window.d2lnTheta_yang        = d2lnTheta_yang;
window.d3lnTheta_yang        = d3lnTheta_yang;
window.yangPotentialDimless  = yangPotentialDimless;
window.yangPotentialXYZ      = yangPotentialXYZ;
window.yangForceCartesian    = yangForceCartesian;
window.yangEquilibriumTheta  = yangEquilibriumTheta;
window.yangEquilibriumRadius = yangEquilibriumRadius;