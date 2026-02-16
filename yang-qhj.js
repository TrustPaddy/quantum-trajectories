"use strict";
// ═══════════════════════════════════════════════════════════════
//  YANG QHJ DYNAMICS MODULE
//  ─────────────────────────
//  Based on: Yang, C.D. (2005)
//  "Solving Quantum Trajectories in Coulomb Potential
//   by Quantum Hamilton–Jacobi Theory"
//  Int. J. Quantum Chem. 106, 1620–1639
//
//  Implements the total potential V̄(ρ,θ) from Eq. (40):
//
//    V̄ = −2Z/ρ  +  (4 + cot²θ)/(4ρ²)
//         − d²ln R_{nl}/dρ²  − (1/ρ²) d²ln Θ_{lm}/dθ²
//
//  Forces are obtained by numerical gradient of V̄ in Cartesian
//  coordinates and converted to atomic units (Hartree/a₀).
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

  // Near poles: cot²θ → ∞, return large barrier
  if (Math.abs(st) < 1e-8) return -1e6;

  const P = assocLegendreP(l, absm, ct);

  // Near angular node: use relative threshold
  let Pm1 = 0;
  if (l - 1 >= absm) Pm1 = assocLegendreP(l - 1, absm, ct);
  const angScale = Math.max(Math.abs(Pm1), 1);
  if (Math.abs(P) < 1e-12 * angScale) return -1e6;

  // u = (dΘ/dθ) / Θ = l·cotθ − (l+|m|)·P_{l-1}^{|m|} / (sinθ · P_l^{|m|})
  let u;
  if (l - 1 >= absm) {
    u = l * ct / st - (l + absm) * Pm1 / (st * P);
  } else {
    // P_{l-1}^{|m|} = 0 when l-1 < |m|
    u = l * ct / st;
  }

  // Cap u to prevent overflow
  const U_CAP = 1e4;
  u = Math.max(-U_CAP, Math.min(U_CAP, u));

  return -ct / st * u - l * (l + 1) + absm * absm / (st * st) - u * u;
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
// Force = −∇V computed by central finite differences.
const YANG_FORCE_CAP_BASE = 200;   // base max force magnitude (Hartree/a₀)

function yangForceCartesian(x, y, z, n, l, m, Z) {
  // Scale force cap with n² — higher shells need stronger barrier forces
  const forceCap = YANG_FORCE_CAP_BASE * Math.max(1, n * n);

  const r = Math.sqrt(x * x + y * y + z * z);
  let h = Math.max(r * 1e-4, 1e-6);

  // First pass: compute gradient with initial h
  let Vxp = yangPotentialXYZ(x + h, y, z, n, l, m, Z);
  let Vxm = yangPotentialXYZ(x - h, y, z, n, l, m, Z);
  let Vyp = yangPotentialXYZ(x, y + h, z, n, l, m, Z);
  let Vym = yangPotentialXYZ(x, y - h, z, n, l, m, Z);
  let Vzp = yangPotentialXYZ(x, y, z + h, n, l, m, Z);
  let Vzm = yangPotentialXYZ(x, y, z - h, n, l, m, Z);

  // If potential variation across h is huge, refine with smaller step
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

  let fx = -(Vxp - Vxm) / (2 * h);
  let fy = -(Vyp - Vym) / (2 * h);
  let fz = -(Vzp - Vzm) / (2 * h);

  // NaN guard
  if (!Number.isFinite(fx)) fx = 0;
  if (!Number.isFinite(fy)) fy = 0;
  if (!Number.isFinite(fz)) fz = 0;

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
window.yangPotentialDimless  = yangPotentialDimless;
window.yangPotentialXYZ      = yangPotentialXYZ;
window.yangForceCartesian    = yangForceCartesian;
window.yangEquilibriumTheta  = yangEquilibriumTheta;
window.yangEquilibriumRadius = yangEquilibriumRadius;