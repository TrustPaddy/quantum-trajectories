"use strict";

// ═════════════════════════════════════════════════════════════
//  ATOMIC UNITS (a.u.)
//  ───────────────────
//  ℏ = mₑ = e = 4πε₀ = 1
//  Length:  1 a.u. = a₀ = 0.529 177 Å
//  Mass:   1 a.u. = mₑ
//  Charge: 1 a.u. = e
//  Energy: 1 a.u. = Eₕ = 27.211 eV  (1 Ry = 0.5 Eₕ)
//  Time:   1 a.u. = ℏ/Eₕ ≈ 2.419 × 10⁻¹⁷ s
//
//  Coulomb force:  F = qi·qj / r²   (no 4πε₀ needed)
//  Kratzer term:   F = qi·qj / r³
//
//  All internal positions, velocities, accelerations in a.u.
//  Conversion to SI / scene units only at boundaries.
// ═════════════════════════════════════════════════════════════

// ─── SI constants (boundary conversions only) ────────────────
const SI = Object.freeze({
  a0:   5.29177210903e-11,
  me:   9.10938370e-31,
  mp:   1.6726219e-27,
  e:    1.602176634e-19,
  c:    2.99792458e8,
  hbar: 1.054571817e-34,
});

const MP_AU = SI.mp / SI.me;                         // proton mass in a.u. ≈ 1836.15
const SCENE_SCALE = SI.a0 * 1e12;                    // 1 a₀ → ~52.917 scene units
const C_AU = SI.c * SI.me * SI.a0 / SI.hbar;         // c in a.u. ≈ 137.036
const RAD_CONST_AU = 2 / (3 * C_AU * C_AU * C_AU);   // radiation prefactor
const V_AU = SI.hbar / (SI.me * SI.a0);              // 1 a.u. velocity in m/s

const AtomModel = Object.freeze({
  HAtom: "H Atom", H2Cation: "H2 Cation", HAnion: "H Anion",
  HECation: "HE Cation", HEAtom: "HE Atom",
});

// ─── Utilities ───────────────────────────────────────────────
const clamp = (n, lo, hi) => Math.max(lo, Math.min(hi, n));

function cartInPolar(x, y, z) {
  const r = Math.sqrt(x * x + y * y + z * z) || 0;
  return [r, r === 0 ? 0 : Math.acos(clamp(z / r, -1, 1)), Math.atan2(y, x)];
}

function createEntity(tag, attrs) {
  const el = document.createElement(tag);
  if (attrs) for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v);
  return el;
}

// ─── Reusable vector pool ────────────────────────────────────
const _v = {
  dv:   new THREE.Vector3(),
  dir:  new THREE.Vector3(),
  perp: new THREE.Vector3(),
};

// ─── Collision helper (a.u.) ─────────────────────────────────
function findNonCollidingPosition(type, generator) {
  const r1 = type === "electron" ? 0.1 : 0.27;
  for (let n = 0; n < 50; n++) {
    const c = generator();
    let ok = true;
    for (const p of STATE.particles) {
      const r2 = p.type === "electron" ? 0.1 : 0.27;
      if (p.pos.distanceTo(c) < 0.2 + r1 + r2) { ok = false; break; }
    }
    if (ok) return c;
  }
  console.warn("findNonCollidingPosition: gave up");
  return null;
}

// ─── DOM cache ───────────────────────────────────────────────
const UI = {
  Distance: document.getElementById("Distance"),
  Energy:   document.getElementById("Energy"),
  posx1: document.getElementById("posx1"), 
  y1: document.getElementById("y1"), 
  z1: document.getElementById("z1"),
  x2: document.getElementById("x2"), 
  y2: document.getElementById("y2"), 
  z2: document.getElementById("z2"),
  nLvl: document.getElementById("nLvl"), 
  lLvl: document.getElementById("lLvl"), 
  mLvl: document.getElementById("mLvl"),
  DispCountElectron: document.getElementById("DispCountElectron"),
  DispCountProton:   document.getElementById("DispCountProton"),
  Xs: document.getElementById("Xs"), 
  Ys: document.getElementById("Ys"), 
  Zs: document.getElementById("Zs"),
  led: document.getElementById("led"), 
  myButton: document.getElementById("myButton"),
  pathButton: document.getElementById("PathButton"),
  fx: document.getElementById("fx"), 
  fy: document.getElementById("fy"), 
  fz: document.getElementById("fz"),
  relaxRead: document.getElementById("relaxRead"),
};
const SCENE = { 
  el: document.querySelector("a-scene"), 
  camera: document.querySelector("#camera") 
};

// ─── Particle (all in a.u.) ─────────────────────────────────
class Particle {
  constructor(mass, charge, pos, type) {
    this.mass   = mass;     // mₑ=1, mp≈1836
    this.charge = charge;   // −1 electron, +Z proton
    this.pos    = pos;      // a₀
    this.vel    = new THREE.Vector3();
    this.acc    = new THREE.Vector3();
    this.type   = type;
  }
}

// ─── State ───────────────────────────────────────────────────
const STATE = {
  atomLabel: AtomModel.HAtom,
  particles: [], spheres: [], bohrRings: [],
  isAnimationRunning: false, isNucleusLocked: false,
  isRecordPath: false, isAngularMomentum: false,
  isLarmor: false, isThetaPhi: false, isRelaxation: false,
  _relaxTimer: 0,
  dt: 0.08,               // a.u. time step
  dtMs: 1000 / 30,        // display interval
  radiation: 1,
  N: 1, L: 0, M: 0, Z: 1,
  indexOfParticle: -1, timePassedPs: 0, countTime: 0, angleCache: 0,
  _tick: 0, _uiEvery: 3,
};

// ─── Scene management ────────────────────────────────────────
function clearScene() {
  for (let i = STATE.spheres.length - 1; i >= 0; i--) try { 
    SCENE.el.removeChild(STATE.spheres[i]); 
  } catch {}
  for (let i = STATE.bohrRings.length - 1; i >= 0; i--) try { 
    SCENE.el.removeChild(STATE.bohrRings[i]); 
  } catch {}
  STATE.particles.length = STATE.spheres.length = STATE.bohrRings.length = 0;
}

function addParticleEntity(particle, opts = {}) {
  const sphere = createEntity("a-sphere", {
    radius: opts.radius ?? (particle.type === "proton" ? 14 : 5),
    color: opts.color ?? (particle.type === "proton" ? "red" : "blue"),
    class: "clickable",
  });
  SCENE.el.appendChild(sphere); 
  STATE.spheres.push(sphere);
  const idx = STATE.spheres.length - 1;
  sphere.addEventListener("click", () => {
    STATE.indexOfParticle = idx;
    const p = STATE.particles[idx];
    if (p) { 
      UI.Xs.value = (p.pos.x * SI.a0 * 1e12).toFixed(2); 
      UI.Ys.value = (p.pos.y * SI.a0 * 1e12).toFixed(2); 
      UI.Zs.value = (p.pos.z * SI.a0 * 1e12).toFixed(2); 
    }
  });
}

function addBohrRing() {
  const r = SCENE_SCALE;
  const ring = createEntity("a-ring", { 
    "radius-inner": r * 0.99, 
    "radius-outer": r * 1.01, 
    color: "grey", 
    opacity: 0.2, 
    "ignore-ray": true }
  );
  SCENE.el.appendChild(ring); 
  STATE.bohrRings.push(ring);
}

function syncEntitiesToParticles() {
  const S = SCENE_SCALE, doUI = (STATE._tick % STATE._uiEvery === 0), pm = SI.a0 * 1e12;
  let ringIdx = 0;
  for (let i = 0; i < STATE.particles.length; i++) {
    const p = STATE.particles[i];
    const sx = p.pos.x * S, 
    sy = p.pos.y * S, 
    sz = p.pos.z * S;
    const sphere = STATE.spheres[i];
    if (sphere) sphere.object3D.position.set(sx, sy, sz);
    if (p.type === "proton") { 
      const ring = STATE.bohrRings[ringIdx++]; 
      if (ring) ring.object3D.position.set(sx, sy, sz); 
    }
    if (doUI && p.type === "electron") {
      if (i === 1) { 
        UI.posx1.textContent = (p.pos.x * pm).toFixed(3); 
        UI.y1.textContent = (p.pos.y * pm).toFixed(3);
        UI.z1.textContent = (p.pos.z * pm).toFixed(3); 
      }
      else if (i === 2) { 
        UI.x2.textContent = (p.pos.x * pm).toFixed(3); 
        UI.y2.textContent = (p.pos.y * pm).toFixed(3); 
        UI.z2.textContent = (p.pos.z * pm).toFixed(3); 
      }
    }
  }
  if (!STATE.particles[1] || STATE.particles[1].type !== "electron") 
    UI.posx1.textContent = UI.y1.textContent = UI.z1.textContent = "-";
  if (!STATE.particles[2] || STATE.particles[2].type !== "electron") 
    UI.x2.textContent = UI.y2.textContent = UI.z2.textContent = "-";
}

// ─── Atom builders (positions in a₀) ─────────────────────────
function pushProton(posAU, radius = 14)  { 
  STATE.particles.push(new Particle(MP_AU, +1, posAU.clone(), "proton"));  
  addParticleEntity(STATE.particles.at(-1), { radius, color: "red" });  
  addBohrRing(); 
}
function pushElectron(posAU, radius = 5) { 
  STATE.particles.push(new Particle(1, -1, posAU.clone(), "electron")); 
  addParticleEntity(STATE.particles.at(-1), { radius, color: "blue" }); 
}

function resetAndClear() { 
  STATE.isAnimationRunning = false; 
  UI.myButton.textContent = "Start"; 
  clearScene(); 
}

function createHAtom() { 
  resetAndClear(); 
  pushProton(new THREE.Vector3(0,0,0)); 
  pushElectron(new THREE.Vector3(-1,0,0)); 
  STATE.Z=1; 
  STATE.atomLabel=AtomModel.HAtom; 
}
function createH2Cation() { 
  resetAndClear(); 
  pushProton(new THREE.Vector3(-1,0,0)); 
  pushProton(new THREE.Vector3(1,0,0)); 
  pushElectron(new THREE.Vector3(0,0,0)); 
  STATE.Z=1; 
  STATE.atomLabel=AtomModel.H2Cation; 
}
function createHAnion() { 
  resetAndClear(); 
  pushProton(new THREE.Vector3(0,0,0)); 
  pushElectron(new THREE.Vector3(-1,0,0)); 
  pushElectron(new THREE.Vector3(1,0,0)); 
  STATE.Z=1; STATE.atomLabel=AtomModel.HAnion; 
}
function createHECation() { 
  resetAndClear(); 
  pushProton(new THREE.Vector3(0,0,0)); 
  STATE.particles[0].charge=2; 
  STATE.spheres[0].setAttribute("radius",18); 
  pushElectron(new THREE.Vector3(0.5,0,0)); 
  STATE.Z=2; 
  STATE.atomLabel=AtomModel.HECation; 
}
function createHEAtom() { 
  resetAndClear(); 
  pushProton(new THREE.Vector3(0,0,0)); 
  STATE.particles[0].charge=2; 
  STATE.spheres[0].setAttribute("radius",18); 
  pushElectron(new THREE.Vector3(0.5,0,0)); 
  pushElectron(new THREE.Vector3(-0.5,0,0)); 
  STATE.Z=2; STATE.atomLabel=AtomModel.HEAtom; 
}

createHAtom();

// ─── Energy / distance readout (all in a₀ / Ry) ─────────────
function updateDistanceAndEnergy() {
  if (STATE._tick % STATE._uiEvery !== 0) return;
  const p = STATE.particles;
  switch (STATE.atomLabel) {
    case AtomModel.HAtom: {
      if (p.length < 2) return;
      const d = p[1].pos.distanceTo(p[0].pos);
      UI.Distance.textContent = d.toFixed(3) + " a₀";
      UI.Energy.textContent = (1/(d*d) - 2/d).toFixed(3) + " Ry";
      break;
    }
    case AtomModel.H2Cation: {
      if (p.length < 3) return;
      const d1 = p[2].pos.distanceTo(p[1].pos), d2 = p[2].pos.distanceTo(p[0].pos);
      UI.Distance.textContent = d1.toFixed(3) + " a₀";
      UI.Energy.textContent = ((1/(d1*d1))-(2/d1)+(1/(d2*d2))-(2/d2)+(2*0.79473/(d1+d2))).toFixed(3) + " Ry";
      break;
    }
    case AtomModel.HAnion: {
      if (p.length < 3) return;
      const d1 = p[1].pos.distanceTo(p[0].pos), d2 = p[2].pos.distanceTo(p[0].pos);
      UI.Distance.textContent = d1.toFixed(3) + " a₀";
      UI.Energy.textContent = ((1/(d1*d1))-(2/d1)+(1/(d2*d2))-(2/d2)+(2*0.72644/(d1+d2))).toFixed(3) + " Ry";
      break;
    }
    case AtomModel.HECation: {
      if (p.length < 2) return;
      const d = p[1].pos.distanceTo(p[0].pos);
      UI.Distance.textContent = d.toFixed(3) + " a₀";
      UI.Energy.textContent = (1/(d*d) - 4/d).toFixed(3) + " Ry";
      break;
    }
    case AtomModel.HEAtom: {
      if (p.length < 3) return;
      const d = p[1].pos.distanceTo(p[0].pos);
      UI.Distance.textContent = d.toFixed(3) + " a₀";
      UI.Energy.textContent = ((2/(d*d))-(4*1.95393/d)+(1/d)).toFixed(3) + " Ry";
      break;
    }
  }
}

// ─── Controls ────────────────────────────────────────────────
function startAnimation() { 
  STATE.isAnimationRunning = !STATE.isAnimationRunning; 
  UI.myButton.textContent = STATE.isAnimationRunning ? "Stop" : "Start"; 
}
function toggleThetaPhi(on) { 
  STATE.isThetaPhi = !!on; 
}
function toggleLamor(on) { 
  STATE.isLarmor   = !!on; 
}
function changeAnsatz() { 
  STATE.isAngularMomentum = !STATE.isAngularMomentum; 
}

function changeAtomModel(v) {
  const m = { 
    HAtom: createHAtom, 
    H2Cation: createH2Cation, 
    HAnion: createHAnion, 
    HECation: createHECation, 
    HEAtom: createHEAtom 
  };
  (m[v] || (() => console.warn("Unknown:", v)))();
}

function setVelocitiesToZero() {
  for (const p of STATE.particles) { 
    p.vel.set(0,0,0); 
    p.acc.set(0,0,0); 
  }
  const toAU = v => v / V_AU;  // SI m/s → a.u. velocity
  const p = STATE.particles;
  switch (STATE.atomLabel) {
    case AtomModel.H2Cation: 
      if (p[2]) p[2].vel.x = toAU(1.14e-13); 
      break;
    case AtomModel.HAnion:   
      if (p[1]) p[1].vel.x = toAU(-6.35e-14); 
      if (p[2]) p[2].vel.x = toAU(+6.35e-14); 
      break;
    case AtomModel.HEAtom:   
      if (p[1]) p[1].vel.x = toAU(+7.35e-13); 
      if (p[2]) p[2].vel.x = toAU(-7.35e-13); 
      break;
  }
}

function lockNucleus() { 
  STATE.isNucleusLocked = !STATE.isNucleusLocked; 
  if (STATE.isNucleusLocked) 
    for (const p of STATE.particles) 
      if (p.type === "proton") p.vel.set(0,0,0); 
}
function startRecord() { 
  STATE.isRecordPath = !STATE.isRecordPath; 
  UI.pathButton.textContent = STATE.isRecordPath ? "⏸" : "⏵"; 
  UI.led.classList.toggle("recording", STATE.isRecordPath); 
}
function deletePath() { 
  SCENE.el.querySelectorAll("a-circle.trajectory-point").forEach(p => p.parentNode?.removeChild(p)); 
}

// ─── Quantum numbers & selection rules ───────────────────────
const LEVEL_CAMERA = { 
  1: { radiation: 1, z: 200 }, 
  2: { radiation: 0.999905, z: 250 }, 
  3: { radiation: 0.99995, z: 400 } 
};
function applyLevelCamera() { 
  const c = LEVEL_CAMERA[STATE.N]; 
  if (c) { STATE.radiation = c.radiation; 
    SCENE.camera.setAttribute("position", { x:0, y:0, z:c.z }); 
  } 
}
function updateQuantumUI() {
   UI.nLvl.textContent = String(STATE.N); 
   UI.lLvl.textContent = String(STATE.L); 
   if (UI.mLvl) UI.mLvl.textContent = String(STATE.M); 
}

function tryTransition(nN, nL, nM) {
  if (nN<1||nL<0||nL>=nN||Math.abs(nM)>nL) return false;
  if (Math.abs(nL-STATE.L)!==1) { 
    console.warn(`Δℓ=${nL-STATE.L} verboten`); 
    return false; 
  }
  if (Math.abs(nM-STATE.M)>1) { 
    console.warn(`Δm=${nM-STATE.M} verboten`); 
    return false; 
  }
  if (STATE.atomLabel===AtomModel.HAtom && STATE.isAngularMomentum && STATE.particles[1] && nN>STATE.N) {
    STATE.particles[1].pos.multiplyScalar(4);
    const conn = STATE.particles[0].pos.clone().sub(STATE.particles[1].pos);
    STATE.particles[1].vel.add(new THREE.Vector3(-conn.y, conn.x, 0).normalize().multiplyScalar(0.02));
  }
  STATE.N=nN; STATE.L=nL; STATE.M=nM; STATE._relaxTimer=0;
  applyLevelCamera(); updateQuantumUI(); return true;
}

function incrementEnergyLevel() { 
  tryTransition(STATE.N+1, STATE.L+1, STATE.M); 
}
function decrementEnergyLevel() { 
  tryTransition(STATE.N-1, STATE.L-1, STATE.M); 
}
function incrementAngularMomentum() { 
  tryTransition(STATE.N, STATE.L+1, STATE.M); 
}
function decrementAngularMomentum() { 
  tryTransition(STATE.N, STATE.L-1, STATE.M); 
}
function incrementM() { 
  if (!tryTransition(STATE.N, STATE.L+1, STATE.M+1)) 
    tryTransition(STATE.N, STATE.L-1, STATE.M+1); 
}
function decrementM() { 
  if (!tryTransition(STATE.N, STATE.L+1, STATE.M-1)) 
    tryTransition(STATE.N, STATE.L-1, STATE.M-1); 
}

// ─── Relaxation ──────────────────────────────────────────────
const RELAX_BASE_TICKS = 300;
function toggleRelaxation(on) { 
  STATE.isRelaxation = on !== undefined ? !!on : !STATE.isRelaxation; 
  STATE._relaxTimer = 0; 
}
function tickRelaxation() {
  if (!STATE.isRelaxation || STATE.N <= 1) return;
  STATE._relaxTimer++;
  if (STATE._relaxTimer >= Math.round(RELAX_BASE_TICKS / (STATE.N * STATE.N))) {
    if (!tryTransition(STATE.N-1, STATE.L-1, STATE.M) && STATE.L===0) STATE._relaxTimer = 0;
  }
}

// ─── Particle count controls ─────────────────────────────────
function changeCountProton(delta) {
  if (delta > 0) { 
    const pos = findNonCollidingPosition("proton", () => 
      new THREE.Vector3(Math.random()*4-2, Math.random()*4-2, 0)
    ); 
    if (pos) pushProton(pos); 
  }
  else if (delta < 0) { 
    for (let i = STATE.particles.length-1; i >= 0; i--) { 
      if (STATE.particles[i].type==="proton") { 
        const ring = STATE.bohrRings.pop(); 
        if (ring) SCENE.el.removeChild(ring); 
        const s = STATE.spheres[i]; 
        if (s) SCENE.el.removeChild(s); 
        STATE.spheres.splice(i,1); 
        STATE.particles.splice(i,1); 
        break; 
      } 
    } 
  }
  UI.DispCountProton.textContent = String(STATE.particles.filter(p => p.type==="proton").length);
}

function changeCountElectron(delta) {
  if (delta > 0) { 
    const pos = findNonCollidingPosition("electron", () => 
      new THREE.Vector3(Math.random()*2-1, Math.random()*2-1, 0)
    ); 
    if (pos) pushElectron(pos, 5); 
  }
  else if (delta < 0) { 
    for (let i = STATE.particles.length-1; i >= 0; i--) { 
      if (STATE.particles[i].type==="electron") { 
        const s=STATE.spheres[i]; 
        if(s)SCENE.el.removeChild(s); 
        STATE.spheres.splice(i,1); 
        STATE.particles.splice(i,1); 
        break; 
      } 
    } 
  }
  UI.DispCountElectron.textContent = String(STATE.particles.filter(p => p.type==="electron").length);
}

// ─── θ/φ quantum terms ───────────────────────────────────────
function d2ThetaLnPsi(n, l, theta) {
  const c2 = Math.cos(theta)**2, s2 = Math.sin(theta)**2;
  if (n===2 && l===1) return -1/(c2||1e-9);
  if (n===3 && l===2) { const d=3*c2-1; return 6*(s2-2)/((d*d)||1e-9); }
  return 0;
}
function d2PhiLnPsi() { return 0; }

// ─── Position input (pm → a.u.) ─────────────────────────────
["Xs","Ys","Zs"].forEach((key, idx) => {
  UI[key].addEventListener("keypress", ev => {
    if (ev.key !== "Enter") return;
    const pm = parseFloat(UI[key].value);
    if (!Number.isFinite(pm)) return;
    const p = STATE.particles[STATE.indexOfParticle];
    if (!p) return;
    p.pos[["x","y","z"][idx]] = (pm * 1e-12) / SI.a0;
    p.acc.set(0,0,0);
  });
});

// ═════════════════════════════════════════════════════════════
//  PHYSICS: Gather-then-Integrate
//  ─────────────────────────────
//  Phase 1: Compute acceleration for ALL particles from their
//           current positions (forces evaluated simultaneously).
//  Phase 2: Update velocity and position for ALL particles.
//
//  In atomic units with 4πε₀=1:
//    F_coulomb = −Z · qi · qj / r²
//    F_kratzer = qi · qj / r³
//    a = F / m,   v += a·dt,   x += v·dt·radiation
// ═════════════════════════════════════════════════════════════

const SOFTENING = 1e-4;  // a.u.

function animation() {
  updateDistanceAndEnergy();
  if (!STATE.isAnimationRunning || STATE.particles.length === 0) return;
  STATE._tick++;
  tickRelaxation();

  const pArr = STATE.particles, dt = STATE.dt;
  const Lterm = STATE.isAngularMomentum ? 0 : STATE.L * (STATE.L + 1);

  // ── Phase 1: Gather forces ──
  for (let i = 0; i < pArr.length; i++) {
    const pi = pArr[i];
    pi.acc.set(0, 0, 0);
    if (STATE.isNucleusLocked && pi.type === "proton") continue;

    for (let j = 0; j < pArr.length; j++) {
      if (i === j) continue;
      const pj = pArr[j];

      _v.dv.copy(pj.pos).sub(pi.pos);
      const dx = _v.dv.x, dy = _v.dv.y, dz = _v.dv.z;
      const r2 = dx*dx + dy*dy + dz*dz + SOFTENING*SOFTENING;
      const r  = Math.sqrt(r2);
      const r3 = r2 * r;

      // Coulomb
      const Fel = (-STATE.Z * pi.charge * pj.charge) / r2;

      // Kratzer + spin (opposite charges only)
      let Fk = 0, Fspin = 0;
      if (pi.charge * pj.charge < 0) {
        Fk = (pi.charge * pj.charge) / r3;
        if (Lterm) Fspin = (pi.charge * pj.charge * Lterm) / r3;
      }

      // Radiation reaction
      let Frad = 0;
      const sgn = pi.pos.x < pj.pos.x ? -1 : pi.pos.x > pj.pos.x ? +1 : 0;
      if (sgn !== 0) Frad = sgn * RAD_CONST_AU * (pi.charge*pj.charge)**2 * pi.vel.x / (r3 * pi.mass);

      // θ/φ correction
      let Ftp = 0;
      if (STATE.isThetaPhi) {
        const [, th] = cartInPolar(dx, dy, dz);
        const s2 = Math.sin(th)**2 || 1e-9;
        Ftp = 2 * (d2ThetaLnPsi(STATE.N, STATE.L, th) + d2PhiLnPsi() / s2) / (r3 || 1e-18);
      }

      // Radial acceleration (ACCUMULATED over all pairs)
      _v.dir.copy(_v.dv).divideScalar(r);
      pi.acc.addScaledVector(_v.dir, (Fel + Fk + Frad + Ftp) / pi.mass);

      // Spin term (perpendicular)
      if (Lterm) {
        _v.perp.set(-dy, dx, 0);
        const pLen = _v.perp.length();
        if (pLen > 1e-12) { _v.perp.divideScalar(pLen); pi.acc.addScaledVector(_v.perp, Fspin / pi.mass); }
      }
    }
  }

  // ── Phase 2: Integrate ──
  for (let i = 0; i < pArr.length; i++) {
    const pi = pArr[i];
    if (STATE.isNucleusLocked && pi.type === "proton") { pi.vel.set(0,0,0); continue; }
    if (STATE.isLarmor && pi.type === "electron") pi.vel.multiplyScalar(Math.exp(-1e-3 * dt));
    pi.vel.addScaledVector(pi.acc, dt);
    pi.pos.addScaledVector(pi.vel, dt * STATE.radiation);
  }

  // ── Force readout (throttled) ──
  if (STATE._tick % STATE._uiEvery === 0) {
    for (const pi of pArr) {
      if (pi.type !== "electron") continue;
      UI.fx.textContent = (pi.acc.x * pi.mass).toExponential(2) + " Eₕ/a₀";
      UI.fy.textContent = (pi.acc.y * pi.mass).toExponential(2) + " Eₕ/a₀";
      UI.fz.textContent = (pi.acc.z * pi.mass).toExponential(2) + " Eₕ/a₀";
      break;
    }
  }
}

// ─── Timer ───────────────────────────────────────────────────
let interval = setInterval(animation, STATE.dtMs);
function changeInterval(newDt) {
  const v = Number(newDt); if (!Number.isFinite(v) || v <= 0) return;
  STATE.dtMs = v; clearInterval(interval); interval = setInterval(animation, v);
}

// ─── Render loop ─────────────────────────────────────────────
function intoRealWorld() {
  if (typeof frames === "number") frames++;
  if (STATE.particles.length > 0) {
    syncEntitiesToParticles();
    if (STATE.isRecordPath && typeof frames === "number" && frames % 5 === 0) {
      const S = SCENE_SCALE;
      for (const pi of STATE.particles) {
        if (pi.type !== "electron") continue;
        const pt = createEntity("a-circle", { radius: 1, color: "grey" });
        pt.classList.add("trajectory-point");
        pt.object3D.position.set(pi.pos.x*S, pi.pos.y*S, pi.pos.z*S);
        SCENE.el.appendChild(pt);
      }
    }
  }
  requestAnimationFrame(intoRealWorld);
}
requestAnimationFrame(intoRealWorld);

// ─── Public API ──────────────────────────────────────────────
window.cartInPolar = cartInPolar;
window.startAnimation = startAnimation;
window.toggleLamor = toggleLamor;
window.toggleThetaPhi = toggleThetaPhi;
window.createH2Cation = createH2Cation;
window.createHAnion = createHAnion;
window.createHAtom = createHAtom;
window.createHECation = createHECation;
window.createHEAtom = createHEAtom;
window.changeAtomModel = changeAtomModel;
window.changeCountProton = changeCountProton;
window.changeCountElectron = changeCountElectron;
window.setVelocitiesToZero = setVelocitiesToZero;
window.lockNucleus = lockNucleus;
window.startRecord = startRecord;
window.deletePath = deletePath;
window.animation = animation;
window.changeInterval = changeInterval;
window.intoRealWorld = intoRealWorld;
window.incrementEnergyLevel = incrementEnergyLevel;
window.decrementEnergyLevel = decrementEnergyLevel;
window.incrementAngularMomentum = incrementAngularMomentum;
window.decrementAngularMomentum = decrementAngularMomentum;
window.changeAnsatz = changeAnsatz;
window.tryTransition = tryTransition;
window.toggleRelaxation = toggleRelaxation;
window.incrementM = incrementM;
window.decrementM = decrementM;