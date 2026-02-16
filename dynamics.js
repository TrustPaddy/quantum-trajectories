"use strict";
// Test comment
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
  HAtom: "H Atom", 
  H2Cation: "H2 Cation", 
  HAnion: "H Anion",
  HECation: "HE Cation", 
  HEAtom: "HE Atom",
});

// ─── Utilities ───────────────────────────────────────────────
const clamp = (n, lo, hi) => Math.max(lo, Math.min(hi, n));

function cartInPolar(x, y, z) {
  const r = Math.sqrt(x * x + y * y + z * z) || 0;
  return [r, r === 0 ? 0 : Math.acos(clamp(z / r, -1, 1)), Math.atan2(y, x)];
}

function createEntity(tag, attrs) {
  const el = document.createElement(tag);
  if (attrs) 
    for (const [k, v] of Object.entries(attrs)) 
      el.setAttribute(k, v);
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
      if (p.pos.distanceTo(c) < 0.2 + r1 + r2) { 
        ok = false; 
        break; 
      }
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
  x1: document.getElementById("x1"),  // FIX: was posx1
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
  myButton: document.getElementById("myButton"),
  pathButton: document.getElementById("PathButton"),
  fx: document.getElementById("fx"),
  fy: document.getElementById("fy"),
  fz: document.getElementById("fz"),
  relaxRead: document.getElementById("relaxRead"),
  modeRead: document.getElementById("modeRead"),
  // Selected particle display
  selectedX: document.getElementById("selectedX"),
  selectedY: document.getElementById("selectedY"),
  selectedZ: document.getElementById("selectedZ"),
  // Velocity display
  vx: document.getElementById("vx"),
  vy: document.getElementById("vy"),
  vz: document.getElementById("vz"),
  vMag: document.getElementById("vMag"),
  // Panel-specific readouts
  dtPanelRead: document.getElementById("dtPanelRead"),
  fpsPanelRead: document.getElementById("fpsPanelRead"),
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
  particles: [], 
  spheres: [], 
  bohrRings: [],
  isAnimationRunning: false, 
  isNucleusLocked: false,
  isRecordPath: false, 
  isAngularMomentum: false,
  isLarmor: false, 
  isThetaPhi: false, 
  isRelaxation: false,
  dynamicsMode: 'classic',  // 'classic' or 'qhj'
  _relaxTimer: 0,
  dt: 0.08,               // a.u. time step
  dtMs: 1000 / 30,        // display interval
  radiation: 1,
  N: 1, L: 0, M: 0, Z: 1,
  indexOfParticle: -1, 
  timePassedPs: 0, 
  countTime: 0, 
  angleCache: 0,
  _tick: 0, _uiEvery: 3,
};

// ─── Scene management ────────────────────────────────────────
function clearScene() {
  for (let i = STATE.spheres.length - 1; i >= 0; i--) 
    try { SCENE.el.removeChild(STATE.spheres[i]); 
  } catch {}
  for (let i = STATE.bohrRings.length - 1; i >= 0; i--) 
    try { SCENE.el.removeChild(STATE.bohrRings[i]); 
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
}

function addBohrRing() {
  const r = SCENE_SCALE;
  const ring = createEntity("a-ring", { 
    "radius-inner": r * 0.99, 
    "radius-outer": r * 1.01, 
    color: "grey", 
    opacity: 0.2, 
    "ignore-ray": true 
  });
  SCENE.el.appendChild(ring); 
  STATE.bohrRings.push(ring);
}

function syncEntitiesToParticles() {
  const S = SCENE_SCALE, pm = SI.a0 * 1e12;
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
    if (p.type === "electron") {
      if (i === 1) { 
        UI.x1.textContent = (p.pos.x * pm).toFixed(3); 
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
    UI.x1.textContent = UI.y1.textContent = UI.z1.textContent = "-";
  if (!STATE.particles[2] || STATE.particles[2].type !== "electron") 
    UI.x2.textContent = UI.y2.textContent = UI.z2.textContent = "-";
}

function updateSelectedParticleDisplay() {
  const idx = STATE.indexOfParticle;
  if (idx < 0 || !STATE.particles[idx]) {
    if (UI.selectedX) UI.selectedX.textContent = "–";
    if (UI.selectedY) UI.selectedY.textContent = "–";
    if (UI.selectedZ) UI.selectedZ.textContent = "–";
    if (UI.vx) UI.vx.textContent = "0";
    if (UI.vy) UI.vy.textContent = "0";
    if (UI.vz) UI.vz.textContent = "0";
    if (UI.vMag) UI.vMag.textContent = "0";
    return;
  }

  const p = STATE.particles[idx];
  const pm = SI.a0 * 1e12;

  // Position
  if (UI.selectedX) UI.selectedX.textContent = (p.pos.x * pm).toFixed(3);
  if (UI.selectedY) UI.selectedY.textContent = (p.pos.y * pm).toFixed(3);
  if (UI.selectedZ) UI.selectedZ.textContent = (p.pos.z * pm).toFixed(3);

  // Velocity
  if (UI.vx) UI.vx.textContent = (p.vel.x * V_AU).toExponential(2);
  if (UI.vy) UI.vy.textContent = (p.vel.y * V_AU).toExponential(2);
  if (UI.vz) UI.vz.textContent = (p.vel.z * V_AU).toExponential(2);

  const vMag = Math.sqrt(p.vel.x**2 + p.vel.y**2 + p.vel.z**2);
  if (UI.vMag) UI.vMag.textContent = (vMag * V_AU).toExponential(2);
}

// ─── Drag & Drop ─────────────────────────────────────────────
// Click+hold on a particle → drag it in the camera-facing plane.
// While dragging, velocity is zeroed so it stays put on release.
// Works by raycasting from the mouse into the scene, then
// projecting onto a plane through the particle position that
// faces the camera.

const _drag = {
  active: false,
  particleIdx: -1,
  plane: new THREE.Plane(),
  intersect: new THREE.Vector3(),
  raycaster: new THREE.Raycaster(),
  mouse: new THREE.Vector2(),
  offset: new THREE.Vector3(),  // grab offset so particle doesn't jump to cursor
};

function _getCanvasEl() {
  return SCENE.el?.canvas || SCENE.el?.querySelector('canvas');
}

function _updateMouse(ev) {
  const canvas = _getCanvasEl();
  if (!canvas) return;
  const rect = canvas.getBoundingClientRect();
  _drag.mouse.x = ((ev.clientX - rect.left) / rect.width) * 2 - 1;
  _drag.mouse.y = -((ev.clientY - rect.top) / rect.height) * 2 + 1;
}

function _getCameraObj() {
  return SCENE.camera?.object3D?.children?.find(c => c.isCamera)
    || SCENE.camera?.object3D;
}

function _onPointerDown(ev) {
  _updateMouse(ev);
  const cam = _getCameraObj();
  if (!cam) return;

  _drag.raycaster.setFromCamera(_drag.mouse, cam);

  // Test against all particle spheres
  const meshes = [];
  for (let i = 0; i < STATE.spheres.length; i++) {
    const obj = STATE.spheres[i]?.object3D;
    if (obj) {
      obj.traverse(child => { 
        if (child.isMesh) { 
          child._particleIdx = i; 
          meshes.push(child); 
        } 
      });
    }
  }

  const hits = _drag.raycaster.intersectObjects(meshes, false);
  if (hits.length === 0) return;

  const hit = hits[0];
  const idx = hit.object._particleIdx;
  if (idx == null || !STATE.particles[idx]) return;

  _drag.active = true;
  _drag.particleIdx = idx;
  STATE.indexOfParticle = idx;

  // Update position input fields immediately
  const particle = STATE.particles[idx];
  const pm = SI.a0 * 1e12;
  UI.Xs.value = (particle.pos.x * pm).toFixed(2);
  UI.Ys.value = (particle.pos.y * pm).toFixed(2);
  UI.Zs.value = (particle.pos.z * pm).toFixed(2);

  // Add dragging cursor class
  const sceneEl = SCENE.el;
  if (sceneEl) sceneEl.classList.add('dragging-particle');

  // Set up drag plane: faces camera, passes through particle
  const worldPos = new THREE.Vector3(
    particle.pos.x * SCENE_SCALE,
    particle.pos.y * SCENE_SCALE,
    particle.pos.z * SCENE_SCALE,
  );

  const camDir = new THREE.Vector3();
  cam.getWorldDirection(camDir);
  _drag.plane.setFromNormalAndCoplanarPoint(camDir, worldPos);

  // Calculate grab offset (so particle doesn't snap to cursor)
  if (_drag.raycaster.ray.intersectPlane(_drag.plane, _drag.intersect)) {
    _drag.offset.copy(worldPos).sub(_drag.intersect);
  }

  // Disable look-controls while dragging
  const lookControls = SCENE.camera?.getAttribute?.('look-controls');
  if (lookControls !== null) SCENE.camera.setAttribute('look-controls', 'enabled', false);

  // Visual highlight: make dragged particle emissive
  const sphere = STATE.spheres[idx];
  if (sphere) sphere.setAttribute('material', 'emissive', '#446');

  ev.preventDefault();
}

function _onPointerMove(ev) {
  _updateMouse(ev);

  // Hover feedback: show grab cursor when hovering over particles
  if (!_drag.active) {
    const cam = _getCameraObj();
    const sceneEl = SCENE.el;
    if (cam && sceneEl) {
      _drag.raycaster.setFromCamera(_drag.mouse, cam);
      const meshes = [];
      for (let i = 0; i < STATE.spheres.length; i++) {
        const obj = STATE.spheres[i]?.object3D;
        if (obj) {
          obj.traverse(child => {
            if (child.isMesh) meshes.push(child);
          });
        }
      }
      const hits = _drag.raycaster.intersectObjects(meshes, false);
      if (hits.length > 0) {
        sceneEl.classList.add('hovering-particle');
      } else {
        sceneEl.classList.remove('hovering-particle');
      }
    }
    return;
  }

  const cam = _getCameraObj();
  if (!cam) return;

  _drag.raycaster.setFromCamera(_drag.mouse, cam);

  if (_drag.raycaster.ray.intersectPlane(_drag.plane, _drag.intersect)) {
    // Apply offset so particle stays where you grabbed it
    _drag.intersect.add(_drag.offset);

    // Convert scene units → a.u.
    const particle = STATE.particles[_drag.particleIdx];
    if (particle) {
      particle.pos.set(
        _drag.intersect.x / SCENE_SCALE,
        _drag.intersect.y / SCENE_SCALE,
        _drag.intersect.z / SCENE_SCALE,
      );
      // Zero velocity while dragging (so it doesn't fly away on release)
      particle.vel.set(0, 0, 0);
      particle.acc.set(0, 0, 0);

      // Update position input fields
      const pm = SI.a0 * 1e12;
      UI.Xs.value = (particle.pos.x * pm).toFixed(2);
      UI.Ys.value = (particle.pos.y * pm).toFixed(2);
      UI.Zs.value = (particle.pos.z * pm).toFixed(2);
    }
  }

  ev.preventDefault();
}

function _onPointerUp(ev) {
  if (!_drag.active) return;

  // Remove highlight
  const sphere = STATE.spheres[_drag.particleIdx];
  if (sphere) sphere.setAttribute('material', 'emissive', '#000');

  _drag.active = false;

  // Remove dragging cursor class
  const sceneEl = SCENE.el;
  if (sceneEl) sceneEl.classList.remove('dragging-particle');

  // Re-enable look-controls
  const lookControls = SCENE.camera?.getAttribute?.('look-controls');
  if (lookControls !== null) SCENE.camera.setAttribute('look-controls', 'enabled', true);
}

// Attach listeners once A-Frame scene is ready
function _initDragDrop() {
  const canvas = _getCanvasEl();
  if (!canvas) {
    // A-Frame not ready yet, retry
    setTimeout(_initDragDrop, 200);
    return;
  }
  canvas.addEventListener('pointerdown', _onPointerDown);
  canvas.addEventListener('pointermove', _onPointerMove);
  canvas.addEventListener('pointerup',   _onPointerUp);
  canvas.addEventListener('pointerleave', _onPointerUp);
  canvas.style.touchAction = 'none';  // prevent scroll on touch devices
}
_initDragDrop();

// ─── Atom builders (positions in a₀) ─────────────────────────
function pushProton(posAU, radius = 14) { 
  STATE.particles.push(new Particle(MP_AU, +1, posAU.clone(), "proton"));  
  addParticleEntity(STATE.particles.at(-1), { radius, color: "red" });  
  addBohrRing(); 
}
function pushElectron(posAU, radius = 5) { 
  STATE.particles.push(new Particle(1,-1, posAU.clone(), "electron")); 
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
  STATE.Z=1; 
  STATE.atomLabel=AtomModel.HAnion; 
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
      if (p.type === "proton") 
        p.vel.set(0,0,0); }
function startRecord() {
  STATE.isRecordPath = !STATE.isRecordPath;
  if (STATE.isRecordPath) {
    UI.pathButton.textContent = "⏹";
    UI.pathButton.classList.add("recording");
  } else {
    UI.pathButton.textContent = "▶";
    UI.pathButton.classList.remove("recording");
  }
}
// deletePath() is defined below in the Trajectory section

// ─── Quantum numbers & selection rules ───────────────────────
const LEVEL_CAMERA = { 
  1: { radiation: 1, z: 200 }, 
  2: { radiation: 0.999905, z: 350 }, 
  3: { radiation: 0.99995, z: 600 }, 
  4: { radiation: 0.99998, z: 1000 }, 
  5: { radiation: 0.99999, z: 1600 }, 
  6: { radiation: 0.999995, z: 2200 } 
};
function applyLevelCamera() { 
  const c = LEVEL_CAMERA[STATE.N]; 
  if (c) { 
    STATE.radiation = c.radiation; 
    SCENE.camera.setAttribute("position", { x:0, y:0, z:c.z }); 
  } 
}
function updateQuantumUI() { 
  UI.nLvl.textContent = String(STATE.N); 
  UI.lLvl.textContent = String(STATE.L); 
  if (UI.mLvl) UI.mLvl.textContent = String(STATE.M); 
}

function tryTransition(nN, nL, nM) {
  if (nN<1||nL<0||nL>=nN||Math.abs(nM)>nL) {
    if (typeof window.showToast === 'function') {
      window.showToast(`Invalid quantum numbers: ℓ must be < n, |m| ≤ ℓ`, 'error');
    }
    return false;
  }
  if (Math.abs(nL-STATE.L)!==1) {
    console.warn(`Δℓ=${nL-STATE.L} verboten`);
    if (typeof window.showToast === 'function') {
      window.showToast(`Selection rule violated: Δℓ must be ±1 (current Δℓ = ${nL-STATE.L})`, 'warning');
    }
    return false;
  }
  if (Math.abs(nM-STATE.M)>1) {
    console.warn(`Δm=${nM-STATE.M} verboten`);
    if (typeof window.showToast === 'function') {
      window.showToast(`Selection rule violated: Δm must be ≤ ±1 (current Δm = ${nM-STATE.M})`, 'warning');
    }
    return false;
  }
  // Rescale electron to new equilibrium radius
  STATE.N=nN; STATE.L=nL; STATE.M=nM; STATE._relaxTimer=0;
  if (STATE.atomLabel === AtomModel.HAtom && STATE.particles[1]) {
    const e = STATE.particles[1];
    const nuc = STATE.particles[0];

    // Choose radius based on dynamics mode
    const newR = (STATE.dynamicsMode === 'qhj' && window.yangEquilibriumRadius)
      ? window.yangEquilibriumRadius(nN, nL, nM, STATE.Z)
      : bohrRadius(nN, nL, STATE.Z);

    // Set velocity depending on mode
    e.vel.set(0, 0, 0);
    if (STATE.dynamicsMode === 'qhj') {
      // QHJ mode: place electron at correct (ρ_eq, θ_eq) in spherical coords
      const thetaEq = window.yangEquilibriumTheta
        ? window.yangEquilibriumTheta(nL, nM)
        : Math.PI / 2;

      // Preserve φ from current position
      const oldDir = e.pos.clone().sub(nuc.pos);
      const phi = Math.atan2(oldDir.y, oldDir.x);

      const sinTh = Math.sin(thetaEq);
      const cosTh = Math.cos(thetaEq);

      e.pos.set(
        nuc.pos.x + newR * sinTh * Math.cos(phi),
        nuc.pos.y + newR * sinTh * Math.sin(phi),
        nuc.pos.z + newR * cosTh
      );

      if (nM !== 0) {
        // v_φ = |m| / (ρ sinθ) from Yang Eq. (43c)
        const vCirc = Math.abs(nM) / (newR * sinTh);
        // φ-hat = (−sinφ, cosφ, 0)
        e.vel.set(
          -Math.sign(nM) * vCirc * Math.sin(phi),
           Math.sign(nM) * vCirc * Math.cos(phi),
          0
        );
      } else if (nL > 0) {
        // m=0: electron is stationary at equilibrium (Yang Eq. 43c: dφ/dτ=0).
        // Add small radial kick so oscillation around equilibrium is visible.
        // Scale kick inversely with n to avoid overshooting for higher shells.
        const radHat = new THREE.Vector3(
          sinTh * Math.cos(phi),
          sinTh * Math.sin(phi),
          cosTh
        );
        const kickFraction = 0.08 / nN;
        e.pos.addScaledVector(radHat, kickFraction * newR);
      }
    } else {
      // Classic mode: scale position along existing direction
      const dir = e.pos.clone().sub(nuc.pos);
      const dirLen = dir.length() || 1e-9;
      dir.divideScalar(dirLen);
      e.pos.copy(nuc.pos).addScaledVector(dir, newR);

      if (nL > 0 && !STATE.isThetaPhi) {
        const vCirc = Math.sqrt(nL * (nL + 1)) / newR;
        const perp = new THREE.Vector3(-dir.y, dir.x, 0);
        const pLen = perp.length();
        if (pLen > 1e-9) { perp.divideScalar(pLen); e.vel.addScaledVector(perp, vCirc); }
      }
    }
    e.acc.set(0, 0, 0);
  }
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
  STATE.isRelaxation = on !== undefined ? 
    !!on : 
    !STATE.isRelaxation; STATE._relaxTimer = 0; 
}
function tickRelaxation() {
  if (!STATE.isRelaxation || STATE.N <= 1) return;
  STATE._relaxTimer++;
  if (STATE._relaxTimer >= Math.round(RELAX_BASE_TICKS / (STATE.N * STATE.N))) {
    if (!tryTransition(STATE.N-1, STATE.L-1, STATE.M) && STATE.L===0) 
      STATE._relaxTimer = 0;
  }
}

// ─── QHJ Mode Toggle ───────────────────────────────────────────
function toggleQHJ(on) {
  STATE.dynamicsMode = (on !== undefined ? !!on : STATE.dynamicsMode !== 'qhj')
    ? 'qhj' : 'classic';
  return STATE.dynamicsMode;
}

// ─── QHJ Force Gathering (Yang 2005, Eq. 40–41) ───────────────
// Electron–proton pairs use the full QHJ total potential.
// All other pairs use classical Coulomb.
function gatherForcesQHJ() {
  const pArr = STATE.particles;

  for (let i = 0; i < pArr.length; i++) {
    const pi = pArr[i];
    pi.acc.set(0, 0, 0);
    if (STATE.isNucleusLocked && pi.type === "proton") continue;

    for (let j = 0; j < pArr.length; j++) {
      if (i === j) continue;
      const pj = pArr[j];

      _v.dv.copy(pj.pos).sub(pi.pos);
      const dx = _v.dv.x, dy = _v.dv.y, dz = _v.dv.z;
      const r2 = dx * dx + dy * dy + dz * dz + SOFTENING * SOFTENING;
      const r  = Math.sqrt(r2);

      if (pi.type === "electron" && pj.type === "proton" && window.yangForceCartesian) {
        // ── QHJ: electron in proton's quantum potential ──
        const rx = pi.pos.x - pj.pos.x;
        const ry = pi.pos.y - pj.pos.y;
        const rz = pi.pos.z - pj.pos.z;
        const Z = Math.abs(pj.charge);
        const f = window.yangForceCartesian(rx, ry, rz, STATE.N, STATE.L, STATE.M, Z);
        pi.acc.x += f.fx;   // m_e = 1 in a.u.
        pi.acc.y += f.fy;
        pi.acc.z += f.fz;
      } else {
        // ── Classical Coulomb for proton–electron (on proton),
        //    electron–electron, and proton–proton pairs ──
        const Fel = (-pi.charge * pj.charge) / r2;
        _v.dir.copy(_v.dv).divideScalar(r);
        pi.acc.addScaledVector(_v.dir, Fel / pi.mass);
      }
    }
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
      if (STATE.particles[i].type === "proton") { 
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
      if (STATE.particles[i].type === "electron") { 
        const s = STATE.spheres[i]; 
        if (s) SCENE.el.removeChild(s); 
        STATE.spheres.splice(i,1); 
        STATE.particles.splice(i,1); 
        break; 
      } 
    } 
  }
  UI.DispCountElectron.textContent = String(STATE.particles.filter(p => p.type==="electron").length);
}

// ─── Radial quantum correction: d²_r ln R_{nl}(r) ────────────
// Generalized Laguerre polynomial L^alpha_k(x) via recurrence
function laguerreL(k, alpha, x) {
  if (k === 0) return 1;
  if (k === 1) return 1 + alpha - x;
  let Lm2 = 1, Lm1 = 1 + alpha - x;
  for (let j = 2; j <= k; j++) {
    const Lj = ((2*j - 1 + alpha - x) * Lm1 - (j - 1 + alpha) * Lm2) / j;
    Lm2 = Lm1; Lm1 = Lj;
  }
  return Lm1;
}

// Second derivative of ln(R_{nl}) w.r.t. r
// R_{nl}(r) ~ rho^l * L^{2l+1}_{n-l-1}(rho) * exp(-rho/2),  rho = 2Zr/n
// d2lnR/dr2 = -l/r^2 + (2Z/n)^2 * (L''*L - L'^2) / L^2
// where L' = -L^{alpha+1}_{k-1}, L'' = L^{alpha+2}_{k-2}
function d2RadialLnR(n, l, r, Z) {
  if (n < 1 || l >= n || r <= 0) return 0;
  const k = n - l - 1;        // degree of Laguerre polynomial
  const alpha = 2 * l + 1;
  const rho = 2 * Z * r / n;
  const c = 2 * Z / n;        // drho/dr

  const L = laguerreL(k, alpha, rho);
  if (Math.abs(L) < 1e-15) return -l / (r * r);  // near a node of R

  const Lp  = (k >= 1) ? -laguerreL(k - 1, alpha + 1, rho) : 0;  // dL/drho
  const Lpp = (k >= 2) ?  laguerreL(k - 2, alpha + 2, rho) : 0;  // d2L/drho2

  return -l / (r * r) + c * c * (Lpp * L - Lp * Lp) / (L * L);
}

// Force from radial correction (code convention: positive = inward)
// V_rad = -(1/2) * d2RadialLnR(r)  (Hartree)
// F_inward = dV/dr = -(1/2) * d/dr[d2RadialLnR]  (numerical derivative)
function radialCorrectionForce(n, l, r, Z) {
  if (n <= 1 && l === 0) return 0;  // (1,0,0) has no correction
  const h = Math.max(r * 1e-4, 1e-6);
  const dp = d2RadialLnR(n, l, r + h, Z);
  const dm = d2RadialLnR(n, l, r - h, Z);
  return -0.5 * (dp - dm) / (2 * h);
}

// Most-probable radius for state (n,l) - used for transition rescaling
// Solves dV_eff/dr = 0 numerically (Newton's method on total radial force)
function bohrRadius(n, l, Z) {
  let r = n * n / Z;
  for (let iter = 0; iter < 40; iter++) {
    if (r < 0.1) r = 0.1;
    const r2 = r * r, r3 = r2 * r;
    // Total radial force: Coulomb + Kratzer + centrifugal + radial correction
    const Fc = Z / r2;                      // Coulomb inward
    const Fk = -1 / r3;                     // Kratzer outward
    const Fcent = -l * (l + 1) / r3;        // centrifugal outward
    const Frad = radialCorrectionForce(n, l, r, Z);
    const F = Fc + Fk + Fcent + Frad;
    // Numerical derivative of F for Newton step
    const h = r * 1e-4;
    const Fp = (() => { const rr=r+h, rr2=rr*rr, rr3=rr2*rr;
      return Z/rr2 - 1/rr3 - l*(l+1)/rr3 + radialCorrectionForce(n,l,rr,Z); })();
    const Fm = (() => { const rr=r-h, rr2=rr*rr, rr3=rr2*rr;
      return Z/rr2 - 1/rr3 - l*(l+1)/rr3 + radialCorrectionForce(n,l,rr,Z); })();
    const dF = (Fp - Fm) / (2 * h);
    if (Math.abs(dF) < 1e-20) break;
    const dr = -F / dF;
    r += clamp(dr, -r * 0.5, r * 0.5);
    if (Math.abs(F) < 1e-10) break;
  }
  return Math.max(r, 0.5);
}

// ─── θ/φ quantum terms ───────────────────────────────────────
// d2_theta ln(Psi) depends only on l (and m), not on n.
// Formulas from QHM Total-potential table (m = 0 only).
function d2ThetaLnPsi(n, l, theta) {
  if (l === 0) return 0;

  const c  = Math.cos(theta);
  const s  = Math.sin(theta);
  const c2 = c * c,  s2 = s * s;
  const c4 = c2*c2,  s4 = s2*s2;
  const c6 = c4*c2,  s6 = s4*s2;
  const EPS = 1e-12;

  if (l === 1) {
    // -1 / cos^2(theta)
    return -1 / (c2 || EPS);
  }
  if (l === 2) {
    // 6(sin^2 - 2) / (3cos^2 - 1)^2
    const d = 3*c2 - 1;
    return 6 * (s2 - 2) / (d*d || EPS);
  }
  if (l === 3) {
    // 3(-25sin^6 + 70sin^4 - 65sin^2 - 25cos^6 + 17)
    //   / ((5sin^2 - 2)^2 * cos^2)
    const num = -25*s6 + 70*s4 - 65*s2 - 25*c6 + 17;
    const d   = 5*s2 - 2;
    const den = d*d * c2;
    return 3 * num / (den || EPS);
  }
  if (l === 4) {
    // 20(35sin^6 - 84sin^4 + 72sin^2 - 32)
    //   / (35cos^4 - 30cos^2 + 3)^2
    const num = 35*s6 - 84*s4 + 72*s2 - 32;
    const d   = 35*c4 - 30*c2 + 3;
    return 20 * num / (d*d || EPS);
  }
  if (l === 5) {
    // 15(-350sin^6 + 980sin^4 - 910sin^2 - 147cos^8 - 182cos^6 + 265)
    //   / ((63sin^4 - 56sin^2 + 8)^2 * cos^2)
    const c8  = c4 * c4;
    const num = -350*s6 + 980*s4 - 910*s2 - 147*c8 - 182*c6 + 265;
    const d   = 63*s4 - 56*s2 + 8;
    const den = d*d * c2;
    return 15 * num / (den || EPS);
  }

  return 0;   // l > 5 not tabulated
}
function d2PhiLnPsi() { return 0; }   // correct for all m = 0 states

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
  if (!STATE.isAnimationRunning || STATE.particles.length === 0) return;
  STATE._tick++;
  tickRelaxation();

  const pArr = STATE.particles, dt = STATE.dt;

  // ── Phase 1: Gather forces (mode-dependent) ──
  if (STATE.dynamicsMode === 'qhj') {
    gatherForcesQHJ();
  } else {
    // ── Classic force gathering ──
    const Lterm = STATE.isAngularMomentum ? 0 : STATE.L * (STATE.L + 1);

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

        // Radial quantum correction (always active for n > 1)
        let Frad_q = 0;
        if (pi.charge * pj.charge < 0 && (STATE.N > 1 || STATE.L > 0)) {
          Frad_q = radialCorrectionForce(STATE.N, STATE.L, r, STATE.Z);
        }

        // θ/φ correction (with force clamping for stability)
        let Ftp = 0;
        if (STATE.isThetaPhi && pi.charge * pj.charge < 0) {
          const [, th] = cartInPolar(dx, dy, dz);
          const s2 = Math.sin(th)**2 || 1e-9;
          const raw = d2ThetaLnPsi(STATE.N, STATE.L, th) + d2PhiLnPsi() / s2;
          const maxAng = 50;
          const clamped = Math.max(-maxAng, Math.min(maxAng, raw));
          Ftp = 2 * clamped / (r3 || 1e-18);
        }

        // Radial acceleration (ACCUMULATED over all pairs)
        _v.dir.copy(_v.dv).divideScalar(r);
        pi.acc.addScaledVector(_v.dir, (Fel + Fk + Frad + Frad_q + Ftp) / pi.mass);

        // Spin term (perpendicular)
        if (Lterm) {
          _v.perp.set(-dy, dx, 0);
          const pLen = _v.perp.length();
          if (pLen > 1e-12) { 
            _v.perp.divideScalar(pLen); pi.acc.addScaledVector(_v.perp, Fspin / pi.mass); 
          }
        }
      }
    }
  }

  // ── Phase 2: Integrate ──
  for (let i = 0; i < pArr.length; i++) {
    const pi = pArr[i];
    if (STATE.isNucleusLocked && pi.type === "proton") { 
      pi.vel.set(0,0,0); 
      continue; 
    }
    if (STATE.isLarmor && pi.type === "electron") pi.vel.multiplyScalar(Math.exp(-1e-3 * dt));
    pi.vel.addScaledVector(pi.acc, dt);
    pi.pos.addScaledVector(pi.vel, dt * STATE.radiation);
    // NaN recovery: reset to safe position if numerical blowup
    if (!Number.isFinite(pi.pos.x + pi.pos.y + pi.pos.z)) {
      console.warn("NaN detected, resetting particle", i);
      const safeR = bohrRadius(STATE.N, STATE.L, STATE.Z);
      pi.pos.set(-safeR, 0, 0);
      pi.vel.set(0, 0, 0);
      pi.acc.set(0, 0, 0);
    }
  }
}

// ═════════════════════════════════════════════════════════════
//  MAIN LOOP: rAF-driven with sub-stepping
//  ─────────────────────────────────────────
//  A single requestAnimationFrame loop handles both physics
//  and rendering.  Each frame:
//    1. Compute how many physics sub-steps fit in dtMs
//    2. Run animation() that many times (stable, no drift)
//    3. Sync 3D entities & record trajectory
//
//  Benefits over setInterval:
//  - No timer drift (rAF is V-synced)
//  - Physics & render always in sync
//  - Sub-stepping: lower dtMs → more steps/frame → smoother
//  - Pauses automatically when tab is hidden (saves CPU)
// ═════════════════════════════════════════════════════════════

let _lastFrameTime = 0;
let _physicsAccumulator = 0;
let _fpsLastTime = 0;
let _fpsFrameCount = 0;

function changeInterval(newDt) {
  const v = Number(newDt);
  if (!Number.isFinite(v) || v <= 0) return;
  STATE.dtMs = v;

  // Update Panel dt display
  if (UI.dtPanelRead) UI.dtPanelRead.textContent = v + ' ms';
}

// ─── Trajectory: InstancedMesh ───────────────────────────────
// Instead of creating hundreds of individual <a-circle> DOM
// nodes, we use a single THREE.InstancedMesh with a shared
// SphereGeometry.  This is orders of magnitude faster for
// large trajectories (1000+ points vs 1000+ DOM elements).

const TRAJ_MAX_POINTS = 4096;    // max trajectory points
const TRAJ_RADIUS = 1;       // scene units

const _trajState = {
  mesh: null,                  // THREE.InstancedMesh (created lazily)
  count: 0,                     // current number of active points
  dummy: new THREE.Object3D(),  // reusable for setMatrixAt
};

/** Get or create the InstancedMesh, attached to the A-Frame scene's three.js object */
function getTrajMesh() {
  if (_trajState.mesh) return _trajState.mesh;

  // Wait for A-Frame scene to be ready
  const sceneObj = SCENE.el?.object3D;
  if (!sceneObj) return null;

  const geo = new THREE.SphereGeometry(TRAJ_RADIUS, 4, 4);  // low-poly = fast
  const mat = new THREE.MeshBasicMaterial({ color: 0x888888, transparent: true, opacity: 0.6 });
  const mesh = new THREE.InstancedMesh(geo, mat, TRAJ_MAX_POINTS);
  mesh.count = 0;         // start with 0 visible instances
  mesh.frustumCulled = false;
  sceneObj.add(mesh);

  _trajState.mesh = mesh;
  return mesh;
}

function addTrajPoint(x, y, z) {
  const mesh = getTrajMesh();
  if (!mesh) return;

  if (_trajState.count >= TRAJ_MAX_POINTS) {
    // Ring buffer: overwrite oldest point
    _trajState.count = 0;
  }

  _trajState.dummy.position.set(x, y, z);
  _trajState.dummy.updateMatrix();
  mesh.setMatrixAt(_trajState.count, _trajState.dummy.matrix);
  _trajState.count++;
  mesh.count = Math.max(mesh.count, _trajState.count);
  mesh.instanceMatrix.needsUpdate = true;
}

function deletePath() {
  // Clear InstancedMesh
  if (_trajState.mesh) {
    _trajState.mesh.count = 0;
    _trajState.count = 0;
    _trajState.mesh.instanceMatrix.needsUpdate = true;
  }
  // Also remove any legacy <a-circle> nodes
  SCENE.el.querySelectorAll("a-circle.trajectory-point")
    .forEach(p => p.parentNode?.removeChild(p));
}

// ─── Unified rAF loop ────────────────────────────────────────
function mainLoop(timestamp) {
  if (typeof frames === "number") frames++;

  // FPS calculation (every second)
  _fpsFrameCount++;
  if (_fpsLastTime === 0) _fpsLastTime = timestamp;
  const fpsElapsed = timestamp - _fpsLastTime;

  if (fpsElapsed >= 1000) {
    const fps = Math.round((_fpsFrameCount * 1000) / fpsElapsed);
    const fpsEl = document.getElementById('fpsRead');
    if (fpsEl) fpsEl.textContent = String(fps);

    // Update Panel FPS
    if (UI.fpsPanelRead) UI.fpsPanelRead.textContent = String(fps);

    _fpsFrameCount = 0;
    _fpsLastTime = timestamp;
  }

  // Compute elapsed time since last frame
  if (_lastFrameTime === 0) _lastFrameTime = timestamp;
  const elapsed = Math.min(timestamp - _lastFrameTime, 200); // cap at 200ms to avoid spiral of death
  _lastFrameTime = timestamp;

  // Sub-step physics
  _physicsAccumulator += elapsed;
  const maxSteps = 10;  // safety cap
  let steps = 0;
  while (_physicsAccumulator >= STATE.dtMs && steps < maxSteps) {
    animation();
    _physicsAccumulator -= STATE.dtMs;
    steps++;
  }

  // Render & UI updates
  if (STATE.particles.length > 0) {
    syncEntitiesToParticles();

    // Force & energy readout (throttled by frame count)
    if (typeof frames === "number" && frames % 3 === 0) {
      updateDistanceAndEnergy();
      updateSelectedParticleDisplay();
      for (const pi of STATE.particles) {
        if (pi.type !== "electron") continue;
        UI.fx.textContent = (pi.acc.x * pi.mass).toExponential(2) + " Eₕ/a₀";
        UI.fy.textContent = (pi.acc.y * pi.mass).toExponential(2) + " Eₕ/a₀";
        UI.fz.textContent = (pi.acc.z * pi.mass).toExponential(2) + " Eₕ/a₀";
        break;
      }
    }

    // Record trajectory using InstancedMesh
    if (STATE.isRecordPath && typeof frames === "number" && frames % 5 === 0) {
      const S = SCENE_SCALE;
      for (const pi of STATE.particles) {
        if (pi.type !== "electron") continue;
        addTrajPoint(pi.pos.x * S, pi.pos.y * S, pi.pos.z * S);
      }
    }
  }

  requestAnimationFrame(mainLoop);
}
requestAnimationFrame(mainLoop);

// Initialize Panel dt display
if (UI.dtPanelRead) UI.dtPanelRead.textContent = STATE.dtMs.toFixed(0) + ' ms';

// ─── Public API ──────────────────────────────────────────────
const intoRealWorld = mainLoop;  // backward compat alias

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
window.mainLoop = mainLoop;
window.incrementEnergyLevel = incrementEnergyLevel;
window.decrementEnergyLevel = decrementEnergyLevel;
window.incrementAngularMomentum = incrementAngularMomentum;
window.decrementAngularMomentum = decrementAngularMomentum;
window.changeAnsatz = changeAnsatz;
window.tryTransition = tryTransition;
window.toggleRelaxation = toggleRelaxation;
window.toggleQHJ = toggleQHJ;
window.incrementM = incrementM;
window.decrementM = decrementM;w