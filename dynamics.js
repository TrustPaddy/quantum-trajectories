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
  selectedR:     document.getElementById("selectedR"),
  selectedPhi:   document.getElementById("selectedPhi"),
  selectedTheta: document.getElementById("selectedTheta"),
  // Spherical inputs
  Rs:     document.getElementById("Rs"),
  Phis:   document.getElementById("Phis"),
  Thetas: document.getElementById("Thetas"),
  // Coordinate toggle
  btnCartesian:  document.getElementById("btnCartesian"),
  btnSpherical:  document.getElementById("btnSpherical"),
  coordCartesian: document.getElementById("coordCartesian"),
  coordSpherical: document.getElementById("coordSpherical"),
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
  dynamicsMode: 'qhj',  // 'classic' or 'qhj'
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
    for (const r of STATE.bohrRings[i])
      try { SCENE.el.removeChild(r); } catch {}
  STATE.particles.length = STATE.spheres.length = STATE.bohrRings.length = 0;
}

function _removeParticleAt(idx) {
  const p = STATE.particles[idx];
  if (!p) return;
  const sphere = STATE.spheres[idx];
  if (sphere) try { SCENE.el.removeChild(sphere); } catch {}
  STATE.spheres.splice(idx, 1);
  if (p.type === 'proton') {
    let ringIdx = 0;
    for (let i = 0; i < idx; i++)
      if (STATE.particles[i].type === 'proton') ringIdx++;
    const rings = STATE.bohrRings[ringIdx];
    if (rings) rings.forEach(r => { try { SCENE.el.removeChild(r); } catch {} });
    STATE.bohrRings.splice(ringIdx, 1);
  }
  STATE.particles.splice(idx, 1);
  STATE.indexOfParticle = -1;
  updateSelectedParticleDisplay();
}

function addParticleEntity(particle, opts = {}) {
  const isProton = particle.type === "proton";
  const baseColor = opts.color ?? (isProton ? "#ff2222" : "#00aaff");
  const sphere = createEntity("a-sphere", {
    radius: opts.radius ?? (isProton ? 14 : 5),
    color: baseColor,
    emissive: baseColor,
    "emissive-intensity": 0.8,
    class: "clickable",
  });
  SCENE.el.appendChild(sphere);
  STATE.spheres.push(sphere);
}

function addBohrRings() {
  const r = SCENE_SCALE;
  const base = { "radius-inner": r * 0.99, "radius-outer": r * 1.01,
                  color: "white", opacity: 0.6, "ignore-ray": true };
  const rXY = createEntity("a-ring", base);
  const rXZ = createEntity("a-ring", { ...base, rotation: "-90 0 0" });
  const rYZ = createEntity("a-ring", { ...base, rotation: "0 90 0" });
  SCENE.el.appendChild(rXY);
  SCENE.el.appendChild(rXZ);
  SCENE.el.appendChild(rYZ);
  STATE.bohrRings.push([rXY, rXZ, rYZ]);
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
      const rings = STATE.bohrRings[ringIdx++];
      if (rings) rings.forEach(r => r.object3D.position.set(sx, sy, sz));
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
    if (UI.selectedR)     UI.selectedR.textContent     = "–";
    if (UI.selectedPhi)   UI.selectedPhi.textContent   = "–";
    if (UI.selectedTheta) UI.selectedTheta.textContent = "–";
    if (UI.vx) UI.vx.textContent = "0";
    if (UI.vy) UI.vy.textContent = "0";
    if (UI.vz) UI.vz.textContent = "0";
    if (UI.vMag) UI.vMag.textContent = "0";
    return;
  }

  const p = STATE.particles[idx];
  const pm = SI.a0 * 1e12;

  // Cartesian position
  if (UI.selectedX) UI.selectedX.textContent = (p.pos.x * pm).toFixed(3);
  if (UI.selectedY) UI.selectedY.textContent = (p.pos.y * pm).toFixed(3);
  if (UI.selectedZ) UI.selectedZ.textContent = (p.pos.z * pm).toFixed(3);

  // Spherical position
  const [r, theta, phi] = cartInPolar(p.pos.x, p.pos.y, p.pos.z);
  if (UI.selectedR)     UI.selectedR.textContent     = (r * pm).toFixed(3);
  if (UI.selectedPhi)   UI.selectedPhi.textContent   = phi.toFixed(4);
  if (UI.selectedTheta) UI.selectedTheta.textContent = theta.toFixed(4);

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

// Shift+Drag velocity arrow (sets initial velocity of selected particle)
const _velDrag = {
  active: false,
  wasRunning: false,
  plane:   new THREE.Plane(),
  lastHit: new THREE.Vector3(),
  arrow:   null,  // THREE.ArrowHelper, added directly to THREE.js scene
};
const VEL_SCALE = 0.5; // a.u. velocity per a.u. of arrow length

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

function _snapToAxisPlane(camDir, point) {
  const ax = Math.abs(camDir.x), ay = Math.abs(camDir.y), az = Math.abs(camDir.z);
  const normal = new THREE.Vector3();
  if (ax >= ay && ax >= az)      normal.set(1, 0, 0);
  else if (ay >= ax && ay >= az) normal.set(0, 1, 0);
  else                           normal.set(0, 0, 1);
  return new THREE.Plane().setFromNormalAndCoplanarPoint(normal, point);
}

function _onPointerDown(ev) {
  // Shift+click → velocity arrow mode for selected particle
  if (ev.shiftKey && STATE.indexOfParticle >= 0) {
    _updateMouse(ev);
    const p = STATE.particles[STATE.indexOfParticle];
    if (p) {
      // Remove any existing arrow before starting a new one
      if (_velDrag.arrow) {
        SCENE.el.object3D.remove(_velDrag.arrow);
        _velDrag.arrow = null;
      }
      const cam = _getCameraObj();
      const origin = p.pos.clone().multiplyScalar(SCENE_SCALE);
      const camDir = new THREE.Vector3();
      cam.getWorldDirection(camDir);
      _velDrag.plane.copy(_snapToAxisPlane(camDir, origin));
      _velDrag.lastHit.copy(origin);
      _velDrag.active = true;
      _velDrag.wasRunning = STATE.isAnimationRunning;
      if (STATE.isAnimationRunning) {
        STATE.isAnimationRunning = false;
        UI.myButton.textContent = "Start";
      }
      _velDrag.arrow = new THREE.ArrowHelper(
        new THREE.Vector3(1, 0, 0), origin, 0, 0x00ffff
      );
      _velDrag.arrow.visible = false;
      SCENE.el.object3D.add(_velDrag.arrow);
      ev.preventDefault();
      return;
    }
  }

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
  const [r0, theta0, phi0] = cartInPolar(particle.pos.x, particle.pos.y, particle.pos.z);
  if (UI.Rs)     UI.Rs.value     = (r0 * pm).toFixed(2);
  if (UI.Phis)   UI.Phis.value   = phi0.toFixed(4);
  if (UI.Thetas) UI.Thetas.value = theta0.toFixed(4);

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
  _drag.plane.copy(_snapToAxisPlane(camDir, worldPos));

  // Calculate grab offset (so particle doesn't snap to cursor)
  if (_drag.raycaster.ray.intersectPlane(_drag.plane, _drag.intersect)) {
    _drag.offset.copy(worldPos).sub(_drag.intersect);
  }


  // Visual highlight: make dragged particle emissive
  const sphere = STATE.spheres[idx];
  if (sphere) sphere.setAttribute('material', 'emissive', '#446');

  ev.preventDefault();
}

function _onPointerMove(ev) {
  _updateMouse(ev);

  // Velocity arrow drag (Shift+drag)
  if (_velDrag.active) {
    const cam = _getCameraObj();
    if (cam) {
      _drag.raycaster.setFromCamera(_drag.mouse, cam);
      const hit = new THREE.Vector3();
      if (_drag.raycaster.ray.intersectPlane(_velDrag.plane, hit)) {
        _velDrag.lastHit.copy(hit);
        const p = STATE.particles[STATE.indexOfParticle];
        if (p) {
          const origin = p.pos.clone().multiplyScalar(SCENE_SCALE);
          const delta = hit.clone().sub(origin);
          const len = delta.length();
          if (len > 1) {
            _velDrag.arrow.visible = true;
            _velDrag.arrow.setDirection(delta.clone().normalize());
            _velDrag.arrow.setLength(len, Math.min(len * 0.25, 30), Math.min(len * 0.12, 15));
          } else {
            _velDrag.arrow.visible = false;
          }
        }
      }
    }
    ev.preventDefault();
    return;
  }

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
      const [rD, thetaD, phiD] = cartInPolar(particle.pos.x, particle.pos.y, particle.pos.z);
      if (UI.Rs)     UI.Rs.value     = (rD * pm).toFixed(2);
      if (UI.Phis)   UI.Phis.value   = phiD.toFixed(4);
      if (UI.Thetas) UI.Thetas.value = thetaD.toFixed(4);
    }
  }

  ev.preventDefault();
}

function _onPointerUp(ev) {
  // Release velocity arrow
  if (_velDrag.active) {
    _velDrag.active = false;
    const p = STATE.particles[STATE.indexOfParticle];
    if (p && _velDrag.arrow?.visible) {
      const origin = p.pos.clone().multiplyScalar(SCENE_SCALE);
      const delta = _velDrag.lastHit.clone().sub(origin);
      const s = VEL_SCALE / SCENE_SCALE;
      p.vel.x += delta.x * s;
      p.vel.y += delta.y * s;
      p.vel.z += delta.z * s;
    }
    // Only remove arrow if simulation was running (it will resume).
    // If simulation was already stopped, keep arrow visible until Start is pressed.
    if (_velDrag.wasRunning && _velDrag.arrow) {
      SCENE.el.object3D.remove(_velDrag.arrow);
      _velDrag.arrow = null;
    }
    // Infer new quantum state from updated position + velocity
    if (p && STATE.particles[0]) {
      const qs = inferQuantumState(p.pos, STATE.particles[0].pos, p.vel);
      if (qs) {
        STATE.N = qs.n; STATE.L = qs.l; STATE.M = qs.m;
        updateQuantumUI();
        _notifyQuantumChange(qs.n, qs.l, qs.m);
      }
    }
    if (_velDrag.wasRunning) {
      STATE.isAnimationRunning = true;
      UI.myButton.textContent = "Stop";
    }
    return;
  }

  if (!_drag.active) return;

  // Remove highlight – restore to particle's base emissive colour
  const sphere = STATE.spheres[_drag.particleIdx];
  if (sphere) {
    const releasedParticle = STATE.particles[_drag.particleIdx];
    const baseEmissive = releasedParticle?.type === 'proton' ? '#ff2222' : '#00aaff';
    sphere.setAttribute('material', 'emissive', baseEmissive);
  }

  _drag.active = false;

  // Infer n from electron-proton distance (Bohr radii: r ≈ n²a₀)
  const draggedP = STATE.particles[_drag.particleIdx];
  if (draggedP?.type === 'electron' && STATE.particles[0]) {
    const r = draggedP.pos.distanceTo(STATE.particles[0].pos);
    const n = Math.max(1, Math.min(Math.round(Math.sqrt(r * STATE.Z)), 7));
    if (n !== STATE.N) {
      const l = Math.min(STATE.L, n - 1);
      const m = Math.max(-l, Math.min(STATE.M, l));
      STATE.N = n; STATE.L = l; STATE.M = m;
      updateQuantumUI();
      _notifyQuantumChange(n, l, m);
    }
  }

  // Remove dragging cursor class
  const sceneEl = SCENE.el;
  if (sceneEl) sceneEl.classList.remove('dragging-particle');

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

// Delete / Backspace → remove selected particle
document.addEventListener('keydown', ev => {
  if (ev.key !== 'Delete' && ev.key !== 'Backspace') return;
  // Don't fire when typing in an input/textarea
  const tag = ev.target?.tagName;
  if (tag === 'INPUT' || tag === 'TEXTAREA') return;
  const idx = STATE.indexOfParticle;
  if (idx < 0 || !STATE.particles[idx]) return;
  _removeParticleAt(idx);
});

// ─── Atom builders (positions in a₀) ─────────────────────────
function pushProton(posAU, radius = 14) { 
  STATE.particles.push(new Particle(MP_AU, +1, posAU.clone(), "proton"));  
  addParticleEntity(STATE.particles.at(-1), { radius, color: "#ff2222" });
  addBohrRings();
}
function pushElectron(posAU, radius = 5) { 
  STATE.particles.push(new Particle(1,-1, posAU.clone(), "electron")); 
  addParticleEntity(STATE.particles.at(-1), { radius, color: "#00aaff" });
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
  // Remove velocity preview arrow when simulation starts
  if (STATE.isAnimationRunning && _velDrag.arrow) {
    SCENE.el.object3D.remove(_velDrag.arrow);
    _velDrag.arrow = null;
  }
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
  if (!c) return;
  STATE.radiation = c.radiation;

  // Protonen-Schwerpunkt in Szeneneinheiten
  const protons = STATE.particles.filter(p => p.type === 'proton');
  const centroid = new THREE.Vector3();
  if (protons.length) {
    protons.forEach(p => centroid.add(p.pos));
    centroid.divideScalar(protons.length).multiplyScalar(SCENE_SCALE);
  }

  const oc = SCENE.camera?.components?.['orbit-controls'];
  if (oc) {
    oc.target.copy(centroid);
    oc.distance = Math.max(oc.data.minDistance, Math.min(oc.data.maxDistance, c.z));
    oc.theta = 0;
    oc.phi   = Math.PI / 2;
  } else {
    SCENE.camera.setAttribute("position", { x: centroid.x, y: centroid.y, z: centroid.z + c.z });
  }
}
function updateQuantumUI() { 
  UI.nLvl.textContent = String(STATE.N); 
  UI.lLvl.textContent = String(STATE.L); 
  if (UI.mLvl) UI.mLvl.textContent = String(STATE.M); 
}

function _notifyQuantumChange(n, l, m) {
  if (typeof window.showToast === 'function') {
    window.showToast(`Quantum state → (${n}, ${l}, ${m})`, 'info');
  }
}

// Infer closest quantum state (n,l,m) from classical position and velocity.
// pos/posNuc: THREE.Vector3 in a₀, vel: THREE.Vector3 in a.u.
// Returns {n,l,m} or null if unbound (E≥0).
function inferQuantumState(pos, posNuc, vel) {
  const rVec = pos.clone().sub(posNuc);
  const r = rVec.length();
  if (r < 1e-6) return null;

  const v2 = vel.lengthSq();
  const E = 0.5 * v2 - STATE.Z / r;   // Hartree: E_n = -Z²/(2n²)
  if (E >= 0) return null;             // ionised

  const n_raw = STATE.Z / Math.sqrt(-2 * E);
  const n = Math.max(1, Math.min(Math.round(n_raw), 7));

  // Angular momentum L = r × v (in ℏ, since ℏ=1 a.u.)
  const Lvec = rVec.clone().cross(vel);
  const L = Lvec.length();
  // l(l+1) = L²  →  l = (−1 + √(1+4L²)) / 2
  const l_raw = (-1 + Math.sqrt(1 + 4 * L * L)) / 2;
  const l = Math.max(0, Math.min(Math.round(l_raw), n - 1));

  // m from z-component of L (quantisation axis = z).
  // If the orbit is mostly around the z-axis (|Lz|/|L| > 0.5), maximize |m| = l
  // so that a horizontal tangential impulse cleanly gives m = ±l.
  const Lz = Lvec.z;
  let m;
  if (l === 0) {
    m = 0;
  } else if (L > 0.1 && Math.abs(Lz) / L > 0.5) {
    m = Math.sign(Lz) * l;
  } else {
    m = Math.max(-l, Math.min(Math.round(Lz), l));
  }

  return { n, l, m };
}

function tryTransition(nN, nL, nM) {
  if (nN<1||nL<0||nL>=nN||Math.abs(nM)>nL) {
    if (typeof window.showToast === 'function') {
      window.showToast(`Invalid quantum numbers: n≥1, 0≤ℓ≤n−1, |m|≤ℓ`, 'error');
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
      // Classic mode: place electron at (ρ_eq, θ_eq) from Yang potential
      const thetaEq = window.yangEquilibriumTheta
        ? window.yangEquilibriumTheta(nL, nM)
        : Math.PI / 2;

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
        e.vel.set(
          -Math.sign(nM) * vCirc * Math.sin(phi),
           Math.sign(nM) * vCirc * Math.cos(phi),
          0
        );
      } else if (nL > 0) {
        // m=0: add small radial kick to show oscillation
        const radHat = new THREE.Vector3(
          sinTh * Math.cos(phi),
          sinTh * Math.sin(phi),
          cosTh
        );
        const kickFraction = 0.08 / nN;
        e.pos.addScaledVector(radHat, kickFraction * newR);
      }
    }
    e.acc.set(0, 0, 0);
  }
  applyLevelCamera(); updateQuantumUI(); return true;
}

function incrementEnergyLevel() {
  const newN = STATE.N + 1;
  tryTransition(newN, newN - 1, newN - 1);
}
function decrementEnergyLevel() {
  const newN = STATE.N - 1;
  if (newN < 1) return;
  const newL = Math.min(STATE.L, newN - 1);
  const newM = Math.max(-newL, Math.min(STATE.M, newL));
  tryTransition(newN, newL, newM);
}
function incrementAngularMomentum() {
  tryTransition(STATE.N, STATE.L+1, STATE.M);
}
function decrementAngularMomentum() {
  const newL = STATE.L - 1;
  const newM = Math.max(-newL, Math.min(STATE.M, newL));
  tryTransition(STATE.N, newL, newM);
}
function incrementM() {
  tryTransition(STATE.N, STATE.L, STATE.M + 1);
}
function decrementM() {
  tryTransition(STATE.N, STATE.L, STATE.M - 1);
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
        const rings = STATE.bohrRings.pop();
        if (rings) rings.forEach(r => SCENE.el.removeChild(r));
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

// Most-probable radius for state (n,l,m) - used for transition rescaling
function bohrRadius(n, l, Z) {
  // Prefer yangEquilibriumRadius (handles all states correctly)
  if (window.yangEquilibriumRadius) {
    return window.yangEquilibriumRadius(n, l, STATE.M || 0, Z);
  }
  // Fallback: Newton's method on the classic radial potential
  let r = n * n / Z;
  for (let iter = 0; iter < 40; iter++) {
    if (r < 0.1) r = 0.1;
    const h = Math.max(r * 1e-4, 1e-6);
    // Use classicPotentialRy for radial force at θ=π/2
    const Vp = classicPotentialRy(r + h, Math.PI/2, n, l, STATE.M || 0, Z);
    const Vm = classicPotentialRy(r - h, Math.PI/2, n, l, STATE.M || 0, Z);
    const F  = -(Vp - Vm) / (2 * h);  // f̄_ρ at θ=π/2
    const Vpp = classicPotentialRy(r + 2*h, Math.PI/2, n, l, STATE.M || 0, Z);
    const Vmm = classicPotentialRy(r - 2*h, Math.PI/2, n, l, STATE.M || 0, Z);
    const Fp  = -(Vpp - Vp) / (2 * h);
    const Fm  = -(Vm - Vmm) / (2 * h);
    const dF  = (Fp - Fm) / (2 * h);
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

// ═════════════════════════════════════════════════════════════
//  CLASSIC MODE: Explicit Yang quantum-potential per state
//  ─────────────────────────────────────────────────────────
//  V̄(ρ,θ)/Ry = -2Z/ρ + A(θ)/ρ² + B(ρ)
//
//  Forces (Ry):  f̄_ρ = -2Z/ρ² + 2A(θ)/ρ³ - B'(ρ)
//                f̄_θ = -A'(θ)/ρ²
//
//  Conversion to Hartree/a₀:  F_r = 0.5·f̄_ρ
//                              F_θ = 0.5·f̄_θ/ρ
//
//  Source: Yang_quantum-potential(4).odt & equations.tex
// ═════════════════════════════════════════════════════════════

// ─── Explicit total potential V̄/Ry for tabulated states ─────
function classicPotentialRy(rho, theta, n, l, m, Z) {
  rho   = Math.max(rho, 1e-6);
  theta = Math.max(1e-6, Math.min(Math.PI - 1e-6, theta));

  const rho2 = rho * rho;
  const c = Math.cos(theta), s = Math.sin(theta);
  const c2 = c * c, s2 = s * s;
  const EPS = 1e-15;
  const cot2 = c2 / (s2 || EPS);             // cot²θ

  // Coulomb (common to all states)
  const vCoul = -2 * Z / rho;

  // Encode state as single key for switch (potential is symmetric in m)
  const key = n * 100 + l * 10 + Math.abs(m);

  let A, B;
  switch (key) {

    case 100:  // (1,0,0)
      A = (4 + cot2) / 4;
      B = 0;
      break;

    case 200:  // (2,0,0)
      A = (4 + cot2) / 4;
      { const d = Z * rho - 2; B = (Z * Z) / (d * d || EPS); }
      break;

    case 210:  // (2,1,0)
      A = (8 + cot2 + 4 / (c2 || EPS)) / 4;   // 4sec²θ = 4/cos²θ
      B = 0;
      break;

    case 211:  // (2,1,1)
      A = (8 + cot2 + 4 / (s2 || EPS)) / 4;   // 4csc²θ = 4/sin²θ
      B = 0;
      break;

    case 300:  // (3,0,0)
      A = (4 + cot2) / 4;
      { const D = 2*Z*Z*rho2 - 18*Z*rho + 27;
        const t = 4*Z*rho - 18;
        B = -4*Z*Z / (D || EPS) + Z*Z*t*t / ((D*D) || EPS);
      }
      break;

    case 310:  // (3,1,0)
      { const tan2 = s2 / (c2 || EPS);
        A = (12 + cot2 + 4 * tan2) / 4;
      }
      { const d = Z * rho - 6; B = (Z * Z) / (d * d || EPS); }
      break;

    case 311:  // (3,1,1)
      A = (12 + 5 * cot2) / 4;
      { const d = Z * rho - 6; B = (Z * Z) / (d * d || EPS); }
      break;

    case 320:  // (3,2,0)
      { const d320 = 3 * c2 - 1;
        A = (12 + cot2) / 4 - 6 * (s2 - 2) / (d320 * d320 || EPS);
      }
      B = 0;
      break;

    case 321:  // (3,2,1)
      A = (12 + (4 - 5 * s2) / ((s2 * c2) || EPS)) / 4;
      B = 0;
      break;

    case 322:  // (3,2,2)
      A = (12 + 9 / (s2 || EPS)) / 4;          // 9csc²θ = 9/sin²θ
      B = 0;
      break;

    default:
      // Not tabulated — fall back to general yang potential
      if (window.yangPotentialDimless)
        return window.yangPotentialDimless(rho, theta, n, l, m, Z);
      // Last resort: Coulomb + universal quantum term only
      A = (4 + cot2) / 4;
      B = 0;
      break;
  }

  let V = vCoul + A / rho2 + B;
  if (!Number.isFinite(V)) V = 2000;
  return Math.max(-2000, Math.min(2000, V));
}

// ─── Explicit forces f̄_ρ, f̄_θ in Ry units ─────────────────
function classicForceRy(rho, theta, n, l, m, Z) {
  rho   = Math.max(rho, 1e-6);
  theta = Math.max(1e-6, Math.min(Math.PI - 1e-6, theta));

  const rho2 = rho * rho, rho3 = rho2 * rho;
  const c = Math.cos(theta), s = Math.sin(theta);
  const c2 = c * c, s2 = s * s;
  const EPS = 1e-15;
  const cot2 = c2 / (s2 || EPS);
  // cotθ·csc²θ = cosθ / sin³θ
  const cotCsc2 = c / ((s2 * s) || EPS);

  const key = n * 100 + l * 10 + Math.abs(m);

  let A, Ap, Bp;  // A(θ), A'(θ), B'(ρ)

  switch (key) {

    case 100:  // (1,0,0)
      A  = (4 + cot2) / 4;
      Ap = -cotCsc2 / 2;
      Bp = 0;
      break;

    case 200:  // (2,0,0)
      A  = (4 + cot2) / 4;
      Ap = -cotCsc2 / 2;
      { const d = Z * rho - 2; Bp = -2*Z*Z*Z / ((d*d*d) || EPS); }
      break;

    case 210:  // (2,1,0)
      { const sec2 = 1 / (c2 || EPS);
        const tan_ = s / (c || EPS);
        A  = (8 + cot2 + 4 * sec2) / 4;
        Ap = (-cotCsc2 + 4 * sec2 * tan_) / 2;
      }
      Bp = 0;
      break;

    case 211:  // (2,1,1)
      { const csc2 = 1 / (s2 || EPS);
        A  = (8 + cot2 + 4 * csc2) / 4;
        Ap = -5 * cotCsc2 / 2;
        // d/dθ(cot²θ + 4csc²θ) = -2cotθcsc²θ - 8csc²θcotθ = -10cotθcsc²θ → /4 = -5cotθcsc²θ/2
      }
      Bp = 0;
      break;

    case 300:  // (3,0,0)
      A  = (4 + cot2) / 4;
      Ap = -cotCsc2 / 2;
      { const D = 2*Z*Z*rho2 - 18*Z*rho + 27;
        const D3 = D * D * D;
        const t = 4*Z*rho - 18;
        // B = -4Z²/D + Z²t²/D²
        // B' = -4Z³·t·(2Z²ρ²-18Zρ+81)/D³
        const inner = 2*Z*Z*rho2 - 18*Z*rho + 81;
        Bp = -4*Z*Z*Z * t * inner / (D3 || EPS);
      }
      break;

    case 310:  // (3,1,0)
      { const tan2 = s2 / (c2 || EPS);
        const tan_ = s / (c || EPS);
        const sec2 = 1 / (c2 || EPS);
        A  = (12 + cot2 + 4 * tan2) / 4;
        Ap = (-cotCsc2 + 4 * tan_ * sec2) / 2;
      }
      { const d = Z * rho - 6; Bp = -2*Z*Z*Z / ((d*d*d) || EPS); }
      break;

    case 311:  // (3,1,1)
      A  = (12 + 5 * cot2) / 4;
      Ap = -5 * cotCsc2 / 2;
      { const d = Z * rho - 6; Bp = -2*Z*Z*Z / ((d*d*d) || EPS); }
      break;

    case 320:  // (3,2,0) — use numerical derivatives
    case 321:  // (3,2,1) — use numerical derivatives
      { const h = Math.max(rho * 1e-4, 1e-6);
        const hth = 1e-5;
        const Vrp = classicPotentialRy(rho + h, theta, n, l, m, Z);
        const Vrm = classicPotentialRy(rho - h, theta, n, l, m, Z);
        const Vtp = classicPotentialRy(rho, theta + hth, n, l, m, Z);
        const Vtm = classicPotentialRy(rho, theta - hth, n, l, m, Z);
        return {
          fr:  -(Vrp - Vrm) / (2 * h),
          fth: -(Vtp - Vtm) / (2 * hth)
        };
      }

    case 322:  // (3,2,2)
      { const csc2 = 1 / (s2 || EPS);
        A  = (12 + 9 * csc2) / 4;
        Ap = -9 * cotCsc2 / 2;
      }
      Bp = 0;
      break;

    default:
      // Fall back to general yang force if available
      if (window.yangForceCartesian) return null;  // signal to use yangForceCartesian
      // Last resort: Coulomb + universal quantum term
      A  = (4 + cot2) / 4;
      Ap = -cotCsc2 / 2;
      Bp = 0;
      break;
  }

  // f̄_ρ = -2Z/ρ² + 2A(θ)/ρ³ - B'(ρ)
  const fr  = -2 * Z / rho2 + 2 * A / rho3 - Bp;
  // f̄_θ = -A'(θ)/ρ²
  const fth = -Ap / rho2;

  return { fr, fth };
}

// ─── Classic force in Cartesian coords [Hartree/a₀] ─────────
const CLASSIC_FORCE_CAP_BASE = 200;

function classicForceCartesian(x, y, z, n, l, m, Z) {
  const forceCap = CLASSIC_FORCE_CAP_BASE * Math.max(1, n * n);

  const r2  = x * x + y * y + z * z;
  const r   = Math.sqrt(r2);
  if (r < 1e-10) return { fx: 0, fy: 0, fz: 0 };

  const rho = r;
  const cosTheta = Math.max(-1, Math.min(1, z / r));
  const theta    = Math.acos(cosTheta);
  const rxy      = Math.sqrt(x * x + y * y);

  let fx, fy, fz;

  // Near z-axis: numerical gradient fallback (avoid 1/sinθ singularity)
  const USE_NUMERICAL = rxy < 1e-8 * r;

  if (!USE_NUMERICAL) {
    const forces = classicForceRy(rho, theta, n, l, m, Z);

    // If null: fallback to yangForceCartesian for unsupported states
    if (forces === null && window.yangForceCartesian) {
      return window.yangForceCartesian(x, y, z, n, l, m, Z);
    }
    if (forces === null) return { fx: 0, fy: 0, fz: 0 };

    const { fr: fbar_r, fth: fbar_th } = forces;

    // Convert Ry → Hartree:  F_r = 0.5·f̄_ρ,  F_θ_eff = 0.5·f̄_θ/r
    const Fr      = 0.5 * fbar_r;
    const Fth_eff = 0.5 * fbar_th / r;

    // Spherical → Cartesian:
    //   r̂  = (x/r, y/r, z/r)
    //   ∂θ/∂(x,y,z) = (xz, yz, −rxy²) / (r²·rxy)
    const geom = 1.0 / (r2 * rxy);
    fx = Fr * (x / r) + Fth_eff * (x * z) * geom;
    fy = Fr * (y / r) + Fth_eff * (y * z) * geom;
    fz = Fr * (z / r) - Fth_eff * (rxy * rxy) * geom;

    if (!Number.isFinite(fx)) fx = 0;
    if (!Number.isFinite(fy)) fy = 0;
    if (!Number.isFinite(fz)) fz = 0;

  } else {
    // Numerical gradient fallback (near z-axis)
    let h = Math.max(r * 1e-4, 1e-6);
    const potXYZ = (xx, yy, zz) => {
      const rr = Math.sqrt(xx*xx + yy*yy + zz*zz);
      if (rr < 1e-6) return 0;
      const th = Math.acos(Math.max(-1, Math.min(1, zz / rr)));
      return 0.5 * classicPotentialRy(rr, th, n, l, m, Z);
    };

    let Vxp = potXYZ(x+h, y, z), Vxm = potXYZ(x-h, y, z);
    let Vyp = potXYZ(x, y+h, z), Vym = potXYZ(x, y-h, z);
    let Vzp = potXYZ(x, y, z+h), Vzm = potXYZ(x, y, z-h);

    const maxDV = Math.max(Math.abs(Vxp-Vxm), Math.abs(Vyp-Vym), Math.abs(Vzp-Vzm));
    if (maxDV > 10) {
      h *= 0.01;
      Vxp = potXYZ(x+h, y, z); Vxm = potXYZ(x-h, y, z);
      Vyp = potXYZ(x, y+h, z); Vym = potXYZ(x, y-h, z);
      Vzp = potXYZ(x, y, z+h); Vzm = potXYZ(x, y, z-h);
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

// ─── Spherical position input (r in pm, φ/θ in rad → a.u.) ──
["Rs","Phis","Thetas"].forEach((key, dimIdx) => {
  UI[key].addEventListener("keypress", ev => {
    if (ev.key !== "Enter") return;
    const val = parseFloat(UI[key].value);
    if (!Number.isFinite(val)) return;
    const p = STATE.particles[STATE.indexOfParticle];
    if (!p) return;
    const [curR, curTheta, curPhi] = cartInPolar(p.pos.x, p.pos.y, p.pos.z);
    let r = curR, theta = curTheta, phi = curPhi;
    if (dimIdx === 0)      r     = (val * 1e-12) / SI.a0; // r: pm → a.u.
    else if (dimIdx === 1) phi   = val;                    // φ in rad
    else                   theta = val;                    // θ in rad
    const sinT = Math.sin(theta);
    p.pos.set(r * sinT * Math.cos(phi), r * sinT * Math.sin(phi), r * Math.cos(theta));
    p.acc.set(0, 0, 0);
  });
});

// ─── Coordinate toggle buttons ───────────────────────────────
if (UI.btnCartesian) {
  UI.btnCartesian.addEventListener("click", () => {
    UI.coordCartesian.style.display = "";
    UI.coordSpherical.style.display = "none";
    UI.btnCartesian.style.borderColor = "rgba(125,211,252,.6)";
    UI.btnSpherical.style.borderColor = "";
  });
  UI.btnSpherical.addEventListener("click", () => {
    UI.coordCartesian.style.display = "none";
    UI.coordSpherical.style.display = "";
    UI.btnSpherical.style.borderColor = "rgba(125,211,252,.6)";
    UI.btnCartesian.style.borderColor = "";
  });
}

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
    // ── Classic: explicit Yang quantum-potential per state ──
    for (let i = 0; i < pArr.length; i++) {
      const pi = pArr[i];
      pi.acc.set(0, 0, 0);
      if (STATE.isNucleusLocked && pi.type === "proton") continue;

      for (let j = 0; j < pArr.length; j++) {
        if (i === j) continue;
        const pj = pArr[j];

        const rx = pi.pos.x - pj.pos.x;
        const ry = pi.pos.y - pj.pos.y;
        const rz = pi.pos.z - pj.pos.z;

        if (pi.type === "electron" && pj.type === "proton") {
          // Full Yang quantum potential (explicit closed-form per state)
          const Z = Math.abs(pj.charge);
          const f = classicForceCartesian(rx, ry, rz, STATE.N, STATE.L, STATE.M, Z);
          pi.acc.x += f.fx;
          pi.acc.y += f.fy;
          pi.acc.z += f.fz;
        } else if (pi.type === "proton" && pj.type === "electron") {
          // Newton's 3rd law: proton feels opposite force (scaled by mass)
          const Z = Math.abs(pi.charge);
          const f = classicForceCartesian(-rx, -ry, -rz, STATE.N, STATE.L, STATE.M, Z);
          pi.acc.x -= f.fx / pi.mass;
          pi.acc.y -= f.fy / pi.mass;
          pi.acc.z -= f.fz / pi.mass;
        } else {
          // e-e or p-p: classical Coulomb only
          const r2 = rx*rx + ry*ry + rz*rz + SOFTENING*SOFTENING;
          const r  = Math.sqrt(r2);
          const Fel = (-pi.charge * pj.charge) / r2;
          _v.dir.set(rx, ry, rz).divideScalar(r);
          pi.acc.addScaledVector(_v.dir, Fel / pi.mass);
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