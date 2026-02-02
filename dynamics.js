"use strict";

// -------------------------------------------------------------
// Constants & utilities
// -------------------------------------------------------------
/**
 * Physical constants used in the simulation.
 * @constant
 * @type {Object}
 * @property {number} e0 - Electric constant (8.8541878128e-12 C^2/N m^2)
 * @property {number} a0 - Bohr radius (5.29177210903e-11 m)
 * @property {number} c - Speed of light (2.99792e8 m/s)
 */
const PHYS = Object.freeze({
  e0: 8.8541878128e-12,
  a0: 5.29177210903e-11,
  c: 2.99792e8,
});

/**
 * Coulomb constant (4 * pi * e0)
 * @constant
 * @type {number}
 */
const K = 4 * Math.PI * PHYS.e0;

/**
 * Radiation constant (2 / (3 * c^3 * c^3 * c^3))
 * @constant
 * @type {number}
 */
const RAD_CONST = 2 / (3 * PHYS.c * PHYS.c * PHYS.c);

/**
 * Enum for atomic models.
 * @constant
 * @type {Object}
 * @property {string} HAtom - H atom
 * @property {string} H2Cation - H2+ cation
 * @property {string} HAnion - H- anion
 * @property {string} HECation - He+ cation
 * @property {string} HEAtom - He atom
 */
const AtomModel = Object.freeze({
  HAtom: "H Atom",
  H2Cation: "H2 Cation",
  HAnion: "H Anion",
  HECation: "HE Cation",
  HEAtom: "HE Atom",
});

/**
 * Clamps a number between a minimum and maximum value.
 * @param {number} n The number to clamp.
 * @param {number} min The minimum value.
 * @param {number} max The maximum value.
 * @returns {number} The clamped number.
 */
function clamp(n, min, max) {
  // Ensures that n is within the range [min, max].
  return Math.max(min, Math.min(max, n));
}

/**
 * Converts a 3D cartesian coordinate (x, y, z) into a 3D polar coordinate (r, theta, phi).
 * @param {number} x The x coordinate of the cartesian point.
 * @param {number} y The y coordinate of the cartesian point.
 * @param {number} z The z coordinate of the cartesian point.
 * @returns {Array<number>} An array containing the polar coordinates [r, theta, phi].
 */
function cartInPolar(x, y, z) {
  // Calculate the radial distance (r) from the origin.
  const r = Math.sqrt(x * x + y * y + z * z) || 0;
  
  // Calculate the polar angle (theta) from the positive z-axis.
  const theta = r === 0 ? 0 : Math.acos(z / r);
  
  // Calculate the azimuthal angle (phi) from the positive x-axis.
  const phi = Math.atan2(y, x);
  
  return [r, theta, phi];
}

/**
 * Creates a new HTML element with the given tag and attributes.
 * @param {string} tag The tag name of the HTML element to create.
 * @param {Object<string, string>} [attrs] An object containing key-value pairs of attributes to set on the element.
 * @returns {HTMLElement} The created HTML element.
 */
function createEntity(tag, attrs) {
  const el = document.createElement(tag);
  if (attrs) {
    // Iterate over the given attributes and set them on the element.
    for (const [k, v] of Object.entries(attrs)) {
      el.setAttribute(k, v);
    }
  }
  return el;
}


const SAFE_DISTANCE = 1e-11; // ~1 Szene-Einheit, weil später mit 1e12 skaliert wird
const MAX_TRIES = 50;
/**
 * Try to find a non-colliding position for a new particle.
 * @param {string} type - either "proton" or "electron"
 * @param {() => THREE.Vector3} generator - function that returns a candidate position
 * @param {number} safetyDistance - minimal allowed distance to any existing particle (in simulation units)
 * @param {number} maxTries - how many random attempts before giving up
 * @returns {THREE.Vector3|null} a valid position or null if none found
 */
function findNonCollidingPosition(type, generator, safetyDistance = SAFE_DISTANCE, maxTries = MAX_TRIES,) {
  const radius1 = type === "electron" ? 5e-12 : 14e-12;
  for (let attempt = 0; attempt < maxTries; attempt++) {
    const candidate = generator();
    let ok = true;

    for (const p of STATE.particles) {
      const radius2 = p.type === "electron" ? 5e-12 : 14e-12;
      if (p.cartPos.distanceTo(candidate) < safetyDistance + radius1 + radius2) {

        console.log(p.cartPos.distanceTo(candidate))
        console.log(p.type)

        ok = false;
        break;
      }
    }

    if (ok) return candidate;
  }

  console.warn("Could not find non-colliding position for new particle.");
  return null;
}

// -------------------------------------------------------------
// DOM & scene cache
// -------------------------------------------------------------
const UI = {
  Distance: document.getElementById("Distance"),
  Energy: document.getElementById("Energy"),
  posx1: document.getElementById("posx1"),
  y1: document.getElementById("y1"),
  z1: document.getElementById("z1"),
  x2: document.getElementById("x2"),
  y2: document.getElementById("y2"),
  z2: document.getElementById("z2"),
  nLvl: document.getElementById("nLvl"),
  lLvl: document.getElementById("lLvl"),
  DispCountElectron: document.getElementById("DispCountElectron"),
  DispCountProton: document.getElementById("DispCountProton"),
  Xs: document.getElementById("Xs"),
  Ys: document.getElementById("Ys"),
  Zs: document.getElementById("Zs"),
  led: document.getElementById("led"),
  myButton: document.getElementById("myButton"),
  pathButton: document.getElementById("PathButton"),
  fx: document.getElementById("fx"),
  fy: document.getElementById("fy"),
  fz: document.getElementById("fz"),
};

const SCENE = {
  sceneEl: document.querySelector("a-scene"),
  cameraEl: document.querySelector("#camera"),
};

// -------------------------------------------------------------
// Types & state
// -------------------------------------------------------------
class Particle {
  constructor(mass, charge, spin, velocity, cartPos, type) {
    this.mass = mass;
    this.charge = charge;
    this.spin = spin;
    this.velocity = velocity; // THREE.Vector3
    this.cartPos = cartPos;   // THREE.Vector3
    this.polarPos = cartInPolar(cartPos.x, cartPos.y, cartPos.z);
    this.type = type;         // "proton" | "electron"
  }
}

const STATE = {
  atomLabel: AtomModel.HAtom, // pretty label used in UI/energy switch
  particles: /** @type {Particle[]} */([]),
  spheres: /** @type {HTMLElement[]} */([]),
  bohrRings: /** @type {HTMLElement[]} */([]),
  isAnimationRunning: false,
  isNucleusLocked: false,
  isRecordPath: false,
  isAngularMomentum: false,
  isLarmor: false,
  isThetaPhi: false,
  dtMs: 1000 / 30, // physics step in ms (changed by changeInterval)
  radiation: 1,
  N: 1,
  L: 0,
  Z: 1, // atomic number for energy/forces where used
  indexOfParticle: -1,
  timePassedPs: 0,
  countTime: 0,
  angleCache: 0,
};

// -------------------------------------------------------------
// Scene helpers
// -------------------------------------------------------------
/**
 * Clears the scene of all particles and Bohr rings.
 * This function is typically called when the user wants to reset the simulation.
 */
function clearScene() {
  // remove all particle entities
  for (let i = STATE.spheres.length - 1; i >= 0; i--) {
    try { SCENE.sceneEl.removeChild(STATE.spheres[i]); } catch {}
  }
  // remove all Bohr rings
  for (let i = STATE.bohrRings.length - 1; i >= 0; i--) {
    try { SCENE.sceneEl.removeChild(STATE.bohrRings[i]); } catch {}
  }
  // reset particle state
  STATE.particles.length = 0;
  STATE.spheres.length = 0;
  STATE.bohrRings.length = 0;
}

/**
 * Adds a particle entity to the scene.
 * @param {Particle} particle - The particle to render
 * @param {Object} [opts] - Optional parameters
 * @param {number} [opts.radius] - The radius of the particle
 * @param {string} [opts.color] - The color of the particle
 */
function addParticleEntity(particle, opts = {}) {
  // Create a new sphere entity
  const sphere = createEntity("a-sphere", {
    radius: opts.radius ?? (particle.type === "proton" ? 14 : 5),
    color: opts.color ?? (particle.type === "proton" ? "red" : "blue"),
    class: "clickable", // wichtig für den Raycaster
  });

  // Add the sphere to the scene
  SCENE.sceneEl.appendChild(sphere);

  // Store the sphere in the state
  STATE.spheres.push(sphere);

  // Add an event listener to the sphere to handle selection
  const idx = STATE.spheres.length - 1;
  sphere.addEventListener("click", () => {
    // dieses Partikel ist jetzt ausgewählt
    STATE.indexOfParticle = idx;

    // aktuelle Position aus der Szene (in pm) lesen
    const pos = sphere.object3D.position;

    UI.Xs.value = String(pos.x.toFixed(2));
    UI.Ys.value = String(pos.y.toFixed(2));
    UI.Zs.value = String(pos.z.toFixed(2));
  });
}

/**
 * Adds a new Bohr ring to the scene.
 * A Bohr ring is a grey ring that represents the orbit of an electron.
 */
function addBohrRing() {
  const ring = createEntity("a-ring", {
    /**
     * The inner radius of the Bohr ring.
     * This is set to the Bohr radius of a hydrogen atom (9.99e11 pm).
     */
    "radius-inner": PHYS.a0 * 9.99e11,

    /**
     * The outer radius of the Bohr ring.
     * This is set to the Bohr radius of a hydrogen atom (1.01e12 pm).
     */
    "radius-outer": PHYS.a0 * 1.01e12,
    color: "grey",
    opacity: 0.2,
    "ignore-ray": true,
  });

  // Add the ring to the scene
  SCENE.sceneEl.appendChild(ring);

  // Store the ring in the state
  STATE.bohrRings.push(ring);
}

/**
 * Synchronizes the position of particles in the state with the position of the sphere entities in the scene.
 * This function is called every frame to ensure that the particles are correctly positioned in the scene.
 */
function syncEntitiesToParticles() {
  let ringIndex = 0;
  for (let i = 0; i < STATE.particles.length; i++) {
    const p = STATE.particles[i].cartPos.clone().multiplyScalar(1e12);
    const sphere = STATE.spheres[i];
    if (sphere) sphere.object3D.position.set(p.x, p.y, p.z);

    // If particle is a proton, update the position of the associated Bohr ring
    if (STATE.particles[i].type === "proton") {
      const ring = STATE.bohrRings[ringIndex++];
      if (ring) ring.object3D.position.set(p.x, p.y, p.z);
    }

    // UI readout for first/second electron
    if (STATE.particles[i].type === "electron") {
      if (i === 1) {
        // Update the UI fields for the first electron
        UI.posx1.textContent = p.x.toFixed(3);
        UI.y1.textContent = p.y.toFixed(3);
        UI.z1.textContent = p.z.toFixed(3);
      } else if (i === 2) {
        // Update the UI fields for the second electron
        UI.x2.textContent = p.x.toFixed(3);
        UI.y2.textContent = p.y.toFixed(3);
        UI.z2.textContent = p.z.toFixed(3);
      }
    }
  }

  // If first/second electron not present, blank the fields
  if (!STATE.particles[1] || STATE.particles[1].type !== "electron") {
    // Blank the UI fields for the first electron if it doesn't exist
    UI.posx1.textContent = UI.y1.textContent = UI.z1.textContent = "NaN";
  }
  if (!STATE.particles[2] || STATE.particles[2].type !== "electron") {
    // Blank the UI fields for the second electron if it doesn't exist
    UI.x2.textContent = UI.y2.textContent = UI.z2.textContent = "NaN";
  }
}

// -------------------------------------------------------------
// Builders for atoms
// -------------------------------------------------------------
/**
 * Pushes a proton to the given position in the scene.
 * @param {THREE.Vector3} pos - The position of the proton in the scene.
 */
function pushProton(pos, radius = 14) {
  // Create a new particle with the given position and properties of a proton
  const particle = new Particle(
    1.6726219e-27, // mass of a proton in kg
    +1.602176634e-19, // charge of a proton in C
    +0.5, // spin of a proton 
    new THREE.Vector3(0, 0, 0), // initial velocity of the proton
    pos.clone(), // position of the proton in the scene
    "proton" // type of the particle
  );

  // Add the proton to the state
  STATE.particles.push(particle);

  // Add a sphere entity to the scene to represent the proton
  addParticleEntity(STATE.particles[STATE.particles.length - 1], {
    radius, // radius of the sphere entity
    color: "red" // color of the sphere entity
  });

  // Add a Bohr ring to the scene
  addBohrRing();
}

/**
 * Pushes an electron to the given position in the scene.
 * @param {THREE.Vector3} pos - The position of the electron in the scene.
 * @param {number} [radius=5] - The radius of the sphere entity that represents the electron.
 */
function pushElectron(pos, radius = 5) {
  // Create a new particle with the given position and properties of an electron
  STATE.particles.push(new Particle(
    9.10938356e-31, // mass of an electron in kg
    -1.602176634e-19, // charge of an electron in C
    +0.5, // spin of an electron
    new THREE.Vector3(0, 0, 0), // initial velocity of the electron
    pos.clone(), // position of the electron in the scene
    "electron" // type of the particle
  ));

  // Add a sphere entity to the scene to represent the electron
  addParticleEntity(STATE.particles[STATE.particles.length - 1], {
    radius, // radius of the sphere entity
    color: "blue" // color of the sphere entity
  });
}

/**
 * Creates a Hydrogen atom in the scene.
 */
function createHAtom() {
  /**
   * Resets the animation state and clears the scene.
   */
  STATE.isAnimationRunning = false;
  UI.myButton.textContent = "Start";
  clearScene();

  /**
   * Pushes a proton to the origin of the scene.
   */
  pushProton(new THREE.Vector3(0, 0, 0));

  /**
   * Pushes an electron to the position (-52.917e-12, 0, 0) in the scene.
   */
  pushElectron(new THREE.Vector3(-52.917e-12, 0, 0));

  /**
   * Sets the atomic number and label of the atom.
   */
  STATE.Z = 1;
  STATE.atomLabel = AtomModel.HAtom;
}

/**
 * Creates a H2+ cation in the scene.
 * @description This function creates a H2+ cation in the scene.
 * @param {none} - No parameters are required.
 * @return {none} - No return value is produced.
 */
function createH2Cation() {
  // Reset the animation state and clear the scene
  STATE.isAnimationRunning = false;
  UI.myButton.textContent = "Start";
  clearScene();

  // Create two protons at the positions (-52.917e-12, 0, 0) and (52.917e-12, 0, 0)
  pushProton(new THREE.Vector3(-52.917e-12, 0, 0));
  pushProton(new THREE.Vector3( 52.917e-12, 0, 0));

  // Create an electron at the position (0, 0, 0)
  pushElectron(new THREE.Vector3(0, 0, 0));

  // Set the atomic number and label of the atom
  STATE.Z = 1;
  STATE.atomLabel = AtomModel.H2Cation;
}

/**
 * Creates a H⁻ anion in the scene.
 * @description This function creates a H⁻ anion in the scene.
 * @param {none} - No parameters are required.
 * @return {none} - No return value is produced.
 */
function createHAnion() {
  // Reset the animation state and clear the scene
  STATE.isAnimationRunning = false;
  UI.myButton.textContent = "Start";
  clearScene();

  // Create a proton at the position (0, 0, 0)
  pushProton(new THREE.Vector3(0, 0, 0));

  // Create two electrons at the positions (-52.917e-12, 0, 0) and (52.917e-12, 0, 0)
  pushElectron(new THREE.Vector3(-52.917e-12, 0, 0));
  pushElectron(new THREE.Vector3( 52.917e-12, 0, 0));

  // Set the atomic number and label of the atom
  STATE.Z = 1;
  STATE.atomLabel = AtomModel.HAnion;
}

/**
 * Creates a HE⁺ cation in the scene.
 * @description This function creates a HE⁺ cation in the scene.
 * @param {none} - No parameters are required.
 * @return {none} - No return value is produced.
 */
function createHECation() {
  // Reset the animation state and clear the scene
  STATE.isAnimationRunning = false;
  UI.myButton.textContent = "Start";
  clearScene();

  // Combine 2p charge on one nucleus (same as original semantics)
  // Create a proton at the origin of the scene with a charge of 2p
  pushProton(new THREE.Vector3(0, 0, 0));
  STATE.particles[0].charge = 2 * 1.602176634e-19;
  STATE.spheres[0].setAttribute("radius", 18); // Increase the radius of the sphere to 18

  // Create an electron at the position (0.5 a0, 0, 0)
  pushElectron(new THREE.Vector3(0.5 * PHYS.a0, 0, 0));

  // Set the atomic number and label of the atom
  STATE.Z = 2;
  STATE.atomLabel = AtomModel.HECation;
}

/**
 * Creates a HE atom in the scene.
 * @description This function creates a HE atom in the scene.
 * @param {none} - No parameters are required.
 * @return {none} - No return value is produced.
 */
function createHEAtom() {
  // Reset the animation state and clear the scene
  STATE.isAnimationRunning = false;
  UI.myButton.textContent = "Start";
  clearScene();

  // Create a proton at the origin of the scene with a charge of 2p
  pushProton(new THREE.Vector3(0, 0, 0));
  STATE.particles[0].charge = 2 * 1.602176634e-19;
  STATE.spheres[0].setAttribute("radius", 18); // Increase the radius of the sphere to 18

  // Create two electrons at the positions (0.5 a0, 0, 0) and (-0.5 a0, 0, 0)
  pushElectron(new THREE.Vector3( 0.5 * PHYS.a0, 0, 0));
  pushElectron(new THREE.Vector3(-0.5 * PHYS.a0, 0, 0));

  // Set the atomic number and label of the atom
  STATE.Z = 2;
  STATE.atomLabel = AtomModel.HEAtom;
}

// default
createHAtom();

// -------------------------------------------------------------
// UI updaters & controls
// -------------------------------------------------------------
/**
 * Updates the distance and energy labels based on the current atom model
 * @description This function computes the distance between the particles and the energy
 * of the system based on the current atom model. It then updates the UI labels
 * with the computed values.
 * @param {none} - No parameters are required.
 * @return {none} - No return value is produced.
 */
function updateDistanceAndEnergy() {
  const p = STATE.particles;
  const a0 = PHYS.a0;
  switch (STATE.atomLabel) {
    /**
     * Case for H atom
     * @description Compute the distance between the electron and the proton and the
     * energy of the system. Then update the UI labels with the computed values.
     */
    case AtomModel.HAtom: {
      if (p.length < 2) return;
      const len = p[1].cartPos.distanceTo(p[0].cartPos);
      UI.Distance.textContent = (len / a0).toFixed(3) + " a0";
      const energy = (a0 ** 2 / (len ** 2)) - (2 * a0 / len);
      UI.Energy.textContent = energy.toFixed(3) + " Ry";
      break;
    }
    /**
     * Case for H2+ cation
     * @description Compute the distance between the electrons and the protons and the
     * energy of the system. Then update the UI labels with the computed values.
     */
    case AtomModel.H2Cation: {
      if (p.length < 3) return;
      const l1 = p[2].cartPos.distanceTo(p[1].cartPos);
      const l2 = p[2].cartPos.distanceTo(p[0].cartPos);
      UI.Distance.textContent = (l1 / a0).toFixed(3) + " a0";
      const energy = (a0 ** 2 / (l1 ** 2)) - (2 * a0 / l1)
                   + (a0 ** 2 / (l2 ** 2)) - (2 * a0 / l2)
                   + (2 * a0 * 0.79473 / (l1 + l2));
      UI.Energy.textContent = energy.toFixed(3) + " Ry";
      break;
    }
    /**
     * Case for H⁻ anion
     * @description Compute the distance between the electrons and the proton and the
     * energy of the system. Then update the UI labels with the computed values.
     */
    case AtomModel.HAnion: {
      if (p.length < 3) return;
      const l3 = p[1].cartPos.distanceTo(p[0].cartPos);
      const l4 = p[2].cartPos.distanceTo(p[0].cartPos);
      UI.Distance.textContent = (l3 / a0).toFixed(3) + " a0";
      const energy = (a0 ** 2 / (l3 ** 2)) - (2 * a0 / l3)
                   + (a0 ** 2 / (l4 ** 2)) - (2 * a0 / l4)
                   + (2 * a0 * 0.72644 / (l3 + l4));
      UI.Energy.textContent = energy.toFixed(3) + " Ry";
      break;
    }
    /**
     * Case for HE⁺ cation
     * @description Compute the distance between the electron and the proton and the
     * energy of the system. Then update the UI labels with the computed values.
     */
    case AtomModel.HECation: {
      if (p.length < 2) return;
      const l = p[1].cartPos.distanceTo(p[0].cartPos);
      UI.Distance.textContent = (l / a0).toFixed(3) + " a0";
      const energy = (a0 ** 2 / (l ** 2)) - (2 * 2 * a0 / l);
      UI.Energy.textContent = energy.toFixed(3) + " Ry";
      break;
    }
    /**
     * Case for HE atom
     * @description Compute the distance between the electrons and the proton and the
     * energy of the system. Then update the UI labels with the computed values.
     */
    case AtomModel.HEAtom: {
      if (p.length < 3) return;
      const l = p[1].cartPos.distanceTo(p[0].cartPos);
      UI.Distance.textContent = (l / a0).toFixed(3) + " a0";
      const energy = (2 * a0 ** 2 / (l ** 2)) - (4 * a0 * 1.95393 / l) + (a0 / l);
      UI.Energy.textContent = energy.toFixed(3) + " Ry";
      break;
    }
  }
}

/**
 * Toggles the animation state and updates the UI accordingly.
 */
function startAnimation() {
  // Toggle the animation state
  STATE.isAnimationRunning = !STATE.isAnimationRunning;

  // Update the UI button text
  UI.myButton.textContent = STATE.isAnimationRunning ? "Stop" : "Start";
}

/**
 * Toggles the theta-phi coordinate system flag.
 * @description When enabled, the animation will use the theta-phi coordinate
 * system instead of the cartesian coordinate system.
 * @param {boolean} on - Toggle state (true: enable, false: disable)
 */
function toggleThetaPhi(on) {
  STATE.isThetaPhi = !!on;
}

/**
 * Toggles the Larmor radiation flag.
 * @param {boolean} a - Toggle state (true: enable, false: disable)
 */
function toggleLamor(a) {
  STATE.isLarmor = !!a;
}

/**
 * Toggles the Angular Momentum flag.
 * @description When enabled, the animation will use the Angular Momentum
 * Ansatz instead of the linear momentum Ansatz.
 * @param {boolean} a - Toggle state (true: enable, false: disable)
 */
function changeAnsatz() {
  STATE.isAngularMomentum = !STATE.isAngularMomentum;
}

/**
 * Changes the atomic model of the simulation.
 * @param {string} modelValue - The value of the selected atomic model.
 * @description This function changes the atomic model of the simulation by calling the appropriate
 * creation function based on the selected model. If the model is unknown, it logs a warning.
 */
function changeAtomModel(modelValue) {
  switch (modelValue) {
    /**
     * Creates a H atom in the scene.
     * @description This function creates a H atom in the scene by creating one proton and one electron.
     */
    case "HAtom": createHAtom(); break;
    /**
     * Creates a H2+ cation in the scene.
     * @description This function creates a H2+ cation in the scene by creating two protons and one electron.
     */
    case "H2Cation": createH2Cation(); break;
    /**
     * Creates a H- anion in the scene.
     * @description This function creates a H- anion in the scene by creating one proton and two electrons.
     */
    case "HAnion": createHAnion(); break;
    /**
     * Creates a He+ cation in the scene.
     * @description This function creates a He+ cation in the scene by creating two protons and one electron.
     */
    case "HECation": createHECation(); break;
    /**
     * Creates a He atom in the scene.
     * @description This function creates a He atom in the scene by creating two protons and two electrons.
     */
    case "HEAtom": createHEAtom(); break;
    default: console.warn("Unknown atom model:", modelValue);
  }
}

/**
 * Sets all particle velocities to zero.
 * @description This function sets all particle velocities to zero and
 * also applies legacy nudges to certain atomic models.
 */
function setVelocitiesToZero() {
  for (let i = 0; i < STATE.particles.length; i++) {
    STATE.particles[i].velocity.set(0, 0, 0);
  }

  // Legacy nudges for certain atomic models
  switch (STATE.atomLabel) {
    case AtomModel.H2Cation:
      // H2+ cation: legacy nudge for the second proton
      if (STATE.particles[2]) {
        // Give the second proton a small velocity in the x direction
        STATE.particles[2].velocity.x = 1.14e-13;
      }
      break;
    case AtomModel.HAnion:
      // H- anion: legacy nudges for the two electrons
      if (STATE.particles[1]) {
        // Give the first electron a negative velocity in the x direction
        STATE.particles[1].velocity.x = -6.35e-14;
      }
      if (STATE.particles[2]) {
        // Give the second electron a positive velocity in the x direction
        STATE.particles[2].velocity.x = +6.35e-14;
      }
      break;
    case AtomModel.HEAtom:
      // He atom: legacy nudges for the two electrons
      if (STATE.particles[1]) {
        // Give the first electron a positive velocity in the x direction
        STATE.particles[1].velocity.x = +7.35e-13;
      }
      if (STATE.particles[2]) {
        // Give the second electron a negative velocity in the x direction
        STATE.particles[2].velocity.x = -7.35e-13;
      }
      break;
    default: break; // H Atom, HE Cation etc.: do nothing
  }
}

/**
 * Toggle the nucleus lock state.
 * When the nucleus is locked, it cannot move due to external forces.
 */
function lockNucleus() {
  /**
   * The current state of the nucleus lock.
   * @type {boolean}
   */
  STATE.isNucleusLocked = !STATE.isNucleusLocked;

  if (STATE.isNucleusLocked) {
    for (const p of STATE.particles) {
      if (p.type === "proton") {
        p.velocity.set(0, 0, 0);
      }
    }
  }
}

/**
 * Toggle path recording.
 * If the path recording is enabled, the positions of all particles will be recorded every 5 frames.
 * If the path recording is disabled, the path recording will be stopped.
 */
function startRecord() {
  /**
   * The current state of the path recording.
   * @type {boolean}
   */
  STATE.isRecordPath = !STATE.isRecordPath;

  // Update the button text
  const pathButton = UI.pathButton;
  pathButton.textContent = STATE.isRecordPath ? "⏸" : "⏵";

    // Toggle the recording LED class
  UI.led.classList.toggle("recording", STATE.isRecordPath);
}

/**
 * Deletes the recorded path.
 * @description This function deletes the recorded path by removing all the
 * trajectory points from the scene.
 */
function deletePath() {
  // only delete our points
  const points = SCENE.sceneEl.querySelectorAll("a-circle.trajectory-point");
  // iterate over all points and remove them from the scene
  points.forEach(p => {
    if (p.parentNode) {
      // remove the point from the scene
      p.parentNode.removeChild(p);
    }
  });
}

/**
 * Increment the energy level of the atom.
 * @description This function increments the energy level of the atom.
 * @returns {void}
 */
function incrementEnergyLevel() {
  if (STATE.atomLabel === AtomModel.HAtom) {
    if (STATE.isAngularMomentum && STATE.particles[1]) {
      // up-scale distance
      STATE.particles[1].cartPos.multiplyScalar(4);
      const conn = STATE.particles[0].cartPos.clone().sub(STATE.particles[1].cartPos);
      const perp = new THREE.Vector3(-conn.y, conn.x, 0).normalize().multiplyScalar(1e-12);
      STATE.particles[1].velocity.add(perp);
    }
    STATE.N++; STATE.L++;
  }
  // radiation & camera based on N
  if (STATE.N === 1) {
    STATE.radiation = 1; SCENE.cameraEl.setAttribute("position", { x: 0, y: 0, z: 200 });
  } else if (STATE.N === 2) {
    STATE.radiation = 0.999905; SCENE.cameraEl.setAttribute("position", { x: 0, y: 0, z: 250 });
  } else if (STATE.N === 3) {
    STATE.radiation = 0.99995; SCENE.cameraEl.setAttribute("position", { x: 0, y: 0, z: 400 });
  }
  UI.nLvl.textContent = String(parseInt(UI.nLvl.textContent, 10) + 1);
  UI.lLvl.textContent = String(parseInt(UI.lLvl.textContent, 10) + 1);
}

/**
 * Decrement the energy level of the atom.
 * @description This function decrements the energy level of the atom.
 * @returns {void}
 */
function decrementEnergyLevel() {
  // get the current energy level
  const currentN = parseInt(UI.nLvl.textContent, 10) || 1;

  // if the energy level is greater than 1, decrement it
  if (currentN > 1) {
    // decrement the energy level
    STATE.N--;

    // if the energy level is equal to the angular momentum level, decrement the angular momentum level
    if (STATE.L === STATE.N) {
      STATE.L--; UI.lLvl.textContent = String(parseInt(UI.lLvl.textContent, 10) - 1);
    }

    // update the UI
    UI.nLvl.textContent = String(currentN - 1);

    // update the camera position and radiation level
    switch (STATE.N) {
      case 1:
        STATE.radiation = 1;
        SCENE.cameraEl.setAttribute("position", { x: 0, y: 0, z: 200 });
        break;
      case 2:
        STATE.radiation = 0.999905;
        SCENE.cameraEl.setAttribute("position", { x: 0, y: 0, z: 250 });
        break;
      case 3:
        STATE.radiation = 0.999999;
        SCENE.cameraEl.setAttribute("position", { x: 0, y: 0, z: 400 });
        break;
      default:
        console.warn("Invalid energy level");
    }
  } else {
    console.warn("Invalid energy level");
  }
}

function incrementAngularMomentum() {
  if (STATE.L < STATE.N - 1) { STATE.L++; UI.lLvl.textContent = String(parseInt(UI.lLvl.textContent, 10) + 1); }
  else { console.warn("Invalid spin"); }
}

/**
 * Decrement the angular momentum level of the atom.
 * @description This function decrements the angular momentum level of the atom.
 * @returns {void}
 */
function decrementAngularMomentum() {
  // check if the angular momentum level is greater than 0
  if (STATE.L > 0) {
    // decrement the angular momentum level
    STATE.L--;
    // update the UI
    UI.lLvl.textContent = String(parseInt(UI.lLvl.textContent, 10) - 1);
  } else {
    // print a warning if the angular momentum level is invalid
    console.warn("Invalid spin");
  }
}

/**
 * Change the count of protons by the given delta.
 * @param {number} delta the amount to change the count by
 * @description This function changes the count of protons by the given delta.
 * If the delta is positive, a new proton is pushed to a free position.
 * If the delta is negative, the last proton is removed from the scene.
 * @returns {void}
 */
function changeCountProton(delta) {
  if (delta > 0) {
    // find a random non-colliding x position for the new proton
    const pos = findNonCollidingPosition(
      "proton",
      () => {
        const newX = Math.random() * 2e-10 - 1e-10;
        const newY = Math.random() * 2e-10 - 1e-10;
        return new THREE.Vector3(newX, newY, 0);
      }
    );

    if (pos) {
      pushProton(pos);
    }
  } else if (delta < 0) {
    // find the last proton and remove it from the scene
    for (let i = STATE.particles.length - 1; i >= 0; i--) {
      if (STATE.particles[i].type === "proton") {
        // remove the sphere and ring aligned with the proton
        const ring = STATE.bohrRings.pop(); if (ring) SCENE.sceneEl.removeChild(ring);
        const sphere = STATE.spheres[i]; if (sphere) SCENE.sceneEl.removeChild(sphere);
        // remove the proton from the particles array
        STATE.spheres.splice(i, 1);
        STATE.particles.splice(i, 1);
        break;
      }
    }
  }

  // update the UI
  UI.DispCountProton.textContent = String(
    STATE.particles.filter(p => p.type === "proton").length
  );
}

/**
 * Change the count of electrons by the given delta.
 * @param {number} delta the amount to change the count by
 * @description This function changes the count of electrons by the given delta.
 * If the delta is positive, a new electron is pushed to a free position.
 * If the delta is negative, the last electron is removed from the scene.
 * @returns {void}
 */
function changeCountElectron(delta) {
  if (delta > 0) {
    // find a random non-colliding x/y position for the new electron
    const pos = findNonCollidingPosition(
      "electron",
      () => {
        const newX = (Math.random() * 2 * PHYS.a0) - PHYS.a0;
        const newY = (Math.random() * 2 * PHYS.a0) - PHYS.a0;
        return new THREE.Vector3(newX, newY, 0);
      }
    );

    if (pos) {
      // nimm hier deinen gewünschten Radius (5, 1.5, …)
      pushElectron(pos, 5);
    }
  } else if (delta < 0) {
    // find the last electron and remove it from the scene
    for (let i = STATE.particles.length - 1; i >= 0; i--) {
      if (STATE.particles[i].type === "electron") {
        // remove the sphere from the scene
        const sphere = STATE.spheres[i]; if (sphere) SCENE.sceneEl.removeChild(sphere);
        // remove the electron from the particles array
        STATE.spheres.splice(i, 1);
        STATE.particles.splice(i, 1);
        break;
      }
    }
  }

  // update the UI
  UI.DispCountElectron.textContent = String(
    STATE.particles.filter(p => p.type === "electron").length
  );
}

function d2ThetaLnPsi(n, l, theta) {
  const c = Math.cos(theta);
  const s = Math.sin(theta);
  const c2 = c*c;
  const s2 = s*s;

  // Beispiele m = 0 (aus Tabelle 2)
  if (n === 2 && l === 1) {
    // 2p:  -1 / cos²θ
    return -1 / (c2 || 1e-9);
  }
  if (n === 3 && l === 2) {
    // aus Tabelle: 6 (sin²θ - 2) / (3 cos²θ - 1)²
    const denom = (3*c2 - 1);
    return 6 * (s2 - 2) / ((denom*denom) || 1e-9);
  }

  return 0; 
}

function d2PhiLnPsi(n, l, m, phi) {
  // in deinen Tabellen ist für m=0 überall 0
  // für m ≠ 0: später ergänzen
  return 0;
}



/**
 * Add event listeners to the position input fields.
 * The event listeners are triggered when the user presses enter in the input fields.
 * If the user presses enter, the position of the corresponding particle is updated.
 * @description This function adds event listeners to the position input fields.
 * @returns {void}
 */
["Xs", "Ys", "Zs"].forEach((key, idx) => {
  /**
   * Event listener for the position input fields.
   * The event listener is triggered when the user presses enter in the input field.
   * @param {KeyboardEvent} ev the event triggered by the user
   * @description This event listener is triggered when the user presses enter in the input field.
   * If the user presses enter, the position of the corresponding particle is updated.
   * @returns {void}
   */
  UI[key].addEventListener("keypress", (ev) => {
    if (ev.key !== "Enter") return;
    const value = parseFloat(UI[key].value);
    if (!Number.isFinite(value)) return;
    const p = STATE.particles[STATE.indexOfParticle];
    if (!p) return;
    if (idx === 0) p.cartPos.x = value * 1e-12;
    if (idx === 1) p.cartPos.y = value * 1e-12;
    if (idx === 2) p.cartPos.z = value * 1e-12;
  });
});

// -------------------------------------------------------------
// Physics step
// -------------------------------------------------------------

/**
 * Animation loop (physics step).
 * This function is called every STATE.dtMs milliseconds.
 * It updates the positions and velocities of all particles in the scene.
 * @description This function updates the positions and velocities of all particles in the scene.
 * @returns {void}
 */
function animation() {
  updateDistanceAndEnergy();

  if (!STATE.isAnimationRunning || STATE.particles.length === 0) return;

  const pArr = STATE.particles;
  const Lterm = STATE.isAngularMomentum ? 0 : (STATE.L * (STATE.L + 1));

  for (let i = 0; i < pArr.length; i++) {
    const pi = pArr[i];

    // --- Nucleus-Lock: Protonen bleiben stehen (v = 0, a = 0) ---
    if (STATE.isNucleusLocked && pi.type === "proton") {
      pi.velocity.set(0, 0, 0);
      // Position bleibt einfach so, wie sie ist
      continue; // keine Kräfte auf Proton integrieren
    }

    // Gesamtbeschleunigung aus allen Wechselwirkungen
    let acc = new THREE.Vector3(0, 0, 0);
    let accSpin = new THREE.Vector3(0, 0, 0);

    // Für das UI / Periodenmessung merken wir uns den letzten dv / phi der Elektronen
    let lastDvForUI = null;
    let lastPhiForUI = 0;

    for (let j = 0; j < pArr.length; j++) {
      if (i === j) continue;
      const pj = pArr[j];

      const dv = new THREE.Vector3().copy(pj.cartPos).sub(pi.cartPos);
      const [r, theta, phi] = cartInPolar(dv.x, dv.y, dv.z);
      if (r === 0) continue;

      const perp = new THREE.Vector3(-dv.y, dv.x, 0);

      // --- Coulomb + Kratzer + evtl. Spin-Anteil ---
      let Fel = (-STATE.Z * pi.charge * pj.charge) / (K * r * r);
      let Fk = 0;
      let Fspin = 0;
      if (pi.charge * pj.charge < 0) {
        Fk = (pi.charge * pj.charge * PHYS.a0) / (K * r * r * r);
        if (Lterm) {
          Fspin = (pi.charge * pj.charge * Lterm * PHYS.a0) / (K * r * r * r);
        }
      }

      // const aDir = dv.clone().normalize();
      // const aPair = aDir.multiplyScalar((Fel + Fk) / pi.mass);
      // acc.add(aPair);

      // if (Lterm) {
      //   const aSpinPair = perp.clone().normalize().multiplyScalar(Fspin / pi.mass);
      //   accSpin.add(aSpinPair);
      // }

      // Richtung vom Kraftvektor
      const aDir = dv.clone().normalize();

      // θ/φ-Quantenterme aus Tabellen
      let FthetaPhi = 0;
      if (STATE.isThetaPhi && pi.type === "electron" && pj.type === "proton") {
        const d2θ = d2ThetaLnPsi(STATE.N, STATE.L, theta);
        const d2φ = d2PhiLnPsi(STATE.N, STATE.L, /*m=*/0, phi); // m=0 bisher

        const sin2 = Math.sin(theta); 
        const sin2sq = sin2*sin2 || 1e-9;
        const A = d2θ + d2φ / sin2sq;

        // Radialkraft aus V_theta,phi:
        // F_r = 2 A / r^3
        FthetaPhi = -2 * A / (r*r*r);
      }

      // Gesamtkraft (radial) + Larmor/Kratzer
      const Ftotal = Fel + Fk + FthetaPhi;

      // wie gehabt:
      acc = aDir.multiplyScalar(Ftotal / pi.mass);
      accSpin = ( Lterm ? 
        perp.clone().normalize().multiplyScalar(Fspin / pi.mass): 
        new THREE.Vector3(0,0,0)
      );


      // Für die Anzeige (nur Elektronen)
      if (pi.charge < 0) {
        lastDvForUI = dv;
        lastPhiForUI = phi;
      }
    }

    // --- Larmor-Strahlung: Energieverlust ~ v² ---
    if (STATE.isLarmor && pi.type === "electron") {
      // "Rate" pro Millisekunde 
      const LARMOR_RATE = 1e-5; // ≈ 0.01 pro 100 ms
      const damping = Math.exp(-LARMOR_RATE * STATE.dtMs);
      pi.velocity.multiplyScalar(damping);
    }

    // --- Integration (ähnliche Größenordnung wie vorher) ---
    pi.velocity.add(acc.multiplyScalar(1e-35)).add(accSpin.multiplyScalar(1e-37));
    pi.cartPos.add(pi.velocity.clone().multiplyScalar(STATE.radiation));

    // --- Kraft-Readout & Periodenmessung für Elektronen ---
    if (pi.charge < 0 && lastDvForUI) {
      const forceVec = acc.clone().multiplyScalar(pi.mass);
      UI.fx.textContent = String(forceVec.x);
      UI.fy.textContent = String(forceVec.y);
      UI.fz.textContent = String(forceVec.z);

      const dv = lastDvForUI;
      const phi = lastPhiForUI;

      if (STATE.countTime === 1) {
        if (Number(STATE.angleCache.toFixed(2)) === Number(phi.toFixed(2))) {
          STATE.countTime = 2;
        }
        STATE.timePassedPs += 2 * 1.6 / 885; // legacy constant
      }
      if (dv.length() >= 4 * PHYS.a0 && STATE.countTime === 0) {
        STATE.countTime = 1;
        STATE.angleCache = phi;
      }
    }
  }
}

// physics timer (keeps original dt/select behavior)
let interval = setInterval(animation, STATE.dtMs);
/**
 * Changes the interval of the animation loop.
 * @param {number} newDt - new interval in milliseconds
 */
function changeInterval(newDt) {
  const v = Number(newDt);
  if (!Number.isFinite(v) || v <= 0) return;

  // update the interval of the animation loop
  STATE.dtMs = v;

  // clear the old interval and set a new one
  clearInterval(interval);
  interval = setInterval(animation, STATE.dtMs);
}

// -------------------------------------------------------------
// Render loop (rAF), path recording & frames counter hook
// -------------------------------------------------------------

/**
 * Called every frame by requestAnimationFrame.
 * Updates the trajectory points and syncs particle positions to entities.
 * @description This function updates the trajectory points and syncs particle positions to entities.
 * @returns {void}
 */
function intoRealWorld() {
  // global `frames` variable from fps.js
  if (typeof frames === "number") frames++; // increment frames counter

  if (STATE.particles.length > 0) {
    // sync particle positions to entities
    syncEntitiesToParticles();

    // record trajectory points
    if (STATE.isRecordPath && (typeof frames === "number") && frames % 5 === 0) {
      for (let i = 0; i < STATE.particles.length; i++) {
        if (STATE.particles[i].type !== "electron") continue;
        const p = STATE.particles[i].cartPos.clone().multiplyScalar(1e12);
        const point = createEntity("a-circle", { radius: 1, color: "grey" });
        point.classList.add("trajectory-point");
        point.object3D.position.set(p.x, p.y, p.z);
        SCENE.sceneEl.appendChild(point);
      }
    }
  }
  requestAnimationFrame(intoRealWorld);
}
requestAnimationFrame(intoRealWorld);

// -------------------------------------------------------------
// Public API (keeps particles.html working)
// -------------------------------------------------------------
window.cartInPolar = cartInPolar;
window.startAnimation = startAnimation;
window.toggleLamor = toggleLamor;
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
window.animation = animation; // still available (interval calls this)
window.changeInterval = changeInterval;
window.intoRealWorld = intoRealWorld;
window.incrementEnergyLevel = incrementEnergyLevel;
window.decrementEnergyLevel = decrementEnergyLevel;
window.incrementAngularMomentum = incrementAngularMomentum;
window.decrementAngularMomentum = decrementAngularMomentum;
window.changeAnsatz = changeAnsatz;
