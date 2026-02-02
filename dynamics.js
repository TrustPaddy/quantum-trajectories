// Parameter-Objekt für alle Controls
const settings = {
  running: false,
  velocitiesZero: () => setVelocitiesToZero(),
  lockNucleus: false,
  recordPath: false,
  deletePath: () => deletePath(),
  energyLevel: 1,
  angularMomentum: 0,
  atom: "HAtom",
  electronCount: 1,
  protonCount: 1,
  larmor: false,
  spin: false,
  dt: 33.33,
  posX: 0,
  posY: 0,
  posZ: 0
};

const gui = new dat.GUI({ width: 350 });

const animFolder = gui.addFolder("Animation");
animFolder.add(settings, "running").name("Running").onChange(value => {
  if (value) {
    startAnimation();
  } else {
    // stopAnimation() 
  }
});
animFolder.add(settings, "velocitiesZero").name("v = 0");


const atomFolder = gui.addFolder("Atom Model");
atomFolder.add(settings, "atom", {
  "H Atom": "HAtom",
  "H2+ Cation": "H2Cation",
  "H- Anion": "HAnion",
  "He+ Cation": "HECation",
  "He Atom": "HEAtom"
}).name("Atom").onChange(value => {
  changeAtomModel(value);
});
atomFolder.add(settings, "electronCount", 1, 10, 1).name("Electrons").onChange(val => {
  changeCountElectron(val - settings.electronCount); 
});
atomFolder.add(settings, "protonCount", 1, 10, 1).name("Protons").onChange(val => {
  changeCountProton(val - settings.protonCount); 
});

const quantumFolder = gui.addFolder("Quantum States");
// quantumFolder.add(settings, "energyLevel", 1, 10, 1).name("Energy Level").onChange(val => {
//   document.getElementById("nLvl").innerText = val;
// });
// quantumFolder.add(settings, "angularMomentum", 0, 10, 1).name("Angular Momentum").onChange(val => {
//   document.getElementById("lLvl").innerText = val;
// });
quantumFolder.add(settings, "energyLevel", 1, 10).step(1).name("Energy Level");
quantumFolder.add(settings, "angularMomentum", 0, 10).step(1).name("Angular Momentum");

const physicsFolder = gui.addFolder("Physics Options");
physicsFolder.add(settings, "lockNucleus").name("Lock Nucleus");
physicsFolder.add(settings, "larmor").name("Larmor Precession");
physicsFolder.add(settings, "spin").name("Spin");

const pathFolder = gui.addFolder("Path Recording");
pathFolder.add(settings, "recordPath").name("Record Path").onChange(value => {
  if (value) {
    startRecord();
  } else {
    // evtl. stopRecord()
  }
});
pathFolder.add(settings, "deletePath").name("Delete Path");

const speedFolder = gui.addFolder("Simulation Speed");
speedFolder.add(settings, "dt", { Normal: 33.33, Fast: 6.67, Slow: 166.67 }).name("dt").onChange(value => {
  changeInterval(value);
});

const posFolder = gui.addFolder("Initial Position");
posFolder.add(settings, "posX", -100, 100, 1).name("X").onChange(val => {
  document.getElementById("Xs").value = val;
});
posFolder.add(settings, "posY", -100, 100, 1).name("Y").onChange(val => {
  document.getElementById("Ys").value = val;
});
posFolder.add(settings, "posZ", -100, 100, 1).name("Z").onChange(val => {
  document.getElementById("Zs").value = val;
});

animFolder.close();
atomFolder.close();
quantumFolder.close();
physicsFolder.close();
pathFolder.close();
speedFolder.close();
posFolder.close();


const sceneEl = document.querySelector('a-scene');
const cameraEl = document.querySelector('#camera');
const button = document.getElementById("myButton");
const pathButton = document.getElementById("PathButton");
const led = document.getElementById('led');
            
const particles = [];

const particleEl = []; 

const bohrRadius = [];

let atom = "H Atom";
    
const e0 = 8.8541878128e-12;
const a0 = 5.29177210903e-11;
const c = 2.99792e8;
const k = 4 * Math.PI * e0;

let isAnimationRunning = false;

let isLamor = false;

const inputField = {
    X: document.getElementById("Xs"),
    Y: document.getElementById("Ys"),
    Z: document.getElementById("Zs"),
    Veff: document.getElementById("Veff"),
}

let indexOfParticle;

let countProton = 1;
let countElectron = 1;

let dt = 1000/30;

let Z = 1;

let isNucleusLocked = false;

const RadConst = 2/(3*c*c*c)

let N = 1;

let L = 0;

let isAngularMomentum = 0;

let radiation = 1;

let timePassed = 0;

let countTime = 0;

let angle = 0;

let isRecordPath = false;

/**
 * Toggles the animation state between running and stopped.
 * Updates the button text and LED recording indicator accordingly.
 */
function startAnimation() {
    isAnimationRunning = !isAnimationRunning;
    button.textContent = isAnimationRunning ? "Stop" : "Start";
    led.classList.toggle('recording', isAnimationRunning);
}



/**
 * Toggles the Larmor precession on or off.
 * @param {boolean} a - true to enable Larmor precession, false to disable it
 */
function toggleLamor(a) {
    isLamor = a;
}



/**
 * Converts a 3D cartesian coordinate (x, y, z) to a 3D polar coordinate (r, θ, φ).
 * @param {Number} x - x-coordinate in the cartesian system
 * @param {Number} y - y-coordinate in the cartesian system
 * @param {Number} z - z-coordinate in the cartesian system
 * @returns {Array<number>} [r, θ, φ] - polar coordinates (radius, theta, phi)
 */
function cartInPolar(x, y, z) {
    const r = Math.sqrt(x * x + y * y + z * z);
    const θ = Math.acos(z/r);
    const φ = Math.atan2(y, x);
    return [r, θ , φ];
}
    


// new class of particles with the attributes of particles
class Particle {
    /**
     * Creates a new particle with the given properties.
     * @param {number} mass - mass of the particle
     * @param {number} charge - charge of the particle
     * @param {number} spin - spin of the particle
     * @param {THREE.Vector3} velocity - velocity of the particle
     * @param {THREE.Vector3} cartPos - position of the particle in cartesian coordinates
     * @param {string} type - type of the particle (proton or electron)
     */
    constructor(mass, charge, spin, velocity, cartPos, type) {
    this.mass = mass;
    this.charge = charge;  
    this.spin = spin;
    this.velocity = velocity;
    this.cartPos = cartPos;
    this.polarPos = cartInPolar(cartPos.x, cartPos.y, cartPos.z);
    this.type = type;
    this.index; // index of the particle
    this.n; 
    this.l;
    this.m_l;
    this.m_s;
    }
}



/**
 * Creates a new atom model for a Hydrogen 2 Cation (H2+).
 * @function createH2Cation
 */
function createH2Cation() {
    isAnimationRunning = false;
    button.textContent = "Start";

    // Reset the particle list
    while (particles.length > 0) {
    particles.pop();
    sceneEl.removeChild(particleEl[particleEl.length-1]);
    particleEl.pop();     
    }

    // Reset the bohr radius ring list
    while (bohrRadius.length > 0) {
    sceneEl.removeChild(bohrRadius[bohrRadius.length-1]);
    bohrRadius.pop();
    }

    // Create the first proton
    particles[0] = new Particle(
        1.6726219e-27, 
        1.602176634e-19, 
        +0.5,
        new THREE.Vector3(0, 0, 0), 
        new THREE.Vector3(-52.917e-12, 0 ,0),
        "proton"
    );
    particleEl[0] = document.createElement('a-sphere'); 
    particleEl[0].setAttribute('radius', 14);
    particleEl[0].setAttribute('color', 'red');
    sceneEl.appendChild(particleEl[0]); 
    bohrRadius[0] = document.createElement('a-ring');   
    bohrRadius[0].setAttribute('radius-inner', a0 * 9.99e11);
    bohrRadius[0].setAttribute('radius-outer', a0 * 1.01e12);
    bohrRadius[0].setAttribute('color', 'grey'); 
    sceneEl.appendChild(bohrRadius[0]);

    // Create the second proton
    particles[1] = new Particle(
        1.6726219e-27, 
        1.602176634e-19, 
        0.5, 
        new THREE.Vector3(0, 0, 0), 
        new THREE.Vector3(52.917e-12, 0, 0),
        "proton"
    ); 
    particleEl[1] = document.createElement('a-sphere'); 
    particleEl[1].setAttribute('radius', 14);
    particleEl[1].setAttribute('color', 'red');
    sceneEl.appendChild(particleEl[1]); 
    bohrRadius[1] = document.createElement('a-ring');   
    bohrRadius[1].setAttribute('radius-inner', a0 * 9.99e11);
    bohrRadius[1].setAttribute('radius-outer', a0 * 1.01e12);
    bohrRadius[1].setAttribute('color', 'grey'); 
    sceneEl.appendChild(bohrRadius[1]);

    // Create the electron
    particles[2] = new Particle(
    9.10938356e-31, 
    -1.602176634e-19, 
    0.5, 
    new THREE.Vector3(0, 0, 0), 
    new THREE.Vector3(0, 0, 0),
    "electron"
    )
    particleEl[2] = document.createElement('a-sphere');
    particleEl[2].setAttribute('radius', 5);
    particleEl[2].setAttribute('color', 'blue');
    sceneEl.appendChild(particleEl[2]);

    // Add click event listeners to the particles
    for (let i in particleEl) {
    particleEl[i].index = i;
    particleEl[i].addEventListener('click', (event) => {
        indexOfParticle = particleEl[i].index;
        document.getElementById("Xs").value = particleEl[i].getAttribute('position').x.toString();
        document.getElementById("Ys").value = particleEl[i].getAttribute('position').y.toString();
        document.getElementById("Zs").value = particleEl[i].getAttribute('position').z.toString();
    });               
    } 

    // Set the atomic number (Z) to 1
    Z = 1;
    // Set the atom name to "H2 Cation"
    atom = "H2 Cation"
}

/**
 * Creates an H anion
 */
function createHAnion() {
    isAnimationRunning = false;
    button.textContent = "Start";

    // Remove all particles and Bohr radius rings from the scene
    while (particles.length > 0) {
    particles.pop();
    sceneEl.removeChild(particleEl[particleEl.length-1]);
    particleEl.pop();     
    }

    while (bohrRadius.length > 0) {
    sceneEl.removeChild(bohrRadius[bohrRadius.length-1]);
    bohrRadius.pop();
    }

    // Create the proton
    particles[0] = new Particle(
    1.6726219e-27, 
    1.602176634e-19, 
    0.5, 
    new THREE.Vector3(0, 0, 0), 
    new THREE.Vector3(0, 0, 0),
    "proton"
    );
    const newProton = document.createElement('a-sphere'); 
    newProton.setAttribute('radius', 14);
    newProton.setAttribute('color', 'red');
    particleEl[0] = newProton;   
    newProton.object3D.position.set(0, 0, 0);   
    sceneEl.appendChild(newProton); 
    const newBohr = document.createElement('a-ring');   
    newBohr.setAttribute('radius-inner', a0 * 9.99e11);
    newBohr.setAttribute('radius-outer', a0 * 1.01e12);
    newBohr.setAttribute('color', 'grey');
    bohrRadius[0] = newBohr;  
    newBohr.object3D.position.set(0, 0, 0);    
    sceneEl.appendChild(newBohr);

    // Create the two electrons
    particles[1] = new Particle(
    9.10938356e-31, 
    -1.602176634e-19, 
    0.5, 
    new THREE.Vector3(0, 0, 0), 
    new THREE.Vector3(-52.917e-12, 0, 0),
    "electron"
    )   
    particleEl[1] = document.createElement('a-sphere');
    particleEl[1].setAttribute('radius', 5);
    particleEl[1].setAttribute('color', 'blue');
    particleEl[1].object3D.position.set(5.2918, 0, 0); 
    sceneEl.appendChild(particleEl[1]);

    particles[2] = new Particle(
    9.10938356e-31, 
    -1.602176634e-19, 
    0.5, 
    new THREE.Vector3(0, 0, 0), 
    new THREE.Vector3(52.917e-12, 0, 0),
    "electron"
    );
    particleEl[2] = document.createElement('a-sphere');
    particleEl[2].setAttribute('radius', 5);
    particleEl[2].setAttribute('color', 'blue');
    particleEl[2].object3D.position.set(5.2918, 0, 0); 
    sceneEl.appendChild(particleEl[2]);

    // Add event listeners to the particles
    for (let i in particleEl) {
    particleEl[i].index = i;
    particleEl[i].addEventListener('click', (event) => {
        indexOfParticle = particleEl[i].index;
        document.getElementById("Xs").value = particleEl[i].getAttribute('position').x.toString();
        document.getElementById("Ys").value = particleEl[i].getAttribute('position').y.toString();
        document.getElementById("Zs").value = particleEl[i].getAttribute('position').z.toString();
    });               
    } 
    // Set the atomic number (Z) to 1
    Z= 1;
    // Set the atom name to "H Anion"
    atom = "H Anion"
}

/**
 * Creates a new atom model for a Hydrogen Atom (H).
 * @function createHAtom
 */
function createHAtom() {
    isAnimationRunning = false;
    button.textContent = "Start";

    // Reset the particle list
    while (particles.length > 0) {
    particles.pop();
    sceneEl.removeChild(particleEl[particleEl.length-1]);
    particleEl.pop();     
    }

    // Reset the bohr radius ring list
    while (bohrRadius.length > 0) {
    sceneEl.removeChild(bohrRadius[bohrRadius.length-1]);
    bohrRadius.pop();
    }

    // Create the proton
    particles[0] = new Particle(
    1.6726219e-27, 
    1.602176634e-19, 
    0.5, 
    new THREE.Vector3(0, 0, 0),
    new THREE.Vector3(0, 0, 0),
    "proton"
    );
    const newProton = document.createElement('a-sphere'); 
    newProton.setAttribute('radius', 14);
    newProton.setAttribute('color', 'red');
    particleEl[0] = newProton;   
    newProton.object3D.position.set(0, 0, 0);   
    sceneEl.appendChild(newProton); 
    const newBohr = document.createElement('a-ring');   
    newBohr.setAttribute('radius-inner', a0 * 9.99e11);
    newBohr.setAttribute('radius-outer', a0 * 1.01e12);
    newBohr.setAttribute('color', 'grey');
    bohrRadius[0] = newBohr;  
    newBohr.object3D.position.set(0, 0, 0);    
    sceneEl.appendChild(newBohr);

    // Create the electron
    particles[1] = new Particle(
    9.10938356e-31, 
    -1.602176634e-19, 
    0.5, 
    new THREE.Vector3(0, 0, 0), 
    new THREE.Vector3(-52.917e-12, 0, 0),
    "electron"
    ) 
    particleEl[1] = document.createElement('a-sphere');
    particleEl[1].setAttribute('radius', 5);
    particleEl[1].setAttribute('color', 'blue');
    particleEl[1].object3D.position.set(5.2918, 0, 0); 
    sceneEl.appendChild(particleEl[1]);

    // Add event listeners to the particles
    for (let i in particleEl) {
        particleEl[i].index = i;
        particleEl[i].addEventListener('click', (event) => {
        indexOfParticle = particleEl[i].index;
        document.getElementById("Xs").value = particleEl[i].getAttribute('position').x.toString();
        document.getElementById("Ys").value = particleEl[i].getAttribute('position').y.toString();
        document.getElementById("Zs").value = particleEl[i].getAttribute('position').z.toString();
        });               
    } 

    // Set the atomic number (Z) to 1
    Z= 1;
    // Set the atom name to "H Atom"
    atom = "H Atom"
}

/**
 * Creates a Helium Cation (HE+)
 */
function createHECation() {
    isAnimationRunning = false;
    button.textContent = "Start";

    // Remove all particles and Bohr radius rings from the scene
    while (particles.length > 0) {
    particles.pop();
    sceneEl.removeChild(particleEl[particleEl.length - 1]);
    particleEl.pop();
    }

    while (bohrRadius.length > 0) {
    sceneEl.removeChild(bohrRadius[bohrRadius.length - 1]);
    bohrRadius.pop();
    }

    // Create the two protons
    particles[0] = new Particle(
    1.6726219e-27, 
    2 * 1.602176634e-19,
    0.5, 
    new THREE.Vector3(0, 0, 0), 
    new THREE.Vector3(0, 0 ,0),
    "proton"
    );
    particles[0].index = 1;   
    particleEl[0] = document.createElement('a-sphere'); 
    particleEl[0].setAttribute('radius', 18);
    //particleEl[0].setAttribute('scale', { x: 1, y: 0.5, z: 0.5 });
    particleEl[0].setAttribute('color', 'red');
    sceneEl.appendChild(particleEl[0]); 
    bohrRadius[0] = document.createElement('a-ring');   
    bohrRadius[0].setAttribute('radius-inner', a0 * 9.99e11);
    bohrRadius[0].setAttribute('radius-outer', a0 * 1.01e12);
    bohrRadius[0].setAttribute('color', 'grey');  
    sceneEl.appendChild(bohrRadius[0]);

    // Create the electron
    particles[1] = new Particle(
    9.10938356e-31, 
    -1.602176634e-19, 
    0.5, 
    new THREE.Vector3(0, 0, 0), 
    new THREE.Vector3(0.5 * a0, 0, 0),
    "electron"
    );
    particles[1].index = 3;
    particleEl[1] = document.createElement('a-sphere');
    particleEl[1].setAttribute('radius', 5);
    particleEl[1].setAttribute('color', 'blue');
    particleEl[1].object3D.position.set(5.2918, 0, 0); 
    sceneEl.appendChild(particleEl[1]);

    // Add event listeners to the particles
    for (let i in particleEl) {
    particleEl[i].index = i;
    particleEl[i].addEventListener('click', (event) => {
        indexOfParticle = particleEl[i].index;
        document.getElementById("Xs").value = particleEl[i].getAttribute('position').x.toString();
        document.getElementById("Ys").value = particleEl[i].getAttribute('position').y.toString();
        document.getElementById("Zs").value = particleEl[i].getAttribute('position').z.toString();
    });
    }

    // Set the atomic number (Z) to 2
    Z = 2;
    // Set the atom name to "HE Cation"
    atom = "HE Cation";
}

/**
 * Creates a Helium Atom (HE)
 * Initializes the scene with two protons and two electrons
 * @function createHEAtom
 */
function createHEAtom() {
    // Stop animation and reset button text
    isAnimationRunning = false;
    button.textContent = "Start";

    // Clear existing particles and Bohr radius rings from the scene
    while (particles.length > 0) {
    particles.pop();
    sceneEl.removeChild(particleEl[particleEl.length - 1]);
    particleEl.pop();
    }

    while (bohrRadius.length > 0) {
    sceneEl.removeChild(bohrRadius[bohrRadius.length - 1]);
    bohrRadius.pop();
    }

    // Create the first proton
    particles[0] = new Particle(
    1.6726219e-27,
    2 * 1.602176634e-19,
    0.5,
    new THREE.Vector3(0, 0, 0),
    new THREE.Vector3(0, 0, 0),
    "proton"
    );
    particles[0].index = 1;
    particleEl[0] = document.createElement('a-sphere');
    particleEl[0].setAttribute('radius', 18);
    particleEl[0].setAttribute('color', 'red');
    sceneEl.appendChild(particleEl[0]);

    // Create Bohr radius ring for the proton
    bohrRadius[0] = document.createElement('a-ring');
    bohrRadius[0].setAttribute('radius-inner', a0 * 9.99e11);
    bohrRadius[0].setAttribute('radius-outer', a0 * 1.01e12);
    bohrRadius[0].setAttribute('color', 'grey');
    sceneEl.appendChild(bohrRadius[0]);

    // Create the first electron
    particles[1] = new Particle(
    9.10938356e-31,
    -1.602176634e-19,
    0.5,
    new THREE.Vector3(0, 0, 0),
    new THREE.Vector3(0.5 * a0, 0, 0),
    "electron"
    );
    particles[1].index = 3;
    particleEl[1] = document.createElement('a-sphere');
    particleEl[1].setAttribute('radius', 5);
    particleEl[1].setAttribute('color', 'blue');
    particleEl[1].object3D.position.set(5.2918, 0, 0);
    sceneEl.appendChild(particleEl[1]);

    // Create the second electron
    particles[2] = new Particle(
    9.10938356e-31,
    -1.602176634e-19,
    0.5,
    new THREE.Vector3(0, 0, 0),
    new THREE.Vector3(-0.5 * a0, 0, 0),
    "electron"
    );
    particles[2].index = 3;
    particleEl[2] = document.createElement('a-sphere');
    particleEl[2].setAttribute('radius', 5);
    particleEl[2].setAttribute('color', 'blue');
    particleEl[2].object3D.position.set(5.2918, 0, 0);
    sceneEl.appendChild(particleEl[2]);

    // Add event listeners to the particles for interaction
    for (let i in particleEl) {
    particleEl[i].index = i;
    particleEl[i].addEventListener('click', (event) => {
        indexOfParticle = particleEl[i].index;
        document.getElementById("Xs").value = particleEl[i].getAttribute('position').x.toString();
        document.getElementById("Ys").value = particleEl[i].getAttribute('position').y.toString();
        document.getElementById("Zs").value = particleEl[i].getAttribute('position').z.toString();
    });
    }

    // Set atomic number and atom name
    Z = 2;
    atom = "HE Atom";
}

// H-atom as default
createHAtom();



/**
 * Change the atom model
 * @param {String} atomModel - HAtom, H2Cation, HAnion, HECation, HEAtom
 */
function changeAtomModel(atomModel) {
    // Change the atom model based on the option selected
    switch (atomModel) {
    case "HAtom":
        createHAtom();
        break;
    case "H2Cation":
        createH2Cation();
        break;
    case "HAnion":
        createHAnion();
        break;
    case "HECation":
        createHECation();
        break;
    case "HEAtom":
        createHEAtom();
        break;
    default:
        console.log("Invalid atom model");
    }
}



/**
 * Change the number of protons in the atom
 * @param {Number} a - the change in the number of protons
 */
function changeCountProton(a) {
    // if the number of protons is increasing
    if (a > 0) {
    // increase the number of protons
    countProton++;
    // generate a random position for the new proton
    let newX = (Math.random() * 1e-11);  
    // create a new particle
    particles.push( 
        new Particle(
        1.6726219e-27, 
        1.602176634e-19, 
        0.5, new THREE.Vector3(0, 0, 0), 
        new THREE.Vector3(0, 0, 0),
        "proton"
        )
    ); 
    // set the position of the new particle
    particles[particles.length-1].cartPos.x = newX;
    // convert the position to polar coordinates
    const p = particles[particles.length-1].cartPos
    particles[particles.length-1].polarPos = cartInPolar(p.x, p.y, p.z);

    // create a new sphere to represent the proton
    const newProton = document.createElement('a-sphere');        
    particleEl.push(newProton);
    sceneEl.appendChild(newProton);
    // set the radius and color of the sphere
    newProton.setAttribute('radius', 6);
    newProton.setAttribute('color', 'red');

    // create a new ring to represent the Bohr radius
    const newBohr = document.createElement('a-ring');        
    bohrRadius.push(newBohr);
    sceneEl.appendChild(newBohr);
    // set the radius and color of the ring
    newBohr.setAttribute('radius-inner', a0 * 9.99e11);
    newBohr.setAttribute('radius-outer', a0 * 1.01e12);
    newBohr.setAttribute('color', 'grey');
    newBohr.setAttribute('opacity', 0.2);
    newBohr.setAttribute('ignore-ray', true);

    // set the position of the sphere and ring
    const pr = p.clone().multiplyScalar(1e12);

    newProton.object3D.position.set(pr.x, pr.y, pr.z);  
    newBohr.object3D.position.set(pr.x, pr.y, pr.z);

    // if the number of protons is decreasing
    } else if (a < 0 && countProton > 0) {
    // decrease the number of protons
    countProton--;
    // remove the last proton
    const removedProton = bohrRadius.pop();
    sceneEl.removeChild(removedProton);
    // remove the corresponding particle
    for (let i = particles.length -1; i >=0; i--){
        if (particles[i].type === "proton"){
        particles.splice(i, 1);
        sceneEl.removeChild(particleEl[i])
        particleEl.splice(i, 1);
        break;
        }
    }
    }

    // update the display
    particleEl.forEach((element, i) => {
        element.index = i;
        element.addEventListener('click', () => {
        indexOfParticle = element.index;
        const position = element.getAttribute('position');
        document.getElementById("Xs").value = position.x.toString();
        document.getElementById("Ys").value = position.y.toString();
        document.getElementById("Zs").value = position.z.toString();
        });               
    });
    
    // update the display of the number of protons
    document.getElementById("DispCountProton").innerHTML = countProton;
}    

/**
 * Changes the number of electrons in the model.
 * @param {number} a Change in the number of electrons.
 */
function changeCountElectron(a) {
    if (a > 0) {
        // increase the number of electrons
        countElectron++;
        // generate a random position for the new electron
        let newX = (Math.random() * 2 *a0) - (a0);  
        let newY = (Math.random() * 2 *a0) - (a0); 

        // create a new electron
        particles.push(
            new Particle(
                9.10938356e-31, 
                -1.602176634e-19, 
                0.5, 
                new THREE.Vector3(0, 0, 0), 
                new THREE.Vector3(3, 3, 0),
                "electron"
            )
        );
        // set the new electron's position
        particles[particles.length-1].cartPos.x = newX;
        particles[particles.length-1].cartPos.y = newY;
        const p = particles[particles.length-1].cartPos;
        // set the new electron's polar position
        particles[particles.length-1].polarPos = cartInPolar(p.x, p.y, p.z);

        // create a new sphere for the electron
        const newElectron = document.createElement('a-sphere');
        particleEl.push(newElectron);
        sceneEl.appendChild(newElectron);
        newElectron.setAttribute('radius', 1.5);
        newElectron.setAttribute('color', 'blue');
        const pr = p.clone().multiplyScalar(1e12)
        newElectron.object3D.position.set(pr.x, pr.y, pr.z);      
    } 

    else if (a < 0 && countElectron > 0) {
        // decrease the number of electrons
        countElectron--;
        // remove the last electron
        for (let i = particles.length -1; i >=0; i--){
            if (particles[i].type === "electron"){
                particles.splice(i, 1);
                sceneEl.removeChild(particleEl[i])
                particleEl.splice(i, 1);
                break;
            }
        }
    }

    // update the display of the particles
    particleEl.forEach((el, i) => {
        el.index = i;
        el.addEventListener('click', () => {
            indexOfParticle = i;
            const position = el.object3D.position;
            document.getElementById("Xs").value = position.x.toString();
            document.getElementById("Ys").value = position.y.toString();
            document.getElementById("Zs").value = position.z.toString();
            console.log(particles[i]);
        });
    });

    // update the display of the number of electrons
    document.getElementById("DispCountElectron").innerHTML = countElectron;
}  



// Updates the x-position of the particle at the current indexp
inputField.X.addEventListener('keypress', (event) => {
    if (event.key === 'Enter') {
        const posx = inputField.X.value;
        particles[indexOfParticle].cartPos.x = posx * 1e-12;
    }
});

// Updates the y-position of the particle at the current indexp
inputField.Y.addEventListener('keypress', (event) => {
    if (event.key === 'Enter') {
        const posy = inputField.Y.value;
        particles[indexOfParticle].cartPos.y = posy * 1e-12;
    }
});

// Updates the z-position of the particle at the current indexp
inputField.Z.addEventListener('keypress', (event) => {
    if (event.key === 'Enter') {
        const posz = inputField.Z.value;
        particles[indexOfParticle].cartPos.z = posz * 1e-12;
    }
});

/*  inputField.Veff.addEventListener('keypress', (event) => {
    if (event.key === 'Enter') {
        const Veff = parseFloat(inputField.Veff.value);

        if(atom === "H Atom"){
            particles[1].cartPos.x = (-a0 * (Math.sqrt(Veff+1)+1))/Veff;
            //particles[1].cartPos.x = (a0 * (Math.sqrt(Veff+1)-1))/Veff;
        }

        else if(atom === "H2 Cation"){
            particles[0].cartPos.x = (-a0 * 1.4141* (Math.sqrt(Veff+1.2842)+1.1332))/Veff;
            particles[1].cartPos.x = (a0 * 1.4141* (Math.sqrt(Veff+1.2842)+1.1332))/Veff;
            //particles[0].cartPos.x = (-a0 * 1.4141* (Math.sqrt(Veff+1.2842)-1.1332))/Veff;
            //particles[1].cartPos.x = (a0 * 1.4141* (Math.sqrt(Veff+1.2842)-1.1332))/Veff;
        }

        else if(atom === "H Anion"){
            particles[2].cartPos.x = (-a0 * 1.4141* (Math.sqrt(Veff+1.3169)+1.1476))/Veff;
            particles[1].cartPos.x = (a0 * 1.4141* (Math.sqrt(Veff+1.3169)+1.1476))/Veff;
            //particles[2].cartPos.x = (-a0 * 1.4141* (Math.sqrt(Veff+1.3169)-1.1476))/Veff;
            //particles[1].cartPos.x = (a0 * 1.4141* (Math.sqrt(Veff+1.3169)-1.1476))/Veff;

        }

        else if(atom === "HE Cation"){
            particles[1].cartPos.x = (-a0 * (Math.sqrt(Veff+4)+2))/Veff;
            //particles[1].cartPos.x = (a0 * (Math.sqrt(Veff+1)-1))/Veff;
        }

        else if(atom === "HE Atom"){
            particles[2].cartPos.x = (-a0 * 1.4141* (Math.sqrt(Veff+5.8068)+2.4097))/Veff;
            particles[1].cartPos.x = (a0 * 1.4141* (Math.sqrt(Veff+5.8068)+2.4097))/Veff;
            //particles[2].cartPos.x = (-a0 * 1.4141* (Math.sqrt(Veff+5.8068)-2.4097))/Veff;
            //particles[1].cartPos.x = (a0 * 1.4141* (Math.sqrt(Veff+5.8068)-2.4097))/Veff;
        }

        else{}

    }
});
*/



/**
 * Toggles the physical "ansatz" of the model. The "ansatz" is the
 * mathematical assumption made in the model about the shape of the
 * atomic orbitals. The two possible "ansatze" are the 1s and 2p
 * orbitals.
 * @function changeAnsatz
 */
function changeAnsatz() {
    isAngularMomentum = !isAngularMomentum;
}

/**
 * Increases the energy-level of the electron in the H atom by one.
 * If the angular momentum is "on", it also increases the angular momentum
 * of the electron by one.
 * @function sp1
 */
function incrementEnergyLevel() {
    if(atom === "H Atom"){  
    if(isAngularMomentum){     
        // increase the distance of the electron from the proton by a factor of 4
        particles[1].cartPos.multiplyScalar(4);       
        // calculate a vector perpendicular to the connection vector between the electron and the proton
        const connectionVector = particles[0].cartPos.clone().sub(particles[1].cartPos);
        const perpendicularVector = new THREE.Vector3(-connectionVector.y, connectionVector.x, 0);   
        
        // normalize the perpendicular vector and multiply it by 1e-12 to get the correct units
        perpendicularVector.normalize().multiplyScalar(1e-12);
        
        // add the perpendicular vector to the velocity of the electron
        particles[1].velocity.add(perpendicularVector);
    }
    // increase the energy level (n) of the electron by one
    N++;
    // increase the angular momentum (l) of the electron by one
    L++;
    }

    /*
    else if(atom === "H2 Cation"){   
    // increase the velocity of the electron in the y direction by 1e-12
    particles[2].velocity.y =  1e-12;
    }

    else if(atom === "H Anion"){  
    // increase the velocity of the electron in the y direction by 1e-12
    particles[1].velocity.y =  1e-12;
    // decrease the velocity of the electron in the y direction by 1e-12
    particles[2].velocity.y =  -1e-12;
    }

    else if(atom === "HE Cation"){    
    // increase the velocity of the electron in the y direction by 1e-12
    particles[1].velocity.y = 1e-12;   
    }

    else if(atom === "HE Atom"){ 
    // increase the velocity of the electron in the y direction by 1e-12
    particles[1].velocity.y =  1e-12;
    // decrease the velocity of the electron in the y direction by 1e-12
    particles[2].velocity.y =  -1e-12;
    }

    else{}
    */

    // change the background color of the scene depending on the energy level
    if (N == 1) {
    radiation = 1;
    } else if (N == 2) { 
    radiation= 0.999905
    cameraEl.setAttribute('position', {x: 0, y: 0, z: 1200});
    } else if (N == 3) { 
    radiation = 0.99995
    cameraEl.setAttribute('position', {x: 0, y: 0, z: 2000});
    }

    // update the energy level and angular momentum in the UI
    document.getElementById("nLvl").innerHTML++;
    document.getElementById("lLvl").innerHTML++;
}

/**
 * Lower the energy level
 * Decrease the energy level by one if it is greater than 1
 * Decrease the angular momentum if the energy level is also decreased
 * Update the UI to reflect the new energy level and angular momentum
 */
function decrementEnergyLevel() {
    // only allow the energy level to be decreased if it is greater than 1
    if (parseInt(document.getElementById("nLvl").innerHTML) > 1) {
    // decrease the energy level
    N--;
    // if the energy level is decreased, also decrease the angular momentum if it is equal to the energy level
    if (L == N) {
        L--;  
        document.getElementById("lLvl").innerHTML--;
    }

    // update the UI to reflect the new energy level and angular momentum
    document.getElementById("nLvl").innerHTML = parseInt(document.getElementById("nLvl").innerHTML) - 1;
    }

    // update the radiation level based on the energy level
    if ( N == 1 ) { 
    radiation = 1
    cameraEl.setAttribute('position', {x: 0, y: 0, z: 800});
    } else if ( N == 2 ) { 
    radiation= 0.999905
    cameraEl.setAttribute('position', {x: 0, y: 0, z: 1200});
    } else if ( N == 3 ) { 
    radiation = 0.999999
    cameraEl.setAttribute('position', {x: 0, y: 0, z: 2000});
    }
}



/**
 * Increase the angular momentum
 * Increase the angular momentum by one if it is less than the energy level
 * Update the UI to reflect the new angular momentum
 */
function incrementAngularMomentum(){
    // only allow the angular momentum to be increased if it is less than the energy level
    if (L < N - 1) {
    // increase the angular momentum
    L++;
    // update the UI to reflect the new angular momentum
    document.getElementById("lLvl").innerHTML++;
    } else {
    // throw an error if the angular momentum is already at the maximum
    throw new Error("Invalid spin");
    }
}

/**
 * Decreases the angular momentum (L) by one if it is greater than zero.
 * Updates the UI to reflect the new angular momentum.
 * Throws an error if the angular momentum is already zero.
 * @function lm1
 */
function decrementAngularMomentum() {
    // check if angular momentum is greater than zero
    if (L > 0) {
    // decrease angular momentum
    L--;
    // update the UI to reflect the new angular momentum
    document.getElementById("lLvl").innerHTML--;
    } else {
    // throw an error if angular momentum is already zero
    throw new Error("Invalid spin");
    }
}



/**
 * Sets the velocity of all particles to zero.
 * @function zeroV
 */
function setVelocitiesToZero(){
    // iterate over all particles and set their velocity to zero
    particles.forEach(particle => particle.velocity.set(0, 0, 0));

    switch (atom) {
    case "H2 Cation":
        particles[2].velocity.x = 1.14e-13;
        break;
    case "H Anion":
        particles[1].velocity.x = -6.35e-14;
        particles[2].velocity.x = 6.35e-14;
        break;
    case "HE Atom":
        particles[1].velocity.x = 7.35e-13;
        particles[2].velocity.x = -7.35e-13;
        break;
    default:
        throw new Error(`Invalid atom: ${atom}`);
    }
}



/**
 * Changes between locked nucleus and movable nucleus
 * Toggles the `isPlocked` variable to switch between a locked nucleus and a movable nucleus
 * @function lockNucleus
 */
function lockNucleus(){
    // toggle the isPlocked variable
    isNucleusLocked = !isNucleusLocked;
}



/**
 * Toggles the path recording on or off.
 * @function startRecord
 */
function startRecord(){
    // toggle the RecordPath variable
    isRecordPath = !isRecordPath;

    // change the button text depending on the state
    if (isRecordPath) {
    pathButton.textContent = "⏸";
    } else {
    pathButton.textContent = "⏵";
    }
}



/**
 * Deletes all recorded paths from the scene.
 * This function is called when the user clicks the "Delete Path" button.
 * @function deletePath
 */
function deletePath(){
    // get all recorded paths from the scene
    const circles = document.querySelectorAll('a-circle');

    // iterate over all paths and remove them from the scene
    circles.forEach(circle => circle.parentNode.removeChild(circle));
}



// calculation of the model
function animation() {
    // calculate and display distance and energy
    const calculateDistanceAndEnergy = (particles, atom) => {
    switch (atom) {
        case "H Atom":
        const length = particles[1].cartPos.distanceTo(particles[0].cartPos);
        const relLen = length / a0;
        document.getElementById("Distance").innerHTML = relLen.toFixed(3) + ' a0';
        const energy = (a0**2 / length**2) - (2 * a0 / length);
        document.getElementById("Energy").innerHTML = energy.toFixed(3) + ' Ry';
        break;
        case "H2 Cation":
        const length1 = particles[2].cartPos.distanceTo(particles[1].cartPos);
        const length2 = particles[2].cartPos.distanceTo(particles[0].cartPos);
        const relLen2 = length1 / a0;
        document.getElementById("Distance").innerHTML = relLen2.toFixed(3) + ' a0';
        const energy2 = (a0**2 / length1**2) - (2 * a0 / length1) + (a0**2 / length2**2) - (2 * a0 / length2) + (2 * a0 * 0.79473 / (length1 + length2));
        document.getElementById("Energy").innerHTML = energy2.toFixed(3) + ' Ry';
        break;
        case "H Anion":
        const length3 = particles[1].cartPos.distanceTo(particles[0].cartPos);
        const length4 = particles[2].cartPos.distanceTo(particles[0].cartPos);
        const relLen4 = length3 / a0;
        document.getElementById("Distance").innerHTML = relLen4.toFixed(3) + ' a0';
        const energy4 = (a0**2 / length3**2) - (2 * a0 / length3) + (a0**2 / length4**2) - (2 * a0 / length4) + (2 * a0 * 0.72644 / (length3 + length4));
        document.getElementById("Energy").innerHTML = energy4.toFixed(3) + ' Ry';
        break;
        case "HE Cation":
        const length5 = particles[1].cartPos.distanceTo(particles[0].cartPos);
        const relLen5 = length5 / a0;
        document.getElementById("Distance").innerHTML = relLen5.toFixed(3) + ' a0';
        const energy5 = (a0**2 / length5**2) - (2 * 2 * a0 / length5);
        document.getElementById("Energy").innerHTML = energy5.toFixed(3) + ' Ry';
        break;
        case "HE Atom":
        const length6 = particles[1].cartPos.distanceTo(particles[0].cartPos);
        const length7 = particles[2].cartPos.distanceTo(particles[0].cartPos);
        const relLen7 = length6 / a0;
        document.getElementById("Distance").innerHTML = relLen7.toFixed(3) + ' a0';
        const energy7 = (2 * a0**2 / length6**2) - (4 * a0 * 1.95393 / length6) + (a0 / length6);
        document.getElementById("Energy").innerHTML = energy7.toFixed(3) + ' Ry';
        break;
        default:
        throw new Error("Invalid atom");
    }
    }

    calculateDistanceAndEnergy(particles, atom);

    // calculation of forces and resulting movements
    if (isAnimationRunning) {
    for (let i in particles) {
        for (let j in particles) {
        if (i !== j) {
            // different vector between two particles
            const distanceVector = new THREE.Vector3().copy(particles[j].cartPos).sub(particles[i].cartPos);
            const distancePolar = cartInPolar(distanceVector.x, distanceVector.y, distanceVector.z); 
            const perpendicularVector = new THREE.Vector3(-distanceVector.y, distanceVector.x, 0);
            
            // electric Force
            let electricForce = ( -Z * particles[i].charge * particles[j].charge) 
                                / (k * distancePolar[0] * distancePolar[0]); 
            let kratzerForce = 0;
            let spinForce = 0;
            let Frad = 0;

            // "kratzer" or repulsive Force
            if(particles[i].charge * particles[j].charge < 0){ 
            kratzerForce = (particles[i].charge * particles[j].charge  * a0) / (k * distancePolar[0] * distancePolar[0] *distancePolar[0]);
            if(!isAngularMomentum){
            spinForce = (particles[i].charge * particles[j].charge * (L * ( L+ 1))* a0) / (k * distancePolar[0] * distancePolar[0] *distancePolar[0]);
            }
            } else {}

            /*
            if(particles[i].charge * particles[j].charge > 0){
            if(atom==="HE Atom"){
                kratzerForce = (particles[i].charge * particles[j].charge * a0 * (-0.836709)) / (k * distancePolar[0] * distancePolar[0] *distancePolar[0]) ;
                electricForce = electricForce * (0.836709);
            }

            else if(atom === "H2 Cation"){
            electricForce = 0;
            }
            
            else{}
            }else{}*/

            // radiation Force
            if (particles[i].cartPos.x < particles[j].cartPos.x) {
            Frad = - RadConst * particles[i].charge * particles[j].charge * particles[i].charge * particles[j].charge * particles[i].velocity.x /(k * distancePolar[0] * distancePolar[0] *distancePolar[0] * particles[i].mass)
            } else if (particles[i].cartPos.x > particles[j].cartPos.x) {
            Frad = RadConst * particles[i].charge * particles[j].charge * particles[i].charge * particles[j].charge * particles[i].velocity.x /(k * distancePolar[0] * distancePolar[0] *distancePolar[0] * particles[i].mass)
            } else{}

            // acceleration of the particle
            const particleAcceleration = distanceVector.clone().normalize().multiplyScalar((electricForce+kratzerForce+(Frad*1e0))/particles[i].mass);
            const spinAcceleration = perpendicularVector.clone().normalize().multiplyScalar(spinForce/(particles[i].mass));

            // velocity of the particle
            particles[i].velocity.add(particleAcceleration.multiplyScalar(1e-35)).add(spinAcceleration.multiplyScalar(1e-37));

            // position of the particle
            particles[i].cartPos.add(particles[i].velocity.multiplyScalar(radiation));
            
            // duration of one circulation on first excited state
            if (particles[i].charge < 0 ){
            const force = particleAcceleration.clone().multiplyScalar(particles[i].mass);
            document.getElementById("fx").innerHTML = force.x;
            document.getElementById("fy").innerHTML = force.y;
            document.getElementById("fz").innerHTML = force.z;

            if (countTime == 1){
                if ( angle.toFixed(2) == distancePolar[2].toFixed(2) ){
                console.log('ja')
                countTime = 2
                }

                timePassed = timePassed + 2 * 1.6 / 885 ;
                console.log(timePassed + " ps")
            }

            if(distanceVector.length() >= 4 * a0){
                if(countTime == 0){
                countTime = 1
                angle = distancePolar[2]
                }
            }
            }    
        }
        }
    }
    } else{}
}



// interval for repeating animation
let interval = setInterval(animation, dt);   

/**
 * Changes interval for repeating animation
 * @param {number} newDt interval time in ms
 */
function changeInterval(newDt) {
    // clear the old interval
    clearInterval(interval);
    // set the new interval
    interval = setInterval(animation, Number(newDt));
    // update the global variable dt
    dt = Number(newDt);
    //console.log(dt + typeof(dt))
}



/**
 * Converts real coordinates into "image" coordinates and creates animation
 * @param {int} frames - number of frames passed
 */
function intoRealWorld() {
    frames++;
    let a = 0;

    if(particles.length > 1){
    for (let i in particles) {      
        const pr = particles[i].cartPos.clone().multiplyScalar(1e12);

        // set positions of particles in the scene
        particleEl[i].object3D.position.set(pr.x, pr.y, pr.z);

        // set positions of protons and electrons in the scene
        if (particles[i].type === 'proton') {
        bohrRadius[a].object3D.position.set(pr.x, pr.y, pr.z);
        a++;
        } else if (particles[i].type === "electron" && i == 1) {
        document.getElementById("posx1").innerHTML = pr.x.toFixed(3);
        document.getElementById("y1").innerHTML = pr.y.toFixed(3);
        document.getElementById("z1").innerHTML = pr.z.toFixed(3);
        } else if (particles[i].type === "electron" && i ==2 ) {
        document.getElementById("x2").innerHTML = pr.x.toFixed(3);
        document.getElementById("y2").innerHTML = pr.y.toFixed(3);
        document.getElementById("z2").innerHTML = pr.z.toFixed(3);
        } else if (particles[i].type === "proton" && i == 1) {
        document.getElementById("posx1").innerHTML = NaN;
        document.getElementById("y1").innerHTML = NaN;
        document.getElementById("z1").innerHTML = NaN;
        } else if(particles[i].type === "proton" && i == 2 || particles.length < 3) {
        document.getElementById("x2").innerHTML = NaN;
        document.getElementById("y2").innerHTML = NaN;
        document.getElementById("z2").innerHTML = NaN;
        } else {}

        // record path of particles if requested
        if(frames % 5 == 0 & isRecordPath == true){
        if(particles[i].type === 'electron'){
            const point = document.createElement('a-circle'); 
            point.setAttribute('radius', 1);
            point.setAttribute('color', 'grey');
            point.object3D.position.set(pr.x, pr.y, pr.z);
            sceneEl.appendChild(point); 
        }
        }
    }
    }
    // call the function again after the frame is rendered
    requestAnimationFrame(intoRealWorld);
}

intoRealWorld();    