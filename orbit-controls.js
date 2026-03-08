/**
 * orbit-controls – CAD-style camera component for A-Frame
 *
 * Controls:
 *   Mouse wheel              → zoom in / out
 *   Middle mouse + drag      → orbit around target
 *   Right mouse  + drag      → pan target point
 *
 * Rotation fix:
 *   obj.lookAt() on a THREE.Group uses the non-camera convention
 *   (makes +Z point toward target, so camera looks AWAY).
 *   Instead we compute Euler angles directly from spherical coords:
 *     pitch (X) = phi   - PI/2
 *     yaw   (Y) = theta
 *   and call obj.rotation.set(pitch, yaw, 0, 'YXZ').
 */
AFRAME.registerComponent('orbit-controls', {
  schema: {
    target:       { type: 'vec3',   default: { x: 0, y: 0, z: 0 } },
    minDistance:  { type: 'number', default: 10 },
    maxDistance:  { type: 'number', default: 15000 },
    zoomSpeed:    { type: 'number', default: 0.12 },
    rotateSpeed:  { type: 'number', default: 0.005 },
    panSpeed:     { type: 'number', default: 0.002 },
  },

  init() {
    this.target   = new THREE.Vector3(
      this.data.target.x, this.data.target.y, this.data.target.z
    );

    // Start at sensible defaults; will be refined on renderstart
    this.distance = 200;
    this.theta    = 0;           // azimuth angle around Y axis
    this.phi      = Math.PI / 2; // polar angle from Y axis (π/2 = equatorial)

    this._mode  = null; // 'orbit' | 'pan'
    this._lastX = 0;
    this._lastY = 0;

    // Bind handlers
    this._onWheel       = this._onWheel.bind(this);
    this._onPointerDown = this._onPointerDown.bind(this);
    this._onPointerMove = this._onPointerMove.bind(this);
    this._onPointerUp   = this._onPointerUp.bind(this);
    this._onContext     = e => e.preventDefault();

    // Defer canvas attachment until the renderer has started
    // (canvas can be null during component init in embedded scenes)
    this._attachListeners = this._attachListeners.bind(this);
    if (this.el.sceneEl.canvas) {
      this._attachListeners();
    } else {
      this.el.sceneEl.addEventListener('renderstart', this._attachListeners, { once: true });
    }

    // Once scene is ready, read initial position to derive spherical coords
    const ready = () => {
      const p = this.el.object3D.position;
      const dx = p.x - this.target.x;
      const dy = p.y - this.target.y;
      const dz = p.z - this.target.z;
      const d  = Math.sqrt(dx * dx + dy * dy + dz * dz);
      if (d > 0) {
        this.distance = d;
        this.phi      = Math.acos(Math.max(-1, Math.min(1, dy / d)));
        this.theta    = Math.atan2(dx, dz);
      }
    };
    if (this.el.sceneEl.hasLoaded) {
      ready();
    } else {
      this.el.sceneEl.addEventListener('loaded', ready, { once: true });
    }
  },

  _attachListeners() {
    const canvas = this.el.sceneEl.canvas;
    canvas.addEventListener('wheel',        this._onWheel,       { passive: false });
    canvas.addEventListener('pointerdown',  this._onPointerDown);
    canvas.addEventListener('pointermove',  this._onPointerMove);
    canvas.addEventListener('pointerup',    this._onPointerUp);
    canvas.addEventListener('pointerleave', this._onPointerUp);
    canvas.addEventListener('contextmenu',  this._onContext);
    this._canvas = canvas;
  },

  remove() {
    const canvas = this._canvas;
    if (!canvas) return;
    canvas.removeEventListener('wheel',        this._onWheel);
    canvas.removeEventListener('pointerdown',  this._onPointerDown);
    canvas.removeEventListener('pointermove',  this._onPointerMove);
    canvas.removeEventListener('pointerup',    this._onPointerUp);
    canvas.removeEventListener('pointerleave', this._onPointerUp);
    canvas.removeEventListener('contextmenu',  this._onContext);
  },

  _onWheel(e) {
    e.preventDefault();
    const dir    = e.deltaY > 0 ? 1 : -1;
    const factor = 1 + dir * this.data.zoomSpeed * 3;
    this.distance = Math.max(this.data.minDistance,
                    Math.min(this.data.maxDistance, this.distance * factor));
  },

  _onPointerDown(e) {
    if (e.button === 1) {         // middle mouse → orbit
      this._mode = 'orbit';
      e.preventDefault();         // suppress browser autoscroll cursor
    } else if (e.button === 2) {  // right mouse → pan
      this._mode = 'pan';
    } else {
      return;
    }
    this._lastX = e.clientX;
    this._lastY = e.clientY;
    try { this._canvas.setPointerCapture(e.pointerId); } catch (_) {}
  },

  _onPointerMove(e) {
    if (!this._mode) return;
    const dx = e.clientX - this._lastX;
    const dy = e.clientY - this._lastY;
    this._lastX = e.clientX;
    this._lastY = e.clientY;

    if (this._mode === 'orbit') {
      this.theta -= dx * this.data.rotateSpeed;
      this.phi    = Math.max(0.05,
                    Math.min(Math.PI - 0.05,
                    this.phi - dy * this.data.rotateSpeed));
    } else if (this._mode === 'pan') {
      const scale = this.distance * this.data.panSpeed;
      const obj   = this.el.object3D;
      const right = new THREE.Vector3(1, 0, 0).applyEuler(obj.rotation);
      const up    = new THREE.Vector3(0, 1, 0).applyEuler(obj.rotation);
      this.target.addScaledVector(right, -dx * scale);
      this.target.addScaledVector(up,     dy * scale);
    }
  },

  _onPointerUp(e) {
    if (!this._mode) return;
    this._mode = null;
    try { this._canvas.releasePointerCapture(e.pointerId); } catch (_) {}
  },


  tick() {
    const { distance, theta, phi } = this;
    const sinPhi = Math.sin(phi);

    // Spherical → Cartesian camera position
    const x = this.target.x + distance * sinPhi * Math.sin(theta);
    const y = this.target.y + distance * Math.cos(phi);
    const z = this.target.z + distance * sinPhi * Math.cos(theta);

    const obj = this.el.object3D;
    obj.position.set(x, y, z);

    // Set rotation via Euler angles (YXZ order) derived from spherical coords.
    // obj.lookAt() must NOT be used here: on a THREE.Group (non-camera),
    // THREE.js uses the reversed convention, making the camera look away from target.
    //   pitch (X rotation) = phi - PI/2
    //   yaw   (Y rotation) = theta
    obj.rotation.set(phi - Math.PI / 2, theta, 0, 'YXZ');
  }
});
