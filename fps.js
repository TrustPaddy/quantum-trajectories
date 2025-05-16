var frames = 0;
var fps;

/**
 * A component to generate data for the FPS counter.
 * Calculates the average number of frames rendered per second.
 * Emits the "event1" event with the FPS value every second.
 */
AFRAME.registerComponent('generate-data', {
    /**
     * Initializes the component.
     * Sets up the timer and binding for the FPS event.
     */
    init() {
    this.startTime = Date.now();
    this.emitFPS = this.emitFPS.bind(this);

    // Set up the timer to emit the FPS value every second
    setInterval(this.emitFPS, 1000);
    },
    /**
     * Calculates the average number of frames rendered per second.
     * Emits the "event1" event with the FPS value.
     */
    emitFPS() {
    const currentTime = Date.now();
    const elapsedTime = (currentTime - this.startTime) / 1000;
    fps = frames / elapsedTime;

    // Emit the FPS value as an event
    this.el.emit("event1", { fps: fps.toFixed(2) });

    // Reset the frame counter
    frames = 0;
    this.startTime = currentTime;
    },
});