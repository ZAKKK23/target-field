/* ===========================================================
   plane-vtp-target-field — page glue
   =========================================================== */

/* ---------- tabs ---------- */
const tabButtons = document.querySelectorAll("nav.tabs button");
const views = document.querySelectorAll(".view");
function showView(name) {
  tabButtons.forEach((b) => b.classList.toggle("active", b.dataset.view === name));
  views.forEach((v) => v.classList.toggle("active", v.id === "view-" + name));
}
tabButtons.forEach((b) => b.addEventListener("click", () => showView(b.dataset.view)));
window.addEventListener("hashchange", () => {
  const h = location.hash.replace("#", "");
  if (h) showView(h);
});
if (location.hash) showView(location.hash.replace("#", ""));

/* ---------- simulation wiring ---------- */
const canvas = document.getElementById("simCanvas");
const sim = new VTPSim(canvas, PARAMS.N_DEFAULT);

const $ = (id) => document.getElementById(id);
const statusLine = $("statusLine");
const numAgents = $("numAgents");
const spfSlider = $("spfSlider"), spfVal = $("spfVal");
const followCb = $("followCb"), headingsCb = $("headingsCb");
const nuSlider = $("nuSlider"), nuVal = $("nuVal");
const spacingSlider = $("spacingSlider"), spacingVal = $("spacingVal");
const halfSlider = $("halfSlider"), halfVal = $("halfVal");
const staticCount = $("staticCount");
const freqSlider = $("freqSlider"), freqNum = $("freqNum"), freqInfo = $("freqInfo");
const ampSlider = $("ampSlider"), ampVal = $("ampVal");
const angSlider = $("angSlider"), angVal = $("angVal");
const cxSlider = $("cxSlider"), cxVal = $("cxVal");
const cySlider = $("cySlider"), cyVal = $("cyVal");
const btnPause = $("btnPause"), btnReset = $("btnReset"), btnRandom = $("btnRandom");

function showFreq() {
  const f = sim.oscFreq;
  freqInfo.textContent = f > 0 ? `cycles/step, period ${(1 / f).toFixed(f >= 0.1 ? 1 : 0)} steps` : "cycles/step (frozen)";
}
function showStaticCount() {
  staticCount.textContent = `${sim.nStatic} static targets`;
}
function setFreq(f) {
  f = clamp(Number.isFinite(f) ? f : 0, 0, 1);
  sim.oscFreq = f;                       // phase accumulates, so changing f never makes the target jump
  freqSlider.value = f;
  if (document.activeElement !== freqNum) freqNum.value = String(f);
  showFreq();
}
function rebuildLattice() {
  sim.buildLattice();
  showStaticCount();
}

spfSlider.addEventListener("input", () => { sim.stepsPerFrame = +spfSlider.value; spfVal.textContent = spfSlider.value; });
followCb.addEventListener("change", () => { sim.followSwarm = followCb.checked; });
headingsCb.addEventListener("change", () => { sim.showHeadings = headingsCb.checked; });
nuSlider.addEventListener("input", () => { sim.nu = +nuSlider.value; nuVal.textContent = sim.nu.toFixed(2); });
spacingSlider.addEventListener("input", () => {
  sim.gridSpacing = +spacingSlider.value; spacingVal.textContent = sim.gridSpacing.toFixed(1); rebuildLattice();
});
halfSlider.addEventListener("input", () => {
  sim.gridHalf = +halfSlider.value; halfVal.textContent = String(sim.gridHalf); rebuildLattice();
});
freqSlider.addEventListener("input", () => setFreq(+freqSlider.value));
freqNum.addEventListener("input", () => { const v = parseFloat(freqNum.value); if (!Number.isNaN(v)) setFreq(v); });
freqNum.addEventListener("change", () => { freqNum.value = String(sim.oscFreq); });
ampSlider.addEventListener("input", () => { sim.oscAmp = +ampSlider.value; ampVal.textContent = sim.oscAmp.toFixed(1); });
angSlider.addEventListener("input", () => { sim.oscAngleDeg = +angSlider.value; angVal.textContent = sim.oscAngleDeg + "\u00b0"; });
cxSlider.addEventListener("input", () => { sim.oscCx = +cxSlider.value; cxVal.textContent = sim.oscCx.toFixed(1); });
cySlider.addEventListener("input", () => { sim.oscCy = +cySlider.value; cyVal.textContent = sim.oscCy.toFixed(1); });

btnPause.addEventListener("click", () => {
  sim.paused = !sim.paused;
  btnPause.textContent = sim.paused ? "Resume" : "Pause";
});

function restoreDefaults() {
  sim.nu = PARAMS.NU_DEFAULT;               nuSlider.value = sim.nu;        nuVal.textContent = sim.nu.toFixed(2);
  sim.gridSpacing = PARAMS.GRID_SPACING;    spacingSlider.value = sim.gridSpacing; spacingVal.textContent = sim.gridSpacing.toFixed(1);
  sim.gridHalf = PARAMS.GRID_HALF;          halfSlider.value = sim.gridHalf; halfVal.textContent = String(sim.gridHalf);
  sim.oscAmp = PARAMS.OSC_AMP;              ampSlider.value = sim.oscAmp;   ampVal.textContent = sim.oscAmp.toFixed(1);
  sim.oscAngleDeg = PARAMS.OSC_ANGLE_DEG;   angSlider.value = sim.oscAngleDeg; angVal.textContent = sim.oscAngleDeg + "\u00b0";
  sim.oscCx = PARAMS.OSC_CX;                cxSlider.value = sim.oscCx;     cxVal.textContent = sim.oscCx.toFixed(1);
  sim.oscCy = PARAMS.OSC_CY;                cySlider.value = sim.oscCy;     cyVal.textContent = sim.oscCy.toFixed(1);
  setFreq(PARAMS.OSC_FREQ);
  rebuildLattice();
}
btnReset.addEventListener("click", () => {
  restoreDefaults();
  sim.setup(+numAgents.value, true);
});
btnRandom.addEventListener("click", () => sim.setup(+numAgents.value, false));
numAgents.addEventListener("change", () => sim.setup(+numAgents.value, true));

showFreq();
showStaticCount();

/* ---------- animation loop ---------- */
function loop() {
  if (!sim.paused) for (let k = 0; k < sim.stepsPerFrame; k++) sim.step();
  sim.draw();
  const f = sim.oscFreq;
  statusLine.textContent =
    `t = ${sim.t} | N = ${sim.N} | ${sim.nStatic} static + 1 oscillating target | ` +
    `f = ${f.toFixed(4)} cycles/step${f > 0 ? ` (period ${(1 / f).toFixed(1)})` : ""} | \u03bd = ${sim.nu.toFixed(2)}`;
  requestAnimationFrame(loop);
}
requestAnimationFrame(loop);

/* ---------- matlab file viewer ---------- */
const FILE_NOTES = {
  "dynamics.m": "Main script. Target-field settings (lattice and oscillation frequency) are at the top; then the VTP loop: Delaunay neighbors, repulsion + alignment + homing, Voronoi speed governor.",
  "Target.m": "Target class: target geometry and the homeToTarget homing-vector method. The lattice and the oscillating target are point components of one Target object.",
  "nearestOnSegment.m": "Geometry helper used internally by Target.m for segment-shaped targets.",
  "neighborhoods.m": "Delaunay-graph neighbors of every agent, plus each agent's nearest neighbor and its distance.",
  "transition.m": "Smooth cutoff functions (the model uses 'expReciprocal'): 1 at 0, exactly 0 from 1 onward.",
  "alignTo.m": "Alignment term: Delaunay neighbors' headings weighted by how well they agree with the agent's own.",
  "voronoiProjectToBoundary.m": "Distance from each agent to its Voronoi cell boundary along its direction of motion (the speed governor).",
  "voronoiForwardArea.m": "Forward Voronoi area, used only when fdim = 2 in dynamics.m.",
  "poly_area.m": "Polygon area helper used by voronoiForwardArea.m.",
};
const FILE_ORDER = ["dynamics.m", "Target.m", "nearestOnSegment.m", "neighborhoods.m", "transition.m", "alignTo.m", "voronoiProjectToBoundary.m", "voronoiForwardArea.m", "poly_area.m"];

const fileList = $("fileList");
const codeView = $("codeView");
const fileHint = $("fileHint");
let MATLAB_SRC = null;

async function loadMatlabSource() {
  try {
    const res = await fetch("assets/matlab_src.json");
    MATLAB_SRC = await res.json();
  } catch (e) {
    MATLAB_SRC = null;
    codeView.textContent = "Could not load assets/matlab_src.json (serve the folder over http, e.g. python3 -m http.server).";
    return;
  }
  let first = true;
  FILE_ORDER.forEach((name) => {
    if (!(name in MATLAB_SRC)) return;
    const btn = document.createElement("button");
    btn.textContent = name;
    btn.addEventListener("click", () => selectFile(name));
    fileList.appendChild(btn);
    if (first) { selectFile(name); first = false; }
  });
}
function selectFile(name) {
  [...fileList.children].forEach((b) => b.classList.toggle("active", b.textContent === name));
  codeView.textContent = MATLAB_SRC[name];
  fileHint.textContent = FILE_NOTES[name] || "";
}
loadMatlabSource();
