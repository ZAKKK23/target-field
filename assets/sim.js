/* ===========================================================
   VTP — Target Field  |  live simulation
   Ported from the MATLAB model (matlab/dynamics.m). This is the
   Voronoi Topological Perception (VTP) model on the plane:
     - Neighbors are TOPOLOGICAL: an agent's neighbors are its
       Delaunay-triangulation neighbors (computed here with a small
       built-in Bowyer-Watson triangulation — no external libraries).
     - Repulsion from the nearest neighbor, fading out with distance
       (expReciprocal transition, length scale L).
     - Alignment with Delaunay neighbors, weighted by how well their
       headings already agree (alignTo.m).
     - Homing toward the nearest target, with weight (1 - s) so that
       repulsion takes priority when an agent is crowded.
     - Speed governor: each agent moves at tanh(l/L) times its desired
       velocity, where l is the distance from the agent to its Voronoi
       cell boundary along the direction of motion (so crowded agents
       slow down, agents on the hull move at full speed).
   TARGET FIELD (the scenario this repo is about):
     - a square lattice of STATIC point targets across the plane
     - ONE oscillating target on a line through osc_center, whose
       frequency (cycles per time step) you choose live.
   Every agent homes toward whichever target (static or oscillating)
   is nearest to it.
   =========================================================== */

const PARAMS = {
  N_DEFAULT: 400,
  L: 1,                 // interaction length scale (MATLAB: L = 1); also M0 for the speed governor
  NU_DEFAULT: 2.5,      // alignment strength (live-editable: sim.nu)
  GRID_HALF: 20,        // static lattice covers [-GRID_HALF, GRID_HALF]^2
  GRID_SPACING: 5,      // one static target every GRID_SPACING units
  OSC_FREQ: 0.01,       // oscillating target frequency, cycles per time step
  OSC_AMP: 8,           // amplitude of the swing
  OSC_ANGLE_DEG: 90,    // direction of the swing (0 = horizontal, 90 = vertical)
  OSC_CX: 0,            // center of the swing
  OSC_CY: 0,
  POS_SEED: 2,          // same seeds the MATLAB script uses (rng(2), rng(18))
  ANG_SEED: 18,
};

function clamp(x, a, b) { return Math.max(a, Math.min(b, x)); }

// small seeded PRNG (mulberry32) so a Reset reproduces the same initial condition
function mulberry32(seed) {
  let a = (seed >>> 0) + 0x6D2B79F5;
  return function () {
    a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/* transition.m, 'expReciprocal': smooth step, 1 at x=0 and exactly 0 for x>=1 */
function expReciprocal(x) {
  if (x >= 1) return 0;
  const a = Math.exp(-1 / (1 - x));
  const b = Math.exp(-1 / x);        // x = 0 -> exp(-Infinity) = 0 -> result 1
  return a / (b + a);
}

/* ---------- Delaunay triangulation (Bowyer-Watson) ----------
   Returns nbr[i] = array of Delaunay neighbors of point i.
   X is a flat [x0,y0,x1,y1,...] array of N points. */
function delaunayNeighbors(X, N) {
  const nbr = new Array(N);
  for (let i = 0; i < N; i++) nbr[i] = [];
  if (N < 2) return nbr;

  let minx = Infinity, maxx = -Infinity, miny = Infinity, maxy = -Infinity;
  for (let i = 0; i < N; i++) {
    const x = X[2 * i], y = X[2 * i + 1];
    if (x < minx) minx = x; if (x > maxx) maxx = x;
    if (y < miny) miny = y; if (y > maxy) maxy = y;
  }
  const ox = (minx + maxx) / 2, oy = (miny + maxy) / 2;
  const ext = Math.max(maxx - minx, maxy - miny, 1);
  const D = 60 * ext;                       // super-triangle size

  const P = new Float64Array(2 * (N + 3));
  for (let i = 0; i < N; i++) { P[2 * i] = X[2 * i] - ox; P[2 * i + 1] = X[2 * i + 1] - oy; }
  P[2 * N] = 0;       P[2 * N + 1] = 3 * D;
  P[2 * N + 2] = -3 * D; P[2 * N + 3] = -3 * D;
  P[2 * N + 4] = 3 * D;  P[2 * N + 5] = -3 * D;

  const makeTri = (a, b, c) => {
    const ax = P[2 * a], ay = P[2 * a + 1], bx = P[2 * b], by = P[2 * b + 1], cx = P[2 * c], cy = P[2 * c + 1];
    const d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    if (Math.abs(d) < 1e-300) return { a, b, c, x: 0, y: 0, r2: 0 };   // degenerate: never "contains" anything
    const a2 = ax * ax + ay * ay, b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
    const ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d;
    const uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d;
    const dx = ax - ux, dy = ay - uy;
    return { a, b, c, x: ux, y: uy, r2: dx * dx + dy * dy };
  };

  let tris = [makeTri(N, N + 1, N + 2)];
  const M = N + 3;
  for (let p = 0; p < N; p++) {
    const px = P[2 * p], py = P[2 * p + 1];
    const keep = [];
    const edgeCount = new Map();
    const addEdge = (u, v) => {
      const k = u < v ? u * M + v : v * M + u;
      const e = edgeCount.get(k);
      if (e) e.n++; else edgeCount.set(k, { u, v, n: 1 });
    };
    for (let t = 0; t < tris.length; t++) {
      const T = tris[t];
      const dx = px - T.x, dy = py - T.y;
      if (dx * dx + dy * dy < T.r2) { addEdge(T.a, T.b); addEdge(T.b, T.c); addEdge(T.c, T.a); }
      else keep.push(T);
    }
    for (const e of edgeCount.values()) if (e.n === 1) keep.push(makeTri(e.u, e.v, p));
    tris = keep;
  }

  const link = (i, j) => { if (nbr[i].indexOf(j) < 0) nbr[i].push(j); if (nbr[j].indexOf(i) < 0) nbr[j].push(i); };
  for (const T of tris) {
    if (T.a >= N || T.b >= N || T.c >= N) continue;   // drop triangles touching the super-triangle
    link(T.a, T.b); link(T.b, T.c); link(T.c, T.a);
  }
  return nbr;
}

class VTPSim {
  constructor(canvas, N = PARAMS.N_DEFAULT) {
    this.canvas = canvas;
    this.ctx = canvas ? canvas.getContext("2d") : null;

    // ---- live-editable parameters ----
    this.nu = PARAMS.NU_DEFAULT;
    this.L = PARAMS.L;
    this.gridHalf = PARAMS.GRID_HALF;
    this.gridSpacing = PARAMS.GRID_SPACING;
    this.oscFreq = PARAMS.OSC_FREQ;
    this.oscAmp = PARAMS.OSC_AMP;
    this.oscAngleDeg = PARAMS.OSC_ANGLE_DEG;
    this.oscCx = PARAMS.OSC_CX;
    this.oscCy = PARAMS.OSC_CY;
    this.followSwarm = true;       // MATLAB: fixframe = false
    this.showHeadings = true;
    this.stepsPerFrame = 1;
    this.paused = false;
    this.cam = null;
    this.onSetupChange = null;

    this.buildLattice();
    this.setup(N);
  }

  /* static lattice: meshgrid(-gridHalf:gridSpacing:gridHalf) */
  buildLattice() {
    const pts = [];
    const sp = Math.max(this.gridSpacing, 0.25);
    const n = Math.floor((2 * this.gridHalf) / sp + 1e-9);
    for (let iy = 0; iy <= n; iy++)
      for (let ix = 0; ix <= n; ix++)
        pts.push(-this.gridHalf + ix * sp, -this.gridHalf + iy * sp);
    this.statics = Float64Array.from(pts);
    this.nStatic = pts.length / 2;
  }

  oscPos() {
    const a = (this.oscAngleDeg * Math.PI) / 180;
    const s = this.oscAmp * Math.sin(this.phase);
    return [this.oscCx + s * Math.cos(a), this.oscCy + s * Math.sin(a)];
  }

  /* new run. seeded = true reproduces the MATLAB-style fixed seeds. */
  setup(N, seeded = true) {
    this.N = N = Math.max(3, Math.round(N));
    const rngP = mulberry32(seeded ? PARAMS.POS_SEED : (Math.random() * 1e9) | 0);
    const rngA = mulberry32(seeded ? PARAMS.ANG_SEED : (Math.random() * 1e9) | 0);
    const icRad = 0.5 * Math.sqrt((N * Math.PI) / 4 / 0.91);
    this.X = new Float64Array(2 * N);
    this.U = new Float64Array(2 * N);     // velocity from the previous step (alignment uses its direction)
    this.U1 = new Float64Array(2 * N);    // desired direction this step (what MATLAB plots)
    this.Unext = new Float64Array(2 * N); // velocity applied at the next commit
    for (let i = 0; i < N; i++) {
      this.X[2 * i] = icRad * (2 * rngP() - 1);
      this.X[2 * i + 1] = icRad * (2 * rngP() - 1);
    }
    for (let i = 0; i < N; i++) {
      const ang = 2 * Math.PI * rngA();
      this.U[2 * i] = Math.cos(ang);
      this.U[2 * i + 1] = Math.sin(ang);
    }
    this.t = 1;                                   // MATLAB loop index of the displayed state
    this.phase = 2 * Math.PI * this.oscFreq * this.t;
    this.cam = null;
    this.compute();
    if (this.onSetupChange) this.onSetupChange();
  }

  /* one MATLAB loop iteration up to (not including) the position update:
     moves the oscillating target, then computes U1 (desired direction)
     and Unext (= tanh(l/L) * U1) for the current positions. */
  compute() {
    const N = this.N, X = this.X, L = this.L, nu = this.nu;
    const [ox, oy] = this.oscPos();
    this.oscX = ox; this.oscY = oy;

    const nbr = delaunayNeighbors(X, N);
    this.nbr = nbr;

    // nearest Delaunay neighbor + distance (neighborhoods.m)
    const s = new Float64Array(N);
    const rx = new Float64Array(N), ry = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      let best = Infinity, bj = -1;
      for (const j of nbr[i]) {
        const dx = X[2 * i] - X[2 * j], dy = X[2 * i + 1] - X[2 * j + 1];
        const d2 = dx * dx + dy * dy;
        if (d2 < best) { best = d2; bj = j; }
      }
      if (bj < 0) { s[i] = 0; continue; }
      const d = Math.sqrt(best);
      s[i] = expReciprocal(d / L);
      if (d > 0) {                         // repulsion: length s, pointing away from the nearest neighbor
        rx[i] = s[i] * (X[2 * i] - X[2 * bj]) / d;
        ry[i] = s[i] * (X[2 * i + 1] - X[2 * bj + 1]) / d;
      }
    }

    // normalized previous headings for alignment
    const nx = new Float64Array(N), ny = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      const m = Math.hypot(this.U[2 * i], this.U[2 * i + 1]);
      if (m > 0) { nx[i] = this.U[2 * i] / m; ny[i] = this.U[2 * i + 1] / m; }
    }

    for (let i = 0; i < N; i++) {
      // alignment (alignTo.m): sum of w(<ui,uj>) uj over Delaunay neighbors, times 1/6
      let ax = 0, ay = 0;
      for (const j of nbr[i]) {
        let dp = nx[i] * nx[j] + ny[i] * ny[j];
        if (dp > 1) dp = 1; else if (dp < -1) dp = -1;
        const w = expReciprocal(Math.acos(dp) / Math.PI);
        ax += w * nx[j]; ay += w * ny[j];
      }
      ax /= 6; ay /= 6;

      // homing: unit vector toward the nearest target, length (1 - s)
      let bestD2 = Infinity, hx = 0, hy = 0;
      const px = X[2 * i], py = X[2 * i + 1];
      const S = this.statics;
      for (let k = 0; k < this.nStatic; k++) {
        const dx = S[2 * k] - px, dy = S[2 * k + 1] - py;
        const d2 = dx * dx + dy * dy;
        if (d2 < bestD2) { bestD2 = d2; hx = dx; hy = dy; }
      }
      {
        const dx = ox - px, dy = oy - py;
        const d2 = dx * dx + dy * dy;
        if (d2 < bestD2) { bestD2 = d2; hx = dx; hy = dy; }
      }
      const hn = Math.hypot(hx, hy);
      if (hn > 0) { hx = ((1 - s[i]) * hx) / hn; hy = ((1 - s[i]) * hy) / hn; } else { hx = 0; hy = 0; }

      // direction
      const u1x = (rx[i] + hx + nu * ax) / (1 + nu);
      const u1y = (ry[i] + hy + nu * ay) / (1 + nu);
      this.U1[2 * i] = u1x; this.U1[2 * i + 1] = u1y;

      // speed governor: distance to the Voronoi cell boundary along u1.
      // The cell is the intersection of the half-planes toward each Delaunay
      // neighbor, so the exit distance is min over neighbors j with u.d_j > 0
      // of |d_j|^2 / (2 u.d_j). No such neighbor -> unbounded cell -> Infinity.
      let m = Math.hypot(u1x, u1y), ux = 1, uy = 0;
      if (m > 0) { ux = u1x / m; uy = u1y / m; }
      let l = Infinity;
      for (const j of nbr[i]) {
        const dx = X[2 * j] - px, dy = X[2 * j + 1] - py;
        const den = ux * dx + uy * dy;
        if (den > 1e-12) {
          const tt = (dx * dx + dy * dy) / (2 * den);
          if (tt < l) l = tt;
        }
      }
      const sp = Math.tanh(l / L);          // tanh(Infinity) = 1
      this.Unext[2 * i] = sp * u1x; this.Unext[2 * i + 1] = sp * u1y;
    }
  }

  /* apply the position update, advance time and the oscillation phase, compute the next step */
  step() {
    const N = this.N;
    for (let i = 0; i < 2 * N; i++) { this.X[i] += this.Unext[i]; this.U[i] = this.Unext[i]; }
    this.t++;
    this.phase += 2 * Math.PI * this.oscFreq;
    this.compute();
  }

  updateCamera() {
    const N = this.N, X = this.X;
    let cx, cy, hw;
    if (this.followSwarm) {
      // MATLAB: frame = center of mass +/- 3 * median radius
      let mx = 0, my = 0;
      for (let i = 0; i < N; i++) { mx += X[2 * i]; my += X[2 * i + 1]; }
      mx /= N; my /= N;
      const r = new Array(N);
      for (let i = 0; i < N; i++) r[i] = Math.hypot(X[2 * i] - mx, X[2 * i + 1] - my);
      r.sort((a, b) => a - b);
      cx = mx; cy = my; hw = Math.max(3 * r[N >> 1], 4);
    } else {
      cx = 0; cy = 0;
      hw = 1.15 * Math.max(this.gridHalf, Math.hypot(this.oscCx, this.oscCy) + this.oscAmp) + 2;
    }
    if (!this.cam) this.cam = { cx, cy, hw };
    const a = 0.08;
    this.cam.cx += (cx - this.cam.cx) * a;
    this.cam.cy += (cy - this.cam.cy) * a;
    this.cam.hw += (hw - this.cam.hw) * a;
  }

  draw() {
    const canvas = this.canvas, ctx = this.ctx;
    const dpr = window.devicePixelRatio || 1;
    const w = canvas.clientWidth, h = canvas.clientHeight;
    if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
      canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
    }
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);

    this.updateCamera();
    const { cx, cy, hw } = this.cam;
    const scale = Math.min(w, h) / (2 * hw);
    const toPx = (x, y) => [w / 2 + (x - cx) * scale, h / 2 - (y - cy) * scale];

    // path of the oscillating target
    {
      const a = (this.oscAngleDeg * Math.PI) / 180;
      const [x0, y0] = toPx(this.oscCx - this.oscAmp * Math.cos(a), this.oscCy - this.oscAmp * Math.sin(a));
      const [x1, y1] = toPx(this.oscCx + this.oscAmp * Math.cos(a), this.oscCy + this.oscAmp * Math.sin(a));
      ctx.strokeStyle = "rgba(242,104,95,0.35)"; ctx.lineWidth = 1; ctx.setLineDash([4, 4]);
      ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
      ctx.setLineDash([]);
    }

    // static lattice targets
    ctx.fillStyle = "#6b7280";
    const sq = clamp(0.35 * scale, 3, 7);
    for (let k = 0; k < this.nStatic; k++) {
      const [px, py] = toPx(this.statics[2 * k], this.statics[2 * k + 1]);
      ctx.fillRect(px - sq / 2, py - sq / 2, sq, sq);
    }

    // oscillating target
    {
      const [px, py] = toPx(this.oscX, this.oscY);
      ctx.fillStyle = "rgba(242,104,95,0.25)";
      ctx.beginPath(); ctx.arc(px, py, 14, 0, Math.PI * 2); ctx.fill();
      ctx.fillStyle = "#f2685f";
      ctx.beginPath(); ctx.arc(px, py, 7, 0, Math.PI * 2); ctx.fill();
    }

    // agents
    const N = this.N;
    if (this.showHeadings) {
      ctx.strokeStyle = "#5aa4f2"; ctx.lineWidth = 1; ctx.globalAlpha = 0.85;
      for (let i = 0; i < N; i++) {
        const x0 = this.X[2 * i], y0 = this.X[2 * i + 1];
        const [px0, py0] = toPx(x0, y0);
        const [px1, py1] = toPx(x0 + this.U1[2 * i], y0 + this.U1[2 * i + 1]);
        ctx.beginPath(); ctx.moveTo(px0, py0); ctx.lineTo(px1, py1); ctx.stroke();
      }
      ctx.globalAlpha = 1;
    }
    ctx.fillStyle = "#f4f5f7";
    for (let i = 0; i < N; i++) {
      const [px, py] = toPx(this.X[2 * i], this.X[2 * i + 1]);
      ctx.beginPath(); ctx.arc(px, py, 2, 0, Math.PI * 2); ctx.fill();
    }
  }
}

if (typeof module !== "undefined") module.exports = { VTPSim, PARAMS, delaunayNeighbors, expReciprocal };
