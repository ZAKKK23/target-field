# plane-vtp-target-field (Voronoi Topological Perception)

The VTP agent model on the plane with a **lattice of static targets across
the whole plane** and **one oscillating target whose frequency you choose**,
as a MATLAB model and a dependency-free browser port of the same dynamics.

**[Live site →](#)** (enable GitHub Pages, see below)

## What it is

Agents follow the Voronoi Topological Perception (VTP) model:

- **Topological neighbors.** An agent's neighbors are its Delaunay-triangulation
  neighbors (the agents whose Voronoi cells touch its own).
- **Repulsion.** Away from the nearest neighbor, with magnitude
  *s* = transition(distance / L): 1 when touching, fading smoothly to exactly
  0 at distance L (`L = 1`).
- **Alignment.** With its Delaunay neighbors, weighted by **&nu;**. Neighbors
  whose heading already agrees with the agent's count for more.
- **Homing.** Toward the nearest target, with magnitude 1 &minus; *s*, so
  repulsion takes priority when an agent is crowded.
- **Speed governor.** Each agent moves at tanh(*l* / L) times its desired
  direction, where *l* is the distance to the boundary of its own Voronoi cell
  along its direction of motion. Crowded agents slow down; agents on the edge of
  the group have an unbounded cell and move at full speed.

### The target field

All targets are point targets, and each agent homes toward whichever one is
currently nearest to it.

- **Static targets** form a square lattice covering `[-H, H]^2`, one target
  every `spacing` units (defaults: `H = 20`, `spacing = 5`, i.e. 81 targets).
- **One oscillating target** swings along a line through a chosen center:

  ```
  position(t) = center + amplitude * sin(2*pi * frequency * t) * (cos(angle), sin(angle))
  ```

  `frequency` is in **cycles per time step** (period = 1 / frequency steps), and
  is the parameter you choose. `frequency = 0` freezes the target at its center.

Because agents choose the *nearest* target, the oscillating target mostly pulls
in agents that are closer to it than to any lattice point. A larger lattice
spacing or swing amplitude gives it a larger catchment area.

## Structure

```
index.html              site shell — tabs for about / live simulation / matlab version
assets/style.css        site styling
assets/sim.js           the VTP simulation engine (Delaunay neighbors, repulsion + alignment + homing, Voronoi speed governor, canvas rendering) — no external libraries
assets/app.js           page wiring — tabs, sliders, MATLAB source viewer
assets/matlab_src.json  bundled MATLAB source (for the in-page code viewer)
matlab/                 MATLAB implementation (dynamics.m and the files it depends on)
```

## Running the web version

No build step. Serve the folder (the MATLAB viewer tab loads
`assets/matlab_src.json` with `fetch`, which needs http rather than `file://`):

```
python3 -m http.server 8000
```

then visit `http://localhost:8000/`. The live controls cover the number of
agents, &nu;, the lattice spacing and half-width, and the oscillating target's
frequency (slider or exact number), amplitude, angle and center.

## Running the MATLAB version

```
cd matlab
matlab -r dynamics
```

or open `matlab/dynamics.m` in the MATLAB editor and run it. The target-field
settings are at the top of the file:

```matlab
grid_half    = 20;      % static lattice covers [-grid_half, grid_half]^2
grid_spacing = 5;       % one static target every grid_spacing units

osc_freq   = 0.01;      % oscillation frequency, cycles per time step  <-- choose this
osc_amp    = 8;         % amplitude of the swing
osc_center = [0 0];     % center of the swing
osc_angle  = pi/2;      % direction of the swing (0 = horizontal)
```

The lattice is plotted as grey squares and the oscillating target as a red
dot. `matlab/` contains `dynamics.m` plus the files it depends on:
`Target.m`, `nearestOnSegment.m`, `neighborhoods.m`, `transition.m`,
`alignTo.m`, `voronoiProjectToBoundary.m`, and (only needed if you set
`fdim = 2`) `voronoiForwardArea.m` and `poly_area.m`. MATLAB needs the
functions these use for geometry (`delaunayTriangulation`, `voronoiDiagram`,
`polyxpoly`, `polyshape`).

## Publishing to GitHub Pages

1. Create a new GitHub repository and push this folder to it (see commands
   below).
2. In the repo, go to **Settings → Pages**.
3. Under **Build and deployment**, set **Source** to `Deploy from a branch`,
   branch `main`, folder `/ (root)`.
4. Save — the site will be published at
   `https://<your-username>.github.io/<repo-name>/` within a minute or two.

```bash
git remote add origin https://github.com/<your-username>/<repo-name>.git
git push -u origin main
```

## Note on the JS port

The browser port implements the same VTP model as `dynamics.m` (with
`fdim = 1`, the default):

- The Delaunay triangulation is a small built-in Bowyer–Watson implementation.
  On a test configuration of 400 agents it produced exactly the same neighbor
  sets as SciPy's Delaunay.
- The speed governor is computed from the Delaunay neighbors rather than by
  intersecting a ray with the Voronoi polygon. A Voronoi cell is the
  intersection of the half-planes toward each Delaunay neighbor, so the
  distance to the boundary along a direction is the minimum of
  |d|² / (2 u·d) over neighbors with u·d > 0, and an unbounded cell gives
  tanh(∞) = 1. This is the same quantity `voronoiProjectToBoundary.m` returns.
  Repulsion, alignment, homing and the speed governor for one step matched an
  independent NumPy/SciPy implementation to machine precision.
- Initial conditions use a seeded generator with the same seeds as the MATLAB
  script (2 and 18), but MATLAB's `rng` cannot be reproduced in JavaScript, so
  the exact starting positions (and therefore trajectories) differ from
  MATLAB's. The behavior is statistically the same.
- Only `fdim = 1` is ported; the forward-area speed option (`fdim = 2`) is
  MATLAB-only.
- The MATLAB script was not run in MATLAB while it was being written.
