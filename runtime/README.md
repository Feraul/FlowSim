# runtime/ — Active-runtime code paths

Files here are called during a live FlowSim simulation. Grouped by role in
the pipeline. All four subdirs are added to the MATLAB path by
`flowsim_init.m` via `genpath`.

## Directory layout

```
runtime/
├── preproc/    mesh + method preprocessing (called once per run)
├── time/       time-stepping drivers (called once per timestep)
├── plug/       PLUG_* callbacks (called every substep — the "hot" plane)
└── util/       helpers reused across the pipeline
```

## preproc/ — one-shot setup

| File | Called by | Role |
|---|---|---|
| `preprocessormod.m` | `main.m` | reads `Start.dat`, builds `env.geometry` from `.msh` |
| `preprocessmethod.m` | `main.m` | per-`pmethod` premethod struct setup |
| `preRE.m` | `hydraulic_RE` | Richards-specific preprocessing |
| `preSaturation.m` | `IMPES` | saturation-specific preprocessing |
| `preconcentration.m` | `IMPEC` | concentration-specific preprocessing |

## time/ — time-driver dispatch + step loops

| File | phasekey | Role |
|---|---|---|
| `setmethod.m` | (any) | dispatch to the phasekey-appropriate driver |
| `hydraulic.m` | 4 | steady groundwater / hydraulic-head solver |
| `hydraulic_RE.m` | 6 | Richards equation transient solver |
| `IMPES.m` | 1 | Implicit Pressure Explicit Saturation |
| `IMPEC.m` | 5 | Implicit Pressure Explicit Concentration |
| `IMHEC.m` | — | variant (hydraulic-head + concentration) |

## plug/ — per-step callbacks (extension points)

These are the physics plugins — replace one to change the physics without
touching the time driver. All accept `(env, parms, tempo)` and return the
scalar / matrix the driver needs for the current substep.

| File | Callback for | Returns |
|---|---|---|
| `PLUG_kfunction.m` | permeability tensor | `env.config.kmap` for current `h` |
| `PLUG_bcfunction.m` | pressure BC | `env.config.nflag`, `nflagface` |
| `PLUG_bcfunction_con.m` | concentration BC | ditto (for tracers) |
| `PLUG_bcfunction_con_mpfa_o_fps.m` | MPFA-O concentration BC | ditto |
| `PLUG_sourcefunction.m` | source term (wells, injection) | element-wise source |
| `PLUG_dfunction.m` | dispersion / diffusion | tensor |
| `PLUG_Gfunction.m` | gravity | element-wise vector |

## util/ — cross-cutting helpers

| File | Role |
|---|---|
| `solver.m` | linear-system solve wrapper (backslash + preconditioner options) |
| `solvePressure.m` / `solvePressure_TPFA.m` | pressure-eq wrapper |
| `solveSaturation.m` | saturation-eq wrapper |
| `addsource.m` | inserts wells/source terms into the assembled system |
| `postprocessor.m` | writes VTK / TecPlot / mat output |
| `plotandwrite.m` / `plotandwrite_pressfield.m` | plot + persist result fields |
| `soil_properties.m` | Richards `dtheta/dt` accumulation term |
| `thetafunction.m` / `theta_n.m` | Richards water-content models |
| `gravitation.m` / `calcnormk.m` | gravity + permeability norm helpers |
| `applyinicialcond.m` / `attribinitialcond.m` / `IC.m` | initial-condition setup |
| `setrestartinicond.m` / `getrestartdata.m` | restart-file handling |
| `defineWells.m` | well-source definition helper |
| `benchmark.m` | benchmark utility (bench identifier lookup) |

## Path-setup order

`flowsim_init.m` adds paths in this order (later = higher precedence):

1. `+fs/`  — vectorized packages (highest precedence — shadow legacy)
2. `base/`, `solvers/`, `factories/`, `simulacoes/`, `benchmarks/` — OOP contract
3. `runtime/**` (this tree, via `genpath`)
4. `legacy/**` — legacy fallbacks (lowest precedence)

If a symbol exists both in `+fs/…` and in `legacy/…`, MATLAB resolves to the
`+fs/` version because it comes earlier on the path.

## Relationship to `legacy/`

Files here are the **live** copies used at runtime. Their vectorized twins
(where they exist) live in `+fs/` and shadow them transparently. Files
that are shadowed and still called from here are:

- `preprocessmethod` → calls `fs.lpew.v2.preLPEW2` (transparent shadow of
  `ferncodes_Pre_LPEW_2_vect`) when `+fs/` is on the path
- Assembly calls hit `fs.assembly.<pmethod>.build` first; for `mpfad` and
  `tpfa` these are fully vectorized, others delegate to legacy

## See also

- `../README.md` — repository-level overview
- `../docs/code-map.md` — function-level "where does X live?" table
- `../docs/vectorization-guide.md` — recipe for extending `+fs/`
- `../legacy/README.md` — index of legacy clusters
