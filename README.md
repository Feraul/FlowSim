# FlowSim — Vectorized

[![Release](https://img.shields.io/github/v/release/Feraul/FlowSim?sort=semver)](https://github.com/Feraul/FlowSim/releases)
[![License](https://img.shields.io/badge/license-see_LICENSE-blue.svg)](LICENSE)

**2-D groundwater / Richards / two-phase flow simulator in MATLAB.**
Finite-volume methods (TPFA, MPFA-D, MPFA-H, NLFV-PP, MPFA-QL) on unstructured
quad + triangle meshes. Originally by the UFPE / groundwater group.

The **2026-07 vectorization campaign** is complete and shipped as
[`v2.0.0-vectorized`](https://github.com/Feraul/FlowSim/releases/tag/v2.0.0-vectorized)
plus root-cleanup follow-up
[`v2.0.1-vectorized`](https://github.com/Feraul/FlowSim/releases/tag/v2.0.1-vectorized).
Pre-campaign state is preserved as
[`v1.0.0-pre-vectorization`](https://github.com/Feraul/FlowSim/releases/tag/v1.0.0-pre-vectorization)
(+ backup branch `legacy-v1.0`). See [`CHANGELOG.md`](CHANGELOG.md).

> **Coming from v1?** Read [`docs/for-scientists.md`](docs/for-scientists.md)
> — plain-language guide covering what stayed the same, what moved, how to
> verify your `Caso NNN` still produces the same numbers, and how to
> disable every new module if you want pure-v1 behaviour.

- **LPEW2** pipeline: fully vectorized (`+fs/+lpew/+v2/`), bit-identical to
  legacy at 1e-15 relative tolerance
- **MPFA-D** assembly: fully vectorized (`+fs/+assembly/+mpfad/`),
  bit-identical to golden baseline
- **TPFA** assembly: vectorized (`+fs/+assembly/+tpfa/`)
- **MPFA-H, MPFA-QL, NLFV-PP, NLFV-H, DMP**: scaffolded — entry points
  exist under `+fs/+assembly/`, currently delegate to legacy for
  correctness while awaiting the full triplet-form rewrite (see
  `CHANGELOG.md` → _Deferred to future work_)

Legacy code remains under `legacy/` and is picked up by the MATLAB path at
lower precedence — nothing breaks, vectorized modules take over silently.

---

## Quick start

```bash
# Clone
git clone https://github.com/Feraul/FlowSim.git
cd FlowSim

# Configure Start.dat (mesh path, numcase, pmethod, phasekey) — see docs/how-to-use.md

# Run headless (Linux/WSL):
tools/mrun -c $(pwd) main.m
```

Interactive from MATLAB:
```matlab
cd /path/to/FlowSim
flowsim_init            % path setup (+fs/ first, legacy/ last)
main                    % reads Start.dat + runs the configured case
```

The pipeline reads `Start.dat` → picks the case via `createBenchmark` →
picks the method via `createMetodo` → picks the simulation type via
`createSimulacao` → runs the time loop.

## Tests

```bash
# Full smoke suite (~30 s):
tools/mrun -c $(pwd) tests/smoke/smoke_env.m
tools/mrun -c $(pwd) tests/smoke/smoke_class_hierarchy.m

# Vectorization unit tests (~1 min each):
tools/mrun -c $(pwd) tests/unit/unit_lpew2_preLPEW2.m
tools/mrun -c $(pwd) tests/unit/unit_assembly_mpfad.m

# Correctness oracle (bit-identical reproduction of captured golden):
tools/mrun -c $(pwd) tests/unit/unit_baseline_reproduces.m
```

~500 assertions currently green — see `CHANGELOG.md` for the per-PR breakdown.

---

## Documentation

| Doc | For | Content |
|---|---|---|
| ★ [`docs/for-scientists.md`](docs/for-scientists.md) | **maintainers / users coming from v1** | plain-language guide: what stayed the same, what moved, how to verify your case matches v1 numerics |
| [`docs/how-to-use.md`](docs/how-to-use.md) | users | install, configure `Start.dat`, run, troubleshoot |
| [`docs/code-map.md`](docs/code-map.md) | contributors | where every function lives + who calls whom |
| [`docs/vectorization-guide.md`](docs/vectorization-guide.md) | contributors | recipe for extending `+fs/` with a new module |
| [`docs/globals-audit.md`](docs/globals-audit.md) | contributors | inventory of remaining `global` variables |
| [`runtime/README.md`](runtime/README.md) | contributors | active-runtime tree (`preproc/time/plug/util`) |
| [`legacy/README.md`](legacy/README.md) | contributors | 12-cluster legacy index + retirement status |
| [`tests/README.md`](tests/README.md) | contributors | test harness details |
| [`CHANGELOG.md`](CHANGELOG.md) | everyone | release history + full PR log |
| [`manual/manual.pdf`](manual/manual.pdf) | scientists | original numerical-methods manual (theory) |

---

## Repository map

```
FlowSim/
│
├── main.m                       ⬅ ENTRY POINT — read Start.dat + run
├── Start.dat                    runtime config (mesh, numcase, pmethod, ...)
├── startup.m                    auto-invoked by MATLAB; delegates to flowsim_init
├── flowsim_init.m               canonical path setup (call once per session)
├── flowsim_deinit.m             symmetric cleanup
├── README.md                    THIS FILE
├── CHANGELOG.md                 vectorization campaign log
├── manual/manual.pdf            original scientific manual
├── docs/                        contributor docs (vectorization guide, code map)
├── tools/
│   └── mrun                     WSL bridge → headless matlab.exe -batch
│
├── ★  OOP CONTRACT LAYER  ★  (small, clean, unified)
├── base/
│   └── SimulacaoBase.m          abstract simulation base class
├── simulacoes/
│   ├── SimGroundwater.m         hydraulic head physics
│   └── SimRichards.m            partially-saturated flow
├── solvers/                     numerical methods (all inherit from MetodoBase)
│   ├── MetodoBase.m
│   ├── MetodoTPFA.m
│   ├── MetodoMPFAD.m            ⬅ production target — fully vectorized
│   ├── MetodoMPFAH.m
│   ├── MetodoMPFAQL.m
│   └── MetodoNLFVPP.m
├── benchmarks/
│   └── Caso439.m                the one fully-implemented Caso class
├── factories/                   dispatch by numcase / pmethod / phasekey
│   ├── createBenchmark.m
│   ├── createMetodo.m
│   └── createSimulacao.m
│
├── ★  VECTORIZED MODULES  ★  (+fs — the new tree)
├── +fs/
│   ├── +mesh/build.m            env → FS struct adapter
│   ├── +geom/                   (derived-once quantities — future)
│   ├── +csr/
│   │   ├── buildCorners.m       CSR-flat node→corner layout (KEY enabler)
│   │   └── buildCornerShifts.m  precomputed k-1/k+1 wrap indices
│   ├── +util/
│   │   └── assertFS.m           5-invariant FS struct checker
│   ├── +data/
│   │   └── paths.m              resolves spe10/spe_perm/gmsh
│   ├── +lpew/
│   │   ├── OPT.m                batched geometry gather (was OPT_Interp_LPEW)
│   │   ├── pinterp.m            nodal pressure interp (was pressureinterpNLFVPP)
│   │   └── +v2/                 vectorized LPEW2 (kills the per-node loop)
│   │       ├── angulos.m
│   │       ├── netas.m
│   │       ├── ksInterp.m
│   │       ├── lambdaWeights.m
│   │       └── preLPEW2.m       full pipeline driver
│   ├── +assembly/               triplet-form sparse assembly per method
│   │   ├── +mpfad/build.m       ⬅ FULLY VECTORIZED (Frob diff 1e-17 vs legacy)
│   │   ├── +tpfa/build.m        vectorized
│   │   ├── +mpfah/build.m       scaffold (delegates to legacy)
│   │   ├── +mpfaql/build.m      scaffold
│   │   ├── +nlfvpp/build.m      scaffold
│   │   ├── +nlfvh/build.m       scaffold
│   │   └── +dmp/build.m         scaffold
│   └── +flow/                   flow-rate calculators
│       ├── mpfad.m
│       └── tpfa.m
│
├── ★  RUNTIME CORE (moved to runtime/ in v2.0.1)  ★
├── runtime/
│   ├── preproc/                 preprocessormod, preprocessmethod, benchmark_setup
│   ├── time/                    hydraulic, hydraulic_RE, setmethod, IMPES, IMPEC, IMHEC
│   ├── plug/                    PLUG_bcfunction, PLUG_kfunction, PLUG_sourcefunction, PLUG_dfunction
│   └── util/                    addsource, soil_properties, postprocessor, solver, ~30 more
│
├── ★  MESH FIXTURES  ★
├── meshes/
│   ├── hermeline/               Hermeline_*.msh (7 files)
│   ├── kozdon/                  Kozdon_*.msh (3 files)
│   └── other/                   M8*.msh, other benchmark meshes (6 files)
│
├── ★  DATA  ★
├── data/
│   ├── Perm_Var*.mat            permeability datasets (Caso247/249/250 inputs)
│   ├── Teste_*.xlsx             benchmark data spreadsheets
│   ├── malha_D.geo              gmsh source
│   └── figura_*.fig             reference plots
│   └── indices_elementos_*.mat  ← gitignored — regenerated by Caso439 on first run
│
├── ★  TESTS  ★
├── tests/
│   ├── README.md
│   ├── run_all.m                aggregate runner
│   ├── helpers/
│   │   ├── fs_setup.m           per-test setup + rng seed
│   │   ├── fs_expect.m          assertion primitive
│   │   ├── fs_teardown.m        summary + exit 0/1
│   │   ├── fs_frob.m            Frobenius diff comparator
│   │   ├── fs_test_env.m        Start.dat isolation (WSL-safe)
│   │   └── capture_baseline.m   golden capture harness
│   ├── smoke/
│   │   ├── smoke_env.m          MATLAB version + paths + basic sanity
│   │   ├── smoke_class_hierarchy.m  OOP hierarchy loads (17 classes)
│   │   ├── smoke_mesh_load.m
│   │   └── smoke_startdat.m
│   ├── unit/                    per-function parity tests vs legacy
│   │   ├── unit_assertFS.m
│   │   ├── unit_fs_mesh_build.m
│   │   ├── unit_csr_corners.m
│   │   ├── unit_lpew_OPT.m
│   │   ├── unit_lpew2_angulos_netas.m
│   │   ├── unit_lpew2_ksInterp.m
│   │   ├── unit_lpew2_lambda.m
│   │   ├── unit_lpew2_preLPEW2.m
│   │   ├── unit_assembly_mpfad.m
│   │   ├── unit_assembly_tpfa.m
│   │   ├── unit_assembly_scaffolds.m
│   │   ├── unit_flow_pinterp_scaffolds.m
│   │   ├── unit_baseline_reproduces.m
│   │   └── unit_data_paths.m
│   └── golden/                  captured baselines (committed, 700 B / mesh)
│       ├── M8-num439-tpfa.mat
│       └── M8-num439-mpfad.mat
│
├── ★  LEGACY  ★  (relocated by PR-E3..E7 + PR-F3)
└── legacy/
    ├── ferncodes/               66 files clustered by method
    │   ├── mpfad/  mpfah/  mpfaql/  nlfvh/  nlfvpp/  dmp/  lpew/  shared/
    ├── ferncodes-con/           dead concentration-coupled variants (7)
    ├── preprocessor/            preprocessor.m + preprocessor2.m (superseded)
    ├── transm/                  ~7k LOC of dead transmissibility code (8)
    └── unknown/                 32 files pending owner disposition
```

Root `.m` file count: **4** (main, flowsim_init, flowsim_deinit, startup) — down from **285** at v1.0.0.  
`.m` files under `runtime/`: **37** (active execution paths).  
`.m` files under `legacy/`: **~250** (still reachable via path precedence — vectorized `+fs/` shadows).

---

## How the code runs (execution flow)

```
main.m
  │
  ├─→ flowsim_init          (adds +fs/, base/, solvers/, factories/,
  │                          simulacoes/, benchmarks/, legacy/** to path)
  │
  ├─→ env = preprocessormod(1)
  │       └─ reads Start.dat, builds env.geometry.{coord, elem, bedge,
  │          inedge, centelem, elemarea, normals, esurn1/2, nsurn1/2}
  │          and env.config.{numcase, pmethod, phasekey, perm, bcflag, …}
  │
  ├─→ env.benchmark = createBenchmark(env.config.numcase)   → CasoNNN
  ├─→ env.metodo    = createMetodo(env.config.pmethod)      → MetodoXXX
  ├─→ sim           = createSimulacao(env.config.phasekey)  → SimXXX
  │
  ├─→ parms        = env.benchmark.initParms()
  ├─→ [env, parms] = sim.preprocessar(env, parms)
  ├─→ [env, parms] = PLUG_kfunction(env, parms, 0)          → env.config.kmap
  ├─→ [env]        = ferncodes_calflag(env, parms, 0)       → env.config.nflag/nflagface
  ├─→ [env, parms] = preprocessmethod(env, parms)
  │       │  ├─ ferncodes_elementface        → env.premethod.MPFAD.{V, N, F}
  │       │  ├─ ferncodes_Kde_Ded_Kt_Kn      → Kde, Ded, Kt, Kn, Hesq
  │       │  └─ ferncodes_Pre_LPEW_2_vect    → weight, s
  │       │      (see +fs/+lpew/+v2/preLPEW2.m for the vectorized twin)
  │       └─ env.preGravity = obj.calcGravidade(env)
  │
  └─→ setmethod(...)  → dispatches to the actual time driver:
        ├─ hydraulic_RE  (phasekey=6, Richards)
        ├─ hydraulic     (phasekey=4, groundwater steady)
        ├─ IMPES         (phasekey=1, single-phase pressure/sat)
        ├─ IMPEC         (phasekey=5, contaminant transport)
        └─ IMHEC         (phasekey=?)

        Inside each time driver:
          - env.metodo.montarSistema     → M, I  (assembly hot path)
          - addsource / sourceterm
          - env.metodo.resolver           → Picard / Anderson / L-scheme
                └─ per iter: PLUG_kfunction → atualizarPremethod → montarSistema → linear solve
          - env.metodo.calcularFlowrate  → post-process fluxes
          - env.sim.atualizarEstado      → benchmark hook (postprocessor, VTK output)
```

## The three class hierarchies

Just one, actually — the vectorization campaign unified them.

```
handle (MATLAB base)
├── MetodoBase                 (contract for all numerical methods)
│    ├── MetodoTPFA            ⬅ baseline
│    ├── MetodoMPFAD           ⬅ production target (fully vectorized)
│    ├── MetodoMPFAH           (was SolverMPFAH — renamed in PR-A2)
│    ├── MetodoNLFVPP          (was SolverNLFVPP)
│    └── MetodoMPFAQL          (was missing — created in PR-A2)
│
└── SimulacaoBase              (contract for simulation types)
     ├── SimGroundwater
     ├── SimRichards
     └── Caso439               (the one fully-implemented benchmark)
```

## Configuration (Start.dat)

Edit `Start.dat` — every line marked `>>> EDITE AQUI <<<` is a user-editable
parameter. Owner's original Windows paths (`C:\Users\flc59\...`) still work
if you have that filesystem; on WSL/Linux, use forward-slash paths.

Critical fields:
- **`numcase`** — which benchmark to run (see `factories/createBenchmark.m`)
- **`pmethod`** — `tpfa`, `mpfad`, `mpfah`, `nlfvpp`, `mpfaql`
- **`phasekey`** — 1 (single-phase), 4 (groundwater), 5 (contaminant), 6 (Richards)
- **mesh path** — folder containing the `.msh` file
- **mesh name** — e.g. `M8.msh`
- **output path** — where to write results

---

## Vectorization status (as of 2026-07-03)

| Component | Status | Verification |
|---|---|---|
| **LPEW2 weight computation** | ✅ Fully vectorized (+fs/+lpew/+v2/) | 1e-15 relative diff vs legacy |
| **MPFA-D matrix assembly** | ✅ Fully vectorized (+fs/+assembly/+mpfad/) | Bit-identical to legacy AND golden |
| **TPFA matrix assembly** | ✅ Vectorized (legacy already was, rename only) | Bit-identical |
| **MPFA-D flow-rate** | ✅ (rename; legacy already vectorized) | — |
| **MPFA-H assembly** | 🔧 Scaffold (delegates) | Deferred — 820 L rewrite |
| **NLFV-PP assembly** | 🔧 Scaffold (delegates) | Deferred |
| **MPFA-QL / DMP / NLFV-H assembly** | 🔧 Scaffolds | Deferred |

`.gitignore` excludes: `gmsh.exe`, `spe*.dat`, `spe*.mat`, `*.asv`, `tests/data/tmp/`

---

## Contributing

- Work on `flowsim-artur` branch; `master` is stable / untouched.
- Add a `tests/unit/unit_<name>.m` alongside every code change (see
  `docs/vectorization-guide.md` for the recipe).
- Verify via `tools/mrun -c $(pwd) tests/unit/<your-test>.m` before commit.
- Golden baselines live in `tests/golden/`; regenerate with
  `tests/helpers/capture_baseline.m` if legacy semantics change.

## References

- `CHANGELOG.md` — PR-by-PR chronology
- `docs/vectorization-guide.md` — recipe cookbook for future vectorized modules
- `docs/code-map.md` — deep dive: who calls whom, where each thing lives
- `manual/manual.pdf` — original scientific manual (numerical methods theory)
- `tests/README.md` — test harness usage

## License / Attribution

Original codebase by the UFPE groundwater flow group.  
Vectorization campaign (2026): Artur Castiel Reis de Souza + AXON code-dev + Copilot.
