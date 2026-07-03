# FlowSim test harness

Run everything from WSL via `tools/mrun`. This folder assumes MATLAB R2024a
on the Windows side (see `../tools/mrun` and the FlowSim vectorization project
AEGIS policy at `<axon>/my-axon/dev-projects/flowsim-vectorize/_policy.md`).

## Layout

```
tests/
├── README.md            # this file
├── run_all.m            # top-level runner — invokes smoke/*, unit/*, prints report
├── run_smoke.m          # only the smoke tests (fast: env sanity, class loads)
├── run_unit.m           # only the unit tests (medium: per-function oracle diffs)
├── helpers/             # test infrastructure — assertions, oracle-diff, fixtures
│   ├── fs_expect.m      # simple assertion (fs_expect(cond, msg))
│   ├── fs_reltol.m      # relative-tolerance comparator
│   ├── fs_frob.m        # Frobenius relative diff on sparse/dense matrices
│   ├── fs_setup.m       # standard test setup: addpath, seed rng, capture t0
│   ├── fs_teardown.m    # cleanup: report elapsed, restore rng
│   └── fs_load_mesh.m   # load a canonical mesh by name (M8, HermelineQuad_12, ...)
├── smoke/               # fast (seconds) — sanity checks, no numerical claims
│   ├── smoke_env.m              # MATLAB version, WSL cwd, paths configured
│   ├── smoke_class_hierarchy.m  # verify OOP hierarchy state (which classes load)
│   ├── smoke_mesh_load.m        # every .msh in meshes/ loads and reports node/elem counts
│   └── smoke_startdat.m         # Start.dat parses; report numcase, pmethod, phasekey
├── unit/                # medium (10-60 s per test) — per-function correctness
│   ├── unit_preprocessormod.m   # mesh build produces expected shapes/counts
│   ├── unit_lpew2_reference.m   # baseline LPEW2 output on a small mesh → captures in golden/
│   ├── unit_assembly_mpfad.m    # MPFA-D assembly nnz, condest, first-5 pressure values
│   └── unit_richards_caso439.m  # one Picard iteration of Caso439 completes + numerical fingerprint
├── golden/              # captured baseline outputs (committed; regenerated with --update)
│   └── (populated by unit tests on first run)
└── data/                # small test fixtures (tiny meshes, canned Start.dat variants)
    └── (empty until tests need them)
```

## Running

```bash
# Full harness (smoke + unit)
tools/mrun -c $(pwd) tests/run_all.m

# Just smoke tests (fast, ~1 min total)
tools/mrun -c $(pwd) tests/run_smoke.m

# One specific test
tools/mrun -c $(pwd) tests/smoke/smoke_env.m

# One unit test with auto-log
tools/mrun -L -c $(pwd) tests/unit/unit_lpew2_reference.m
```

Each test script:
1. Calls `fs_setup` at the start (paths, rng seed, timer).
2. Uses `fs_expect(cond, msg)` for assertions — first failure prints then continues; final summary reports pass/fail count.
3. Ends with `fs_teardown` which prints `TEST OK` or `TEST FAIL: N failures` and exits with 0 or 1.
4. Numerical unit tests read/write baselines in `golden/<test-name>.mat`. On first run they capture. On subsequent runs they compare within tolerance.

## Tolerances

Default per `fs_frob.m`:
- Sparse-matrix Frobenius relative: `1e-12`
- Dense-vector L2 relative: `1e-10`
- Scalar relative: `1e-10`
- `condest` allowed drift: `10× ratio`

Override per-call via `fs_frob(A, B, tol)`.

## Design choices

- **No third-party test framework** — MATLAB has `matlab.unittest` but adds
  overhead; for a code-dev study harness, plain scripts + `fs_expect` are
  simpler and diff cleanly under version control.
- **Golden files are `.mat`** (binary) — deterministic, small, git-friendly for
  scalars but ignored for large fields via `.gitignore` (see next section).
- **Tests are runnable one-at-a-time** so a failing unit doesn't halt the batch.
- **No side-effects outside `tests/`** — every test creates its scratch in
  `tests/data/tmp/` and cleans up in `fs_teardown`.

## Gitignore

Add to `.gitignore`:
```
tests/data/tmp/
tests/golden/*.big.mat   # over-100KB baselines stay local
```
