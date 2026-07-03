# How to use FlowSim

_User-facing guide. For contributors, see `README.md` + `docs/code-map.md` + `docs/vectorization-guide.md`._

## Prerequisites

- **MATLAB R2019b or later** (tested with R2024a on Windows via WSL)
- **Gmsh** (system install: `apt install gmsh` on Linux, `choco install gmsh` on Windows) — used only for editing mesh files, not required to RUN the simulator

## First-time setup

1. **Clone the repo** and check out the working branch:
   ```bash
   git clone https://github.com/Feraul/FlowSim.git
   cd FlowSim
   git checkout flowsim-artur
   ```

2. **Edit `Start.dat`** — the runtime config file. Every line marked
   `>>> EDITE AQUI <<<` is a user parameter. Most important:
   - **Output path** (line 33) — where results go. Use forward slashes on
     Linux/WSL, backslashes on Windows.
   - **Mesh path** (line 218) — folder containing `.msh` files.
   - **Mesh filename** (line 224) — e.g. `M8.msh`.
   - **`numcase`** — which benchmark case (see `factories/createBenchmark.m`
     for the full list).
   - **`pmethod`** — one of: `tpfa` / `mpfad` / `mpfah` / `nlfvpp` / `mpfaql`.
   - **`phasekey`** — physics family:
     - `1` = single-phase pressure/saturation
     - `4` = groundwater hydraulic head
     - `5` = contaminant + hydraulic head
     - `6` = Richards (partially saturated)

3. **Confirm SPE10 permeability data** (optional — only for SPE10 cases). By
   default they've been relocated out of the repo (`spe10.mat`, `spe_perm.dat`
   ~90 MB total). If a benchmark needs them:
   - Symlink from the location `fs.data.paths('spe10')` reports, OR
   - Set the environment variable `FS_DATA_DIR` to where you keep them.

## Running from MATLAB (interactive)

```matlab
cd /path/to/FlowSim
flowsim_init          % sets up all paths (only needed once per session)
main                  % reads Start.dat, runs the configured case
```

## Running from WSL / terminal (headless)

```bash
cd /path/to/FlowSim

# One-shot run of main.m:
tools/mrun -c $(pwd) main.m

# Or a specific script:
tools/mrun -c $(pwd) /path/to/myscript.m

# One-liner:
tools/mrun -e "disp('hi'); disp(version)"

# With timeout (default 30 min):
tools/mrun -t 3600 -c $(pwd) main.m

# With auto-log to /tmp/mrun-logs/:
tools/mrun -L -c $(pwd) main.m
```

`mrun` wraps `matlab.exe -batch` and:
- Filters harmless UNC change-notification warnings
- Auto-cd's into the script's directory (or a `-c` override)
- Propagates the exit code (0 = ok, non-0 = MATLAB error, 124 = timeout)

## Running the test suite

```bash
# Full smoke pass (fast, ~30s):
tools/mrun -c $(pwd) tests/smoke/smoke_env.m

# All smoke + unit (~5 min):
tools/mrun -c $(pwd) tests/run_all.m

# One specific test (with verbose progress):
tools/mrun -c $(pwd) tests/unit/unit_baseline_reproduces.m
```

All tests use `tests/helpers/fs_test_env.m` which isolates
Start.dat into a temp dir — your real `Start.dat` is never mutated by tests.

## Adding a new benchmark case

1. **Create `benchmarks/CasoNNN.m`** — inherit from `SimulacaoBase`, implement
   the abstract methods (see `base/SimulacaoBase.m` for the contract).
   Model after the existing `benchmarks/Caso439.m` — it's the only fully
   implemented one.

2. **Register in `factories/createBenchmark.m`**:
   ```matlab
   case NNN,  bench = CasoNNN();
   ```

3. **Update `Start.dat`**: set `numcase = NNN`.

4. **Optionally add a smoke test** in `tests/smoke/` that verifies the class
   loads (`meta.class.fromName('CasoNNN')`).

## Adding a new numerical method

1. **Create `solvers/MetodoXXX.m`** — inherit from `MetodoBase`, implement:
   - `preprocessar(env, parms)` — precompute mesh-invariant quantities
   - `atualizarPremethod(env, parms)` — refresh on new kmap (Richards)
   - `montarSistema(env, parms, dt)` — return `[M, I]`
   - `resolver(M, I, parms, env, tempo, dt, source_wells)` — solve (linear or Picard)
   - `calcularFlowrate(p, env, parms)` — postprocess

2. **Register in `factories/createMetodo.m`**:
   ```matlab
   case 'xxx',  metodo = MetodoXXX();
   ```

3. **Update `Start.dat`**: `pmethod = xxx`.

4. **If you want the vectorized package structure**: create
   `+fs/+assembly/+xxx/build.m` (delegate to legacy at first; incrementally
   port to triplet form).

## Common gotchas

- **"Undefined function or variable"** on some `Caso...` or `Metodo...` — the
  factory expects a class file that doesn't exist. Only `Caso439` + the 5
  `MetodoXXX` classes are wired in this branch.
- **"mkdir Access is denied"** during preprocessormod — Start.dat's output
  path is Windows-style but MATLAB can't create it (path doesn't exist or
  no write permission). Fix the output path.
- **"Cannot open file: ...\some_mesh.msh"** — Start.dat's mesh dir + mesh
  filename don't match a real file. Meshes now live at repo root (see
  `README.md` structure section).
- **Baseline test failing** — probably a numerical change slipped in. Compare
  `norm(M_new - M_leg, 'fro') / norm(M_leg, 'fro')`. If < 1e-12, adjust
  tolerance; if > 1e-12, there's a real bug.

## Where to look for problems

| Symptom | Look at |
|---|---|
| MATLAB can't find a function | `flowsim_init.m` — is legacy tree on path? |
| Wrong pressures / flowrates | Golden baseline test + `tests/golden/*.mat` |
| Class won't instantiate | `tests/smoke/smoke_class_hierarchy.m` |
| Path setup issues | Run `flowsim_init('verbose', true)` to see the addpath order |
| WSL path issues | Add `-v` to mrun to see the actual MATLAB_CMD |

## Cleaning up

- **Test scratch files** live in `/tmp/flowsim_test_env/` (recreated per test).
- **Auto-logs** from `mrun -L` live in `/tmp/mrun-logs/`.
- **MATLAB `.asv` autosaves** — gitignored, safe to delete.

## Getting help

- `manual/manual.pdf` — original scientific reference (numerical methods theory)
- `CHANGELOG.md` — what changed and when
- `docs/vectorization-guide.md` — recipe for extending the `+fs/` tree
- `docs/code-map.md` — where every function lives + who calls whom
- `tests/README.md` — test harness details
