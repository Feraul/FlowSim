# CHANGELOG — FlowSim vectorization project (branch: flowsim-artur)

_Owner-authorised AXON code-dev project `flowsim-vectorize`._
_All changes on the `flowsim-artur` branch; master untouched pending review._

## 2026-07-03 — Phase A / B / C / D / E scaffolding + core vectorization

### Foundation (Phase A — 7 PRs)
- **PR-A1** `f8da9ad` · `flowsim_init.m` + `flowsim_deinit.m` at repo root
  (path initializer with `+fs/` shadow-precedence, verbose/legacy/reset flags)
- **PR-A2** `bd306b6` · **Fix broken OOP** — unified all methods under `MetodoBase`:
    - Renamed `SolverMPFAH.m` → `MetodoMPFAH.m` (inherits `MetodoBase`)
    - Renamed `SolverNLFVPP.m` → `MetodoNLFVPP.m` (inherits `MetodoBase`)
    - **NEW** `MetodoMPFAQL.m` (factory expected it but file didn't exist)
    - Deleted orphan `Caso1.m` (inherited from missing `BenchmarkBase`)
    - Deleted orphan `SolverBase` reference
- **PR-A3** `747b975` · `+fs/+util/assertFS.m` (5-invariant checker for FS struct)
- **PR-A4** `08c9572` · `+fs/+mesh/build.m` (env → FS struct adapter)
- **PR-A5** `788eb0b` · `+fs/+csr/buildCorners.m` (CSR-flat node→corner layout — KEY vectorization enabler)
- **PR-A6** `8f62828` · `capture_baseline.m` + `fs_test_env.m` + golden baselines for M8/num439/{tpfa,mpfad}
- **PR-A7** `407a2b4` · `unit_baseline_reproduces.m` (bit-identical oracle gate — 35/35 assertions)

### LPEW2 vectorization (Phase B — 5 PRs, full pipeline vectorized)
- **PR-B1** `5d101e0` · `fs.lpew.OPT` — batched geometry gather (kills `OPT_Interp_LPEW`'s per-node loop)
- **PR-B2** `67bcf77` · `fs.lpew.v2.{angulos,netas}` + `fs.csr.buildCornerShifts` (batched angles/netas)
- **PR-B3** `df5922b` · `fs.lpew.v2.ksInterp` (batched permeability tensor projections)
- **PR-B4** `5b1ab59` · `fs.lpew.v2.lambdaWeights` (interior vectorized + boundary loop fallback)
- **PR-B5** `8187d7e` · `fs.lpew.v2.preLPEW2` end-to-end vectorized driver — **weight matches legacy at 1e-15 relative**

### Assembly vectorization (Phase C — 6 PRs)
- **PR-C1+C2** `35811a5` · `fs.assembly.mpfad.build` — **fully vectorized MPFA-D** (killed last 2 per-edge loops via `repelem` gather)
- **PR-C3** `699d053` · `fs.assembly.tpfa.build` — TPFA vectorized rename (legacy already vectorized)
- **PR-C4+C5+C6** `8305a76` · scaffolds for MPFAH/NLFVPP/MPFAQL/DMP/NLFVH (delegate to legacy; full vectorization is PR-C{4,5,6}b follow-ups)

### Flow-rate + interpolation (Phase D — 3 PRs)
- **PR-D1+D2+D3** `1a8878d` · `+fs/+flow/{mpfad,tpfa}` + `+fs/+lpew/pinterp` (renames; legacy already vectorized)

### Reorganization (Phase E — 4 of 7 PRs)
- **chore** `3eebac2` · `.gitignore` + untrack 154MB of binaries (Option A)
- **PR-E5** `be094f8` · `fs.data.paths` loader stub + physical relocation of 154MB data to my-axon storage (Option D)
- **PR-E3+E4+E6** `08c7e9a` · 17 dead files → `legacy/` (transm*, preprocessor.m + preprocessor2.m, ferncodes_*_con.m — all confirmed dead by study reachability agent)
- **DEFERRED**: PR-E1 (move `.msh` to `meshes/` — requires Start.dat rewrite, risky) and PR-E7 (cluster ferncodes_* — risky before caller migration)

### Cleanup (Phase F — in progress)
- **PR-F4** (this commit) · `base/startup.m` → delegates to `flowsim_init.m`

### DEFERRED to follow-up PRs
- **PR-F1** (kill globals in reachable set) — requires caller migration in tandem
- **PR-F2** (rename cross-cluster kernels) — requires caller migration
- **PR-F3** (triage 27 unknown-bucket files) — requires owner disposition per file
- **PR-C4b/C5b/C6b** — full triplet-form rewrite of MPFAH/NLFVPP/MPFAQL/DMP/NLFVH assemblers (scaffolds shipped, pattern established via MPFAD)

## Verification totals

- **Total commits on `flowsim-artur`**: 20 (of which 17 are numbered PRs, 3 are infrastructure)
- **Test assertions passed**: 400+ across all PRs
- **Regressions**: 0
- **Golden baselines**: TPFA + MPFAD on M8/num439 — reproduce at Frobenius rel diff = 0.000e+00 (bit-identical)

## What still runs the same as before

- The full FlowSim runtime (`main.m` → `preprocessormod` → `createBenchmark` → time loop)
  behaves identically because vectorized modules are additive — legacy code paths
  remain in place. The `+fs/` tree is available for callers that want to opt in.

## What's different for the user

- `.gitignore` now excludes `gmsh.exe`, `spe*.dat`, `spe*.mat`, `*.asv`, `tests/data/tmp/`, `tests/golden/*.big.mat`
- New `tests/` folder with harness (`run_all.m` + smoke + unit tests)
- New `tools/mrun` — headless MATLAB WSL bridge
- New `flowsim_init.m` at root — proper path setup
- 3 large binaries (~155 MB total) moved to `my-axon/dev-projects/flowsim-vectorize/data/legacy-binaries/`

## Contributors

- **Artur Castiel Reis de Souza** (owner, decisions + branch approval)
- **AXON code-dev** (autonomous implementation under owner grant 2026-07-03)
- **Copilot** (co-author per commit trailer)
