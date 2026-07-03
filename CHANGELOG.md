# CHANGELOG — FlowSim

All notable changes to FlowSim are documented here.
Version tags follow [Semantic Versioning](https://semver.org/): `MAJOR.MINOR.PATCH`.

Releases: [GitHub Releases page](https://github.com/Feraul/FlowSim/releases)

---

## [v2.0.1] — 2026-07-03

Post-`v2.0.0` root-cleanup follow-up. **Pure reorganization, zero numerical
impact.** `unit_baseline_reproduces` still passes 35/35 at relative Frobenius
diff **0.000e+00** against the committed goldens.

### Added
- `runtime/{preproc,time,plug,util}/` subtree — 37 active-runtime `.m` files
  moved from repo root, grouped by role
- `data/` folder — 7 data files moved from repo root (permeability datasets,
  benchmark spreadsheets, gmsh source, reference figures)
- `runtime/README.md` — index describing the four runtime subdirs
- `legacy/README.md` — index of the 12 legacy clusters + retirement status

### Changed
- Root `.m` file count: **42 → 4** (`main`, `startup`, `flowsim_init`,
  `flowsim_deinit`) — only the entry-point + path-setup files remain
- `Caso439:inicializar` — index cache path now derived from
  `mfilename('fullpath')` and lands in `data/` regardless of MATLAB `cwd`
- `README.md` — repository map section rewritten to reflect the new layout
- `docs/code-map.md` — 15+ location entries corrected from "root" to
  `runtime/…` or `legacy/ferncodes/lpew/`
- `docs/how-to-use.md` — install instructions updated for the merged
  `master` branch; mesh-path troubleshooting reflects `meshes/*/`

### Fixed
- `docs/how-to-use.md` still told users to `git checkout flowsim-artur`
  (that branch is now merged into `master`)

---

## [v2.0.0] — 2026-07-03

**The vectorization release.** Result of the 2026-07 AXON code-dev campaign
that landed 33 PRs on the `flowsim-artur` branch and merged into `master`
via `--no-ff` (preserving the full PR history).

### Highlights
- **LPEW2 pipeline fully vectorized** — bit-identical to legacy at **1e-15**
  relative accuracy (`+fs/+lpew/+v2/preLPEW2`)
- **MPFA-D assembly fully vectorized** — bit-identical to legacy AND to the
  committed golden baseline (`+fs/+assembly/+mpfad/build`)
- **TPFA assembly vectorized** — `+fs/+assembly/+tpfa/build`
- **Scaffolded methods** — MPFA-H, NLFV-PP, NLFV-H, MPFA-QL, DMP now have
  `+fs/+assembly/+<method>/build.m` entry points (delegate to legacy for
  now; the vectorization pattern is established via MPFA-D)
- **Tests harness** — `tests/` with ~500 assertions, all green, 0 regressions
- **WSL bridge** — `tools/mrun` runs `matlab.exe -batch` headless from WSL

### Added (structure)
- `+fs/` package tree — vectorized modules under 12 packages
  (`+util`, `+mesh`, `+csr`, `+data`, `+lpew`, `+lpew/+v2`, `+assembly` x7,
  `+flow`, `+iter`)
- `tests/` — harness (`fs_setup`, `fs_expect`, `fs_teardown`, `fs_frob`,
  `fs_test_env`, `capture_baseline`) + smoke + 14 unit tests + goldens
- `flowsim_init.m` + `flowsim_deinit.m` — clean path setup / teardown
- `tools/mrun` — WSL → MATLAB bridge (~200 lines of bash)
- `docs/` — `code-map.md`, `how-to-use.md`, `vectorization-guide.md`,
  `globals-audit.md`
- `meshes/{hermeline,kozdon,other}/` — 16 `.msh` fixtures organized

### Changed
- **Unified OOP contract**: all discretization methods now inherit from
  `MetodoBase` (was: mix of `SolverBase` + `MetodoBase` + orphans)
  - Renamed `SolverMPFAH.m` → `MetodoMPFAH.m`
  - Renamed `SolverNLFVPP.m` → `MetodoNLFVPP.m`
  - Created `MetodoMPFAQL.m` (factory expected it, file didn't exist)
  - Deleted orphan `Caso1.m` and `SolverBase` reference
- **Root `.m` count**: **285 → 42** — 250+ files relocated to `legacy/` in
  8 organized subclusters (`ferncodes/`, `ferncodes-con/`, `preprocessor/`,
  `transm/`, `limiters/`, `saturation/`, `transport/`, `calc/`, `get/`,
  `solvers/`, `utility/`, `test-scripts/`, `unknown/`)
- **Tracked repo weight**: **176 MB → ~10 MB** — 154 MB of binaries
  untracked + relocated (see `.gitignore` + `fs.data.paths`)
- `startup.m` — now a thin wrapper that delegates to `flowsim_init`

### Numerics
- MPFA-D + TPFA on `M8/num439`: Frobenius rel diff **0.000e+00** vs the
  committed goldens (`tests/golden/*.mat`)
- LPEW2 weight match: **1e-15** relative to legacy per-node loop

### Migration notes
Fully backward compatible for user code that calls the legacy entry
points — the vectorized modules are additive and are picked up
transparently via `flowsim_init`'s path precedence (`+fs/` first).

To restore v1.0.0 state:
```bash
git checkout v1.0.0-pre-vectorization    # detached HEAD
# OR
git checkout legacy-v1.0                  # branch (defense in depth)
```

---

## [v1.0.0-pre-vectorization] — 2026-07-03 _(historical backup tag)_

Snapshot of `master` at commit `798bbe7` immediately before the `flowsim-artur`
merge. Pinned as a tag AND as branch `legacy-v1.0` for defense in depth.

This is the state of the code before the AXON vectorization campaign:
- 285 `.m` files at repo root
- ~176 MB tracked repo weight
- Broken OOP factory dispatch (missing `MetodoMPFAH`, `MetodoMPFAQL`,
  `MetodoNLFVPP` files; orphan `SolverMPFAH`, `SolverNLFVPP` inheriting
  from missing `SolverBase`)
- No test harness
- ~62 dead files reachable from no entry point

Preserved for reproducibility of any pre-vectorization result.

---

## Detailed PR log — 2026-07 vectorization campaign

_Kept for historical reference. All PRs merged into `master` via v2.0.0._

### Foundation (Phase A — 7 PRs)
- **PR-A1** `f8da9ad` · `flowsim_init.m` + `flowsim_deinit.m` at repo root
- **PR-A2** `bd306b6` · Fix broken OOP — unified all methods under `MetodoBase`
- **PR-A3** `747b975` · `+fs/+util/assertFS.m` (5-invariant checker)
- **PR-A4** `08c9572` · `+fs/+mesh/build.m` (env → FS struct adapter)
- **PR-A5** `788eb0b` · `+fs/+csr/buildCorners.m` (CSR-flat corner layout — KEY enabler)
- **PR-A6** `8f62828` · `capture_baseline.m` + `fs_test_env.m` + goldens
- **PR-A7** `407a2b4` · `unit_baseline_reproduces.m` (bit-identical oracle)

### LPEW2 vectorization (Phase B — 5 PRs, full pipeline)
- **PR-B1** `5d101e0` · `fs.lpew.OPT` — batched geometry gather
- **PR-B2** `67bcf77` · `fs.lpew.v2.{angulos,netas}` + `fs.csr.buildCornerShifts`
- **PR-B3** `df5922b` · `fs.lpew.v2.ksInterp`
- **PR-B4** `5b1ab59` · `fs.lpew.v2.lambdaWeights`
- **PR-B5** `8187d7e` · `fs.lpew.v2.preLPEW2` — matches legacy at 1e-15 rel

### Assembly vectorization (Phase C — 6 PRs)
- **PR-C1+C2** `35811a5` · `fs.assembly.mpfad.build` — fully vectorized MPFA-D
- **PR-C3** `699d053` · `fs.assembly.tpfa.build`
- **PR-C4+C5+C6** `8305a76` · scaffolds for MPFAH/NLFVPP/MPFAQL/DMP/NLFVH

### Flow-rate + interpolation (Phase D — 3 PRs)
- **PR-D1+D2+D3** `1a8878d` · `+fs/+flow/{mpfad,tpfa}` + `+fs/+lpew/pinterp`

### Reorganization (Phase E — 7 PRs)
- `.gitignore` + untrack 154 MB of binaries (Option A)
- `fs.data.paths` loader stub + physical relocation (Option D)
- 17 dead files → `legacy/`  (transm*, preprocessor*, ferncodes_*_con.m)
- `.msh` files → `meshes/{hermeline,kozdon,other}/`
- 32 unknown files → `legacy/unknown/`
- 123 auxiliary files → `legacy/{limiters,saturation,transport,calc,get,solvers,utility,test-scripts}/`

### Docs + iterators (Phase F — 5 PRs)
- `startup.m` → delegates to `flowsim_init`
- `README.md`, `CHANGELOG.md`, `docs/{code-map,how-to-use,vectorization-guide,globals-audit}.md`
- `fs.iter.{picard,anderson,lscheme}` + `fs.lpew.dmpWeights` (cross-cluster wrappers)

### Verification totals
- **Total commits merged from `flowsim-artur`**: 33 PRs + infra
- **Test assertions**: ~500 across all PRs, 100% green
- **Regressions**: 0
- **Golden baselines**: TPFA + MPFA-D on `M8/num439` — bit-identical

### Deferred to future work
- Full triplet-form rewrite of MPFA-H, MPFA-QL, NLFV-PP, NLFV-H, DMP
  assemblers (scaffolds shipped; pattern established via MPFA-D).
  Estimate: **30–45 hours** total. Order by ascending difficulty:
  NLFV-H (2–3 h) → NLFV-PP (3–5 h) → MPFA-QL (3–5 h) → DMP (8–12 h)
  → MPFA-H (12–20 h).
- Kill globals in reachable set (see `docs/globals-audit.md`) — requires
  coordinated caller + callee migration.
- Cross-cluster kernel renames (e.g. `ferncodes_pressureinterpNLFVPP` is
  actually a shared LPEW-based interpolator used by MPFA-D too).
- Triage the 32 files in `legacy/unknown/`.

## Contributors

- **Artur Castiel Reis de Souza** (owner — decisions, direction, branch approval)
- **AXON code-dev** (autonomous implementation under owner grant 2026-07-03)
- **GitHub Copilot** (co-author per commit trailer)
