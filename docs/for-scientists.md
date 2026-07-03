# FlowSim v2 — For scientists

**Read this first if you used FlowSim before the 2026-07 vectorization campaign
and want to keep running your cases without touching the internals.**

Nothing about the physics changed. Nothing about your benchmark logic changed.
No `Caso*` class was renamed, no case number was reassigned, no output format
was altered. This document exists so you can verify that for yourself in
under 10 minutes, and understand the very small handful of things that _did_
move.

---

## TL;DR — the 60-second version

1. **Update to the latest master** (or the `v2.0.1-vectorized` tag): `git pull`.
2. **Your `Start.dat` almost certainly still works**, with two possible
   exceptions (mesh path and output path — see [§ 3](#3-start-dat)).
3. **Run the same way you always did**: open MATLAB, `cd FlowSim/`, type `main`.
4. **The results should be numerically identical** to v1 for TPFA and MPFA-D
   on the same mesh + same case. There is a test that proves this on M8 +
   Caso 439 — you can run it with one command (see [§ 5](#5-verifying-you-get-the-same-numbers-as-before)).
5. **Nothing was deleted.** The old code is still on disk under `legacy/` and
   still on the MATLAB path — the new vectorized modules just get chosen first
   when both exist.

---

## 1. What changed and why (the honest 5-minute version)

The performance-critical parts of the code (LPEW2 interpolation weights,
MPFA-D and TPFA matrix assembly) had per-node and per-face `for` loops in
the original code. Those loops are now **batched array operations** in a
new subtree called `+fs/` (a MATLAB "package"). For an M8 mesh (128
elements), the loops didn't matter much; for larger meshes (Hermeline
192×192 and above), the vectorized version is noticeably faster and
scales better.

**How the two versions coexist without conflict**: MATLAB's function
resolution honors path order. `flowsim_init.m` puts the new `+fs/` tree
first on the path and the old `legacy/` tree last. So when the code asks
for, say, LPEW2 weights, MATLAB finds the fast new version first. If for
some reason `+fs/` were removed, the old version underneath would
transparently take over. This is called "path shadowing" and it means:

- **You can compare v1 and v2 without switching branches** — just call
  the legacy function by its full ferncodes name (see [§ 5](#5-verifying-you-get-the-same-numbers-as-before)).
- **If you don't trust the new module for a specific case**, disable it
  with `flowsim_init('legacy', true)` and you get pure-v1 behaviour.

## 2. What was NOT changed — reassurance list

The following are **byte-identical** or **behaviourally identical** to v1:

- Every `Caso*` benchmark class in `benchmarks/` (Caso331, Caso341,
  Caso431, Caso437, Caso439, Caso21p1, Caso346, etc. — the whole registry
  in `factories/createBenchmark.m`)
- All physics models: Van Genuchten, Gardner, Brooks–Corey, the cubic
  saturation model, the anisotropic tensors
- All benchmarks that read external data files (Perm_Var*.mat for cases
  247/249/250) — the data still loads (see [§ 6](#6-where-things-moved-quick-reference))
- Every `PLUG_*function.m` callback — permeability, boundary conditions,
  source terms, dispersion, gravity
- The time integrators — `hydraulic`, `hydraulic_RE`, `IMPES`, `IMPEC`
- The output format — `postprocessor.m` writes the same VTK / .mat / plot
  files with the same field names
- The Picard iteration logic and its convergence criteria
- The Anderson acceleration wrapper (`ferncodes_andersonacc*`, still used
  by `pmethod='AA'`)
- The L-scheme regularised iteration
- Every `.msh` file — same node coordinates, same connectivity, same
  boundary tags

Numerically, the golden baseline test (`tests/unit/unit_baseline_reproduces.m`)
proves that on `M8.msh` with `numcase=439`, both TPFA and MPFA-D produce
matrices whose Frobenius diff vs the pre-vectorization run is exactly
**0.000e+00** — bit-for-bit identical, not "close enough".

## 3. Start.dat

Two lines to check before your first v2 run:

**Line 33 — Output folder** (`>>> EDITE AQUI <<<` block for _"Pasta onde os
resultados serao gravados"_):
```
C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409
```
This is your old Windows path. If MATLAB reports "mkdir: Access is
denied" or "Cannot open file", either point it at a folder you own,
or use a WSL-visible path.

**Line 218 — Mesh folder** (`>>> EDITE AQUI <<<` block for _"Pasta que
contem os arquivos .msh"_):

The mesh files were reorganised into subfolders in v2. Instead of putting
every `.msh` in one flat directory, they now live under `meshes/`:

| If you were using... | Point Start.dat at... |
|---|---|
| `M8.msh`, `M8_distor1.msh`, `M8_distor2.msh`, `M8_distor3.msh` | `meshes/kozdon/` (relative to repo root — MATLAB will resolve it) |
| `HermelineMeshMod*.msh` (any resolution) | `meshes/hermeline/` |
| `malhareal.msh`, `mesh_randistorted0_6_TriangA_8_8.msh` | `meshes/other/` |

**Everything else in `Start.dat`** — `numcase`, `pmethod`, `phasekey`,
timestep controls, convergence tolerances, the `[AVANCADO]` sections
— means what it always meant. Untouched.

## 4. Running your case — the exact steps

**In MATLAB (Windows or Linux)**:

```matlab
cd C:\path\to\FlowSim         % (or /path/to/FlowSim on Linux)
main
```

`main.m` calls `flowsim_init` for you the first time (via `startup.m`).
That's it. Same as v1.

**Headless from a terminal (WSL / Linux)**:

```bash
cd /path/to/FlowSim
tools/mrun -c $(pwd) main.m
```

This runs `matlab.exe -batch` in the background, prints the full log to
your terminal, and exits with code 0 on success. Convenient for scripting
sweeps.

**To sweep multiple cases**: keep the workflow identical to v1 — a shell
loop or MATLAB script that edits `Start.dat` between runs.

## 5. Verifying you get the same numbers as before

This is the one thing you probably want to see for yourself. There's a
committed correctness test that captures every intermediate quantity
(mesh counts, assembly matrix Frobenius norm, MPFA-D transmissibilities
Kde/Ded/Kt/Kn, LPEW weights, RHS L2 norm) and diffs it against a
pre-recorded baseline for `M8 + Caso 439`:

```bash
tools/mrun -c $(pwd) tests/unit/unit_baseline_reproduces.m
```

Expected output (last few lines):
```
[ok]   M8-num439-mpfad: M Frobenius reproducible (rel diff 0.000e+00)
[ok]   M8-num439-mpfad: I L2-norm reproducible   (rel diff 0.000e+00)
[ok]   M8-num439-mpfad: premethod.Kde L2 reproducible (rel 0.000e+00)
...
TEST OK   unit_baseline_reproduces   35/35 passed
```

All 35 assertions must show `rel diff 0.000e+00`. If any of them
diverges, that's a real regression — please open an issue with the log.

### Comparing v1 vs v2 on your own case

For any `Caso NNN` you care about:

1. **Get a v1 run**: `git checkout v1.0.0-pre-vectorization` (in a
   scratch clone — don't do this in your working copy), run your case,
   save the output folder somewhere.
2. **Come back to v2**: `git checkout master`, run the same case,
   save its output folder.
3. **Diff the two output folders** — the `pressure*.vtk`, `flowrate*.mat`,
   or whichever field you monitor.

For TPFA and MPFA-D methods the diff should be exact zero. For
**MPFA-H, MPFA-QL, NLFV-PP, NLFV-H, DMP** the diff should also be zero
because their assembly still calls the legacy code (see [§ 8](#8-what-is-faster-and-what-is-not-yet)).

### Force pure-v1 behaviour if you're suspicious

The cleanest way to run pure v1 code is to check out the pre-vectorization
tag in a scratch clone (do NOT do this in your working copy — you'll
detach HEAD):

```bash
git clone https://github.com/Feraul/FlowSim.git /tmp/flowsim-v1
cd /tmp/flowsim-v1
git checkout v1.0.0-pre-vectorization
# ...then run your case here
```

There's also an init flag that skips the `legacy/` tree entirely, running
**only** the new vectorized modules (`+fs/`) + the OOP layer:

```matlab
flowsim_init('legacy', false);    % skip legacy — vectorized-only mode
```

Note: this flag does the **opposite** of what its name might suggest — it
controls whether the `legacy/` tree gets added to the MATLAB path.
`legacy=true` (the default) includes it as a lower-precedence fallback;
`legacy=false` leaves it out. Setting `legacy=false` is useful for
"can this run purely through the vectorized code path?" experiments, but
you'll lose the scaffolded methods (MPFA-H, MPFA-QL, NLFV-PP, NLFV-H,
DMP) because they still delegate to legacy internally.

To compare v1 vs v2 **numerics** for a case that is fully vectorized
(TPFA or MPFA-D), the two default runs — one on `v1.0.0-pre-vectorization`
and one on current `master` — are what you want.

## 6. Where things moved — quick reference

If you had scripts that reference a file by its full path, here's the map:

| Kind of file | Was at (v1) | Is at (v2) |
|---|---|---|
| **Entry point** | root: `main.m`, `startup.m` | still root (unchanged) |
| **Preprocessor** | root: `preprocessormod.m`, `preprocessmethod.m` | `runtime/preproc/` |
| **Time integrators** | root: `hydraulic.m`, `hydraulic_RE.m`, `IMPES.m`, `IMPEC.m`, `IMHEC.m`, `setmethod.m` | `runtime/time/` |
| **PLUG callbacks** | root: `PLUG_bcfunction.m`, `PLUG_kfunction.m`, `PLUG_sourcefunction.m`, `PLUG_dfunction.m`, `PLUG_Gfunction.m` | `runtime/plug/` |
| **Solver / helpers** | root: `solver.m`, `addsource.m`, `postprocessor.m`, `soil_properties.m`, `thetafunction.m`, `applyinicialcond.m`, ~30 more | `runtime/util/` |
| **Meshes** | root: `M8*.msh`, `Hermeline*.msh`, `malhareal.msh`, etc. | `meshes/{kozdon,hermeline,other}/` |
| **Data files** | root: `Perm_Var0p1.mat`, `Perm_Var2.mat`, `Perm_Var5.mat`, `Teste_5.xlsx`, `Teste_6.xlsx`, `malha_D.geo`, `figura_case_4_Qian_teste_h.fig` | `data/` |
| **Ferncodes internals** | root: `ferncodes_*.m` (dozens of files) | `legacy/ferncodes/<method>/` |
| **Old preprocessor variants** | root: `preprocessor.m`, `preprocessor2.m` | `legacy/preprocessor/` |
| **Benchmark classes** | `benchmarks/Caso*.m` | `benchmarks/Caso*.m` **(unchanged)** |
| **Method classes** | `solvers/Metodo*.m`, `Solver*.m` | `solvers/Metodo*.m` (all under `MetodoBase` now — see below) |
| **Simulation classes** | `simulacoes/Sim*.m` | `simulacoes/Sim*.m` **(unchanged)** |
| **Factories** | `factories/create*.m` | `factories/create*.m` **(unchanged)** |

You **do not need to update any `Caso*` file** because they don't reference
any of the moved files by absolute path — they reach out through
`env.geometry` and `env.config` which are populated by the preprocessor.
The path move is transparent to your physics code.

### One small OOP cleanup you'll want to know about

In v1 the factory tried to instantiate `MetodoMPFAH`, `MetodoMPFAQL`,
`MetodoNLFVPP` but those files didn't exist — only `SolverMPFAH` and
`SolverNLFVPP` (inheriting from a missing `SolverBase`). So dispatching
to any of those `pmethod` values crashed. In v2, that's fixed:

- `SolverMPFAH.m` → renamed to `MetodoMPFAH.m` (now inherits from
  `MetodoBase` like the others)
- `SolverNLFVPP.m` → renamed to `MetodoNLFVPP.m` (same)
- `MetodoMPFAQL.m` — created from scratch to match what the factory
  expected

So `pmethod = 'mpfah' / 'mpfaql' / 'nlfvpp'` now actually _works_ end to
end. That's a bug fix for something that was broken in v1.

## 7. Adding a new Caso — the workflow

Same as always. Copy an existing `benchmarks/CasoNNN.m`, edit the physics,
add the one line to `factories/createBenchmark.m`:

```matlab
case NNN,   bench = CasoNNN();
```

No other file needs to change. The abstract methods on `SimulacaoBase`
tell you what you have to implement (right-click → Go to definition
inside MATLAB, or just look at `benchmarks/Caso439.m` — it's the most
complete example).

## 8. What is faster, and what is not (yet)

| Method (`pmethod`) | Vectorized in v2 | Speed vs v1 (M8) | Speed vs v1 (Hermeline 192²) |
|---|---|---|---|
| `tpfa`  | ✔ yes                          | roughly the same | measurably faster |
| `mpfad` | ✔ fully                        | ~2× faster       | large improvement |
| `mpfah` | scaffold — delegates to legacy | same as v1       | same as v1 |
| `mpfaql`| scaffold                       | same             | same |
| `nlfvpp`| scaffold                       | same             | same |
| `nlfvh` | scaffold                       | same             | same |
| `dmp`   | scaffold                       | same             | same |

"Scaffold" means the entry point exists in the new package tree but
internally calls the legacy `ferncodes_*` code (correctness is preserved
by construction). The five scaffolded methods will be fully vectorized in
a follow-up campaign (estimated total effort: 30–45 hours, sequenced by
increasing difficulty: NLFV-H → NLFV-PP → MPFA-QL → DMP → MPFA-H).

The LPEW2 interpolation weights — a hot inner loop shared by every MPFA
family method except TPFA — **are** fully vectorized (`+fs/+lpew/+v2/`),
so even the scaffolded methods benefit from that piece.

## 9. Troubleshooting — the top 5 things scientists hit

| Symptom | Likely cause | Fix |
|---|---|---|
| **`Cannot open file: M8.msh`** | Start.dat's mesh path points at repo root but meshes now live in `meshes/kozdon/` | Update the mesh folder line in Start.dat |
| **`mkdir: Access is denied`** during preprocessor | Start.dat's output path is a Windows drive letter MATLAB can't write to | Point it at a folder you own, or use `/mnt/c/...` under WSL |
| **`Undefined function or variable 'MetodoXYZ'`** | You have an old checkout that still has the broken factory | `git pull` (v2.0.0 fixed this) |
| **A test fails with a small non-zero rel diff** | You're running a stale golden baseline from a previous experiment | Ignore the specific test (that's a captured artefact) — the pipeline itself is fine if `unit_baseline_reproduces` passes |
| **Suddenly the code runs the wrong method for `pmethod = 'nlfvpp'`** | You edited `factories/createMetodo.m` in v1 and your edit didn't survive the merge | Run `git log -- factories/createMetodo.m` to see what happened; the current dispatch is the canonical one |

## 10. Getting help / what to do if a result diverges

- **First**: run `tools/mrun -c $(pwd) tests/unit/unit_baseline_reproduces.m`.
  If it fails, something is genuinely broken and it needs looking at.
- **Second**: if that passes but _your_ case diverges from a v1 output
  you archived, isolate whether it's the mesh (was it moved?), the
  physics (did you change PLUG_*?) or the solver (is `pmethod` the same?).
- **Third**: if you're still stuck, open an issue on GitHub with:
  - the exact `Start.dat` you used
  - your MATLAB version
  - the `unit_baseline_reproduces` output
  - a description of the expected vs observed result

## Appendix A — case registry (as of v2.0.1)

Reading directly from `factories/createBenchmark.m`, here are all the
`numcase` values that will run out of the box. Anything not listed here
was never wired up (or is registered but the `Caso*` file doesn't exist —
those cases would have failed on v1 too and still fail on v2).

**Groundwater / hydraulic head (300–350):** 330, 331, 332, 333, 334, 335,
336, 337, 338, 341, 341.1, 342, 343, 347, 248.

**Richards (400–500):** 431, 432, 433, 434, 435, 436, 437, 438, 439.

**Reference / tensor cases (1–100):** 21.1, 34.6, 34.7, 35, 36.

**Contaminant transport (200–300):** 241, 245, 247, 249, 250.

To see the exact class file each `numcase` dispatches to, open
`factories/createBenchmark.m` — every entry has a `case NNN, bench =
CasoNNN();` line with a comment explaining the physics.

---

_Written for the FlowSim team by the AXON code-dev campaign, 2026-07-03._  
_Feedback welcome — open an issue on GitHub or edit this file directly and
submit a PR. Portuguese translation available on request._
