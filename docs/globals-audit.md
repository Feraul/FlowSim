# Globals inventory + migration path

_Documents the top globals still in use in legacy code + their FS struct equivalents._
_Refreshed 2026-07-03 post-v2.0.1 (paths reflect the `runtime/` + `legacy/` layout)._

## Why not "just kill them"?

The `global` declarations in `legacy/ferncodes/**` cannot be removed piecemeal
without breaking the runtime — each `global` in a callee is paired with an
implicit expectation that the same global was populated by some caller (usually
`runtime/preproc/preprocessormod` at startup). Killing them requires
**coordinated caller + callee migration** in tandem, which is a per-file
mechanical task better done as its own campaign.

This document is the **migration map**: for each global still in use, the FS
struct field that replaces it and the migration pattern.

## Scale of the problem (post-v2.0.1)

- **~183 `.m` files** still declare at least one `global` (counted via
  `grep -rl "^\s*global\s" --include='*.m' .`)
- Almost all live under `legacy/ferncodes/**` (67 files) + `runtime/**` +
  a few remaining under `benchmarks/`, `factories/`, `simulacoes/`

## Top-13 globals (by files-using count — refreshed)

| Global | Files | Kind | FS struct equivalent | Notes |
|---|---:|---|---|---|
| `bedge`     | 107 | mesh face table | `FS.mesh.bedge` | boundary faces |
| `inedge`    | ~100 | mesh face table | `FS.mesh.inedge` | interior faces |
| `coord`     | ~100 | mesh coords    | `FS.mesh.coord` | node coordinates |
| `elem`      |  ~85 | mesh conn      | `FS.mesh.elem` | element vertex tables |
| `centelem`  |  ~85 | derived geom   | `FS.geom.centElem` | element centroids |
| `numcase`   |  58 | physics tag    | `FS.cfg.numcase` | benchmark selector |
| `bcflag`    |  ~35 | BC map         | `FS.bc.bcflag` | flag→value dispatch |
| `normals`   |  ~29 | derived geom   | `FS.geom.normalBnd` / `.normalInt` | face normals |
| `esurn2`    |  ~25 | CSR ptr        | `FS.mesh.esurn2` (alias `FS.csr.nodePtr`) | elements-around-node row-pointer |
| `nsurn2`    |  ~22 | CSR ptr        | `FS.mesh.nsurn2` | nodes-around-node row-pointer |
| `phasekey`  |  22 | physics tag    | `FS.cfg.phasekey` | 1/4/5/6 dispatch |
| `elemarea`  |  ~19 | derived geom   | `FS.geom.elemArea` | per-element areas |
| `bcflagc`   |  ~19 | BC map (conc)  | `FS.bc.bcflagc` | concentration BCs |

_(Counts are approximate — bedge is exact. See `grep -rl "^\s*global.*<name>"
--include='*.m' .` to refresh.)_

## Migration pattern (per-file mechanical)

For each legacy `ferncodes_*.m` file to be migrated:

**Before**:
```matlab
function [out] = ferncodes_foo(other, args)
    global bedge inedge coord elem numcase bcflag;
    % ... uses bedge(:,5), inedge(:,3), etc.
end
```

**After**:
```matlab
function [out] = ferncodes_foo(FS, other, args)   % FS added as first arg
    bedge      = FS.mesh.bedge;
    inedge     = FS.mesh.inedge;
    coord      = FS.mesh.coord;
    elem       = FS.mesh.elem;
    numcase    = FS.cfg.numcase;
    bcflag     = FS.bc.bcflag;
    % ... unchanged body
end
```

Every caller then passes `FS` as an extra arg. That's the tandem-migration
constraint: caller AND callee change together per PR.

## What's already migrated (v2.0.x)

- **`+fs/+mesh/build.m`** — reads `env.geometry.*` + `env.config.*`, packs into FS
- **`+fs/+csr/**`** — builds CSR corners + shifts, all FS-based
- **`+fs/+lpew/**`** — reads `FS.mesh.*`, `FS.geom.*`, `FS.csr.*`, `FS.perm.*`
- **`+fs/+assembly/+mpfad/build.m`** — reads env-style directly (no globals)
- **`+fs/+assembly/+tpfa/build.m`** — same
- **`+fs/+iter/{picard,anderson,lscheme}`** — env-based wrappers over legacy iterators
- **`+fs/+lpew/dmpWeights`** — env-based wrapper over `ferncodes_weightnlfvDMP`
- **`+fs/+flow/{mpfad,tpfa}`** — env-based

The **modern OOP layer** (`Metodo*` / `Sim*` / `Caso*` classes under `solvers/`,
`simulacoes/`, `benchmarks/`) is already env-based — they don't touch globals
directly. `MetodoBase` provides `env` as method context.

## What still uses globals

- **All `legacy/ferncodes/**` files** (~67 `.m` files) — pending per-file migration
- **`runtime/time/hydraulic.m`** — the big steady driver, still using globals
- **`runtime/time/hydraulic_RE.m`** — Richards transient driver
- **`runtime/time/IMPES.m` / `IMPEC.m` / `IMHEC.m`** — time-loop drivers
- **`runtime/plug/PLUG_*.m`** — some read `numcase`/`phasekey` for benchmark dispatch
- **`runtime/util/solver.m`** / `solvePressure*.m` — solver wrappers
- **`benchmarks/Caso439.m`** — reads `centelem`, `elem` from globals

## Priority for future migration

1. **Files on the Richards runtime path** (highest value):
   `runtime/time/hydraulic_RE`, `ferncodes_iterpicard`, `ferncodes_calflag`,
   `ferncodes_Kde_Ded_Kt_Kn`
2. **Assembly variants when full-vectorization lands**:
   `ferncodes_assemblematrix{MPFAH,NLFVPP,MPFAQL,DMP,NLFVH}` — already
   scaffolded under `+fs/+assembly/`, but their internal ferncodes callees
   still use globals
3. **Post-processing**:
   `ferncodes_flowrate*`, `ferncodes_pressureinterp*`
4. **PLUG_ callbacks** — smallest surface, could be batched
5. **Everything else** — lowest priority, code path frequency unknown

## Not migrating

- **`legacy/transm/**`** — 8 files, all dead. Delete if desired.
- **`legacy/preprocessor/`** — 2 files, dead. Superseded by `runtime/preproc/preprocessormod`.
- **`legacy/ferncodes-con/`** — 7 dead concentration-coupled variants.
- **`legacy/solvers/`** — 2 dead solver variants.
- **`legacy/unknown/**`** — 32 files pending owner triage; no migration until
  known-live status.

## Verification

To verify a migrated file, run its unit test (`tests/unit/unit_<name>.m`) which
diffs the new signature against the legacy signature on M8 mesh. Expected diff:

- **Matrices** — `< 1e-12` relative Frobenius
- **Vectors** — `< 1e-10` relative L2 norm
- **The oracle** — `unit_baseline_reproduces` must still hit rel diff
  `0.000e+00` on the committed golden (M8/num439/mpfad and M8/num439/tpfa).

## Estimated cost

A **single-file** migration is 15-30 min of mechanical work (declare param,
substitute globals, update all callers). But because migrations must be tandem
(caller + callee), each cascade can pull in 3-10 files. Realistic estimate for
killing 90% of globals in the reachable set: **40-60 hours** as a dedicated
follow-up campaign, likely best sequenced right after the remaining 5-method
vectorization work (which itself already touches the reachable assembly path).
