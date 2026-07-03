# Globals inventory + migration path

_PR-F1 artifact. Documents the top globals still in use in legacy code + their FS struct equivalents._

## Why not "just kill them"?

The `global` declarations in `legacy/ferncodes/**` cannot be removed piecemeal
without breaking the runtime — each `global` in a callee is paired with an
implicit expectation that the same global was populated by some caller (usually
`preprocessormod` at startup). Killing them requires **coordinated caller +
callee migration** in tandem, which is a per-file mechanical task better done
as its own campaign.

This document is the **migration map**: for each global still in use, the FS
struct field that replaces it and the migration pattern.

## Top-13 globals (by files-using count)

| Global | Files | Kind | FS struct equivalent | Notes |
|---|---:|---|---|---|
| `bedge`     | 116 | mesh face table | `FS.mesh.bedge` | boundary faces |
| `inedge`    | 107 | mesh face table | `FS.mesh.inedge` | interior faces |
| `coord`     | 107 | mesh coords    | `FS.mesh.coord` | node coordinates |
| `elem`      |  90 | mesh conn      | `FS.mesh.elem` | element vertex tables |
| `centelem`  |  90 | derived geom   | `FS.geom.centElem` | element centroids |
| `numcase`   |  62 | physics tag    | `FS.cfg.numcase` | benchmark selector |
| `bcflag`    |  38 | BC map         | `FS.bc.bcflag` | flag→value dispatch |
| `normals`   |  29 | derived geom   | `FS.geom.normalBnd` / `.normalInt` | face normals |
| `esurn2`    |  25 | CSR ptr        | `FS.mesh.esurn2` (alias `FS.csr.nodePtr`) | elements-around-node row-pointer |
| `nsurn2`    |  22 | CSR ptr        | `FS.mesh.nsurn2` | nodes-around-node row-pointer |
| `phasekey`  |  21 | physics tag    | `FS.cfg.phasekey` | 1/4/5/6 dispatch |
| `elemarea`  |  19 | derived geom   | `FS.geom.elemArea` | per-element areas |
| `bcflagc`   |  19 | BC map (conc)  | `FS.bc.bcflagc` | concentration BCs |

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

## What's already migrated

- **`+fs/+mesh/build.m`** — reads env.geometry.* + env.config.*, packs into FS
- **`+fs/+csr/**`** — builds CSR corners + shifts, all FS-based
- **`+fs/+lpew/**`** — reads FS.mesh.*, FS.geom.*, FS.csr.*, FS.perm.*
- **`+fs/+assembly/+mpfad/build.m`** — reads env-style directly (no globals)
- **`+fs/+assembly/+tpfa/build.m`** — same

The **modern OOP layer** (Metodo*/Sim*/Caso* classes) is already env-based —
they don't touch globals directly.

## What still uses globals

- **All `legacy/ferncodes/**` files** (66 files) — pending per-file migration
- **`hydraulic.m`** — one big driver still using globals
- **`IMPES` / `IMPEC` / `IMHEC`** — time-loop drivers
- **`PLUG_*.m`** — some read globals for numcase dispatch

## Priority for future migration

1. **Files on the Richards runtime path** (highest value):
   `ferncodes_iterpicard`, `ferncodes_calflag`, `ferncodes_Kde_Ded_Kt_Kn`
2. **Assembly variants when full-vectorization lands**:
   `ferncodes_assemblematrix{MPFAH,NLFVPP,MPFAQL,DMP,NLFVH}` (already
   scaffolded — see `+fs/+assembly/`)
3. **Post-processing**:
   `ferncodes_flowrate*`, `ferncodes_pressureinterp*`
4. **Everything else** — lowest priority, code path frequency unknown

## Not migrating

- **`transm*.m`** — all dead (see study `call-graph-and-reachability.md`),
  already in `legacy/transm/`. Delete if desired.
- **`legacy/unknown/**`** — 32 files pending owner triage; no migration until
  known-live status.
- **`ferncodes_*_con.m`** — dead concentration-coupled variants, in
  `legacy/ferncodes-con/`.

## Verification

To verify a migrated file, run its unit test (`tests/unit/unit_<name>.m`) which
diffs the new signature against the legacy signature on M8 mesh. Expected diff:
< 1e-12 relative Frobenius (matrices) or < 1e-10 L2 (vectors).
