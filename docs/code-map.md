# Code map — FlowSim

_Deep-dive companion to `README.md`. Answers "where does X live?" and "who calls Y?"._

## By concern

### Mesh & connectivity (once per mesh load)
| Function | Location | Purpose |
|---|---|---|
| `preprocessormod` | `runtime/preproc/` | reads Start.dat, builds env.geometry from .msh |
| `preprocessmethod` | `runtime/preproc/` | per-method premethod setup |
| `ferncodes_elementface` | `legacy/ferncodes/shared/` | builds V, N, F element-face maps |
| `fs.mesh.build` | `+fs/+mesh/` | wraps env.geometry → FS.mesh + FS.geom |
| `fs.csr.buildCorners` | `+fs/+csr/` | CSR-flat node→corner layout |
| `fs.csr.buildCornerShifts` | `+fs/+csr/` | k-1 / k+1 wrap indices for LPEW |

### Boundary conditions (per-timestep)
| Function | Location | Purpose |
|---|---|---|
| `PLUG_bcfunction` | `runtime/plug/` | evaluates BC via benchmark |
| `PLUG_bcfunction_con` | `runtime/plug/` | concentration BC |
| `PLUG_kfunction` | `runtime/plug/` | evaluates permeability tensor (per-h) |
| `PLUG_sourcefunction` | `runtime/plug/` | source term |
| `PLUG_dfunction` | `runtime/plug/` | dispersion |
| `PLUG_Gfunction` | `runtime/plug/` | gravitational source |
| `ferncodes_calflag` | `legacy/ferncodes/shared/` | builds nflag + nflagface |

### Permeability & tensors
| Function | Location | Purpose |
|---|---|---|
| `ferncodes_Kde_Ded_Kt_Kn` | `legacy/ferncodes/mpfad/` | MPFA-D transmissibilities |
| `ferncodes_Kde_Ded_Kt_Kn_TPFA` | `legacy/ferncodes/mpfad/` | TPFA transmissibilities |
| `ferncodes_coefficient` | `legacy/ferncodes/shared/` | NLFV-PP coefficients |
| `ferncodes_coefficientmpfaH` | `legacy/ferncodes/mpfah/` | MPFA-H coefficients |
| `ferncodes_harmonicopoint` | `legacy/ferncodes/mpfah/` | harmonic points for MPFA-H |
| `ferncodes_weightnlfvDMP` | `legacy/ferncodes/dmp/` | DMP weights (used by MPFA-H, DMP) |

### LPEW2 (linearity-preserving interpolation weights)
| Function | Location | Purpose |
|---|---|---|
| `OPT_Interp_LPEW` | `legacy/ferncodes/lpew/` | per-node geometry gather (legacy) |
| `angulos_Interp_LPEW2` | `legacy/ferncodes/lpew/` | corner angles (legacy) |
| `netas_Interp_LPEW` | `legacy/ferncodes/lpew/` | netas ratios (legacy) |
| `Lamdas_Weights_LPEW2` | `legacy/ferncodes/lpew/` | lambda weights (legacy) |
| `ferncodes_Ks_Interp_LPEW2` | `legacy/ferncodes/lpew/` | permeability projections (legacy) |
| `ferncodes_Pre_LPEW_2_vect` | `legacy/ferncodes/lpew/` | ★ legacy driver (partial vect) |
| `fs.lpew.OPT` | `+fs/+lpew/` | batched geometry gather (vect) |
| `fs.lpew.v2.angulos` | `+fs/+lpew/+v2/` | batched angles (vect) |
| `fs.lpew.v2.netas` | `+fs/+lpew/+v2/` | batched netas (vect) |
| `fs.lpew.v2.ksInterp` | `+fs/+lpew/+v2/` | batched permeability projections (vect) |
| `fs.lpew.v2.lambdaWeights` | `+fs/+lpew/+v2/` | interior-vectorized lambda (vect) |
| **`fs.lpew.v2.preLPEW2`** | `+fs/+lpew/+v2/` | ★ full pipeline driver (vect) |

### Assembly (matrix construction)
| Function | Location | Purpose |
|---|---|---|
| `ferncodes_globalmatrix_MPFAD` | `legacy/ferncodes/mpfad/` | MPFA-D legacy assembler |
| `ferncodes_globalmatrix_TPFA` | `legacy/ferncodes/mpfad/` | TPFA legacy assembler |
| `ferncodes_globalmatrix` | `legacy/ferncodes/shared/` | steady/groundwater assembler |
| `ferncodes_assemblematrixMPFAH` | `legacy/ferncodes/mpfah/` | 820 L MPFA-H legacy assembler |
| `ferncodes_assemblematrixMPFAQL` | `legacy/ferncodes/mpfaql/` | MPFA-QL |
| `ferncodes_assemblematrixNLFVPP` | `legacy/ferncodes/nlfvpp/` | NLFV-PP |
| `ferncodes_assemblematrixNLFVH` | `legacy/ferncodes/nlfvh/` | NLFV-H |
| `ferncodes_assemblematrixDMP` | `legacy/ferncodes/dmp/` | DMP-preserving |
| **`fs.assembly.mpfad.build`** | `+fs/+assembly/+mpfad/` | ★ fully vectorized MPFA-D |
| `fs.assembly.tpfa.build` | `+fs/+assembly/+tpfa/` | vectorized TPFA |
| `fs.assembly.{mpfah,mpfaql,nlfvpp,nlfvh,dmp}.build` | `+fs/+assembly/+<x>/` | scaffolds (delegate to legacy) |

### Nonlinear iterators (called inside resolver)
| Function | Location | Purpose |
|---|---|---|
| `ferncodes_iterpicard` | `legacy/ferncodes/shared/` | Picard fixed-point |
| `ferncodes_iterpicardANLFVPP2` | `legacy/ferncodes/nlfvpp/` | Anderson accel (NLFV-PP shape) |
| `ferncodes_andersonacc` | `legacy/ferncodes/shared/` | generic Anderson |
| `ferncodes_andersonacc2` | `legacy/ferncodes/shared/` | Anderson variant |
| `L_scheme` | `legacy/unknown/` | L-scheme regularized iteration |

### Flow rate & pressure interpolation (post-solve)
| Function | Location | Purpose |
|---|---|---|
| `ferncodes_flowrate` | `legacy/ferncodes/mpfad/` | MPFA-D flow-rate (already vect) |
| `ferncodes_flowrateTPFA` | `legacy/ferncodes/mpfad/` | TPFA flow-rate |
| `ferncodes_flowratelfvHP` | `legacy/ferncodes/mpfah/` | MPFA-H flow-rate |
| `ferncodes_flowratelfvMPFAQL` | `legacy/ferncodes/mpfaql/` | MPFA-QL flow-rate |
| `ferncodes_pressureinterpNLFVPP` | `legacy/ferncodes/nlfvpp/` | nodal pressure interp (shared) |
| `ferncodes_pressureinterpMPFAQL` | `legacy/ferncodes/mpfaql/` | MPFA-QL nodal pressure |
| `ferncodes_pressureinterpHP` | `legacy/ferncodes/mpfah/` | MPFA-H nodal pressure |
| `fs.flow.mpfad` | `+fs/+flow/` | MPFA-D flow-rate (renamed) |
| `fs.flow.tpfa` | `+fs/+flow/` | TPFA flow-rate (renamed) |
| `fs.lpew.pinterp` | `+fs/+lpew/` | shared nodal pressure interp (renamed) |

### Time drivers
| Function | Location | Purpose |
|---|---|---|
| `main` | root | ENTRY POINT |
| `hydraulic` | `runtime/time/` | steady hydraulic solver |
| `hydraulic_RE` | `runtime/time/` | Richards transient solver |
| `IMPES` | `runtime/time/` | Implicit Pressure Explicit Saturation |
| `IMPEC` | `runtime/time/` | Implicit Pressure Explicit Concentration |
| `IMHEC` | `runtime/time/` | (variant) |
| `preRE` | `runtime/preproc/` | Richards preprocessor |
| `setmethod` | `runtime/time/` | time-driver dispatch |
| `soil_properties` | `runtime/util/` | Richards accumulation term (dtheta/dt) |
| `addsource` | `runtime/util/` | inserts wells into system |
| `solver` | `runtime/util/` | linear solve wrapper |
| `postprocessor` | `runtime/util/` | writes results |

### Test harness
| Function | Location | Purpose |
|---|---|---|
| `fs_setup / fs_expect / fs_teardown / fs_frob` | `tests/helpers/` | primitives |
| `fs_test_env` | `tests/helpers/` | WSL-safe Start.dat isolation |
| `capture_baseline` | `tests/helpers/` | golden baseline capture |
| `run_all` | `tests/` | aggregate runner |

## By call graph (Richards runtime path)

```
main
└─ hydraulic_RE (via setmethod, phasekey=6)
   └─ ferncodes_solver (per time step)
      ├─ env.metodo.montarSistema
      │  └─ ferncodes_globalmatrix_MPFAD  (or fs.assembly.mpfad.build)
      │     ├─ ferncodes_Kde_Ded_Kt_Kn (via env.premethod.MPFAD)
      │     └─ ferncodes_Pre_LPEW_2_vect (via env.premethod.MPFAD)
      │        (or fs.lpew.v2.preLPEW2 — vectorized twin)
      ├─ addsource
      └─ env.metodo.resolver
         └─ ferncodes_iterpicard (or L_scheme, or ANLFVPP2)
            └─ per Picard iter:
               ├─ PLUG_kfunction        → new kmap
               ├─ metodo.atualizarPremethod  → refresh Kde/Ded/Kt/Kn + LPEW2
               ├─ metodo.montarSistema  → rebuild M, I
               ├─ solver(M, I)          → new p
               └─ residual check
```

Per Picard iteration, the two hot spots are:
1. **LPEW2 weight recompute** — now vectorized (`fs.lpew.v2.preLPEW2`)
2. **Matrix assembly** — now vectorized (`fs.assembly.mpfad.build`)

## By dead-code inventory (see study `call-graph-and-reachability.md`)

- **~62 files** were confirmed dead by grep (no callers from any entry point)
- **~7000 LOC** of `transm*` code is dead → moved to `legacy/transm/`
- **7 `*_con.m` files** (concentration-coupled variants) are dead → `legacy/ferncodes-con/`
- **32 files** in the "unknown" bucket → `legacy/unknown/` pending owner triage
- **4 files** were data tables masquerading as .m (`parametros*.m`, `conduchidraulica.m`, `getchue.m`) — should move to `.mat` (owner decision pending)

## Naming discipline gotchas

- **`ferncodes_pressureinterpNLFVPP`** is NOT NLFVPP-specific — it's a shared
  LPEW-based nodal pressure interpolator used by MPFA-D too. Renamed to
  `fs.lpew.pinterp`.
- **`ferncodes_iterpicardANLFVPP2`** — Anderson acceleration wrapper, used
  by MPFA-D case `'AA'`. Also misnamed for its cross-method role.
- **`ferncodes_weightnlfvDMP`** — shared by MPFA-H AND DMP.

If you see a function name suggesting one method but it's called by others,
that's a legacy naming artefact. Cross-cluster kernels are true shared
infrastructure — renaming them belongs to future PR-F2.
