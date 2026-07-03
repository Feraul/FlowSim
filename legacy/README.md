# legacy/ — Pre-vectorization code

Everything under `legacy/` is code that predates the 2026-07 vectorization
campaign (v2.0.0). It is **still on the MATLAB path** — added by
`flowsim_init.m` via `genpath` at the lowest precedence — so any caller
that hasn't migrated to `+fs/` continues to work transparently.

For each legacy symbol that has a vectorized replacement, MATLAB resolves
to the `+fs/…` version because `+fs/` is added first (higher precedence).
This is the **shadow** pattern — legacy stays as the correctness reference,
vectorized code takes over silently.

> 🇧🇷 **Versao em portugues**: [`legacy/LEIAME.md`](LEIAME.md)

---

## Cluster index

| Cluster | Files (.m) | What it is | Status | Shadowed by |
|---|---:|---|---|---|
| `ferncodes/` | 72 | The methods themselves — MPFA-D, MPFA-H, MPFA-QL, NLFV-PP, NLFV-H, DMP, LPEW1/2, shared kernels. Organized by method. | **live** | partially (`+fs/+lpew/+v2/`, `+fs/+assembly/{+mpfad,+tpfa}`) |
| `ferncodes-con/` | 7 | Concentration-coupled variants (`_con.m`) — dead code in this branch | **dead** | none |
| `preprocessor/` | 2 | `preprocessor.m` + `preprocessor2.m` — superseded by `runtime/preproc/preprocessormod.m` | **dead** | `runtime/preproc/preprocessormod` |
| `transm/` | 8 | Transmissibility computations from an older attempt (~7 kLOC of unreached code) | **dead** | none |
| `limiters/` | 5 | Flux/slope limiters — used only by the transport family | **live** (opt-in) | none |
| `saturation/` | 9 | IMPES-side saturation update helpers | **live** | none |
| `transport/` | 11 | Advection/dispersion transport, hyperbolic sub-solvers | **live** (opt-in) | none |
| `calc/` | 28 | `calc*` helpers — mass fluxes, viscosities, spectral fluxes | **partial** | some unreached |
| `get/` | 31 | `get*` helpers — mass extraction, boundary states, analytical solutions | **partial** | some unreached |
| `solvers/` | 2 | Old `solveEnriched` + `solvefluxadvecdispersive` — coupled solvers | **dead** | `runtime/util/solver*` |
| `utility/` | 25 | Grab-bag utilities (`Durlofsk`, `L_scheme`, `MalhaKozdon`, `SPEField`, etc.) | **partial** | `+fs/+iter/{picard,anderson,lscheme}` wraps `L_scheme` |
| `test-scripts/` | 6 | Pre-existing test/plot scripts (`Teste02.m`, `buckey_levalidation.m`, `plot_cilamce2023.m`, `testeBuckley.m`) | **live** (opt-in) | none |
| `unknown/` | 32 | Files that couldn't be classified during triage; owner disposition pending | **pending** | none |
| **Total** | **~238** | | | |

## `ferncodes/` subclusters (the meat of legacy)

`ferncodes/` — the pre-existing "Fernando's codes" tree — is organized by
discretization method to make it clear which files any given method needs:

| Subdir | Files | Method | Vectorized twin |
|---|---:|---|---|
| `ferncodes/shared/` | 18 | cross-method kernels (`ferncodes_globalmatrix`, `ferncodes_calflag`, `ferncodes_elementface`, `ferncodes_andersonacc*`, …) | none — most are already vectorized |
| `ferncodes/lpew/` | 10 | LPEW1 + LPEW2 (`OPT_Interp_LPEW`, `angulos_Interp_LPEW2`, `netas_Interp_LPEW`, `Lamdas_Weights_LPEW2`, `ferncodes_Ks_Interp_LPEW2`, `ferncodes_Pre_LPEW_2_vect`, …) | ★ **`+fs/+lpew/+v2/`** — full pipeline, 1e-15 rel |
| `ferncodes/mpfad/` | 10 | MPFA-D (`ferncodes_globalmatrix_MPFAD`, `ferncodes_Kde_Ded_Kt_Kn`, `ferncodes_flowrate`, …) | ★ **`+fs/+assembly/+mpfad/build`** + `+fs/+flow/mpfad` |
| `ferncodes/mpfah/` | 8 | MPFA-H (`ferncodes_assemblematrixMPFAH` — 820 lines, `ferncodes_coefficientmpfaH`, `ferncodes_flowratelfvHP`, `ferncodes_harmonicopoint`, …) | scaffold only (`+fs/+assembly/+mpfah/build` delegates) |
| `ferncodes/mpfaql/` | 6 | MPFA-QL (`ferncodes_assemblematrixMPFAQL`, `ferncodes_pressureinterpMPFAQL`, `ferncodes_flowratelfvMPFAQL`, …) | scaffold only |
| `ferncodes/nlfvpp/` | 8 | NLFV-PP (`ferncodes_assemblematrixNLFVPP`, `ferncodes_pressureinterpNLFVPP`, `ferncodes_iterpicardANLFVPP2`, `ferncodes_coefficient`, …) | scaffold only |
| `ferncodes/nlfvh/` | 6 | NLFV-H (`ferncodes_assemblematrixNLFVH`, …) | scaffold only |
| `ferncodes/dmp/` | 6 | DMP-preserving (`ferncodes_assemblematrixDMP`, `ferncodes_weightnlfvDMP`) | scaffold only |

## Retirement path

- **dead** clusters — safe to physically delete once someone confirms no
  external tooling depends on them. Current best-guess candidates for
  removal: `ferncodes-con/`, `preprocessor/`, `transm/`, `solvers/`.
- **live** clusters — cannot be deleted; retire only after the shadowing
  `+fs/` module lands + all callers migrate.
- **partial** clusters — mixed. Individual files can be deleted once
  reachability from `main.m` → factory → benchmark is confirmed empty.
  See `docs/code-map.md` § "By dead-code inventory".
- **pending** cluster (`unknown/`) — needs owner disposition. Each file
  is either dead (delete) or missed during initial triage (recategorize).

## Cross-cutting naming gotchas

Some legacy names suggest a single method but the file is actually shared
infrastructure. These are candidates for future renames (deferred as PR-F2):

| Misleading name | Reality |
|---|---|
| `ferncodes_pressureinterpNLFVPP` | LPEW-based nodal-pressure interp used by MPFA-D too. Vectorized as `+fs/+lpew/pinterp`. |
| `ferncodes_iterpicardANLFVPP2` | Anderson-acceleration wrapper used by MPFA-D case `'AA'`. |
| `ferncodes_weightnlfvDMP` | Shared by MPFA-H **and** DMP. |

## See also

- `../runtime/README.md` — active-runtime tree (what shadows what)
- `../docs/code-map.md` — function-level "where does X live?"
- `../docs/vectorization-guide.md` — how to add a new `+fs/` shadow
- `../CHANGELOG.md` — full campaign history
