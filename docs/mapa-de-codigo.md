# Mapa do codigo — FlowSim

_Companheiro deep-dive do `LEIAME.md`. Responde "onde vive X?" e "quem chama Y?"._

> **English version**: [`docs/code-map.md`](code-map.md)

## Arquitetura de uma olhada

O pipeline ponta-a-ponta, mostrando onde cada peca vive no disco e qual
camada substitui qual:

```
                                ┌──────────────────────────────────────┐
   USUARIO                      │  Start.dat  ← configure numcase,    │
    │                           │              pmethod, phasekey,     │
    │                           │              path da malha, arquivo │
    │                           └──────────────────────────────────────┘
    │
    │                     ╔═══════════════════════════════════════════╗
    ▼                     ║              (raiz)                       ║
  ┌─────────┐             ║  main.m ─▶ flowsim_init ─▶ setup do path: ║
  │ main.m  │──────────┐  ║    (+fs/) > (base|solvers|factories|      ║
  └─────────┘          │  ║     simulacoes|benchmarks) > (runtime/**) ║
                       │  ║     > (legacy/**)                          ║
                       │  ╚═══════════════════════════════════════════╝
                       │
                       ▼
      ┌───────────────────────────────────────────────────────────────┐
      │  runtime/preproc/preprocessormod                               │
      │    • le Start.dat                                              │
      │    • parseia .msh (de meshes/{hermeline,kozdon,other}/)        │
      │    • constroi env.geometry.{coord,elem,bedge,inedge,esurn,nsurn}│
      │    • constroi env.config.{numcase,pmethod,phasekey,perm,bcflag}│
      └───────────────────────────────────────────────────────────────┘
                       │
       ┌───────────────┼────────────────┬──────────────────┐
       ▼               ▼                ▼                  ▼
    factory:       factory:         factory:         (benchmark
    createBench    createMetodo     createSimulacao   initParms)
    ─────────      ─────────        ─────────
       │               │                │
       ▼               ▼                ▼
   benchmarks/     solvers/         simulacoes/
   Caso439.m       Metodo{TPFA,     Sim{Groundwater,
   (o unico        MPFAD,MPFAH,     Richards}
    totalmente     MPFAQL,NLFVPP}    (SimulacaoBase)
    impl)          (MetodoBase)
                       │
                       ▼
      ┌───────────────────────────────────────────────────────────────┐
      │  runtime/preproc/preprocessmethod   ▲ struct premethod por-metodo│
      │    ├── ferncodes_elementface  (legacy/ferncodes/shared/)      │
      │    ├── ferncodes_Kde_Ded_Kt_Kn   (legacy/ferncodes/mpfad/)    │
      │    └── ferncodes_Pre_LPEW_2_vect (legacy/ferncodes/lpew/)     │
      │           │                                                    │
      │           └──▶ ★ SOMBREADO por fs.lpew.v2.preLPEW2 (+fs/+lpew/)│
      └───────────────────────────────────────────────────────────────┘
                       │
                       ▼
      ┌───────────────────────────────────────────────────────────────┐
      │  runtime/time/setmethod   dispatch por phasekey               │
      │    ├── phasekey=1 → IMPES              (monofasico)           │
      │    ├── phasekey=4 → hydraulic          (groundwater estacionario)│
      │    ├── phasekey=5 → IMPEC              (concentracao)         │
      │    └── phasekey=6 → hydraulic_RE       (Richards transiente)  │
      └───────────────────────────────────────────────────────────────┘
                       │
                       ▼
      ┌───────────────────────────────────────────────────────────────┐
      │  LOOP DE TEMPO  (por passo, por iteracao Picard)              │
      │    ┌─────────────────────────────────────────────────────┐   │
      │    │  runtime/plug/PLUG_kfunction    → novo kmap         │   │
      │    │  metodo.atualizarPremethod      → refresh Kde/etc   │   │
      │    │                                                      │   │
      │    │  metodo.montarSistema           → montar [M, I]:    │   │
      │    │    ┌──────────────────────────────────────────────┐│   │
      │    │    │ ★ NOVO (v2.0): +fs/+assembly/<metodo>/build  ││   │
      │    │    │   ├── mpfad: VETORIZADO TOTAL (bit-identico) ││   │
      │    │    │   ├── tpfa:  vetorizado                       ││   │
      │    │    │   └── {mpfah,mpfaql,nlfvpp,nlfvh,dmp}:        ││   │
      │    │    │        scaffold → delega pro legado           ││   │
      │    │    │                                                ││   │
      │    │    │ ANTIGO: ferncodes_globalmatrix_<METODO>       ││   │
      │    │    │      (legacy/ferncodes/<metodo>/)             ││   │
      │    │    └──────────────────────────────────────────────┘│   │
      │    │                                                      │   │
      │    │  addsource                        (runtime/util/)   │   │
      │    │  metodo.resolver                                     │   │
      │    │    └── fs.iter.{picard,anderson,lscheme}            │   │
      │    │        (+fs/+iter/) — wrapper sobre iters legadas    │   │
      │    │                                                      │   │
      │    │  metodo.calcularFlowrate         (+fs/+flow/ ou     │   │
      │    │                                    legacy/ferncodes)│   │
      │    │  runtime/util/postprocessor      → saida VTK/mat    │   │
      │    └─────────────────────────────────────────────────────┘   │
      └───────────────────────────────────────────────────────────────┘
```

### Shadow de path em duas camadas (o padrao chave de correcao)

```
    +fs/…      ← maior precedencia  (modulos vetorizados)
    ─────
    base/, solvers/, factories/, simulacoes/, benchmarks/
    ─────
    runtime/{preproc,time,plug,util}/
    ─────
    legacy/…   ← menor precedencia   (referencia de correcao legada)
```

Quando um simbolo existe em ambos, o MATLAB resolve para `+fs/`. O legado
fica como oraculo de correcao — todo modulo `+fs/` tem um
`tests/unit/unit_*.m` que faz o diff contra o gemeo legado com
tolerancia relativa `< 1e-12`.

---

## Por preocupacao

### Malha & conectividade (uma vez por carregamento)
| Funcao | Localizacao | Proposito |
|---|---|---|
| `preprocessormod` | `runtime/preproc/` | le Start.dat, constroi env.geometry a partir do .msh |
| `preprocessmethod` | `runtime/preproc/` | setup do premethod por-metodo |
| `ferncodes_elementface` | `legacy/ferncodes/shared/` | constroi mapas V, N, F de element-face |
| `fs.mesh.build` | `+fs/+mesh/` | wrappea env.geometry → FS.mesh + FS.geom |
| `fs.csr.buildCorners` | `+fs/+csr/` | layout CSR-flat de no→canto |
| `fs.csr.buildCornerShifts` | `+fs/+csr/` | indices k-1 / k+1 pro LPEW |

### Condicoes de contorno (por passo)
| Funcao | Localizacao | Proposito |
|---|---|---|
| `PLUG_bcfunction` | `runtime/plug/` | avalia BC pelo benchmark |
| `PLUG_bcfunction_con` | `runtime/plug/` | BC de concentracao |
| `PLUG_kfunction` | `runtime/plug/` | avalia tensor de permeabilidade (por-h) |
| `PLUG_sourcefunction` | `runtime/plug/` | termo-fonte |
| `PLUG_dfunction` | `runtime/plug/` | dispersao |
| `PLUG_Gfunction` | `runtime/plug/` | fonte gravitacional |
| `ferncodes_calflag` | `legacy/ferncodes/shared/` | constroi nflag + nflagface |

### Permeabilidade & tensores
| Funcao | Localizacao | Proposito |
|---|---|---|
| `ferncodes_Kde_Ded_Kt_Kn` | `legacy/ferncodes/mpfad/` | transmissibilidades MPFA-D |
| `ferncodes_Kde_Ded_Kt_Kn_TPFA` | `legacy/ferncodes/mpfad/` | transmissibilidades TPFA |
| `ferncodes_coefficient` | `legacy/ferncodes/shared/` | coeficientes NLFV-PP |
| `ferncodes_coefficientmpfaH` | `legacy/ferncodes/mpfah/` | coeficientes MPFA-H |
| `ferncodes_harmonicopoint` | `legacy/ferncodes/mpfah/` | pontos harmonicos pro MPFA-H |
| `ferncodes_weightnlfvDMP` | `legacy/ferncodes/dmp/` | pesos DMP (usado por MPFA-H, DMP) |

### LPEW2 (pesos de interpolacao lineares-preservadores)
| Funcao | Localizacao | Proposito |
|---|---|---|
| `OPT_Interp_LPEW` | `legacy/ferncodes/lpew/` | gather de geometria por-no (legado) |
| `angulos_Interp_LPEW2` | `legacy/ferncodes/lpew/` | angulos de canto (legado) |
| `netas_Interp_LPEW` | `legacy/ferncodes/lpew/` | razoes netas (legado) |
| `Lamdas_Weights_LPEW2` | `legacy/ferncodes/lpew/` | pesos lambda (legado) |
| `ferncodes_Ks_Interp_LPEW2` | `legacy/ferncodes/lpew/` | projecoes de permeabilidade (legado) |
| `ferncodes_Pre_LPEW_2_vect` | `legacy/ferncodes/lpew/` | ★ driver legado (parcialmente vet) |
| `fs.lpew.OPT` | `+fs/+lpew/` | gather de geometria em batelada (vet) |
| `fs.lpew.v2.angulos` | `+fs/+lpew/+v2/` | angulos em batelada (vet) |
| `fs.lpew.v2.netas` | `+fs/+lpew/+v2/` | netas em batelada (vet) |
| `fs.lpew.v2.ksInterp` | `+fs/+lpew/+v2/` | projecoes de permeabilidade em batelada (vet) |
| `fs.lpew.v2.lambdaWeights` | `+fs/+lpew/+v2/` | lambda vetorizado por dentro (vet) |
| **`fs.lpew.v2.preLPEW2`** | `+fs/+lpew/+v2/` | ★ driver de pipeline completo (vet) |

### Assembly (construcao da matriz)
| Funcao | Localizacao | Proposito |
|---|---|---|
| `ferncodes_globalmatrix_MPFAD` | `legacy/ferncodes/mpfad/` | assembler legado MPFA-D |
| `ferncodes_globalmatrix_TPFA` | `legacy/ferncodes/mpfad/` | assembler legado TPFA |
| `ferncodes_globalmatrix` | `legacy/ferncodes/shared/` | assembler steady/groundwater |
| `ferncodes_assemblematrixMPFAH` | `legacy/ferncodes/mpfah/` | assembler legado MPFA-H (820 L) |
| `ferncodes_assemblematrixMPFAQL` | `legacy/ferncodes/mpfaql/` | MPFA-QL |
| `ferncodes_assemblematrixNLFVPP` | `legacy/ferncodes/nlfvpp/` | NLFV-PP |
| `ferncodes_assemblematrixNLFVH` | `legacy/ferncodes/nlfvh/` | NLFV-H |
| `ferncodes_assemblematrixDMP` | `legacy/ferncodes/dmp/` | preservador DMP |
| **`fs.assembly.mpfad.build`** | `+fs/+assembly/+mpfad/` | ★ MPFA-D totalmente vetorizado |
| `fs.assembly.tpfa.build` | `+fs/+assembly/+tpfa/` | TPFA vetorizado |
| `fs.assembly.{mpfah,mpfaql,nlfvpp,nlfvh,dmp}.build` | `+fs/+assembly/+<x>/` | scaffolds (delegam pro legado) |

### Iteradores nao-lineares (chamados dentro do resolver)
| Funcao | Localizacao | Proposito |
|---|---|---|
| `ferncodes_iterpicard` | `legacy/ferncodes/shared/` | ponto-fixo de Picard |
| `ferncodes_iterpicardANLFVPP2` | `legacy/ferncodes/nlfvpp/` | acel Anderson (forma NLFV-PP) |
| `ferncodes_andersonacc` | `legacy/ferncodes/shared/` | Anderson generico |
| `ferncodes_andersonacc2` | `legacy/ferncodes/shared/` | variante Anderson |
| `L_scheme` | `legacy/unknown/` | iteracao regularizada L-scheme |

### Flow rate & interpolacao de pressao (pos-solve)
| Funcao | Localizacao | Proposito |
|---|---|---|
| `ferncodes_flowrate` | `legacy/ferncodes/mpfad/` | flow-rate MPFA-D (ja vet) |
| `ferncodes_flowrateTPFA` | `legacy/ferncodes/mpfad/` | flow-rate TPFA |
| `ferncodes_flowratelfvHP` | `legacy/ferncodes/mpfah/` | flow-rate MPFA-H |
| `ferncodes_flowratelfvMPFAQL` | `legacy/ferncodes/mpfaql/` | flow-rate MPFA-QL |
| `ferncodes_pressureinterpNLFVPP` | `legacy/ferncodes/nlfvpp/` | interp de pressao nodal (compartilhado) |
| `ferncodes_pressureinterpMPFAQL` | `legacy/ferncodes/mpfaql/` | pressao nodal MPFA-QL |
| `ferncodes_pressureinterpHP` | `legacy/ferncodes/mpfah/` | pressao nodal MPFA-H |
| `fs.flow.mpfad` | `+fs/+flow/` | flow-rate MPFA-D (renomeado) |
| `fs.flow.tpfa` | `+fs/+flow/` | flow-rate TPFA (renomeado) |
| `fs.lpew.pinterp` | `+fs/+lpew/` | interp de pressao nodal compartilhada (renomeada) |

### Drivers de tempo
| Funcao | Localizacao | Proposito |
|---|---|---|
| `main` | raiz | PONTO DE ENTRADA |
| `hydraulic` | `runtime/time/` | solver estacionario de hidraulica |
| `hydraulic_RE` | `runtime/time/` | solver transiente de Richards |
| `IMPES` | `runtime/time/` | Implicit Pressure Explicit Saturation |
| `IMPEC` | `runtime/time/` | Implicit Pressure Explicit Concentration |
| `IMHEC` | `runtime/time/` | (variante) |
| `preRE` | `runtime/preproc/` | preprocessador de Richards |
| `setmethod` | `runtime/time/` | dispatch de time-driver |
| `soil_properties` | `runtime/util/` | termo de acumulacao (dtheta/dt) de Richards |
| `addsource` | `runtime/util/` | insere pocos no sistema |
| `solver` | `runtime/util/` | wrapper de solve linear |
| `postprocessor` | `runtime/util/` | escreve resultados |

### Harness de teste
| Funcao | Localizacao | Proposito |
|---|---|---|
| `fs_setup / fs_expect / fs_teardown / fs_frob` | `tests/helpers/` | primitivas |
| `fs_test_env` | `tests/helpers/` | isolacao WSL-safe do Start.dat |
| `capture_baseline` | `tests/helpers/` | captura de baseline dourado |
| `run_all` | `tests/` | runner agregado |

## Por call graph (caminho de runtime do Richards)

```
main
└─ hydraulic_RE (via setmethod, phasekey=6)
   └─ ferncodes_solver (por passo de tempo)
      ├─ env.metodo.montarSistema
      │  └─ ferncodes_globalmatrix_MPFAD  (ou fs.assembly.mpfad.build)
      │     ├─ ferncodes_Kde_Ded_Kt_Kn (via env.premethod.MPFAD)
      │     └─ ferncodes_Pre_LPEW_2_vect (via env.premethod.MPFAD)
      │        (ou fs.lpew.v2.preLPEW2 — gemeo vetorizado)
      ├─ addsource
      └─ env.metodo.resolver
         └─ ferncodes_iterpicard (ou L_scheme, ou ANLFVPP2)
            └─ por iter Picard:
               ├─ PLUG_kfunction        → novo kmap
               ├─ metodo.atualizarPremethod  → refresh Kde/Ded/Kt/Kn + LPEW2
               ├─ metodo.montarSistema  → rebuild M, I
               ├─ solver(M, I)          → novo p
               └─ check de residuo
```

Por iteracao Picard, os dois hot spots sao:
1. **Recomputacao dos pesos LPEW2** — agora vetorizada (`fs.lpew.v2.preLPEW2`)
2. **Montagem da matriz** — agora vetorizada (`fs.assembly.mpfad.build`)

## Por inventario de codigo morto (veja o estudo `call-graph-and-reachability.md`)

- **~62 arquivos** foram confirmados mortos por grep (nenhum caller de nenhum ponto de entrada)
- **~7000 LOC** de codigo `transm*` esta morto → movido pra `legacy/transm/`
- **7 arquivos `*_con.m`** (variantes com acoplamento de concentracao) mortos → `legacy/ferncodes-con/`
- **32 arquivos** no bucket "unknown" → `legacy/unknown/` pendente de triagem pelo owner
- **4 arquivos** eram tabelas de dados fingindo ser .m (`parametros*.m`, `conduchidraulica.m`, `getchue.m`) — deveriam virar `.mat` (decisao do owner pendente)

## Pegadinhas de disciplina de nomes

- **`ferncodes_pressureinterpNLFVPP`** NAO e especifica de NLFVPP — e uma
  interp de pressao nodal baseada em LPEW usada tambem por MPFA-D.
  Renomeada pra `fs.lpew.pinterp`.
- **`ferncodes_iterpicardANLFVPP2`** — wrapper de aceleracao Anderson,
  usado pelo caso `'AA'` do MPFA-D. Tambem tem nome errado pro seu papel
  transversal.
- **`ferncodes_weightnlfvDMP`** — compartilhado por MPFA-H E DMP.

Se voce ve um nome de funcao sugerindo um metodo mas ela e chamada por
outros, esse e um artefato legado de nomenclatura. Kernels transversais
sao infra compartilhada de verdade — renomea-los pertence a PR-F2
futuro.
