# legacy/ — Codigo pre-vetorizacao

Tudo sob `legacy/` e codigo que antecede a campanha de vetorizacao de
2026-07 (v2.0.0). Ele **ainda esta no path do MATLAB** — adicionado pelo
`flowsim_init.m` via `genpath` na menor precedencia — entao qualquer
caller que ainda nao migrou pra `+fs/` continua funcionando de forma
transparente.

Para cada simbolo legado que tenha uma substituicao vetorizada, o MATLAB
resolve para a versao `+fs/…` porque `+fs/` e adicionado primeiro
(maior precedencia). Esse e o padrao **shadow** — o legado permanece
como referencia de correcao, o codigo vetorizado assume silenciosamente.

> **English version**: [`legacy/README.md`](README.md)

---

## Indice de clusters

| Cluster | Arquivos (.m) | O que e | Status | Sombreado por |
|---|---:|---|---|---|
| `ferncodes/` | 72 | Os proprios metodos — MPFA-D, MPFA-H, MPFA-QL, NLFV-PP, NLFV-H, DMP, LPEW1/2, kernels compartilhados. Organizados por metodo. | **vivo** | parcialmente (`+fs/+lpew/+v2/`, `+fs/+assembly/{+mpfad,+tpfa}`) |
| `ferncodes-con/` | 7 | Variantes com acoplamento de concentracao (`_con.m`) — codigo morto neste branch | **morto** | nenhum |
| `preprocessor/` | 2 | `preprocessor.m` + `preprocessor2.m` — superados por `runtime/preproc/preprocessormod.m` | **morto** | `runtime/preproc/preprocessormod` |
| `transm/` | 8 | Computacoes de transmissibilidade de tentativa anterior (~7 kLOC de codigo nao-alcancado) | **morto** | nenhum |
| `limiters/` | 5 | Limitadores de fluxo/inclinacao — usados so pela familia de transporte | **vivo** (opt-in) | nenhum |
| `saturation/` | 9 | Helpers de update de saturacao do lado IMPES | **vivo** | nenhum |
| `transport/` | 11 | Advecao/dispersao, sub-solvers hiperbolicos de transporte | **vivo** (opt-in) | nenhum |
| `calc/` | 28 | Helpers `calc*` — fluxos de massa, viscosidades, fluxos espectrais | **parcial** | alguns nao-alcancados |
| `get/` | 31 | Helpers `get*` — extracao de massa, estados de contorno, solucoes analiticas | **parcial** | alguns nao-alcancados |
| `solvers/` | 2 | `solveEnriched` + `solvefluxadvecdispersive` antigos — solvers acoplados | **morto** | `runtime/util/solver*` |
| `utility/` | 25 | Sacola de utilitarios (`Durlofsk`, `L_scheme`, `MalhaKozdon`, `SPEField`, etc.) | **parcial** | `+fs/+iter/{picard,anderson,lscheme}` envolve `L_scheme` |
| `test-scripts/` | 6 | Scripts de teste/plot pre-existentes (`Teste02.m`, `buckey_levalidation.m`, `plot_cilamce2023.m`, `testeBuckley.m`) | **vivo** (opt-in) | nenhum |
| `unknown/` | 32 | Arquivos que nao puderam ser classificados durante a triagem; disposicao pelo owner pendente | **pendente** | nenhum |
| **Total** | **~238** | | | |

## Subclusters de `ferncodes/` (o miolo do legado)

`ferncodes/` — a arvore pre-existente "codigos do Fernando" — e
organizada por metodo de discretizacao para deixar claro quais arquivos
qualquer metodo dado precisa:

| Subdir | Arquivos | Metodo | Gemeo vetorizado |
|---|---:|---|---|
| `ferncodes/shared/` | 18 | kernels transversais entre metodos (`ferncodes_globalmatrix`, `ferncodes_calflag`, `ferncodes_elementface`, `ferncodes_andersonacc*`, …) | nenhum — a maioria ja e vetorizada |
| `ferncodes/lpew/` | 10 | LPEW1 + LPEW2 (`OPT_Interp_LPEW`, `angulos_Interp_LPEW2`, `netas_Interp_LPEW`, `Lamdas_Weights_LPEW2`, `ferncodes_Ks_Interp_LPEW2`, `ferncodes_Pre_LPEW_2_vect`, …) | ★ **`+fs/+lpew/+v2/`** — pipeline completo, 1e-15 rel |
| `ferncodes/mpfad/` | 10 | MPFA-D (`ferncodes_globalmatrix_MPFAD`, `ferncodes_Kde_Ded_Kt_Kn`, `ferncodes_flowrate`, …) | ★ **`+fs/+assembly/+mpfad/build`** + `+fs/+flow/mpfad` |
| `ferncodes/mpfah/` | 8 | MPFA-H (`ferncodes_assemblematrixMPFAH` — 820 linhas, `ferncodes_coefficientmpfaH`, `ferncodes_flowratelfvHP`, `ferncodes_harmonicopoint`, …) | so scaffold (`+fs/+assembly/+mpfah/build` delega) |
| `ferncodes/mpfaql/` | 6 | MPFA-QL (`ferncodes_assemblematrixMPFAQL`, `ferncodes_pressureinterpMPFAQL`, `ferncodes_flowratelfvMPFAQL`, …) | so scaffold |
| `ferncodes/nlfvpp/` | 8 | NLFV-PP (`ferncodes_assemblematrixNLFVPP`, `ferncodes_pressureinterpNLFVPP`, `ferncodes_iterpicardANLFVPP2`, `ferncodes_coefficient`, …) | so scaffold |
| `ferncodes/nlfvh/` | 6 | NLFV-H (`ferncodes_assemblematrixNLFVH`, …) | so scaffold |
| `ferncodes/dmp/` | 6 | Preservador DMP (`ferncodes_assemblematrixDMP`, `ferncodes_weightnlfvDMP`) | so scaffold |

## Caminho de aposentadoria

- Clusters **mortos** — seguro deletar fisicamente uma vez que alguem
  confirme que nenhuma tooling externa depende deles. Candidatos atuais
  para remocao: `ferncodes-con/`, `preprocessor/`, `transm/`, `solvers/`.
- Clusters **vivos** — nao podem ser deletados; aposentar so depois que o
  modulo `+fs/` que sombreia aterrissar + todos os callers migrarem.
- Clusters **parciais** — misturado. Arquivos individuais podem ser
  deletados quando a alcancabilidade a partir de `main.m` → factory →
  benchmark for confirmada vazia. Veja `docs/mapa-de-codigo.md`
  § "Inventario de codigo morto".
- Cluster **pendente** (`unknown/`) — precisa de disposicao do owner. Cada
  arquivo e ou morto (delete) ou perdido na triagem inicial (recategorize).

## Pegadinhas de nomes transversais

Alguns nomes legados sugerem um metodo unico mas o arquivo e infra
compartilhada de verdade. Sao candidatos a rename futuro (adiado como
PR-F2):

| Nome enganoso | Realidade |
|---|---|
| `ferncodes_pressureinterpNLFVPP` | Interp de pressao nodal baseada em LPEW, usada por MPFA-D tambem. Vetorizada como `+fs/+lpew/pinterp`. |
| `ferncodes_iterpicardANLFVPP2` | Wrapper de aceleracao Anderson usado pelo caso `'AA'` do MPFA-D. |
| `ferncodes_weightnlfvDMP` | Compartilhado por MPFA-H **e** DMP. |

## Veja tambem

- `../runtime/LEIAME.md` — arvore ativa de runtime (o que sombreia o que)
- `../docs/mapa-de-codigo.md` — "onde vive a funcao X?" no nivel de funcoes
- `../docs/guia-de-vetorizacao.md` — como adicionar um novo shadow `+fs/`
- `../CHANGELOG.md` — historico completo da campanha
