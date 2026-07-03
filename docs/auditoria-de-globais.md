# Inventario de globais + caminho de migracao

_Documenta as principais variaveis globais ainda em uso no codigo legado + os equivalentes na struct FS._
_Atualizado em 2026-07-03 pos-v2.0.1 (paths refletem o layout de `runtime/` + `legacy/`)._

> **English version**: [`docs/globals-audit.md`](globals-audit.md)

## Por que nao "simplesmente matar todas"?

As declaracoes `global` dentro de `legacy/ferncodes/**` nao podem ser removidas
uma-a-uma sem quebrar o runtime — cada `global` num callee esta pareado com uma
expectativa implicita de que a mesma global tenha sido populada por algum
caller (geralmente `runtime/preproc/preprocessormod` no startup). Mata-las
exige **migracao coordenada caller + callee** em conjunto, o que e uma
tarefa mecanica arquivo-a-arquivo melhor feita como uma campanha propria.

Este documento e o **mapa de migracao**: para cada global ainda em uso,
o campo da struct FS que a substitui, e o padrao de migracao.

## Escala do problema (pos-v2.0.1)

- **~183 arquivos `.m`** ainda declaram pelo menos uma `global` (contado via
  `grep -rl "^\s*global\s" --include='*.m' .`)
- Quase todos vivem sob `legacy/ferncodes/**` (67 arquivos) + `runtime/**` +
  alguns poucos remanescentes em `benchmarks/`, `factories/`, `simulacoes/`

## Top-13 globais (por contagem de arquivos-que-usam — atualizada)

| Global | Arquivos | Tipo | Equivalente na struct FS | Notas |
|---|---:|---|---|---|
| `bedge`     | 107 | tabela de faces da malha | `FS.mesh.bedge` | faces de contorno |
| `inedge`    | ~100 | tabela de faces da malha | `FS.mesh.inedge` | faces internas |
| `coord`     | ~100 | coordenadas | `FS.mesh.coord` | coordenadas dos nos |
| `elem`      |  ~85 | conectividade | `FS.mesh.elem` | tabelas de vertices por elemento |
| `centelem`  |  ~85 | geom derivada | `FS.geom.centElem` | centroides de elementos |
| `numcase`   |  58 | tag de fisica | `FS.cfg.numcase` | seletor de benchmark |
| `bcflag`    |  ~35 | mapa de BC | `FS.bc.bcflag` | dispatch flag→valor |
| `normals`   |  ~29 | geom derivada | `FS.geom.normalBnd` / `.normalInt` | normais de face |
| `esurn2`    |  ~25 | ponteiro CSR | `FS.mesh.esurn2` (alias `FS.csr.nodePtr`) | row-pointer de elementos-em-torno-de-no |
| `nsurn2`    |  ~22 | ponteiro CSR | `FS.mesh.nsurn2` | row-pointer de nos-em-torno-de-no |
| `phasekey`  |  22 | tag de fisica | `FS.cfg.phasekey` | dispatch 1/4/5/6 |
| `elemarea`  |  ~19 | geom derivada | `FS.geom.elemArea` | areas por elemento |
| `bcflagc`   |  ~19 | mapa de BC (conc) | `FS.bc.bcflagc` | BCs de concentracao |

_(Contagens sao aproximadas — bedge e exata. Rode `grep -rl "^\s*global.*<nome>"
--include='*.m' .` para atualizar.)_

## Padrao de migracao (mecanico arquivo-a-arquivo)

Para cada arquivo legado `ferncodes_*.m` a ser migrado:

**Antes**:
```matlab
function [out] = ferncodes_foo(other, args)
    global bedge inedge coord elem numcase bcflag;
    % ... usa bedge(:,5), inedge(:,3), etc.
end
```

**Depois**:
```matlab
function [out] = ferncodes_foo(FS, other, args)   % FS adicionado como primeiro arg
    bedge      = FS.mesh.bedge;
    inedge     = FS.mesh.inedge;
    coord      = FS.mesh.coord;
    elem       = FS.mesh.elem;
    numcase    = FS.cfg.numcase;
    bcflag     = FS.bc.bcflag;
    % ... corpo inalterado
end
```

Cada caller entao passa `FS` como um arg extra. Essa e a restricao de
migracao-conjunta: caller E callee mudam juntos por PR.

## O que ja foi migrado (v2.0.x)

- **`+fs/+mesh/build.m`** — le `env.geometry.*` + `env.config.*`, empacota em FS
- **`+fs/+csr/**`** — constroi cantos CSR + shifts, tudo baseado em FS
- **`+fs/+lpew/**`** — le `FS.mesh.*`, `FS.geom.*`, `FS.csr.*`, `FS.perm.*`
- **`+fs/+assembly/+mpfad/build.m`** — le estilo-env diretamente (sem globais)
- **`+fs/+assembly/+tpfa/build.m`** — idem
- **`+fs/+iter/{picard,anderson,lscheme}`** — wrappers baseados em env sobre iteradores legados
- **`+fs/+lpew/dmpWeights`** — wrapper baseado em env sobre `ferncodes_weightnlfvDMP`
- **`+fs/+flow/{mpfad,tpfa}`** — baseado em env

A **camada OOP moderna** (classes `Metodo*` / `Sim*` / `Caso*` em `solvers/`,
`simulacoes/`, `benchmarks/`) ja e baseada em env — nao toca em globais
diretamente. `MetodoBase` fornece `env` como contexto do metodo.

## O que ainda usa globais

- **Todos os arquivos de `legacy/ferncodes/**`** (~67 arquivos `.m`) — pendentes de migracao arquivo-a-arquivo
- **`runtime/time/hydraulic.m`** — o grande driver estacionario, ainda usa globais
- **`runtime/time/hydraulic_RE.m`** — driver transiente de Richards
- **`runtime/time/IMPES.m` / `IMPEC.m` / `IMHEC.m`** — drivers do loop de tempo
- **`runtime/plug/PLUG_*.m`** — alguns leem `numcase`/`phasekey` para dispatch de benchmark
- **`runtime/util/solver.m`** / `solvePressure*.m` — wrappers do solver
- **`benchmarks/Caso439.m`** — le `centelem`, `elem` de globais

## Prioridade para migracao futura

1. **Arquivos no caminho de runtime de Richards** (maior valor):
   `runtime/time/hydraulic_RE`, `ferncodes_iterpicard`, `ferncodes_calflag`,
   `ferncodes_Kde_Ded_Kt_Kn`
2. **Variantes de montagem quando a vetorizacao completa aterrissar**:
   `ferncodes_assemblematrix{MPFAH,NLFVPP,MPFAQL,DMP,NLFVH}` — ja estao
   scaffolded sob `+fs/+assembly/`, mas seus callees ferncodes internos
   ainda usam globais
3. **Pos-processamento**:
   `ferncodes_flowrate*`, `ferncodes_pressureinterp*`
4. **Callbacks PLUG_** — menor superficie, podiam ser feitos em batelada
5. **Todo o resto** — prioridade minima, frequencia do caminho desconhecida

## Nao migrar

- **`legacy/transm/**`** — 8 arquivos, todos mortos. Delete se quiser.
- **`legacy/preprocessor/`** — 2 arquivos, mortos. Superados por `runtime/preproc/preprocessormod`.
- **`legacy/ferncodes-con/`** — 7 variantes mortas com acoplamento de concentracao.
- **`legacy/solvers/`** — 2 variantes de solver mortas.
- **`legacy/unknown/**`** — 32 arquivos pendentes de triagem do owner; sem migracao ate status conhecido.

## Verificacao

Para verificar um arquivo migrado, rode seu teste unitario
(`tests/unit/unit_<nome>.m`) que faz o diff da assinatura nova contra a
assinatura legada na malha M8. Diff esperado:

- **Matrizes** — `< 1e-12` Frobenius relativo
- **Vetores** — `< 1e-10` norma L2 relativa
- **O oraculo** — `unit_baseline_reproduces` deve continuar batendo rel diff
  `0.000e+00` no golden commitado (M8/num439/mpfad e M8/num439/tpfa).

## Custo estimado

Uma migracao de **arquivo unico** e 15-30 min de trabalho mecanico
(declarar param, substituir globais, atualizar todos os callers). Mas como
as migracoes precisam ser em conjunto (caller + callee), cada cascata pode
puxar 3-10 arquivos. Estimativa realista para matar 90% das globais no
conjunto alcancavel: **40-60 horas** como campanha de follow-up dedicada,
provavelmente melhor sequenciada logo apos o trabalho remanescente de
vetorizacao dos 5 metodos (que ja toca no caminho de montagem alcancavel).
