# runtime/ — Codigo dos caminhos ativos de runtime

Os arquivos aqui sao chamados durante uma simulacao FlowSim ao vivo.
Agrupados por papel no pipeline. Todos os quatro subdiretorios sao
adicionados ao path do MATLAB pelo `flowsim_init.m` via `genpath`.

> **English version**: [`runtime/README.md`](README.md)

## Layout dos diretorios

```
runtime/
├── preproc/    preprocessamento de malha + metodo (chamado uma vez por rodada)
├── time/       drivers de time-stepping (chamados uma vez por passo de tempo)
├── plug/       callbacks PLUG_* (chamados a cada substep — o plano "quente")
└── util/       helpers reutilizados pelo pipeline
```

## preproc/ — setup unico

| Arquivo | Chamado por | Papel |
|---|---|---|
| `preprocessormod.m` | `main.m` | le `Start.dat`, constroi `env.geometry` a partir do `.msh` |
| `preprocessmethod.m` | `main.m` | setup por-`pmethod` da struct premethod |
| `preRE.m` | `hydraulic_RE` | preprocessamento especifico de Richards |
| `preSaturation.m` | `IMPES` | preprocessamento especifico de saturacao |
| `preconcentration.m` | `IMPEC` | preprocessamento especifico de concentracao |

## time/ — dispatch de time-driver + loops de passo

| Arquivo | phasekey | Papel |
|---|---|---|
| `setmethod.m` | (qualquer) | dispatch pro driver adequado ao phasekey |
| `hydraulic.m` | 4 | solver estacionario de groundwater / carga hidraulica |
| `hydraulic_RE.m` | 6 | solver transiente da equacao de Richards |
| `IMPES.m` | 1 | Implicit Pressure Explicit Saturation |
| `IMPEC.m` | 5 | Implicit Pressure Explicit Concentration |
| `IMHEC.m` | — | variante (carga hidraulica + concentracao) |

## plug/ — callbacks por-passo (pontos de extensao)

Esses sao os plugins de fisica — troque um para mudar a fisica sem tocar
no driver de tempo. Todos aceitam `(env, parms, tempo)` e retornam o
escalar / matriz que o driver precisa pro substep atual.

| Arquivo | Callback pra | Retorna |
|---|---|---|
| `PLUG_kfunction.m` | tensor de permeabilidade | `env.config.kmap` pro `h` atual |
| `PLUG_bcfunction.m` | BC de pressao | `env.config.nflag`, `nflagface` |
| `PLUG_bcfunction_con.m` | BC de concentracao | idem (pra tracers) |
| `PLUG_bcfunction_con_mpfa_o_fps.m` | BC MPFA-O de concentracao | idem |
| `PLUG_sourcefunction.m` | termo-fonte (pocos, injecao) | fonte por elemento |
| `PLUG_dfunction.m` | dispersao / difusao | tensor |
| `PLUG_Gfunction.m` | gravidade | vetor por elemento |

## util/ — helpers transversais

| Arquivo | Papel |
|---|---|
| `solver.m` | wrapper de solve de sistema linear (backslash + opcoes de precondicionador) |
| `solvePressure.m` / `solvePressure_TPFA.m` | wrapper da equacao de pressao |
| `solveSaturation.m` | wrapper da equacao de saturacao |
| `addsource.m` | insere pocos / termos-fonte no sistema montado |
| `postprocessor.m` | escreve saida VTK / TecPlot / mat |
| `plotandwrite.m` / `plotandwrite_pressfield.m` | plota + persiste campos de resultado |
| `soil_properties.m` | termo de acumulacao `dtheta/dt` de Richards |
| `thetafunction.m` / `theta_n.m` | modelos de conteudo de agua de Richards |
| `gravitation.m` / `calcnormk.m` | helpers de gravidade + norma de permeabilidade |
| `applyinicialcond.m` / `attribinitialcond.m` / `IC.m` | setup de condicao inicial |
| `setrestartinicond.m` / `getrestartdata.m` | manipulacao de arquivo de restart |
| `defineWells.m` | helper de definicao de pocos |
| `benchmark.m` | utilitario de benchmark (lookup de identificador) |

## Ordem de setup do path

`flowsim_init.m` adiciona paths nesta ordem (mais tarde = maior precedencia):

1. `+fs/`  — packages vetorizados (maior precedencia — sombreiam legado)
2. `base/`, `solvers/`, `factories/`, `simulacoes/`, `benchmarks/` — contrato OOP
3. `runtime/**` (esta arvore, via `genpath`)
4. `legacy/**` — fallbacks legados (menor precedencia)

Se um simbolo existe em ambos `+fs/…` e `legacy/…`, o MATLAB resolve pra
versao `+fs/` porque ela vem antes no path.

## Relacao com `legacy/`

Os arquivos aqui sao as copias **vivas** usadas em runtime. Seus gemeos
vetorizados (quando existem) vivem em `+fs/` e os sombreiam de forma
transparente. Arquivos que sao sombreados e ainda chamados a partir daqui
sao:

- `preprocessmethod` → chama `fs.lpew.v2.preLPEW2` (shadow transparente de
  `ferncodes_Pre_LPEW_2_vect`) quando `+fs/` esta no path
- Chamadas de assembly caem primeiro em `fs.assembly.<pmethod>.build`;
  pra `mpfad` e `tpfa` estao totalmente vetorizadas, outras delegam pro legado

## Veja tambem

- `../LEIAME.md` — visao geral do repositorio
- `../docs/mapa-de-codigo.md` — tabela "onde vive a funcao X?" no nivel de funcoes
- `../docs/guia-de-vetorizacao.md` — receita pra estender `+fs/`
- `../legacy/LEIAME.md` — indice de clusters legados
