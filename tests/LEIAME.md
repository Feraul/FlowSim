# Harness de testes do FlowSim

Rode tudo pelo WSL via `tools/mrun`. Esta pasta assume MATLAB R2024a no
lado Windows (veja `../tools/mrun` e a policy AEGIS do projeto de
vetorizacao FlowSim em
`<axon>/my-axon/dev-projects/flowsim-vectorize/_policy.md`).

> **English version**: [`tests/README.md`](README.md)

## Layout

```
tests/
├── README.md            # arquivo em ingles
├── LEIAME.md            # este arquivo (pt-BR)
├── run_all.m            # runner de topo — invoca smoke/*, unit/*, imprime relatorio
├── run_smoke.m          # so os testes de smoke (rapido: sanidade de env, class loads)
├── run_unit.m           # so os testes de unit (medio: oracle-diffs por funcao)
├── helpers/             # infra de teste — assercoes, oracle-diff, fixtures
│   ├── fs_expect.m      # assercao simples (fs_expect(cond, msg))
│   ├── fs_reltol.m      # comparador com tolerancia relativa
│   ├── fs_frob.m        # diff Frobenius relativo em matrizes sparse/dense
│   ├── fs_setup.m       # setup padrao de teste: addpath, seed rng, captura t0
│   ├── fs_teardown.m    # cleanup: reporta tempo, restaura rng
│   └── fs_load_mesh.m   # carrega uma malha canonica por nome (M8, HermelineQuad_12, ...)
├── smoke/               # rapido (segundos) — sanity checks, sem afirmacoes numericas
│   ├── smoke_env.m              # versao MATLAB, cwd WSL, paths configurados
│   ├── smoke_class_hierarchy.m  # verifica estado da hierarquia OOP (quais classes carregam)
│   ├── smoke_mesh_load.m        # todo .msh em meshes/ carrega e reporta contagens de no/elem
│   └── smoke_startdat.m         # Start.dat parseia; reporta numcase, pmethod, phasekey
├── unit/                # medio (10-60 s por teste) — correcao por funcao
│   ├── unit_preprocessormod.m   # build de malha produz formas/contagens esperadas
│   ├── unit_lpew2_reference.m   # saida baseline do LPEW2 numa malha pequena → captura em golden/
│   ├── unit_assembly_mpfad.m    # nnz da montagem MPFA-D, condest, primeiros-5 valores de pressao
│   └── unit_richards_caso439.m  # uma iteracao Picard do Caso439 completa + fingerprint numerico
├── golden/              # baseline outputs capturados (commitados; regenerados com --update)
│   └── (populado por testes de unit na primeira rodada)
└── data/                # fixtures pequenas de teste (malhas minusculas, variantes de Start.dat)
    └── (vazio ate testes precisarem)
```

## Rodando

```bash
# Harness completo (smoke + unit)
tools/mrun -c $(pwd) tests/run_all.m

# So testes de smoke (rapido, ~1 min total)
tools/mrun -c $(pwd) tests/run_smoke.m

# Um teste especifico
tools/mrun -c $(pwd) tests/smoke/smoke_env.m

# Um teste unit com auto-log
tools/mrun -L -c $(pwd) tests/unit/unit_lpew2_reference.m
```

Cada script de teste:
1. Chama `fs_setup` no comeco (paths, seed rng, timer).
2. Usa `fs_expect(cond, msg)` para assercoes — primeira falha imprime mas continua; sumario final reporta contagem pass/fail.
3. Termina com `fs_teardown` que imprime `TEST OK` ou `TEST FAIL: N failures` e sai com 0 ou 1.
4. Testes unit numericos leem/escrevem baselines em `golden/<nome-do-teste>.mat`. Na primeira rodada capturam. Nas rodadas seguintes comparam dentro da tolerancia.

## Tolerancias

Padrao em `fs_frob.m`:
- Matriz esparsa Frobenius relativa: `1e-12`
- Vetor denso L2 relativo: `1e-10`
- Escalar relativo: `1e-10`
- Deriva permitida em `condest`: `10× ratio`

Override por-chamada via `fs_frob(A, B, tol)`.

## Escolhas de design

- **Sem framework de teste third-party** — MATLAB tem `matlab.unittest`
  mas adiciona overhead; para um harness de estudo code-dev, scripts
  planos + `fs_expect` sao mais simples e fazem diff limpo sob version
  control.
- **Arquivos golden sao `.mat`** (binario) — deterministicos, pequenos,
  amigaveis-a-git para escalares mas ignorados para campos grandes via
  `.gitignore` (veja proxima secao).
- **Testes sao roda-veis um-de-cada-vez** entao uma unit que falha nao
  para o batch.
- **Sem side-effects fora de `tests/`** — todo teste cria seu scratch em
  `tests/data/tmp/` e limpa em `fs_teardown`.

## Gitignore

Adicionar ao `.gitignore`:
```
tests/data/tmp/
tests/golden/*.big.mat   # baselines maiores que 100 KB ficam so locais
```
