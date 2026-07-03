# FlowSim — Vetorizado

[![Release](https://img.shields.io/github/v/release/Feraul/FlowSim?sort=semver)](https://github.com/Feraul/FlowSim/releases)
[![License](https://img.shields.io/badge/license-veja_LICENSE-blue.svg)](LICENSE)

> **English version**: [`README.md`](README.md)

**Simulador 2-D de groundwater / Richards / escoamento bifasico em MATLAB.**
Metodos de volumes finitos (TPFA, MPFA-D, MPFA-H, NLFV-PP, MPFA-QL) em
malhas nao-estruturadas de quadrilateros + triangulos. Originalmente
pelo grupo de groundwater / UFPE.

A **campanha de vetorizacao de 2026-07** esta completa e shipada como
[`v2.0.0-vectorized`](https://github.com/Feraul/FlowSim/releases/tag/v2.0.0-vectorized)
mais um follow-up de limpeza de raiz
[`v2.0.1-vectorized`](https://github.com/Feraul/FlowSim/releases/tag/v2.0.1-vectorized).
O estado pre-campanha esta preservado como
[`v1.0.0-pre-vectorization`](https://github.com/Feraul/FlowSim/releases/tag/v1.0.0-pre-vectorization)
(+ branch de backup `legacy-v1.0`). Veja [`CHANGELOG.md`](CHANGELOG.md).

> **Vindo da v1?** Leia [`docs/para-cientistas.md`](docs/para-cientistas.md)
> (ou a versao em ingles [`docs/for-scientists.md`](docs/for-scientists.md))
> — guia em linguagem simples cobrindo o que continuou igual, o que se
> moveu, como verificar que seu `Caso NNN` continua produzindo os mesmos
> numeros, e como desabilitar todo modulo novo se voce quiser
> comportamento puro-v1.

- **LPEW2** pipeline: totalmente vetorizado (`+fs/+lpew/+v2/`), bit-identico ao legado em tolerancia relativa 1e-15
- **MPFA-D** assembly: totalmente vetorizado (`+fs/+assembly/+mpfad/`), bit-identico ao baseline dourado
- **TPFA** assembly: vetorizado (`+fs/+assembly/+tpfa/`)
- **MPFA-H, MPFA-QL, NLFV-PP, NLFV-H, DMP**: scaffolded — pontos de entrada
  existem sob `+fs/+assembly/`, no momento delegam pro legado por
  correcao enquanto o rewrite total pra forma de triplas nao chega
  (veja `CHANGELOG.md` → _Adiado para trabalho futuro_)

O codigo legado permanece sob `legacy/` e e pego pelo path do MATLAB em
menor precedencia — nada quebra, os modulos vetorizados assumem
silenciosamente.

---

## Quick start

```bash
# Clone
git clone https://github.com/Feraul/FlowSim.git
cd FlowSim

# Configure o Start.dat (mesh path, numcase, pmethod, phasekey) — veja docs/como-usar.md

# Rode headless (Linux/WSL):
tools/mrun -c $(pwd) main.m
```

Interativo pelo MATLAB:
```matlab
cd /caminho/para/FlowSim
flowsim_init            % setup do path (+fs/ primeiro, legacy/ ultimo)
main                    % le Start.dat + roda o caso configurado
```

O pipeline le `Start.dat` → pega o caso via `createBenchmark` → pega o
metodo via `createMetodo` → pega o tipo de simulacao via
`createSimulacao` → roda o loop de tempo.

## Testes

```bash
# Passada smoke completa (~30 s):
tools/mrun -c $(pwd) tests/smoke/smoke_env.m
tools/mrun -c $(pwd) tests/smoke/smoke_class_hierarchy.m

# Testes unit de vetorizacao (~1 min cada):
tools/mrun -c $(pwd) tests/unit/unit_lpew2_preLPEW2.m
tools/mrun -c $(pwd) tests/unit/unit_assembly_mpfad.m

# Oraculo de correcao (reproducao bit-identica do golden capturado):
tools/mrun -c $(pwd) tests/unit/unit_baseline_reproduces.m
```

~500 assercoes atualmente verdes — veja `CHANGELOG.md` pro breakdown por-PR.

---

## Documentacao

| Doc | Para | Conteudo |
|---|---|---|
| ★ [`docs/para-cientistas.md`](docs/para-cientistas.md) | **mantenedores / usuarios vindo da v1** | guia em linguagem simples: o que continuou igual, o que se moveu, como verificar que seu caso bate com a numerica da v1 |
| ★ [`docs/for-scientists.md`](docs/for-scientists.md) | maintainers / users coming from v1 | (versao em ingles do acima) |
| [`docs/como-usar.md`](docs/como-usar.md) | usuarios | instalar, configurar `Start.dat`, rodar, troubleshoot |
| [`docs/mapa-de-codigo.md`](docs/mapa-de-codigo.md) | contribuidores | onde cada funcao vive + quem chama quem |
| [`docs/guia-de-vetorizacao.md`](docs/guia-de-vetorizacao.md) | contribuidores | receita pra estender `+fs/` com um novo modulo |
| [`docs/auditoria-de-globais.md`](docs/auditoria-de-globais.md) | contribuidores | inventario das variaveis `global` remanescentes |
| [`runtime/LEIAME.md`](runtime/LEIAME.md) | contribuidores | arvore de runtime ativa (`preproc/time/plug/util`) |
| [`legacy/LEIAME.md`](legacy/LEIAME.md) | contribuidores | indice dos 12 clusters legados + status de aposentadoria |
| [`tests/LEIAME.md`](tests/LEIAME.md) | contribuidores | detalhes do harness de teste |
| [`CHANGELOG.md`](CHANGELOG.md) | todos | historico de release + log completo de PRs |
| [`manual/manual.pdf`](manual/manual.pdf) | cientistas | manual cientifico original de metodos numericos (teoria) |

_Versoes em ingles de todos os docs acima estao lado-a-lado no mesmo
diretorio: `docs/*.md` sem prefixo pt (ex.: `docs/how-to-use.md`,
`docs/code-map.md`), ou `README.md` em vez de `LEIAME.md`._

---

## Mapa do repositorio

Veja `docs/mapa-de-codigo.md` para o mapa detalhado com diagrama de
arquitetura. Contagem enxuta:

- **Raiz:** 4 arquivos `.m` (main, flowsim_init, flowsim_deinit,
  startup) + Start.dat + config/docs. Era 285 na v1.0.0.
- **`+fs/`:** ~30 arquivos em 12 packages (`+util`, `+mesh`, `+csr`,
  `+data`, `+lpew`, `+lpew/+v2`, `+assembly` x7, `+flow`, `+iter`) —
  novo codigo vetorizado
- **`runtime/`:** 37 arquivos em 4 subdirs (`preproc`, `time`, `plug`,
  `util`) — codigo ativo de runtime
- **`legacy/`:** ~250 arquivos em 12 clusters — codigo legado ainda
  alcancavel via precedencia de path
- **`meshes/`:** 16 arquivos `.msh` em `{hermeline,kozdon,other}/`
- **`data/`:** arquivos de dados de benchmark (`Perm_Var*.mat`,
  `Teste_*.xlsx`, etc.)
- **`tests/`:** harness + 4 smoke + 14 unit + 2 goldens (`M8-num439-*`)
- **`docs/`:** 8 arquivos de docs (4 pares EN/PT)
- **`tools/`:** `mrun` (bridge WSL → MATLAB batch)

## Executando o codigo (fluxo de execucao)

Veja `docs/mapa-de-codigo.md` § "Arquitetura de uma olhada" para o
diagrama ASCII completo do fluxo `main.m` → `preprocessormod` →
`createBenchmark` → `preprocessmethod` → loop de tempo → assembly →
solver → postprocessor.

## As tres hierarquias de classe

So uma, na verdade — a campanha de vetorizacao unificou.

```
handle (base do MATLAB)
├── MetodoBase                 (contrato pra todos os metodos numericos)
│    ├── MetodoTPFA            ⬅ baseline
│    ├── MetodoMPFAD           ⬅ target de producao (totalmente vetorizado)
│    ├── MetodoMPFAH           (era SolverMPFAH — renomeado no PR-A2)
│    ├── MetodoNLFVPP          (era SolverNLFVPP)
│    └── MetodoMPFAQL          (nao existia — criado no PR-A2)
│
└── SimulacaoBase              (contrato pra tipos de simulacao)
     ├── SimGroundwater
     ├── SimRichards
     └── Caso439               (o unico benchmark totalmente implementado)
```

## Configuracao (Start.dat)

Edite o `Start.dat` — toda linha marcada com `>>> EDITE AQUI <<<` e um
parametro editavel pelo usuario. Os paths Windows originais do owner
(`C:\Users\flc59\...`) ainda funcionam se voce tiver esse filesystem;
no WSL/Linux, use paths com barras pra frente.

Campos criticos:
- **`numcase`** — qual benchmark rodar (veja `factories/createBenchmark.m`)
- **`pmethod`** — `tpfa`, `mpfad`, `mpfah`, `nlfvpp`, `mpfaql`
- **`phasekey`** — 1 (monofasico), 4 (groundwater), 5 (contaminante), 6 (Richards)
- **pasta da malha** — pasta contendo o arquivo `.msh`
- **nome da malha** — ex.: `M8.msh`
- **pasta de saida** — onde escrever resultados

---

## Status de vetorizacao (em 2026-07-03)

| Componente | Status | Verificacao |
|---|---|---|
| **Computacao de pesos LPEW2** | ✅ Totalmente vetorizada (+fs/+lpew/+v2/) | diff relativo 1e-15 vs legado |
| **Assembly de matriz MPFA-D** | ✅ Totalmente vetorizado (+fs/+assembly/+mpfad/) | Bit-identico ao legado E golden |
| **Assembly de matriz TPFA** | ✅ Vetorizado (legado ja era, so rename) | Bit-identico |
| **Flow-rate MPFA-D** | ✅ (rename; legado ja era vetorizado) | — |
| **Assembly MPFA-H** | 🔧 Scaffold (delega) | Adiado — 820 L de rewrite |
| **Assembly NLFV-PP** | 🔧 Scaffold (delega) | Adiado |
| **Assembly MPFA-QL / DMP / NLFV-H** | 🔧 Scaffolds | Adiado |

O `.gitignore` exclui: `gmsh.exe`, `spe*.dat`, `spe*.mat`, `*.asv`,
`tests/data/tmp/`.

---

## Contribuindo

- Trabalhe no branch `flowsim-artur`; `master` esta estavel / merged.
- Adicione um `tests/unit/unit_<nome>.m` junto de toda mudanca de codigo
  (veja `docs/guia-de-vetorizacao.md` pra receita).
- Verifique via `tools/mrun -c $(pwd) tests/unit/<seu-teste>.m` antes de
  commitar.
- Baselines dourados vivem em `tests/golden/`; regenere com
  `tests/helpers/capture_baseline.m` se a semantica legada mudar.

## Referencias

- `CHANGELOG.md` — cronologia PR-por-PR
- `docs/guia-de-vetorizacao.md` — cookbook de receitas pra futuros modulos vetorizados
- `docs/mapa-de-codigo.md` — deep dive: quem chama quem, onde cada coisa vive
- `manual/manual.pdf` — manual cientifico original (teoria dos metodos numericos)
- `tests/LEIAME.md` — uso do harness de teste

## Licenca / Atribuicao

Codebase original pelo grupo de fluxo em meios porosos da UFPE.
Campanha de vetorizacao (2026): Artur Castiel Reis de Souza + AXON code-dev + Copilot.
