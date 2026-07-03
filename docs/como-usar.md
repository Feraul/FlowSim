# Como usar o FlowSim

_Guia voltado ao usuario. Para contribuidores, veja `LEIAME.md` + `docs/mapa-de-codigo.md` + `docs/guia-de-vetorizacao.md`._

> **English version**: [`docs/how-to-use.md`](how-to-use.md)

## Pre-requisitos

- **MATLAB R2019b ou mais novo** (testado com R2024a no Windows via WSL)
- **Gmsh** (install pelo sistema: `apt install gmsh` no Linux, `choco install gmsh` no Windows) — usado so pra editar arquivos de malha, nao e requerido pra RODAR o simulador

## Setup inicial

1. **Clone o repo** (o release fica no `master`):
   ```bash
   git clone https://github.com/Feraul/FlowSim.git
   cd FlowSim
   # master = release atual v2.0.1-vectorized
   # pro backup pre-vetorizacao, use: git checkout v1.0.0-pre-vectorization
   ```

2. **Edite o `Start.dat`** — o arquivo de configuracao de runtime. Toda
   linha marcada com `>>> EDITE AQUI <<<` e um parametro do usuario. Os
   mais importantes:
   - **Pasta de saida** (linha 33) — onde os resultados vao. Use barras
     pra frente no Linux/WSL, barras invertidas no Windows.
   - **Pasta da malha** (linha 218) — pasta contendo os arquivos `.msh`.
     Todas as malhas agora vivem sob `meshes/{hermeline,kozdon,other}/`
     — aponte o Start.dat pra subpasta certa.
   - **Nome do arquivo de malha** (linha 224) — ex.: `M8.msh`.
   - **`numcase`** — qual caso de benchmark (veja
     `factories/createBenchmark.m` pra lista completa).
   - **`pmethod`** — um de: `tpfa` / `mpfad` / `mpfah` / `nlfvpp` / `mpfaql`.
   - **`phasekey`** — familia de fisica:
     - `1` = monofasico pressao/saturacao
     - `4` = carga hidraulica de groundwater
     - `5` = contaminante + carga hidraulica
     - `6` = Richards (parcialmente saturado)

3. **Confirme os dados de permeabilidade SPE10** (opcional — so pra
   casos SPE10). Por padrao eles foram relocados pra fora do repo
   (`spe10.mat`, `spe_perm.dat`, ~90 MB no total). Se um benchmark
   precisar deles:
   - Faca um symlink a partir da localizacao que `fs.data.paths('spe10')`
     reportar, OU
   - Defina a variavel de ambiente `FS_DATA_DIR` onde voce guarda os
     arquivos.

## Rodando pelo MATLAB (interativo)

```matlab
cd /caminho/para/FlowSim
flowsim_init          % configura todos os paths (so precisa uma vez por sessao)
main                  % le o Start.dat, roda o caso configurado
```

## Rodando pelo WSL / terminal (headless)

```bash
cd /caminho/para/FlowSim

# Rodada unica do main.m:
tools/mrun -c $(pwd) main.m

# Ou um script especifico:
tools/mrun -c $(pwd) /caminho/pro/meuscript.m

# One-liner:
tools/mrun -e "disp('oi'); disp(version)"

# Com timeout (padrao 30 min):
tools/mrun -t 3600 -c $(pwd) main.m

# Com auto-log em /tmp/mrun-logs/:
tools/mrun -L -c $(pwd) main.m
```

O `mrun` embrulha `matlab.exe -batch` e:
- Filtra os warnings inocuos de change-notification UNC
- Auto-cd na pasta do script (ou um override `-c`)
- Propaga o exit code (0 = ok, non-0 = erro do MATLAB, 124 = timeout)

## Rodando a suite de testes

```bash
# Passada smoke completa (rapido, ~30 s):
tools/mrun -c $(pwd) tests/smoke/smoke_env.m

# Todos os smoke + unit (~5 min):
tools/mrun -c $(pwd) tests/run_all.m

# Um teste especifico (com progresso verbose):
tools/mrun -c $(pwd) tests/unit/unit_baseline_reproduces.m
```

Todos os testes usam `tests/helpers/fs_test_env.m` que isola o
Start.dat num diretorio temporario — seu `Start.dat` real nunca e
mutado pelos testes.

## Adicionando um novo caso de benchmark

1. **Crie `benchmarks/CasoNNN.m`** — herde de `SimulacaoBase`, implemente
   os metodos abstratos (veja `base/SimulacaoBase.m` pro contrato).
   Modele a partir do `benchmarks/Caso439.m` existente — e o unico
   totalmente implementado.

2. **Registre em `factories/createBenchmark.m`**:
   ```matlab
   case NNN,  bench = CasoNNN();
   ```

3. **Atualize o `Start.dat`**: setar `numcase = NNN`.

4. **Opcionalmente adicione um smoke test** em `tests/smoke/` que
   verifica que a classe carrega
   (`meta.class.fromName('CasoNNN')`).

## Adicionando um novo metodo numerico

1. **Crie `solvers/MetodoXXX.m`** — herde de `MetodoBase`, implemente:
   - `preprocessar(env, parms)` — precompute quantidades mesh-invariant
   - `atualizarPremethod(env, parms)` — refresh em novo kmap (Richards)
   - `montarSistema(env, parms, dt)` — retorna `[M, I]`
   - `resolver(M, I, parms, env, tempo, dt, source_wells)` — resolve (linear ou Picard)
   - `calcularFlowrate(p, env, parms)` — postprocess

2. **Registre em `factories/createMetodo.m`**:
   ```matlab
   case 'xxx',  metodo = MetodoXXX();
   ```

3. **Atualize o `Start.dat`**: `pmethod = xxx`.

4. **Se quiser a estrutura de package vetorizada**: crie
   `+fs/+assembly/+xxx/build.m` (delegate pro legado a principio;
   incrementalmente porte pra forma de triplas).

## Pegadinhas comuns

- **"Undefined function or variable"** em alguma `Caso...` ou `Metodo...`
  — o factory espera um arquivo de classe que nao existe. So `Caso439` +
  as 5 classes `MetodoXXX` estao wired neste branch.
- **"mkdir Access is denied"** durante o preprocessormod — o caminho de
  saida do Start.dat e Windows-style mas o MATLAB nao consegue criar
  (path nao existe ou sem permissao de escrita). Corrija o caminho de
  saida.
- **"Cannot open file: ...\some_mesh.msh"** — a pasta de malhas e o
  nome de arquivo de malha do Start.dat nao batem com um arquivo real.
  As malhas agora vivem sob `meshes/{hermeline,kozdon,other}/` —
  atualize o caminho de malha do Start.dat.
- **Teste baseline falhando** — provavelmente uma mudanca numerica se
  infiltrou. Compare `norm(M_new - M_leg, 'fro') / norm(M_leg, 'fro')`.
  Se < 1e-12, ajuste a tolerancia; se > 1e-12, tem um bug de verdade.

## Onde olhar em busca de problemas

| Sintoma | Olhe em |
|---|---|
| MATLAB nao acha uma funcao | `flowsim_init.m` — a arvore legacy esta no path? |
| Pressoes / flowrates errados | Teste de baseline dourado + `tests/golden/*.mat` |
| Classe nao instancia | `tests/smoke/smoke_class_hierarchy.m` |
| Problemas de setup de path | Rode `flowsim_init('verbose', true)` pra ver a ordem addpath |
| Problemas de path WSL | Adicione `-v` ao mrun pra ver o MATLAB_CMD real |

## Limpeza

- **Arquivos scratch de teste** vivem em `/tmp/flowsim_test_env/`
  (recriados por teste).
- **Auto-logs** do `mrun -L` vivem em `/tmp/mrun-logs/`.
- **Autosaves `.asv` do MATLAB** — gitignored, seguro deletar.

## Pedindo ajuda

- `manual/manual.pdf` — referencia cientifica original (teoria dos metodos numericos)
- `CHANGELOG.md` — o que mudou e quando
- `docs/guia-de-vetorizacao.md` — receita pra estender a arvore `+fs/`
- `docs/mapa-de-codigo.md` — onde cada funcao vive + quem chama quem
- `tests/LEIAME.md` — detalhes do harness de teste
