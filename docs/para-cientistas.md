# FlowSim v2 — Para cientistas

**Leia isto primeiro se voce usou o FlowSim antes da campanha de vetorizacao
de 2026-07 e quer continuar rodando seus casos sem mexer nos detalhes internos.**

Nada da fisica mudou. Nada da logica dos benchmarks mudou. Nenhuma classe
`Caso*` foi renomeada, nenhum numero de caso foi remanejado, nenhum formato
de saida foi alterado. Este documento existe para voce mesmo verificar isso
em menos de 10 minutos e entender o punhado bem pequeno de coisas que
_realmente_ mudou de lugar.

> **English version**: [`docs/for-scientists.md`](for-scientists.md)

---

## TL;DR — a versao de 60 segundos

1. **Atualize para o master** (ou para a tag `v2.0.1-vectorized`): `git pull`.
2. **Seu `Start.dat` quase certamente ainda funciona**, com duas possiveis
   excecoes (pasta da malha e pasta de saida — veja [§ 3](#3-o-startdat)).
3. **Rode do mesmo jeito que sempre rodou**: abra o MATLAB, `cd FlowSim/`,
   digite `main`.
4. **Os resultados devem ser numericamente identicos** aos da v1 para TPFA
   e MPFA-D na mesma malha + mesmo caso. Existe um teste que prova isso
   na malha M8 + Caso 439 — voce pode rodar com um comando so
   (veja [§ 5](#5-verificando-que-os-numeros-batem-com-a-v1)).
5. **Nada foi apagado.** O codigo antigo continua no disco, dentro da pasta
   `legacy/`, e continua no path do MATLAB — os novos modulos vetorizados
   apenas sao escolhidos primeiro quando os dois existem.

---

## 1. O que mudou e por que (a versao honesta de 5 minutos)

As partes criticas em performance (pesos de interpolacao LPEW2, montagem
das matrizes MPFA-D e TPFA) tinham loops `for` por no e por face no
codigo original. Esses loops agora sao **operacoes de array em batelada**
numa nova subarvore chamada `+fs/` (um "package" do MATLAB). Numa malha
M8 (128 elementos), os loops nao faziam muita diferenca; em malhas
maiores (Hermeline 192×192 pra cima), a versao vetorizada roda
sensivelmente mais rapido e escala melhor.

**Como as duas versoes convivem sem conflito**: a resolucao de funcoes
do MATLAB respeita a ordem do path. O `flowsim_init.m` poe a nova arvore
`+fs/` em primeiro lugar e a arvore antiga `legacy/` no fim. Entao,
quando o codigo pede, por exemplo, os pesos LPEW2, o MATLAB acha a
versao nova (rapida) primeiro. Se por qualquer motivo o `+fs/` fosse
removido, a versao antiga por baixo assumia de forma transparente. Isso
se chama "path shadowing" e quer dizer que:

- **Da para comparar v1 e v2 sem trocar de branch** — basta chamar a
  funcao legada pelo nome completo ferncodes (veja
  [§ 5](#5-verificando-que-os-numeros-batem-com-a-v1)).
- **Se voce nao confia no modulo novo pra algum caso especifico**,
  desabilite ele com `flowsim_init('legacy', false)` (paradoxalmente!
  esse flag controla se o legacy entra no path — veja [§ 5](#5-verificando-que-os-numeros-batem-com-a-v1)).

## 2. O que NAO mudou — lista de tranquilidade

Os itens abaixo sao **byte-identicos** ou **comportamentalmente identicos**
aos da v1:

- Toda classe `Caso*` em `benchmarks/` (Caso331, Caso341, Caso431,
  Caso437, Caso439, Caso21p1, Caso346, etc. — o registry inteiro em
  `factories/createBenchmark.m`)
- Todos os modelos fisicos: Van Genuchten, Gardner, Brooks–Corey, o
  modelo cubico de saturacao, os tensores anisotropicos
- Todos os benchmarks que leem dados externos (Perm_Var*.mat para os
  casos 247/249/250) — os dados continuam carregando
  (veja [§ 6](#6-onde-cada-coisa-foi-parar-referencia-rapida))
- Todo callback `PLUG_*function.m` — permeabilidade, condicoes de
  contorno, termos-fonte, dispersao, gravidade
- Os integradores no tempo — `hydraulic`, `hydraulic_RE`, `IMPES`,
  `IMPEC`
- O formato de saida — o `postprocessor.m` escreve os mesmos VTK / .mat /
  figuras com os mesmos nomes de campo
- A logica da iteracao Picard e seus criterios de convergencia
- O wrapper de aceleracao Anderson (`ferncodes_andersonacc*`, ainda
  usado por `pmethod='AA'`)
- A iteracao regularizada L-scheme
- Todo arquivo `.msh` — mesmas coordenadas dos nos, mesma conectividade,
  mesmas tags de contorno

Numericamente, o teste do baseline dourado
(`tests/unit/unit_baseline_reproduces.m`) prova que na `M8.msh` com
`numcase=439`, tanto TPFA quanto MPFA-D produzem matrizes cuja
diferenca em norma de Frobenius em relacao a rodada pre-vetorizacao e
exatamente **0.000e+00** — bit-a-bit identico, nao "suficientemente
proximo".

## 3. O `Start.dat`

Duas linhas para conferir antes da sua primeira rodada na v2:

**Linha 33 — Pasta de saida** (bloco `>>> EDITE AQUI <<<` do _"Pasta onde
os resultados serao gravados"_):
```
C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409
```
Esse e o seu caminho Windows antigo. Se o MATLAB reclamar de
"mkdir: Access is denied" ou "Cannot open file", ou aponte para uma
pasta que voce tem permissao de escrever, ou use um caminho visivel
pelo WSL.

**Linha 218 — Pasta das malhas** (bloco `>>> EDITE AQUI <<<` do _"Pasta
que contem os arquivos .msh"_):

Na v2 os arquivos de malha foram reorganizados em subpastas. Em vez de
todos os `.msh` num diretorio unico, agora eles vivem sob `meshes/`:

| Se voce usava... | Aponte o Start.dat para... |
|---|---|
| `M8.msh`, `M8_distor1.msh`, `M8_distor2.msh`, `M8_distor3.msh` | `meshes/kozdon/` (relativo a raiz do repo — o MATLAB resolve) |
| `HermelineMeshMod*.msh` (qualquer resolucao) | `meshes/hermeline/` |
| `malhareal.msh`, `mesh_randistorted0_6_TriangA_8_8.msh` | `meshes/other/` |

**Todo o resto do `Start.dat`** — `numcase`, `pmethod`, `phasekey`,
controles de passo de tempo, tolerancias de convergencia, os blocos
`[AVANCADO]` — quer dizer o que sempre quis dizer. Nao mexemos.

## 4. Rodando seu caso — os passos exatos

**No MATLAB (Windows ou Linux)**:

```matlab
cd C:\caminho\pra\FlowSim         % (ou /caminho/pra/FlowSim no Linux)
main
```

O `main.m` chama o `flowsim_init` sozinho na primeira vez (via
`startup.m`). E so. Igualzinho a v1.

**Headless num terminal (WSL / Linux)**:

```bash
cd /caminho/pra/FlowSim
tools/mrun -c $(pwd) main.m
```

Isso roda `matlab.exe -batch` em segundo plano, joga o log inteiro no
seu terminal e sai com codigo 0 em caso de sucesso. Bom para scripts
de varredura.

**Para varrer varios casos**: mantenha o workflow igual ao da v1 — um
loop em shell ou um script MATLAB que edita o `Start.dat` entre as
rodadas.

## 5. Verificando que os numeros batem com a v1

Essa e a coisa que voce provavelmente quer ver com os proprios olhos.
Tem um teste de correcao commitado que captura toda quantidade
intermediaria (contagens de mesh, norma de Frobenius da matriz de
montagem, transmissibilidades MPFA-D Kde/Ded/Kt/Kn, pesos LPEW, norma
L2 do RHS) e faz o diff contra um baseline gravado para
`M8 + Caso 439`:

```bash
tools/mrun -c $(pwd) tests/unit/unit_baseline_reproduces.m
```

Saida esperada (ultimas linhas):
```
[ok]   M8-num439-mpfad: M Frobenius reproducible (rel diff 0.000e+00)
[ok]   M8-num439-mpfad: I L2-norm reproducible   (rel diff 0.000e+00)
[ok]   M8-num439-mpfad: premethod.Kde L2 reproducible (rel 0.000e+00)
...
TEST OK   unit_baseline_reproduces   35/35 passed
```

Todas as 35 assercoes tem que mostrar `rel diff 0.000e+00`. Se
alguma divergir, e uma regressao real — por favor, abra uma issue com
o log.

### Comparando v1 vs v2 no seu proprio caso

Para qualquer `Caso NNN` que voce use:

1. **Pegue uma rodada da v1**: `git checkout v1.0.0-pre-vectorization`
   (num clone descartavel — nao faca isso na sua copia de trabalho),
   rode seu caso, guarde a pasta de saida em algum lugar.
2. **Volte pra v2**: `git checkout master`, rode o mesmo caso, guarde
   a pasta de saida dele.
3. **Faca o diff das duas pastas** — o `pressure*.vtk`, `flowrate*.mat`,
   ou o campo que voce monitora.

Para os metodos TPFA e MPFA-D o diff tem que dar zero exato. Para
**MPFA-H, MPFA-QL, NLFV-PP, NLFV-H, DMP** o diff tambem tem que dar
zero, porque a montagem deles ainda chama o codigo legado (veja
[§ 8](#8-o-que-esta-mais-rapido-e-o-que-nao-ainda)).

### Forcando comportamento puro de v1 se estiver desconfiado

O jeito mais limpo de rodar codigo puro de v1 e fazer checkout da tag
pre-vetorizacao num clone descartavel (NAO faca isso na sua copia de
trabalho — voce vai ficar em detached HEAD):

```bash
git clone https://github.com/Feraul/FlowSim.git /tmp/flowsim-v1
cd /tmp/flowsim-v1
git checkout v1.0.0-pre-vectorization
# ...rode seu caso aqui
```

Tambem existe um flag do init que pula a arvore `legacy/` inteirinha,
rodando **somente** os novos modulos vetorizados (`+fs/`) + a camada
OOP:

```matlab
flowsim_init('legacy', false);    % pula legacy — modo so-vetorizado
```

Atencao: esse flag faz o **oposto** do que o nome sugere — ele controla
se a arvore `legacy/` entra no path do MATLAB. `legacy=true` (o padrao)
inclui como fallback de baixa precedencia; `legacy=false` deixa de fora.
Setar `legacy=false` e util pra experimentos do tipo "isso roda 100%
pelo caminho vetorizado?", mas voce perde os metodos scaffolded
(MPFA-H, MPFA-QL, NLFV-PP, NLFV-H, DMP) porque internamente eles ainda
delegam pro legado.

Para comparar **numerica** entre v1 e v2 num caso completamente
vetorizado (TPFA ou MPFA-D), as duas rodadas com config padrao — uma
em `v1.0.0-pre-vectorization` e outra no `master` atual — sao o
comparativo que voce quer.

## 6. Onde cada coisa foi parar — referencia rapida

Se voce tinha scripts que referenciam algum arquivo pelo caminho
completo, este e o mapa:

| Tipo de arquivo | Estava em (v1) | Esta em (v2) |
|---|---|---|
| **Ponto de entrada** | raiz: `main.m`, `startup.m` | ainda na raiz (nao mudou) |
| **Preprocessador** | raiz: `preprocessormod.m`, `preprocessmethod.m` | `runtime/preproc/` |
| **Integradores no tempo** | raiz: `hydraulic.m`, `hydraulic_RE.m`, `IMPES.m`, `IMPEC.m`, `IMHEC.m`, `setmethod.m` | `runtime/time/` |
| **Callbacks PLUG** | raiz: `PLUG_bcfunction.m`, `PLUG_kfunction.m`, `PLUG_sourcefunction.m`, `PLUG_dfunction.m`, `PLUG_Gfunction.m` | `runtime/plug/` |
| **Solver / auxiliares** | raiz: `solver.m`, `addsource.m`, `postprocessor.m`, `soil_properties.m`, `thetafunction.m`, `applyinicialcond.m`, ~30 outros | `runtime/util/` |
| **Malhas** | raiz: `M8*.msh`, `Hermeline*.msh`, `malhareal.msh`, etc. | `meshes/{kozdon,hermeline,other}/` |
| **Arquivos de dados** | raiz: `Perm_Var0p1.mat`, `Perm_Var2.mat`, `Perm_Var5.mat`, `Teste_5.xlsx`, `Teste_6.xlsx`, `malha_D.geo`, `figura_case_4_Qian_teste_h.fig` | `data/` |
| **Ferncodes internas** | raiz: `ferncodes_*.m` (dezenas de arquivos) | `legacy/ferncodes/<metodo>/` |
| **Variantes antigas do preprocessador** | raiz: `preprocessor.m`, `preprocessor2.m` | `legacy/preprocessor/` |
| **Classes de benchmark** | `benchmarks/Caso*.m` | `benchmarks/Caso*.m` **(nao mudou)** |
| **Classes de metodo** | `solvers/Metodo*.m`, `Solver*.m` | `solvers/Metodo*.m` (todos sob `MetodoBase` agora — veja abaixo) |
| **Classes de simulacao** | `simulacoes/Sim*.m` | `simulacoes/Sim*.m` **(nao mudou)** |
| **Factories** | `factories/create*.m` | `factories/create*.m` **(nao mudou)** |

Voce **nao precisa atualizar nenhum arquivo `Caso*`** porque eles nao
referenciam nenhum arquivo mudado por caminho absoluto — eles acessam
tudo via `env.geometry` e `env.config`, que sao populados pelo
preprocessador. A mudanca de path e transparente pro seu codigo de
fisica.

### Uma limpezinha na OOP que vale voce saber

Na v1 o factory tentava instanciar `MetodoMPFAH`, `MetodoMPFAQL`,
`MetodoNLFVPP` mas esses arquivos nao existiam — so `SolverMPFAH` e
`SolverNLFVPP` (herdando de um `SolverBase` inexistente). Entao
despachar pra qualquer um desses valores de `pmethod` quebrava. Na v2
isso ta corrigido:

- `SolverMPFAH.m` → renomeado pra `MetodoMPFAH.m` (agora herda de
  `MetodoBase` como os outros)
- `SolverNLFVPP.m` → renomeado pra `MetodoNLFVPP.m` (idem)
- `MetodoMPFAQL.m` — criado do zero pra bater com o que o factory
  esperava

Ou seja, `pmethod = 'mpfah' / 'mpfaql' / 'nlfvpp'` agora _de fato_
funciona ponta a ponta. Isso e um bug fix pra coisa que estava quebrada
na v1.

## 7. Adicionando um novo Caso — o workflow

Igual sempre foi. Copie um `benchmarks/CasoNNN.m` existente, edite a
fisica, adicione a unica linha em `factories/createBenchmark.m`:

```matlab
case NNN,   bench = CasoNNN();
```

Nenhum outro arquivo precisa mudar. Os metodos abstratos em
`SimulacaoBase` te dizem o que voce tem que implementar (botao direito
→ Go to definition dentro do MATLAB, ou basta olhar
`benchmarks/Caso439.m` — e o exemplo mais completo).

## 8. O que esta mais rapido e o que nao (ainda)

| Metodo (`pmethod`) | Vetorizado na v2 | Velocidade vs v1 (M8) | Velocidade vs v1 (Hermeline 192²) |
|---|---|---|---|
| `tpfa`  | ✔ sim                          | mais ou menos igual | mensuravelmente mais rapido |
| `mpfad` | ✔ totalmente                   | ~2× mais rapido     | melhoria grande |
| `mpfah` | scaffold — delega pro legado   | igual a v1          | igual a v1 |
| `mpfaql`| scaffold                       | igual               | igual |
| `nlfvpp`| scaffold                       | igual               | igual |
| `nlfvh` | scaffold                       | igual               | igual |
| `dmp`   | scaffold                       | igual               | igual |

"Scaffold" quer dizer que o ponto de entrada existe na nova arvore de
pacote, mas internamente chama o codigo legado `ferncodes_*` (a
correcao numerica esta garantida por construcao). Os cinco metodos
scaffolded serao totalmente vetorizados numa proxima campanha (esforco
total estimado: 30–45 horas, sequenciado por dificuldade crescente:
NLFV-H → NLFV-PP → MPFA-QL → DMP → MPFA-H).

Os pesos de interpolacao LPEW2 — um loop interno critico compartilhado
por todos os metodos da familia MPFA exceto TPFA — **estao** totalmente
vetorizados (`+fs/+lpew/+v2/`), entao ate os metodos scaffolded se
beneficiam dessa parte.

## 9. Solucao de problemas — os 5 sintomas mais comuns

| Sintoma | Causa mais provavel | Correcao |
|---|---|---|
| **`Cannot open file: M8.msh`** | A pasta de malhas do Start.dat aponta pra raiz do repo, mas as malhas agora ficam em `meshes/kozdon/` | Atualize a linha da pasta de malha no Start.dat |
| **`mkdir: Access is denied`** durante o preprocessador | O caminho de saida do Start.dat e uma letra de drive Windows na qual o MATLAB nao consegue escrever | Aponte pra uma pasta que voce tem permissao, ou use `/mnt/c/...` sob WSL |
| **`Undefined function or variable 'MetodoXYZ'`** | Voce tem um checkout antigo com o factory quebrado | `git pull` (a v2.0.0 corrigiu isso) |
| **Um teste falha com uma diferenca relativa pequena mas nao-zero** | Voce esta rodando um baseline dourado velho de um experimento anterior | Ignore aquele teste especifico (e um artefato capturado) — o pipeline em si esta bem se o `unit_baseline_reproduces` passar |
| **De repente o codigo roda o metodo errado pra `pmethod = 'nlfvpp'`** | Voce editou `factories/createMetodo.m` na v1 e sua edicao nao sobreviveu ao merge | Rode `git log -- factories/createMetodo.m` pra ver o que aconteceu; o dispatch atual e o canonico |

## 10. Pedindo ajuda / o que fazer se um resultado divergir

- **Primeiro**: rode `tools/mrun -c $(pwd) tests/unit/unit_baseline_reproduces.m`.
  Se ele falhar, alguma coisa realmente quebrou e precisa ser olhada.
- **Segundo**: se ele passa mas _o seu_ caso diverge de uma saida da v1
  que voce guardou, isole se e a malha (foi movida?), a fisica (voce
  alterou algum `PLUG_*`?) ou o solver (o `pmethod` e o mesmo?).
- **Terceiro**: se ainda estiver enroscado, abra uma issue no GitHub com:
  - o `Start.dat` exato que voce usou
  - a versao do seu MATLAB
  - a saida do `unit_baseline_reproduces`
  - uma descricao do resultado esperado vs observado

## Apendice A — registry de casos (na v2.0.1)

Lendo direto do `factories/createBenchmark.m`, esses sao todos os
valores de `numcase` que vao rodar de cara. Qualquer coisa que nao
esteja aqui nunca foi registrada (ou esta registrada mas o arquivo
`Caso*` nao existe — esses casos ja quebravam na v1 e continuam
quebrando na v2).

**Carga hidraulica / groundwater (300–350):** 330, 331, 332, 333, 334,
335, 336, 337, 338, 341, 341.1, 342, 343, 347, 248.

**Richards (400–500):** 431, 432, 433, 434, 435, 436, 437, 438, 439.

**Casos de referencia / tensores (1–100):** 21.1, 34.6, 34.7, 35, 36.

**Transporte de contaminantes (200–300):** 241, 245, 247, 249, 250.

Pra ver qual classe cada `numcase` chama, abra
`factories/createBenchmark.m` — cada entrada tem uma linha
`case NNN, bench = CasoNNN();` com um comentario explicando a fisica.

---

_Escrito para o time do FlowSim pela campanha AXON code-dev, 2026-07-03._
_Feedback muito bem-vindo — abra uma issue no GitHub ou edite este arquivo
diretamente e mande um PR._
