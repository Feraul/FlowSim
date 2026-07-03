# Guia de vetorizacao — FlowSim

_Companheiro do `CHANGELOG.md`. Escrito para quem for estender a arvore `+fs/`._

> **English version**: [`docs/vectorization-guide.md`](vectorization-guide.md)

## A struct FS

O habilitador da vetorizacao e a **struct `FS`** (construida por `fs.mesh.build`):

```
FS
├── mesh       conectividade bruta (coord, elem, bedge, inedge, CSR nsurn, CSR esurn)
├── geom       derivadas uma vez (centElem, elemArea, normalInt, faceMidBnd/Int, ...)
├── csr        ★ O HABILITADOR CHAVE DA VETORIZACAO — layout achatado de cantos para LPEW2
├── perm       tensor / kmap / material
├── bc         nflag, nflagFace, bcflag
├── cfg        numcase, phasekey, pmethod, ...
├── state      p, h, Sw, Con
└── workspace  buffers de triplas pre-alocados (rows/cols/vals/rhs)
```

## Layout CSR de cantos — o habilitador da vetorizacao

O LPEW2 legado faz loop por no porque `nec` (cantos por no) e uma quantidade
ragged (irregular). O **layout achatado CSR de cantos** (construido por
`fs.csr.buildCorners`) achata isso em arrays onde cada linha e um canto:

```
FS.csr.nCorners       escalar
FS.csr.cornerNode     [nCorners x 1]  qual no possui cada canto
FS.csr.cornerElem     [nCorners x 1]  qual elemento (= esurn1)
FS.csr.nodePtr        [nNodes+1 x 1]  row-pointer (= esurn2)
FS.csr.nodeNec        [nNodes x 1]    cantos por no
FS.csr.maxNec         escalar         maximo de cantos num unico no
```

## Receita: loop por no → arrays em batelada

Legado:
```matlab
for No = 1:nNodes
    for k = 1:nec_No
        result(k) = K11*v(1)^2 + K22*v(2)^2;
    end
end
```

Vetorizado:
```matlab
[~, ~, O, Qc] = fs.lpew.OPT(FS);
v = O - Qc;
K11 = perm.tensor(FS.csr.cornerElem, 2);
result = K11 .* v(:,1).^2 + K22 .* v(:,2).^2;
per_node = accumarray(FS.csr.cornerNode, result);
```

## Receita: montagem esparsa → forma de triplas

```matlab
rows = [rowsD; rowsI; rowsN];
cols = [colsD; colsI; colsN];
vals = [valsD; valsI; valsN];
M = sparse(rows, cols, vals, N, N);   % coalesce duplicatas automaticamente
```

## Receita: scatter de pesos LPEW (repelem + cumsum)

Para cada face qualificada, espalhar contribuicoes pelos cantos CSR:
```matlab
ncQ   = nec_per(vArr(qedges));
startQ= esurn2(vArr(qedges)) + 1;
csum  = cumsum([0; ncQ]);
posInGroup = (1:csum(end))' - repelem(csum(1:end-1), ncQ) - 1;
cornerFlat = repelem(startQ, ncQ) + posInGroup;
```

Ver `+fs/+assembly/+mpfad/build.m` para o padrao completo.

## Testes

Toda funcao vetorizada tem um `tests/unit/unit_<nome>.m` que faz o diff
contra a chamada legada por-no com tolerancia relativa de 1e-12.

Baselines dourados em `tests/golden/`. Capture com:
```matlab
capture_baseline('mesh', 'M8.msh', 'method', 'mpfad', 'numcase', 439);
```
