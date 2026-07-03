# Vectorization guide — FlowSim

_Companion to `CHANGELOG.md`. Written for developers extending the `+fs/` tree._

## The FS struct

The vectorization enabler is the **`FS` struct** (built by `fs.mesh.build`):

```
FS
├── mesh       raw connectivity (coord, elem, bedge, inedge, nsurn CSR, esurn CSR)
├── geom       derived-once (centElem, elemArea, normalInt, faceMidBnd/Int, ...)
├── csr        ★ THE KEY VECT ENABLER — flat corner layout for LPEW2
├── perm       tensor / kmap / material
├── bc         nflag, nflagFace, bcflag
├── cfg        numcase, phasekey, pmethod, ...
├── state      p, h, Sw, Con
└── workspace  preallocated triplet buffers (rows/cols/vals/rhs)
```

## CSR corner layout — the vectorization enabler

Legacy LPEW2 loops per node because `nec` (corners per node) is ragged. The **CSR-flat
corner layout** (built by `fs.csr.buildCorners`) flattens this into arrays where each
row is one corner:

```
FS.csr.nCorners       scalar
FS.csr.cornerNode     [nCorners x 1]  which node owns each corner
FS.csr.cornerElem     [nCorners x 1]  which element (= esurn1)
FS.csr.nodePtr        [nNodes+1 x 1]  row-pointer (= esurn2)
FS.csr.nodeNec        [nNodes x 1]    corners per node
FS.csr.maxNec         scalar          max corners on any single node
```

## Recipe: per-node loop → batched arrays

Legacy:
```matlab
for No = 1:nNodes
    for k = 1:nec_No
        result(k) = K11*v(1)^2 + K22*v(2)^2;
    end
end
```

Vectorized:
```matlab
[~, ~, O, Qc] = fs.lpew.OPT(FS);
v = O - Qc;
K11 = perm.tensor(FS.csr.cornerElem, 2);
result = K11 .* v(:,1).^2 + K22 .* v(:,2).^2;
per_node = accumarray(FS.csr.cornerNode, result);
```

## Recipe: sparse assembly → triplet form

```matlab
rows = [rowsD; rowsI; rowsN];
cols = [colsD; colsI; colsN];
vals = [valsD; valsI; valsN];
M = sparse(rows, cols, vals, N, N);   % coalesces duplicates
```

## Recipe: LPEW weight scatter (repelem + cumsum)

For each qualifying edge, scatter contributions across CSR corners:
```matlab
ncQ   = nec_per(vArr(qedges));
startQ= esurn2(vArr(qedges)) + 1;
csum  = cumsum([0; ncQ]);
posInGroup = (1:csum(end))' - repelem(csum(1:end-1), ncQ) - 1;
cornerFlat = repelem(startQ, ncQ) + posInGroup;
```

See `+fs/+assembly/+mpfad/build.m` for the full pattern.

## Testing

Every vectorized function has `tests/unit/unit_<name>.m` that diffs against
the legacy per-node call at 1e-12 relative tolerance.

Golden baselines in `tests/golden/`. Capture with:
```matlab
capture_baseline('mesh', 'M8.msh', 'method', 'mpfad', 'numcase', 439);
```
