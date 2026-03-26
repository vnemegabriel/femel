# femel — Mindlin Plate FEM in MATLAB

Finite element analysis of thick plates using the **Mindlin–Reissner** theory (first-order shear deformation). The element is a four-node quadrilateral (Q4) with three DOFs per node.

---

## DOF Convention

Each node carries **three degrees of freedom**, ordered as:

| Local index | Meaning |
|-------------|---------|
| 1 | `w` — transverse deflection |
| 2 | `φ_x` — rotation in the x–z plane |
| 3 | `φ_y` — rotation in the y–z plane |

Node `i` occupies global DOFs `[3i-2, 3i-1, 3i]`.

---

## Step-by-Step Workflow

```
1. Define geometry, material, and mesh
        ↓
2. getElementConstitutive   →  Cb, Cs
        ↓
3. getStiffnessMatrix_Plates  (loop over elements)  →  cell arrays {Kele}, {rowIdx}, {colIdx}
        ↓
4. assembleSystem            →  K  (sparse global stiffness)
        ↓
5. applyLoadOnElement        (loop over elements)   →  assemble R
        ↓
6. fixBoundary               →  fixedDofs, freeDofs
        ↓
7. solveSystem               →  D  (global displacement vector)
        ↓
8. getCharacteristicStresses →  M (moments), Q (shear forces)
```

---

## Function Reference

### `getElementConstitutive(material, t)`
Builds the bending and shear constitutive matrices for a homogeneous isotropic plate.

```matlab
material.E  = 210e9;   % Young's modulus [Pa]
material.nu = 0.3;     % Poisson's ratio
material.G  = material.E / (2*(1 + material.nu));
t = 0.1;               % plate thickness [m]

[Cb, Cs] = getElementConstitutive(material, t);
% Cb: (3x3) bending  ->  [Mx; My; Mxy] = Cb * [kappa_x; kappa_y; kappa_xy]
% Cs: (2x2) shear    ->  [Qx; Qy]      = Cs * [gamma_xz; gamma_yz]
```

---

### `getStiffnessMatrix_Plates(nodesEle, eleDofs, Cb, Cs, integrationMethod)`
Computes the element stiffness matrix using the Mindlin plate formulation.
Calls `gauss` once per quadrature scheme, precomputes shape functions outside
the integration loop, then integrates with a **single loop** over all Gauss points.

| `integrationMethod` | Bending pts | Shear pts | Use case |
|---|---|---|---|
| `'full'` | 2×2 | 2×2 | General / thick plates |
| `'selective'` | 2×2 | 1×1 | Avoids shear locking (recommended) |
| `'reduced'` | 1×1 | 1×1 | Thin-plate limit / patch tests |

```matlab
integrationMethod = 'selective';   % recommended for most problems

elementStiffness = cell(nEle, 1);
rowIndex         = cell(nEle, 1);
columnIndex      = cell(nEle, 1);

for iEle = 1:nEle
    eleDofs  = elements.dof(iEle,:)';
    nodesEle = nodes(elements.connectivity(iEle,:), :);
    [elementStiffness{iEle}, ~, ~, rowIndex{iEle}, columnIndex{iEle}] = ...
        getStiffnessMatrix_Plates(nodesEle, eleDofs, Cb, Cs, integrationMethod);
end
```

---

### `assembleSystem(elementStiffness, rowIndex, columnIndex)`
Assembles the sparse global stiffness matrix from element contributions using
`vertcat` + `sparse`.

```matlab
K = assembleSystem(elementStiffness, rowIndex, columnIndex);
```

---

### `applyLoadOnElement(nodesEle, eleDofs, fz, mx, my, npg)`
Computes the element load vector for a distributed transverse load and/or
distributed moment loads, all given at element nodes.

```matlab
R = zeros(nDofTot, 1);

for iEle = 1:nEle
    eleDofs  = elements.dof(iEle,:)';
    nodesEle = nodes(elements.connectivity(iEle,:), :);
    nNodEle  = size(nodesEle, 1);

    fz = q0 * ones(nNodEle, 1);   % uniform transverse load [N/m^2]
    mx = zeros(nNodEle, 1);        % no distributed moment about x
    my = zeros(nNodEle, 1);        % no distributed moment about y

    [Rele, ~] = applyLoadOnElement(nodesEle, eleDofs, fz, mx, my, 2);
    R(eleDofs) = R(eleDofs) + Rele;
end
```

---

### `fixBoundary(nodes, coordIdx, coordVal, dofs2fix, nDofNod)`
Finds all nodes on a boundary edge defined by `nodes(:, coordIdx) == coordVal`
and returns their constrained global DOF indices.

| Boundary condition | `dofs2fix` |
|---|---|
| Simply supported (soft) — `w = 0` only | `1` |
| Simply supported (hard) — `w = 0, φ_n = 0` | `[1 2]` or `[1 3]` |
| Clamped — `w = φ_x = φ_y = 0` | `1:3` |

```matlab
nDofNod = 3;
[fd1, ~] = fixBoundary(nodes, 1,  0, 1, nDofNod);  % x = 0 edge, w = 0
[fd2, ~] = fixBoundary(nodes, 1,  a, 1, nDofNod);  % x = a edge, w = 0
[fd3, ~] = fixBoundary(nodes, 2,  0, 1, nDofNod);  % y = 0 edge, w = 0
[fd4, ~] = fixBoundary(nodes, 2,  a, 1, nDofNod);  % y = a edge, w = 0

fixedDofs = unique([fd1; fd2; fd3; fd4]);
freeDofs  = setdiff((1:nDofTot)', fixedDofs);
```

---

### `solveSystem(K, R, freeDofs)`
Solves the reduced system `K_ff * D_f = R_f` and returns the full displacement
vector (zeros at constrained DOFs).

```matlab
D = solveSystem(K, R, freeDofs);
% D(3*i - 2) = w     at node i
% D(3*i - 1) = phi_x at node i
% D(3*i    ) = phi_y at node i
```

---

### `getCharacteristicStresses(nodes, elements, properties, D, evaluationPoints)`
Computes moment and shear-force resultants at specified isoparametric coordinates.
Shape functions are precomputed at all evaluation points before the element loop.

```matlab
% Evaluate at element centroids (xi=0, eta=0)
evaluationPoints = [0, 0];

[M, Q] = getCharacteristicStresses(nodes, elements, properties, D, evaluationPoints);
% M: (3 x nEvalPoints x nEle)  ->  [Mx; My; Mxy]  [N·m/m]
% Q: (2 x nEvalPoints x nEle)  ->  [Qx; Qy]       [N/m]

Mx_centroid = squeeze(M(1, 1, :));   % Mx at each element centroid
```

The `elements` struct requires fields:
- `elements.connectivity` — `(nEle × nNodEle)` node IDs
- `elements.property`     — `(nEle × 1)` property index per element
- `elements.dof`          — `(nEle × nDofEle)` global DOF indices per element

The `properties` struct array requires fields `.Cb` and `.Cs` per entry.

---

## Low-Level Utilities

| Function | Purpose |
|---|---|
| `gauss(n)` | Tensor-product 2D Gauss points. `n = [n_xi; n_eta]`. Returns `w (npg×1)`, `gp (npg×2)`, `npg` |
| `gauss1D(n)` | 1D Gauss–Legendre points and weights (`n` = 1, 2, or 3) |
| `shapefuns(pts, eleType)` | Shape function values. Returns `(1 × nNodEle × nPts)` for multiple points, `(1 × nNodEle)` for one |
| `shapefunsDer(pts, eleType)` | Shape function derivatives. Returns `(2 × nNodEle × nPts)` for multiple points, `(2 × nNodEle)` for one |
| `meshplot(conn, nodes, color, showtext)` | Quick mesh visualization |

---

## Integration Scheme Recommendation

- **Thick plates** (`t/a > 0.1`): `'full'` integration is appropriate.
- **Thin plates** (`t/a < 0.1`): use `'selective'` to avoid shear locking.
- `'reduced'` is useful for patch tests or the thin-plate limit check.

---

## Running the Example

```matlab
cd matlab
example_simply_supported_plate
```

See `matlab/example_simply_supported_plate.m` for a complete worked problem:
a simply-supported square plate under uniform load, compared against the
Navier series analytical solution (Timoshenko & Woinowsky-Krieger).
