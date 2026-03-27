# Root Cause Analysis: QLQL Deviation from QLLL Benchmark

## Problem Summary
The QLQL element shows significant deviation from the QLLL reference solution in the Scordelis-Lo roof benchmark, particularly at coarser mesh densities (2×2, 4×4, 8×8). QLQL should provide comparable or improved convergence through hierarchical enrichment, but instead it shows worse accuracy.

## Root Cause: **Shear Integration Discrepancy**

### The Fundamental Difference

**QLLL (Selective Integration - Correct)**
- Membrane: 2×2 full Gauss integration
- Flexion (bending): 2×2 full Gauss integration
- Shear: **1-point (reduced/selective) integration** ← Prevents shear locking
- Shear integration call: `getGaussPoints(1)` for 1 Gauss point

**QLQL (Full Integration - Problem)**
- Membrane: 2×2 full Gauss integration
- Flexion: 2×2 full Gauss integration
- Shear: **2×2 full Gauss integration** ← Induces artificial stiffness

### Why This Matters

The Scordelis-Lo roof is a **thin shell structure** (t = 0.25 in, R = 25 in, aspect ratio ≈ 1:100). Thin structures are prone to **shear locking** — a numerical artifact where:

1. **Full integration of shear strain** (2 Gauss points per direction) overstiffens the element
2. The element cannot properly represent pure bending modes without inducing parasitic shear strains
3. This artificially increases the vertical stiffness, reducing vertical displacement w_B

The selective (reduced) 1-point integration used in QLLL naturally alleviates this locking by underintegrating the shear term, allowing more flexible bending deformation.

### Evidence from the Code

**CalcularRigidezQLLL.m (lines 49-65):**
```matlab
% INTEGRACIÓN DE CORTE (Selectiva o Full)
if strcmpi(integrationType, 'full')
    nGaussS = 2;
else
    nGaussS = 1;  % ← Called with 'selective' in main2.m line 53
end
[gpS, gwS] = getGaussPoints(nGaussS);
```

**main2.m (line 53):**
```matlab
Ke = CalcularRigidezQLLL(cords_nodos, geometria, material, 'selective');
```
QLLL receives 'selective' → **nGaussS = 1**

**CalcularRigidezQLQL.m (lines 41-76):**
```matlab
% Integración 2×2 Gauss (no argument for integration type)
[gp, gw] = getGaussPoints(2);
...
% Corte: 2×2 full integration
idxS = [3,4,5, 9,10,11, 15,16,17, 21,22,23];
Kee(idxS, idxS) = Kee(idxS, idxS) + Bs' * Ds * Bs * wJ;
```
QLQL **always uses 2×2 Gauss points** — No option for selective integration

### Design Comment Indicates False Assumption

**CalcularRigidezQLQL.m (lines 12-13):**
```matlab
% Signatura compatible con CalcularRigidezQLLL (integrationType omitido
% porque QLQL siempre usa 2×2 full y no necesita integración selectiva).
```
**Translation:** "...because QLQL always uses 2×2 full and **doesn't need selective integration**."

**This assumption is incorrect for thin shell analysis!** The hierarchical DOFs alone cannot overcome shear locking induced by full integration. The hierarchical enhancement helps with bending representation, but it doesn't suppress the numerical locking mechanism.

## Impact on Convergence

Looking at the w_B convergence graph:
- **QLLL** (with selective shear): w_B → -0.3024 in (reference value) ✓
- **QLQL** (with full shear): w_B significantly underestimates displacement ✗

The shear locking effect is most pronounced at coarser meshes because:
- Coarse elements have larger aspect ratios
- Fewer elements mean the shear locking effect is less diluted by averaging
- As mesh refines, the effect becomes smaller (relative to overall deformation)

## The Fix

**Modify `CalcularRigidezQLQL.m` to support selective integration for shear:**

1. Add an optional `integrationType` parameter (default: 'selective')
2. Apply the same selective integration logic as QLLL for the shear term
3. Update `main2.m` to pass 'selective' to QLQL, or keep the default

This allows QLQL to benefit from:
- Hierarchical enrichment (better bending modes through side DOFs)
- Reduced shear integration (prevents locking through selective integration)

## Verification Strategy

After applying the fix:
1. QLQL should show convergence pattern similar to QLLL
2. w_B values should track QLLL solution across all mesh sizes
3. Stress distributions (Nx', My', Qy') should align better with QLLL
4. Coarser meshes should show improved accuracy relative to current implementation

## References

- Oñate, E. (2009). "Structural Analysis with the Finite Element Method" §6.7.1 (QLLL selective integration) and §6.7.4 (QLQL hierarchical enrichment)
- The selective integration (1-point for shear vs 2×2 for flexion) is a standard technique to suppress shear locking in Mindlin-Reissner plate/shell elements
- Hierarchical enrichment complements but does not replace selective integration for thin structures
