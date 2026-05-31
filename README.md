# lau-sheaf-cohomology

Cellular sheaf cohomology in Rust — coboundary operators, sheaf Laplacians, spectral analysis, and cohomology computation.

## What is this?

A **cellular sheaf** F on a cell complex X assigns a vector space F(σ) (the *stalk*) to each cell σ, and linear *restriction maps* F(σ⪯τ): F(τ) → F(σ) for each face relation, satisfying functoriality conditions.

This library computes:
- **Sheaf cohomology** H⁰(F) and H¹(F) via Hodge theory (kernel of the sheaf Laplacian)
- **Sheaf Laplacian** eigenvalues, spectral gap, and harmonic sections
- **Cellular (co)homology** — Betti numbers, Euler characteristic, boundary/coboundary matrices
- **Functoriality verification** of restriction maps

## Core Types

| Type | Description |
|------|-------------|
| `CellComplex` | Cell complex with boundary/coboundary operators, Betti numbers, Euler characteristic |
| `CellularSheaf` | Sheaf with stalks and restriction maps |
| `SparseMatrix` | Sparse matrix for boundary/coboundary operators |
| `DenseMatrix` | Dense matrix for stalk computations |
| `CoboundaryOperator` | Sheaf coboundary δₖᶠ |
| `SheafCohomology` | Computed H⁰, H¹, Laplacians |
| `SheafLaplacian` | Spectral analysis (eigenvalues, spectral gap, energy) |

## Examples

```rust
use lau_sheaf_cohomology::*;

// Triangle complex (contractible)
let cx = triangle_complex();
println!("Betti: {:?}", cx.betti_numbers()); // [1, 0]
println!("Euler: {}", cx.euler_characteristic()); // 1

// Constant sheaf
let sheaf = constant_sheaf(&cx);
let cohom = SheafCohomology::compute(&sheaf);
println!("H⁰ = {}, H¹ = {}", cohom.h0_dimension, cohom.h1_dimension); // 1, 0

// Spectral analysis
let lap = cohom.sheaf_laplacian_0();
println!("Eigenvalues: {:?}", lap.eigenvalues());
println!("Spectral gap: {}", lap.spectral_gap());

// Circle S¹
let circ = circle_complex(5);
let sheaf = constant_sheaf(&circ);
let cohom = SheafCohomology::compute(&sheaf);
println!("H⁰ = {}, H¹ = {}", cohom.h0_dimension, cohom.h1_dimension); // 1, 1
```

## Theorems Verified

1. **Euler characteristic**: χ = Σ(-1)^k |Cₖ|
2. **Betti numbers**: Hₖ = ker(∂ₖ)/im(∂ₖ₊₁)
3. **∂² = 0**: boundary of boundary is zero
4. **Sheaf cohomology = cellular cohomology** for constant sheaf (stalks = ℝ, restrictions = id)
5. **H⁰(F) = global sections** (compatible assignments)
6. **Functoriality**: F(σ⪯σ) = id and F(σ⪯ρ) = F(σ⪯τ) ∘ F(τ⪯ρ)
7. **Spectral gap > 0 ⟹ unique global section** on connected complexes
8. **Energy of harmonic section = 0**
9. **Laplacian is symmetric positive semidefinite**

## License

MIT
