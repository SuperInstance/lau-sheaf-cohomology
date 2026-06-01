# lau-sheaf-cohomology

Cellular sheaf cohomology in Rust — coboundary operators, sheaf Laplacians, spectral analysis, and cohomology computation.

## Table of Contents

1. [Overview](#overview)
2. [Mathematical Background](#mathematical-background)
3. [Architecture](#architecture)
4. [API Reference](#api-reference)
5. [Examples](#examples)
6. [Theorems Verified](#theorems-verified)
7. [Installation & Usage](#installation--usage)
8. [License](#license)

---

## Overview

A **cellular sheaf** F on a cell complex X assigns a vector space F(σ) (the *stalk*) to each cell σ, and linear *restriction maps* F(σ⪯τ): F(τ) → F(σ) for each face relation, satisfying functoriality conditions.

This library computes:
- **Sheaf cohomology** H⁰(F) and H¹(F) via Hodge theory (kernel of the sheaf Laplacian)
- **Sheaf Laplacian** eigenvalues, spectral gap, and harmonic sections
- **Cellular (co)homology** — Betti numbers, Euler characteristic, boundary/coboundary matrices
- **Functoriality verification** of restriction maps

All results are verified against 73 property-based tests covering every theorem and invariant.

---

## Mathematical Background

### Cell Complexes

A **cell complex** X is a collection of cells (vertices, edges, faces, …) with boundary relations. This library works with simplicial complexes built from vertices, edges, and triangles.

The **boundary operator** ∂ₖ: Cₖ → Cₖ₋₁ maps each k-cell to a signed sum of its (k-1)-dimensional faces. The fundamental property is:

> **∂² = 0**: the boundary of a boundary is always zero.

This gives a chain complex … → C₂ → C₁ → C₀ → 0, and the **homology** groups are:

> Hₖ(X) = ker(∂ₖ) / im(∂ₖ₊₁)

The **Betti numbers** βₖ = dim Hₖ(X) count topological features: β₀ = number of connected components, β₁ = number of "holes", β₂ = number of "voids".

The **Euler characteristic** is the alternating sum:

> χ(X) = Σₖ (-1)ᵏ · |Cₖ| = Σₖ (-1)ᵏ · βₖ

### Cellular Sheaves

A **cellular sheaf** F on a cell complex X assigns:
- A **stalk** (vector space) F(σ) to each cell σ
- A **restriction map** F(σ⪯τ): F(τ) → F(σ) for each face relation σ ⪯ τ

These must satisfy the **functoriality** axioms:
1. F(σ⪯σ) = id (identity on stalks)
2. F(σ⪯ρ) = F(σ⪯τ) ∘ F(τ⪯ρ) for σ ⪯ τ ⪯ ρ (composition)

### Sheaf Cohomology

The **sheaf coboundary** δₖᶠ: Cᵏ(F) → Cᵏ⁺¹(F) lifts the cellular coboundary to the sheaf setting using restriction maps. The **sheaf cohomology** groups are:

> Hᵏ(F) = ker(δₖᶠ) / im(δₖ₋₁ᶠ)

**Key theorem**: For the constant sheaf (all stalks = ℝ, all restrictions = id), sheaf cohomology recovers ordinary cellular cohomology.

### Hodge Theory and the Sheaf Laplacian

The **sheaf Laplacian** Lₖ = (δₖᶠ)ᵀ δₖᶠ + δₖ₋₁ᶠ (δₖ₋₁ᶠ)ᵀ is a symmetric positive semidefinite operator on k-cochains.

**Hodge decomposition**: every cochain splits orthogonally as:

> Cᵏ(F) = im(δₖ₋₁ᶠ) ⊕ Hᵏ(F) ⊕ im((δₖᶠ)ᵀ)

where the **harmonic space** Hᵏ(F) = ker(Lₖ) is isomorphic to the sheaf cohomology Hᵏ(F).

**Consequences**:
- H⁰(F) = ker(L₀) = space of **global sections** (compatible assignments across all stalks)
- The **spectral gap** (smallest nonzero eigenvalue of L₀) quantifies how tightly the sheaf constrains local-to-global consistency. Spectral gap > 0 on a connected complex implies a unique global section.
- Harmonic sections have **zero Dirichlet energy**: vᵀL₀v = 0.

---

## Architecture

```
src/
├── lib.rs          — All types, algorithms, and example builders
└── main.rs         — Demonstration program
```

The library is structured as a layered pipeline:

```
Cell Complex (topology)
     │
     ▼
Cellular Sheaf (algebraic data on topology)
     │
     ▼
Coboundary Operator (sheaf coboundary matrix)
     │
     ▼
Sheaf Laplacian (spectral analysis)
     │
     ▼
Sheaf Cohomology (H⁰, H¹ dimensions and bases)
```

Each layer is independent and inspectable — you can stop at any level and examine the matrices, ranks, kernels, etc.

### Core Types

| Type | Role |
|------|------|
| `CellComplex` | Cell complex with boundary/coboundary operators, Betti numbers, Euler characteristic |
| `CellularSheaf` | Sheaf with stalks and restriction maps; verifies functoriality |
| `SparseMatrix` | Sparse matrix for boundary/coboundary operators (HashMap-backed) |
| `DenseMatrix` | Dense matrix for stalk computations and Laplacians |
| `CoboundaryOperator` | Wrapper around the sheaf coboundary δₖᶠ |
| `SheafCohomology` | Computed H⁰, H¹, Laplacians, and global sections |
| `SheafLaplacian` | Spectral analysis — eigenvalues (Jacobi iteration), spectral gap, harmonic space, energy |

### Serialization

All major types implement `Serialize`/`Deserialize` via serde + serde_json. Complexes, sheaves, and cohomology results can be persisted to JSON and reloaded.

---

## API Reference

### CellComplex

```rust
// Construction
let mut cx = CellComplex::new();
cx.add_cell(id, dimension, vertices, face_ids);

// Topology
cx.num_cells(k);                // count of k-cells
cx.boundary_matrix(k);          // ∂ₖ: Cₖ → Cₖ₋₁
cx.coboundary_matrix(k);        // δₖ = ∂ₖ₊₁ᵀ: Cᵏ → Cᵏ⁺¹
cx.betti_numbers();             // Vec<usize> of β₀, β₁, …
cx.euler_characteristic();      // χ = Σ(-1)ᵏ|Cₖ|
cx.is_connected();              // connectivity of 1-skeleton
```

### CellularSheaf

```rust
let mut sheaf = CellularSheaf::new(complex);
sheaf.set_stalk(cell_id, dimension);
sheaf.set_restriction(face_id, cell_id, matrix);  // F(σ⪯τ): F(τ) → F(σ)

// Queries
sheaf.stalk_dimension(cell_id);
sheaf.total_stalk_dimension(k);       // total dim of all k-cell stalks
sheaf.verify_functoriality();         // checks identity + composition axioms
sheaf.coboundary_matrix(k);           // sheaf coboundary δₖᶠ
```

### SheafCohomology

```rust
let cohom = SheafCohomology::compute(&sheaf);
cohom.h0_dimension;              // dim H⁰(F)
cohom.h1_dimension;              // dim H¹(F)
cohom.betti;                     // cellular Betti numbers for comparison
cohom.global_sections();         // basis for H⁰ (harmonic 0-cochains)
cohom.laplacian_0();             // L₀ matrix (DenseMatrix)
cohom.laplacian_1();             // L₁ matrix (DenseMatrix)
cohom.sheaf_laplacian_0();       // SheafLaplacian with spectral methods
cohom.sheaf_laplacian_1();       // SheafLaplacian for 1-cochains
```

### SheafLaplacian

```rust
let lap = cohom.sheaf_laplacian_0();
lap.eigenvalues();               // sorted eigenvalues (Jacobi algorithm)
lap.spectral_gap();              // smallest nonzero eigenvalue
lap.harmonic_space();            // basis for ker(L)
lap.is_harmonic(&v);             // check if v is in kernel
lap.energy(&v);                  // vᵀLv (Dirichlet energy)
```

### Pre-built Complexes

```rust
let tri = triangle_complex();           // 2-simplex (3V, 3E, 1F)
let tet = tetrahedron_complex();        // boundary of 3-simplex (4V, 6E, 4F)
let circ = circle_complex(n);           // S¹ with n vertices

let sheaf = constant_sheaf(&cx);        // stalks = ℝ, restrictions = id
let ori = orientation_sheaf(&cx);       // sign-twisted restrictions
```

### Matrix Operations

Both `SparseMatrix` and `DenseMatrix` support:

```rust
// Arithmetic
m.multiply(&other);        // matrix multiplication
m.multiply_vec(&v);        // matrix-vector product
m.transpose();
m.add(&other);
m.rank();                  // via Gaussian elimination
m.nullity();               // cols - rank
m.kernel_basis();          // null space via RREF
```

`DenseMatrix` additionally provides `image_basis()`, `approx_eq()`, `is_zero()`, and `scale()`.

---

## Examples

### Basic: Betti Numbers and Euler Characteristic

```rust
use lau_sheaf_cohomology::*;

// Triangle (contractible → β₀=1, β₁=0)
let tri = triangle_complex();
println!("Betti: {:?}", tri.betti_numbers());     // [1, 0]
println!("Euler: {}", tri.euler_characteristic()); // 1

// Circle S¹ (β₀=1, β₁=1)
let circ = circle_complex(5);
println!("Betti: {:?}", circ.betti_numbers());     // [1, 1]
println!("Euler: {}", circ.euler_characteristic()); // 0

// Tetrahedron surface = S² (β₀=1, β₁=0, β₂=1)
let tet = tetrahedron_complex();
println!("Betti: {:?}", tet.betti_numbers());     // [1, 0, 1]
println!("Euler: {}", tet.euler_characteristic()); // 2
```

### Constant Sheaf Cohomology

```rust
let circ = circle_complex(5);
let sheaf = constant_sheaf(&circ);
let cohom = SheafCohomology::compute(&sheaf);

// Constant sheaf on S¹: H⁰ = ℝ (one global section), H¹ = ℝ (one "twist")
println!("H⁰ = {}, H¹ = {}", cohom.h0_dimension, cohom.h1_dimension); // 1, 1
```

### Spectral Analysis

```rust
let tri = triangle_complex();
let sheaf = constant_sheaf(&tri);
let cohom = SheafCohomology::compute(&sheaf);
let lap = cohom.sheaf_laplacian_0();

println!("Eigenvalues: {:?}", lap.eigenvalues());  // [0.0, λ₁, λ₂]
println!("Spectral gap: {:.4}", lap.spectral_gap()); // λ₁ > 0

// Harmonic sections have zero energy
for v in &lap.harmonic_space() {
    assert!(lap.energy(v).abs() < 1e-10);
}
```

### Custom Sheaf with Projection Maps

```rust
let tri = triangle_complex();
let mut sheaf = CellularSheaf::new(tri.clone());

// Vertices get 1-dim stalks, edges get 2-dim stalks
for cell in tri.cells_of_dimension(0) {
    sheaf.set_stalk(cell.id, 1);
}
for cell in tri.cells_of_dimension(1) {
    sheaf.set_stalk(cell.id, 2);
}

// Restriction: edge → vertex projects onto first coordinate
let proj = DenseMatrix::from_vec(1, 2, vec![1.0, 0.0]);
for cell in &tri.cells {
    for &face_id in &cell.faces {
        if tri.cell(face_id).unwrap().dimension == 0 && cell.dimension == 1 {
            sheaf.set_restriction(face_id, cell.id, proj.clone());
        }
    }
}

let cohom = SheafCohomology::compute(&sheaf);
println!("H⁰ = {}", cohom.h0_dimension);
```

### Higher-Dimensional Stalks

```rust
// ℝ²-valued constant sheaf on triangle → H⁰ = 2, H¹ = 0
let cx = triangle_complex();
let mut sheaf = CellularSheaf::new(cx.clone());
let id2 = DenseMatrix::identity(2);
for cell in &cx.cells {
    sheaf.set_stalk(cell.id, 2);
    sheaf.set_restriction(cell.id, cell.id, id2.clone());
    for &face_id in &cell.faces {
        sheaf.set_restriction(face_id, cell.id, id2.clone());
    }
}
let cohom = SheafCohomology::compute(&sheaf);
assert_eq!(cohom.h0_dimension, 2);
assert_eq!(cohom.h1_dimension, 0);
```

### Orientation Sheaf

```rust
let circ = circle_complex(5);
let sheaf = orientation_sheaf(&circ);
println!("Functoriality: {}", sheaf.verify_functoriality()); // true

let cohom = SheafCohomology::compute(&sheaf);
// Orientation sheaf on orientable S¹: H⁰ = 0 (no compatible global section)
println!("H⁰ = {}", cohom.h0_dimension); // 0
```

### Serialization

```rust
let cx = triangle_complex();
let sheaf = constant_sheaf(&cx);
let cohom = SheafCohomology::compute(&sheaf);

// Serialize to JSON
let json = serde_json::to_string(&cohom).unwrap();
let restored: SheafCohomology = serde_json::from_str(&json).unwrap();
assert_eq!(restored.h0_dimension, cohom.h0_dimension);
```

---

## Theorems Verified

All 73 tests serve as machine-checked proofs of the following theorems:

### Topological Invariants

1. **∂² = 0**: boundary of boundary is zero — `test_boundary_squared_zero`
2. **Euler characteristic formula**: χ = Σ(-1)ᵏ|Cₖ| — `test_euler_characteristic_*`
3. **Euler–Poincaré theorem**: χ = Σ(-1)ᵏβₖ — `test_euler_characteristic_equals_alternating_betti_sum`
4. **Betti numbers**: βₖ = dim ker(∂ₖ) − rank(∂ₖ₊₁) — `test_betti_*`

### Known Spaces

5. **Triangle** (contractible): β₀=1, β₁=0, χ=1
6. **Circle Sⁿ** (any n≥3 vertices): β₀=1, β₁=1, χ=0 — `test_constant_sheaf_different_circles`
7. **Tetrahedron** (S²): β₀=1, β₁=0, β₂=1, χ=2

### Sheaf Theory

8. **Functoriality**: constant and orientation sheaves satisfy F(σ⪯σ)=id and composition law — `test_*_functoriality`
9. **Sheaf cohomology = cellular cohomology** for constant sheaf — `test_sheaf_cohomology_specializes_*`
10. **H⁰(F) = global sections**: dimension matches harmonic space — `test_h0_equals_global_sections`
11. **Higher-dimensional stalks**: ℝⁿ-valued constant sheaf has H⁰ = n — `test_higher_dimensional_stalks`

### Spectral / Hodge Theory

12. **Laplacian is symmetric** — `test_laplacian_symmetry`, `test_laplacian_1_symmetry`
13. **Laplacian is positive semidefinite** — `test_laplacian_positive_semidefinite`
14. **Harmonic sections have zero energy** — `test_harmonic_section_energy_zero`
15. **Spectral gap > 0 ⟹ unique global section** on connected complexes — `test_spectral_gap_implies_unique_global_section`
16. **Eigenvalue count matches harmonic dimension** — `test_laplacian_eigenvalues_*`

---

## Installation & Usage

### Prerequisites

- Rust 1.56+ (2021 edition)

### Add to your project

```toml
[dependencies]
lau-sheaf-cohomology = { git = "https://github.com/SuperInstance/lau-sheaf-cohomology" }
```

Or clone and use locally:

```bash
git clone https://github.com/SuperInstance/lau-sheaf-cohomology.git
cd lau-sheaf-cohomology
cargo run    # Run the demo
cargo test   # Run all 73 tests
```

### Dependencies

| Crate | Purpose |
|-------|---------|
| `serde` + `serde_json` | JSON serialization of complexes, sheaves, and cohomology results |

No external linear algebra or numerics crates — all algorithms (Gaussian elimination, Jacobi eigenvalue, matrix arithmetic) are implemented from scratch for transparency and educational value.

---

## License

MIT
