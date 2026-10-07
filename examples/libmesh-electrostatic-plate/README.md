# Surface charge density on a unit square plate (libMesh backend)

This example computes the surface charge density $\sigma$ on a thin, perfectly conducting unit square plate $\Gamma = [0,1]^2 \times \{0\}$ held at a constant potential $U$ in free space. It uses the libMesh backend of HierBEM with triangular elements.

The total charge density (sum of both sides of the plate) satisfies the first kind integral equation

$$
\frac{1}{4\pi\varepsilon_0} \int_{\Gamma} \frac{\sigma(y)}{|x-y|} \, \mathrm{d}s_y = U, \quad x \in \Gamma,
$$

which is solved by a Galerkin BEM with piecewise constant (P0), continuous linear Lagrange (P1), or continuous quadratic Lagrange (P2) basis functions on triangles. P1 is the default. The Galerkin matrix of the single layer operator is computed with the Sauter–Schwab quadrature for triangle pairs (`include/libmesh_bem/sauter_quadrature_triangle.h`), the triangular counterpart of the deal.II based quadrilateral quadrature in `include/quadrature/sauter_quadrature_tools.h`). All orders use a **dense matrix and dense LU solve**, not an H-matrix.

## Mesh

`unit_square.msh` is a triangular mesh (242 triangles, mesh size 0.1) in the **Gmsh MSH 4.1 binary** format. It is read by `HierBEM::LibMeshBEM::read_gmsh41_binary`, since libMesh's `GmshIO` only supports ASCII files. Meshes with other mesh sizes can be generated with the Gmsh Python API:

```bash
pip install gmsh
python3 generate_unit_square_mesh.py unit_square_fine.msh 0.05
```

For P2, the example promotes the input mesh to libMesh TRI6 elements, sharing edge-midpoint DOFs between adjacent triangles. Geometry remains flat and affine; curved quadratic triangles and hanging-node constraints are not supported.

## Build and run

Build HierBEM with `-DHBEM_USE_LIBMESH=ON` (see `BUILD.md`), then:

```bash
cd <build_dir>/examples/libmesh-electrostatic-plate
./solve-electrostatic-plate [--mesh unit_square.msh] [--potential 1.0] \
  [--same-panel 5] [--common-edge 4] [--common-vertex 4] [--regular 3] \
  [--order 0|1|2] [--output-prefix plate]
./solve-electrostatic-plate --order 1 --output-prefix plate-p1
./solve-electrostatic-plate --order 2 --output-prefix plate-p2
```

The quadrature orders are the numbers of Gauss points per direction in the 4D unit cube for the singular cases, and a triangle Gauss rule of the same polynomial exactness for the regular case, matching the deal.II based `SauterQuadOrder`.

Use `--order 0` to reproduce the original piecewise-constant solve. `--potential` must be finite and nonzero. Higher quadrature orders can be selected to check integration accuracy, particularly for P2.

## Output

* `plate-charge-density.csv`: element centroids, areas and charge densities evaluated at the centroids in C/m².
* `plate-charge-density.msh`: mesh with the charge density as element data (P0) or element-node data (P1/P2) in the Gmsh MSH 2.2 ASCII format, which can be opened in Gmsh directly. P2 output retains the six-node triangles.
* The total charge and the capacitance $C = Q/U$, compared with the reference value $C \approx 0.3667874 \cdot 4\pi\varepsilon_0 \cdot 1\,\mathrm{m} \approx 40.81\,\mathrm{pF}$:

Total charge is integrated using the basis-function integrals, rather than treating nodal coefficients as elementwise constants. The printed density range is the range of coefficients; quadratic interpolants may have extrema between nodes.

Dense results with unit potential and default quadrature:

| Basis | Mesh size | Triangles | DOFs | $C/(4\pi\varepsilon_0)$ | Relative error |
|-------|-----------|-----------|------|-------------------------|----------------|
| P0    | 0.1       | 242       | 242  | 0.36205                 | 1.3 %          |
| P1    | 0.1       | 242       | 142  | 0.36406                 | 0.74 %         |
| P2    | 0.1       | 242       | 525  | 0.36551                 | 0.35 %         |
| P0    | 0.05      | 944       | 944  | 0.36437                 | 0.66 %         |

The charge density is singular at the edges and corners of the plate, which limits the convergence rate on quasi-uniform meshes.
