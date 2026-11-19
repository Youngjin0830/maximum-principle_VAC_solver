# maximum-principle_VAC_solver

MATLAB implementation of the explicit Euler method with a central difference Laplacian
for the vector-valued Allen--Cahn (VAC) equation with $N$ components. The time step
condition keeps every component in $[0,1]$ and is independent of $N$.

## Requirements

* MATLAB (no additional toolboxes required)
* The 3D spatial convergence test with $M=1024$ needs a large amount of memory (about 35 GB for one array of size $1026^3\times4$)

## Folders

* `Fig1_Stability_test_1D/` — Boundedness and discrete energy test in one dimension (Figure 1).
* `Fig2_Stability_test_2D/` — Boundedness and discrete energy test in two dimensions (Figure 2).
* `Fig3_Stability_test_3D/` — Boundedness and discrete energy test in three dimensions (Figure 3).
* `Fig4_Components_test/` — Maximum and minimum values for $N=2,3,\ldots,7$ in two dimensions (Figure 4).
* `Tab1_Temporal_convergence_test/` — Temporal convergence test in 1D, 2D, and 3D (Table 1).
* `Tab2_Spatial_convergence_test/` — Spatial convergence test in 1D, 2D, and 3D (Table 2).

## Numerical method

The VAC equation is discretized using a central difference Laplacian on a cell-centered
grid with the homogeneous Neumann boundary condition and the explicit Euler method in
time. The Lagrange multiplier keeps the sum of the components equal to one.

The time step condition used to keep every component in $[0,1]$ is

$$
\Delta t \leq \frac{2\epsilon^2h^2}{h^2+4d\epsilon^2},
$$

where $h$ is the spatial grid size, $\epsilon$ is the interfacial-width parameter, and
$d$ is the spatial dimension.

## Usage

Open MATLAB, move to the folder of the desired test, and run the script. For example,

```
cd Fig1_Stability_test_1D
Fig1_Stability_test_1D
```

The figures are saved as EPS files in the same folder. The convergence tests print the
errors and the convergence rates.

The parameters, initial conditions, and output settings can be changed near the
beginning of each script.

## Reference

Youngjin Hwang, "Gibbs simplex preservation of the explicit Euler method for the
multi-component Allen--Cahn equation," manuscript.

## Author

Youngjin Hwang
