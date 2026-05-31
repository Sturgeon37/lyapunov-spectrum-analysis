# Lyapunov Spectrum Analysis

Pure Python implementation of the Benettin–Wolf algorithm (tangent-space shadow trajectory method) for computing the full Lyapunov spectrum of continuous-time dynamical systems without explicit construction of the Jacobian matrix.

The implementation is validated on the 3D Lorenz system using a 4th-order Runge–Kutta integrator and a QR/Gram–Schmidt orthogonalization scheme.

---

## Overview

Lyapunov exponents quantify the asymptotic exponential rates of divergence or convergence of infinitesimally close trajectories in phase space. They characterize the stability structure of dynamical systems and provide a standard diagnostic for deterministic chaos.

For a dynamical system:

$$
\frac{d\mathbf{y}}{dt} = f(\mathbf{y}, t)
$$

the Lyapunov spectrum $\{\lambda_i\}$ describes growth rates of perturbations along principal directions of the tangent space.

Key interpretations:
- $\lambda_1 > 0$: sensitive dependence on initial conditions (chaos)
- $\lambda_i < 0$: contraction along corresponding directions
- one exponent in continuous flows is expected to converge to zero (tangent to trajectory)
- $\sum_i \lambda_i$: phase space volume contraction rate

---

## Implementation Status

This repository is currently a reference-grade, research-oriented implementation intended for algorithmic transparency and reproducibility rather than production performance.

The current codebase prioritizes:
- explicit step-by-step numerical procedures
- educational clarity over vectorization
- minimal abstraction layers
- direct implementation of Benettin–Wolf workflow

As a result, the architecture is intentionally not optimized for performance and does not yet integrate with the `numtoolkit` framework or other vectorized solver backends.

Future versions will migrate core numerical routines toward a modular solver architecture while preserving algorithmic equivalence.

## Methodological Idea

Instead of explicitly constructing the Jacobian $J = \partial f / \partial y$, the algorithm approximates tangent-space evolution using a set of perturbation (shadow) trajectories evolved alongside the reference trajectory.

This yields a numerical approximation of the variational equation:

$$
\dot{\delta \mathbf{y}} = J(\mathbf{y}) \, \delta \mathbf{y}
$$

without explicit evaluation of $J$.

---

## Benettin–Wolf Algorithm

The method evolves a reference trajectory and a set of perturbation vectors, periodically re-orthonormalizing them to prevent alignment along the dominant expanding direction.

Let $M$ be the number of exponents.

### Step 1 — Initialization
Initialize orthonormal perturbation vectors:

$$
\mathbf{p}_j(0) = \epsilon \mathbf{e}_j, \quad j = 1, \dots, M
$$

where $\epsilon \ll 1$.

---

### Step 2 — Time Evolution
Integrate:
- reference trajectory $\mathbf{y}(t)$
- perturbed trajectories $\mathbf{y}(t) + \mathbf{p}_j(t)$

using a fixed-step RK4 scheme.

---

### Step 3 — Deviation Construction
Compute deviation vectors:

$$
\mathbf{d}_j = \mathbf{y}_j^{\text{pert}} - \mathbf{y}^{\text{ref}}
$$

These approximate tangent-space directions.

---

### Step 4 — QR / Gram–Schmidt Orthogonalization
Construct an orthogonal basis:

$$
\mathbf{d}_j^\perp = \mathbf{d}_j - \sum_{i < j} \mathrm{proj}_{\mathbf{d}_i^\perp}(\mathbf{d}_j)
$$

This step is equivalent to a QR decomposition of the tangent matrix.

---

### Step 5 — Renormalization and Accumulation
Compute scaling factors:

$$
r_j^{(k)} = \|\mathbf{d}_j^\perp\|
$$

Accumulate:

$$
\lambda_j \approx \frac{1}{N \Delta t} \sum_{k=1}^{N} \ln r_j^{(k)}
$$

Renormalize perturbations:

$$
\mathbf{p}_j \leftarrow \epsilon \frac{\mathbf{d}_j^\perp}{r_j^{(k)}}
$$

---

## Lorenz System

The implementation is tested on the Lorenz attractor:

$$
\frac{dx}{dt} = \sigma (y - x)
$$
$$
\frac{dy}{dt} = x(\rho - z) - y
$$
$$
\frac{dz}{dt} = xy - \beta z
$$

with standard chaotic parameters:
- $\sigma = 10$
- $\rho = 28$
- $\beta = 8/3$

---

## Numerical Setup

- Integrator: RK4 (fixed step)
- Step size: $\Delta t = 0.005$
- Initial perturbation: $\epsilon = 10^{-8}$
- Orthogonalization: Gram–Schmidt (QR-equivalent)
- Spectrum size: 3 exponents (Lorenz system)

---

## Results

Finite-time Lyapunov exponents converge as integration time increases.

| Iterations (N) | Time (T) | λ₁ | λ₂ | λ₃ | Interpretation |
|---:|---:|---:|---:|---|---|
| 1,000 | 5.0 | 0.549 | -0.123 | -14.092 | transient regime |
| 5,000 | 25.0 | 0.951 | -0.078 | -14.540 | partial convergence |
| 100,000 | 500.0 | 0.907 | -0.003 | -14.571 | numerically stable regime |

---

## Consistency Checks

### 1. Zero exponent (flow direction)
The second exponent is expected to converge to zero for continuous-time flows due to invariance along the trajectory direction. Finite-time estimates exhibit small drift due to numerical discretization and finite sampling.

### 2. Volume contraction
The sum of Lyapunov exponents converges to:

$$
\sum \lambda_i \approx -13.666
$$

which is consistent with the divergence of the Lorenz vector field:

$$
-(\sigma + 1 + \beta) = -(10 + 1 + 8/3)
$$

---

## Convergence Visualization

![Lyapunov Convergence](figures/lorenz_lyapunov_convergence.png)

---

## Repository Structure

```text
lyapunov-spectrum-analysis/
├── systems/
│   └── lorenz.py              # Lorenz system + RK4 integrator
├── methods/
│   └── benettin.py            # Benettin–Wolf algorithm implementation
├── figures/
│   └── lorenz_lyapunov_convergence.png
├── examples/
│   └── lorenz_lyapunov.py     # execution script
├── requirements.txt
└── README.md
