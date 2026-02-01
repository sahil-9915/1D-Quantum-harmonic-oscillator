#  Time-Dependent Schrödinger Equation Solver for Harmonic Oscillator

This repository provides a numerical framework for solving the TDSE for a harmonic oscillator potential. The implementation uses the second-order product formula that approximates the time-step operator: $$e^{-i\tau H} \approx e^{-i\tau K_1/2} e^{-i\tau K_2/2} e^{-i\tau V} e^{-i\tau K_2/2} e^{-i\tau K_1/2}$$. Then the time-evolution operator is obtained by concatenating the time-step operator N (= time/tau) times.

# QuantumSimulator.py Methods Documentation

## Constructor
- `__init__(L, x_min, x_max)` - Initialize simulator with spatial grid parameters (L grid points from x_min to x_max)

## Setup & Configuration
- `discretization()` - Create spatial grid and compute step size dx
- `Initial_wavefunction(sigma, x0)` - Set initial Gaussian wavepacket: ψ(x) = (πσ²)^{-1/4} exp(-(x-x₀)²/(2σ²))
- `discrete_Potential(omega)` - Define harmonic oscillator potential V(x) = ½ω²x²

## Time Evolution
- `approximate_evolution(tau, time)` - Evolve wavefunction using splitting method
  - Implements: e^{-iτH} ≈ e^{-iτK₁/2} e^{-iτK₂/2} e^{-iτV} e^{-iτK₂/2} e^{-iτK₁/2}
  - τ = time step, time = total evolution time
  - Returns evolved wavefunction and updates probability density/expectation values

## Analysis & Computation
- `compute_norm()` - Calculate and return wavefunction normalization ∫|ψ|²dx
- `plot_probability_density()` - Visualize |ψ(x)|² with shaded probability region
- `plot_expectation_x()` - Plot ⟨x⟩ vs time (0 to 10 with default parameters)
- `plot_expectation_var()` - Plot variance Δx² = ⟨x²⟩ - ⟨x⟩² vs time

## Attributes (Available After Methods)
- `pdf` - Probability density |ψ(x)|² (after evolution)
- `expectation_x` - Position expectation value ⟨x⟩
- `expectation_x2` - ⟨x²⟩ expectation value
- `norm` - Wavefunction normalization
- `evolved_wavefunction` - Current wavefunction state
- `x_discretization` - Spatial grid array
- `V` - Potential energy array
- `dx` - Spatial step size


---
## 📄 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👤 Author
**Sahil Sarbadhikary**  
- GitHub: [@sahil-9915](https://github.com/sahil-9915)

