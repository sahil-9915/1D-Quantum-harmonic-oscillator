#  Time-Dependent Schrödinger Equation Solver for Harmonic Oscillator

This repository provides a numerical framework for solving the TDSE for a harmonic oscillator potential. The implementation uses the second-order product formula that approximates the time-step operator: $$e^{-i\tau H} \approx e^{-i\tau K_1/2} e^{-i\tau K_2/2} e^{-i\tau V} e^{-i\tau K_2/2} e^{-i\tau K_1/2}$$. Then the time-evolution operator is obtained by concatenating the time-step operator N (= time/tau) times.

QuantumSimulator.py Methods Documentation
Constructor
__init__(L, x_min, x_max) - Initialize simulator with spatial grid parameters (L grid points from x_min to x_max)

Setup & Configuration
discretization() - Create spatial grid and compute step size dx

Initial_wavefunction(sigma, x0) - Set initial Gaussian wavepacket: ψ(x) = (πσ²)^{-1/4} exp(-(x-x₀)²/(2σ²))

discrete_Potential(omega) - Define harmonic oscillator potential V(x) = ½ω²x²
