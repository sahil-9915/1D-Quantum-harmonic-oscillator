#  Time-Dependent Schrödinger Equation Solver for Harmonic Oscillator

This repository provides a numerical framework for solving the TDSE for a harmonic oscillator potential. The implementation uses the second-order product formula that approximates the time-step operator: $$e^{-i\tau H} \approx e^{-i\tau K_1/2} e^{-i\tau K_2/2} e^{-i\tau V} e^{-i\tau K_2/2} e^{-i\tau K_1/2}$$. Then the time-evolution operator is obtained by concatenating the time-step operator N (= time/tau) times.


