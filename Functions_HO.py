import numpy as np

def discretization(L, x_min, x_max):
    dx = (x_max - x_min) / (L - 1)
    x_discretization = np.linspace(x_min, x_max, L)
    
    # Function body
    print(f"Spacing dx: {dx}")
    return x_discretization, dx



def Initial_wavefunction(sigma,x0, x_discretization):
    init_wavefunction=np.zeros(len(x_discretization))
    norm=0
    n=0
    dx = x_discretization[1] - x_discretization[0]
    for i in x_discretization:
        init_wavefunction[n]= (np.pi * sigma**2)**(-0.25) * np.exp(-(i - x0)**2 / (2 * sigma**2))
        n=n+1
    
    psi_squared = np.abs(init_wavefunction)**2
    norm = np.trapz(psi_squared, x_discretization)
    return init_wavefunction, norm

def discrete_Potential(omega, x_discretization):
    #define the discretized Harmonic Oscillator potential
    V=  0.5 * (omega**2) * (x_discretization**2 )
    return V



def approximate_time_step(tau, V, dx):
    #For time step operator U(t) = exp(-iτH)
    # e^{-iτH} ≈ e^{-iτK₁/2} e^{-iτK₂/2} e^{-iτV} e^{-iτK₂/2} e^{-iτK₁/2}
    # where H = K₁ + K₂ + V is the Hamiltonian and K1, K2 and V are hermitian matrices to ensure approx time step operator is unitary.
    L=len(V)
    # defining e^{-iτV}
    V_diag = (1 / dx**2) * (1 + dx**2 * V)
    exp_V_diag_1D= np.exp(-1j * tau * V_diag)
    exp_V_diag=np.diag(exp_V_diag_1D)
    # defining e^{-iτK₂/2} and e^{-iτK₁/2}
    a = tau / (4 * dx**2)
    c = np.cos(a)
    s = np.sin(a)
    is_val = 1j * s
    exp_K1_half = np.zeros((L, L), dtype=complex)
    exp_K2_half = np.zeros((L, L), dtype=complex)
    for i in range(0, L-1, 2):
        # Only fill if we have a complete 2×2 block
        if i+1 < L:
            exp_K1_half[i, i] = c
            exp_K1_half[i, i+1] = is_val
            exp_K1_half[i+1, i] = is_val
            exp_K1_half[i+1, i+1] = c
    
    # If L is odd, last element is just c
    if L % 2 == 1:
        exp_K1_half[L-1, L-1] = 1

    exp_K2_half[0, 0] = 1
    
    for i in range(1, L-1, 2):
        # Only fill if we have a complete 2×2 block
        if i+1 < L:
            exp_K2_half[i, i] = c
            exp_K2_half[i, i+1] = is_val
            exp_K2_half[i+1, i] = is_val
            exp_K2_half[i+1, i+1] = c
    
    # If L is even, need to handle last element
    if L % 2 == 0 and L > 0:
        exp_K2_half[L-1, L-1] = c

    approx_time_step= exp_K1_half @ exp_K2_half @ exp_V_diag @ exp_K2_half @ exp_K1_half

    return approx_time_step




    
    


               
    




    

        






