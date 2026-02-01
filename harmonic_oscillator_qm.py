import numpy as np
import matplotlib.pyplot as plt

class QuantumSimulator:
     
     def __init__(self, L, x_min, x_max):
         
         self.L=L
         self.x_min = x_min
         self.x_max = x_max
         
         self.approx_time_step=None
         self.approx_evolution=None
         self.pdf= None   #probability density
         self.sigma=None
         self.x0=None
         self.omega=None
         self.tau=None
         self.timesteps=None
         self.x_discretization = None
         self.dx = None
         self.V = None
         self.init_wavefunction = None
         self.evolved_wavefunction = None
         self.norm = None
         
         self.discretization()
        

     def discretization(self):

        self.dx = (self.x_max - self.x_min) / (self.L - 1)
        self.x_discretization = np.linspace(self.x_min, self.x_max, self.L)
        


     def Initial_wavefunction(self, sigma, x0):
        self.init_wavefunction=np.zeros(len(self.x_discretization))
        self.norm=0
        n=0
        for i in self.x_discretization:
            self.init_wavefunction[n]= (np.pi * sigma**2)**(-0.25) * np.exp(-(i - x0)**2 / (2 * sigma**2))
            n=n+1
        self.evolved_wavefunction=self.init_wavefunction


     def compute_norm(self):
        psi_squared = np.abs(self.evolved_wavefunction)**2
        norm_squared = np.trapz(psi_squared, self.x_discretization)
        self.norm = np.sqrt(norm_squared)
        return self.norm
     
     def probability_density(self):
         self.pdf = np.abs(self.evolved_wavefunction)**2
         return self.pdf
     
     def plot_probability_density(self):
         prob_density = self.probability_density()
         fig, ax = plt.subplots(figsize=(12, 7))
         ax.plot(self.x_discretization, prob_density, color='darkblue', linewidth=2.5, label='|ψ(x)|²')
         ax.fill_between (self.x_discretization, 0, prob_density, alpha=0.3, color='skyblue', label='Probability')
         
         ax.set_xlim(self.x_min, self.x_max)
         ax.set_ylim(0, 1.1 * np.max(prob_density))

         ax.set_xlabel('Position (x)', fontsize=14, fontweight='bold')
         ax.set_ylabel('Probability Density', fontsize=14, fontweight='bold')
         ax.grid(True, alpha=0.3, linestyle='--')
         ax.legend(fontsize=12, loc='best')
         plt.show()
    

     def discrete_Potential(self, omega):
        #define the discretized Harmonic Oscillator potential
        self.V=  0.5 * (omega**2) * (self.x_discretization**2 )


     def approximate_evolution(self,tau,time):
        #For time step operator U(t) = exp(-iτH)
        # e^{-iτH} ≈ e^{-iτK₁/2} e^{-iτK₂/2} e^{-iτV} e^{-iτK₂/2} e^{-iτK₁/2}
        # where H = K₁ + K₂ + V is the Hamiltonian and K1, K2 and V are hermitian matrices to ensure approx time step operator is unitary.
        # defining e^{-iτV}
        self.tau=tau
        V_diag = (1 / self.dx**2) * (1 + self.dx**2 * self.V)
        exp_V_diag_1D= np.exp(-1j * tau * V_diag)
        exp_V_diag=np.diag(exp_V_diag_1D)
        # defining e^{-iτK₂/2} and e^{-iτK₁/2}
        a = tau / (4 * self.dx**2)
        c = np.cos(a)
        s = np.sin(a)
        is_val = 1j * s
        exp_K1_half = np.zeros((self.L, self.L), dtype=complex)
        exp_K2_half = np.zeros((self.L, self.L), dtype=complex)
        for i in range(0, self.L-1, 2):
            # Only fill if we have a complete 2×2 block
            if i+1 < self.L:
                exp_K1_half[i, i] = c
                exp_K1_half[i, i+1] = is_val
                exp_K1_half[i+1, i] = is_val
                exp_K1_half[i+1, i+1] = c
        
        # If L is odd, last element is just c
        if self.L % 2 == 1:
            exp_K1_half[self.L-1, self.L-1] = 1

        exp_K2_half[0, 0] = 1
        
        for i in range(1, self.L-1, 2):
            # Only fill if we have a complete 2×2 block
            if i+1 < self.L:
                exp_K2_half[i, i] = c
                exp_K2_half[i, i+1] = is_val
                exp_K2_half[i+1, i] = is_val
                exp_K2_half[i+1, i+1] = c
        
        # If L is even, need to handle last element
        if self.L % 2 == 0 and self.L > 0:
            exp_K2_half[self.L-1, self.L-1] = c

        approx_time_step= exp_K1_half @ exp_K2_half @ exp_V_diag @ exp_K2_half @ exp_K1_half
        self.approx_time_step=approx_time_step
        self.approx_evolution = np.eye(self.approx_time_step.shape[0])
        
        self.timesteps= int(time/tau)
        self.approx_evolution = np.linalg.matrix_power(self.approx_time_step, self.timesteps)    
        self.evolved_wavefunction = self.approx_evolution @ self.init_wavefunction

        return self.evolved_wavefunction
        
        
        
     

 
             
    



    


               
    




    

        






