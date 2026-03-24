import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve

def box_function(x, width=4.0):
    """Defines a simple rectangular top-hat function."""
    return np.where(np.abs(x) <= width/2, 1.0, 0.0)

def gaussian_kernel(x, t, D=1.0):
    """Defines the Gaussian diffusion kernel G(x, t)."""
    if t <= 0: return np.zeros_like(x)
    # Variance grows as 2Dt
    sigma = np.sqrt(2 * D * t)
    return (1.0 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-x**2 / (2 * sigma**2))

L=100
# 1. Setup spatial grid
x = np.linspace(-L, L, 2000)
dx = x[1] - x[0]
box_width = 100
D = 1.0  # Diffusion coefficient
times = [0.05, 1.0, 5.0, 20.0]  # From very early to 'late' t

# 2. Create the initial source (box function)
f_box = box_function(x, width=box_width)

# 3. Plotting
plt.figure(figsize=(10, 6))
plt.plot(x, f_box, 'k--', label=r'Initial Box Source ($t \approx 0$)', lw=2)

for t in times:
    # Generate kernel for current time
    kernel = gaussian_kernel(x, t, D)
    
    # Perform the spatial convolution
    # We multiply by dx to maintain the physical normalization of the integral
    convolved = convolve(f_box, kernel, mode='same') * dx
    
    plt.plot(x, convolved, label=f'Time $t$ = {t}')

plt.title('Diffusion Evolution: Box Function Convolved with Gaussian Kernel')
plt.xlabel('Position ($r$)')
plt.ylabel('Temperature Change ($\Delta T$)')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.6)
plt.savefig('diffusion_convolution.png')


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve

def spatial_box(x, width=2.0):
    """The spatial distribution of the heat source."""
    return np.where(np.abs(x) <= width/2, 1.0, 0.0)

def diffusion_kernel_3d(x, tau, D=1.0):
    """The Green's function kernel for 3D diffusion at time lag tau."""
    if tau <= 0: return np.zeros_like(x)
    
    # The prefactor from your image: 1 / [4*pi*D*tau]^(3/2)
    prefactor = 1.0 / (4 * np.pi * D * tau)**(1.5)
    # The exponential part: exp(-|r|^2 / 4D*tau)
    exponent = - (x**2) / (4 * D * tau)
    
    return prefactor * np.exp(exponent)

# Parameters
N_x = 1000          # Spatial resolution
x = np.linspace(-L, L, N_x)
dx = x[1] - x[0]

T_final = 1000       # The "late t" we want to observe
dt = 0.1           # Time step for the integral
t_primes = np.arange(dt, T_final, dt) # Integration samples

# Initialize the total temperature change
delta_T = np.zeros_like(x)

# The spatial source (the box)
source_spatial = spatial_box(x, width=box_width)

# Numerical integration: Summing (Source * Kernel) * dt
for tp in t_primes:
    tau = T_final - tp  # The time elapsed since this heat was released
    kernel = diffusion_kernel_3d(x, tau, D)
    
    # Perform spatial convolution
    # mode='same' keeps the output size equal to x
    step_contribution = convolve(source_spatial, kernel, mode='same') * dx
    
    # Add to the integral (multiplied by the size of the timestep dt)
    delta_T += step_contribution * dt

# Plotting the result
plt.figure(figsize=(10, 6))
plt.plot(x, delta_T, color='firebrick', lw=2, label=f'Total $\Delta T$ at $t={T_final}$')
plt.fill_between(x, delta_T, color='orange', alpha=0.2)
# plt.plot(x, source_spatial * 0.1, 'k--', alpha=0.5, label='Source Location (Scaled)')

plt.title(r'Numerical Integration of $\Delta T(\mathbf{r}, t)$')
plt.xlabel(r'Position $\mathbf{r}$')
plt.ylabel(r'Temperature Change $\Delta T$')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()