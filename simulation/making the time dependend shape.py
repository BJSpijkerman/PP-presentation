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

L = 100
N_x = 1000  # Set this once
x = np.linspace(-L, L, N_x)
dx = x[1] - x[0]
box_width = 75
D = 3.0  # Diffusion coefficient
times = [5, 20, 50, 200, 1000]  # From very early to 'late' t

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
plt.plot(x, source_spatial, 'k--', alpha=0.5, label='Source Location (Scaled)')

plt.title(r'Numerical Integration of $\Delta T(\mathbf{r}, t)$')
plt.xlabel(r'Position $\mathbf{r}$')
plt.ylabel(r'Temperature Change $\Delta T$')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

#%%
import pandas as pd

# 1. Recalculate the initial source on the CURRENT x grid to avoid length mismatches
f_box_current = box_function(x, width=box_width)

# 2. Prepare Data for Plot 1 (Diffusion Snapshots)
# We calculate these here to ensure they match the current grid length
t_05 = convolve(f_box_current, gaussian_kernel(x, 0.05, D), mode='same') * dx
t_5  = convolve(f_box_current, gaussian_kernel(x, 5.0, D), mode='same') * dx
t_20 = convolve(f_box_current, gaussian_kernel(x, 20.0, D), mode='same') * dx

df1 = pd.DataFrame({
    'x': x,
    'initial': f_box_current,
    't_0.05': t_05,
    't_5.0': t_5,
    't_20.0': t_20
})

# 3. Prepare Data for Plot 2 (Integrated Total)
df2 = pd.DataFrame({
    'x': x, 
    'delta_T': delta_T,
    'initial': f_box_current  # We include 'initial' here so the dashed line works in LaTeX
})

# 4. Downsample and Save
# Using .copy() avoids a common Pandas warning; iloc[::10] keeps the file small for LaTeX
df1.iloc[::10, :].copy().to_csv('diffusion_evolution.csv', index=False)
df2.iloc[::10, :].copy().to_csv('total_delta_t.csv', index=False)

print(f"Success! CSVs saved with {len(df1)//10} points each.")

#%%

L_0 = 8.5E-6
alpha = 0.55E-6
dndT = 1.1E-5
n_0 = 1.4682

def L(dT):
    dl_1 = alpha*dT
    dl_2 = dndT * dT
    
    dl = L_0*(n_0 * dl_1 + dl_1 + dl_1 *dl_2)
    l = n_0 *L_0
    return l,dl

l,dl = L(delta_T)

plt.plot(x,l+dl)
plt.ylim([0.0, max(l+dl)*1.1])
plt.show()











