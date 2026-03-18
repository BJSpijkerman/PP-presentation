import numpy as np


"""
This system is analogous to the finite potential well in 2D figure out how to solve that system using the finite
diferences method
"""

# Simulation parameters
k_x_outside = #some number		# Wave vector outside the core
k_x_inside = #some number		# Wave vector inside the core
f = #some number				# Vacuum frequency of the light
epsilon_inside = #some number	# Refractive index inside the core
epsilon_outside = #some number	# Refractive index outside the core
step_size = #some small number	# Simulation step size
L = #some number				# Diameter of the fiber core

# Simulation constants
c = 2.99 * 10 **8

# Calculate numer of steps and create matrix
N_steps = L//step_size
system_matrix = np.zeros((N_steps, N_steps))

# Potential and kinetic analogue terms
potentail_inside = k_x**2 - (2*np.pi*f / c)**2 * epsilon_inside
potential_outside = k_x**2 - (2*np.pi*f / c)**2 * epsilon_outside
kinetic_term = 1 / (2 * step_size)

# Creating the matrix
......

# Solving for the Eigen values and functions
......

# Visualizing
......
