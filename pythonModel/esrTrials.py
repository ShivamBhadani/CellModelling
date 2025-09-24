# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 21:59:59 2025

@author: iitm9
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Define the function 'i' as a function of time
def i(t):
    return np.sin(t)  # Example: sinusoidal function

# Define the coupled ODE system
def coupled_ode(t, y, tau):
    i1 = y[0]
    di1_dt = (i(t) - i1) / tau
    return [di1_dt]

# Parameters
tau = 1  # Relaxation time constant
i1_initial = 0.0  # Initial condition for i1
t_span = (0, 10)  # Time range
t_eval = np.linspace(t_span[0], t_span[1], 20)  # Time points for evaluation
dt=t_eval[1]-t_eval[0]
# Solve the ODE using the original method
original_solution = solve_ivp(coupled_ode, t_span, [i1_initial], t_eval=t_eval, args=(tau,), method='RK45')

# Time-step approach
#dt = 0.5  # Time increment for the step-by-step solution
#time_steps = np.arange(t_span[0], t_span[1] + dt, dt)  # Generate time steps

i1_timestep = [i1_initial]  # Initialize i1 with the initial condition
i2_timestep = [i1_initial]  # Initialize i1 with the initial condition

# Iterative computation for the time-step approach
for t in t_eval[:-1]:  # Iterate through the time steps
    previous_i1 = i1_timestep[-1]
    previous_i2 = i2_timestep[-1]
    new_i = i(t + dt)  # Evaluate i(t) for the next time step
    new_i1 = (previous_i1 + new_i * dt / tau)/(1+dt/tau) # Calculate i1 for the next step
    new_i2 = previous_i2 + (new_i - previous_i2) * dt / tau
    i1_timestep.append(new_i1)
    i2_timestep.append(new_i2)

# Convert i1_timestep to a numpy array for plotting
i1_timestep = np.array(i1_timestep)
i2_timestep = np.array(i2_timestep)

# Plot both solutions for comparison
plt.plot(original_solution.t, original_solution.y[0], label='Original Solution (solve_ivp)', linestyle='--')
plt.plot(t_eval, i1_timestep, label='Time-Step Solution', marker='o', markersize=4, alpha=0.7)
plt.plot(t_eval, i2_timestep, label='Time-Step Solution old', marker='x', markersize=4, alpha=0.7)
plt.plot(original_solution.t, i(original_solution.t), label='i(t) (Input Function)', linestyle=':')
plt.xlabel('Time (t)')
plt.ylabel('i1(t)')
plt.title('Comparison of Original and Time-Step Solutions')
plt.legend()
plt.grid()
plt.show()