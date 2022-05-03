# -*- coding: utf-8 -*-
"""
Arbeidskrav 4

@authors: Sigmund Grimsrud, sigmund.grimsrud@nmbu.no
          Muhammed Idris Omar, mohammed.idris.omar@nmbu.no 
          Liibaan Hassan Osman, liibaan.hassan.osman@nmbu.no
"""
import numpy as np
import matplotlib.pyplot as plt

A = 2000         # [m^2]
# delta_Fin = 1.0  # [m^3/s]
# delta_Hmax = 0.5 # [m]

Tc = 1000   # [s] Tc = A*dh_max/dF_in
Kc = -2     # [m^3/s^2] Kc= -A/Tc
Kp = Kc
Ti = 2000   # [s] Ti = 2*Tc
u_max = 100
u_min = 0
u_i_km1 = 3

f_out = u_i_km1

# Height Set Point
h_sp = 1   # [m]
h_min = 0  # [m]
h_max = 2  # [m]
h_k = h_sp


t_step = 1       # [s]
t_start = 0      # [s]
t_stop = 10000   # [s]
N_sim = int((t_stop - t_start)/t_step) + 1  # Number of time-steps

"""
output skal være en sinusbølge med:
    Amplitude 1 [m^3/s]
    Period 1200 [s]
"""
period = 1200  # [s]
amplitude = 1  # [m^3/s]

#%% Preallocation of arrays for plotting:
t_array = np.zeros(N_sim)
h_array = np.zeros(N_sim)
f_in_array = np.zeros(N_sim)
f_out_array = np.zeros(N_sim)
delay_array = np.zeros(Tc)

#%% Sinusoid F_in
radiant_list = np.linspace(0, 2*np.pi*N_sim/period, N_sim)
sin_array = amplitude * np.sin(radiant_list) + 3


# %% Function defining a PI controller:

def pi_contr(h_sp, h_k, u_i_km1, Kp, Ti, u_min, u_max, t_step):
    
    e_k = h_sp - h_k  # Control error
    u_p_k = Kp*e_k  # P term
    u_i_k = u_i_km1 + (Kp*t_step/Ti)*e_k  # PI term
    u_i_min = u_min
    u_i_max = u_max
    u_i_k = np.clip(u_i_k, u_i_min, u_i_max)  # Limit ui
    u_k = u_p_k + u_i_k  # Total control signal
    u_k = np.clip(u_k, u_min, u_max)  # Limit of control
    
    return (u_k, u_i_k)

#%% Simulation

for k in range(0, N_sim):
    # State limitation:
    h_k = np.clip(h_k, h_min, h_max)
    t_k = k*t_step  # Current time [s]
    
    f_in_k = sin_array[k]
    
    # PI
    u_k, u_i_k = pi_contr(h_sp, h_k, u_i_km1, Kp, Ti, u_min, u_max, t_step)
    
    f_out_k = u_k
    
    # Hight Time-derivative:
    dh_dt_k = ((f_in_k/A))*(f_in_k - f_out_k)
    
    # State updates using the Euler method:
    h_kp1 = h_k + dh_dt_k*t_step
    
    # Arrays for plotting:
    t_array[k] = t_k
    h_array[k] = h_k
    f_in_array[k] = f_in_k
    f_out_array[k] = f_out_k

    # Time shift:
    u_i_km1 = u_i_k
    h_k = h_kp1
    
#%% Plotting
plt.close('all')  # Closes all figures before plotting
plt.figure(num=1, figsize=(12, 9))

h_sp_array = np.zeros(N_sim)
for k in range(0, N_sim):
    h_sp_array[k] = h_sp
plt.subplot(2, 1, 1)
plt.plot(t_array, h_array, 'c')
plt.plot(t_array, h_sp_array, 'g')
plt.grid()
plt.ylim(h_min, h_max)
plt.xlabel('t [s]')
plt.ylabel('[m]')
plt.legend(labels=('Level h', 'h set point'))

plt.subplot(2, 1, 2)
plt.plot(t_array, f_in_array, 'c')
plt.plot(t_array, f_out_array, 'y')
plt.grid()
plt.ylim(0, 6)
plt.xlabel('t [s]')
plt.ylabel('[m^3/s]')
plt.legend(labels=('Inflow f_in', 'Outflow f_out = u'))

plt.show()
