import numpy as np
import matplotlib.pyplot as plt

from SEIR import SEIR


# DEFINING PARAMETERS

# Time in days
t_max = 100
dt = 0.1
t = np.linspace(0, t_max, int(t_max/dt) + 1)

# Initial Values
N = 10000
S0 = 1 - 1/N
E0 = 1/N
I0 = 0
R0 = 0

# Model Parameters: COVID_19

# 1/t_incupation. t_incubation = 5days
alpha = 0.2

# the invearse of the mean infectious period (1/t_infectious).
# t_infectious = 2 days
gamma = 0.5

# average contact rate in the population,
# R0 = beta/gamma. R0 represents how quickly the disease spreads.
# For COVID-19 this value is around 3.5.
beta = 1.75

# FITTING MODEL
model = SEIR(beta, alpha, gamma, S0, E0, I0, R0)
U = model.compile(t)


# Model Without Social Distancing
plt.figure()
S_t = plt.plot(t, U[:, 0])
E_t = plt.plot(t, U[:, 1])
I_t = plt.plot(t, U[:, 2])
R_t = plt.plot(t, U[:, 3])
plt.legend(['Susceptible', 'Exposed', 'Infected', 'Recovered'])
plt.xlabel('Days')
plt.ylabel('Fraction of Population')

plt.title('SEIR model for COVID-19')
plt.show()
