import numpy as np


class SEIR:

    def __init__(self, beta, alpha, gamma, S0, E0, I0, R0):
        """Initalize initial values and model parameters

        Args:
            beta (function): parameter of the model SEIR
            alpha (function): parameter of the model SEIR
            gamma (function): parameter of the model SEIR
            S0 (int): initial suceptible population
            E0 (int): initial exposed population
            I0 (int): initial infected population
            R0 (int): initial recovered/dead population

        Returns:
            None
        """

        self.initial_conditions = [S0, E0, I0, R0]
        # beta, alpha and gamma can, maybe, vary over time.
        # If a constant value is passed, it will be transforme to a function
        if (isinstance(beta, (int, float))):
            self.beta = lambda t: beta
        else:
            self.beta = beta
        if (isinstance(alpha, (int, float))):
            self.alpha = lambda t: alpha
        else:
            self.alpha = alpha
        if (isinstance(gamma, (int, float))):
            self.gamma = lambda t: gamma
        else:
            self.gamma = gamma

    def compile(self, timestemp):
        """Apply the SEIR model to the given timestemp

        Args:
            timestemp (array): timestemp array

        Returns:
            u (matrix): matrix containing results of each poplutation through
            time. The matrix follows the format:
            (len(timestemp) X [S, E, I, R])
        """
        # In order to solve the differential equations the semi-implicit
        # Euler method will be use.

        n = len(timestemp)

        # Matrix of values: n X number_of_equations
        # timestemp x [S, E, I , R]
        u = np.zeros((n, len(self.initial_conditions)))
        u[0, :] = self.initial_conditions

        # Timestemp resolution
        dt = timestemp[1] - timestemp[0]

        for t in range(n - 1):
            # S{t + 1} = S{t} - beta{t}*S{t}*I{t}*dt
            u[t + 1, 0] = u[t, 0] - (self.beta(t)*u[t, 0]*u[t, 2]) * dt
            # E{t + 1} = E{t} + (beta{t}*S{t}*I{t} - alpha(t)*E{t})*dt
            u[t + 1, 1] = u[t, 1] + (self.beta(t)*u[t, 0]*u[t, 2]
                                     - self.alpha(t)*u[t, 1]) * dt
            # I{t + 1} = I{t} + (alpha(t)*E{t} - gamma(t)*I{t})*dt
            u[t + 1, 2] = u[t, 2] + (self.alpha(t)*u[t, 1]
                                     - self.gamma(t)*u[t, 2]) * dt
            # R{t + 1} = R{t} + gamma(t)*I{t}*dt
            u[t + 1, 3] = u[t, 3] + self.gamma(t)*u[t, 2] * dt

        return u
