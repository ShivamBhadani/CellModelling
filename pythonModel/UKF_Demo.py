# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 12:21:39 2025

@author: iitm9
"""

import numpy as np


class Cell:
    def __init__(self, capacity, eta, Id):
        self.capacity = capacity
        self.eta = eta
        self.Id = Id
        self.SoC = [100.0]
        self.voltage = []
        self.current = []
        self.Temperature = [25.0]

    def step(self, I_k, delta_t):
        I_eff = self.eta * I_k if I_k > 0 else I_k
        soc_new = self.SoC[-1] + (I_eff * delta_t) / self.capacity - (self.Id * delta_t) / self.capacity
        self.SoC.append(soc_new)
        self.current.append(I_k)
        self.voltage.append(self.get_voltage(soc_new))
        self.Temperature.append(25.0)

    def get_voltage(self, soc):
        return 3.0 + 0.01 * soc  # Simplified OCV model

    def get_noisy_observation(self, voltage_noise_std=0.01, current_noise_std=0.05):
        noisy_voltage = self.voltage[-1] + np.random.normal(0, voltage_noise_std)
        noisy_current = self.current[-1] + np.random.normal(0, current_noise_std)
        temperature = self.Temperature[-1]
        return noisy_voltage, noisy_current, temperature
    
class SoCUKFEstimator:
    def __init__(self, initial_soc, capacity, eta, Id, Q, R):
        self.x_est = initial_soc
        self.P = 1.0
        self.C = capacity
        self.eta = eta
        self.Id = Id
        self.Q = Q
        self.R = R
        self.alpha = 1e-1
        self.beta = 2
        self.kappa = 0

    def _generate_sigma_points(self, x, P):
        n = 1
        lambda_ = self.alpha**2 * (n + self.kappa) - n
        sigma_points = [x]
        sqrt_P = np.sqrt((n + lambda_) * P)
        sigma_points.append(x + sqrt_P)
        sigma_points.append(x - sqrt_P)
        return np.array(sigma_points), lambda_

    def predict(self, I_k, delta_t):
        I_eff = self.eta * I_k if I_k > 0 else I_k
        self.x_pred = self.x_est + (I_eff * delta_t) / self.C - (self.Id * delta_t) / self.C
        self.P_pred = self.P + self.Q

    def update(self, V_k):
        sigma_points, lambda_ = self._generate_sigma_points(self.x_pred, self.P_pred)
        weights_mean = [lambda_ / (1 + lambda_), 0.5 / (1 + lambda_), 0.5 / (1 + lambda_)]
        weights_cov = weights_mean.copy()

        predicted_measurements = [3.0 + 0.01 * sp for sp in sigma_points]
        z_pred = np.dot(weights_mean, predicted_measurements)
        P_zz = self.R + sum(weights_cov[i] * (predicted_measurements[i] - z_pred)**2 for i in range(3))
        P_xz = sum(weights_cov[i] * (sigma_points[i] - self.x_pred) * (predicted_measurements[i] - z_pred) for i in range(3))

        K = P_xz / P_zz
        self.x_est = self.x_pred + K * (V_k - z_pred)
        self.P = self.P_pred - K * P_zz * K

    def get_estimate(self):
        return self.x_est
    

#R0=2e-3 #ohms
#RZlist=[[1e-3,10,0],[1.5e-3,1e4,0],[8,5e-7,1]]


import matplotlib.pyplot as plt

cell = Cell(capacity=100.0, eta=0.97, Id=1e-6)
estimator = SoCUKFEstimator(initial_soc=50.0, capacity=100.0, eta=0.99, Id=1e-6, Q=0.1, R=0.05)

soc_true = []
soc_est = []
P_values = []

for t in range(1000):
    I_k = -1.0  # Discharge
    delta_t = 1.0
    cell.step(I_k, delta_t)
    V_k, _, _ = cell.get_noisy_observation()
    estimator.predict(I_k, delta_t)
    estimator.update(V_k)
    soc_true.append(cell.SoC[-1])
    soc_est.append(estimator.get_estimate())
    P_values.append(estimator.P)

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot(soc_true, label='True SoC')
plt.plot(soc_est, label='Estimated SoC')
plt.legend()
plt.title('SoC Estimation')

plt.subplot(1, 2, 2)
plt.plot(P_values)
plt.title('Covariance P Over Time')
plt.xlabel('Time Step')
plt.tight_layout()
plt.show()
    
