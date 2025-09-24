# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 19:55:21 2025

@author: iitm9
"""

import numpy as np

class StatisticsEstimator:
    def __init__(self, variable_stats):
        self.stats = variable_stats


    def __call__(self, common_factors, inner_term):
        mean_c = self.product_mean(common_factors)
        var_c  = self.product_variance(common_factors)
        mean_inner = sum(self.product_mean(t) for t in inner_term)
        var_inner  = sum(self.product_variance(t) for t in inner_term)
        return mean_c*mean_inner, (mean_c**2 * var_inner) + (var_c * mean_inner**2)

    # Utility functions to extract mean and variance
    def get_mean(self, var):
        if '^' in var:
            base, exp = var.split('^')
            return self.stats[base]['mean'] ** float(exp)
        elif var.startswith('-'):
            return -self.get_mean(var[1:])
        else:
            return self.stats[var]['mean']

    def get_var(self,var):
        if '^' in var:
            base, power = var.split('^')
            m = self.stats[base]['mean']
            v = self.stats[base]['var']
            exp = int(power)
            #print(var," =",np.sign(exp)*exp*((m)**(exp-1))*v)
            return ((np.sign(exp)*exp*((m)**(exp-1)))**2)*v
        elif var.startswith('-'):
            return self.get_var(var[1:])
        else:
            return self.stats[var]['var']

    # Mean and variance of a product
    def product_mean(self, vars):
        result = 1
        for v in vars:
            result *= self.get_mean(v)
        return result
    
    def product_variance(self, vars):
        total = 0
        for i, v in enumerate(vars):
            v_i = self.get_var(v)
            print(v," get_var=",v_i)
            other_means_squared = 1
            for j, u in enumerate(vars):
                if i != j:
                    other_means_squared *= self.get_mean(u)**2                    
            total += v_i * other_means_squared
        #ignore higher order variance terms
        return total


class SoCUKFEstimator:
    def __init__(self, battery_stats, cell_stats, RZlist, ocvLookup):
        #Initialize Model parameters: initial_soc, capacity, eta, Id, Q, R
        self.battery_stats=battery_stats #passed by reference
        self.cell_stats = cell_stats.copy()
        self.variable_stats=self.battery_stats|self.cell_stats
        self.stats=StatisticsEstimator(self.variable_stats)
        self.x_est = 0.9  # initial_soc
        self.cellESR=cellESR(cell_stats['R0']['mean'],RZlist)
        self.socEstimator=LookupTableEstimator(ocvLookup)
        self.P = 1.0  # initial variance on soc estimate
        #self.C = capacity
        #self.eta = eta
        #self.Id = Id
        self.Q = 0.1  #variance in SoC update lets calculate this
        self.R = R  #variance in voltage/ocv measurements from lookup and noise/errors, lets 
        self.lambda_=0.1


    def _generate_sigma_points(self, x, P):
        n = 1
        #lambda_ = self.alpha**2 * (n + self.kappa) - n
        sigma_points = [x]
        sqrt_P = np.sqrt((n + self.lambda_) * P)
        sigma_points.append(min(x + sqrt_P,1))
        sigma_points.append(max(x - sqrt_P,0))
        weights=[y/x for y in sigma_points]
        return np.array(sigma_points), np.array(weights)

    def predict(self, I_k, delta_t):
        I_eff = self.eta * I_k if I_k > 0 else I_k
        
        self.x_pred = self.x_est + (I_eff * delta_t) / self.C - (self.Id * delta_t) / self.C
        self.P_pred = self.P + self.Q

    def update(self, V_k):
        sigma_points, weights = self._generate_sigma_points(self.x_pred, self.P_pred)
        weights_mean = np.mean(weights)
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