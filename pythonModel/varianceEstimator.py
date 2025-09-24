# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 19:41:45 2025

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


# Define variable statistics
battery_stats = {
    'dt':   {'mean': 2.0, 'var': 0.1},
    'I':    {'mean': 2.5, 'var': 0.3},
    'p':    {'mean': 1, 'var': 0.0}
}

cell_stats = {
    'R0':   {'mean':0.1, 'var':0.005},
    'Qc':   {'mean': 3.0, 'var': 0.2},
    'eta':  {'mean': 1.5, 'var': 0.05},
    'etaD': {'mean': 1.5, 'var': 0.05},
    'Id':   {'mean': 1.2, 'var': 0.1}
}
variable_stats=battery_stats|cell_stats
stats=StatisticsEstimator(variable_stats)

# Define the parsed expression
# Expression: (dt / Qc) * (eta * I - Id)
common_factors = ['dt', 'Qc^-1']
inner_term = [
    ['eta', 'I'],
    ['-Id'] ]

final_mean, final_variance = stats(common_factors, inner_term)
# Output
print(f"Estimated mean: {final_mean:.4f} and Estimated variance: {final_variance:.4f} ")