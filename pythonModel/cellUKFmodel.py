# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 19:55:21 2025

@author: iitm9
"""

import numpy as np
from cellESR import cellESR, LookupTableEstimator

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
            #print(v," get_var=",v_i)
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
        self.vNoise=0.001 #1mV
        self.dOCV= 0.001 #1mV
        self.dZ=0.05 #5%
        #self.C = capacity
        #self.eta = eta
        #self.Id = Id
        self.x_pred=0.9
        self.common_factors = ['dt', 'Qc^-1']
        self.inner_term = [
            ['eta^-1', 'I'],
            ['-Id'] ]


        self.Q = 0.1  #variance in SoC update lets calculate this
        #self.R = R  #variance in voltage/ocv measurements from lookup and noise/errors, lets 
        self.lambda_=0.1


    def _generate_sigma_points(self, x, P):
        n = 1
        #lambda_ = self.alpha**2 * (n + self.kappa) - n
        sigma_points = [x]
        sqrt_P = np.sqrt((self.lambda_) * P)
        sigma_points.append(min(x + sqrt_P,1))
        sigma_points.append(max(x - sqrt_P,0))
        weights=[y/x for y in sigma_points]
        #print("p=",P," sqrt_P=",sqrt_P," sigma points=",sigma_points)
        #print("weights=",weights)
        return np.array(sigma_points), np.array(weights)

    def predict(self):
        #I_eff = self.eta * I_k if I_k > 0 else I_k
        if self.variable_stats['I']['mean']<0:
            self.inner_term = [
                ['eta^-1', 'I'],
                ['-Id'] ]
        else:
            self.inner_term = [
                ['I'],
                ['-Id'] ]
                           
        self.Qcharge, self.Q = self.stats(self.common_factors, self.inner_term)
        self.x_pred = self.x_est - self.Qcharge
        self.P_pred = self.P + self.Q
#        print("dQ=",self.Qcharge," post predict x_pred=",self.x_pred)

    def update(self, V_k):
        
        sigma_points, weights = self._generate_sigma_points(self.x_pred, self.P_pred)
        weights_mean = np.mean(weights)
        #weights_cov = weights_mean.copy()
        #Variance is v_noise+d_ocv+dI*z+I*dZ, since dI/I is much smaller than dZ/Z, ignore dI
        deltaV=self.cellESR.calculateESR(self.variable_stats['I']['mean'],self.variable_stats['dt']['mean'])*self.variable_stats['I']['mean']
        #print("deltaV=",deltaV," time=",self.cellESR.time)
        self.R=self.vNoise**2+self.dOCV**2+(self.dZ**2)*abs(deltaV)
        predicted_measurements = [self.socEstimator.output(100*(1-sp),self.variable_stats['T']['mean'])+deltaV for sp in sigma_points]
        
        if (max(sigma_points)<1) and (min(sigma_points)>0):
            z_pred = np.dot(weights, predicted_measurements)/np.sum(weights)
            P_xz = sum(weights[i] * (sigma_points[i] - self.x_pred) * (predicted_measurements[i] - z_pred) for i in range(3))
            P_zz = sum(weights[i] * (predicted_measurements[i] - z_pred)**2 for i in range(3))
            slope=(P_xz / P_zz)
        else:
            z_pred=predicted_measurements[0]
            P_zz = sum(weights[i] * (predicted_measurements[i] - z_pred)**2 for i in range(3))
            if (max(sigma_points)<1):
                slope=(predicted_measurements[0]-predicted_measurements[2])/(sigma_points[0]-sigma_points[2])                
            else:
                slope=(predicted_measurements[0]-predicted_measurements[1])/(sigma_points[0]-sigma_points[1])

       # print("z_pred=",z_pred," predicted measurements=",predicted_measurements)
        
        

        K = slope*(self.P_pred/(self.P_pred+self.R))
        self.x_est=min(1,max(self.x_pred + K * (V_k - z_pred),0))
        self.P = max(self.P_pred - K * P_zz * K, self.Q)
        
        x_est=self.x_pred + K * (V_k - z_pred)
        P=self.P_pred - K * P_zz * K
        if x_est>0 and x_est<1:
            if self.x_est != x_est:
                print("lets investigate x_est=",x_est," self.x_est=",self.x_est)
            if P!= self.P:
                print("lets investigate p=",P," self.P=",self.P)
                print("self.P_pred=",self.P_pred," K=",K, "Pzz=",P_zz)
        elif x_est>=1:
            if self.x_est !=1:
                print("lets investigate why x_est not clipped x_est=",x_est," self.x_est=",self.x_est)
            #don't update self.P
        else:
            if self.x_est !=0:
                print("lets investigate why x_est not clipped x_est=",x_est," self.x_est=",self.x_est)
            
            
            
        #print("K=",K," Slope=",slope)
        #print("v_k=",V_k," Zpred=",z_pred)
        #print("post update x_est=",self.x_est," x_pred=",self.x_pred)
    def get_estimate(self):
        return self.x_est