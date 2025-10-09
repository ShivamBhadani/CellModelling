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
        self.Innovation=[]
        self.Vcmin=3.5  #at 25C
        self.Vcmax=4.25

        """
        Model Udate part-1: Update Qc Cell capacity
        RLS estimator for battery capacity Qc.
        :param initial_Qc: Initial guess for capacity in mAh
        :param QcRLSmemory: Forgetting factor (close to 1 for slow aging)
        :param sigma_Q: Std deviation of charge measurement noise (mAh)

        self.Qc = self.variable_stats['Qc']['mean']  # Initial estimate of capacity
        self.QcUncertainity = 16*self.variable_stats['Qc']['var']  # Initial uncertainty (scalar)
        self.QcRLSmemory = 0.995  # may be optimized based on how fast/slow Qc ages/changes
        self.Qvariance = 0  # Variance of Q(t) measurement
        self.coulombCount=0
        self.absCoulombCount=0
        self.x_est_prev=self.x_est
        self.coulombCount_prev=0
        """
        # RLS-based capacity estimation
        self.Qc_charge = self.variable_stats['Qc']['mean']
        self.Qc_discharge = self.variable_stats['Qc']['mean']
        self.P_charge = 1e6
        self.P_discharge = 1e6
        self.QcRLSmemory = 0.99
        self.sigma_Q = 5.0  # mAh noise in coulomb count
        self.SoC_variance_threshold = 0.002
        self.SoC_start = None
        self.Q_accumulated = 0.0
        self.mode = None  # 'charge' or 'discharge'
        self.rls_active = False
        self.min_delta_soc_factor = 5.0  # Minimum delta SoC = factor × sqrt(P)
        self.min_delta_soc_forQcUpdate=0.01
        self.qcUpdateMode= 1 # use min_delta, if 0 use the delta_soc_factor
        self.chargeQcUpdateCount=0
        self.dischargeQcUpdateCount=0
        self.charge_soc=[]
        self.discharge_soc=[]
        self.initialCalibrationComplete=False
        
    def checkFuse(self, discharge):
        #dt= time - self.time[-1]
        #discharge=i*dt
        if discharge==0:
            return True
        if discharge>0:
            # We are discharging
            dSoC=discharge/(self.variable_stats['eta']['mean']*self.variable_stats['Qc']['mean'])
                        
        else:
            # We are charging
            dSoC=discharge/self.variable_stats['Qc']['mean']

        if (self.x_est-dSoC)>=0.99 and discharge<0:
            self.charge=False
            #print("Cannot charge anymore battery full. SoC",self.SoC[-1]," SoH=",self.SoH[-1]," time=",time)
            return False
        elif (self.x_est-dSoC)<=0.01 and discharge>0:
            #print("Cannot discharge anymore battery empty. SoC",self.SoC[-1]," SoH=",self.SoH[-1]," time=",time)
            self.discharge=False
            return False
        else:
            self.charge=True
            self.discharge=True
            return True


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

    def update(self, V_k, Qaccumulated):
         sigma_points, weights = self._generate_sigma_points(self.x_pred, self.P_pred)
         #weights_mean = np.mean(weights)
         #weights_cov = weights_mean.copy()
         #Variance is v_noise+d_ocv+dI*z+I*dZ, since dI/I is much smaller than dZ/Z, ignore dI
         deltaV=self.cellESR.calculateESR(self.variable_stats['I']['mean'],self.variable_stats['dt']['mean'])*self.variable_stats['I']['mean']
         #print("deltaV=",deltaV," time=",self.cellESR.time)
         self.R=self.vNoise**2+self.dOCV**2+(self.dZ**2)*abs(deltaV)
         predicted_measurements = [self.socEstimator.get_ocv(100*(1-sp),self.variable_stats['T']['mean'])+deltaV for sp in sigma_points]
         #print("x0=",sigma_points[0]," OCV=",self.socEstimator.get_ocv(100*(1-sigma_points[0]),self.variable_stats['T']['mean']), " DeltaV=",deltaV, "Vpred=",predicted_measurements)
         
         
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
         #print("slope=",slope)
         
 
         K = slope*(self.P_pred/(self.P_pred+self.R))
         
         #K = slope
         self.x_est=min(1,max(self.x_pred + K * (V_k - z_pred),0))
         self.P = max(self.P_pred - K * P_zz * K, self.Q)
                
         self.Innovation.append([self.x_pred,self.x_est,slope,K,V_k,z_pred])    
         self.update_capacity_estimate(Qaccumulated)
             
         #print("K=",K," Slope=",slope)
         #print("v_k=",V_k," Zpred=",z_pred)
         #print("post update x_est=",self.x_est," x_pred=",self.x_pred)
    def get_estimate(self):
        return self.x_est
    
    def start_cycle(self, mode, Q):
        self.SoC_start = self.x_est
        self.Q_offset = Q
        self.mode = mode  # 'charge' or 'discharge'

        self.rls_active = True
        print("started Qc RLS soc=",self.SoC_start," Qoffset=",self.Q_offset)
    
    def reset_cycle(self):
        self.SoC_start = None
        self.rls_active = False
        self.updateCount=None
        print("stopped Qc RLS")
        
    def update_capacity_estimate(self, Qaccumulated):
        # Check if SoC variance is low enough to trust
   
        if self.rls_active and self.SoC_start is not None:
            if ((self.mode=="charge") and (self.variable_stats['I']['mean']>0)) or ((self.mode=="discharge") and (self.variable_stats['I']['mean']<0)):
                self.reset_cycle()
            else:            
                delta_SoC = self.x_est - self.SoC_start
                if self.qcUpdateMode:
                    dsocThreshold = self.min_delta_soc_forQcUpdate
                else:
                    dsocThreshold = self.min_delta_soc_factor * np.sqrt(self.P)
                #print("dsocthreshold=",dsocThreshold)
                if abs(delta_SoC) > dsocThreshold:  # Minimum observable change
                    phi = delta_SoC
                    y = Qaccumulated-self.Q_offset
                    if self.mode == 'charge':
                        denom = self.P_charge * phi**2 + self.sigma_Q**2
                        K = (self.P_charge * phi) / denom
                        error = y - phi * self.Qc_charge
                        self.Qc_charge += K * error
                        self.P_charge = (self.P_charge - K * phi * self.P_charge) / self.QcRLSmemory
                        self.chargeQcUpdateCount+=1
                        self.charge_soc.append([phi, y])
                        #print(self.mode,"SoC=",self.x_est,"dQ=",Qaccumulated," Qc=",self.Qc_discharge)
                        #self.variable_stats['eta']['mean']=-self.Qc_charge
                        #self.variable_stats['Qc']['var']=-self.P_charge                        print(self.mode,"SoC=",self.x_est," Qc=",self.Qc_discharge)
                    elif self.mode == 'discharge':
                        denom = self.P_discharge * phi**2 + self.sigma_Q**2
                        K = (self.P_discharge * phi) / denom
                        error = y - phi * self.Qc_discharge
                        self.Qc_discharge += K * error
                        self.P_discharge = (self.P_discharge - K * phi * self.P_discharge) / self.QcRLSmemory
                        self.dischargeQcUpdateCount+=1
                        self.discharge_soc.append([phi, y])
                        #print(self.mode,"SoC=",self.x_est,"dQ=",Qaccumulated," Qc=",self.Qc_discharge)
                        #self.variable_stats['Qc']['mean']=-self.Qc_discharge
                        #self.variable_stats['Qc']['var']=self.P_discharge
                        #print(self.mode,"SoC=",self.x_est," Qc=",self.Qc_discharge)
                    self.SoC_start=self.x_est
                    self.Q_offset=Qaccumulated

                    self.update_capacity_estimates()

    def update_capacity_estimates(self):
        efficiency =  self.Qc_discharge / self.Qc_charge if self.Qc_charge != 0 else None
        if (self.dischargeQcUpdateCount>3):
            self.variable_stats['Qc']['mean']=-self.Qc_discharge
            self.variable_stats['Qc']['var']=self.P_discharge
            print("Updated Qc=",self.variable_stats['Qc']['mean'])
            if (self.chargeQcUpdateCount>3):
                self.variable_stats['eta']['mean']=efficiency
                print("Updated Eta=",self.variable_stats['eta']['mean'])
        return {
            "Qc_charge": self.Qc_charge,
            "Qc_discharge": self.Qc_discharge,
            "efficiency": efficiency
        }

    def start_calibration_cycle(self, V, T, mode):
        self.calibration_soc_start = (100-self.socEstimator.get_dod(V, T))/100
        self.calibration_Q = 0
        self.calibration_cycle_data = []
        self.calibration_mode = mode
        self.calibration_active = True
        self.slopeCheck_soc= self.calibration_soc_start
        self.slopeCheck_ocv= V
        self.scale_factor=1
        self.slopes=[]
 
                    
    def record_calibration_point(self, I, V, T, dt, slope_threshold=2, min_soc_span=0.01):
        """
        Record a calibration observation and check if slope-based anchor condition is met.
        calibrationping completes only when dOCV/dSoC exceeds threshold.
        """
        dQ=I*dt
        self.calibration_Q +=dQ
        self.deltaV = self.cellESR.calculateDeltaV(I,dt)
        V_ocv = V - self.deltaV
        SoC = self.calibration_soc_start - (self.calibration_Q/self.variable_stats['Qc']['mean'])
        self.calibration_cycle_data.append((SoC, V_ocv, T))
        #if V_ocv>self.Vcmax:
            #print("Reached Max Voltage")
    
        if abs(self.slopeCheck_soc-SoC)>min_soc_span:
            # Sort by SoC
            slope=(self.slopeCheck_ocv-V_ocv)/(self.slopeCheck_soc-SoC)
            self.slopes.append(slope)
            self.slopeCheck_ocv=V_ocv
            self.slopeCheck_soc=SoC
            soc_now=max(min(1,(100-self.socEstimator.get_dod(V_ocv,T))/100),0)
            if len(self.slopes)>4:
                if ((np.mean(self.slopes[-3:])>slope_threshold) and abs(self.calibration_soc_start-soc_now)>0.5)or(self.calibration_mode=="discharge" and V_ocv<self.Vcmin) or (self.calibration_mode=="charge" and V_ocv>self.Vcmax):               
                    self.scale_factor = (self.calibration_soc_start - soc_now) / (self.calibration_soc_start - SoC)
                    print(self.calibration_mode," rescaling factor=",self.scale_factor)
                    rescaled_data = [
                        (self.calibration_soc_start + (s - self.calibration_soc_start) * self.scale_factor, v, t)
                        for s, v, t in self.calibration_cycle_data
                    ]
                    self.calibration_cycle_data=rescaled_data
                    self.socEstimator.update_from_observations(rescaled_data)
                    #self.socEstimator.update_ocv_grid_from_calibration(rescaled_data)
                    self.x_est=soc_now
                    print("Previous Qc=",self.variable_stats['Qc']['mean'])
                    self.variable_stats['Qc']['mean']=abs((self.calibration_Q)/(self.calibration_soc_start - soc_now)) 
                    print("Updated Qc=",self.variable_stats['Qc']['mean'])
                    self.calibration_active = False
                    print(f"✅ calibrationping complete: max slope = {slope:.3f} V/SoC, SoC swing = {(soc_now-self.calibration_soc_start):.2f}, scale = {self.scale_factor:.2f}, points = {len(rescaled_data)}")
                