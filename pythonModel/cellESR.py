# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 21:17:26 2025

@author: iitm9
"""

import pandas as pd
from scipy.interpolate import interp1d
import numpy as np

class cellESR:
    def __init__(self, R0, RZlist):
        self.R0=R0
        self.RTlist=[]
        self.previousIx=[]
        self.esrDC=self.R0
        for RZ in RZlist:
            if RZ[2]==0:
                self.RTlist.append([RZ[0],RZ[0]*RZ[1],0,0])
                self.esrDC+=RZ[0] 
            else:
                self.RTlist.append([RZ[0],RZ[1]/RZ[0],1,0])
        self.esr=R0
        self.time=0
            
    def ix_next(self,ix_previous, tau, i_new, dt):
        # Using the equation: i1_next = previous_i1 + (new_i - previous_i1) * dt / tau
       ix_next = (ix_previous + i_new * dt / tau)/(1+dt/tau)
       if abs(ix_next)>abs(i_new):
           ix_next=i_new
       if abs(ix_next/i_new)>1:
           print("inext=",ix_next," i_new=",i_new)
       return ix_next

    def calculateESR(self,i_new, dt):
        self.time+=dt
        ESRout=self.R0
        if i_new!=0:
            for index,[r,tau,ind,ix_previous] in enumerate(self.RTlist):
                self.RTlist[index][3]=self.ix_next(ix_previous,tau,i_new,dt)
                if ind:
                    vz=r*(i_new-self.RTlist[index][3])
                else:
                    vz=r*self.RTlist[index][3]
                ESRout+=(vz/i_new)
        return ESRout

    def calculateDeltaV(self,i_new, dt):
        self.time+=dt
        ESRout=self.R0
        dV=i_new*ESRout

        for index,[r,tau,ind,ix_previous] in enumerate(self.RTlist):
            self.RTlist[index][3]=self.ix_next(ix_previous,tau,i_new,dt)
            if ind:
                dV+=r*(i_new-self.RTlist[index][3])
            else:
                dV+=r*self.RTlist[index][3]
        return dV
    
class LookupTableEstimator: #TODO older implementation is a lot faster
    def __init__(self, file_path):
        # Load the lookup table once during initialization
        self.lookup_table = pd.read_csv(file_path)
        # Extract column headers
        headers = self.lookup_table.columns.tolist()
        # First column: input values
        self.inputs = self.lookup_table[headers[0]].values
        # Remaining columns: output data
        self.temperatures = [int(header.split(" ")[-1][:-1]) for header in headers[1:]]
        self.data = self.lookup_table.iloc[:, 1:].values

    def estimateOCV(self, input_value, temperature):
        input_func = interp1d(self.inputs, self.data, axis=0, fill_value="extrapolate")
        outputs_at_temperatures = input_func(input_value)
        temp_func = interp1d(self.temperatures, outputs_at_temperatures, fill_value="extrapolate")
        return temp_func(temperature)

# class LookupTableEstimator:
#     def __init__(self, file_path):
#         self.lookup_table = pd.read_csv(file_path)

#         headers = self.lookup_table.columns.tolist()
#         self.DoD_vector = self.lookup_table[headers[0]].values
#         self.temperature_vector = [int(header.split(" ")[-1][:-1]) for header in headers[1:]]
#         self.OCV_grid = self.lookup_table.iloc[:, 1:].values

#     def estimateOCV(self, dod, temperature):
#         # Interpolate OCV at each DoD for the given temperature
#         ocv_at_dod = []
#         for row in self.OCV_grid:
#             temp_interp = interp1d(self.temperature_vector, row, bounds_error=False, fill_value="extrapolate")
#             ocv_at_dod.append(temp_interp(temperature))

#         # Interpolate OCV from DoD
#         dod_interp = interp1d(self.DoD_vector, ocv_at_dod, bounds_error=False, fill_value="extrapolate")
#         return float(dod_interp(dod))

#     def get_dod(self, ocv, temperature):
#         # Interpolate OCV at each DoD for the given temperature
#         ocv_at_dod = []
#         for row in self.OCV_grid:
#             temp_interp = interp1d(self.temperature_vector, row, bounds_error=False, fill_value="extrapolate")
#             ocv_at_dod.append(temp_interp(temperature))

#         # Interpolate DoD from OCV
#         ocv_interp = interp1d(ocv_at_dod, self.DoD_vector, bounds_error=False, fill_value="extrapolate")
#         return float(ocv_interp(ocv))