# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 19:54:10 2025

@author: iitm9
"""



import numpy as np
from batteryCell import Cell
from cellUKFmodel import SoCUKFEstimator
import pandas as pd
import matplotlib.pyplot as plt


def fuseCheck(battery, current):
    fuse=True
    for cell in battery:
        fuse = fuse and cell.checkFuse(time,current)
    return fuse
    

#initial_soc, capacity, eta, Id, Q, R

battery_stats = {
    'dt':   {'mean': 1, 'var': 1e-8},
    'I':    {'mean': 0, 'var': 0}
# period of IV measurements, jitter in IV mesurements
#I and (0.0005*I)**2
}
    
cell_stats = {
    'Qc':   {'mean': 65*3600, 'var': 0.1*0.1*65},
    'eta':  {'mean': 0.98, 'var': 0.01*0.01},
    'Id':   {'mean': 1e-3, 'var': 0.25e-6},
    'T' :   {'mean': 25, 'var': 0.01},
    'R0' :   {'mean': 0.002 , 'var': 0.25e-6},    
}
R0=cell_stats['R0']['mean'] #ohms
RZlist=[[1e-3,10,0],[1.5e-3,1e4,0],[8,5e-7,1]]
ocvLookup="Sample_OCV_DoD.csv"

dCurrent=0.0005 # sigma for current measurement
dVoltage=5/64000  # 15bit enob for 1 sigma
dTemperature=0.1 # 0.1deg sigma for temperature measurement
cell_stats['T']['var']=dTemperature**2


#corners=["nom","low","high","mc",'mc',"mc",'mc',"mc"]
corners=['low']
battery=[]
battery_model=[]
soc_true=[]
soc_est=[]
P_values=[]
I_values=[]
for i, corner in enumerate(corners):
    battery.append(Cell("Li_Ion", R0, RZlist, corner=corner))
    battery_model.append(SoCUKFEstimator(battery_stats, cell_stats, RZlist, ocvLookup))
    soc_true.append([])
    soc_est.append([])
    P_values.append([])
    I_values.append([])
           
input_csv = 'test_sequence.csv'
#input_csv = 'la92shortdds.csv'
data = pd.read_csv(input_csv)
offset=0
time_measured = data['Time'].values[offset:]
current = data['Current'].values[offset:]
temperature= data['Temperature'].values[offset:]


max_loop=10000
loop_count=0
charge=0
previous_t=0
t_model=0
previous_fuse=True
for i, time in enumerate(time_measured):
    if i!=0:
        start_t=time_measured[i-1]
        for nclk in range(int((time-previous_t)/battery_stats['dt']['mean'])):
            loop_count=loop_count+1
            if loop_count>max_loop:
                break            
            t=time_measured[i-1] + nclk*battery_stats['dt']['mean']+np.random.normal(loc=0.0, scale=np.sqrt(battery_stats['dt']['var']), size=None)
            current_now=current[i-1]+(current[i]-current[i-1])*nclk*battery_stats['dt']['mean']/(time-previous_t)
            current_measured=current_now*(1+np.random.normal(loc=0.0, scale=dCurrent, size=None))
            if not(fuseCheck(battery,current_now)):
                #print("protection fuse blown, run simulation with zero current since protection FET is OFF, current=", current_now, " first cell SoC=", battery[0].SoC[-1])
                current_now=0
                current_measured=0
                if previous_fuse:
                    print("fuse tripped at time=",t,"current=", current_now, " first cell SoC=", battery[0].SoC[-1] )
                    previous_fuse=False
                break
            else:
                if not(previous_fuse):
                    print("fuse recovered at time=",t,"current=", current_now, " first cell SoC=", battery[0].SoC[-1])           
                previous_fuse=True
                
            temperature_now=temperature[i-1]+(temperature[i]-temperature[i-1])*nclk*battery_stats['dt']['mean']/(time-previous_t)
            temperature_measured=temperature_now+np.random.normal(loc=0.0, scale=dTemperature, size=None)         
            for j, (cell, cell_model) in enumerate(zip(battery,battery_model)):
                cell(t,current_now,temperature_now)                        
                voltage_measured=cell.voltage[-1] + np.random.normal(loc=0.0, scale=dVoltage, size=None)
                battery_stats['I']['mean']=current_measured
                #battery_stats['I']['mean']=current_now                
                battery_stats['I']['var']=(1e-3*current_measured)**2
                cell_stats['T']['mean']=temperature_measured                    
                # Run the UKF steps
                cell_model.predict()
                cell_model.update(voltage_measured)
                soc_true[j].append(cell.SoC[-1])
                soc_est[j].append(100*cell_model.get_estimate())
                P_values[j].append(cell_model.P)
                I_values[j].append(current_measured)
    previous_t=time    
    if loop_count>max_loop:
        break


for i, cell in enumerate(battery):
    #print(cell)
    plt.figure(figsize=(12, 4))
    plt.subplot(1, 3, 1)
    plt.plot(soc_true[i], label='True SoC')
    plt.plot(soc_est[i], label='Estimated SoC')
    plt.legend()
    plt.title('SoC Estimation')
    
    plt.subplot(1, 3, 2)
    plt.plot(P_values[i])
    plt.title('Covariance P Over Time')
    plt.xlabel('Time Step')
