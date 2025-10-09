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


def fuseCheck(battery,battery_model,current,time):
    fuse=True
    for cell, cell_model in zip(battery,battery_model):
        fuse = fuse and cell.checkFuse(time*current) and cell_model.checkFuse(time*current) 
    return fuse

def calibrate(battery,battery_model,current,Temperature):
    if current>0:
        mode="discharge"
    elif current <0:
        mode="charge"
    else:
        print("Error, calibrating with zero current!!")
        return
    for cell, cell_model in zip(battery,battery_model):
        cell(0,0,Temperature)                       
        voltage_measured=cell.voltage[-1] + np.random.normal(loc=0.0, scale=dVoltage/10, size=None)  #10 since we expect atleast 100x averaging due to 100Hz measurements
        cell_model.start_calibration_cycle(voltage_measured,Temperature,mode)
    coulombCounter_discharge=0
    battery_calibrate=True
    deltaT=1  #time step is 1 sec to reduce simulation time
    t=0
    while (battery_calibrate):
        current_measured=current*(1+np.random.normal(loc=0.0, scale=dCurrent, size=None))
        coulombCounter_discharge += current_measured * deltaT
        t+=deltaT
        battery_calibrate=False
        for j, (cell, cell_model) in enumerate(zip(battery,battery_model)):
            cell(deltaT,current,Temperature)                       
            voltage_measured=cell.voltage[-1] + np.random.normal(loc=0.0, scale=dVoltage, size=None)
            temperature_measured=Temperature + np.random.normal(loc=0.0, scale=dTemperature, size=None)
            if cell_model.calibration_active:
                cell_model.record_calibration_point(current_measured, voltage_measured,temperature_measured, deltaT)
            battery_calibrate=battery_calibrate or cell_model.calibration_active
    return t



#initial_soc, capacity, eta, Id, Q, R

battery_stats = {
    'dt':   {'mean': 1, 'var': 0},
    'I':    {'mean': 0, 'var': 0}
# period of IV measurements, jitter in IV mesurements
#I and (0.0005*I)**2
}
    
cell_stats = {
    'Qc':   {'mean': 65*3600, 'var': 9*0.1*0.1*65*65*3600*3600},
    'eta':  {'mean': 0.98, 'var': 0.01*0.01},
    'Id':   {'mean': 1e-3, 'var': 0.25e-6},
    'T' :   {'mean': 25, 'var': 0.01},
    'R0' :   {'mean': 0.002 , 'var': 0.25e-6},    
}
R0=cell_stats['R0']['mean'] #ohms
RZlist=[[1e-3,10,0],[1.5e-3,1e4,0],[8,5e-7,1]]
ocvLookup="Sample_OCV_DoD.csv"
ocvLookup="lowCornerCalib_ocvTable.csv"
calibrate=False

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
SoC_LastQcUpdate=[]
for i, corner in enumerate(corners):
    battery.append(Cell("Li_Ion", R0, RZlist, corner=corner))
    battery[i].ocvLookup=ocvLookup
    battery_model.append(SoCUKFEstimator(battery_stats, cell_stats, RZlist, ocvLookup))
    soc_true.append([])
    soc_est.append([])
    P_values.append([])
    SoC_LastQcUpdate.append(battery_model[-1].x_est)
           

if calibrate:
    Temperature=25
    #2C discharge
    tsimulated=calibrate(battery,battery_model,130,Temperature)
    #1C charge
    tsimulated+=calibrate(battery,battery_model,-65,Temperature)
else:
    tsimulated=0
    
    #battery_model[0].socEstimator.export_updated_grid_to_csv("lowCornerCalib_ocvTable.csv")

input_csv = 'test_sequence.csv'
#input_csv = 'la92shortdds.csv'
data = pd.read_csv(input_csv)
offset=0
time_measured = data['Time'].values[offset:]
current = data['Current'].values[offset:]
temperature= data['Temperature'].values[offset:]

max_loop=30000
loop_count=0
charge=0
previous_t=tsimulated
t_model=0
previous_fuse=True
dSoC_QcUpdate=0.1
t_measured_prev=0
coulombCount=0
absCoulombCount=0
started_QcRLS=False
toffset=tsimulated
for i, time in enumerate(time_measured):
    if i!=0:
        for nclk in range(int((time+toffset-previous_t)/battery_stats['dt']['mean'])):
            loop_count=loop_count+1
            if loop_count>max_loop:
                break            
            t=toffset+time_measured[i-1] + nclk*battery_stats['dt']['mean']+np.random.normal(loc=0.0, scale=np.sqrt(battery_stats['dt']['var']), size=None)
            t_measured=toffset+time_measured[i-1] + nclk*battery_stats['dt']['mean']          
            current_now=current[i-1]+(current[i]-current[i-1])*nclk*battery_stats['dt']['mean']/(time-previous_t)
            current_measured=current_now*(1+np.random.normal(loc=0.0, scale=dCurrent, size=None))
            fuse_blown=not(fuseCheck(battery,battery_model,current_now,battery_stats['dt']['mean']))
            if fuse_blown:
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
                    t_measured_prev=t_measured-battery_stats['dt']['mean'] 
                previous_fuse=True
                
            temperature_now=temperature[i-1]+(temperature[i]-temperature[i-1])*nclk*battery_stats['dt']['mean']/(time-previous_t)
            temperature_measured=temperature_now+np.random.normal(loc=0.0, scale=dTemperature, size=None)         
            coulombCount+=(t_measured-t_measured_prev)*current_measured
            absCoulombCount+=abs((t_measured-t_measured_prev)*current_measured)
            t_measured_prev=t_measured
            I_values.append(current_measured)
            #if current_measured<0:
                #print("debug now")
            for j, (cell, cell_model) in enumerate(zip(battery,battery_model)):
                cell(battery_stats['dt']['mean'],current_now,temperature_now)                       
                voltage_measured=cell.voltage[-1] + np.random.normal(loc=0.0, scale=dVoltage, size=None)
                battery_stats['I']['mean']=current_measured
                #battery_stats['I']['mean']=current_now                
                battery_stats['I']['var']=(1e-3*current_measured)**2
                cell_stats['T']['mean']=temperature_measured                    
                # Run the UKF steps
                cell_model.predict()
                cell_model.update(voltage_measured,coulombCount)
                soc_true[j].append(cell.SoC[-1])
                soc_est[j].append(100*cell_model.get_estimate())
                P_values[j].append(cell_model.P)

                if not(cell_model.rls_active):
                    if len(P_values[j])>10:
                        if max(P_values[j][-10:])<1e-3:
                            if  all(i>0 for i in I_values[-10:]):
                                cell_model.start_cycle("discharge",coulombCount)
                            elif all(i<0 for i in I_values[-10:]):
                                cell_model.start_cycle("charge",coulombCount)

    previous_t=time    
    if loop_count>max_loop:
        break
column_names=['x_pred','x_est','slope','K','V_observed','v_pred']
df = pd.DataFrame(battery_model[0].Innovation, columns=column_names)
df.to_csv('Innovation_output.csv', index=False)

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
