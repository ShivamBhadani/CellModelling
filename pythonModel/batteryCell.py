# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 20:35:09 2025

@author: Rajesh N Gupta
"""
import numpy as np
from cellESR import cellESR, LookupTableEstimator

class Cell:
    def __init__(self, chemID, R0, RZlist, corner="nom", seed=None):
        self.cellName='generic'
        self.SoC=[100.0]   # in percentage
        self.SoH=[100.0]   # percentage of max capacity
        self.MaxCapacity=65000   #units of milli-amp hour at time zero
        self.capacityStdev=300   # sigma is 3%
        self.charging=True #current charging status
        self.discharge=True #current discharging status
        self.cellESR=cellESR(R0,RZlist)
        self.charge=0 # denotes readiness to charge
        ocvLookup="Sample_OCV_DoD.csv"  # Path to the lookup table file
        self.ocvEstimator=LookupTableEstimator(ocvLookup)
        self.Temperature=[25]  # cell temperature in degree centigrade
        self.current=[0]   # cell current positive means discharge negative means charging
        self.esr=[self.cellESR.esrDC] 
        self.ocv=[self.ocvEstimator.estimateOCV(100-self.SoC[-1],self.Temperature[-1])] # initial OCV given temp and DoD
        self.voltage=[self.ocv[-1]+self.current[-1]*self.esr[-1]] # initial voltage
        self.time=[0]          # keep track of time
        self.coulombEffeciency=[0.95]
        self.SoCLmax=0
        self.SoCLmin=100
        self.SoCnoiseThreshold=0.5 # threshold to determine if we are in charging or discharging mode
        self.SoCextremaList=[]
        self.SoCmaxDegradation=0
        self.DoDmaxDegradation=0
        if seed is not None: #seed for random number generator
            np.random.seed(seed)

        if corner == "low":
            self.MaxCapacity *= 0.9
            esrFactor = 1.1           
            self.SoH = [90.0]
            self.SoC = [85.0]
            self.coulombEffeciency=[0.90]
        elif corner == "high":
            self.MaxCapacity *= 1.1
            esrFactor = 0.9
            self.SoH = [100.0]
            self.SoC = [100.0]
            self.coulombEffeciency=[0.98]
        elif corner == "mc":
            self.MaxCapacity = np.random.normal(self.MaxCapacity, self.capacityStdev)
            esrFactor = np.random.normal(1, 0.03)
            self.SoH = [min(100,np.random.normal(95.0, 1))]
            self.SoC = [min(100,np.random.normal(95.0, 1))]
            self.coulombEffeciency=[min(95,np.random.normal(90.0, 1))]
        else:
            #we are at nominal
            esrFactor = 1
            
        R0*=esrFactor
        for index, rc in enumerate(RZlist):
            RZlist[index][0]=rc[0]*esrFactor

    def updateSoH(self):
        dSoH=self.DoDmaxDegradation+self.SoCmaxDegradation
        self.SoH.append(self.SoH[-1]-dSoH)
        
    def checkFuse(self, time, i):
        dt= time - self.time[-1]
        discharge=i*dt
        if discharge>0:
            # We are discharging
            dSoC=100*discharge/(self.coulombEffeciency[-1]*(3.6*self.SoH[-1]*self.MaxCapacity/100))          
        else:
            # We are charging
            dSoC=100*discharge/(3.6*self.SoH[-1]*self.MaxCapacity/100)

        if (self.SoC[-1]-dSoC)>=99 and i<0:
            self.charge=False
            return False
        elif (self.SoC[-1]-dSoC)<=1 and i>0:
            self.discharge=False
            return False
        else:
            self.charge=True
            self.discharge=True
            return True
        
    def SoCdegradation(self):
        print("SoCdegradation called")
        self.DoDmaxDegradation=0
        self.SoCmaxDegradation=0
        if self.charging and self.SoCLmin<80:
            self.DoDmaxDegradation=(1/1000)*(50/(90*90))*((100-self.SoCLmin)**2)
        elif not self.charging and self.SoCLmax>20:
            #SoCMax just got calculated
            self.SoCmaxDegradation=(1/1000)*(50/(90*90))*((self.SoCLmax)**2)
        self.updateSoH()
           
    def ocvSoH_gain(self):
        return 1+((100-self.SoH[-1])/15000)*(self.SoC[-1]-50)
    
    def __call__(self,time,i, Temperature):
        dt= time - self.time[-1]
        discharge=i*dt
        if discharge>0:
            # We are discharging
            dSoC=100*discharge/(self.coulombEffeciency[-1]*(3.6*self.SoH[-1]*self.MaxCapacity/100))            
        else:
            # We are charging
            dSoC=100*discharge/(3.6*self.SoH[-1]*self.MaxCapacity/100)

        if(not self.checkFuse(time,i)):
            print("protection fuse blown, run simulation with zero current , current=", i, " SoC=", self.SoC[-1])
            i=0
            return
        
        self.Temperature.append(Temperature)
        self.time.append(time)
        self.esr.append(self.cellESR.calculateESR(i,dt))
        self.SoC.append(self.SoC[-1]-dSoC)
        self.ocv.append((self.ocvEstimator.estimateOCV(100-self.SoC[-1],self.Temperature[-1]))*self.ocvSoH_gain())
        self.voltage.append(self.ocv[-1]+i*self.esr[-1])
        self.current.append(i)
        if self.charging:
            if (self.SoCLmax-self.SoC[-1])>self.SoCnoiseThreshold:
                self.charging=0
                self.SoCextremaList.append(self.SoCLmax)
                self.SoCLmin=self.SoC[-1]
                self.SoCdegradation()
            else:
                self.SoCLmax=max(self.SoCLmax,self.SoC[-1])
        else:
            if (self.SoC[-1]-self.SoCLmin)>self.SoCnoiseThreshold:
                self.charging=1
                self.SoCextremaList.append(self.SoCLmin)
                self.SoCLmax=self.SoC[-1]
                self.SoCdegradation()
            else:
                self.SoCLmin=min(self.SoCLmin,self.SoC[-1])