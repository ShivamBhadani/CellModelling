# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 19:54:10 2025

@author: iitm9
"""


#initial_soc, capacity, eta, Id, Q, R
variable_stats = {
    'dt':   {'mean': 2.0, 'var': 0.1},
    'Qc':   {'mean': 3.0, 'var': 0.2},
    'eta':  {'mean': 1.5, 'var': 0.05},
    'I':    {'mean': 0, 'var': 0.3},
    'Id':   {'mean': 1.2, 'var': 0.1}
}


from batteryCell import Cell







corners=["nom","low","high","mc",'mc',"mc",'mc',"mc"]
#corners=['nom']
battery=[]
battery_model=[]
for corner in corners:
    battery.append(Cell("Li_Ion", corner=corner))
                
input_csv = 'test_sequence.csv'
#input_csv = 'la92shortdds.csv'
data = pd.read_csv(input_csv)
offset=0
time = data['Time'].values[offset:]
current = data['Current'].values[offset:]
temperature= data['Temperature'].values[offset:]
charge=0
for i, t in enumerate(time):
    fuse=True
    if i!=0:
        for cell in battery:
            fuse = fuse and cell.checkFuse(t,current[i])

        if fuse:
            for cell in battery:
                cell(t,current[i],temperature[i])
                
        else:
            print("protection fuse blown")
        
    


for cell in battery:
    print(cell)
