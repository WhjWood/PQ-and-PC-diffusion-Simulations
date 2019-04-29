import numpy as np
import os



FULL_DATA =  []
for filename in os.listdir(os.getcwd()):
    if filename[8:12] != "dark":
        #print(filename)
        
        DATA = open(filename, 'r').readlines()
        T_data = []
        for line in DATA:
            f = line.split(" ")
            
            if len(f) == 5:
                T_data.append(float(f[4]))
        T = np.array(T_data)
        FULL_DATA.append([int(filename[-9:-6].strip("_")),T.mean()])
        print(filename, T.min(),T.mean(),T.max())
            
        

