import numpy as np
import random


endfile = open("end.txt",'r').readlines()
END = []
for e in endfile:
    END.append(float(e))

darkfile = open("Dark.txt",'r').readlines()
dark = []
for d in darkfile:
    dark.append(float(d))

lightfile = open("Light.txt", 'r').readlines()
light = []
for l in lightfile:
    light.append(float(l))
    
current = light
Pend = 2.0/10.6
# Dark = 2.0/17.16
#Light = 2.0/10.6

for i in range(1000):
    if random.random() < Pend:
        X = random.choice(END)
    else:
        X = random.choice(current)
    
    print(X)