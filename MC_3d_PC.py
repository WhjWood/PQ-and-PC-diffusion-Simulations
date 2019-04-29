import random,math
import numpy as np
import time
import matplotlib.pyplot as plt


T0 = time.time()

## wwood 25/10/16
# simulating the diffusion of plastocyanin in 3D
# uses brownian motion and monte carlo for a single particle
### The -350 centering has been removed

## wwood 25/10/17
# changed to smulate end membranes
# upper grana and lower stromal lamellae removed and
# upper sl replaced with end membrane
# no junctional slit

## Model constants ( in nm )
PSII_radius = 6.7
PSI_radius = 5
grana_radius = 170 
stroma_radius = 250
boundary = 700 #total width
grana_edge = 2*grana_radius*np.pi

# Diffusion constants
# constants
dx = 1
dt = 1 #/(0.9*10**8)
D = 0.9
RateConst=D*dt/(dx**2) # must be less than one half



def dist(P1,P2):
    if len(P1) == 2 or len(P2) == 2:
        return np.sqrt((P1[0]-P2[0])**2 + (P1[1]-P2[1])**2)
    else:
        return np.sqrt((P1[0]-P2[0])**2 + (P1[1]-P2[1])**2 + (P1[2]-P2[2])**2)

# download PSII/PSII/b6f particles from previous MCMC
# change to integer (nm grid at a later stage

PSII_COORDINATES = open('grana_particles_22_07_16_experiment_1','r').readlines()# this is double layer overlapped
UPPER_PSII_POPULATION = []
for i in PSII_COORDINATES:
    i.strip("\n")
    F = i.split(" ")
    if dist([int(float(F[0]))-350,int(float(F[1]))-350],[0,0]) < grana_radius-PSII_radius:
        UPPER_PSII_POPULATION.append([int(float(F[0]))-350,int(float(F[1]))-350]) ## - 350 is needed to normalise


PSII_COORDINATES = open('grana_particles_22_07_16_experiment_2','r').readlines()# this is double layer overlapped
LOWER_PSII_POPULATION = []
for i in PSII_COORDINATES:
    i.strip("\n")
    F = i.split(" ")
    if dist([int(float(F[0]))-350,int(float(F[1]))-350],[0,0]) < grana_radius-PSII_radius:
        LOWER_PSII_POPULATION.append([int(float(F[0]))-350,int(float(F[1]))-350]) ## - 350 is needed to normalise



PSII_POPULATION = [UPPER_PSII_POPULATION,LOWER_PSII_POPULATION]


# x2 PSI coordinates representing 2 layers
PSI_COORDINATES = open('end_particles_25_07_17_experiment_','r').readlines()
UPPER_PSI_POPULATION = []
for i in PSI_COORDINATES:
    i.strip("\n")
    F = i.split(" ")
    if (dist([int(float(F[0]))-350,int(float(F[1]))-350],[0,0]) > grana_radius + PSI_radius) & (dist([int(float(F[0]))-350,int(float(F[1]))-350],[0,0]) < grana_radius + PSI_radius+100):
        UPPER_PSI_POPULATION.append([int(float(F[0]))-350,int(float(F[1]))-350,10])

PSI_COORDINATES = open('end_particles_25_07_17_experiment_','r').readlines()
LOWER_PSI_POPULATION = []
for i in PSI_COORDINATES:
    i.strip("\n")
    F = i.split(" ")
    if (dist([int(float(F[0]))-350,int(float(F[1]))-350],[0,0]) > grana_radius + PSI_radius) & (dist([int(float(F[0]))-350,int(float(F[1]))-350],[0,0]) < grana_radius + PSI_radius+100):
        LOWER_PSI_POPULATION.append([-(int(float(F[0]))-350),-(int(float(F[1]))-350),0])



PSI_POPULATION = [UPPER_PSI_POPULATION,LOWER_PSI_POPULATION]
## create boundaries no. 2 commented out for now


# outer boundary
BOUNDARIES2 = [(int(stroma_radius*np.cos(T*(2*np.pi/1000.0))-np.sin(T*(2*np.pi/1000.0))),int(stroma_radius*np.sin(T*(2*np.pi/1000.0))+np.cos(T*(2*np.pi/1000.0)))) for T in range(1000)]



# create b6f population
"""
Nb6f = int(len(PSII_POPULATION)/3.0) # this is where the fraction of particles that are b6f is set
B6F_POPULATION = []
while len(B6F_POPULATION) <= Nb6f:
    NEW = random.choice(PSII_POPULATION)
    NEW2 = [int(NEW[0]),int(NEW[1])]
    if NEW2 not in B6F_POPULATION:
        B6F_POPULATION.append(NEW)
        PSII_POPULATION.pop(PSII_POPULATION.index(NEW))
"""



def step(P):
    S = random.choice([[-1,0,0],[0,-1,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
    
    if P[2]+S[2] <= 0 or P[2]+S[2] >= 10 or (P[0]+S[0],P[1]+S[1]) in BOUNDARIES:

        return P

    
    if P[2]+S[2] >=5:
        for psii in PSII_POPULATION[0]:
            
            if dist([P[0]+S[0],P[1]+S[1]],psii)<PSII_radius:

                return P
   
    elif P[2]+S[2] <=5:
        for psii in PSII_POPULATION[1]:
            if dist([P[0]+S[0],P[1]+S[1]],psii)<PSII_radius:

                 return P

    return [P[0]+S[0],P[1]+S[1],P[2]+S[2]]
   


def display_coordinates(pop):
    for p in pop:
        print(p[0],p[1])
      
#main
DATA = {}

for B in [5000]: #range(1000,10000,1000): # no junctional slit
    print("Current B:"+str(B)+" Current time" +str(time.time()-T0))
    Swidth = B
    #3 equal size slits starting from 2pi/3 radians apart
    BOUNDARIES1a = [(int(grana_radius*np.cos(T*(2*np.pi/10000.0))-np.sin(T*(2*np.pi/10000.0))),int(grana_radius*np.sin(T*(2*np.pi/10000.0))+np.cos(T*(2*np.pi/10000.0)))) for T in range(int(Swidth/3.0))]
    BOUNDARIES1b = [(int(grana_radius*np.cos(T*(2*np.pi/10000.0 + 2*np.pi/3))-np.sin(T*(2*np.pi/10000.0 + 2*np.pi/3))),int(grana_radius*np.sin(T*(2*np.pi/10000.0 + 2*np.pi/3))+np.cos(T*(2*np.pi/10000.0 + 2*np.pi/3)))) for T in range(int(Swidth/3.0))]
    BOUNDARIES1c = [(int(grana_radius*np.cos(T*(2*np.pi/10000.0 + 4*np.pi/3))-np.sin(T*(2*np.pi/10000.0 + 4*np.pi/3))),int(grana_radius*np.sin(T*(2*np.pi/10000.0 + 4*np.pi/3))+np.cos(T*(2*np.pi/10000.0 + 4*np.pi/3)))) for T in range(int(Swidth/3.0))]
    #BOUNDARIES1 = [(int(grana_radius*np.cos(T*(2*np.pi/1000.0))-np.sin(T*(2*np.pi/1000.0))),int(grana_radius*np.sin(T*(2*np.pi/1000.0))+np.cos(T*(2*np.pi/1000.0)))) for T in range(1000)]
    
    #Bshow = set(BOUNDARIES1a[:]) | set(BOUNDARIES1b[:]) | set(BOUNDARIES1c[:])
    
    BOUNDARIES = set(BOUNDARIES1a[:]) | set(BOUNDARIES1b[:]) | set(BOUNDARIES1c[:]) | set(BOUNDARIES2[:])
    
    """
    #the following is for displaying the slits
    Bx = []
    By = []
    for b in Bshow:
        Bx.append(b[0])
        By.append(b[1])
    plt.scatter(np.array(Bx),np.array(By))
    plt.show()
    """
    T = []
    #datafile=open("results_light_grana_01_09_17_slit_"+str(Swidth)+"graph1", 'w')
    
    
    for i in range(1):
        print(i)
        L = random.choice([1,0]) # the layer used
        
        P = random.choice(PSII_POPULATION[L])
        #P = random.choice(PSII_POPULATION)
        
        PSII_POPULATION[L].pop(PSII_POPULATION[L].index(P)) # remove the current particle from the population
        #PSII_POPULATION.pop(PSII_POPULATION.index(P)) # remove the current particle from the population
        B6F = P
        
        
        if L == 0:
            P.append(9) # the z component
        elif L == 1:
            P.append(1)
        
        
        
        psi_collision = False
        t = 0
    
    
    
        while psi_collision == False and t <= 10**8: # the 10**8 is meant to signify trapped particles
            P = step(P)
            
            t += dt
        #X.append(P[0])
        #Y.append(P[1])
            if dist(P[0:2],[0,0]) >= grana_radius:
                if P[2] > 5:
                    for psi in PSI_POPULATION[0]:
                    
                        if dist(P,psi)<PSI_radius:
                            psi_collision = True
                            break
                elif P[2] < 5:
                        
                    for psi in PSI_POPULATION[1]:
                        if dist(P,psi)<PSI_radius:
                            psi_collision = True
                            break
                    
            
        T.append(t)
        #datafile.write(str(B6F[0])+" "+str(B6F[1])+" "+str(P[0])+" "+str(P[1])+" "+str(t)+"\n")
        PSII_POPULATION[L].append(B6F)
        #PSII_POPULATION.append(B6F[0:2])
        #print(i)
    DATA[B] = T
    
    
    #datafile.close()

print(time.time()-T0)
