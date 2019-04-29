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


def Save_system(SYSTEM, filename):
            # saves population in CD 
    SYS_FILE = open(filename, 'w')
    for i in SYSTEM:
        SYS_FILE.writelines(str(i[0])+" "+str(i[1])+"\n")

for grana_radius in [190]:    
    for Nos in [0.3]: #[0.25, 0.5, 0.75]: #the number of b6fs in the G and SL
    
        ## Model constants ( in nm )
        PSII_radius = 6.7
        PSI_radius = 5
        LHCII_radius = 3
        #grana_radius = 190 
        stroma_radius = grana_radius+50
        boundary = 700 #total width
        grana_edge = 2*grana_radius*np.pi
        
        # Diffusion constants
        # constants
        dx = 1
        dt = 1 #/(0.9*10**8)
        D = 0.9
        RateConst=D*dt/(dx**2) # must be less than one half
        
        def Plot_coords(POP,r):
            X = []
            Y = []
            for xy in POP:
                X.append(xy[0])
                Y.append(xy[1])
            plt.scatter(np.array(X),np.array(Y),s=r)
            plt.show()
        
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
        
        
        
        PSII_POPULATION = UPPER_PSII_POPULATION #,LOWER_PSII_POPULATION]
        """
        LHCII_POPULATION = []
        # the following creates the LHCII population to be commented out once saved
        # code adapted from hard discs 2
        def COLLISION(P1,P2,R):
            return dist(P1,P2) < R
        
        def SYST_COLL(P,SYSTEM,R):
            
            for P1 in SYSTEM:
                if COLLISION(P,P1,R) == True and P != P1:
                    return True
            return False
         
        while len(LHCII_POPULATION) <= 4*len(PSII_POPULATION): #N LHCII
            P = [int(random.random()*2*grana_radius) - grana_radius,int(random.random()*2*grana_radius)-grana_radius]
            while dist(P,[0,0]) > grana_radius-LHCII_radius:
                P = [int(random.random()*2*grana_radius) - grana_radius,int(random.random()*2*grana_radius)-grana_radius]
            del_x = 0
            del_y = 0
            if len(LHCII_POPULATION)> 0:
                C = SYST_COLL(P,PSII_POPULATION,PSII_radius+LHCII_radius)+SYST_COLL(P,LHCII_POPULATION,LHCII_radius+LHCII_radius)
            else:    
                C = SYST_COLL(P,PSII_POPULATION,PSII_radius+LHCII_radius)
            while C == True:      
                del_x = random.choice([-1,1])*dx
                del_y = random.choice([-1,1])*dx
                P1 = [P[0]+del_x,P[1]+del_y]
                if dist(P1,[0,0]) < grana_radius-LHCII_radius:
                    P = P1
                if len(LHCII_POPULATION)> 0:
                    C = SYST_COLL(P,PSII_POPULATION,PSII_radius+LHCII_radius)+SYST_COLL(P,LHCII_POPULATION,LHCII_radius+LHCII_radius)
                else:    
                    C = SYST_COLL(P,PSII_POPULATION,PSII_radius+LHCII_radius)
            LHCII_POPULATION.append(P)
            #print(P)
        
        def Save_system(SYSTEM, filename):
            # saves population in CD 
            SYS_FILE = open(filename, 'w')
            for i in SYSTEM:
                SYS_FILE.writelines(str(i[0])+" "+str(i[1])+"\n")
        
        
        #Save_system(LHCII_POPULATION,"lhcii_particles_06_09_17")
        #end of LHCII creation
        
        
        #load the LHCII file
        """
        LHCII_COORDINATES = open("lhcii_particles_06_09_17_"+str(Nos)+" "+str(grana_radius),'r').readlines()
        LHCII_POPULATION = []
        for i in LHCII_COORDINATES:
            i.strip("\n")
            F = i.split(" ")
            if dist([int(float(F[0])),int(float(F[1]))],[0,0]) < grana_radius-LHCII_radius:
                LHCII_POPULATION.append([int(float(F[0])),int(float(F[1]))])
        
        # x2 PSI coordinates representing 2 layers
        PSI_COORDINATES = open('end_particles_25_07_17_experiment_','r').readlines()
        UPPER_PSI_POPULATION = []
        for i in PSI_COORDINATES:
            i.strip("\n")
            F = i.split(" ")
            if (dist([int(float(F[0]))-350,int(float(F[1]))-350],[0,0]) > grana_radius + PSI_radius) & (dist([int(float(F[0]))-350,int(float(F[1]))-350],[0,0]) < grana_radius + PSI_radius+50):
                UPPER_PSI_POPULATION.append([int(float(F[0]))-350,int(float(F[1]))-350])
        
        """
        PSI_COORDINATES = open('stroma_particles_22_07_16_experiment_1','r').readlines()
        LOWER_PSI_POPULATION = []
        for i in PSI_COORDINATES:
            i.strip("\n")
            F = i.split(" ")
            if (dist([int(float(F[0]))-350,int(float(F[1]))-350],[0,0]) > grana_radius + PSI_radius) & (dist([int(float(F[0]))-350,int(float(F[1]))-350],[0,0]) < grana_radius + PSI_radius+50):
                LOWER_PSI_POPULATION.append([-(int(float(F[0]))-350),-(int(float(F[1]))-350)])
        """
        
        
        PSI_POPULATION = UPPER_PSI_POPULATION#+LOWER_PSI_POPULATION
        ## create boundaries no. 2 commented out for now
        
        
        # outer boundary
        BOUNDARIES2 = [(int(stroma_radius*np.cos(T*(2*np.pi/1000.0))-np.sin(T*(2*np.pi/1000.0))),int(stroma_radius*np.sin(T*(2*np.pi/1000.0))+np.cos(T*(2*np.pi/1000.0)))) for T in range(1000)]
        
        
       
        # create b6f population
        print(len(PSII_POPULATION))
        Nb6f = 0.6*len(PSI_POPULATION)
        #Nb6f = int(len(PSII_POPULATION)*0.4) # this is where the fraction of particles that are b6f is set
        
        B6F_POPULATION = []
        while len(B6F_POPULATION) <= Nb6f*Nos:
            #print(len(B6F_POPULATION))
            NEW = random.choice(PSII_POPULATION)
            NEW2 = [int(NEW[0]),int(NEW[1])]
            if NEW2 not in B6F_POPULATION:
                B6F_POPULATION.append(NEW)
                PSII_POPULATION.pop(PSII_POPULATION.index(NEW))
        
         #int(len(PSI_POPULATION)/4) # this is where the fraction of particles that are b6f is set
        B6F_POPULATION_stroma = []
        while len(B6F_POPULATION_stroma) <= Nb6f*(1-Nos):
            NEW = random.choice(PSI_POPULATION)
            NEW2 = [int(NEW[0]),int(NEW[1])]
            if NEW2 not in B6F_POPULATION_stroma:
                B6F_POPULATION_stroma.append(NEW)
                PSI_POPULATION.pop(PSI_POPULATION.index(NEW))
        
        Save_system(B6F_POPULATION, "b6f_grana_particles "+str(Nos)+" "+str(grana_radius))
        Save_system(B6F_POPULATION_stroma, "b6f_stroma_particles "+str(Nos)+" "+str(grana_radius))
        
        """
        B6F_POPULATION_stroma = []
        while len(B6F_POPULATION_stroma) <= Nos*len(B6F_POPULATION):
            print(len(B6F_POPULATION_stroma)) 
            P = [int(random.random()*2*stroma_radius) - stroma_radius,int(random.random()*2*stroma_radius)-stroma_radius]
            while dist(P,[0,0]) < grana_radius+PSII_radius:
                P = [int(random.random()*2*stroma_radius) - stroma_radius,int(random.random()*2*stroma_radius)-stroma_radius]
            del_x = 0
            del_y = 0
            if len(B6F_POPULATION_stroma)> 0:
                C = SYST_COLL(P,PSI_POPULATION,PSII_radius+PSI_radius)+SYST_COLL(P,B6F_POPULATION_stroma,PSII_radius+PSII_radius)
            else:    
                C = SYST_COLL(P,PSI_POPULATION,PSII_radius+PSI_radius)
            while C == True:      
                del_x = random.choice([-1,1])*dx
                del_y = random.choice([-1,1])*dx
                P1 = [P[0]+del_x,P[1]+del_y]
                if (dist(P1,[0,0]) > grana_radius+PSII_radius) and (dist(P1,[0,0]) < stroma_radius - PSII_radius):
                    P = P1
                if len(B6F_POPULATION_stroma)> 0:
                    C = SYST_COLL(P,PSI_POPULATION,PSII_radius+PSI_radius)+SYST_COLL(P,B6F_POPULATION_stroma,PSII_radius+PSII_radius)
                else:    
                    C = SYST_COLL(P,PSI_POPULATION,PSII_radius+PSI_radius)
            B6F_POPULATION_stroma.append(P)
        
        
        B6F_POPULATION = []
        
        #B6F_COORDINATES = open("b6f_grana_particles "+str(Nos)+" "+str(grana_radius),'r').readlines()
        
        for i in B6F_COORDINATES:
            B6F_POPULATION
            i.strip("\n")
            F = i.split(" ")
            if dist([int(float(F[0])),int(float(F[1]))],[0,0]) < grana_radius-PSII_radius:
                B6F_POPULATION.append([int(float(F[0])),int(float(F[1]))])
        
        
        #remove B6Fs from PSII pop
        for B in B6F_POPULATION:
            if B in PSII_POPULATION:
                PSII_POPULATION.pop(PSII_POPULATION.index(B))
                
        B6F_POPULATION_stroma = []
        
        
        B6F_COORDINATES_stroma = open("b6f_stroma_particles "+str(Nos)+" "+str(grana_radius),'r').readlines()
        
        for i in B6F_COORDINATES_stroma:
            i.strip("\n")
            F = i.split(" ")
            #if dist([int(float(F[0])),int(float(F[1]))],[0,0]) < grana_radius-PSII_radius:
            B6F_POPULATION_stroma.append([int(float(F[0])),int(float(F[1]))])
        """
        def step(P):
            S = random.choice([[-1,0],[0,-1],[1,0],[0,1]])
            
            if (P[0]+S[0],P[1]+S[1]) in BOUNDARIES:
        
                return P
        
            
        
            for psii in PSII_POPULATION:
                if dist([P[0]+S[0],P[1]+S[1]],psii)<PSII_radius:
                    return P
            
            for psi in PSI_POPULATION:
                if dist([P[0]+S[0],P[1]+S[1]],psi)<PSI_radius:
                    return P
        
        
            return [P[0]+S[0],P[1]+S[1]]
        
        
        
        def display_coordinates(pop):
            for p in pop:
                print(p[0],p[1])
            
        #main
        DATA = {}
        G_to_SL = []
        
        print("PSII: ",len(PSII_POPULATION))
        print("PSI: ",len(PSI_POPULATION))
        print("B6F_GRANA: ",len(B6F_POPULATION))
        print("B6F_SL: ",len(B6F_POPULATION_stroma))
        
        for B in range(5000,6000,1000): # no junctional slit
            print("Current B:"+str(B)+" "+str(grana_radius)+" "+str(Nos)+" Current time" +str(time.time()-T0))
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
            LOC = []
            NOS_STRING = str(10*Nos)[0]
            #datafile=open("results_"+str(grana_radius)+"NosG"+NOS_STRING+"_grana_PQ_25_09_17_slit_"+str(Swidth)+"graph1", 'w')
            
            
            for i in range(1):
                #print(i)
        
                
                P = random.choice(PSII_POPULATION)
        
                
                PSII_POPULATION.pop(PSII_POPULATION.index(P)) # remove the current particle from the population
                
                Pstart = P
                
                
                
                
                
                
                b6f_collision = False
                t = 0
            
            
            
                while b6f_collision == False and t <= 35*10**4: #10 s trapped particles
                    P = step(P)
                     #print(P[0],P[1])
                    t += dt
                #X.append(P[0])
                #Y.append(P[1])
                    if dist(P[0:2],[0,0]) <= grana_radius:
                        for b6f in B6F_POPULATION:
                            
                            if dist(P,b6f)<PSII_radius:
                                b6f_collision = True
                                break
                        
                    else:            
                        for b6f in B6F_POPULATION_stroma:
                            
                            if dist(P,b6f)<PSII_radius:
                                b6f_collision = True
                                break
                            
                    
                T.append(t)
                #datafile.write(str(Pstart[0])+" "+str(Pstart[1])+" "+str(P[0])+" "+str(P[1])+" "+str(t/4)+"\n")
                # the t/4 in file is time in ms
                PSII_POPULATION.append(Pstart)
                LOC.append([P[0],P[1]])
                #PSII_POPULATION.append(B6F[0:2])
                #print(i)
            PrSL = 0
            for l in LOC:
                if dist(l,[0,0])>grana_radius:
                    PrSL += 1
            print([100-B/100,PrSL])
            DATA[B] = T
            
            
            #datafile.close()


print(time.time()-T0)
Plot_coords(PSII_POPULATION,PSII_radius)
Plot_coords(LHCII_POPULATION,LHCII_radius)
Plot_coords(PSI_POPULATION,PSI_radius)
Plot_coords(B6F_POPULATION,PSII_radius)
Plot_coords(B6F_POPULATION_stroma,PSII_radius)

