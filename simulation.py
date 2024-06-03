'''
Add more import commands below as required. 
This is all that is required for the given code.
'''
import random
import math
import numpy as np
import time
import matplotlib.pyplot as plt
from numba import jit ,njit
from mpl_toolkits.mplot3d import Axes3D
from numba.typed import List

    

@njit
def elasticCS(energy: float) -> float:
    '''
    Parameters
    ----------
    energy : float
        electron energy [eV].

    Returns
    -------
    float
        elastic cross section [cm^2]. (Prob of elastic collision)

    Piece-wise cross section approximation.
    '''
    if energy <= 2.0:
        return 7.0e-16 # [cm^2]
    if energy <= 10.0:
        return energy * (-0.27231e-16) + 6.455375e-16
    if energy <= 20.0:
        return energy * (-0.20977e-16) + 6.8801e-16
    if energy <= 100.0:
        return energy * (-0.02981e-16) + 3.2809e-16
    return 3.0e-17 # [cm^2]


@njit
def inelasticCS(energy: float) -> float:
    '''
    Parameters
    ----------
    energy : float
        electron energy [eV].

    Returns
    -------
    float
        inelastic cross section [cm^2]. (Prob of elastic collision)

    Approximate the inelastic cross section as linear from
    the threshold at 24.6 eV to 100 eV, constant after that.
    Cross sections return in units of cm^2.
    '''
    slope = 0.047745e-17
    intercept = -1.174527e-17
    threshold = 24.6 # [eV]
    if energy < threshold:
        return 0.0
    if energy > 100.0:
        return 3.6e-17 # constant max [cm^2]
    return energy * slope + intercept # [cm^2]



@njit
def cross_section(energy: float) -> tuple[bool, float]:

    '''
    Parameters
    ----------
    energy : float
        electron energy [eV].

    Returns
    -------
    inelFlag : bool
        signal inelastic scattering.
    elastic  : float
        elastic scattering cross section [m^2] or zero.

    '''
    inelFlag = False
    elastic = elasticCS(energy) * 1.0e-4 # convert to [m^2]
    inel = inelasticCS(energy) * 1.0e-4 # convert to [m^2]
    if inel <= 0.0: # only elastic
        return inelFlag, elastic
    ratio = inel / (inel + elastic)
    if random.random() < ratio:
        inelFlag = True
        return inelFlag, 0.0
    return inelFlag, elastic



@njit
def run(charges, voltage ):
    interactionPos = []
    transTimes = []
    kmax = 2.0e-12 
    av=6.02214076e23
    tau = 1.0 / ( 0.1664 * 250 * av * kmax) 
    count=0
    while count<len(charges):
        stopFlag = False
        t = 0
        v = np.array([0.0, 0.0, 0.0])
        pos = charges[count]
        print("count= ", count)
        print("num of inelastic= ",len(charges))

        while not stopFlag:
            time_update = -tau * math.log(random.random())                       
            E = electric_field(pos,voltage)               
            v = accelerate(E,time_update,v)
            energy=calcEnergy(v)
            inelFlag , kv = cross_section(energy)    
            kv=kv*np.linalg.norm(v)                      
            #inelastic scattering
            if inelFlag == True:
                kv = 0
                v=setSpeed(v)               
                charges.append(pos)
                interactionPos.append(pos)                
            # elastic scattering
            elif random.random() <= (kv / kmax):
                v = randDirection(v,pos)
                if random.randint(1, 1000) == 1:
                    interactionPos.append(pos)
                
            if math.hypot(pos[0], pos[1]) <  5.0e-4:    
                stopFlag = True
            pos = pos + v*time_update     
            t += time_update       
        # add transport time 
        transTimes.append(t)
        count=count+1
    return interactionPos, transTimes


@njit
def setSpeed(v: np.ndarray) -> np.ndarray:
    m = 9.10938356e-31
    e = 1.602176634e-19  # electron charge
    speed = (0.025* 2  * e / m) ** 0.5
    ratio = speed / np.linalg.norm(v)
    v = ratio * v
    return v

@njit
def calcEnergy(v: np.ndarray) -> np.ndarray:
    m = 9.10938356e-31
    e = 1.602176634e-19
    en=0.5* m * (np.linalg.norm(v))**2 / e           
    return en


@njit
def electric_field(pos , voltage):
    pos2dd =(np.array([pos[0], pos[1]]))
    pos2d = pos2dd/ math.hypot(pos[0], pos[1])
    E = voltage / (math.hypot(pos[0], pos[1]) * math.log(   1.0e-2 /8.0e-5  )) * pos2d          
    return np.array([E[0],E[1], 0])



@njit
def accelerate(electric_field_strength,time_update,v):
    e = -1.60217662e-19 # Coulombs
    m = 9.10938356e-31 # kilograms
    acceleration = electric_field_strength * e / m
    v=acceleration*time_update + v
    return v


@njit
def randDirection(v: float, pos: np.ndarray) -> float:
    found = False
    distance = np.linalg.norm(pos[:2])
    while not found:
        RandomVel = np.array([random.uniform(-1,1),random.uniform(-1,1),random.uniform(-1,1)])
        subtraction = pos - RandomVel
        addition = pos + RandomVel
        if np.linalg.norm(addition[:2]) < np.linalg.norm(subtraction[:2]):
            found = True
    RandomVel = (np.linalg.norm(v) / np.linalg.norm(RandomVel)) * RandomVel
    return RandomVel


positionStart=np.array([3.0e-3 , 0 ,0 ])
lis = List()
lis.append(positionStart)
start = time.process_time()
res, times = run(lis, 1000) # charges, voltages
stop = time.process_time()
print('Process time: ', (stop-start))


for r in res:
    x = r[0]
    z = r[2]
    plt.plot(x, z, 'ro')
plt.title("2D Plot of Transport Locations (x/z)")
plt.xlabel("X Axis")
plt.ylabel("Z Axis")
plt.show()


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=10, azim=90)

for r in res:
    x = r[0]
    y = r[1]
    z = r[2]
    ax.scatter(x,y,z,c='r',marker='o')

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
ax.set_title('Plot of Transport Locations in 3D')

plt.show()




