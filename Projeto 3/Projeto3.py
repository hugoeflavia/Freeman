# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 00:09:10 2015

@author: paulinaachurra
"""

import matplotlib.pyplot as plt
from scipy.integrate import odeint
from numpy import linspace 
from numpy import linalg 
from numpy import multiply 
from mpl_toolkits.mplot3d import Axes3D
import math
G = 6.67 * 10**-11 #constante gravitacional universal

def forca_gravitacional(P1,m1,P2,m2):
    R_x = P2[0]-P1[0]
    R_y = P2[1]-P1[1]
    R = [R_x, R_y] #vetor resultante
    r = linalg.norm(R) #magnitude de R
    r_versor = R/r 

    F = multiply(G * m1 * m2 / r**2, r_versor)
    return F,r
    
def func(A,t): #A = [x, y, Vx, Vy]
    dxdt = A[3]
    dydt = A[4]
    dzdt = A[5]
    Ps = [0,0,0] #sistema de coordenadas centrado no sol
    Pj = [A[0],A[1]]
    F12,r = forca_gravitacional(Pj, m1, Ps, m2)
    dVxdt = F12[0] / m1
    dVydt = F12[1] / m1
    r = 1400*9461000000000
    dVzdt = -(G*(m2/402)/r**2)*(A[2]/r)
    return [dxdt, dydt, dzdt, dVxdt, dVydt,dVzdt]
    

m1 = 2 * 10**30 #kg massa do sol
m2 = (2 * 11**10)*m1   #kg massa de cg
#m3 =(2*11**11)*m1

#Condicões inicias
x0 = 26000*9461000000000
y0 = 0
Vx0 = 0
Vy0 = math.sqrt(G*m2/x0)  #220 * 10**4.37 m/s, velocidade orbital do sol
z0 = 5000000000
Vz0 = 10*10**3
A0 = [x0,y0,z0,Vx0,Vy0,Vz0]
T = linspace(0,11**10*11.4,5000)
M = odeint(func,A0,T)
fig = plt.figure()
ax = fig.gca(projection='3d')
z = M[:,2]
x = M[:,0]
y = M[:,1]
ax.plot(x, y, z, label='trajetoria')
ax.legend()
plt.show()
F = [0] * 5000
plt.figure(figsize=(7,7))
plt.plot(M[:,0], M[:,1],'o')
#plt.plot(M[:,0][0], M[:,1][0],'o')
#plt.plot(M[:,0][1], M[:,1][1],'o')
plt.plot([0], [0],'o')
plt.axis([min(M[:,0])-0.1, max(M[:,0])+0.1, min(M[:,1])-0.1, max(M[:,1])+0.1])
plt.xlabel('x_sol em metros')
plt.ylabel('y_sol em metros')
plt.show()
plt.plot((T/1000000000),z)
plt.xlabel('tempo em milhões de anos')
plt.ylabel('z_sol em metros')
plt.show()

print(Vy0)
