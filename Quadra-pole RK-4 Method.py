#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import math 
import pandas as pd
import matplotlib.pyplot as plt


# In[ ]:




#Initial Value
r_0=input("Initial value of the distance of the charged particle from origin: ")
theta_0=input("Initial value of polar angle: ")
phi_0=input("Initial value of the azimuthal angle: ")
p_r_0=input("Initial value of the radial momentum: ")
p_theta_0=input("Initial value of the polar angular momentum: ")
p_phi_0=input("Initial value of the azimuthal angular momentum: ")
lamb=input("Value of lambda: ")
t=0


# In[ ]:


r_0= float(r_0) 
theta_0= float(theta_0) 
phi_0= float(phi_0) 
p_r_0= float(p_r_0) 
p_theta_0= float(p_theta_0) 
p_phi_0= float(p_phi_0)
ini=np.array([float(r_0), float(theta_0), float(phi_0), float(p_r_0), float(p_theta_0), float(p_phi_0)])
ls=math.sin(float(lamb))
lc=math.cos(float(lamb))
kqS=1
h=0.01

"""
r_dot= p_r_0
theta_dot= p_theta_0 / r_0**2
phi_dot= p_phi_0/(r_0**2 * (math.sin(p_theta_0))**2)
p_r_dot= (p_theta_0**2/r_0**3)+(p_phi_0**2 /(r_0**3 * math.sin(p_theta_0)**2)) + (kqS/2)*(1/r_0**4)*(3*ls*math.sin(2*theta_0)*math.sin(phi_0)+3*lc*math.cos(2*theta_0)+lc)
p_theta_dot= (p_phi_0**2 * math.cos(theta_0)/(r_0**2 * math.sin(p_theta_0)**3)) + (kqS/r_0**3)*(lc*math.sin(2*theta_0)-ls*math.cos(2*theta_0)*math.sin(phi_0))

#x'=f(x)"""


# In[ ]:


def Quadrapole(ini, h):
    r_dot= ini[3]
    theta_dot= ini[4] / ini[0]**2
    phi_dot= ini[5]/(ini[0]**2 * (math.sin(ini[4]))**2)
    p_r_dot= (ini[4]**2/ini[0]**3)+(ini[5]**2 /(ini[0]**3 * math.sin(ini[4])**2)) + (kqS/2)*(1/ini[0]**4)*(3*ls*math.sin(2*ini[1])*math.sin(ini[2])+3*lc*math.cos(2*ini[1])+lc)
    p_theta_dot= (ini[5]**2 * math.cos(ini[1])/(ini[0]**2 * math.sin(ini[4])**3)) + (kqS/ini[0]**3)*(lc*math.sin(2*ini[1])-ls*math.cos(2*ini[1])*math.sin(ini[2]))
    p_phi_dot= -(kqS/(2*ini[0]**3))*ls*math.sin(2*ini[1])*math.cos(ini[2])
    
    
    a=np.array([r_dot, theta_dot, phi_dot, p_r_dot, p_theta_dot, p_phi_dot])
    ini1=ini+(h/2)*np.array([r_dot, theta_dot, phi_dot, p_r_dot, p_theta_dot, p_phi_dot])
    
    
    r_dot= ini1[3]
    theta_dot= ini1[4] / ini1[0]**2
    phi_dot= ini1[5]/(ini1[0]**2 * (math.sin(ini1[4]))**2)
    p_r_dot= (ini1[4]**2/ini1[0]**3)+(ini1[5]**2 /(ini1[0]**3 * math.sin(ini1[4])**2)) + (kqS/2)*(1/ini1[0]**4)*(3*ls*math.sin(2*ini1[1])*math.sin(ini1[2])+3*lc*math.cos(2*ini1[1])+lc)
    p_theta_dot= (ini1[5]**2 * math.cos(ini1[1])/(ini1[0]**2 * math.sin(ini1[4])**3)) + (kqS/ini1[0]**3)*(lc*math.sin(2*ini1[1])-ls*math.cos(2*ini1[1])*math.sin(ini1[2]))
    p_phi_dot= -(kqS/(2*ini1[0]**3))*ls*math.sin(2*ini1[1])*math.cos(ini1[2])
    
    
    ini2=ini+(h/2)*np.array([r_dot, theta_dot, phi_dot, p_r_dot, p_theta_dot, p_phi_dot])
    
    
    r_dot= ini2[3]
    theta_dot= ini2[4] / ini2[0]**2
    phi_dot= ini2[5]/(ini2[0]**2 * (math.sin(ini2[4]))**2)
    p_r_dot= (ini2[4]**2/ini2[0]**3)+(ini2[5]**2 /(ini2[0]**3 * math.sin(ini2[4])**2)) + (kqS/2)*(1/ini2[0]**4)*(3*ls*math.sin(2*ini2[1])*math.sin(ini2[2])+3*lc*math.cos(2*ini2[1])+lc)
    p_theta_dot= (ini2[5]**2 * math.cos(ini2[1])/(ini2[0]**2 * math.sin(ini2[4])**3)) + (kqS/ini2[0]**3)*(lc*math.sin(2*ini2[1])-ls*math.cos(2*ini2[1])*math.sin(ini2[2]))
    p_phi_dot= -(kqS/(2*ini2[0]**3))*ls*math.sin(2*ini2[1])*math.cos(ini2[2])
    
    
    ini3=ini+(h/2)*np.array([r_dot, theta_dot, phi_dot, p_r_dot, p_theta_dot, p_phi_dot])
    
    
    ini=ini+(h/6)*(a+2*ini1+2*ini2+ini3)
    
    
    return list(ini)



trajectory=list(ini)
table=list([[float(r_0), float(theta_0), float(phi_0), float(p_r_0), float(p_theta_0), float(p_phi_0)]])
for t in range(0, 100, 0.0001):
    
    trajectory=Quadrapole(trajectory, .1)
    table.append(trajectory)
    
table_final=np.array(table)


# In[ ]:





# In[ ]:


type(table_final)


# In[ ]:


np.shape(table_final)


# In[ ]:


df = pd.DataFrame(table_final, columns = ['Distance from the origin', 'Polar Angle', 'Azimuthal Angle', 'Radial Momentum', 'Polar Angle Momentum', 'Azimuthal Momentum'])


# In[ ]:


df


# In[ ]:


plt.polar(df['Distance from the origin'], df['Azimuthal Angle'])


# In[ ]:




