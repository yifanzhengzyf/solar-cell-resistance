# -*- coding: utf-8 -*-
"""
@author: Sameer Kesava

In version 2, final data also has efficiency and cell length added.

Run RintrinsicDetermination_v3.py to obtain Rintrinsic"""


Rintrinsic = 1.9 #Obtained from RintrinsicDetermination_v3.py
#Rseries = 2.52248
Rshunt = 579.563 #Ohms-cm2
Jph = 10.6305*10**(-3) #A/cm2
n = 3.14219 #ideality factor
Jo = 4.33538*10**(-8)  #A/cm2
Rsheet = 20 #Ohms/sq
Vth = 25.6825 * 10**(-3) #Thermal voltage-eV
#Voc = 1.0013 #Volt
cellwidth = 0.3 #measured cell width (cm). Starting value
l = 0.3 #measured cell length (cm)
Psun = 0.1 #W/cm2
#%%

""" importing packages"""
import math
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
#from scipy.integrate import odeint
#from scipy.special import lambertw
#import sys
from scipy import optimize
import time  
#%%
def Vocfunc(x):
    
    return Jph-Jo*(math.exp(x/(n*Vth))-1)- x/(Rshunt)

soln = optimize.root(Vocfunc, 0.98, method = 'lm')

Voc=soln.x[0]

#%%
def fun(x,*argv):
    """Function representing diode-equation with input voltage and cell
    intrinsic series resistance as arguments"""
    return x-Jph+Jo*(math.exp((argv[0]+argv[1]*x)/(n*Vth))-1)\
    +(argv[0]+argv[1]*x)/(Rshunt)

def Jx(v,i,r):
    """ Returns Current density J by solving the diode equation while passing
    input voltage and cell intrinsic series resistance as arguments to the 
    diode-equation function"""
    sol=optimize.root(fun, i, args = (v,r), method = 'lm')
    return sol.x[0]

#%%
def VIlist(W_V_input, Rsheet, Rseries, dx):
    """ Inputting voltage and current values at x=0 along cell width.
    Giving sheet resistance of electrode and series resistance of active layer.
    This function returns values of I and V showing variation along the width of the cell
    for a starting V at x = 0"""
   
    Vlist = np.array([W_V_input[1]]) 
    
    Ilist = np.array([0])
    
    num = int(round(W_V_input[0]/dx,0))
    
    for i in range(num):   
        Inext = Ilist[i]+ Jx(Vlist[i],Ilist[i], Rseries)*dx*l
        Vnext = Vlist[i] - Rsheet*dx/l*Inext
        Ilist = np.append(Ilist, Inext)
        Vlist = np.append(Vlist, Vnext)
    
    return [Vlist, Ilist]
    

#%%
"""Calculating the Vmp and Imp of experimental cell using the above function VIlist"""
exp_cell_vilist = np.empty([0,3])

dx = 0.001

for i in np.arange(0, Voc, 0.01):
    calc = VIlist([cellwidth,i], Rsheet, Rintrinsic, dx)
    exp_cell_vilist = np.append(exp_cell_vilist, [[i,calc[0][-1], calc[1][-1]]], axis = 0)
    
""" exp_cell_vilist contains an array of input V, output V and output I"""    
    
exp_cell_maxpw_indx = np.argmax(exp_cell_vilist[:,1]*exp_cell_vilist[:,2])

exp_cell_maxV_x_zero = exp_cell_vilist[:,0][exp_cell_maxpw_indx]
exp_cell_maxV_x_width = exp_cell_vilist[:,1][exp_cell_maxpw_indx]
    

#%%
"""" This code determines the minimum V at x=0 for which the I is positive for x=width"""
dx = 0.01
#Starting integration step

w = 12
#Required maximum width

Num = int(round(w/dx,0))    
#Number of data points for integration

Vstart=0      
#count = 1
Vstartlist = np.empty([0,2])
loop=True

#%%    
while loop:
    
    
    if Jx(Vstart, Jph, Rintrinsic)<0 or Vstart > float(format(Voc,'0.3f')):
        #print('J negative')
        print('\nVoc reached')
        #Vstartlist = np.append(Vstartlist,[[dx*count, Vstart]], axis=0)
        break
    
    else:
        Ilist = np.array([0])
        Vlist = np.array([Vstart])
        print(dx)
    
    for j in range(Num):
        
        Inext = Ilist[j]+ Jx(Vlist[j], Ilist[j], Rintrinsic)*dx*l
        Vnext = Vlist[j] - Rsheet*dx/l*Inext
               
        if Vnext < 0:
            #prevV=Vnext
            Vstartlist = np.append(Vstartlist,[[dx*(j), Vstart]], axis=0)
            print(dx*(j), Vstart, '\n') 
            Vstart = Vstart + abs(Vnext)
            if Vstart > float(format(Voc, '0.3f')) and dx >= 0.005:
                dx = dx*0.5
                Num = int(round(w/dx,0))  
                Vstart = Vstart - abs(Vnext)
                break
            else:
                break
       
             
        else:
            Ilist = np.append(Ilist, Inext)
            Vlist = np.append(Vlist, Vnext)
#%%
np.savetxt(r"VvsW_x_0_Rsheet-%d_Rseries-%0.1f.csv" %(Rsheet, Rintrinsic), Vstartlist, delimiter=",",\
           header = "cell width (cm), Minimum V @ x=0 for which device is at short-circuit",\
           comments = "")
#%%
"""The maximum width above which increase in width would lead to no added value
to the solar cell, calculated as V approximately equal to Voc, upto 3 digits of precision, at x = 0"""
Maxwidth1 = Vstartlist[[-1]][0,0]
print (Maxwidth1)

"""The maximum width above which increase in width would lead to no added value
to the solar cell, calculated as V approximately equal to Voc, upto 2 digits of precision, at x = 0"""
Maxwidth2 = Vstartlist[[np.where(Vstartlist[:,1]<=float(format(Voc, '0.2f'))-0.01)[0][-1]]][0,0]
print (Maxwidth2)
#%%
widthlist = np.arange(cellwidth, w+cellwidth, cellwidth)
Vstartlistcopy = Vstartlist
#%%
widthvoltlist = np.empty([0,2])
#%%

"""Creating an array with elements as multiples of actual cell width and appending
the corresponding starting voltage for IV plotting and Max power determination"""
for i in widthlist:
    for j in Vstartlistcopy:
        
        if i>j[0]:
            continue        
        
        elif i<=j[0]:
            
            widthvoltlist = np.append(widthvoltlist, [j], axis=0)
        
            Vstartlistcopy = np.delete(Vstartlistcopy, np.where(Vstartlistcopy[:,0]<=j[0]),0)
            
            break
    
    if i==widthlist[-1]:
        
        widthvoltlist = np.append(widthvoltlist, [j], axis=0)
        
#%%
def IVdata(wv, dx):
    """ returns V-I list for a given width, and voltage at x=0 such that V at x=width is zero.
    dx is the integration step"""
    #global step   
    """looks like when you have global variables, parallel processing will not work"""
    step = (Voc - wv[1])/100
    
    bias_x_zero = wv[1]
    
    bias_x_width = np.array([])
    I_x_width = np.array([])
    
    v_out_diff = 0
    
    while bias_x_zero < Voc:
        
        ivdata = VIlist([wv[0], bias_x_zero], Rsheet, Rintrinsic, dx)
        
        if abs(v_out_diff - ivdata[0][-1]) > 0.1:
            bias_x_zero = bias_x_zero - 0.5 * step
            step = step*0.5
            continue
        
        else:
            v_out_diff = ivdata[0][-1]
            bias_x_width = np.append(bias_x_width, ivdata[0][-1])
            I_x_width = np.append(I_x_width, ivdata[1][-1])
            bias_x_zero = bias_x_zero + step
        
    return [bias_x_width, I_x_width, wv[0]]
#%%
def JVdataplot(data_input):            
    fig1, ax1 = plt.subplots()
    ax1.plot(data_input[0], data_input[1]*1000/(data_input[2]*l), 'r-', lw=2,\
    label = r'Cell width = %0.1f cm, R sheet = %d $\Omega$/sq' %(data_input[2], Rsheet))
    #ax1.plot(Voutlist, Jlist*1000, 'b--', lw=2, label = r'Diode equation,\
     #Rseries: %0.1f $\Omega$-cm$^2$'%Rseries)
    ax1.legend(loc ='lower left', prop={'size':8})
    ax1.axis([np.min(data_input[0]), np.max(data_input[0]),\
    np.min(data_input[1])*1000/(data_input[2]*l), np.max(data_input[1])*1000/(data_input[2]*l)])
    plt.xlabel('Voltage (V) ')
    plt.ylabel(r'Current density (mA/cm$^2$)')
    #plt.savefig('IV_RintrinsicDetermination v2.jpeg', format = 'jpeg', dpi=600)
#%% 
def Max_Power(wv, dx):
    """ returns Vmp (V), Imp (A), Pmp (W), Efficiency(%), cell width and cell length
    which is constant"""
    data_input = IVdata(wv, dx)
    index = np.argmax(data_input[0]*data_input[1])
    return [data_input[0][index],data_input[1][index],\
    data_input[0][index] * data_input[1][index],\
    (data_input[0][index] * data_input[1][index])/(data_input[2]*l*Psun)*100, data_input[2],l]
    

#%%    
dx2 = 0.01 
max_power_list = np.empty([0,6])

start_time = time.time()
for i in widthvoltlist:
    #temp_ivdata = IVdata(i, dx2)cr
    #temp_max_pw = Max_Power(temp_ivdata)
    #temp_max_pw.append(temp_ivdata[2],l)
    max_power_list = np.append(max_power_list, [Max_Power(i, dx2)], axis=0)

end_time = time.time()
#%%

"""NOT WORKING! Need to fix this at some point"""

"""Using parallel processing"""
"""
from joblib import Parallel, delayed


Parallel(n_jobs = 4)(delayed(Jx)(i, Jph, Rintrinsic) for i in np.arange(0, Voc,0.01))



start_time2 = time.time()

listofMax_Pwr = Parallel(n_jobs = 4)(delayed(Max_Power)(i, dx2) for i in widthvoltlist)

end_time2 = time.time()

print('Time taken for execution is %f s' %(end_time - start_time))
   
max_power_list2 = np.empty([0,6])
for i in listofMax_Pwr:
    max_power_list2 = np.append(max_power_list2, [i], axis=0)
"""    
#%%

#temp_max_pw.remove(widthvoltlist[-1][0]) 
indx = np.where(widthlist>max_power_list[-1][4])
for i in widthlist[indx]:
    #temp_max_pw.append(i)
    max_power_list = np.append(max_power_list, [max_power_list[-1]], axis=0)
    #temp_max_pw.remove(i)
    np.put(max_power_list[-1], [3, 4], [max_power_list[-1][2]/(i*l*Psun)*100, i])
    

   
#%%
np.savetxt(r"MaxPvsW_Rsheet-%d_Rseries-%0.1f.csv" %(Rsheet, Rintrinsic), max_power_list, delimiter=",", header = "Max V, Max I (A),\
Max Power (W), Efficiency (%), Cell Width (cm), Cell length (cm)", comments = "")
#%%
fig1, ax1 = plt.subplots()
ax1.plot(max_power_list[:, 4], np.round(max_power_list[:,0],3), 'r-', lw=2,\
label = r'R Intrinsic = %0.1f $\Omega$-cm$^2$, R sheet = %d $\Omega$/sq' %(Rintrinsic, Rsheet))
#ax1.plot(Voutlist, Jlist*1000, 'b--', lw=2, label = r'Diode equation,\
 #Rseries: %0.1f $\Omega$-cm$^2$'%Rseries)
ax1.legend(loc ='upper right', prop={'size':10})
ax1.axis([0, np.max(max_power_list[:,4]),0, Voc])
plt.xlabel('Cell Width (cm) ')
plt.ylabel('Voltage at Maximum Power Point (V)')
#plt.savefig('VmaxvsCellwidth', format = 'jpeg', dpi = 600)

fig2, ax2 = plt.subplots()
ax2.plot(max_power_list[:, 4], np.round(max_power_list[:,1]*1000,3), 'b-', lw=2,\
label = r'R Intrinsic = %0.1f $\Omega$-cm$^2$, R sheet = %d $\Omega$/sq' %(Rintrinsic, Rsheet))
#ax1.plot(Voutlist, Jlist*1000, 'b--', lw=2, label = r'Diode equation,\
 #Rseries: %0.1f $\Omega$-cm$^2$'%Rseries)
ax2.legend(loc ='upper right', prop={'size':10})
ax2.axis([0, np.max(max_power_list[:,4]),0, np.max(max_power_list[:,1])*1000*2])
plt.xlabel('Cell Width (cm) ')
plt.ylabel('Current at Maximum Power Point (mA)')

fig3, ax3 = plt.subplots()
ax3.plot(max_power_list[:, 4], max_power_list[:,2]*1000, 'g-', lw=2,\
label = r'R Intrinsic = %0.1f $\Omega$-cm$^2$, R sheet = %d $\Omega$/sq' %(Rintrinsic, Rsheet))
#ax1.plot(Voutlist, Jlist*1000, 'b--', lw=2, label = r'Diode equation,\
 #Rseries: %0.1f $\Omega$-cm$^2$'%Rseries)
ax3.legend(loc ='upper right', prop={'size':10})
ax3.axis([0, np.max(max_power_list[:,4]),0, np.max(max_power_list[:,2])*1000*2])
plt.xlabel('Cell Width (cm) ')
plt.ylabel('Maxium Power (mW)')

fig4, ax4 = plt.subplots()
ax4.plot(max_power_list[:, 4], max_power_list[:,3], 'k-', lw=2,\
label = r'R Intrinsic = %0.1f $\Omega$-cm$^2$, R sheet = %d $\Omega$/sq' %(Rintrinsic, Rsheet))
#ax1.plot(Voutlist, Jlist*1000, 'b--', lw=2, label = r'Diode equation,\
 #Rseries: %0.1f $\Omega$-cm$^2$'%Rseries)
ax4.legend(loc ='upper right', prop={'size':10})
ax4.axis([0, np.max(max_power_list[:,4]),0, 1.5*np.max(max_power_list[:,3])])
plt.xlabel('Cell Width (cm) ')
plt.ylabel('Efficiency (%)')
#%%

fig1.savefig(r'VmaxvsCellwidth_Rsheet_%d_ohms.png' %(Rsheet), format = 'png', dpi = 100)
fig2.savefig(r'ImaxvsCellwidth_Rsheet_%d_ohms.png'%(Rsheet), format = 'png', dpi = 100)
fig3.savefig(r'PmaxvsCellwidth_Rsheet_%d_ohms.png'%(Rsheet), format = 'png', dpi = 100)
fig4.savefig(r'EffvsCellwidth_Rsheet_%d_ohms.png'%(Rsheet), format = 'png', dpi = 100)

#%%

    
    




    

    
    
    
    
    
    
    
    
    

     
         
            
