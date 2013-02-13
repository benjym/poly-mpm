from numpy import *
import matplotlib.pyplot as plt


E_0 = 1e7 # from simulation
n=array((2,3,5,10,15,25,50,75,150,250))

T_1=array((220,270,310,340,360,370,380,385,390,390))*(10**-4) # rho = 2650
#T_2=array((9,11,13,15,17,16,17))*(10**-3) # rho = 500

f_1=1./T_1
#f_2=1/T_2

k_1=.5*.5*1.*2650.*4*(pi**2)*(f_1**2)
#k_2=.5*.5*1.*500.*4*(pi**2)*(f_2**2)

k_parallel = E_0*n
k_series = 1./(n/E_0)


plt.plot(n,k_1/E_0,'o',label='measured, rho=2650')
#plt.plot(n,k_2/E_0,'o',label='measured, rho=500')
#plt.plot(n,k_series/E_0,'x',label='series')
#plt.plot(n,k_parallel/E_0,'v',label='parallel')
#plt.legend()


plt.xlabel('number of material points in loading direction')
plt.ylabel('spring constant/E')
plt.savefig('../im/spring.png')

print k_1[-1]/E_0
