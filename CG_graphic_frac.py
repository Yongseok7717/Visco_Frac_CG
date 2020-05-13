#!/usr/bin/python
from __future__ import print_function
from dolfin import *
import numpy as np
import math
import getopt, sys
import matplotlib.pyplot as plt
from mpltools import annotation

'''
p1Eu=np.loadtxt("energy_error_linear.txt")
p1Lu=np.loadtxt("L2_error_linear.txt")
p2Eu=np.loadtxt("energy_error_quad.txt")
p2Lu=np.loadtxt("L2_error_quad.txt")
'''


p1Eu=np.loadtxt("energy_error_linear_smooth.txt")
p1Lu=np.loadtxt("L2_error_linear_smooth.txt")
p2Eu=np.loadtxt("energy_error_quad_smooth.txt")
p2Lu=np.loadtxt("L2_error_quad_smooth.txt")



xh=p1Eu[1::,0]

'''
plt.figure(1)
line_D1=plt.plot(np.log(xh),np.log(p1Eu.diagonal()[1::]),label=r'$H^1(k=1)$',ls='-.',marker="D")
line_D2=plt.plot(np.log(xh),np.log(p1Lu.diagonal()[1::]),label=r'$L_2(k=1)$',ls='--',marker="o")
line_D3=plt.plot(np.log(xh),np.log(p2Eu.diagonal()[1::]),label=r'$H^1(k=2)$',marker='^')
line_D4=plt.plot(np.log(xh),np.log(p2Lu.diagonal()[1::]),label=r'$L_2(k=2)$',marker="s")


annotation.slope_marker((4.5,-5.4),-1,invert=True)
annotation.slope_marker((4.5,-11.5),-1.5,invert=True)

plt.xlabel(r'$\log(1/h)$')
plt.legend(loc='best')

plt.ylabel('log(error)')
plt.title('Numerical results')
plt.show()
'''


plt.figure(2)
line_D1=plt.plot(np.log(xh),np.log(p1Eu.diagonal()[1::]),label=r'$H^1(k=1)$',ls='-.',marker="D")
line_D2=plt.plot(np.log(xh),np.log(p1Lu.diagonal()[1::]),label=r'$L_2(k=1)$',ls='--',marker="o")
line_D3=plt.plot(np.log(xh),np.log(p2Eu.diagonal()[1::]),label=r'$H^1(k=2)$',marker='^')
line_D4=plt.plot(np.log(xh),np.log(p2Lu.diagonal()[1::]),label=r'$L_2(k=2)$',marker="s")


annotation.slope_marker((4.5,-10.8),-1,invert=True)
annotation.slope_marker((4.5,-16.6),-2,invert=True)

plt.xlabel(r'$\log(1/h)$')
plt.legend(loc='best')

plt.ylabel('log(error)')
plt.title('Numerical results')
plt.show()
