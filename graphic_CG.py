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
p1Eu_DG=np.loadtxt("energy_error_linear_DG.txt")
p1Lu_DG=np.loadtxt("L2_error_linear_DG.txt")
'''

'''
p1Eu=np.loadtxt("energy_error_quad.txt")
p1Lu=np.loadtxt("L2_error_quad.txt")
p1Eu_DG=np.loadtxt("energy_error_quad_DG.txt")
p1Lu_DG=np.loadtxt("L2_error_quad_DG.txt")
'''

'''
p1Eu=np.loadtxt("energy_error_linear_smooth.txt")
p1Lu=np.loadtxt("L2_error_linear_smooth.txt")
p1Eu_DG=np.loadtxt("energy_error_linear_DG_smooth.txt")
p1Lu_DG=np.loadtxt("L2_error_linear_DG_smooth.txt")
'''


p1Eu=np.loadtxt("energy_error_quad_smooth.txt")
p1Lu=np.loadtxt("L2_error_quad_smooth.txt")
p1Eu_DG=np.loadtxt("energy_error_quad_DG_smooth.txt")
p1Lu_DG=np.loadtxt("L2_error_quad_DG_smooth.txt")


iMax=len(p1Eu)
jMax=len(p1Eu[0])

iMin=1
jMin= 1
#print(N,N2,M,M2)
xh=p1Eu[1::,0]
#xh2=p2Eu[1::,0]
#print(np.shape(xh2),np.shape(p2Eu_DG.diagonal()[1::]))


print('------------------------- CG ------------------')
print ('H1 error')
print('\\begin{tabular}{|c|',end="")
for j in range(jMin,jMax): print('c',end="")
print('|}\\hline')
for j in range(jMin,jMax):
  if j==jMin: print('    ', end="")
  print(' & 1/%d ' % p1Eu[0,j-jMin+1], end="")
print('  \\\\\\hline')
for i in range(iMin,iMax):
  print('1/%d ' % p1Eu[i-iMin+1,0], end="")
  for j in range(jMin,jMax):
      print('& %11.3le ' % p1Eu[i-iMin+1,j-jMin+1], end="")
  print(' \\\\')
print('\\hline\\end{tabular}')

print('This material is obsolete')
print ('L2_error')
print('\\begin{tabular}{|c|',end="")
for j in range(jMin,jMax): print('c',end="")
print('|}\\hline')
for j in range(jMin,jMax):
  if j==jMin: print('    ', end="")
  print(' &1/%d ' % p1Lu[0,j-jMin+1], end="")
print('  \\\\\\hline')
for i in range(iMin,iMax):
  print('1/%d ' % p1Lu[i-iMin+1,0], end="")
  for j in range(jMin,jMax):
      print('& %11.3le ' % p1Lu[i-iMin+1,j-jMin+1], end="")
  print(' \\\\')
print('\\hline\\end{tabular}')
print('------------------------ DG ------------------')
print ('H1 error')
print('\\begin{tabular}{|c|',end="")
for j in range(jMin,jMax): print('c',end="")
print('|}\\hline')
for j in range(jMin,jMax):
  if j==jMin: print('    ', end="")
  print(' & 1/%d ' % p1Eu_DG[0,j-jMin+1], end="")
print('  \\\\\\hline')
for i in range(iMin,iMax):
  print('1/%d ' % p1Eu_DG[i-iMin+1,0], end="")
  for j in range(jMin,jMax):
      print('& %11.3le ' % p1Eu_DG[i-iMin+1,j-jMin+1], end="")
  print(' \\\\')
print('\\hline\\end{tabular}')

print('This material is obsolete')
print ('L2_error')
print('\\begin{tabular}{|c|',end="")
for j in range(jMin,jMax): print('c',end="")
print('|}\\hline')
for j in range(jMin,jMax):
  if j==jMin: print('    ', end="")
  print(' & 1/%d ' % p1Lu_DG[0,j-jMin+1], end="")
print('  \\\\\\hline')
for i in range(iMin,iMax):
  print('1/%d ' % p1Lu_DG[i-iMin+1,0], end="")
  for j in range(jMin,jMax):
      print('& %11.3le ' % p1Lu_DG[i-iMin+1,j-jMin+1], end="")
  print(' \\\\')
print('\\hline\\end{tabular}')

'''
plt.figure(1)
line_D1=plt.plot(np.log(xh),np.log(p1Eu.diagonal()[1::]),label=r'(CG)$H^1$',ls='-.',marker="D")
line_D2=plt.plot(np.log(xh),np.log(p1Lu.diagonal()[1::]),label=r'(CG)$L_2$',ls='--',marker="o")
line_D3=plt.plot(np.log(xh),np.log(p1Eu_DG.diagonal()[1::]),label=r'(DG)$H^1$',marker='^')
line_D4=plt.plot(np.log(xh),np.log(p1Lu_DG.diagonal()[1::]),label=r'(DG)$L_2$',marker="s")


annotation.slope_marker((4.5,-10.8),-1,invert=True)
annotation.slope_marker((4.5,-16.3),-1.5,invert=True)

plt.xlabel(r'$\log(1/h)$')
plt.legend(loc='best')

plt.ylabel('log(error)')
plt.title('Numerical results of linear basis')
plt.show()
'''

plt.figure(2)
line_D1=plt.plot(np.log(xh),np.log(p1Eu.diagonal()[1::]),label=r'(CG)$H^1$',ls='-.',marker="D")
line_D2=plt.plot(np.log(xh),np.log(p1Lu.diagonal()[1::]),label=r'(CG)$L_2$',ls='--',marker="o")
line_D3=plt.plot(np.log(xh),np.log(p1Eu_DG.diagonal()[1::]),label=r'(DG)$H^1$',marker='^')
line_D4=plt.plot(np.log(xh),np.log(p1Lu_DG.diagonal()[1::]),label=r'(DG)$L_2$',marker="s")


annotation.slope_marker((4.5,-15),-2,invert=True)
annotation.slope_marker((4.5,-16.7),-2,invert=True)

plt.xlabel(r'$\log(1/h)$')
plt.legend(loc='best')

plt.ylabel('log(error)')
plt.title('Numerical results of quadratic basis')
plt.show()
