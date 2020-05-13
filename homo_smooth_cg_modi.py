#!/usr/bin/python

'''
We solve a fractional order viscoelastic problem with spatially continuous Galerking finite element method and Crank-Nicolson type finite difference method in time. We apply a linear interpolation technique to the part of fractional integration.

The reduced model is given as a vector-valued Volterra integral equation of the second kind.

To simplify the model problem, we set the fourth order identity tensor as the Young modulus for relaxation.

#Set H3 regular solution
'''
from __future__ import print_function
from dolfin import *
import numpy as np
import math
import getopt, sys
import matplotlib.pyplot as plt
from mpltools import annotation

# SS added
parameters['form_compiler']['cpp_optimize'] = True
parameters['form_compiler']['optimize'] = True
parameters["ghost_mode"] = "shared_facet"
set_log_active(False)

# SS added
iMin=2; iMax = 7
jMin=0; jMax = 2

# Define parameters
k=1 # degree of poly
T=1.0   # final time
c=4.0/3.0/sqrt(pi) # 1 over the gamma (5/2)

def usage():
  print("-h   or --help")
  print("-k k      to specify k")
  print("-i i or --iMin  i       to specify iMin")
  print("-j j or --jMin  j       to specify jMin")
  print("-I i or --iMax  i       to specify iMax")
  print("-J j or --jMax  j       to specify jMax")
  print(" ")
  os.system('date +%Y_%m_%d_%H-%M-%S')
  print (time.strftime("%d/%m/%Y at %H:%M:%S"))

# parse the command line
try:
  opts, args = getopt.getopt(sys.argv[1:], "hk:i:I:j:J:",
                   [
                    "help",           # obvious
                    "k=",         # alpha
                    "iMin=",          # iMin
                    "iMax=",          # iMax
                    "jMin=",          # jMin
                    "jMax=",          # jMax
                    ])

except getopt.GetoptError as err:
  # print help information and exit:
  print(err) # will print something like "option -a not recognized"
  usage()
  sys.exit(2)

for o, a in opts:
  if o in ("-h", "--help"):
    usage()
    sys.exit()
  elif o in ("-k"):
    k = int(a)
    print('setting: k = %f;' % k),
  elif o in ("-i", "--iMin"):
    iMin = int(a)
    print('setting:  iMin = %f;' % iMin),
  elif o in ("-I", "--iMax"):
    iMax = int(a)
    print('setting:  iMax = %f;' % iMax),
  elif o in ("-j", "--jMin"):
    jMin = int(a)
    print('setting:  jMin = %f;' % jMin),
  elif o in ("-J", "--jMax"):
    jMax = int(a)
    print('setting:  jMax = %f;' % jMax),
  else:
    assert False, "unhandled option"



#Save data for error
V_error=np.zeros((iMax-iMin+1,jMax-jMin+1), dtype=np.float64)
L2_error=np.zeros((iMax-iMin+1,jMax-jMin+1), dtype=np.float64)

#Define coefficients on the quadrature rule
def Qw(n):
    Bn=[]
    if n==0:
        Bn.append(0.0)
    else:
        Bn.append(n**0.5*(1.5-n)+(n-1.0)**1.5)
        for i in range (1,n):
            Bn.append((n-i-1.0)**1.5+(n-i+1.0)**1.5-2.0*(n-i)**1.5)
        Bn.append(1.0)
    return Bn

#Define exact solution and data terms
f=Expression(("(sin(pi*x[0])*sin(pi*x[1])*ftn1)+ftn2*(3.0/2*pi*pi*sin(pi*x[0])*sin(pi*x[1])-0.5*(x[0]-1)*(x[1]-1)-0.5*(x[0])*(x[1]-1)-0.5*(x[0]-1)*(x[1])-0.5*(x[0])*(x[1]))",\
    "(x[0]*x[1]*(1.0-x[0])*(1.0-x[1])*ftn1)\
+ftn2*(-0.5*pi*pi*cos(pi*x[0])*cos(pi*x[1])-2*(x[0]-1)*(x[0])-(x[1])*(x[1]-1))"),ftn1=0,ftn2=0,degree=5)
u0=Expression(("0.0","0.0"),degree=5)
ux=Expression(("utn*sin(pi*x[0])*sin(pi*x[1])","utn*x[0]*x[1]*(1.0-x[0])*(1.0-x[1])"),utn=1.0,degree=5)

#Define time function values for average
def fTime1(tn):
    return 7.0/2*sqrt(tn)*tn*tn

def fTime2(tn):
    return math.gamma(9.0/2)/24.0*(tn**4)

def uTime(tn):
    return (tn**3.5)


#Strain tensor(Cauchy infinitesimal tensor)
def epsilon(v):
    Dv=nabla_grad(v)
    return 0.5*(Dv+Dv.T)

#Main task
for i in range(iMin,iMax):
    for j in range(jMin,jMax):
        Nh=pow(2,i)
        Nt=2**j
        dt=T/Nt #time step size
        n0=(T/8/dt)
        #print(n0)
        print(' degree of poly = %d, i = %d (to %d), ; j = %d (to %d)' % (k, i, iMax-1,  j, jMax-1))
        mesh=UnitSquareMesh(Nh, Nh) #unit square
        V=VectorFunctionSpace(mesh,'CG',k)  #Define test function space
        dx=Measure('dx')
        ux.utn=0.0
        
        #Define homogeneous Dirichlet boundary condition
        def boundary(x, on_boundary):
            return on_boundary
        bc = DirichletBC(V, Constant((0.0,0.0)), boundary)
    
        #Define linear and bilinear forms and global matrix
        u=TrialFunction(V)
        v=TestFunction(V)
        mass=inner(u,v)*dx
        stiffness=inner(epsilon(u),epsilon(v))*dx

        M= assemble(mass)   #mass matrix
        A= assemble(stiffness)  #stiffness matrix
        B=1.0/dt*M+0.5*(dt**0.5)*c*A    #global matrix
        
        #Approximate an initial condition
        uh=Function(V)
        uh=project(u0,V)    #current solution
        
        oldu=Function(V)
        oldu.assign(uh) #past solution
        U=[]
        U.extend(uh.vector())#Save all solutions in U
        numDof=len(U)   #the number of degree of freedoms
        
        #Apply the quadrature rule
        def Quad(n):
            S=Function(V)
            S.vector()[:]=(Qw(n+1)[0]+Qw(n)[0])*np.array(U[0:numDof])
            for i in range(1,n+1):
                S.vector()[:]+=(Qw(n+1)[i]+Qw(n)[i])*np.array(U[numDof*(i):numDof*(i+1)])
            return S
        
        b = None
        for nt in range(0,int(n0)):            
            #Update exact solution and source terms with respect to given time
            tn=(nt+1.0)*dt;th=(nt)*dt;
            ux.utn=0.5*(uTime(tn)+uTime(th));
            f.ftn1=0.5*(fTime1(tn)+fTime1(th));
            f.ftn2=0.5*(fTime2(tn)+fTime2(th));
            
            #Assemble the right hand side
            L=inner(f,v)*dx+1.0/dt*inner(oldu,v)*dx-0.5*(dt**0.5)*c*inner(epsilon(Quad(nt)),epsilon(v))*dx
            solver = KrylovSolver('cg') #use conjugate gradient method for the linear solver
            prm = solver.parameters
            prm.absolute_tolerance = 1E-14  # from 10
            prm.relative_tolerance = 1E-14
            prm.maximum_iterations = 1000000000
            b = assemble(L, tensor=b)
            bc.apply(B,b)      #apply Dirichlet boundary condition
            solve(B, uh.vector(), b)   
            U.extend(uh.vector())   #restore degrees of freedoms        
            oldu.assign(uh) 
            
      
        #Compute numerical errors
        ux.utn=uTime(tn)  
        err0 = errornorm(ux,uh,'H1')
        err1 = errornorm(ux,uh,'L2')
        V_error[0,j-jMin+1]=Nt; V_error[i-iMin+1,0]=Nh; V_error[i-iMin+1,j-jMin+1]=err0;
        L2_error[0,j-jMin+1]=Nt; L2_error[i-iMin+1,0]=Nh; L2_error[i-iMin+1,j-jMin+1]=err1;
        
#Save errors in text files     
if k==1:        
    np.savetxt("energy_error_linear_smooth.txt",V_error,fmt="%2.3e")
    np.savetxt("L2_error_linear_smooth.txt",L2_error,fmt="%2.3e")
    
elif k==2:
    np.savetxt("energy_error_quad_smooth.txt",V_error,fmt="%2.3e")
    np.savetxt("L2_error_quad_smooth.txt",L2_error,fmt="%2.3e")

#Print error tables    
'''      
# SS altered the following loop limits
#print(V_error)
print ('H1 error')
print('\\begin{tabular}{|c|',end="")
for j in range(jMin,jMax): print('c',end="")
print('|}\\hline')
for j in range(jMin,jMax):
  if j==jMin: print('    ', end="")
  print(' & %10d ' % V_error[0,j-jMin+1], end="")
print('  \\\\\\hline')
for i in range(iMin,iMax):
  print('%4d ' % V_error[i-iMin+1,0], end="")
  for j in range(jMin,jMax):
      print('& %11.4le ' % V_error[i-iMin+1,j-jMin+1], end="")
  print(' \\\\')
print('\\hline\\end{tabular}')

print('This material is obsolete')
print ('L2_error')
print('\\begin{tabular}{|c|',end="")
for j in range(jMin,jMax): print('c',end="")
print('|}\\hline')
for j in range(jMin,jMax):
  if j==jMin: print('    ', end="")
  print(' & %10d ' % L2_error[0,j-jMin+1], end="")
print('  \\\\\\hline')
for i in range(iMin,iMax):
  print('%4d ' % L2_error[i-iMin+1,0], end="")
  for j in range(jMin,jMax):
      print('& %11.4le ' % L2_error[i-iMin+1,j-jMin+1], end="")
  print(' \\\\')
print('\\hline\\end{tabular}')

np.savetxt("energy_error.txt",V_error,fmt="%2.3e")
np.savetxt("L2_error_linear.txt",L2_error,fmt="%2.3e")

l2Diag=[]
h1Diag=[]
m= min(iMax-iMin,jMax-jMin)
#print(m)
for i0 in range(1,m+1):
    l2Diag.append(L2_error[i0,i0])
    h1Diag.append(V_error[i0,i0])
#print((l2Diag))
v1=np.array(l2Diag)
t1=np.log(v1[0:m-2]/v1[1:m-1])
d0=np.mean(t1/np.log(2))

v2=np.array(h1Diag)
t2=np.log(v2[0:m-2]/v2[1:m-1])
d1=np.mean(t2/np.log(2))

print(t1/np.log(2),t2/np.log(2))
print('Numeical convergent order when h=dt: L2 error = %5.4f,  H1 error = %5.4f' %(d0,d1))
'''