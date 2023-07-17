Hello, I am a FEniCS user interested in *continuous/discontinuous Galerkin finite element methods*.  
Mainly, during my PhD research, I have contributed to code implementations for solving **viscoelastic problems**.  
My works focused on numerical solutions to **dynamic viscoelastic problems**. 
Code implementations deal with simple wave equations, elastic problems, vicoelastic problems modelled by *generalised Maxwell solid* or *fractional order model*.  
Numerical schemes are based on *spatial finite element methods* (CG/DG) and *Crank-Nicolson method* for time discretisation.


In this numerical experiment, we consider two examples to observe *suboptimal* errors and *optimal* errors with respect to time (we take into account $H^2$ and $H^3$ exact solutions). For numerical approaches of fractional integration, a **linear interpolation technique** is employed and gives second order accurracy.


For more details of the model problem, see **ModelProblem.pdf** or **ModelProblem.html**.

All codes are constructed in Python 2.7.12 and FEniCS(dolfin) 2016.2.0.
- **homo_nonsmooth_cg.py**: Solve the less smooth problem to show suboptimal results.
- **homo_smooth_cg_modi.py**: Solve the $H^3$ regular problem to show optimal results.
- **CG_graphic_frac.py**: Illustrate graphs of numerical errors on *log-log* scale when $h\approx\Delta t$.
- **main.sh**: Main task to run (consider both examples with linear polynomial basis as well as quadratic).


If you have any inquiries, please contact me at email yongseok.jang@brunel.ac.uk or yongseok20007717@gmail.com.
If you are interested in my research, please visit https://yongseokmath.wordpress.com/.
