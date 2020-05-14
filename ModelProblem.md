# A fractional model of viscoelasticity



The viscoelastic model problem is described as a Volterra integral equation of the second kind. A choice of *power law* type kernel yield the **fractional order model of viscoelasticity**. We solve the fractional order model by using spatial finite element method and *Crank-Nicolson* type finite difference schemes in time. 



### Model Problem

Let $\Omega\subset\mathbb{R}^d$ be a spatial domain for $d=2,3$ and $[0,T]$ be a time interval for $T>0$.  We  consider the following model problem such that

$$
\begin{align*}
\rho\dot{\boldsymbol{w}}(t)-\nabla\cdot{_0}I^{1-\alpha}_t({\underline{\boldsymbol{D}}}\ \underline{{\boldsymbol{\varepsilon}}}(\boldsymbol{w}(t)))&=\boldsymbol{f}(t),&\textrm{ on }(0,T]\times\Omega,\\
{_0}I^{1-\alpha}_t({\underline{\boldsymbol{D}}}\ \underline{{\boldsymbol{\varepsilon}}}(\boldsymbol{w}(t)))\cdot\boldsymbol{n}&=\boldsymbol{g}_N(t),&\textrm{ on }[0,T]\times\Gamma_N,\\
\boldsymbol{w}(t)&=\boldsymbol{0},&\textrm{ on }[0,T]\times\Gamma_D,\\
\boldsymbol{w}(0)&=\boldsymbol{w}_0,&\textrm{ on }\Omega,
\end{align*}
$$


where $\boldsymbol{w}$ is velocity, $_{0}I^{1-\alpha}_t$ is fractional integral of order $1-\alpha$ for $\alpha\in(0,1)$, $\underline{\boldsymbol D}$ is a symmetric positive definite fourth order tensor, $\boldsymbol{\varepsilon}$ is the Cauchy infinitesimal tensor, $\boldsymbol{f}$ is the external body force, $\boldsymbol{n}$ is an outward unit normal vector, $\boldsymbol{g}_n$ is the traction, $\Gamma_N$ is Neumann boundary and $\Gamma_D$ is Dirichlet boundary.  We assume that $\Gamma_D$ and $\Gamma_N$ are disjoint with $\partial \Omega=\Gamma_D\cup\Gamma_N$ and $\Gamma_D$ is of positive measured.



### Weak Formulation

Let $\boldsymbol{V}$ be a test space such that $\boldsymbol{V}=\left\{\boldsymbol{v}\in [H^1(\Omega)]^d\ |\ \boldsymbol v=\boldsymbol 0\ \textrm{ on } {\Gamma_D}\right\}$. Integration by parts yields the following weak formulation.

Find a mapping $\boldsymbol{w}:[0,T]\mapsto\boldsymbol{V}$ such that satisfies 
$$
\begin{align*}
(\rho\dot{\boldsymbol{w}}(t),\boldsymbol{v})+a({_0}I^{1-\alpha}_t\boldsymbol{w}(t),\boldsymbol{v})&=F(t;\boldsymbol{v}),\ \forall t\in(0,T],\\
a(\boldsymbol{w}(0),{\boldsymbol v})&=a(\boldsymbol{w}_0,{\boldsymbol v}),
\end{align*}
$$


 for any $\boldsymbol{v}\in\boldsymbol{V}$ where $(\cdot,\cdot)$ is $L_2$ inner product over $\Omega$ $\left((\cdot,\cdot)_{A}\text{ is $L_2$ inner product over $A\subset\Omega$}\right)$,  $a({\cdot},{\cdot})$ and $F$ are defined by 
$$
a({\boldsymbol{v}},{\boldsymbol w})=\int_\Omega{\underline{\boldsymbol{D}}~ \underline{\boldsymbol{\varepsilon}}(\boldsymbol{v})}:{\underline{\boldsymbol{\varepsilon}}(\boldsymbol{w})}d\Omega\text{ and }F(t;\boldsymbol v)=({\boldsymbol{f}(t)},{\boldsymbol{v}})+({\boldsymbol{g}_N(t)},{\boldsymbol{v}})_{\Gamma_N}.
$$



### Fully Discrete Formulation

Let us consider the finite element space. Suppose that we have conformal meshes with respect to $\Omega$ and let $\boldsymbol{V}^h$ be the Lagrange finite element space of polynomial of degree $k$ in $\boldsymbol{V}$.  Then we define the time step size and introduce time discretisation. 

Let $\Delta t=T/N$ for some $N\in\mathbb{N}$. Define $t_n=n\Delta t$ for $n=0,\ldots,N$ and denote our fully discrete solution by $\boldsymbol{W}_h^{n}$ for $n=0,\ldots,N$. For numerical approximations of fractional integration, we use linear interpolation technique such that
$$
{_0}I^{1-\alpha}_{t_n}\boldsymbol{w}(t)=\frac{\Delta t^{1-\alpha}}{\Gamma(3-\alpha)}\sum_{i=0}^{n}B_{n,i}\boldsymbol{w}(t_i)+\boldsymbol O(\Delta t^2)=\boldsymbol{q}_n(\boldsymbol{w})+\boldsymbol O(\Delta t^2),
$$
where $\Gamma$ is the Gamma function and
$$
\begin{align*}
B_{n,i}=\left\{ \begin{array}{cc}
n^{1-\alpha}({2-\alpha}-n)+(n-1)^{{2-\alpha}},&i=0,\\
(n-i-1)^{{2-\alpha}}+(n-i+1)^{{2-\alpha}}-2(n-i)^{{2-\alpha}},&i=1,\ldots,n-1,\\
1,&i=n.
\end{array}\right.
\end{align*}
$$
Consequently, use of Crank-Nicolson method and numerical integration leads us to obtain a fully discrete formulation as follows:
Find $\boldsymbol{W}_h^n\in \boldsymbol{V}^h$ for $n=0,\ldots,N$ such that satisfying for any $\boldsymbol{v}\in \boldsymbol{V}^h$,
$$
\begin{gather}
\left({\rho\frac{\boldsymbol{W}_h^{n+1}-\boldsymbol{W}_h^n}{\Delta t}},{\boldsymbol{v}}\right)+a\left({\frac{\boldsymbol{q}_{n+1}(\boldsymbol{W}_h)+\boldsymbol{q}_n(\boldsymbol{W}_h)}{2}},{\boldsymbol{v}}\right)=\frac{1}{2}(F(t_{n+1};\boldsymbol{v})+F(t_n;\boldsymbol{v})),\\
a({\boldsymbol{W}_h^0},{\boldsymbol v})=a({\boldsymbol{w}_0},{\boldsymbol v}).
\end{gather}
$$


### Stability and Error Analysis

Our numerical scheme has well-posedness as well as convergence. *A priori* error estimates give us suboptimal and optimal convergence rates. For more details, see this paper and author's thesis (the link will appear here soon).



### Numerical Experiments

Let $\Omega$ be the unit square, $\Gamma_D=\partial\Omega$, $\rho=1$ and $\alpha=1/2$ for simplicity.



We consider two examples; one is of $H^2$ regular in time and the other is of $H^3$. The former example does not satisfy the condition of second order schemes hence it fails to give convergence of the second order. But it could show advantage of the second order method rather than first order methods. On the other hand, the latter is sufficiently smooth with respect to time. We are able to observe optimal error rates.



1. Set 
   $$
   \begin{align*}\boldsymbol{w}(x,y,t)=(t+t^{1.5})\left[\begin{array}{c}
   \sin(\pi x)\sin(\pi y)  \\
   xy(1-x)(1-y)
   \end{array}
   \right].
   \end{align*}
   $$
   Then $\boldsymbol{w}\in C^2(0,T;[C^\infty(\Omega)]^2)\cap W^2_1(0,T;[C^\infty(\Omega)]^2)$ with homogeneous Dirichlet boundary condition.

   

2. Set
   $$
   \begin{align*}\boldsymbol{w}(x,y,t)=t^{2.5}\left[\begin{array}{c}
   \sin(\pi x)\sin(\pi y)  \\
   xy(1-x)(1-y)
   \end{array}
   \right].
   \end{align*}
   $$
   

According to the fully discrete formulation, code implementation has been carried out in **FEniCS**. 

Numerical solutions of Example 1 can be obtained by running **homo_nonsmooth_cg.py**. We can specify mesh sizes with respect to the space as well as time. Once you run it, $H^1$ norm errors and $L_2$ norm errors are stored in *txt* files. In this manner, **homo_smooth_cg_modi.py** deals with numerical approximations for Example 2.



- Run **main.sh** to get numerical solutions and their numerical errors in *txt* files.
- Run  **CG_graphic_frac.py** to illustrate graphs of error convergence rates for linear and quadratic polynomial bases when $h\approx\Delta t$.