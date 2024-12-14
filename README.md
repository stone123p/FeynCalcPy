# FeynCalcPy
This is a Python library to do the Feynman rules calculation semi-automatically. 

## Install
Download the folder of feyncalcpy into your program folder. The library is based on the [SymPy](https://github.com/sympy/sympy). Please make sure that the [SymPy](https://github.com/sympy/sympy) also is installed and imported into your program.
Therefore you can run Python:
~~~ python
from feyncalcpy.feynmanruls import *
~~~
To import all the functions in the library.
## Usage
### Some functions in SymPy 
The FeynCalcPy supports the gamma matrix calculations. You would also need some functions from the SymPy:
~~~ python
from sympy import symbols,S,diag,integrate,cos,sin,diag,pi,nsimplify,factor
from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensor_heads
~~~
Define the Tensor index type and Lorentz index :
~~~ python
Lorentz = TensorIndexType('Lorentz', dummy_name='L', dim=4)
mu,nu,alpha,beta,rho,sigma = tensor_indices('mu,nu,alpha,beta,rho,sigma', Lorentz)
~~~

### Gamma Matrix 5, $\gamma^5$ 
It's convenient to define $\gamma^5$ and chiral projection operators:
~~~ python
G5=GammaMatrix5();
PL=(S.One-G5)/2
PR=(S.One+G5)/2
~~~

### Spinors, GammaSlash
To define a Dirac spinor with mass $m1$ and momentum $p1$ and polarization $n$, 
~~~ python
p_1,m_1,n=symbols('p_1,m_1,n')
u1=Spinor(p_1,m_1,GammaSlash(n))
~~~
where we use SymPy.symbols to define p_mu,m_mu,n to do the algebra. And GammaSlash(n) is meant $\gamma^\alpha n_\alpha$. 
~~~ python
>>> GammaSlash(n)
gs(n)
~~~

### GammaMatrix
To define the gamma matrix with a specific superscript index, $\gamma^\mu$, you can 
~~~ python
GammaMatrix(mu) # for gamma^mu
~~~
on the other hand, you can also define the gamma matrix with subscript, $\gamma_\mu$
~~~ python
GammaMatrix(-mu) # for gamma_mu
~~~

### Gamma Matrix calculator, GC 
Using GC is a object to represent $\bar u(p_1)\gamma^\mu\gamma^\nu(\gamma^5+1)u(p_2)$, 
~~~ python
>>> mu,nu = tensor_indices('mu,nu', Lorentz)
>>> G5=GammaMatrix5();
>>> GC(u(p_1),GammaMatrix(mu)*GammaMatrix(nu)*(G5+S.One),u(p_2))
GC(u(p_1),gamma^mu*gamma^nu*(gamma^5+1),u(p_2))
~~~
you can do the complex conjugation
~~~ python
>>> y=GC(u(p_1),GammaMatrix(mu)*GammaMatrix(nu)*(G5+S.One),u(p_2)).conj()
>>> y
GC(u(p_2),(G5-S.One)*GammaMatrix(nu)*GammaMatrix(mu),u(p_1))
~~~

It supports the multiplication between GC. It returns different type of results depending on the given spinors. For example, If one of the momentum of the spinor are the same and contracts to each other, it will automatically do projections concerning the given spinor
~~~ python
>>> GC(u(p_1),H1,u(p_2))*GC(u(p_2),H2,u(p_3))
GC(u(p_1),H1*(gs(p_2)+m2)*H2,u(p_3))
~~~
If two momentum contracts respectively, it can do the trace automatically
~~~ python
>>> GC(u(p_1),H1,u(p_2))*GC(u(p_2),H2,u(p_1))
trace((gs(p_1)+m_1)*H1*(gs(p_2)+m_2)*H2)
~~~
### Trace 
The function trace(...) is also a function in FeynCalcPy.
~~~ python
>>> trace(GammaMatrix(mu)*GammaMatrix(nu))
~~~
The output of the above snippet is as follows -
~~~ math
 4 metric^{\mu\nu}
~~~

~~~ python
>>> trace(GammaMatrix(mu)*GammaMatrix(nu)*GammaMatrix(alpha)*GammaMatrix(beta))
~~~
The output of the above snippet is as follows 
~~~ math
4 ((-1)metric^{\mu\alpha}metric^{\nu\beta}+metric^{\mu\beta}metric^{\nu\alpha}+metric^{\mu\nu}metric^{\beta\alpha})
~~~
~~~
>>> trace(GammaMatrix5*GammaMatrix(mu)*GammaMatrix(nu)*GammaMatrix(alpha)*GammaMatrix(beta))
~~~
The output of the above snippet is as follows 
~~~ math
4i\epsilon^{\mu\nu\alpha\beta}
~~~
The GammaSlash object can also participate with the trace function, for example,
~~~
>>> p = symbols('p')
>>> trace(GammaMatrix(mu)*GammaSlash(p))
~~~
The output of the above snippet is as follows 
~~~ math
 4 metric^{\mu L_0}p_{L_0}
~~~
You can use contract_metric to simplify further
~~~
>>> p = symbols('p')
>>> Lorentz = TensorIndexType('Lorentz', dummy_name='L', dim=4)
>>> Lg = Lorentz.metric
>>> trace(GammaMatrix(mu)*GammaSlash(p)).contract_metric(Lg)
~~~
The output of the above snippet is as follows 
~~~ math
 4p^\mu
~~~
where $p^\mu$ is a [sympy.tensor.Tensor](https://docs.sympy.org/latest/modules/tensor/tensor.html) object.
### Levi-Civita Symbol, Kronecker Delta and their symmetry
You define a four dimensional LeviCivita tensor by runing 
~~~ python
 >>> m0,m1,m2,m3,m4 = tensor_indices('m0,m1,m2,m3,m4', Lorentz)
 >>> LC=LeviCivitaT('LeviCivita',self_indices=(m0,m1,m2,m3))
~~~
~~~ math
 \epsilon^{m_0m_1m_2m_3}
~~~
It has some useful contraction, for examples
~~~python
>>> m0,m1,m2,m3,m4,m5,m6 = tensor_indices('m0,m1,m2,m3,m4,m5,m6', Lorentz)
>>> LeviCivitaT(m0,m1,m2,m3)*LeviCivitaT(-m0,-m1,-m4,-m5)
~~~
~~~ math
 -2((-1)metric_{m5}^{m2}metric_{m4}^{m3}+metric_{m5}^{m3}metric_{m4}^{m2}))
~~~
or 
~~~python
>>> m0,m1,m2,m3,m4,m5,m6 = tensor_indices('m0,m1,m2,m3,m4,m5,m6', Lorentz)
>>> LeviCivitaT(m0,m1,m3,m4)*LeviCivitaT(-m0,-m1,-m3,-m5)
~~~
~~~ math
 -6(metric_{m5}^{m4})
~~~
Because of that LeviCivitaT is a anti-symmetry tensor, the contraction with a symmetry tensor leads to 0, forexample
~~~ python
>>>T=tensor_heads('T', [Lorentz]*2,TensorSymmetry.fully_symmetric(2))
>>>t=(LeviCivitaT('LeviCivita',self_indices=(-m0,-m1,-m3,-m5))*T(m0,m1))
>>>t
~~~
~~~ math
\epsilon_{L0L1m3,m5}T^{L0L1}
~~~
you can use, [.canon_bp()](https://docs.sympy.org/latest/modules/tensor/tensor.html), canonicalization under monoterm symmetrie to canonicalize:
~~~ python
>>>t.canon_bp()
~~~
0
### Effective Hamiltonian 
You can follow your requirement to write done the effective Hamiltonian, for instance, for a four fermion interaction, the Hamiltonian
~~~ math
 H\sim (\bar u(p_1)\gamma^\mu L u(p_2))(\bar u(p_3)\gamma_\mu L u(p_4))+h.c. .
~~~ 
In order to calculate amplitude squared, you can follow
~~~ python
g1=GammaMatrix(mu)
g2=GammaMatrix(nu)
conj_g1=GammaMatrix(-mu)
conj_g2=GammaMatrix(-nu)

u1=Spinor(p_1,m_1,0)
u2=Spinor(p_2,m_2,0)
u3=Spinor(p_3,m_3,0)
u4=Spinor(p_4,m_4,0)
Hv1=[GC(u1,g1*PL,u2),GC(u3,conj_g1*PL,u4)]
Hv2=[GC(u1,g2*PL,u2),GC(u3,conj_g2*PL,u4)]
~~~
Then amplitude squared can be obtained by
~~~ python
Z1=Hs[0].conj()*Hv[0]
Z2=Hs[1].conj()*Hv[1]
Z3=Z1*Z2
~~~
