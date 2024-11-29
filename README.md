# FeynCalcPy
This is a Python library to do the Feynman ruls calculation semi-automatically. 

## Install
Download the folder of feyncalcpy into your program folder. The library is based on the [SymPy](https://github.com/sympy/sympy). Please make sure that the [SymPy](https://github.com/sympy/sympy) also is installed and import into your program.
Therefore you can run Python:
~~~ python
from feyncalcpy.feynmanruls import *
~~~
to import all the function in the library.
## Usage
### Some functions in SymPy 
The FeynCalcPy support the gamma matrix calculations. You would also need some functions from the SymPy:
~~~ python
from sympy import symbols,S,diag,integrate,cos,sin,diag,pi,nsimplify,factor
from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensor_heads
~~~
define the Tensor index type and Lorentz index :
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
### Spinors
To define a Dirac spinor with mass $m$ and momentum $p$ and polarization $n$, 
~~~ python
p_mu,m_mu,n=symbols('p_mu,m_mu,n')
u=Spinor(p_mu,m_mu,GammaSlash(n))
~~~
where we use SymPy.symbols to define p_mu,m_mu,n to do the algebra. And GammaSlash(n) is one of object in the FeynCalcPy defined by
$
  \text{GammaSlash(n)}=\gamma_\nu n^\nu
$
