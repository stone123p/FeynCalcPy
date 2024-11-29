from sympy import symbols,diag,integrate,cos,sin,diag,pi,nsimplify,factor
from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensor_heads

from feyncalcpy.feynmanruls import *

Lorentz = TensorIndexType('Lorentz', dummy_name='L', dim=4)

mu,nu,alpha,beta,rho,sigma = tensor_indices('mu,nu,alpha,beta,rho,sigma', Lorentz)

n1,n2,n3,n=symbols("n1,n2,n3,n")

E_nub,p_nub,theta_nub,phi_nub=symbols(r"E_{\bar\nu} p_{\bar\nu} \theta_{\bar\nu} \phi_{\bar\nu}")
E_nu,p_nu,theta_nu,phi_nu=symbols(r"E_\nu p_\nu \theta_\nu \phi_\nu")
E_e,p_e,m_e,theta_e,phi_e=symbols(r"E_e p_e m_e \theta_e \phi_e")
E_mu,p_mu,m_mu=symbols(r"E_\mu p_\mu m_mu")

MatrixReps={
     r'p_e':[E_e,p_e*cos(phi_e)*sin(theta_e),p_e*sin(phi_e)*sin(theta_e),p_e*cos(theta_e)]
    ,r'n':[0,n1,n2,n3]
    ,r'p_{\bar\nu}':[E_nub,p_nub*cos(phi_nub)*sin(theta_nub),p_nub*sin(phi_nub)*sin(theta_nub),p_nub*cos(theta_nub)]
    ,r'p_\mu':[m_mu,0,0,0]
    ,r'p_\nu':[E_nu,p_nu*cos(phi_nu)*sin(theta_nu),p_nu*sin(phi_nu)*sin(theta_nu),p_nu*cos(theta_nu)]
}

G5=GammaMatrix5()
PL=(S.One-G5)/2
PR=(S.One+G5)/2

u1=Spinor(p_e,m_e,0)
u2=Spinor(p_nu,0,0)
u3=Spinor(p_nub,0,0)
u4=Spinor(p_mu,m_mu,GammaSlash(n))
momentum_conservations_signs={r'p_\mu':+1,r'p_\nu':-1,r'p_{\bar\nu}':-1,r'p_e':-1}

T1=(1j/2)*(GammaMatrix(alpha)*GammaMatrix(beta)-GammaMatrix(beta)*GammaMatrix(alpha))
conj_T1=(1j/2)*(GammaMatrix(-alpha)*GammaMatrix(-beta)-GammaMatrix(-beta)*GammaMatrix(-alpha))
T2=(1j/2)*(GammaMatrix(rho)*GammaMatrix(sigma)-GammaMatrix(sigma)*GammaMatrix(rho))
conj_T2=(1j/2)*(GammaMatrix(-rho)*GammaMatrix(-sigma)-GammaMatrix(-sigma)*GammaMatrix(-rho))

g1=GammaMatrix(mu)
g2=GammaMatrix(nu)
conj_g1=GammaMatrix(-mu)
conj_g2=GammaMatrix(-nu)


"""
define the effective Hamitonian 	
"""
Hs=[GC(u1,PR*S.One,u2),GC(u3,S.One*PR,u4)]
Hv=[GC(u1,PL*g1,u2),GC(u3,conj_g1*PL,u4)]



Z1=Hs[0].conj()*Hv[0]
Z2=Hs[1].conj()*Hv[1]
Z3=Z1*Z2
print('ampiltude square:Y1=',Evl(Z3))
Y1=Evl(Z3)
print('\n')
print('Y1.args[3]:',Y1.args[3])

"""
in order to reduce the number of variable in angles as less as possible, We replace expansion in terms of p_mu at rest frame which has vec{p}_mu=0
so we can use tensorInnerproductSubs

Y2=tensorInnerproductSubs(momentum_conservations_signs, Y1.args[3], r'p_\nu', r'p_{\bar\nu}')

"""


Y2=Y1.args[3]
print('\nY2:',Y2)

if  isinstance(Y2.expand(),Mul):
    Y3=Add(* [replaceTensor(x,MatrixReps) for x in Y2.expand().args])
else:
    Y3=replaceTensor(Y2.expand(),MatrixReps)

print('\nY3:',Y3)


#intergrate the phase space and substitute the momentum in term of known energy and mass for later intergration
Y4=integrate(integrate(Y3,(phi_nu,0,2*pi)),(phi_e,0,2*pi)).subs(p_e,E_e).subs(p_nub,E_nub).subs(p_nu,E_nu).subs(E_nub,m_mu-E_nu-E_e)
print('\nY4:',Y4)


"""
in case there is the factor of cos(theta_nu), we can use
the momentum conservation to instead of the terms of energy and mass
"""
Y5=Y4.subs(cos(theta_nu),(m_mu**2-2*m_mu*(E_e+E_nu)+2*E_e*E_nu)/(2*E_e*E_nu))
print('\nY5:',Y5)


Y6=integrate(Y5,(E_nu,m_mu/2-E_e,m_mu/2))
Yf=factor(nsimplify(Y6.simplify()))

print('Yf:',Yf)

