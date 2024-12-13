from sympy import Symbol, S,symbols,Mul,Function,conjugate,Add,diag,I,integrate,cos,sin,diag,pi,nsimplify,factor
from sympy.tensor.tensor import TensorHead,Tensor,TensorIndexType, tensor_indices, tensor_heads,TensorIndex,TensorSymmetry
from sympy.core.expr import Expr

Lorentz = TensorIndexType('Lorentz', dummy_name='L', dim=4)
Lg = Lorentz.metric

def g(id1,id2):
    m0, m1 = tensor_indices('m0,m1', Lorentz)
    p,q=1,1;
    
    if not isinstance(id1,TensorIndex):
        p = tensor_heads(id1.name, [Lorentz])(-m0)
        id1=m0
        
    if not isinstance(id2,TensorIndex):
        q = tensor_heads(id2.name, [Lorentz])(-m1)
        id2=m1
        
    return (Lg(id1,id2)*p*q)
    
class Spinor(Symbol):
    def __new__(cls,momentum,mass,polarization,**kwargs):
        obj = super().__new__(cls,momentum.name,**kwargs,commutative=False)
        obj.momentum = momentum
        obj.mass = mass
        obj.polarization=polarization
        return obj
    def __init__(self,momentum,mass,polarization):
        self.name =f"u({momentum})" 



class GammaSlash(Symbol):  
    def __new__(cls,name,*args,**kwargs):
        obj_name=name,
        obj = super().__new__(cls,name.name,*args,**kwargs,commutative=False)
        obj.index = name
        return obj
    def __init__(self,name):
        self.name =f" gs({self.index.name})"
        self.is_gamma5=False
        self.is_gs=True
         
class GammaMatrix(Symbol):
    def __new__(cls,name,*args,**kwargs):
        obj_name=f"gamma^{name.name}" if name.is_up else f"gamma_{name.name}"
        obj = super().__new__(cls,obj_name,*args,**kwargs,commutative=False)
        obj.index = name
        return obj
        
    def __init__(self,name):
        self.is_gamma5=False
        self.is_gs=False
            
class GammaMatrix5(Symbol):
    def __new__(cls,*args,**kwargs):
        obj = super().__new__(cls,name='gamma^5',*args,**kwargs,commutative=False)
        return obj
    
    def __init__(self):
        self.name =f"gamma^5"
        self.is_gamma5=True
        self.is_gs=False
        
    def __pow__(self, power):
        if power == 2:
            return S.One
        elif power == 3:
            return self  # G5^3 = -G5
        else:
            raise ValueError(f"Power {power} not implemented for {self.name}")
    def __repr__(self):
        return self.name

    def __mul__(self, other):
        
        if other.__class__.__name__=='GammaMatrix5':
            return S.One
        return super().__mul__(other)

class GC:         
    def __init__(self,left_spinor,matrix,right_spinor):
        self.left_spinor = left_spinor
        self.matrix = matrix
        self.right_spinor = right_spinor
        
    def __repr__(self):
        return f"GC({self.left_spinor}{self.matrix}{self.right_spinor})"
    
    def __mul__(self, other):
        L_contract=self.right_spinor.momentum==other.left_spinor.momentum;
        R_contract=self.left_spinor.momentum==other.right_spinor.momentum;
        if L_contract and (not R_contract):
            p1=(GammaSlash(self.right_spinor.momentum)+self.right_spinor.mass)*(1+self.right_spinor.polarization*GammaMatrix5())
            return GC(self.left_spinor,self.matrix*p1*other.matrix,other.right_spinor)
        if (not L_contract) and  R_contract:
            p1=(GammaSlash(self.left_spinor.momentum)+self.left_spinor.mass)
            return GC(other.left_spinor,other.matrix*p1*self.matrix,self.right_spinor)*(1+self.right_spinor.polarization*GammaMatrix5())
        if L_contract and R_contract:
            p1=(GammaSlash(other.right_spinor.momentum)+S.One*other.right_spinor.mass)*(1+other.right_spinor.polarization*GammaMatrix5())
            p2=(GammaSlash(self.right_spinor.momentum)+S.One*self.right_spinor.mass)*(1+self.right_spinor.polarization*GammaMatrix5())
            return trace(other.matrix*p1*self.matrix*p2)
            
        return Mul(self, other, evaluate=False)
    def conj(self):
        m=(self.matrix) if self.matrix.is_symbol else Mul(* reversed(list(self.matrix.args)))
        return GC(self.right_spinor,m.subs(GammaMatrix5(),-GammaMatrix5()),self.left_spinor)
    
def reorder_with_eps(expr):    
    terms=expr.args
    reordered_terms = []
    if len(terms)==0 and not isinstance(expr,Symbol):
        return expr

    for index, term in enumerate(terms):
        if term.__class__.__name__=='LeviCivitaT':
            reordered_terms.insert(0, term)
        else:
            reordered_terms.append(term)
    # Combine terms back into an expression
    reordered_expr=1
    for x in reordered_terms:
        reordered_expr*=x
    return reordered_expr
    
def reorder_with_sign(expr):   
    
    cof=Mul(*list(filter(lambda x:not hasattr(x,'is_gs'),  expr.args)))
    terms=list(filter(lambda x:hasattr(x,'is_gs'),  expr.args))
    if len(terms)==0 and not isinstance(expr,Symbol):
        return expr
    reordered_terms = []
    # Iterate through terms to move `target` to the front
    n=0;
    for index, term in enumerate(terms):
        if term.__class__.__name__=='GammaMatrix5':
            # If target is found, prepend it and adjust the sign
            reordered_terms.insert(0, term)
            
            cof *= (-1)**(index-n)  # Flip the sign for each swap
            n+=1;
        else:
            reordered_terms.append(term)
    
    # Combine terms back into an expression
    reordered_expr=1
    for x in reordered_terms:
        reordered_expr *=x
    return cof * reordered_expr
    
def GMSimplify(GMs):
    expandsion=(GMs*S.One).expand()
    if not isinstance(expandsion,Mul):
        Z=Add(*[reorder_with_sign(gms_arg) for gms_arg in list(expandsion.args)])
        return Z
    else :
        return reorder_with_sign(expandsion)

def Trace_Sigle(expandsion):
    n= sum([1 for x in list(expandsion.args) if hasattr(x,'index') ])
    m= sum([1 for x in list(expandsion.args) if hasattr(x,'is_gamma5') and x.is_gamma5])
    coeff=(Mul(*filter(lambda x: not hasattr(x,'index')  and not (hasattr(x,'is_gamma5') and x.is_gamma5), expandsion.args) ))
    if(m==0):
        if n==2:
            indies=list(map(lambda x:x.index,filter(lambda x:hasattr(x,'index'),list(expandsion.args))))
            
            return 4*coeff*g(indies[0],indies[1])
        if n==4:
            indies=list(map(lambda x:x.index,filter(lambda x:hasattr(x,'index'),list(expandsion.args))))
            return 4*coeff*(g(indies[0],indies[1])*g(indies[2],indies[3])+g(indies[0],indies[3])*g(indies[1],indies[2])-g(indies[0],indies[2])*g(indies[1],indies[3]))
    else:
        if n==4:
            indies=list(map(lambda x:x.index,filter(lambda x:hasattr(x,'index'),list(expandsion.args))))
            #print(LeviCivitaT(1,self_indices=(indies[0],indies[1],indies[2],indies[3])))
            m0,m1,m2, m3= tensor_indices('m0,m1,m2,m3', Lorentz)
            p0,p1,p2,p3=1,1,1,1;
            id0,id1,id2,id3=indies
            
            if not isinstance(id0,TensorIndex):
               
                p0 = tensor_heads(id0.name, [Lorentz])(-m0)
                id0=m0
            if not isinstance(id1,TensorIndex):
                p1 = tensor_heads(id1.name, [Lorentz])(-m1)
                id1=m1
            if not isinstance(id2,TensorIndex):
                p2 = tensor_heads(id2.name, [Lorentz])(-m2)
                id2=m2
            if not isinstance(id3,TensorIndex):
                p3 = tensor_heads(id3.name, [Lorentz])(-m3)
                id3=m3
            
            return 4j*coeff*LeviCivitaT(1,self_indices=(id0,id1,id2,id3))*p0*p1*p2*p3
        
        if n==6:
            indies=list(map(lambda x:x.index,filter(lambda x:hasattr(x,'index'),list(expandsion.args))))
            m0,m1,m2,m3,m4,m5= tensor_indices('m0,m1,m2,m3,m4,m5', Lorentz)
            p0,p1,p2,p3,p4,p5=1,1,1,1,1,1
            id0,id1,id2,id3,id4,id5=indies
            
            if not isinstance(id0,TensorIndex):
                p0 = tensor_heads(id0.name, [Lorentz])(-m0)
                id1=m0
            if not isinstance(id1,TensorIndex):
                p1 = tensor_heads(id1.name, [Lorentz])(-m1)
                id1=m1
            if not isinstance(id2,TensorIndex):
                p2 = tensor_heads(id2.name, [Lorentz])(-m2)
                id2=m2
            if not isinstance(id3,TensorIndex):
                p3 = tensor_heads(id3.name, [Lorentz])(-m3)
                id3=m3
            if not isinstance(id4,TensorIndex):
                p4 = tensor_heads(id4.name, [Lorentz])(-m4)
                id4=m4
            if not isinstance(id5,TensorIndex):
                p5 = tensor_heads(id5.name, [Lorentz])(-m5)
                id5=m5
                
            S=g(id0,id1)*LeviCivitaT(1,self_indices=(id2,id3,id4,id5))
            S+=-g(id0,id2)*LeviCivitaT(1,self_indices=(id1,id3,id4,id5))
            S+=g(id1,id2)*LeviCivitaT(1,self_indices=(id0,id3,id4,id5))
            
            S+=g(id3,id4)*LeviCivitaT(1,self_indices=(id0,id1,id2,id5))
            S+=-g(id3,id5)*LeviCivitaT(1,self_indices=(id0,id1,id2,id4))
            S+=g(id4,id5)*LeviCivitaT(1,self_indices=(id0,id1,id2,id3))
            return -4j*coeff*S*p0*p1*p2*p3*p4*p5
    return 0

def trace(GMs):
    expandsion=GMSimplify(GMs)
    if not isinstance(expandsion,Mul):
       return Add(*list(map(lambda x:Trace_Sigle(x),list(expandsion.args))))
    else :
        return Trace_Sigle(expandsion)

def delta(id1,id2):
    m0 = tensor_indices('m0', Lorentz)
    return (Lg(id1,m0)*Lg(id2,-m0)).contract_metric(Lg)
    
class LeviCivitaT(Tensor):
    def __init__(self,name,self_indices,**is_canon_bp):
        self.name=1;
    def __new__(cls,name,self_indices,**is_canon_bp):
        instance = super().__new__(cls,tensor_head=TensorHead('epsilon', [Lorentz]*4,TensorSymmetry.fully_symmetric(-4)),indices=self_indices)
        return instance
    
    def __mul__(self, other):
        
        if other.__class__.__name__=='LeviCivitaT':
            
            self_indices=list(map(lambda x:-x,self.get_indices()))
            other_indices=other.get_indices()
            combined_indices =  self_indices+ other_indices
            self_free_indices = [idx for idx in self_indices if combined_indices.count(idx) == 1]
            other_free_indices = [idx for idx in other_indices if combined_indices.count(idx) == 1]
            self_dummy_indices = [idx for idx in self_indices if combined_indices.count(idx) == 2]
            other_dummy_indices = [idx for idx in other_indices if combined_indices.count(idx) == 2]

            if(len(self_free_indices)==2 and len(other_free_indices)==2):
                
                osfids=list(map(lambda x:-x,self_free_indices))
                f1=(-1)**(self_indices.index(self_dummy_indices[0])-0)*(-1)**(self_indices.index(self_dummy_indices[1])-1)
                f2=(-1)**(other_indices.index(other_dummy_indices[0])-0)*(-1)**(other_indices.index(other_dummy_indices[1])-1)
                dd1=delta(osfids[0],other_free_indices[0])*delta(osfids[1],other_free_indices[1])
                dd2=delta(osfids[0],other_free_indices[1])*delta(osfids[1],other_free_indices[0])
                return f1*f2*(-2)*(dd1-dd2)
            
            if(len(self_free_indices)==1 and len(other_free_indices)==1):
                osfids=list(map(lambda x:-x,self_free_indices))
                f1=(-1)**(self_indices.index(self_dummy_indices[0])-0)*(-1)**(self_indices.index(self_dummy_indices[1])-1)
                f2=(-1)**(other_indices.index(other_dummy_indices[0])-0)*(-1)**(other_indices.index(other_dummy_indices[1])-1)
                return f1*f2*(-6)*delta(osfids[0],other_free_indices[0])
        return super().__mul__(other)#.canon_bp()

def Evl(epx):
  if epx==0:
        return 0
    return Add(*[reorder_with_eps(x) for x in epx.expand().args]).canon_bp().contract_metric(Lg)    

def replaceTensor(term,matrix_reps):
    tensors=list(filter(lambda x:isinstance(x,Tensor),term.args))
    my_dict={}
    my_dict[Lorentz] = diag(1, -1, -1, -1)
    for x in tensors:
        my_dict.setdefault(x,matrix_reps[x.head.name])
    return term.replace_with_arrays(my_dict)

    
def tensor_subs(momentum_conservations_signs,exper,in1_a):
    momentum_conservations=lambda i: Add(*[ sign*tensor_heads(key, [Lorentz])(i) for key, sign in momentum_conservations_signs.items()])
    #i0 = tensor_indices('i0', Lorentz)
    sing=1
    for x in exper.args:
        if isinstance(x,Tensor) and x.head.name==in1_a:
            i0=x.indices[0]
            sing=momentum_conservations_signs[in1_a]*-1
            term=momentum_conservations(i0).subs(tensor_heads(in1_a,[Lorentz])(i0),0)
            exper=exper.subs(x,term)
    
    return exper*sing
    
def tensorInnerproductSubs(momentum_conservations_signs,exper,in1_a,in1_b):
    momentum_conservations=lambda i: Add(*[ sign*tensor_heads(key, [Lorentz])(i) for key, sign in momentum_conservations_signs.items()])
    i0 = tensor_indices('i0', Lorentz)
    sing=1
    term=(1/2)*momentum_conservations(i0).subs(tensor_heads(in1_a,[Lorentz])(i0),0).subs(tensor_heads(in1_b,[Lorentz])(i0),0)*momentum_conservations(-i0).subs(tensor_heads(in1_a,[Lorentz])(-i0),0).subs(tensor_heads(in1_b,[Lorentz])(-i0),0)
    term+=-tensor_heads(in1_a,[Lorentz])(i0)*tensor_heads(in1_a,[Lorentz])(-i0)/2
    term+=-tensor_heads(in1_b,[Lorentz])(i0)*tensor_heads(in1_b,[Lorentz])(-i0)/2
    
    sing*=momentum_conservations_signs[in1_a]
    sing*=momentum_conservations_signs[in1_b]
    
    for x in exper.args:
        if isinstance(x,Tensor) and (x.head.name==in1_a or x.head.name==in1_b):
            exper=exper.subs(x,1)  
            
    return exper*sing*term



