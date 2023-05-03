import groupTensor_pack as gt
import numpy as np
import math as mh

def directProd(ten1,ten2):
    n=ten1.shape[0]
    res=[]
    for i in range(n):
        res.append(np.kron(ten1[i,:,:],ten2[i,:,:]))
    return np.array(res)
def prod4d(a,b):
    return np.einsum('aij,akl->aijkl',a,b)

def getCG(ten1,ten2,target):
    d_lam=target.shape[1]
    n=target.shape[0]
    I=ten1.shape[1]
    J=ten2.shape[1]
    tencg=np.zeros((I,J,d_lam))
    for i in range(I):
        for j in range(J):
            for k in range(d_lam):
                s=np.sum(ten1[:,i,i]*ten2[:,j,j]*np.conj(target[:,k,k]))
                tencg[i,j,k]=np.sqrt(abs(s*d_lam/n)) #这样开方相位未必正确
    return np.einsum('ijk->kij',tencg)
    return tencg

def _getFold(prod,target):
    n=prod.shape[0]
    chi12=np.einsum('kii->k',prod)
    chi3=np.einsum('kii->k',target)
    s=np.sum(chi12*np.conj(chi3))/n
    return s
def getFolds(prod,reps):
    res=[_getFold(prod, i) for i in reps]
    return np.array(res)

reps=gt.getReps(3)
a=reps[0]
b=reps[2]
t=getCG(b,b,b)
'''
c=directProd(directProd(b,b),directProd(b,b))
t=getFolds(c,reps)
'''
B=prod4d(b, b)
t=np.array([[[mh.sqrt(1/2),0],[0,-mh.sqrt(1/2)]],[[0,-mh.sqrt(1/2)],[-mh.sqrt(1/2),0]]])
res=np.einsum('amnrl,mri,nlj->aij',B,t,t)
print(res)