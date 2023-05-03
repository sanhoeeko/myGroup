import collections as cl
import groupLinear as gl
import newYoung as ny
import oddEven as oe
import numpy as np

class tgroup:
    def __init__(self,g:gl.group):
        mat=g.getTable(numerical=True)
        n=g.n;self.n=n
        ten=np.zeros((n,n,n)).astype(int)
        for i in range(n):
            for j in range(n):
                ten[i,j,mat[i,j]]=1
        self.tensor=ten
        self.g=g
    def mul(self,a,b):
        return np.einsum('ijk,i,j->k',self.tensor,a,b)
    def vec(self,x):
        v=np.zeros((self.n,)).astype(int)
        v[x]=1
        return v
    def search(self,v):
        m=num(v)
        return self.g.eles[m]
    
def num(vec):
    for i in range(len(vec)):
        if vec[i]==1:return i
def getTGroupS(num):
    return tgroup(gl.S(num))
        
def getGroupVector(num,tg,x:ny.atlas): #list<list<vec>>
    tg.parities=np.zeros((tg.n,)).astype(int)
    def allToVec(x):
        if type(x)==tuple:
            num=oe.link2group(x,tg.g)
            tg.parities[num]=oe.parity(x,tg.g.eles[0].mat.shape[0])
            return tg.vec(num)
        else:
            return [allToVec(i) for i in x]
    x.getPerms()
    p=x.perms
    y=allToVec(p)
    return y
        
def S(lis):
    return sum(lis)
def A(lis,tg):
    for i in range(len(lis)):
        m=num(lis[i])
        lis[i]*=tg.parities[m]
    return sum(lis)

#经验证，产生了各行的对称和反对称算符

def Y(num):
    tg=tgroup(gl.S(num))
    x=ny.atlas(num)
    lin=getGroupVector(num,tg,x)
    y=ny.dagger(x)
    col=getGroupVector(num,tg,y)
    Ss=[S(i) for i in lin]
    As=[A(i,tg) for i in col]
    n=len(Ss)
    res=[]
    for i in range(n):
        res.append(tg.mul(Ss[i],As[i]))
    return res

def Schmidt(vecs):
    res=[]
    for i in vecs:
        for j in res:
            i=i-np.einsum('i,i->',i,j)*j #j已经归一化
        nn=np.linalg.norm(i)
        if abs(nn-0)<1e-6:continue
        i=i/nn
        res.append(i)
    return np.array(res)
def noRep(vecs):
    res=[]
    for i in vecs:
        msg=0
        for j in res:
            if (i==j).all():
                msg=1;break
        if msg==0:res.append(i)
    return res

#-------------------main---------------------#
def getReps(num): #->list<ten3>
    ys=Y(num)
    tg=tgroup(gl.S(num))
    res=[]
    for i in range(len(ys)):
        mat=np.einsum('ijk,k->ij',tg.tensor,ys[i])
        r=np.linalg.matrix_rank(mat)
        basis=[];cnt=0;R=0;Rsave=0
        while(True):
            basis.append(mat[cnt]);cnt+=1
            R=np.linalg.matrix_rank(np.array(basis))
            if R==r:break;
            if R==Rsave:
                basis=basis[:-1]
            else:
                Rsave=R
        basis=Schmidt(basis)
        rep=np.einsum('ijk,nj,mk->inm',tg.tensor,basis,basis) #求表示
        rep=np.around(rep,6)
        res.append(rep)
    return res

def getCharTen2d(reps,noSort=False): #->ten2
    #直接对张量求迹得所有元素的特征标χ(g)
    chars=[]
    for ten in reps: #ten->N*dim*dim
        chi=np.einsum('inn->i',ten)
        chars.append(chi)
    if noSort:return np.array(chars)
    mat=np.array(chars).T
    chars=list(mat)
    ma=np.max(mat);mi=np.min(mat)
    m=ma-mi
    chars.sort(key=lambda x:vecHash(x,m),reverse=True)
    mat=np.array(chars).T
    return mat

def getCharChart(chars):
    #去重得到特征标表
    chart=noRep(chars)
    chart=np.array(chart).T
    chart=noRep(chart)
    chart=np.array(chart).T
    return chart

def vecHash(vec,m):
    arr=[m**i for i in range(len(vec))]
    arr=np.array(arr)
    return np.dot(arr,vec)

if __name__=='__main__':
    reps=getReps(4)
    chars=getCharTen2d(reps)
    chart=getCharChart(chars)
    print(chart)