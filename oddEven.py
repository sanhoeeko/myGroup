import groupLinear as gl
import copy as cp
import numpy as np

#判断轮换的奇偶性
def parity(tuptup,nn):
    return int(np.linalg.det(tup2mat_perm(tuptup, nn)))

def P(i,j,n):
    p=np.eye(n)
    w=cp.deepcopy(p[i])
    p[i]=cp.deepcopy(p[j])
    p[j]=w
    return p

def tup2mat_alter(tuptup,nn):
    lis=[]
    for i in tuptup:
        if len(i)>1:
            lis.append(i)
    mats=[]
    for i in lis:
        n=len(i)
        for j in range(n-1):
            mats.append(P(i[j]-1,i[j+1]-1,nn))
    if mats==[]:return np.eye(nn)
    m=mats[0]
    for i in range(1,len(mats)):
        m=mats[i]@m #左右乘勿搞反了！
    return m

def tup2mat_perm(tuptup,nn):
    lis=[]
    for i in tuptup:
        if len(i)>1:
            lis.append(i)
    mats=[]
    for i in lis:
        a=np.eye(nn)
        b=np.eye(nn)
        z=list(i)
        z.sort()
        cnt=0
        for j in i:
            b[z[cnt]-1,:]=a[j-1];cnt+=1
        mats.append(b)
    if mats==[]:return np.eye(nn)
    m=mats[0]
    #print(tuptup)
    #print(mats)
    for i in range(1,len(mats)):
        m=mats[i]@m #左右乘勿搞反了！
    return m

def link2group(tup,g:gl.group):
    mat=tup2mat_perm(tup,g.eles[0].mat.shape[0])
    for i in range(g.n):
        #print(g.eles[i].mat)
        #print(g.eles[i].perm)
        if (g.eles[i].mat==mat).all():
            return i
    return g.n+1

if __name__=='__main__':
    a=((4,3,2,1),)
    s=gl.S(4)
    print(link2group(a,s))