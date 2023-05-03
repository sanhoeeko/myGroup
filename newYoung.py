#此模块用来做杨图的计算，输入阶数，返回所有计算结果

import numpy as np
import copy as cp
import itertools as it
import collections as cl
import groupLinear as gl

class heap:
    def __init__(self,height):
        self.height=height
        self.eles=np.zeros((height,)).astype(int)
        self.full=False
        self.ptr=0
    def add(self,x):
        if self.full==False:
            self.eles[self.ptr]=x
            self.ptr+=1
            if self.ptr>=self.height:
                self.full=True
    def getPerms(self):
        lis=list(it.permutations(self.eles))
        #res=[Perm(i) for i in lis]
        return lis

class youngDiagram:
    def __init__(self,tup,n):
        self.n=n
        self.tup=tup
        self.heaps=[]
        for i in tup:
            self.heaps.append(heap(i))
        self.heaps[0].add(1)
        self.goal=2
    def add(self):
        res=[]
        for i in range(len(self.heaps)):
            x=cp.deepcopy(self)
            if (i==0 or self.heaps[i-1].ptr>self.heaps[i].ptr)and 1-self.heaps[i].full: #可访问的堆
                x.heaps[i].add(self.goal)
                x.goal+=1
                res.append(x)
        return res
    def getPerms(self): #将所有堆贡献的置换做幂集
        lis=[i.getPerms() for i in self.heaps]
        self.perms=list(it.product(*lis))
    def getMat(self):
        mat=np.zeros((self.n,self.n)).astype(int)
        for i in range(len(self.heaps)):
            for j in range(len(self.heaps[i].eles)):
                mat[i,j]=self.heaps[i].eles[j]
        self.mat=mat
    def __str__(self):
        return str(self.mat)
    def __repr__(self):
        return str(self)
    def findDual(self,dad):
        for diagram in dad.diagrams:
            if (diagram.mat==self.mat.T).all():
                self.T=diagram
                return
        print(self.mat)
        self.T=None
    def __eq__(self,o):
        return (self.mat==o.mat).all()

class atlas:
    def __init__(self,n):
        self.n=n
        self.tups=divide(n)
        self.diagrams=[youngDiagram(i,n) for i in self.tups]
        for i in range(self.n-1): #生长杨表
            self.add()
        for k in self.diagrams: #只有杨表生长完才能做的事
            k.getMat()
        #self.diagrams=clearRepeat_unHashable(self.diagrams)
        for k in self.diagrams:
            k.getPerms()
        for k in self.diagrams:
            k.findDual(self)
    def add(self):
        res=[]
        for i in self.diagrams:
            res.extend(i.add())
        self.diagrams=res
    def __str__(self):
        return '\n\n'.join([str(i) for i in self.diagrams])
    def __repr__(self):
        return str(self)
    def getPerms(self):
        res=[]
        for i in self.diagrams:
            i.getPerms()
            res.append(i.perms)
        #self.perms=allToTup(res)
        self.perms=res

def divide(num): #用深度优先搜索获取整数的所有分支表达
    attempt=[[1]]
    for t in range(num-1): #可以算出迭代次数，故代替递归
        mid=[]
        for arr in attempt:
            if sum(arr)<num:
                x=cp.deepcopy(arr)
                x.append(1)
                mid.append(x)
                y=cp.deepcopy(arr)
                y[-1]+=1
                mid.append(y)
        attempt=mid
    for i in attempt:
        i.sort(reverse=True)
    return clearRepeat_unHashable(attempt)

def clearRepeat_unHashable(lis):
    res=[]
    for i in lis:
        if i not in res:
            res.append(i)
    return res

def dagger(atl):
    res=cp.deepcopy(atl)
    for i in range(len(res.diagrams)):
        res.diagrams[i]=res.diagrams[i].T
    return res

'''   
class Perm:
    def __init__(self,tup):
        self.tup=tup
    def __repr__(self):
        return str(self.tup)
    def __str__(self):
        return repr(self)
    
def allToTup(x):
    if not isinstance(x,cl.Iterable):
        return x.tup
    else:
        return [allToTup(i) for i in x]
    '''

if __name__=='__main__':
    x=atlas(4)
    x.getPerms() #list<list<tuple<tuple>>>>
    print(x.perms)