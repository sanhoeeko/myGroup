#基于线性表示的群的实现
import numpy as np
import itertools as it
import math as mh

class element:
    def __init__(self,msg):
        if type(msg)==tuple:
            self.perm=msg #perm->置换表达
            self.n=len(self.perm)
            self.linearRepresentation()
        if type(msg)==np.ndarray:
            self.mat=msg
            self.n=self.mat.shape[0]
            self.getPerm()
        self.alt=self.toAlternation() #alt->轮换表达
        self.char=None
    def __str__(self):
        return ''.join([str(i) for i in self.alt])
    def __repr__(self):
        return str(self)
    def linearRepresentation(self): #行向量表示
        self.mat=np.zeros((self.n,self.n))
        for i in range(self.n):
            self.mat[i,self.perm[i]-1]=1
        #auto return
    def getPerm(self):
        vec=np.arange(1,self.n+1)
        self.perm=tuple(self.mat@vec)
    def inv(self):
        return element(self.mat.T)
    def __mul__(self,o):
        return element(self.mat@o.mat)
    def toAlternation(self):
        res=[]
        def reduce(perm):
            lis=[int(perm[0])]
            s=self.perm[int(lis[0]-1)]
            while s!=perm[0]:
                lis.append(int(s))
                s=self.perm[int(s-1)]
            res.append(tuple(lis))
            perm=tuple(set(perm)-set(lis))
            if perm!=():
                reduce(perm)
        reduce(self.perm)
        resres=[]
        for i in res:
            if len(i)>1:resres.append(i)
        if len(resres)==0:resres.append('e')
        return resres

class coset:
    def __init__(self,eles):
        self.eles=eles
    def leftMul(self,o):
        el=[o*x for x in self.eles]
        return coset(el)
    def rightMul(self,o):
        el=[x*o for x in self.eles]
        return coset(el)
    def reArrange(self):
        self.eles.sort(key=lambda x:binary(x.perm))
    def __eq__(self,o):
        self.reArrange();o.reArrange()
        n=len(self.eles)
        for i in range(n):
            if not self.eles[i].perm==o.eles[i].perm:
                return False
        return True
    def __str__(self):
        return '\n'.join([str(x) for x in self.eles])
    def __repr__(self):
        return str(self)

#群元的排序：按照permutation自带的字母升序排序        
class group(coset):
    def __init__(self,eles):
        self.eles=eles
        self.n=len(self.eles)
        self.classes=None
    def getTable(self,numerical=False):
        res=np.zeros((self.n,self.n)).astype(element)
        for i in range(self.n):
            for j in range(self.n):
                res[i,j]=self.eles[i]*self.eles[j]
        if numerical==True:
            names=[repr(x) for x in self.eles]
            for i in range(self.n):
                for j in range(self.n):
                    res[i,j]=names.index(repr(res[i,j]))
            res=res.astype(int)
        return res
    def outputTable(self,numerical=False):
        t=self.getTable(numerical)
        with open('output.txt','w') as w:
            for i in t:
                for j in i:
                    w.write(str(j)+';')
                w.write('\n')
        print('written.')
    def getClasses(self):
        if self.classes==None:
            for i in range(len(self.eles)):
                char,vec=np.linalg.eig(self.eles[i].mat)
                char=list(char)
                char.sort(key=posiArg)
                self.eles[i].char=np.array(char)
            res=[]
            for i in self.eles:
                msg=0
                for j in res:
                    if vecEq(i.char,j.char):
                        msg=1
                        j.eles.append(i);break
                if msg==0:
                    res.append(conjClass([i],i.char))
            self.classes=res
    def ifInvariantOf(self,dad):
        for i in dad.eles:
            a=self.leftMul(i)
            b=self.rightMul(i)
            if not a==b:
                return False
        return True
        
class conjClass:
    def __init__(self,eles,char):
        self.eles=eles
        self.char=char
        self.charRepr=self.getCharRepr()
    def getCharRepr(self):
        res=[]
        for i in self.char:
            res.append(cplxRepr(i))
        return res
    def __str__(self):
        return '['+','.join(self.charRepr)+']&'+str(len(self.eles))
    def __repr__(self):
        return str(self)
        
def cplxRepr(x):
    ag=mh.atan2(x.imag,x.real)
    arg=ag/mh.pi
    for i in range(1,100):
        argi=arg*i
        if abs(argi-round(argi))<1e-6:
            rargi=round(argi)%(2*i)
            if rargi==0:return '0'
            if i==1:return 'pi'
            return '('+str(rargi)+'/'+str(i)+')pi'
def posiArg(x):
    ag=mh.atan2(x.imag,x.real)
    if ag<0:return ag+mh.pi*2
    return ag
def vecEq(vec1,vec2):
    x=abs(vec1-vec2)
    y=np.ones((vec1.shape[0]))*1e-6
    if (x<y).all():
        return True
    return False
        
def S(num):
    perms=list(it.permutations(np.arange(1,num+1)))
    eles=[]
    for i in perms:
        eles.append(element(i))
    return group(eles)

def generate(gens): #输入几个生成元，输出群
    n=len(gens[0].perm)
    e=element(np.eye(n))
    edge=[e]
    res=[e]
    def recursion(edge):
        if edge==[]:return
        newEdge=[]
        for i in gens:
            for j in edge:
                y=[i*j,j*i,i.inv()*j,j.inv()*i]
                for x in y:
                    lis=[u.mat for u in res]
                    if not In(x.mat,lis):
                        res.append(x)
                        newEdge.append(x)
        recursion(newEdge)
    recursion(edge)
    g=group(res)
    g.reArrange()
    return g
def In(mat,matlist):
    for i in matlist:
        if (mat==i).all():
            return True
    return False
def binary(tup):
    n=len(tup)
    s=0
    for i in range(n):
        s+=n**(n-i-1)*tup[i]
    return s

def findSubGroups(g):
    g.getClasses()
    res=[]
    #获取全组合，叫号抽取共轭类
    n=len(g.classes)
    call=[]
    for k in range(1,n-1):
        call.extend(list(it.combinations(range(n),k)))
    for i in call:
        gens=[]
        for j in i:
            gens.append(g.classes[j].eles[0])
        u=generate(gens)
        res.append(u)
    return res

#-------------------main---------------------#
if __name__=='__main__':
    g=S(4)
    '''
    res=findSubGroups(g)
    res=filter(lambda x:len(x.eles)==4,res)
    for i in res:
        print(i.getTable(True))
        print(i.ifInvariantOf(g))
        print()
        '''
    a=element((2,1,4,3))
    b=element((4,3,2,1))
    e=element((1,2,3,4))
    c=element((3,4,1,2))
    s=group([e,a,b,c])
    print(s.ifInvariantOf(g))