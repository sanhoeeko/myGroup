# MyGroup

把群论课上讲的一些繁琐的计算交给计算机完成。其中，

`groupLinear.py` 实现了最基本的群置换表示、群的乘法运算、求陪集、求共轭类，以及找到一些正规子群。

`groupTensor_pack.py` 实现了通过画杨表的方式计算置换群的各种表示，并据此求特征标表。

`cg.py` 在此基础上求CG系数。

### 依赖关系

`cg`
 └ `groupTensor_pack`
　├ `newYoung`
　 |　└ `groupLinear`
　├ `oddEven`
　 |　└ `groupLinear`
　└ `groupLinear`