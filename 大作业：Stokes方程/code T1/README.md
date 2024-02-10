# 大作业第一题
## 各函数用途
- discrete：计算网格点上真解。输入N为网格规模，输出f，g，u_true,v_true,p_true分别为网格点上的f、g、u，v，p的真实值。
- Au、Av：计算系数矩阵与乘法乘法。
- ru、rv、rp：计算残量
- DGS：DGS磨光子
- restrict：限制算子
- lift：提升算子
- vcycle：实现一次V-cycle多重网格方法

## 复现上机报告
运行FW1.m即可。改变nu1、nu2、bottom可得到不同参数下的数值结果。
  