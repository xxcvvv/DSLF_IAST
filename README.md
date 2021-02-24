# 用于计算基于DSLF的IAST曲线及其他各项参数,包括吸附热、吸附选择性
## 所用公式
DSLF拟合公式：
![image](https://github.com/xxcvvv/DSLF_IAST/blob/main/images/DSLF.jpg)
Qst计算公式：
![image](https://github.com/xxcvvv/DSLF_IAST/blob/main/images/Qst.jpg)
## Fit_GA
DSLF中的六个参数在实际拟合过程中易出现数值异常问题，影响之后的计算，而利用遗传算法可以在无需初值的情况下得到理想的参数，这里使用Matlab中的GA工具箱，通过提供上下限得到较理想的结果。
## Qst
使用Fortran进行计算，详细步骤见注释
## Selectivity
使用Fortran进行计算，详细步骤见注释
参考自'http://muchong.com/t-5109332-1'
## Demo
所有数据均为matlab随机生成，仅供演示参考