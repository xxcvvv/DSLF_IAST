# 用于计算基于DSLF的IAST曲线及其他各项参数,包括Qst、吸附选择性
## Fit
DSLF中的六个参数在实际拟合过程中易出现数值异常问题，利用遗传算法可以在无需初值的情况下得到理想的参数，这里使用Matlab中的GA工具箱，通过提供上下限得到较理想的结果。推荐初试时设置q为0到20，b为0到0.01，n为0到2.
## Qst
使用Fortran进行计算，详细步骤见注释
## Selectivity
使用Fortran进行计算，详细步骤见注释
参考自'http://muchong.com/t-5109332-1'