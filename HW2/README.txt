一、函数说明
主函数为Poisson_GMRES.cpp，其他函数均为draft

二、参数设置
h：网格尺度，h2=h*h，n=1/h，N=网格总点数；
m：单次GMRES进行m次Arnoldi迭代；
j_num：数组存储每一行的列数；
A，B，PI：存储常用常数
residual：残量

三、结果复现
直接运行Poisson_GMRES.cpp
