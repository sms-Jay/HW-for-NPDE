一、函数说明
主文件为Poisson_GMRES.cpp，实现Poisson方程的GS-GMRES(m)求解器.子函数说明：
1.f,alpha,beta,g,real_u：为方程相应函数
2.vec_add,vec_subtract,num_vec,vec_dot：实现向量加、减、数乘、点乘
3.LS：用基于Givens旋转的QR分解算法求解最小二乘问题argmin \|Hy-d\|_2
4.Arnoldi：Arnoldi迭代，生成V与y
5.RHS：生成方程右端项b
6.GS：单次Gauss-Seidel迭代，用于找初值
7.Ax：矩阵-向量乘法A*x
8.Vy：矩阵-向量乘法V*y
9.GMRES：广义极小剩余法的主迭代

二、参数设置
1.h：网格尺度，h2=h*h，n=1/h，N=网格总点数；
2.m：单次GMRES至多进行m次Arnoldi迭代；
3.j_num：数组存储每一行的列数；V：存储矩阵V
4.A，B，PI：存储常用常数
5.residual：残量
6.epsilon：GMRES停机标准为residual<epsilon
7.max_iter：最大迭代次数

三、结果复现
设置相应的alpha beta，直接运行Poisson_GMRES.cpp
