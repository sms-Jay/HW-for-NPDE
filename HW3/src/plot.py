import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False



# 手动输入数据
X = np.log(np.array([0.1, 0.05, 0.025,0.0125,0.006250, 0.003125]))
Y = np.log(np.array([2.55e-3, 8.15e-4, 3.15e-4, 2.80e-4, 8.90e-5, 4.42e-5])) #FDM forward
Z = np.log(np.array([2.99e-3, 1.41e-3, 7.17e-4, 4.37e-4, 1.95e-4, 9.82e-5])) #FDM backward h
A = np.log(np.array([2.80e-3, 8.76e-4, 3.27e-4, 2.81e-4, 9.00e-5])) #FDM backward h^2
W = np.log(np.array([2.09e-3, 7.71e-4, 2.54e-4, 3.48e-4, 1.02e-4, 6.04e-5])) # FVM forward
V = np.log(np.array([3.24e-3, 1.51e-3, 6.97e-4, 4.91e-4, 2.02e-4, 1.06e-4])) # FVM backward h
B = np.log(np.array([1.96e-3, 7.62e-4, 2.68e-4, 3.47e-4, 1.02e-4])) # FVM backward h^2
# 使用polyfit进行线性回归
slope, intercept = np.polyfit(X, Y, 1)

print(f"回归方程: Y = {intercept:.4f} + {slope:.4f}X")
print(f"斜率 (b): {slope:.4f}")
print(f"截距 (a): {intercept:.4f}")


Y_pred = intercept + slope * X
# 可视化
plt.scatter(X, Y, color='blue', label='data')
plt.plot(X, Y_pred, color='red', label='regression')
plt.xlabel('log h')
plt.ylabel('log e_h')
plt.legend()
plt.title(f'log e_h = {intercept:.4f} + {slope:.4f} log h')
plt.grid(True)
plt.savefig('FDM-forward-dt=0.25h^2.png', dpi=300)