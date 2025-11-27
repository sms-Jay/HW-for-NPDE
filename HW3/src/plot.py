import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False



# 手动输入数据
X = np.log(np.array([0.1, 0.05, 0.025]))
Y = np.log(np.array([2.55e-3, 8.22e-4, 3.10e-4, 4.36e-5])) #FDM forward
Z = np.log(np.array([2.80e-3, 8.82e-4, 3.22e-4, 4.70e-5])) #FDM backward
W = np.log(np.array([3.24e-3, 1.51e-3, 7.37e-4])) # FVM backward
V = np.log(np.array([2.09e-3, 7.71e-4, 3.17e-4]))
# 使用polyfit进行线性回归
slope, intercept = np.polyfit(X, W, 1)

print(f"回归方程: Y = {intercept:.4f} + {slope:.4f}X")
print(f"斜率 (b): {slope:.4f}")
print(f"截距 (a): {intercept:.4f}")


Y_pred = intercept + slope * X
# 可视化
plt.scatter(X, W, color='blue', label='data')
plt.plot(X, Y_pred, color='red', label='regression')
plt.xlabel('log h')
plt.ylabel('log e_h')
plt.legend()
plt.title(f'log e_h = {intercept:.4f} + {slope:.4f} log h')
plt.grid(True)
plt.savefig('FVM-backward-Euler.png', dpi=300)