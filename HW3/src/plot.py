import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False



# 手动输入数据
X = np.log(np.array([0.25, 0.125, 0.0625, 0.03125]))
Y = np.log(np.array([0.0469123,0.0207536 ,0.0106329,0.00782855 ])) 
Z = np.log(np.array([0.216301, 0.0979067 , 0.0375688,0.0217699 ,0.00706033,0.000533549 ])) 

# 使用polyfit进行线性回归（1表示1次多项式，即直线）
slope, intercept = np.polyfit(X, Y, 1)

print(f"回归方程: Y = {intercept:.4f} + {slope:.4f}X")
print(f"斜率 (b): {slope:.4f}")
print(f"截距 (a): {intercept:.4f}")

# 计算R²
Y_pred = intercept + slope * X
ss_res = np.sum((Y - Y_pred) ** 2)
ss_tot = np.sum((Y - np.mean(Y)) ** 2)
r_squared = 1 - (ss_res / ss_tot)
print(f"R² 决定系数: {r_squared:.4f}")

# 可视化
plt.scatter(X, Y, color='blue', label='实际数据')
plt.plot(X, Y_pred, color='red', label='回归线')
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.title('线性回归结果 - 使用numpy.polyfit')
plt.grid(True)
plt.show()