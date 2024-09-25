import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

from matplotlib.font_manager import FontProperties
from qiskit import QuantumCircuit, Aer, execute

from Opertion import Operation
from QAOA_circuit import QAOA_circuit
from QAOA_black_box import QAOA_black_box
from QAOA_init import QAOA_init
from QAOA_cost import QAOA_cost
from QAOA_XY_mixer import QAOA_XY_mixer

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

K = 3
T = 4

op1 = Operation(1, 1, 1, 0, 0)  # timespan,job,machine,th
op2 = Operation(1, 2, 1, 1, 0)
op3 = Operation(1, 1, 2, 2, 0)
# op4 = Operation(2, 3, 1, 3,0)

# ops = []
ops = [op1, op2, op3]

p = 1
QAOA_black_box_obj = QAOA_black_box(K, T, ops, p)


def obj_function(gamma, beta):
    QAOA_circuit = QuantumCircuit(K * T)

    init_circuit_object = QAOA_init(K, T)
    QAOA_circuit.append(init_circuit_object.get_standerd_init_circuit(), range(K * T))
    cost_obj = QAOA_cost(K, T, gamma, ops)
    cost_circuit1 = cost_obj.get_cost_circuit1()
    cost_circuit2 = cost_obj.get_cost_circuit2()
    cost_circuit3 = cost_obj.get_cost_circuit3()
    QAOA_circuit.append(cost_circuit1, range(K * T))
    QAOA_circuit.append(cost_circuit2, range(K * T))
    QAOA_circuit.append(cost_circuit3, range(K * T))
    mixer_obj = QAOA_mixer(K, T, beta)
    mixer_circuit = mixer_obj.get_mixer_circuit()
    QAOA_circuit.append(mixer_circuit, range(K * T))

    QAOA_circuit.measure_all(K * T, K * T)
    backend = Aer.get_backend('qasm_simulator')
    counts = execute(QAOA_circuit, backend, seed_simulator=10, shots=10000).result().get_counts()

    return QAOA_black_box_obj.compute_port_energy(QAOA_black_box_obj.invert_counts(counts))


gamma = np.linspace(0, np.pi, 100)
beta = np.linspace(0, 2 * np.pi, 100)
prefix = 'data_'
suffix = '_portfolio_qaoa_p1.dat'
preprefix = ".\data\\"
filename = preprefix + prefix + suffix

# for x1 in gamma:
#      for x2 in beta:
#          z = (obj_function(x1, x2))
#          f = open(filename, 'a')
#          f.write(str(x1) + ' ' + str(x2)+ ' ' + str(z)+ "\n")
#          f.write("\n")
#          f.close()

# # 绘制三维图-获取数据
X1, X2 = np.meshgrid(gamma, beta)
data = np.loadtxt(filename)
z = data[:, 2]
# print(z)
m = int(np.sqrt(len(z)))
Z = matrix = z.reshape((m, m))

# 1. 创建三维图像
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
# 绘制三维图像
ax.plot_surface(X1, X2, Z, cmap='viridis')
# 添加轴标签
font = FontProperties(size=16)
ax.set_xlabel('gamma', fontproperties=font)
ax.set_ylabel('beta', fontproperties=font)
ax.set_zlabel('能量值', fontproperties=font)
plt.savefig('.\data\\energy_3D_2.png', format='png')
plt.show()
# 设置视角
ax.view_init(elev=30, azim=45)
fig.set_size_inches(12 / 2.54, 12 / 2.54)

# 2. 绘制填充等高线图
CS = plt.contourf(X1, X2, Z)
# CS = plt.contour(X1, X2, Z)
plt.colorbar(CS)
plt.xlabel('gamma')
plt.ylabel('beta')
plt.savefig('.\data\\energy_contour_2.png', format='png')
plt.show()
