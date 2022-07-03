from math import log, exp, sin, cos, pi
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np

def progonka(matrix):
    n = len(matrix)
    a = [-1 * matrix[0][1] / matrix[0][0]]
    b = [matrix[0][n] / matrix[0][0]]
    for i in range(1, n):
        y = matrix[i][i] + matrix[i][i - 1] * a[i - 1]
        a.append(-1 * matrix[i][i + 1] / y)
        b.append((matrix[i][n] - matrix[i][i - 1] * b[i - 1]) / y)
    b.append((matrix[n-1][n] - matrix[n-1][n-2] * b[n-2]) / (matrix[n-1][n-1] + matrix[n-1][n-2] * a[n-2]))
    x = [0 for i in range(len(matrix))]
    x[-1] = b[-1]
    for i in range(len(matrix)-2, -1, -1):
        x[i] = a[i] * x[i+1] + b[i]
    return x

def main():
    p_functons = [
        lambda x: -4,
        lambda x: -1,
        lambda x: -1
    ]
    f_functons = [
        lambda x: 0,
        lambda x: 4*sin(x),
        lambda x: 4*x*exp(x)
    ]
    an_functons = [
        lambda x: (-2/5)*cos(2*x)-(1/5)*sin(2*x),
        lambda x: sin(x)-cos(x)*(1+2*x),
        lambda x: exp(pi/2)*(cos(x)+sin(x))+2*(x-1)*exp(x)
    ]
    koshi_data = [
        (0, pi/4),
        (0, pi/2),
        (0, pi/2)
    ]
    alpha_data = [
        (1, 0),
        (1, 0),
        (1, 0)
    ]
    beta_data = [
        (1, 1),
        (1, 0),
        (1, 0)
    ]
    # wb = xl.Workbook()
    for i in range(len(koshi_data)):
        eps = 0.01
        n = int((koshi_data[i][1]-koshi_data[i][0])/eps)+1
        matrix = np.zeros((n, n+1))
        matrix[0][0] = -1-alpha_data[i][0]*eps
        matrix[0][1] = 1
        matrix[0][n] = alpha_data[i][1]*eps
        for j in range(1, n-1):
            matrix[j][j-1] = 1
            matrix[j][j] = -2-p_functons[i](koshi_data[i][0]+j*eps)*(eps**2)
            matrix[j][j+1] = 1
            matrix[j][n] = f_functons[i](koshi_data[i][0]+j*eps)*(eps**2)
        matrix[n-1][n-2] = -1
        matrix[n-1][n-1] = 1 - beta_data[i][0]*eps
        matrix[n-1][n] = beta_data[i][1]*eps
        a = progonka(matrix)
        coords = []
        graph_coords = []
        for j in range(len(a)):
            graph_coords.append((koshi_data[i][0] + j * eps, a[j]))
            coords.append((koshi_data[i][0] + j * eps, a[j]))
            coords.pop(0)

        eps /= 2
        n = int((koshi_data[i][1] - koshi_data[i][0]) / eps)+1
        matrix2 = np.zeros((n, n + 1))
        matrix2[0][0] = -1 - alpha_data[i][0] * eps
        matrix2[0][1] = 1
        matrix2[0][n] = alpha_data[i][1] * eps
        for j in range(1, n - 1):
            matrix2[j][j - 1] = 1
            matrix2[j][j] = -2 - p_functons[i](koshi_data[i][0] + j * eps) * (eps**2)
            matrix2[j][j + 1] = 1
            matrix2[j][n] = f_functons[i](koshi_data[i][0] + j * eps) * (eps**2)
        matrix2[n - 1][n - 2] = -1
        matrix2[n - 1][n - 1] = 1 - beta_data[i][0] * eps
        matrix2[n - 1][n] = beta_data[i][1] * eps
        a = progonka(matrix2)
        coords = []
        graph_coords_any = []
        for j in range(len(a)):
            graph_coords_any.append((koshi_data[i][0] + j * eps, a[j]))
            coords.append((koshi_data[i][0] + j * eps, a[j]))
            coords.pop(0)

        plt.plot(np.array([deepcopy(a[0]) for a in graph_coords]),
                 np.array([deepcopy(a[1]) for a in graph_coords]), label=f'h = {eps * 2}')
        plt.plot(np.array([deepcopy(a[0]) for a in graph_coords_any],),
                 np.array([deepcopy(a[1]) for a in graph_coords_any]), label=f'h = {eps}')
        plt.plot(np.array([deepcopy(a[0]) for a in graph_coords_any]),
                 np.array([an_functons[i](deepcopy(a[0])) for a in graph_coords_any]), label='Точное')
        plt.legend()
        plt.show()
        mn = 0
        for x in graph_coords:
            if abs(an_functons[i](x[0]) - x[1]) > mn:
                mn = abs(an_functons[i](x[0]) - x[1])
        n1 = mn
        mn = 0
        for x in graph_coords_any:
            if abs(an_functons[i](x[0]) - x[1]) > mn:
                mn = abs(an_functons[i](x[0]) - x[1])
        print(f'''Максимальная невязка на 1 решении = {n1
        }, Максимальная невязка на 2 решении = {mn}, Отношение невязок = {n1 / mn}''')

if __name__ == '__main__':
    main()