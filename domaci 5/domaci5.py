import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

L1=[684, 620, 1]
L2=[678, 499, 1]
L3=[516, 565, 1]
L4=[622, 506, 1]
L5=[972, 537, 1]
L6=[1072, 337, 1]
L7=[619, 221, 1] #uz pomoc prvog domaceg
L8=[499, 433, 1]
L9=[994, 401, 1]
L10=[1107, 186, 1]
L11=[613, 95, 1]
L12=[464, 289, 1]
L13=[816, 309, 1]
L14=[944, 285, 1]
L15=[898, 214, 1] #uz pomoc prvog domaceg
L16=[775, 237, 1]
L17=[819, 253, 1]
L18=[952, 227, 1]
L19=[905, 155, 1]
L20=[774, 181, 1]

D1=[454, 552, 1]
D2=[583, 473, 1]
D3=[408, 452, 1]
D4=[483, 415, 1]
D5=[745, 622, 1]
D6=[1055, 535, 1]
D7=[823, 293, 1] #uz pomoc prvog domaceg
D8=[538, 364, 1]
D9=[761, 435, 1]
D10=[1105, 342, 1]
D11=[849, 101, 1]
D12=[532, 169, 1]
D13=[740, 304, 1]
D14=[859, 337, 1]
D15=[912, 277, 1]
D16=[792, 243, 1] #uz pomoc prvog domaceg
D17=[748, 231, 1]
D18=[872, 270, 1]
D19=[923, 201, 1]
D20=[803, 171, 1]

L = [L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11, L12, L13, L14, L15, L16, L17, L18, L19, L20]
D = [D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20]

LL = np.array([L1, L2, L3, L4, L7, L8, L9, L11])
DD = np.array([D1, D2, D3, D4, D7, D8, D9, D11])

def oneEquation(l, r):

    x1 = l[0]
    y1 = l[1]
    z1 = l[2]

    x2 = r[0]
    y2 = r[1]
    z2 = l[2]

    return np.array([x1*x2, y1*x2, z1*x2, 
                     x1*y2, y1*y2, z1*y2, 
                     x1*z2, y1*z2, z1*z2])

# Matrica formata 8x9 - 8 jednacina dobijenih iz korespodencija
equations8 = []
for i in range(8): 
    equations8.append(oneEquation(LL[i], DD[i]))

# print("Jednacine:")
# for i in range(8):
    # print(equations[i])

# SVD dekompozicija prethodne matrice
U, S, V = np.linalg.svd(equations8)

# Fundamentalna matrica FF
FF = np.zeros((3, 3))
V = V[:][-1]

for i in range(3):
    for j in range(3):
        FF[i][j] = V[3*i+j]

print("1. Fundamentalna matrica:")
print(FF)

# Provera da li je R.T * FF * L = 0
C = np.ones((8, 1))
for i in range(8):
    C[i][0] = round(DD[i].T.dot(FF).dot(LL[i]), 5)

# print(C.T)

print(np.linalg.det(FF)) # treba da bude blizu 0

"""***Trazimo epipolove***"""

U, S, V = np.linalg.svd(FF)
# print("U:\n", U)
# print("S:\n", S)
# print("V:\n", V.T)

print("e1:")
e1 = V[:][-1] # treca vrsta matrice V 
print(e1)
e1 = (1/e1[2]) * e1 # afine koordinate
print(e1)

print("\ne2:")
e2 = U.T[:][-1] # treca kolona matrice U (tj. treca vrsta od U.T)
print(e2)
e2 = (1/e2[2]) * e2 # afine koordinate
print(e2)

"""***Postizanje uslova det(FF) = 0***"""

S = np.array([[S[0], 0, 0],
              [0, S[1], 0],
              [0, 0, 0]])
print(S)

FF1 = U.dot(S).dot(V)
print(FF1)

print(np.linalg.det(FF))
print(np.linalg.det(FF1))
# det(FF1) < det(FF) pa cemo koristiti FF1



"""***Trijangulacija***"""

# Kanonska matrica kamere
T1 = np.array([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1, 0]])
print("T1:")
print(T1)

# Matrica vektorskog mnozenja
E2 = np.array([[0, -e2[2], e2[1]],
               [e2[2], 0, -e2[0]],
               [-e2[1], e2[0], 0]])
print("E2:")
print(E2)

# Matrica druge kamere
tmp = E2.dot(FF1)
T2 = np.array([[tmp[0][0], tmp[0][1], tmp[0][2], e2[0]],
               [tmp[1][0], tmp[1][1], tmp[1][2], e2[1]],
               [tmp[2][0], tmp[2][1], tmp[2][2], e2[2]]])
print("T2:")
print(T2)

# Za svaku tacku dobijamo sistem od 4 jednacine sa 4 homogene nepoznate
def equations(l, d):
    return np.array([l[1]*T1[2] - l[2]*T1[1],
                    -l[0]*T1[2] + l[2]*T1[0],
                    d[1]*T2[2] - d[2]*T2[1],
                    -d[0]*T2[2] + d[2]*T2[0]])

# Vracamo 3D koordinate
def TriD(l, r):
    U, S, V = np.linalg.svd(equations(l, r))
    P = V[-1]
    P = P / P[3]
    return P[:-1]

reconstructed = []
for i in range(len(L)):
    reconstructed.append(TriD(L[i], D[i]))
print("--------------------")

# Mnozimo z koordinatu sa nekoliko stotina
tmp = np.eye(3)
tmp[2][2] = 400
reconstructed_400 = np.zeros((len(L),3))
for i in range(len(L)):
    reconstructed_400[i] = tmp.dot(reconstructed[i])

print(reconstructed_400)
X = reconstructed_400

# ---------------------- ISCRTAVANJE -----------------------------

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim([300, 600])
ax.set_ylim([100, 400])
ax.set_zlim([100, 400])

Z = X[:4]
verts = [[Z[0],Z[1],Z[2]],
         [Z[0],Z[1],Z[3]], 
         [Z[0],Z[2],Z[3]], 
         [Z[1],Z[2],Z[3]] ]
ax.add_collection3d(Poly3DCollection(verts, facecolors='blue', linewidths=1, edgecolors='r', alpha=.25))

Z = X[4:12]
verts = [[Z[0],Z[1],Z[2],Z[3]],
         [Z[4],Z[5],Z[6],Z[7]], 
         [Z[0],Z[1],Z[5],Z[4]], 
         [Z[2],Z[3],Z[7],Z[6]], 
         [Z[1],Z[2],Z[6],Z[5]],
         [Z[4],Z[7],Z[3],Z[0]]]
ax.add_collection3d(Poly3DCollection(verts, facecolors='green', linewidths=1, edgecolors='r', alpha=.25))

Z = X[12:20]
verts = [[Z[0],Z[1],Z[2],Z[3]],
         [Z[4],Z[5],Z[6],Z[7]], 
         [Z[0],Z[1],Z[5],Z[4]], 
         [Z[2],Z[3],Z[7],Z[6]], 
         [Z[1],Z[2],Z[6],Z[5]],
         [Z[4],Z[7],Z[3],Z[0]]]
ax.add_collection3d(Poly3DCollection(verts, facecolors='yellow', linewidths=1, edgecolors='r', alpha=.25))


ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
    
plt.gca().invert_yaxis()
plt.show()