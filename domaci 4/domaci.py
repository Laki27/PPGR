import numpy as np

def ParametriKamere(T):
    T0=T[:, :-1]
    if np.linalg.det(T0)<0:
        T=-T

    c1=np.linalg.det(T[:, 1:])
    c2=-np.linalg.det(np.array([T[:, 0], T[:, 2], T[:, 3]]).T)
    c3=np.linalg.det(np.array([T[:, 0], T[:, 1], T[:, 3]]).T) 
    c4=-np.linalg.det(T[:, :-1])

    c1=c1/c4
    c2=c2/c4
    c3=c3/c4

    C=np.array([c1, c2, c3])
    Q, R = np.linalg.qr(np.linalg.inv(T0))

    if R[0][0] < 0:
        R[0] = -R[0]
        Q[:,0] = -Q[:,0]
    if R[1][1] < 0:
        R[1] = -R[1]
        Q[:,1] = -Q[:,1]
    if R[2][2] < 0:
        R[2] = -R[2]
        Q[:,2] = -Q[:,2]

    A = Q.T
    K = T0.dot(Q) 

    if K[2][2] != 1:
        K = K / K[2][2]

    return K, A, C  

n=4

T = np.array([[5, -1-2*n, 3, 18-3*n],
             [0, -1, 5, 21],
             [0, -1, 0, 1]]) 

K, A, C=ParametriKamere(T)
print("Parametar K:\n", K.round(5))
print("Parametar A:\n", A.round(5))
print("Parametar C:\n", C.round(5))
print("\n")


M1 = np.array([460, 280, 250, 1])
M2 = np.array([50, 380, 350, 1])
M3 = np.array([470, 500, 100, 1])
M4 = np.array([380, 630, 50*n, 1])
M5 = np.array([30*n, 290, 0, 1])
M6 = np.array([580, 0, 130, 1])

M1p = np.array([288, 251, 1])
M2p = np.array([79, 510, 1])
M3p = np.array([470, 440, 1])
M4p = np.array([520, 590, 1])
M5p = np.array([365, 388, 1])
M6p = np.array([365, 20, 1])

originali = np.array([M1, M2, M3, M4, M5, M6])
projekcije = np.array([M1p, M2p, M3p, M4p, M5p, M6p])

def CameraDLP(originali, projekcije):
    x = originali[0][0]
    y = originali[0][1]
    z = originali[0][2]
    t = originali[0][3]

    xp = projekcije[0][0]
    yp = projekcije[0][1]
    zp = projekcije[0][2]

    A = np.array([
        [0, 0, 0, 0, -zp*x, -zp*y, -zp*z, -zp*t, yp*x, yp*y, yp*z, yp*t],
        [zp*x, zp*y, zp*z, zp*t, 0, 0, 0, 0, -xp*x, -xp*y, -xp*z, -xp*t]
    ])

    for i in range(1, len(originali)):
        x = originali[i][0]
        y = originali[i][1]
        z = originali[i][2]
        t = originali[i][3]


        xp = projekcije[i][0]
        yp = projekcije[i][1]
        zp = projekcije[i][2]

        row1 = np.array([0, 0, 0, 0, -zp*x, -zp*y, -zp*z, -zp*t, yp*x, yp*y, yp*z, yp*t])
        row2 = np.array([zp*x, zp*y, zp*z, zp*t, 0, 0, 0, 0, -xp*x, -xp*y, -xp*z, -xp*t])

        A = np.vstack((A, row1))
        A = np.vstack((A, row2))

    U,S,V = np.linalg.svd(A)

    T = V[-1].reshape(3,4)
    
    return T

T = CameraDLP(originali, projekcije)
T = T / T[0][0]
print("Matrica kamere T:\n", T)
print("\n")

T1=np.array([450, 255, 0, 1])
T2=np.array([550, 45, 0, 1])
T3=np.array([350, 290, 85, 1])
T4=np.array([245, 52, 0, 1])
T5=np.array([245, 52, 14, 1])
T6=np.array([9, 30, 25, 1])

t1p=np.array([269, 421, 1])
t2p=np.array([125, 361, 1])
t3p=np.array([337, 328, 1])
t4p=np.array([592, 390, 1])
t5p=np.array([598, 293, 1])
t6p=np.array([474, 100, 1])

originali = np.array([T1, T2, T3, T4, T5, T6])
projekcije = np.array([t1p, t2p, t3p, t4p, t5p, t6p])

T = CameraDLP(originali, projekcije)
print("Matrica moje kamere T: \n", T)

K, A, C = ParametriKamere(T)
print("Parametar K:\n", K.round(5))
print("Parametar A:\n", A.round(5))
print("Parametar C:\n", C.round(5))
print("\n")