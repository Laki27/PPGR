import numpy as np
import math
import cv2
original_tacke=[[-3, -1, 1],
                [3, -1, 1], 
                [1, 1, 1], 
                [-1, 1, 1]]

slika_tacke=[[-2, -1, 1], 
             [2, -1, 1], 
             [2, 1, 1], 
             [-2, 1, 1]]

def pronadji_matricu(x):
    matrica=np.array([
        [x[0][0], x[1][0], x[2][0]],
        [x[0][1], x[1][1], x[2][1]],
        [x[0][2], x[1][2], x[2][2]]
    ])

    D = np.array([x[3][0], x[3][1], x[3][2]])

    result = np.linalg.solve(matrica, D)

    lambda1=result[0]
    lambda2=result[1]
    lambda3=result[2]

    kolona1=np.array([lambda1*x[0][0], lambda1*x[0][1], lambda1*x[0][2]])
    kolona2=np.array([lambda2*x[1][0], lambda2*x[1][1], lambda2*x[1][2]])
    kolona3=np.array([lambda3*x[2][0], lambda3*x[2][1], lambda3*x[2][2]])

    P = np.column_stack([kolona1, kolona2, kolona3])

    return P

def naivni_algoritam(original, slika):
    P1=pronadji_matricu(original)
    P2=pronadji_matricu(slika)

    P=np.dot(P2, np.linalg.inv(P1))
    return P.round(5)

print("Zadati test primer")
print("Naivni algoritam: ")
matrica_naivni=naivni_algoritam(original_tacke, slika_tacke)
print(matrica_naivni) 

def dlt(original, product):

    x = original[0][0]
    y = original[0][1]
    z = original[0][2]

    xp = product[0][0]
    yp = product[0][1]
    zp = product[0][2]

    A = np.array([
        [   0,    0,    0, -zp*x, -zp*y, -zp*z,  yp*x,  yp*y,  yp*z],
        [zp*x, zp*y, zp*z,     0,     0,     0, -xp*x, -xp*y, -xp*z]
    ])

    for i in range(1, len(original)):
        x = original[i][0]
        y = original[i][1]
        z = original[i][2]

        xp = product[i][0]
        yp = product[i][1]
        zp = product[i][2]

        row1 = np.array([0, 0, 0, -zp*x, -zp*y, -zp*z, yp*x, yp*y, yp*z])
        row2 = np.array([zp*x, zp*y, zp*z, 0, 0, 0, -xp*x, -xp*y, -xp*z])

        A = np.vstack((A, row1))
        A = np.vstack((A, row2))

    U, D, Vt = np.linalg.svd(A)

    P = Vt[-1].reshape(3,3)
    
    return P.round(5)

def poredjenje_matrica(matrica_dlt, matrica_naivni):

    koeficijent=(matrica_naivni[0][0]/matrica_dlt[0][0]).round(5)
    matrica_dlt=(matrica_dlt/matrica_dlt[0][0])*matrica_naivni[0][0]
    epsilon=0.001
    razlika_veca_od_epsilon=False
    for i in range(len(matrica_dlt)):
        for j in range(len(matrica_dlt[0])):
            if abs(matrica_dlt[i][j]-matrica_naivni[i][j])>epsilon:
                razlika_veca_od_epsilon=True
                if razlika_veca_od_epsilon:
                    print("Matrice nisu jednake")
                    break


    if razlika_veca_od_epsilon==False:
        print(f"Matrice se razlikuju za koeficijent {koeficijent}\n")

print("\nDLT algoritam: ")
matrica_dlt=dlt(original_tacke, slika_tacke)
print(matrica_dlt)

print("\nPoredjenje naivnog i dlt algoritma: ")
poredjenje_matrica(matrica_dlt, matrica_naivni)

def normalizacija(original):
    x = sum([p[0]/p[2] for p in original]) / len(original)
    y = sum([p[1]/p[2] for p in original]) / len(original)
    p = 0.0

    for i in range(len(original)):
        tmp1 = float(original[i][0]/original[i][2]) - x
        tmp2 = float(original[i][1]/original[i][2]) - y

        p = p + math.sqrt(tmp1**2 + tmp2**2)

    p = p / float(len(original))

    S = float(math.sqrt(2)) / p

    return np.array([[S, 0, -S*x], [0, S, -S*y], [0, 0, 1]])

def dlt_normalizacija(original,product):

    
    T = normalizacija(original) #matrica normalizacije 
    Tp = normalizacija(product)
    

    M_line = T.dot(np.transpose(original)) #normalizacija original tacaka
    Mp = Tp.dot(np.transpose(product))  #normalizacija slike tacaka

    M_line = np.transpose(M_line)
    Mp = np.transpose(Mp)

    P_line = dlt(M_line, Mp)

    P = (np.linalg.inv(Tp)).dot(P_line).dot(T)

    return P.round(5)

print("Matrica normalizacije dlt-a:")
matrica_dlt=dlt_normalizacija(original_tacke, slika_tacke).round(5)
print(matrica_dlt)

print("\nPoredjenje naivnog i normalizovanog dlt algoritma: ")
poredjenje_matrica(matrica_dlt, matrica_naivni)


print("\n\n\n-------Testprimeri asistenta----------")

test1=[[2, 1, 1], [1, 2, 1], [3, 4, 1], [-1, -3, 1]]
slika_test1=[[0, 1, 1], [5, 0, 1], [2, -5, 1], [-1, -1, 1]]

A=naivni_algoritam(test1, slika_test1)
print("1)")
print("Naivni:")
print((A/A[0][0]).round(5))
A_dlt=dlt(test1, slika_test1)
print("DLT: ")
print((A_dlt/A_dlt[0][0]).round(5))

A_dlt_norm=dlt_normalizacija(test1, slika_test1)
print("DLT modifikovan")
print((A_dlt_norm/A_dlt_norm[0][0]).round(5))


test2=[[2, 1, 1], [1, 2, 1], [3, 4, 1], [-1, -3, 1], [-2, 5, 1]]
slika_test2=[[0, 1, 1], [5, 0, 1], [2, -5, 1], [-1, -1, 1], [4, 1, 2]]

print()
print("2)")
A_dlt=dlt(test2, slika_test2)
A_dlt_norm=dlt_normalizacija(test2, slika_test2)
print("DLT 5 tacaka")
print((A_dlt/A_dlt[0][0]).round(5))
print("DLT norm 5 tacaka")
print((A_dlt_norm/A_dlt_norm[0][0]).round(5))

print()
print("3)")

test3=[[0, -3, 1], [0, -1, 1], [4, -1, 1], [-7, -4, 1], [0, 5, 1]]
slika_test3=[[3, -1, 1], [4, 4, 1], [9, 1, 1], [5, -2, 1], [7, 2, 2]]

A_dlt_norm=dlt_normalizacija(test3, slika_test3)
print("Dlt norm")
print((A_dlt_norm/A_dlt_norm[0][0]).round(5))
print("-------------------------------------------------\n\n")
print("Kliknite na 4 tacke sa fotografije koje zelite da ispravite")
print("Redosled je bitan: gore-levo, gore-desno, dole-desno, dole-levo")

br=0
A=[[0, 0, 0], 
    [0, 0, 0], 
    [0, 0, 0], 
    [0, 0, 0]]

def click_event(event, x, y, flags, params):
    global br, A, img
    if event==cv2.EVENT_LBUTTONDOWN:
        if br<4:
            A[br][0]=x
            A[br][1]=y
            A[br][2]=1
            br=br+1

            print(x, y)
            cv2.circle(img, (x, y), 4, (25, 255, 25), -1)
            cv2.imshow("zgrada", img)
        if br==4:
            Ap=[[40, 65, 1], [400, 65, 1], [400, 460, 1], [40, 460, 1]]   
            M=naivni_algoritam(A, Ap)
            print("Matrica transformacije M, naivni algoritam")
            print(M)
            out=cv2.warpPerspective(img, M, (700, 960), flags=cv2.INTER_LINEAR)
            
            M_dlt=dlt(A, Ap)
            print("\nMatrica transformacije M, dlt algoritam")
            print(M_dlt)

            poredjenje_matrica(M, M_dlt)
            cv2.imshow("zgrada", out)
            cv2.waitKey(0)
         
img = cv2.imread("zgrada.jpg", cv2.IMREAD_COLOR)
cv2.imshow("zgrada", img)
cv2.setMouseCallback('zgrada', click_event)
cv2.waitKey(0)

