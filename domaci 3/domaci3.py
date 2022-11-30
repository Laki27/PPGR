import numpy as np
import math

def Euler2A(fi, teta, psi):
    
    Rx = np.array([[1, 0, 0],
                    [0, math.cos(fi), -math.sin(fi)],
                    [0, math.sin(fi), math.cos(fi)]])

    Ry = np.array([[math.cos(teta), 0, math.sin(teta)],
                    [0, 1, 0],
                    [-math.sin(teta), 0, math.cos(teta)]])

    Rz = np.array([[math.cos(psi), -math.sin(psi), 0],
                    [math.sin(psi), math.cos(psi), 0],
                    [0, 0, 1]])


    A = (Rz.dot(Ry)).dot(Rx)
    return A

n=224
fi=(2*math.pi/6)*(n%5+1)
teta=(math.pi/17)*(n%8+1)
psi=(2*math.pi/8)*(n%7+1)

A=Euler2A(fi, teta, psi)
print("Euler2A")
print(A)

def AxisAngle(A):

    eps=0.001
    E=np.eye(3)
    prvi_uslov=np.dot(A, A.transpose()).round(5)-E
    if prvi_uslov.all()<eps:
        prvi_uslov=True
    drugi_uslov=(np.linalg.det(A)-1.0)<=eps
    if not (prvi_uslov and drugi_uslov):
        print("Matrica A nije matrica kretanja")
        return -1, -1    
    X = A - E

    r = np.array([X[0][0], X[0][1], X[0][2]])
    q = np.array([X[1][0], X[1][1], X[1][2]])

    p=np.cross(r, q)
    if p[0]==0 and p[1]==0 and p[2]==0:
        q=np.array([X[2][0], X[2][1], X[2][2]])
        p=np.cross(r, q)

    p_norm = math.sqrt(p[0]**2 + p[1]**2 + p[2]**2)
    p = p * (1/p_norm)

    
    r_norm = math.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
    u = r * (1/r_norm) 

    up = A.dot(u)

    cos_fi = np.dot(u, up)
    fi = math.acos(cos_fi)

    return p, fi

p, fi = AxisAngle(A)
print("\nAxisAngle\nPrava p i ugao fi")
print(p, fi)

p=np.array([p])

def rodrigez(p, fi):

    E=np.eye(3)
    
    R=(1+math.cos(fi))*p.transpose().dot(p)
    R=R+math.cos(fi)*E

    p_x = np.array([[0, -p[0][2], p[0][1]],
                    [p[0][2], 0, -p[0][0]],
                    [-p[0][1], p[0][0], 0]])

    R = R + math.sin(fi)*p_x

    return R

print("\nRodrigez")
print(rodrigez(p, fi))

def A2Euler(A):
    eps=0.001
    uslov = np.linalg.det(A)
    if uslov-1.0>eps:
        print("Matrica A nije ortogonalna!")
        return -1, -1, -1
    
    psi=0 
    teta=0 
    fi=0
    if A[2][0]<1:
        if A[2][0]>-1:
            psi=math.atan2(A[1][0], A[0][0])
            teta=math.asin(-A[2][0])
            fi=math.atan2(A[2][1], A[2][2])
        else:
            psi=math.atan2(-A[0][1], A[1][1])
            teta=math.pi/2
            fi=0
    else:
        psi=math.atan2(-A[0][1], A[1][1])
        teta=-math.pi/2
        fi=0
    return fi, teta, psi         

print("\nA2Euler")
print(A2Euler(A))

print("\nAxisAngle2Q")

def AxisAngle2Q(p, fi):
    
    w = math.cos(fi/2)
    p=p*(1/(math.sqrt(p[0][0]**2+p[0][1]**2+p[0][2]**2)))
    r = math.sin(fi/2) * p
    q = np.array([r[0][0],r[0][1], r[0][2], w])

    return q

q=AxisAngle2Q(p, fi)
print(q)

print("\nQ2AngleAxis")
def Q2AngleAxis(q):
    q=q*(1/(math.sqrt(q[0]**2+q[1]**2+q[2]**2+q[3]**2)))
    if(q[3]<0):
        q=-q
    fi=2*math.acos(q[3])
    print(q[3])
    if(q[3]==1):
        p=np.array([1, 0, 0])
    else:
        p=np.array([q[0], q[1], q[2]])
        p=p*(1/(math.sqrt(p[0]**2+p[1]**2+p[2]**2)))

    return p, fi   
p, fi=Q2AngleAxis(q)
print(p, fi)
