import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as anim
from random import *
from scipy.optimize import fsolve
from scipy.integrate import odeint


i_comp=complex(0,1)

def Newton(xn,F,Jac,eps):    
    x=xn
    while np.linalg.norm(F(x,xn))>eps:
        x=x-np.dot(np.linalg.inv(Jac(x,xn)),F(x,xn)) #x_(n+1)=x_n-inv(J(xn))*F(xn)
    return x



            
## Odre cv kdv

alpha=0.5381301488732363
beta= 0.066633190123881123
a=1.367577724399269
b=0.8234281701082790
c=0.018520783486686603
# 
# alpha=0
# beta=0
# a=3./2.
# b=-3./5.
# c=1./10.

# vitesse=1.
# # 
# u=lambda x,t: 3*vitesse*(1/np.cosh((np.sqrt(vitesse)/2)*((x-vitesse*t)%L-L/2)))**2
# 
# E=[]
# DX=[]
# 
# """Descrétisation spatiale"""
# L=25
# T=2
# J=20
# for m in range(6):
#     dx=L/J
#     DX.append(dx)
#     X=np.linspace(0,L,J+1) 
# 
#     """Descrétisation temporelle"""
#     dt=dx
#     N=int(T/dt)
#     temps=np.linspace(0,T,N+1)
# 
# 
# 
#     A=np.diag(np.ones(J))+np.diag(np.ones(J-1)*alpha,-1)+np.diag(np.ones(J-1)*alpha,1)+np.diag(np.ones(J-2)*beta,-2)+np.diag(np.ones(J-2)*beta,2)
# 
#     A[0][J-1]=alpha
#     A[0][J-2]=beta
#     A[1][J-1]=beta
#     
#     A[J-1][0]=alpha
#     A[J-1][1]=beta
#     A[J-2][0]=beta
# 
#     
#     B=np.diag(np.ones(J-1)*(-a/(2*dx)),-1)+np.diag(np.ones(J-1)*(a/(2*dx)),1)+np.diag(np.ones(J-2)*(-b/(4*dx)),-2)+np.diag(np.ones(J-2)*(b/(4*dx)),2)+np.diag(np.ones(J-3)*(-c/(6*dx)),-3)+np.diag(np.ones(J-3)*(c/(6*dx)),3)
#     
#     B[0][J-1]=-a/(2*dx)
#     B[0][J-2]=-b/(4*dx)
#     B[0][J-3]=-c/(6*dx)
#     B[1][J-1]=-b/(4*dx)
#     B[1][J-2]=-c/(6*dx)
#     B[2][J-1]=-c/(6*dx)
#     
#     B[J-1][0]=a/(2*dx)
#     B[J-1][1]=b/(4*dx)
#     B[J-1][2]=c/(6*dx)
#     B[J-2][0]=b/(4*dx)
#     B[J-2][1]=c/(6*dx)
#     B[J-3][0]=c/(6*dx)
#     
#     
#     A_inv_B=np.matmul(np.linalg.inv(A),B)
#     A_inv_B3=np.matmul(np.matmul(A_inv_B,A_inv_B),A_inv_B)    
#     
#     def F(x,xn):   #application tq F(v^(n+1))=0 + les opération + et * se font points par points sur les vecteurs 
#         return x+(dt/6)*np.dot(A_inv_B,(x*x+x*xn+xn*xn))+(dt/2)*np.dot(A_inv_B3,(x+xn))-xn
#     
#     def Jac(x,xn):  #Jacobienne de F
#         return np.diag(np.ones(J))+dt/2*A_inv_B3+dt/6*np.matmul(A_inv_B,np.diag(2*x+xn))
#     
#     
#     
#     
#     eps=10**(-12)
#     v0=u(X[:J],0)
#     v=v0
#     v0_graphe=np.zeros(J+1)
#     v0_graphe[:J]=v
#     v0_graphe[J]=v[0]
#     V=[v0_graphe]
#     # Hamil.append(H(v0))
#     for n in range(np.size(temps)-1):
#         v=Newton(v,F,Jac,eps)
#         v_graphe=np.zeros(J+1)
#         v_graphe[:J]=v 
#         v_graphe[J]=v[0]
#         V.append(v_graphe)  
#         # Hamil.append(H(v))
#     V_trace=np.array(V)
# 
#     """Erreur de convergence"""
#     erreur=np.linalg.norm(V_trace[0,:]-u(X,0),np.inf)
#     for n in range(1,N+1):
#         if erreur<np.linalg.norm(V_trace[n,:]-u(X,n*dt),np.inf):
#             erreur=np.linalg.norm(V_trace[n,:]-u(X,n*dt),np.inf)
#     E.append(erreur)
#     
#     print('ok',J)
#     J=J+20
#    
# 
# """Tracés"""
# E=np.array(E)
# DX=np.array(DX)
# 
# plt.figure(3)
# plt.clf()
# plt.loglog(DX,E,'r',label='erreur')
# plt.loglog(DX,DX,label='dx')
# plt.loglog(DX,DX**2,label='dx^2')
# plt.loglog(DX,DX**3,label='dx^3')
# plt.loglog(DX,DX**4,label='dx^4')
# plt.loglog(DX,DX**5,label='dx^5')
# plt.loglog(DX,DX**6,label='dx^6')
# plt.loglog(DX,DX**7,label='dx^7')
# plt.xlabel('J')
# # plt.title('Erreur en fonction de dx')
# plt.legend(loc='best')
# plt.show()
# 








##KdV

"""Modélisation du soliton de KdV pour différentes valeurs de C
   ------------------------------------------------------------
"""

# """Déscritisation spatiale"""
# L = 10
# J = 1000
# dx = L/J
# X = np.linspace(-L,L,J+1)
# 
# 
# """Paramétrages des différents solitons"""
# c1 = 1
# c2 = 1./2.
# c3 = 2
# 
# u1 = lambda x,t: 3*c1*(1/np.cosh((np.sqrt(c1)/2)*((x-c1*t))))**2
# u2 = lambda x,t: 3*c2*(1/np.cosh((np.sqrt(c2)/2)*((x-c2*t))))**2
# u3 = lambda x,t: 3*c3*(1/np.cosh((np.sqrt(c3)/2)*((x-c3*t))))**2
# 
# 
# """Tracés""" 
# plt.figure(1)
# plt.clf()
# plt.title('Exemples de solitons pour différentes vitesses')
# plt.plot(X,u1(X,0),'r',label='c=1')
# plt.plot(X,u2(X,0),'b--',label='c=1/2')
# plt.plot(X,u3(X,0),'g-.',label='c=1')
# plt.xlabel('X')
# plt.ylabel('U',rotation=0)
# plt.legend()
# plt.show()



"""Portrait de phase du soliton de KdV
   -----------------------------------
"""

# 
# """Définition des variables"""
# x = np.linspace(-1,3,20)
# y = np.linspace(-2,2,20)
# X,Y = np.meshgrid(x,y)
# 
# vitesse=1
# 
# 
# """Fonction vectorielle définissant l'EDO"""
# def G(X, x=0) : #avec G la fonction vectorielle définissant l'EDO suivante : [u,u']' = [u',c*u-u**2/2] = G([u,u'],x) 
#     return np.array([X[1],vitesse*X[0]-X[0]**2/2.])
# 
# 
# """Tracés"""
# u,v=G([X,Y])
# fig,ax=plt.subplots(figsize=(7,7))
# ax.quiver(X,Y,u,v)
# plt.title('Portrait de phase KdV avec une vitesse égale à 1')
# # set_aspect('equal')
# plt.grid()
# plt.show()


"""Schéma Compact pour KdV
-----------------------------
"""

"""Descrétisation spatiale"""
L = 25
J = 100
dx = L/J
X = np.linspace(0,L,J+1) 


"""Descrétisation temporelle"""
N = 1000
dt = 1/50
T = int(N*dt)
temps = np.linspace(0,T,N+1) 


# 
# """Différents jeux de variables"""
# alpha = 0                         #(C1)
# beta = 0
# a = 3./2.
# b = -3./5.
# c = 1./10.
# 
alpha = 0.5381301488732363        #(C2)
beta = 0.066633190123881123
a = 1.367577724399269
b = 0.8234281701082790
c = 0.018520783486686603


"""Autres variables"""
vitesse = 1.
 
u = lambda x,t: 3*vitesse*(1/np.cosh((np.sqrt(vitesse)/2)*((x-vitesse*t)%L-L/2)))**2

A = np.diag(np.ones(J))+np.diag(np.ones(J-1)*alpha,-1)+np.diag(np.ones(J-1)*alpha,1)+np.diag(np.ones(J-2)*beta,-2)+np.diag(np.ones(J-2)*beta,2)
A[0][J-1] =alpha
A[0][J-2] = beta
A[1][J-1] = beta
A[J-1][0] = alpha
A[J-1][1] = beta
A[J-2][0] = beta

B = np.diag(np.ones(J-1)*(-a/(2*dx)),-1)+np.diag(np.ones(J-1)*(a/(2*dx)),1)+np.diag(np.ones(J-2)*(-b/(4*dx)),-2)+np.diag(np.ones(J-2)*(b/(4*dx)),2)+np.diag(np.ones(J-3)*(-c/(6*dx)),-3)+np.diag(np.ones(J-3)*(c/(6*dx)),3)
B[0][J-1] = -a/(2*dx)
B[0][J-2] = -b/(4*dx)
B[0][J-3] = -c/(6*dx)
B[1][J-1] = -b/(4*dx)
B[1][J-2] = -c/(6*dx)
B[2][J-1] = -c/(6*dx)
B[J-1][0] = a/(2*dx)
B[J-1][1] = b/(4*dx)
B[J-1][2] = c/(6*dx)
B[J-2][0] = b/(4*dx)
B[J-2][1] = c/(6*dx)
B[J-3][0] = c/(6*dx)
 
A_inv_B = np.matmul(np.linalg.inv(A),B)                      #A^{-1}B
A_inv_B3 = np.matmul(np.matmul(A_inv_B,A_inv_B),A_inv_B)     #(A^{-1}B)^{3}

def F(x,xn):   #F_{n}
    return x+(dt/6)*np.dot(A_inv_B,(x*x+x*xn+xn*xn))+(dt/2)*np.dot(A_inv_B3,(x+xn))-xn

def Jac(x,xn):  #J_{n}
    return np.diag(np.ones(J))+dt/2*A_inv_B3+dt/6*np.matmul(A_inv_B,np.diag(2*x+xn))


# """Hamiltonien"""
# def H(v):
#     s=0
#     for i in range(np.size(v)):
#         s+=-1/6*(v[i])**3+1/2*(np.dot(A_inv_B,v)[i])**2
#     return s*dx
# 
# Hamil=[]
# 
# 

# 
# 
"""Schéma"""
eps = 10**(-12)                              #Précision de Newton-Raphson
v0 = u(X[:J],0)
v = v0
v0_graphe = np.zeros(J+1)
v0_graphe[:J] = v                           
v0_graphe[J] = v[0]                          #Condition périodique
V = [v0_graphe]
# Hamil.append(H(v0))

for n in range(np.size(temps)-1):          #Il reste N termes à calculer car le premier terme est connu
    v  =Newton(v,F,Jac,eps)
    v_graphe = np.zeros(J+1)
    v_graphe[:J] = v 
    v_graphe[J] = v[0]                       #Condition périodique
    V.append(v_graphe)  
    # Hamil.append(H(v))
V_trace = np.array(V)
# 
# 
"""Animation"""
fig=plt.figure(2)
plt.clf()
ax=fig.gca()
line1, =ax.plot(X, u(X,0), 'r')
line2, =ax.plot(X, V[0], 'b--')
ax.set_title('Temps initial')
plt.legend(['sol exacte','sol approchée'], loc='best')
# 
# plt.plot(temps,np.array(Hamil))
# 
# 
# """Mise à jour de chaque plot"""
# def runanimate1(n):
#     line1.set_data(X, u(X,dt*n))
#     ax.set_title('Temps t=' + str('{:.3f}'.format(n*dt)))
# 
# def runanimate2(n):
#     line2.set_data(X, V[n])
#     ax.set_title('Temps t=' + str('{:.3f}'.format(n*dt)))
# 
# 
# ani1=anim.FuncAnimation(fig,runanimate1,frames=np.arange(N+1),interval=1,repeat=False)
# ani2=anim.FuncAnimation(fig,runanimate2,frames=np.arange(N+1),interval=1,repeat=False)  
# plt.show()
# 
# 
"""Afficher espace temps"""
fig=plt.figure()
X2,Y2=np.meshgrid(X,temps)
ax=fig.gca(projection='3d')
ax.plot_wireframe(X2,Y2,V_trace,color='r',rstride=100, cstride=0)
plt.title("Tracé filaire KdV")
ax.set_xlabel('X')
ax.set_ylabel('T')
ax.set_zlabel('u')
plt.tight_layout()
plt.show()


# """Heun"""
# def G(x,X):
#     return -X*np.dot(A_inv_B,X)-np.dot(A_inv_B3,X)  
# 
# t = 0   
# w = u(X[:J],0)
# W = [w]
# Hamil.append(H(w))

# for n in range(np.size(temps)-1):
#     k1 = G(t+n*dt,w)
#     k2 = G(t+(n+1)*dt,w+dt*k1)
#     w = w+dt/2*(k1+k2)
#     W.append(w)   #rajouter w[J]=w[0] ?
#     Hamil.append(H(w))


# """Hamiltonien"""
# plt.figure(1)
# plt.clf()
# plt.plot(temps, np.array(Hamil))


# """Animation"""
# fig=plt.figure(2)
# plt.clf()
# ax=fig.gca()
# line1, =ax.plot(X[:J], u(X[:J],0), 'r')
# line2, =ax.plot(X[:J], W[0], 'b--')
# ax.set_title('Temps initial')
# plt.legend(['sol exacte','sol approchée'], loc='best')


# 
# """Mise à jour de chaque plot"""
# def runanimate1(n):
#     line1.set_data(X[:J], u(X[:J],dt*n))
#     ax.set_title('Temps t=' + str('{:.3f}'.format(n*dt)))
# 
# def runanimate2(n):
#     line2.set_data(X[:J], W[n])
#     ax.set_title('Temps t=' + str('{:.3f}'.format(n*dt)))
# 
# 
# ani1=anim.FuncAnimation(fig,runanimate1,frames=np.arange(N+1),interval=1,repeat=False)
# ani2=anim.FuncAnimation(fig,runanimate2,frames=np.arange(N+1),interval=1,repeat=False)  

# """Afficher espace temps"""
# fig=plt.figure()
# X2,Y2=np.meshgrid(X[:J],temps)
# ax=fig.gca(projection='3d')
# ax.plot_wireframe(X2,Y2,np.array(W),color='r',rstride=100, cstride=0)
# ax.set_xlabel('X')
# ax.set_ylabel('T')
# ax.set_zlabel('u')
# plt.tight_layout()
# plt.show()

# 
# """Tracés KdV C1 C2 t=10"""
# 
# """Différents jeux de variables"""
# alpha2=0                         #(C1)
# beta2=0
# a2=3./2.
# b2=-3./5.
# c2=1./10.
# 
# alpha=0.5381301488732363        #(C2)
# beta= 0.066633190123881123
# a=1.367577724399269
# b=0.8234281701082790
# c=0.018520783486686603
# 
# 
# """Descrétisation spatiale"""
# L=30
# J=100
# dx=L/J
# X=np.linspace(0,L,J+1) 
# 
# 
# """Descrétisation temporelle"""
# N=1000
# dt=1/50
# T=int(N*dt)
# temps=np.linspace(0,T,N+1)
# 
# 
# """Autres variables"""
# vitesse=1.
#  
# u=lambda x,t: 3*vitesse*(1/np.cosh((np.sqrt(vitesse)/2)*((x-vitesse*t)%L-L/2)))**2
# 
# A=np.diag(np.ones(J))+np.diag(np.ones(J-1)*alpha,-1)+np.diag(np.ones(J-1)*alpha,1)+np.diag(np.ones(J-2)*beta,-2)+np.diag(np.ones(J-2)*beta,2)
# A[0][J-1]=alpha
# A[0][J-2]=beta
# A[1][J-1]=beta
# A[J-1][0]=alpha
# A[J-1][1]=beta
# A[J-2][0]=beta
# 
# B=np.diag(np.ones(J-1)*(-a/(2*dx)),-1)+np.diag(np.ones(J-1)*(a/(2*dx)),1)+np.diag(np.ones(J-2)*(-b/(4*dx)),-2)+np.diag(np.ones(J-2)*(b/(4*dx)),2)+np.diag(np.ones(J-3)*(-c/(6*dx)),-3)+np.diag(np.ones(J-3)*(c/(6*dx)),3)
# B[0][J-1]=-a/(2*dx)
# B[0][J-2]=-b/(4*dx)
# B[0][J-3]=-c/(6*dx)
# B[1][J-1]=-b/(4*dx)
# B[1][J-2]=-c/(6*dx)
# B[2][J-1]=-c/(6*dx)
# B[J-1][0]=a/(2*dx)
# B[J-1][1]=b/(4*dx)
# B[J-1][2]=c/(6*dx)
# B[J-2][0]=b/(4*dx)
# B[J-2][1]=c/(6*dx)
# B[J-3][0]=c/(6*dx)
#  
# A_inv_B=np.matmul(np.linalg.inv(A),B)                      #A^{-1}B
# A_inv_B3=np.matmul(np.matmul(A_inv_B,A_inv_B),A_inv_B)     #(A^{-1}B)^{3}
# 
# A2=np.diag(np.ones(J))+np.diag(np.ones(J-1)*alpha2,-1)+np.diag(np.ones(J-1)*alpha2,1)+np.diag(np.ones(J-2)*beta2,-2)+np.diag(np.ones(J-2)*beta2,2)
# A2[0][J-1]=alpha2
# A2[0][J-2]=beta2
# A2[1][J-1]=beta2
# A2[J-1][0]=alpha2
# A2[J-1][1]=beta2
# A2[J-2][0]=beta2
# 
# B2=np.diag(np.ones(J-1)*(-a2/(2*dx)),-1)+np.diag(np.ones(J-1)*(a2/(2*dx)),1)+np.diag(np.ones(J-2)*(-b2/(4*dx)),-2)+np.diag(np.ones(J-2)*(b2/(4*dx)),2)+np.diag(np.ones(J-3)*(-c2/(6*dx)),-3)+np.diag(np.ones(J-3)*(c2/(6*dx)),3)
# B2[0][J-1]=-a/(2*dx)
# B2[0][J-2]=-b/(4*dx)
# B2[0][J-3]=-c/(6*dx)
# B2[1][J-1]=-b/(4*dx)
# B2[1][J-2]=-c/(6*dx)
# B2[2][J-1]=-c/(6*dx)
# B2[J-1][0]=a/(2*dx)
# B2[J-1][1]=b/(4*dx)
# B2[J-1][2]=c/(6*dx)
# B2[J-2][0]=b/(4*dx)
# B2[J-2][1]=c/(6*dx)
# B2[J-3][0]=c/(6*dx)
#  
# A_inv_B2=np.matmul(np.linalg.inv(A2),B2)                      #A^{-1}B
# A_inv_B32=np.matmul(np.matmul(A_inv_B2,A_inv_B2),A_inv_B2)     #(A^{-1}B)^{3}
# 
# 
# def F(x,xn):   #F_{n}
#     return x+(dt/6)*np.dot(A_inv_B,(x*x+x*xn+xn*xn))+(dt/2)*np.dot(A_inv_B3,(x+xn))-xn
# 
# def Jac(x,xn):  #J_{n}
#     return np.diag(np.ones(J))+dt/2*A_inv_B3+dt/6*np.matmul(A_inv_B,np.diag(2*x+xn))
# 
# 
# """Schéma"""
# eps=10**(-12)                              #Précision de Newton-Raphson
# v0=u(X[:J],0)
# v=v0
# v0_graphe=np.zeros(J+1)
# v0_graphe[:J]=v                           
# v0_graphe[J]=v[0]                          #Condition périodique
# V=[v0_graphe]
# Hamil.append(H(v0))
# for n in range(np.size(temps)-1):          #Il reste N termes à calculer car le premier terme est connu
#     v=Newton(v,F,Jac,eps)
#     v_graphe=np.zeros(J+1)
#     v_graphe[:J]=v 
#     v_graphe[J]=v[0]                       #Condition périodique
#     V.append(v_graphe)  
#     Hamil.append(H(v))
# V_trace=np.array(V)
# 
# def F2(x,xn):   #F_{n}
#     return x+(dt/6)*np.dot(A_inv_B2,(x*x+x*xn+xn*xn))+(dt/2)*np.dot(A_inv_B32,(x+xn))-xn
# 
# def Jac2(x,xn):  #J_{n}
#     return np.diag(np.ones(J))+dt/2*A_inv_B32+dt/6*np.matmul(A_inv_B2,np.diag(2*x+xn))
# 
# 
# Hamil2=[]
# """Schéma"""
# eps=10**(-12)                              #Précision de Newton-Raphson
# v02=u(X[:J],0)
# v2=v02
# v0_graphe2=np.zeros(J+1)
# v0_graphe2[:J]=v2                          
# v0_graphe2[J]=v2[0]                          #Condition périodique
# V2=[v0_graphe2]
# Hamil2.append(H(v02))
# for n in range(np.size(temps)-1):          #Il reste N termes à calculer car le premier terme est connu
#     v2=Newton(v2,F2,Jac2,eps)
#     v_graphe2=np.zeros(J+1)
#     v_graphe2[:J]=v2 
#     v_graphe2[J]=v2[0]                       #Condition périodique
#     V2.append(v_graphe2)  
#     Hamil2.append(H(v2))
# V_trace2=np.array(V2)

# 
# """Tracés"""
# plt.figure(1)
# plt.clf()
# plt.plot(X,V2[600],color='r',label='C1')
# plt.plot(X,V[600],color='b',label='C2')
# plt.ylim(-0.1,1)
# plt.xlim(4,25)
# plt.xlabel('x')
# plt.ylabel('u',rotation=0)
# # plt.title('C1 et C2 à t=14')
# plt.legend()
# plt.show()
# 
# plt.figure(2)
# plt.clf()
# plt.plot(temps,np.array(Hamil),color='b',label='C2')
# plt.plot(temps,np.array(Hamil2),color='r',label='C1')
# # plt.title('Hamiltonien en fonction de C1 et C2')
# ax = fig.add_subplot(111)
# plt.xlabel('t')
# plt.ylabel('H',rotation=0)
# ax.set_aspect('equal',adjustable='box')
# plt.legend()
# plt.show()


# """Heun"""
# def G(x,X):
#     return -X*np.dot(A_inv_B,X)-np.dot(A_inv_B3,X)  
# 
# t=0   
# w=u(X[:J],0)
# W=[w]
# Hamil.append(H(w))
# for n in range(np.size(temps)-1):
#     k1=G(t+n*dt,w)
#     k2=G(t+(n+1)*dt,w+dt*k1)
#     w=w+dt/2*(k1+k2)
#     W.append(w)   #rajouter w[J]=w[0] ?
#     Hamil.append(H(w))
# 
# 
# """Hamiltonien"""
# plt.figure(1)
# plt.clf()
# plt.plot(temps, np.array(Hamil),label='Heun avec C1')
# plt.xlabel('t')
# plt.ylabel('H',rotation=0)
# plt.legend()
# plt.show()

# 
# """Animation"""
# fig=plt.figure(2)
# plt.clf()
# ax=fig.gca()
# line1, =ax.plot(X[:J], u(X[:J],0), 'r')
# line2, =ax.plot(X[:J], W[0], 'b--')
# ax.set_title('Temps initial')
# plt.legend(['sol exacte','sol approchée'], loc='best')
# 
# 
# 
# """Mise à jour de chaque plot"""
# def runanimate1(n):
#     line1.set_data(X[:J], u(X[:J],dt*n))
#     ax.set_title('Temps t=' + str('{:.3f}'.format(n*dt)))
# 
# def runanimate2(n):
#     line2.set_data(X[:J], W[n])
#     ax.set_title('Temps t=' + str('{:.3f}'.format(n*dt)))
# 
# 
# ani1=anim.FuncAnimation(fig,runanimate1,frames=np.arange(N+1),interval=1,repeat=False)
# ani2=anim.FuncAnimation(fig,runanimate2,frames=np.arange(N+1),interval=1,repeat=False)  
# 
# """Afficher espace temps"""
# fig=plt.figure()
# X2,Y2=np.meshgrid(X[:J],temps)
# ax=fig.gca(projection='3d')
# ax.plot_wireframe(X2,Y2,np.array(W),color='r',rstride=100, cstride=0)
# ax.set_xlabel('X')
# ax.set_ylabel('T')
# ax.set_zlabel('u')
# plt.tight_layout()
# plt.show()



"""Méthode de Tir pour KdV
--------------------------
"""
# 
# vitesse=1.
# 
# u=lambda y: 3*vitesse*(1/np.cosh((np.sqrt(vitesse)/2)*((y))))**2
# # 
# 
# """Descrétisation spatiale""" 
# L=50
# J=40000
# # J2=80000
# # J3=140000
# dx=L/J
# # dx2=L/J2
# # dx3=L/J3
# X=np.linspace(0,L,J+1) 
# # X2=np.linspace(0,L,J2+1) 
# # X3=np.linspace(0,L,J3+1) 
# # alpha=3.0
# # 
# # 
# # # 
# 
# """Définitions des variables"""
# 
# eps=10**(-8)
# a=2
# b=5
# 
# def F(x,X):  
#     # Z=np.zeros(3)
#     # Z[0]=X[1]
#     # Z[1]=X[2]
#     # Z[2]=vitesse*X[1]-X[0]*X[1]
#     Z=np.zeros(2)
#     Z[0]=X[1]
#     Z[1]=vitesse*X[0]-X[0]**2/2.  #- vitesse*alpha +alpha**2/2. +3/2.
#     return Z
#     
# # def G(X,x):  
# #     # Z=np.zeros(3)
# #     # Z[0]=X[1]
# #     # Z[1]=X[2]
# #     # Z[2]=vitesse*X[1]-X[0]*X[1]
# #     Z=np.zeros(2)
# #     Z[0]=X[1]
# #     Z[1]=vitesse*X[0]-X[0]**2/2.  #- vitesse*alpha +alpha**2/2. +3/2.
# #     return Z
# # 
# 
# n=1
# while n>0:
#     alpha=(a+b)/2.
#     Y=np.zeros(2)
#     Y[0]=alpha      
#     # Y[2]=-3./2.
#     y=0              
#     Y_liste=[Y[0]]     #u_tilde
#     Ydx_liste=[Y[1]]   
#     # Ydx2_liste=[Y[2]]
#     for i in range (J): 
#         k1=K(y,Y)
#         k2=K(y+dx/2.,Y+dx/2.*k1)
#         k3=K(y+dx/2.,Y+dx/2.*k2)
#         k4=K(y+dx,Y+dx*k3)
#         Y=Y+dx/6.*(k1+2.*k2+2.*k3+k4)
#         Y_liste.append(Y[0])
#         Ydx_liste.append(Y[1])
#         # Ydx2_liste.append(Y[2])
#         y=y+dx
#     plt.figure(1)
#     plt.clf()
#     Y_tableau=np.array(Y_liste)
#     # plt.plot(X,Y_tableau)
#     # plt.plot(Y_tableau, np.array(Ydx_liste))
#     # plt.axis([0,0.1,0,50])
#     # plt.show()
#     c=0
#     for i in range (len(Y_liste)):
#         # print("Y",Y_liste[i])
#         if Y_liste[i]<0 :   #conteur c 
#             c+=1
#             b=alpha
#             
#     if c==0:
#         for j in range(int(2.*len(Y_liste)/3.),len(Y_liste)):
#             if Y_liste[j]>eps :
#                 # print("aled",Y_liste[j])
#                 a=alpha
#                 n=1
#             else :
#                 n=0
#     print("c",c)
#     print("b",b)
#     print("a",a)
#     print(n)
#     print(alpha)
# 
# 
# 
#  
#  
#  
#         
# 
# plt.figure(1)
# plt.clf()
# Y_tableau=np.array(list(reversed(Y_liste))[:len(Y_liste)-1]+Y_liste)
# X2=np.linspace(-L,L,2*J+1)
# plt.plot(X2,Y_tableau,'b--')
# # plt.plot(X2,u(X2),'r')                         
# plt.show()
  




# 
# eps=10**(-5)
# 
# alpha=3.0
# Y=np.zeros(2)
# Y[0]=alpha
# y=0
# Y_liste=[Y[0]]
# Ydx_liste=[Y[1]]
# for i in range (J): 
#     k1=K(y,Y)
#     k2=K(y+dx/2.,Y+dx/2.*k1)
#     k3=K(y+dx/2.,Y+dx/2.*k2)
#     k4=K(y+dx,Y+dx*k3)
#     Y=Y+dx/6.*(k1+2.*k2+2.*k3+k4)
#     Y_liste.append(Y[0])
#     Ydx_liste.append(Y[1])
#     y=y+dx
# Y_tableau=np.array(Y_liste)
# 
# Y2=np.zeros(2)
# Y2[0]=alpha
# y=0
# Y_liste2=[Y2[0]]
# Ydx_liste2=[Y2[1]]
# for i in range (J2): 
#     k1=K(y,Y2)
#     k2=K(y+dx2/2.,Y2+dx2/2.*k1)
#     k3=K(y+dx2/2.,Y2+dx2/2.*k2)
#     k4=K(y+dx2,Y2+dx2*k3)
#     Y2=Y2+dx2/6.*(k1+2.*k2+2.*k3+k4)
#     Y_liste2.append(Y2[0])
#     Ydx_liste2.append(Y2[1])
#     y=y+dx
# Y_tableau2=np.array(Y_liste2)
# 
# Y3=np.zeros(2)
# Y3[0]=alpha
# y=0
# Y_liste3=[Y3[0]]
# Ydx_liste3=[Y3[1]]
# for i in range (J3): 
#     k1=K(y,Y3)
#     k2=K(y+dx3/2.,Y3+dx3/2.*k1)
#     k3=K(y+dx3/2.,Y3+dx3/2.*k2)
#     k4=K(y+dx3,Y3+dx3*k3)
#     Y3=Y3+dx3/6.*(k1+2.*k2+2.*k3+k4)
#     Y_liste3.append(Y3[0])
#     Ydx_liste3.append(Y3[1])
#     y=y+dx
# Y_tableau3=np.array(Y_liste3)
# 
# 
# 
# plt.figure(5)
# plt.clf()
# # plt.plot(X,odeint(G,[3.0,0],X,rtol=1.e-14,atol=1.e-14)[:,0],'r-',label='odeint')
# # plt.plot(X,Y_tableau,'b--', label='J=40 000')
# # plt.plot(X2,Y_tableau2,'r-',label='J=80 000')
# plt.plot(X3,Y_tableau3,label='J=140 000')
# plt.xlabel('y')
# plt.ylabel('u_tilde')
# plt.legend()
# plt.show()





 ##NLS


"""Modélisation du soliton de NLS pour différentes valeurs de sigma, c et v
   ------------------------------------------------------------------------
"""


# """Déscritisation spatiale"""
# L = 10
# J = 1000
# dx = L/J
# X = np.linspace(-L,L,J+1)
# 
# 
# """Définition des variables"""
# sigma1 = 1.
# c1 = 2.
# v1 = 0.5
# 
# sigma2 = 2.
# c2 = 4.
# v2 = 0.5
# 
# u1 = lambda x,t : (((c1**2-2*c1*v1)/4.*(sigma1+1))/np.cosh(sigma1*np.sqrt(((c1**2-2*c1*v1)/4.))*(x-c1*t))**2)**(1/(2*sigma1))*np.exp(i_comp*c1/2*(x-v1*t))
# u2 = lambda x,t : (((c1**2-2*c1*v1)/4.*(sigma2+1))/np.cosh(sigma2*np.sqrt(((c1**2-2*c1*v1)/4.))*(x-c1*t))**2)**(1/(2*sigma2))*np.exp(i_comp*c1/2*(x-v1*t))
# u3 = lambda x,t : (((c2**2-2*c2*v1)/4.*(sigma1+1))/np.cosh(sigma1*np.sqrt(((c2**2-2*c2*v1)/4.))*(x-c2*t))**2)**(1/(2*sigma1))*np.exp(i_comp*c2/2*(x-v1*t))
# u4 = lambda x,t : (((c1**2-2*c1*v2)/4.*(sigma1+1))/np.cosh(sigma1*np.sqrt(((c1**2-2*c1*v2)/4.))*(x-c1*t))**2)**(1/(2*sigma1))*np.exp(i_comp*c1/2*(x-v2*t))
# u5 = lambda x,t : (((c2**2-2*c2*v2)/4.*(sigma1+1))/np.cosh(sigma1*np.sqrt(((c2**2-2*c2*v2)/4.))*(x-c2*t))**2)**(1/(2*sigma1))*np.exp(i_comp*c2/2*(x-v2*t))
# 
# 
# """Tracés""" 
# plt.figure(1)
# plt.clf()
# # plt.title('Exemples de solitons pour différentes vitesses, sigma')
# plt.plot(X,np.abs(u1(X,0)),'r',label='sigma=1,c=2, v=0.5')
# plt.plot(X,np.abs(u2(X,0)),'b--',label='sigma=2, c=2, v=0.5')
# plt.plot(X,np.abs(u5(X,0)),'g-.',label='sigma=1, c=4, v=1.5')
# plt.xlabel('X')
# plt.ylabel('U',rotation=0)
# plt.legend()
# plt.show()



"""Modélisation filaire du soliton de NLS 
   --------------------------------------
"""


# """Déscritisation spatiale"""
# L = 20
# J = 1000
# dx = L/J
# X = np.linspace(-L,L,J+1)
# 
# 
# """Déscritisation temporelle"""
# N = 500 
# dt = 1/50
# T = int(N*dt)
# temps = np.linspace(0,T,N+1)
# 
# 
# """Définition des variables"""
# sigma1 = 1.
# c1 = 1.
# v1 = 0.2
# 
# u = lambda x,t : (((c1**2-2*c1*v1)/4.*(sigma1+1))/np.cosh(sigma1*np.sqrt(((c1**2-2*c1*v1)/4.))*(x-c1*t))**2)**(1/(2*sigma1))*np.exp(i_comp*c1/2*(x-v1*t))
# 
# 
# """Afficher espace temps"""
# fig=plt.figure()
# X2,Y2=np.meshgrid(X,temps)
# ax=fig.gca(projection='3d')
# # ax.plot_surface(X2,Y2,np.array(W),cmap=cm.coolwarm,linewidth=0)
# ax.plot_wireframe(X2,Y2,np.abs(u(X2,Y2)),color='r',rstride=50, cstride=0)
# plt.title("Tracé filaire NLS")
# ax.set_xlabel('X')
# ax.set_ylabel('T')
# ax.set_zlabel('u')
# plt.tight_layout()
# plt.show()



"""Portrait de phase du soliton de NlS
   -----------------------------------
"""


"""Définition des variables"""
# x = np.linspace(0,2,20)
# y = np.linspace(-0.5,0.5,20)
# X,Y = np.meshgrid(x,y)
# 
# sigma = 1
# c = 2
# v = 0.5 
# 
# 
# X_f1=np.array([0.5,0])  #point critique
# t0 = 0
# tmax= 20
# nbpoints = 1000
# t = np.linspace(t0, tmax,nbpoints) 
# 
# 
# """Fonction vectorielle définissant l'EDO"""
# def G(X, x=0) : #avec G la fonction vectorielle définissant l'EDO suivante : [u,u']' = [u',c*u-u**2/2] = G([u,u'],x) 
#     return np.array([X[1],(c**2-2*c*v)/4*X[0]-X[0]**(2*sigma+1)])
# 
# 
# """Tracés"""
# # f2 = plt.figure()
# # values  = np.linspace(-0.5,0.5 , 10)
# # for v in values:
# #     X0 = v * X_f1                             
# #     X_p = odeint(G,X0,t)
# #     plt.plot( X_p[:,0], X_p[:,1])
# u,v=G([X,Y])
# fig,ax=plt.subplots(figsize=(7,7))
# ax.quiver(X,Y,u,v)
# # plt.title('Portrait de phase KdV avec une vitesse égale à 1')
# ax.set_aspect('equal')
# plt.grid()
# plt.show()



"""Méthode de Tir NLS
---------------------
"""


# """Descrétisation spatiale"""
# L = 20
# J = 2000
# dx = L/J
# X = np.linspace(0,L,J+1) 
# 
# 
# """Définition des variables"""
# eps=10**(-8)
# sigma=1
# 
# a=0.9
# b=2
# 
# u=lambda x,t: ((sigma+1.)*(1./np.cosh(sigma*x))**2)**(1./(2.*sigma))*np.exp(i_comp*t) 
# 
# 
# 
# """Hamiltonien"""
# def Hamil(Y,Ydx):
#     s=0
#     for i in range(len(Y)):
#         s+=(1./2.)*abs(Ydx[i])**2-(1./2.)*(1./(sigma+1.))*abs(Y[i])**(2*sigma+2)
#     return s*dx
# 
# H=[]
# 
# 
# """Fonction F"""
# def F(x,X):
#     Z=np.zeros(2)
#     Z[0]=X[1]
#     Z[1]=(1-(abs(X[0]))**(2*sigma))*X[0]
#     return Z
# 
# 
# """Schéma"""
# n=1
# while n>0:
#     alpha=(a+b)/2.
#     Y=np.zeros(2)
#     Y[0]=alpha           #Conditions initiales              
#     y=0
#     Y_liste=[Y[0]]       #phi_tilde
#     Ydx_liste=[Y[1]]     #dérivée de phi_tilde
#     for i in range (J): 
#         k1=F(y,Y)
#         k2=F(y+dx/2.,Y+dx/2.*k1)
#         k3=F(y+dx/2.,Y+dx/2.*k2)
#         k4=F(y+dx,Y+dx*k3)
#         Y=Y+dx/6.*(k1+2.*k2+2.*k3+k4)
#         Y_liste.append(Y[0])
#         Ydx_liste.append(Y[1])
#         y=y+dx
#     c=0                  #conteur c pour savoir si il existe une valeur <0
#     H.append(Hamil(Y_liste,Ydx_liste))
#     for i in range (len(Y_liste)):
#         if Y_liste[i]<0 :   
#             c+=1
#             b=alpha
#             
#     if c==0:
#         for j in range(int(100.*len(Y_liste)/101.),len(Y_liste)): #dernières valeurs
#             if Y_liste[j]>eps :
#                 a=alpha
#                 n=1
#             else :
#                 n=0     #arrêt car les deux critères ont été satisfaits
# 
# 
# # 
# # """Schéma"""
# # m=2
# # n=1
# # alpha=(a+b)/2
# # Y=np.zeros(2)
# # Y[0]=a
# # y=0
# # Y_liste=[Y[0]]
# # Ydx_liste=[Y[1]]
# # for i in range (J): 
# #     k1=F(y,Y)
# #     k2=F(y+dx/2.,Y+dx/2.*k1)
# #     k3=F(y+dx/2.,Y+dx/2.*k2)
# #     k4=F(y+dx,Y+dx*k3)
# #     Y=Y+dx/6.*(k1+2.*k2+2.*k3+k4)
# #     Y_liste.append(Y[0])
# #     Ydx_liste.append(Y[1])
# #     y=y+dx
# # plt.figure(1)
# # Y_tableau=np.array(Y_liste)
# # plt.plot(X,Y_tableau,'r',label='a_n')
# # c=0
# # 
# # Y=np.zeros(2)
# # Y[0]=b
# # y=0
# # Y_liste=[Y[0]]
# # Ydx_liste=[Y[1]]
# # for i in range (J): 
# #     k1=F(y,Y)
# #     k2=F(y+dx/2.,Y+dx/2.*k1)
# #     k3=F(y+dx/2.,Y+dx/2.*k2)
# #     k4=F(y+dx,Y+dx*k3)
# #     Y=Y+dx/6.*(k1+2.*k2+2.*k3+k4)
# #     Y_liste.append(Y[0])
# #     Ydx_liste.append(Y[1])
# #     y=y+dx
# # Y_tableau=np.array(Y_liste)
# # plt.plot(X,Y_tableau,'b--',label='b_n')
# # c=0
# # 
# # Y=np.zeros(2)
# # Y[0]=alpha
# # y=0
# # Y_liste=[Y[0]]
# # Ydx_liste=[Y[1]]
# # for i in range (J): 
# #     k1=F(y,Y)
# #     k2=F(y+dx/2.,Y+dx/2.*k1)
# #     k3=F(y+dx/2.,Y+dx/2.*k2)
# #     k4=F(y+dx,Y+dx*k3)
# #     Y=Y+dx/6.*(k1+2.*k2+2.*k3+k4)
# #     Y_liste.append(Y[0])
# #     Ydx_liste.append(Y[1])
# #     y=y+dx
# # plt.figure(1)
# # Y_tableau=np.array(Y_liste)
# # plt.plot(X,Y_tableau,'g-.',label='(a_n+b_n)/2')
# # c=0
# # plt.xlabel('x')
# # plt.ylabel('phi_tilde')
# # plt.xlim(0,8)
# # plt.ylim(-2,2)
# # plt.legend()
# # plt.show()
# 
#         
# 
# plt.figure(1)
# plt.clf()
# # Y_tableau=np.array(list(reversed(Y_liste))[:len(Y_liste)-1]+Y_liste)
# # X2=np.linspace(-L,L,2*J+1)
# plt.plot(X,Y_liste,'b--',label='solution approchée')
# plt.plot(X,u(X,0).real,'r-',label='solution exacte')                         
# # Xhamil=np.linspace(0,len(H),len(H))
# # plt.plot(Xhamil,np.array(H))
# # plt.xlabel('Etape de la dichotomie')
# # plt.ylabel('H')
# plt.xlabel('x')
# plt.ylabel('phi_tilde')
# plt.legend()
# plt.show()



"""Schéma compact NLS
---------------------
"""

# 
# """Descrétisation spatiale"""
# L=30
# # dx=0.2
# # J=int(L/dx)
# J=200
# dx=L/J
# X=np.linspace(0,L,J+1) 
# 
# 
# """Descrétisation temporelle"""
# # T=5
# # dt=0.02
# # N=int(T/dt)
# N=2000
# dt=1/50
# T=int(N*dt)
# temps=np.linspace(0,T,N+1)
# 
# 
# """Hamiltonien"""
# def H(v):
#     s=0
#     dxv=np.dot(mat3,v)
#     for i in range(J):
#         s+=1./2.*(dxv[i]**2+dxv[i+J]**2)-1./(2.*sigma+2.)*(v[i]**2+v[i+J]**2)**(sigma+1)
#     return s*dx
#     
# Hamil=[]
# 
# 
# """Norme L_2"""
# def Norme_L(v):
#     s=0
#     for i in range(J):
#         s+=v[i]**2+v[i+J]**2
#     return s*dx
# 
# L_2=[]
# 
# 
# """Moment"""
# def M(v):
#     s=0
#     dxv=np.dot(mat3,v)
#     for i in range(J):
#         s+=(v[i]-i_comp*v[i+J])*(dxv[i]+i_comp*dxv[i+J])
#     return (s*dx).imag
# 
# Mom=[]
# 
# 
# """Variables"""
# alpha=0
# beta=0
# a=3./2.
# b=-3./5.
# c=1./10.
# 
# c_u=2
# v_u=0.5
# sigma=1
# 
# u=lambda x,t: ((((c_u**2-2*c_u*v_u)/4)*(sigma+1))/np.cosh(sigma*np.sqrt((c_u**2-2*c_u*v_u)/4)*(x-c_u*t)%L-L/2)**2)**(1/(2*sigma))*np.exp(i_comp*(c_u*(x-v_u*t))/2)
# 
# A=np.diag(np.ones(J))+np.diag(np.ones(J-1)*alpha,-1)+np.diag(np.ones(J-1)*alpha,1)+np.diag(np.ones(J-2)*beta,-2)+np.diag(np.ones(J-2)*beta,2)
# A[0][J-1]=alpha
# A[0][J-2]=beta
# A[1][J-1]=beta
# A[J-1][0]=alpha
# A[J-1][1]=beta
# A[J-2][0]=beta
# 
# B=np.diag(np.ones(J-1)*(-a/(2*dx)),-1)+np.diag(np.ones(J-1)*(a/(2*dx)),1)+np.diag(np.ones(J-2)*(-b/(4*dx)),-2)+np.diag(np.ones(J-2)*(b/(4*dx)),2)+np.diag(np.ones(J-3)*(-c/(6*dx)),-3)+np.diag(np.ones(J-3)*(c/(6*dx)),3)
# B[0][J-1]=-a/(2*dx)
# B[0][J-2]=-b/(4*dx)
# B[0][J-3]=-c/(6*dx)
# B[1][J-1]=-b/(4*dx)
# B[1][J-2]=-c/(6*dx)
# B[2][J-1]=-c/(6*dx)
# B[J-1][0]=a/(2*dx)
# B[J-1][1]=b/(4*dx)
# B[J-1][2]=c/(6*dx)
# B[J-2][0]=b/(4*dx)
# B[J-2][1]=c/(6*dx)
# B[J-3][0]=c/(6*dx)
# 
# A_inv_B=np.matmul(np.linalg.inv(A),B)
# A_inv_B2=np.matmul(A_inv_B,A_inv_B)  
# 
# mat_echange=np.diag(np.ones(J),J)+np.diag(np.ones(J),-J) #première matrice de F_n
# 
# mat_somme=np.concatenate((np.eye(J),np.eye(J)),axis=1)  #matrice dans la somme de F_n
# 
# mat=np.diag(np.ones(2*J))       #matrice possèdant des -1 sur la moitiée de sa diagonale
# for i in range(J):
#     mat[i+J][i+J]=-1
# 
# mat_opc2=np.zeros((2*J,2*J))  #matrice avec l'operateur compact au carré
# for i in range(J):
#     for j in range(J):
#         mat_opc2[i][j]=A_inv_B2[i][j]
#         mat_opc2[i+J][j+J]=A_inv_B2[i][j]
#         
# mat_opc=np.zeros((2*J,2*J))
# for i in range(J):
#     for j in range(J):
#         mat_opc[i][j]=A_inv_B[i][j]
#         mat_opc[i+J][j+J]=A_inv_B[i][j]
#     
# def somme(sigma,x,xn):   #somme dans F_n
#     s=0
#     for k in range(sigma+1):
#         s+=(np.dot(mat_somme,x*x))**(sigma-k)*(np.dot(mat_somme,xn*xn))**k
#     return s
# 
# def Somme_liste_complex(L1,L2):  #Fonction auxiliaire créant une liste L=L1+iL2 terme par terme
#     n=len(L1)
#     L=[]
#     for i in range(n):
#         L.append(L1[i]+i_comp*L2[i])
#     return L
# 
# def tableau_complex(T1,T2): #Fonction auxiliaire créant le tableau T=T1+iT2 terme par terme
#     n=np.size(T1)
#     T=np.zeros(J+1,dtype=complex)
#     for i in range(n):
#         T[i]=np.complex(T1[i],T2[i])
#     return T
#     
#         
# def F(x,xn):   #application tq F(v^(n+1))=0 
#     return np.dot(mat,dt/2.*np.dot(np.matmul(mat_opc2,mat_echange),x+xn)+dt/(2.*sigma+2.)*np.dot(mat_echange,x+xn)*np.concatenate((somme(sigma,x,xn),somme(sigma,x,xn)),axis=0))+(x-xn)
# 
# def Poly(x,xn):  #P
#     P=np.zeros(J)
#     for j in range(J):
#         s=0
#         for k in range(sigma+1):
#             s+=(x[j]**2+x[j+J]**2)**(sigma-k)*(xn[j]**2+xn[j+J]**2)**k
#         P[j]=s
#     return P
# 
# def da_P(x,xn):  #dérivée partielle de P suivant a
#     da_P=np.zeros(J)
#     for j in range(J):
#         s=0
#         for k in range(sigma):
#             s+=((x[j]**2+x[j+J]**2)**(sigma-k-1))*((xn[j]**2+xn[j+J]**2)**k)*2*x[j]*(sigma-k)
#         da_P[j]=s
#     return da_P
#     
# def db_P(x,xn):   #dérivée partielle de P suivant b
#     db_P=np.zeros(J)
#     for j in range(J):
#         s=0
#         for k in range(sigma):
#             s+=((x[j]**2+x[j+J]**2)**(sigma-k-1))*((xn[j]**2+xn[j+J]**2)**k)*2*x[j+J]*(sigma-k)
#         db_P[j]=s
#     return db_P
# 
# # def Poly(x,xn):
# #     P=np.zeros(J)
# #     for j in range(J):
# #         P[j]=x[j]**2+x[j+J]**2+xn[j]**2+xn[j+J]**2
# #     return P
# # 
# # def da_P(x,xn):
# #     da_P=np.zeros(J)
# #     for j in range(J):
# #         da_P[j]=2*x[j]
# #     return da_P
# #     
# # def db_P(x,xn):
# #     db_P=np.zeros(J)
# #     for j in range(J):
# #         db_P[j]=2*x[j+J]
# #     return db_P
#     
# 
# def Jac(x,xn):  #Jacobienne de F
#     Jac=np.zeros((2*J,2*J))
#     M_Jac=np.zeros((J,J))
#     N_Jac=np.zeros((J,J))
#     E_Jac=np.zeros((J,J))
#     G_Jac=np.zeros((J,J))
#     for i in range(J):
#         M_Jac[i][i]=(xn[i+J]+x[i+J])*da_P(x,xn)[i]
#         N_Jac[i][i]=Poly(x,xn)[i]+(x[i+J]+xn[i+J])*db_P(x,xn)[i]
#         E_Jac[i][i]=Poly(x,xn)[i]+(x[i]+xn[i])*da_P(x,xn)[i]
#         G_Jac[i][i]=(xn[i]+x[i])*db_P(x,xn)[i]
#     Jac[:J,:J]=np.eye(J)+dt/(2*sigma+2)*M_Jac
#     Jac[:J,J:]=dt/2*A_inv_B2+dt/(2*sigma+2)*N_Jac
#     Jac[J:,:J]=-dt/2*A_inv_B2-dt/(2*sigma+2)*E_Jac
#     Jac[J:,J:]=np.eye(J)-dt/(2*sigma+2)*G_Jac
#     return Jac
# 
# 
# 
# 
# V_trace=np.zeros((N+1,J+1),dtype=complex)
# eps=10**(-12)
# v0=u(X[:J],0)
# v=v0
# va=v.real       #séparation des parties imaginaires et réelle
# vb=v.imag
# va0_graphe=np.zeros(J+1) 
# vb0_graphe=np.zeros(J+1)
# va0_graphe[:J]=va   
# va0_graphe[J]=va[0]     #Condition périodique
# vb0_graphe[:J]=vb
# va0_graphe[J]=vb[0]     #Condition périodique
# Va=[va0_graphe]
# Vb=[vb0_graphe]
# V_trace[0,:]=tableau_complex(va0_graphe,vb0_graphe)
# L_2.append(Norme_L(np.concatenate((va,vb),axis=0)))
# Hamil.append(H(np.concatenate((va,vb),axis=0)))
# Mom.append(M(np.concatenate((va,vb),axis=0)))
# 
# for n in range(1,np.size(temps)):
#     xn=np.zeros(2*J)
#     xn[:J]=Va[-1][:J]
#     xn[J:]=Vb[-1][:J]
#     # def F(x):   #application tq F(v^(n+1))=0 + les opération + et * se font points par points sur les vecteurs 
#     #     return np.dot(mat,-dt/2.*np.dot(np.matmul(mat2,mat_echange),x+xn)-dt/(2.*sigma+2.)*np.dot(mat_echange,x+xn)*np.concatenate((somme(sigma,x,xn),somme(sigma,x,xn)),axis=0))+(xn-x)
#     # # print('somme',somme(sigma,xn,xn))
#     # v=fsolve(F,xn,xtol=1.e-13)
#     v=Newton(xn,F,Jac,eps)
#     va_graphe=np.zeros(J+1)
#     vb_graphe=np.zeros(J+1)
#     va_graphe[:J]=v[:J]
#     vb_graphe[:J]=v[J:]
#     va_graphe[J]=va_graphe[0] #Condition périodique
#     vb_graphe[J]=vb_graphe[0] #Condition périodique
#     Va.append(va_graphe)
#     Vb.append(vb_graphe)
#     V_trace[n,:]=tableau_complex(va_graphe,vb_graphe)
#     Hamil.append(H(v))
#     L_2.append(Norme_L(v))
#     Mom.append(M(v))
# V_trace=np.array(Somme_liste_complex(Va,Vb),dtype=complex)
# V_trace=np.array(list_complex(Va,Vb),dtype=complex)
# 
# """Animation"""
# fig=plt.figure(2)
# plt.clf()
# ax=fig.gca()
# line1, =ax.plot(X, np.abs(u(X,0)), 'r')
# line2, =ax.plot(X, np.abs(V_trace[:,0]), 'b--')
# ax.set_title('Temps initial')
# plt.legend(['sol exacte','sol approchée'], loc='best')
# 
# 
# 
# 
# plt.plot(temps,np.array(Hamil))
# plt.plot(temps,np.array(L_2))
# plt.plot(temps,np.array(Mom))
# 
# 
# 
# """Mise à jour de chaque plot"""
# def runanimate1(n):
#     line1.set_data(X, np.abs(u(X,dt*n)))
#     ax.set_title('Temps t=' + str('{:.3f}'.format(n*dt)))
# 
# def runanimate2(n):
#     line2.set_data(X, np.abs(V_trace[:,n]))
#     ax.set_title('Temps t=' + str('{:.3f}'.format(n*dt)))
# 
# 
# ani1=anim.FuncAnimation(fig,runanimate1,frames=np.arange(N+1),interval=1,repeat=False)
# ani2=anim.FuncAnimation(fig,runanimate2,frames=np.arange(N+1),interval=1,repeat=False)  
# 
# plt.show()
# 
# 
# V_trace_abs=np.abs(V_trace)
# 
# """Afficher espace temps"""
# fig=plt.figure()
# X2,Y2=np.meshgrid(X,temps)
# ax=fig.gca(projection='3d')
# ax.plot_wireframe(X2,Y2,V_trace_abs,color='r',rstride=100, cstride=0)
# # ax.plot_wireframe(X2,Y2,np.abs(u(X2,Y2)),color='r',rstride=100, cstride=0)
# ax.set_xlabel('X')
# ax.set_ylabel('T')
# ax.set_zlabel('|u|')
# plt.tight_layout()
# plt.show()


