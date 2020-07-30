from matplotlib.pyplot import *
from numpy import *
import csv

pois=csv.reader(open('plot_LBM.csv','r'),delimiter=',')

ul=[]
yl=[]
k=0

for i in pois:
    if k>0 and k%20==1:
        yl.append(float(i[7]))
        ul.append(float(i[1]))
    k+=1
        
F=1e-7
rho=1
nu=0.1/3
mu=nu/rho
h=99

def u(y):
    return F/(2*mu)*y*(h-y)

y=linspace(0,99,50)

umaxLBM=max(ul)
umax=u(h/2)

ul2=[]
for i in ul:
    i=i/umaxLBM
    ul2.append(i)
u2=[]
for i in u(y):
    i=i/umax
    u2.append(i)
y2=[]  
for i in y:
    i=i/h
    y2.append(i)
yl2=[]  
for i in yl:
    i=i/h
    yl2.append(i)

figure(1)
xlabel('y/h')
ylabel('u(y)/umax')
plot(y2,u2,'r-',label='analytical')
plot(yl2,ul2,'kx',markerfacecolor='none',label='simulation')
legend()
tight_layout()
savefig('poiseuille_profile_scaled.png')
savefig('poiseuille_profile_scaled.eps')
figure(2)
xlabel('y')
ylabel('velocity')
plot(y,u(y),'r-',label='analytical')
plot(yl,ul,'kx',label='simulation')
legend()
tight_layout()
savefig('poiseuille_profile.png')
savefig('poiseuille_profile.eps')
show()
