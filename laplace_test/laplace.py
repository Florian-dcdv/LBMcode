from matplotlib.pyplot import *
from numpy import *
from scipy.optimize import curve_fit

tau='0.6'   #For the Laplace test
#tau='0.6a' #For surface tension adjustment with parameter 'a' from CS EoS

dp=[]
R=[]
st=[]
Tr=['0.50','0.70','0.95']   #For the Laplace test
#Tr=['a0.25','a0.40','a0.45','a0.50']   #For surface tension adjustment with parameter 'a' from CS EoS
text='_tau0.6'
def func(x,a):
    return a*x

for j in Tr:
    file=open("tau"+tau+"/laplaceTest_"+j+".dat","r")
    line=file.readline()
    #line=file.readline()
    for i in range(4):
        line=file.readline()
        listp=line.split()
        dp.append(float(listp[1])-float(listp[2]))
        R.append(1/float(listp[0]))
        #print(R,dp)
    file.close()
    gamma,err=curve_fit(func,R,dp)
    st.append(gamma[0])
    print("Tr=",j, "gamma=","%.3e"%gamma[0], "error=",err[0][0])
    x=linspace(0,max(R),1000)
    fig=figure(1)
    ax=fig.add_subplot(111)
    #title("Laplace test")
    ax.set_xlabel("1/R")
    ax.set_ylabel("Pressure difference")
    
    if j=='0.50' or j=='a0.25':
        plot(R,dp,'ro',label="Tr="+j,markersize=8,markerfacecolor='none')
        plot(x,gamma*x,'r-',label="gamma="+"%.3e"%gamma[0])
    elif j=='0.70' or j=='a0.40':
        plot(R,dp,'bs',label="Tr="+j,markersize=8,markerfacecolor='none')
        plot(x,gamma*x,'b-',label="gamma="+"%.3e"%gamma[0])
    elif j=='0.95' or j=='a0.45':
        plot(R,dp,'g^',label="Tr="+j,markersize=8,markerfacecolor='none')
        plot(x,gamma*x,'g-',label="gamma="+"%.3e"%gamma[0])
    elif j=='0.55' or j=='a0.50':
        plot(R,dp,'kx',label="Tr="+j,markersize=8,markerfacecolor='none')
        plot(x,gamma*x,'k-',label="gamma="+"%.3e"%gamma[0])
    R=[]
    dp=[]
legend()
tight_layout()
savefig("tau"+tau+"/Laplace test"+text+".eps")
savefig("tau"+tau+"/Laplace test"+text+".png")
show()
if tau=='0.6a': #Surface tension adjustment
    Tr2=[]
    for i in Tr:
        i=i.replace('a','')
        i=float(i)
        Tr2.append(i)
    
    def func2(x,a,b):
        return a*x+b
    res,error=curve_fit(func2,Tr2,st)
    print("a=",res[0],"b=",res[1],"error=",error[0][0])
    x=linspace(min(Tr2),max(Tr2),1000)
    
    figure(2)
    xlabel("a")
    ylabel("gamma")
    plot(Tr2,st,'ro',markersize=8,markerfacecolor='none')
    plot(x,func2(x,res[0],res[1]),'r-',label="fitting curve \n y=ax+b \n a="+"%.3e"%res[0]+"\n b="+"%.3e"%res[1])
    legend()
    tight_layout()
    savefig("tau0.6a/surface tension adjustment.eps")
    savefig("tau0.6a/surface tension adjustment.png")
    show()
