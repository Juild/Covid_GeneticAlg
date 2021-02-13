#GAplot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

colora="red"
colori="blue"
colory="green"
colorr="purple"
colord="black"
truedata=np.loadtxt("truedata.txt",delimiter=",")
At=truedata[:,0]
It=truedata[:,1]
Yt=truedata[:,2]
Rt=truedata[:,3]
Dt=truedata[:,4]
t=np.arange(len(At))


expdata=pd.read_csv("fituni.txt")
print(expdata)
expdata["Pt"]=expdata.sum(axis=1)
expdata["P0"]=1e6
aex=expdata["A_d"]
iex=expdata["I_2"]
yex=expdata["Y"]
rex=expdata["R"]
dex=expdata["D"]
ae=aex.to_numpy()
ie=iex.to_numpy()
ye=yex.to_numpy()
re=rex.to_numpy()
de=dex.to_numpy()
print(expdata)
for i in range(len(ae)):
	print(de[i],Dt[i+1])

t1=np.arange(1,len(At))

plt.plot(t,At,color=colora,label="asymptomatic",alpha=0.8)
plt.scatter(t1,ae,color=colora,marker="+")

plt.plot(t,It,color=colori,label="infected",alpha=0.8)
plt.scatter(t1,ie,color=colori,marker="+")

plt.plot(t,Yt,color=colory,label="strong symptoms",alpha=0.8)
plt.scatter(t1,ye,color=colory,marker="+")


plt.plot(t,Rt,color=colorr,label="recovered",alpha=0.8)
plt.scatter(t1,re,color=colorr,marker="+")

plt.plot(t,Dt,color=colord,label="deaths",alpha=0.8)
plt.scatter(t1,de,color=colord,marker="+")
plt.yscale("log")
plt.legend()
plt.grid()
plt.title("uniform")
plt.xlabel("time(days)")
plt.ylabel("number of individuals")
plt.show()
