import utilities as li
import math as m
import numpy as nm
import matplotlib.pyplot as plt
i=0;X=[];Y=[];Z=[]
#reading the matrix
A,col=li.read_sqmatrix("A","r","esem_table.dat")
n=len(A)
for i in range(n):
    X.append(A[i][0])
    Y.append(A[i][1])
n=len(Y)
def func1(x):
    pass

f1=func1
for i in range(n):
    Z.append(m.log(Y[i],m.e))
# fitting for the ist and 2nd

plt.plot(X,Y,'r')

a,b,r=li.Leastsquarefitting(X,Y)
print("values for(i):omega0,omegac and Pearson r resptively:",a,b,r)

a,b,r=li.Leastsquarefitting(X,Z)
print("values for (ii):omega0,omegac and Pearson r resptively",a,b,r)

#output
#     values for(i):omega0,omegac and Pearson r resptively: 4.154487179487182 -1.762820512820514 13.383563201254486
#     values for (ii):omega0,omegac and Pearson r resptively 0.7028121255682788 -0.3425868389257727 0.7486359806957483
