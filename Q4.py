import utilities as li
import math as m
import matplotlib.pyplot as pl
g=9.8;Z=[]
#function writen as less order
def func(x,y,z):
    s=z
    return(s)
def funct_z(x,y,z):
    s=-g
    return(s)
def func3(x):
    y=2*m.exp(x)-x-1
    retyrn(y)
func1=func
func2=funct_z
func3=func3
h=[0.2]
#calll for shooting method
X,Y,z = (li.shooting_method(func2,func1,0.02,0,2,5,45,0,100))
print("The initial velocity",Z)
#output

# The initial velocity=29.999
