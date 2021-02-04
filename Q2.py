#importing library
import utilities as li
import math as m
L=1
g=9.8
def func(x):
    s=x/m.sqrt(1-((m.sin(m.pi/8)**2)*(m.sin(x)**2)))
    return(s)
eq=func
k=li.Simpson_rule(eq,0,m.pi/2,10)
# find the integral
print("value from simpson method",k)
#find the Time
T=4*m.sqrt(L/g)*k
print("Value of T",T)

#output
# value from simpson method 1.303592370560151
# Value of T 1.665669231727196
