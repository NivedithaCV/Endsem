#import library
import utilities as li
import math as m

#function to be solved
def func(x):
    s=(x-5)*m.exp(x)+5
    return(s)
h=6.626*10**(-34)
k=1.381*10**(-23)
c=3*10**(8)
f=func
# calling for function to find the root
x0,n,abs=li.Newton_Raphson(f,1.5,0.0001,'Q1.txt')
print("root from Newton_Raphson ",x0)
b=x0*k/h*c
print("Wein'sconstant b=",b)

# output:
#     root from Newton_Raphson  2.3356212975945758e-17
# Wein'sconstant b= 146.03801744543205
