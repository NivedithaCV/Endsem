
# importing modules
import math as m
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import random
import matplotlib.animation as anim

#read  and write
def read_sqmatrix(A,z,name):
    with open(name,z) as fhand:
        M=[]
        N=[]
        # reading through line
        for line in fhand:
          line=line.rstrip()
          li=list(line.split(","))
          c=len(li)
          M.append(li)
        r=len(M)
        col=len(M[0])
        A=[[0 for y in range(c)]for x in range(r)]

        for i in range(r):
          for j in range(c):
              A[i][j]=float(M[i][j])
        return(A,col)


def write_smatrix(A,name):
    with open(name,'w') as ghand:

        for listitems in A:
            for list in listitems:
                ghand.write('%s,'%list)
            ghand.write('\n')
def read_write(A,z,name):
    if z=="r":
        with open(A,z) as fhand:
            M=[]
            N=[]
            for line in fhand:
              line=line.rstrip()
              li=list(line.split(","))
              c=len(li)
              M.append(li)
            r=len(M)
            c=len(M[0])
            A=[[0 for y in range(c)]for x in range(r)]

            for i in range(r):
              for j in range(c):
                  A[i][j]=int(M[i][j])
            return(A,c)
    if z=="w" or z=="a":
        file1=open(name,z)
        file1.write(A)
        file1.close()
#plotting witha constant function
def plotting_a(n,e,x,y_,u,a):
    plt.plot(n,e,'r',marker=".",label='function')
    plt.axhline(y=a,color='r',linestyle='--',label="analytical value")
    plt.xlabel(x)
    plt.ylabel(y_)
    plt.title(u)
    plt.legend()
    plt.show()
def write_table(a,z):
    with open('your_file.txt', z) as f:
        p=0
        x=' '
        for item in a:
            item=str(item)
            p=30-len(item)
            f.write(item+p*x)
            #f.write("%s" % item)
        f.write("\n")

#partial pivoting
def partial_pivot(a,b,col):
    """This function does partial pivoting of passed matrices"""
    r=len(a)
    for i in range(r):
        if a[i][i]==0:
            for k in range(i,col):
                if k==i or a[i][i]!=0:
                    continue
                else:
                    if abs(a[k][i])>abs(a[i][i]):
                        # c=b[i][col-1]
                        # b[i][col-1]=b[k][col-1]
                        # b[k][col-1]=c
                        for j in range(r):
                            pivot=a[i][j]
                            a[i][j]=a[k][j]
                            a[k][j]=pivot
                        for z in range(col):
                            c=b[i][z]
                            b[i][z]=b[k][z]
                            b[k][z]=c
    return a,b


#Gauss_Jordan elemination
def Gauss_Jordan(a,b,col):
    """Gauss Jordan method of decomposition"""
    for q in range(r):
        pivot=a[q][q]
        for l in range(q,r):
            a[q][l]= a[q][l]/pivot
            b[q][col]=b[q][col]/pivot
        for w in range(r):
            if a[w][q]==0 or q==w:
                continue
            else:
                factor=a[w][q]
                b[w][col]=b[w][col]-factor*b[q][col]
                for c in range(q,r):
                    a[w][c]=a[w][c]-factor*a[q][c]

    return(a,b)


#multiplication of matrice M x N
def mult(M,N):
    E=[[0 for y in range(len(N[0]))]for x in range(len(M))]
    I=[[1,0,0],[0,1,0],[0,0,1]]
    for i in range(len(M)):
        for j in range(len(N[0])):
            for k in range(len(N)):
                E[i][j]+=M[i][k]*N[k][j]
    return(E)


def L_Udec(A):
    for j in range(c_A):
        for i in range(len(A)):

            #diagonal
            if i==j:
                sum=0
                for u in range(i):
                    sum=sum+A[i][u]*A[u][i]
                A[i][i]=A[i][i]-sum

                #elements of upper triangle
            if i<j:
                sum=0
                for k in range(i):
                    sum=sum+A[i][k]*A[k][j]
                A[i][j]=A[i][j]-sum

                #elements of lower triangle
            if i>j:
                sum=0
                for z in range(j):
                    sum=sum+A[i][z]*A[z][j]
                A[i][j]=(A[i][j]-sum)/A[j][j]
    return(A)


def forw_backw(A,B,col):
    r=len(A)
    Y=[[0 for x in range(col)] for y in range(r)]
    X=[[0 for t in range(col)] for w in range(r)]
    for g in range(col):
        #forwardd substitution
        for i in range(r):
            sum=0
            for k in range(i):
                sum=sum+A[i][k]*Y[k][g]
            Y[i][g]=B[i][g]-sum

        #backward substitution
        for l in range(r-1,-1,-1):
            sum=0
            for m in range(l+1,r):
                sum=sum+A[l][m]*X[m][g]
            X[l][g]=(Y[l][g]-sum)/A[l][l]
            X[l][g]=round(X[l][g],4)
    #matrix Y and inverse matrix is X
    return(Y,X)

# bracketing of roots
def bracketing(a,b,equation):
    f_a=equation(a)
    f_b=equation(b)
    if equation(a)*equation(b)<0:
        return(a,b)

    if f_a*f_b>0:
        i=0
        while equation(a)*equation(b)>0 and  i<15:
            if abs(f_a)<abs(f_b):
                    a=a-0.5*(b-a)
                    i=i+1
            else:
                b=b+0.5*(b-a)
                i=i+1

        if equation(a)*equation(b)<0:
            return(a,b)
        if i>14:
            return("please provide another range")


# bisection method
def bisection(a,b,func,A):
    num=[]
    error=[]
    s="{:<16s}{:s}\n".format("iteration","absolute error")
    read_write(s,'w',A)
    i=0;
    while abs(a-b)>=0.000001 and i<200 :
        c=(a+b)/2
        if func(a)*func(c)<0:
            b=c
            i=i+1
            num.append(i)
            error.append(abs(a-b))
            s="{:<16d}{}\n".format(i,abs(a-b))
            read_write(s,'a',A)
        if func(b)*func(c)<0:
            a=c
            i=i+1
            num.append(i)
            error.append(abs(a-b))
            s="{:<16d}{}\n".format(i,abs(a-b))
            read_write(s,'a',A)
    if i<200:
        return(a,b,num,error)
    else:
        return("Please provide another range")


def Regula_falsi(a,b,func,A):
    n=[]
    e=[]
    #starting bracketing
    a,b=bracketing(a,b,func)
    i=0;k=0;error=1
    s="{:<16s}{:s}\n".format("iteration","absolute error")
    read_write(s,'w',A)
    #comparing values
    while abs(a-b)>=0.000001 and i<200 and abs(error)>=0.000001:
        c=b-(((b-a)*func(b))/(func(b)-func(a)))
        error=k-c
        if func(a)*func(c)<0:
            b=c
            k=c
            i=i+1
            if i!=1:
                n.append(i-1)
                e.append(abs(error))
                s="{:<16d}{}\n".format(i-1,abs(error))
                read_write(s,'a',A)
        if func(a)*func(c)>0:
            a=c
            k=c
            i=i+1
            if i!=1:
                n.append(i-1)
                e.append(abs(error))
                s="{:<16d}{}\n".format(i-1,abs(error))
                read_write(s,'a',A)
        if func(c)==0:
            return(c)
    if abs(a-b)<0.000001:
        return(a,n,e)
    if i<200:
        return(c,n,e)
    else:
        print("Please provide another range")


# derivations
def first_derivative(func,x0,h):
    f_=(func(x0+h)-func(x0-h))/(2*h)
    return(f_)

def derivative1(a,func,x0,n,h):
    f_=(func(a,x0+h,n)-func(a,x0-h,n))/(2*h)
    return(f_)

# sec_derivative function
def sec_derivative(func,x0,h):
    f2=(func(x0+h)-2*func(x0)+func(x0-h))/(h*h)
    return(f2)
def derivative2(a,func,x0,n,h):
    f2=(func(a,x0+h,n)-2*func(a,x0,n)+func(a,x0-h,n))/(h*h)
    return(f2)

# Newton_Raphson method
def Newton_Raphson(func,x0,h,A):
    i=0;error=abs(x0);
    N=[]
    ab_er=[]
    s="{:<16s}{:s}\n".format("iteration","absolute error")
    read_write(s,'w',A)
    while error>=0.000001 and i<=200 and func(x0)!=0:
        k=x0
        f_=first_derivative(func,x0,h)
        x0=x0-func(x0)/f_
        error=abs(x0-k)
        i=i+1
        N.append(i)
        ab_er.append(error)
        s="{:<16d}{}\n".format(i,error)
        read_write(s,'a',A)
    if i<200:
        return(x0,N,ab_er)
    else:
        print("Please provide another range")

#polynomial calculator
def polynomial(a,x,n):
    j=0;f=0
    while j<=(n):
        t=a[j]*x**(n-j)
        f=f+t
        j=j+1
    return(f)


# deflection using synthetic division
def deflection(poly,x0,n):
    k=0; d=[0 for y in range(n+1)]
    while k<=n:
        if k==0:
            d[k]=poly[k]
            k=k+1
        else:
            d[k]=(x0*d[k-1])+poly[k]
            k=k+1
    return(d,n-1)

#Laguerres_method function
def Laguerres_method(c,x0,h,n):
    X=[0 for m in range(n)]; z=1;g=n
    while z<=4:
        func=polynomial
        f=func(c,x0,n)
        error=1;i=0
        while error>=0.000001 and i<=200 and func(c,x0,n)!=0:
            G=derivative1(c,func,x0,n,h)/func(c,x0,n)
            H=G**2-(derivative2(c,func,x0,n,h)/func(c,x0,n))
            if G<0:
                x_=x0
                a=n/(G-m.sqrt((n-1)*(n*H-G**2)))
                x0=x_-a
                error=abs(x0-x_)
                i=i+1
            else:
                x_=x0
                a=n/(G+m.sqrt((n-1)*(n*H-G**2)))
                x0=x_-a
                error=abs(x0-x_)
                i=i+1
        x1=x0
        X[z-1]=x0
        c,n=deflection(c,x0,n)
        z=z+1
    return(X)

# code for Midpoint method
def Rectangle_meth(func,a,b,N):
    h=(b-a)/N
    j=1;M=0;i=0
    while j<=N and i<=(N-1):
        x=((a+i*h)+(a+j*h))/2
        M=M+h*func(x)
        j=j+1
        i=i+1
    return(M)


#function for Trapezoidal method
def Trapzoidal_meth(func,a,b,N):
    h=(b-a)/N
    j=0;T=0;
    while j<=N:
        x=a+(j)*h
        if j==0 or j==N:
            T=T+(h/2)*func(a+(j)*h)
        else:
            T=T+(2*(h/2)*func(a+(j)*h))
        j=j+1
    return(T)

#function for simpson technique
def Simpson_rule(func,a,b,N):
    if N%2!=0:
        N=N+1
        print("N made even",N)
    h=(b-a)/N
    j=0;S=0
    while j<=N:
        x=a+(j)*h
        if j==0 or j==N:
            S=S+(h/3)*func(x)
        else:
            if j%2==0:
                S=S+(2*(h/3)*func(a+(j)*h))
            else:
                S=S+4*(h/3)*func(a+(j)*h)

        j=j+1
    return(S)

def N_value(a,b,d,s):
    if s=="M":
        N=m.sqrt((((b-a)**3)/(24*0.001))*abs(d))
    if s=="T":
        N=m.sqrt((((b-a)**3)/(12*0.001))*abs(d))
    if s=="S":
        N=((((b-a)**5)/(180*0.001))*abs(d))**(1/4)

    return(m.ceil(N))



def random_list(N):
    Xi=[]
    for i in range(N):
            n=random.random()
            Xi.append(n)
    return(Xi)

#function for Monte_Carlo
def Monte_Carlo(func,a,b,N):
    S=0;S_2=0
    Xi=random_list(N)
    for x in Xi:
        S=S+func(x)
        S_2=S_2+(func(x)**2)
    Fn=((b-a)/N)*S
    variance=((1/N)*S_2)-((1/N)*S)**2
    return(Fn,variance)

#Assignment 7
# solving differential equation
#function for Euler method
def Explicit_euler(x0,y_x0,h,func,fun,xn):
    X=[[0 for c in range(250)] for r in h]
    Y=[[0 for c_ in range(250)] for r_ in h]
    Z=[[0 for c__ in range(250)] for r_ in h]
    z=[[0 for c__ in range(250)]for r_ in h]
    i=0
    w=0
    # different values of h
    for k in h:
        y_x=y_x0
        x=x0
        while i<250:
            X[w][i]=x
            Y[w][i]=y_x
            z[w][i]=x
            y_xn=y_x+k*func(y_x,x)
            z_xn=fun(x)
            Z[w][i]=z_xn
            i=i+1
            y_x=y_xn
            x=x+k

        i=0
        w=w+1
    for d in h:
        plt.plot(X[i],Y[i],"--",label=str(d)+'solution')
        i=i+1
    plt.plot(z[0],Z[0],label='analytical')
    plt.xlim([X[0][0],X[0][249]])
    plt.ylim(Y[0][0],Y[0][249])
    plt.legend()
    plt.show()

# # function for Runge-Kutta 4 oder
# def RK4(x0,y_x0,v_x0,h,f,g,analy,xn):
#     X=[[0 for c in range(250)] for r in h]
#     Y=[[0 for c_ in range(250)] for r_ in h]
#     Z=[[0 for c__ in range(250)] for r_ in h]
#     z=[[0 for c__ in range(250)]for r_ in h]
#     #a=x,b=y_x0
#     i=0
#     w=0
#     #different vvalues of h
#     for k in h:
#         y=y_x0
#         x=x0
#         v=v_x0
#         while i<250 :
#             X[w][i]=x
#             Y[w][i]=y
#             z[w][i]=x
#             k1y=k*f(x,y,v)
#             k1v=k*g(x,y,v)
#             k2y=k*f(x+k/2,y+k/2,v+k1v/2)
#             k2v=k*g(x+k/2,y+k/2,v+k1v/2)
#             k3y=k*f(x+k/2,y+k/2,v+k2v/2)
#             k3v=k*g(x+k/2,y+k/2,v+k2v/2)
#             k4y=k*f(x+k/2,y+k/2,v+k3v/2)
#             k4v=k*g(x+k/2,y+k/2,v+k3v/2)
#             z_xn=analy(x)
#             Z[w][i]=z_xn
#             i=i+1
#             y=y+1/6*(k1y+(2*k2y)+(2*k3y)+k4y)
#             v=v+1/6*(k1v+(2*k2v)+(2*k3v)+k4v)
#             x=x+k
#
#         i=0
#         w=w+1
#
#     j=0
#     for d in h:
#         plt.plot(X[j],Y[j],"--",label=str(d)+'solution')
#         j=j+1
#     i=0
#
#     plt.plot(z[0],Z[0],'r',label='analyticaal solution')
#     plt.plot(z[1],Z[1],'r')
#     plt.xlim([-6,6])
#     plt.ylim([-5,5])
#     plt.legend()
#     plt.show()
#     return(X,Y,Z)


#
# Assignmen 8 functions
# three D plot
def threeDplot(X,Y,Z):
    fig=plt.figure()
    ax=plt.axes(projection='3d')
    ax.set_xlabel('X axis')
    ax.set_ylabel("Y axis")
    ax.set_zlabel("z axis")
    plt.title("3D plot of ellipsoid Monte Carlo data")
    ax.scatter3D(X,Y,Z,c=X,cmap='rainbow');
    plt.show()
# potting function
def plotting(n,e,x,y,u):
    plt.plot(n,e,'r',marker=".",label='function')
    #plt.axhline(m.pi,label='theoretical pi')
    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(u)
    plt.legend()

# finding volume of ellipsoid using monte carlo method
def M_C_volume(n,l,f):
    M=[];X=[];Y=[];Z=[]
    N=[];Fractional=[]
    j=0
    list=l
    n_list=[-x for x in list]
    # Number of steps
    for j in range(100,20000,2000):
        ellipsoid_pts=0
        i=0;
        while i<=j:
            x=random.uniform(n_list[0],list[0])
            y=random.uniform(n_list[1],list[1])
            z=random.uniform(-2,2)
            pts=f(x,y,z)
            if j==10100 and pts<=1:
                X.append(x)
                Y.append(y)
                Z.append(z)
            if pts<=1:
                ellipsoid_pts=ellipsoid_pts+1
            i=i+1

            # plotting values
        # if j==10100:
        #     threeDplot(X,Y,Z)
        ellipse=(ellipsoid_pts/j)*(2*3*4)
        Fractional.append(abs((ellipse-12.57)/ellipse))
        M.append(ellipse)
        N.append(j)

    # plt.plot(N,M,label="result")
    # plt.axhline(y=12.57,color='r',linestyle='--',label="analytical value")
    # plt.xlabel("Step number")
    # plt.ylabel("result")
    # plt.title("Compare result to analytical value")
    # plt.legend()
    # plt.show()
    return(M,N,Fractional,X,Y,Z)

# function on random walk simulation on a 2 D plane
def random_walk(n):
    Rs=0;R=0
    r=0
    dX=0;dY=0
    X=[[0 for k in range(n+1)] for l in range(5)]
    Y=[[0 for k in range(n+1)] for l in range(5)]
    R2=[0 for i in range (100)]
    # taking 100 different random walks
    while r<100:
        x=0
        y=0
        random.seed(None)
        # takiing n steps
        for i in range(n):
            step =random.random()
            x=x+m.cos(2*m.pi*step)
            y=y+m.sin(2*m.pi*step)
            if r<5:
                X[r][i+1]=x
                Y[r][i+1]=y
        R=R+m.sqrt(x**2+y**2)
        R2[r]=x**2+y**2
        dX=dX+x
        dY=dY+y
        Rs=Rs+R2[r]
        if r<5:
            plt.plot(X[r],Y[r],label="random walk {}".format(r+1))
        r=r+1
    # plotting and calculating average values
    plt.xlabel("x axis")
    plt.ylabel("y axis")
    plt.legend()
    plt.title("Random walks of {} steps".format(n))
    plt.show()
    Rs=Rs/100
    Dx=dX/100
    Dy=dY/100
    R_av=R/100
    Rrms=m.sqrt(Rs)

    return(X,Y,Rrms,Dx,Dy,R_av)
# least square fittiing withoutweight
def Leastsquarefitting(X,Y):
    n=len(X)
    x_=0;y_=0;x_2=0;y_2=0
    for i in range(len(X)):
        x_=x_+X[i]
        y_=y_+Y[i]
        x_2=x_2+X[i]**2
        y_2=y_2+Y[i]**2
        x_y=X[i]*Y[i]
    # find the parameters
    x_av=x_/n
    y_av=y_/n
    a=((y_av*x_2)-(x_av*x_y))/(x_2-(n*(x_av**2)))
    b=(x_y-(n*x_av*y_av))/(x_2-n*(x_av**2))
    Sxx=x_2-n*(x_av**2)
    Syy=y_2-(n*y_av**2)
    Sxy=x_y-(n*x_av*y_av)
    sigma_x2=Sxx/n
    sigma_y2=Syy/n
    r=(Sxy**2)/(Sxx*Syy)
    corr=Sxy/n
    return(a,b,r)

def Rungekutta(x0,y0,z0,h,f,g,x_n):
    X=[x0];Y=[y0];Z=[z0]
    x=x0;y=y0;z=z0
    N = int(abs((x_n - x0))/h)
    for n in range(N):
        k1y = h*f(x,y,z)
        k1z = h*g(x,y,z)
        k2y = h*f(x+h/2,y,z + k1z/2)
        k2z = h*g(x + h/2,y + k1y/2,z+ k1z/2)
        k3y = h*f(x,y,z + k2z/2)
        k3z = h*g((x + h/2),(y + k2y/2),(z+k2z/2))
        k4y = h*f(x,y,z + k3z)
        k4z = h*g((x + h),(y + k2y),(z+k3z))
        x = x + h
        y = y + 1/6*(k1y+(2*k2y)+(2*k3y)+k4y)
        z = z + 1/6*(k1z+(2*k2z)+(2*k3z)+k4z)
        X.append(x)
        Y.append(y)
        Z.append(z)
    return(X,Y,Z)


def fit(X,func):
    Y=[]
    for x in X:
        Y.append(func(x))
    return(Y)

def shooting_method(u_f,v_f,h,a,b,y_a,y_b,func1,func2):
    # takiing value s from Rengakutta function to check the guess
    M,N1,Z1=Rungekutta(a,y_a,func1,h,v_f,u_f,b)
    while(round(N1[-1]-y_b,5)>0):
        return(0,0,0,0)
    if round(N1[-1]-y_b,5)!=0:
        M,N2,Z1=Rungekutta(a,y_a,func2,h,v_f,u_f,b)
        while(round(N2[-1]-y_b,5)<0):
            return(0,0,0,0)
        #intrapolation
        if(round(N2[-1]-y_b,5)==0):
            return(M,N2,func2)
        z=func1+(func2-func1)*(y_b-N1[-1])/(N2[-1]-N1[-1])
        M,N1,Z1=Rungekutta(a,y_a,z,h,v_f,u_f,b)
        return(M,N1,z,Z1[0])
    else:
        return(M,N1,func1,Z1[0])
