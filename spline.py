import numpy as np
import matplotlib.pyplot as plt

############# set x and y values for interpolation points (x,y); afterwards simply run the script
xval = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]).astype(float)
yval = np.array([5, 5, 6, 1, 4, 3, 1, 9, 4, 6, 4]).astype(float)
#############

n    = int(xval.size)
N    = 1000      # number of spline points (more points = smoother spline graph)

mainDiag  = np.array([ 2*(xval[i+2]-xval[i]) for i in range(0,n-3+1)])
upperDiag = np.array([ (xval[k+1]-xval[k])   for k in range(1,n-3+1)])
RHS       = np.array([ 6*((yval[k+1]-yval[k])/(xval[k+1]-xval[k]) - (yval[k]-yval[k-1])/(xval[k]-xval[k-1])) for k in range(1,n-2+1)])

# solve system of linear equations with tridiagonal symmetric matrix (a mainDiag, u upperDiag, b RHS)
def solve(a, u, b):
    k = b.size-1
    c = np.zeros(k+1)
    # produce upper triangle matrix
    for i in range(1,k+1):
        dummy  = u[i-1]/a[i-1]
        a[i]   = a[i] - u[i-1] * dummy
        b[i]   = b[i] - b[i-1] * dummy
    # "backwards substitution"
    c[k] = b[k]/a[k]
    for i in range(k-1,0-1,-1):
        c[i] = (b[i]-c[i+1]*u[i])/a[i]
    return c

c = np.append(0., np.append(solve(mainDiag, upperDiag, RHS), 0.))

def spline(leftindex, x):
    j = leftindex + 1
    return yval[j-1]*(xval[j]-x)/(xval[j]-xval[j-1]) + yval[j]*(x-xval[j-1])/(xval[j]-xval[j-1]) - c[j-1]/6*( (xval[j]-x)*(xval[j]-xval[j-1]) - (xval[j]-x)**3/(xval[j]-xval[j-1]) ) - c[j]/6*( (x-xval[j-1])*(xval[j]-xval[j-1]) - (x-xval[j-1])**3/(xval[j]-xval[j-1]) )

intervall = np.linspace(np.min(xval),np.max(xval), N)
splineval = np.zeros(N)

for i in range(N):
    leftindex = 0
    for m in range(n-1):
        if intervall[i] >= xval[m]:
            leftindex = m
        else:
            break

    splineval[i] = spline(leftindex, intervall[i])

########################## alternatively just output intervall and splinevall to get the spline points (uncomment next line and comment all the plot lines)
#print(intervall)
#print(splineval)
plt.plot(xval, yval, '.')
plt.plot(intervall, splineval)
plt.show()