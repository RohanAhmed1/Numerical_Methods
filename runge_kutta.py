# Required Libraries
from math import sqrt
import matplotlib.pyplot as plt

# Global Values

# values according to my roll #
# only these changes according to your roll #, 
a = 91
b = 1.23

# It can be guessable 
# (c1+c2+c3+c4) == 1

# guess only the w1, w2, t1, t2 , you can change the weightagesa also(ci)

# I assume, all weightages of ki is same
c1, c2, c3, c4 = 0.25, 0.25, 0.25, 0.25 
# location of slope computation
w0, w1, w2, w3 = 0, 0.5, 0.5 , 1
t0, t1, t2, t3 = 0, 0.5, 0.5, 1

# no. of intervals -> 40 or 80
n = 40

# initial values
x = 0
y = round( (b-1) / b, 4 )

# y required at xn  = 1.3162
xn = round( (a*b)/(a+b) , 4) 
print("xn =",xn)

# step size
h = (xn - x) / n
print("h = ", h)





# There is a values defined in the global space. Therefore, these values can access by each functions.


# define dy/dx
def dydx(x,y):
    return ( (b-1)/(a+b) ) * x*y




# Standard Four Stage Runge Kutta Method
def S4RK():
    # initial values
    xi = x
    yi = y
    # this list has ans at every steps in domain i.e. 0 <= x <= 10 ---  
    SRK_domain = []
    for i in range(n):
        # k values
        k1 = h * dydx(xi, yi)     # w0 = 0, that's why the k0 is not present
        k2 = h * dydx(xi + (1/2) * h, yi + (1/2) * k1)
        k3 = h * dydx(xi + (1/2) * h, yi +(1/2) * k2)
        k4 = h * dydx(xi + h, yi + k3)
    
        # xi and yi at every iteration starting from (x1,y1) to (xn, yn)
        xi = xi + h
        yi = yi + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)

        SRK_domain.append(yi)

    # value at last value of domain i.e. at x = 10, calculate y [ y(10) = ? -> Generally it is required ]
    return SRK_domain[-1], SRK_domain



# Modified Four Stage Runge Kutta Method
def M4RK():
    # initial values
    xi = x
    yi = y
    # this list has ans at every steps in domain i.e. 0 <= x <= 10 ---  
    MRK_domain = []
    for i in range(n):
        # k values
        k1 = h * dydx(xi + t0*h, yi + w0)     # w0 = 0, that's why the k0 is not present
        k2 = h * dydx(xi + t1*h, yi + w1*k1)
        k3 = h * dydx(xi + t2*h, yi + w2*k2)
        k4 = h * dydx(xi + t3*h, yi + w3*k3)

        # xi and yi at every iteration starting from (x1,y1) to (xn, yn)
        xi = xi + h
        yi = yi + (k1*c1 + k2*c2 + k3*c3 + k4*c4)

        MRK_domain.append(yi)

    # value at last value of domain i.e. at x = 10, calculate y [ y(10) = ? -> Generally it is required ]
    return MRK_domain[-1], MRK_domain

def rms(values_SRK, values_MRK):
    sum_square_diff = 0
    for i in range(len(values_SRK)):
        square_diff = (values_SRK[i] - values_MRK[i] )**2
        sum_square_diff += square_diff
    
    rms = sqrt( sum_square_diff/len(values_SRK) )
    return rms


# for making graph

step_size = []
rms_at_defined_h = []



# Scalability at defined h by h, h/2, h/4, h/8, h/16
divisor = [1, 2, 4, 8, 16]
for i in range(5):

    print(f'This iteration is for h/{divisor[i]}')
    h = h /divisor[i]
    n = int ( (xn - x) / h )
    
    # main
    result_SRK = S4RK()
    result_MRK = M4RK()
    result_rms = rms(result_SRK[1], result_MRK[1])
    print("Standard Runge Kutta= ", result_SRK[0])
    print("Modified Runge Kutta= ", result_MRK[0])
    print("Root Mean Square = ", result_rms)

    # appending
    step_size.append(h)
    rms_at_defined_h.append(result_rms)



# plot figure
plt.plot(step_size, rms_at_defined_h)
plt.xlabel("Step Size( h )")
plt.ylabel("Root Mean Square Error")
plt.show()
