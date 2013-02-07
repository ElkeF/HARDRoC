#! /opt/local/bin/python

# Purpose to fit a number of input parameters to an exponential
# function of the form: f(x) = factor * exp(-alpha*x) + const

import numpy as np
from scipy import *
from scipy.optimize import leastsq

def residuals(p, y, x):
        err = y-peval(x,p)
        return err

def peval(x, p):
        return p[0]*(exp(-p[1]*x)) + p[2]

filename=('ArXe.txt')
data = np.loadtxt(filename)

y = data[:,1]
x = data[:,0]

#y = y.astype(np.float64)
#x = x.astype(np.float64)

factor = 1E-4
alpha  = 0.1
const  = 0
pname = (['factor','alpha','const'])
p0 = array([factor , alpha, const])
plsq = leastsq(residuals, p0, args=(y, x), maxfev=2000, full_output=1)

print "Final parameters"
for i in range(len(pname)):
        print "%s = %.4f " % (pname[i], p0[i])


#import numpy as np
#from scipy.optimize import curve_fit

#fitfunc = lambda p, x: p[0]*np.exp(-p[1]*x) + p[2]


# Here comes the input data
#x  = [4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.6,6.0,6.5,7.0,7.5,8.0,8.5,9.0]
#y  = [0.00528753,0.00380589,0.002719340,0.001930730,0.001363320,0.000957980,\
#      0.000670090,0.000323350,0.000153080,0.000058960,0.000022870,0.000009290,\
#      0.000003990,0.000001760,0.000000770]
#yn = y + 0.2*np.random.normal(size=len(x))

#popt, pcov = curve_fit(func,x,y)
