import cmath
from cmath import phase
import numpy as np
import matplotlib.pyplot as plt
import sympy
from scipy import signal
from sympy import Symbol, simplify
from sympy.solvers import solve

print "\nElectrical Circuits\n"
vo = Symbol("vo")
vr = Symbol("vr")
ir = Symbol("ir")
ic = Symbol("ic")
il = Symbol("il")
r = Symbol("r")
omega = Symbol("omega")
c = Symbol("c")
l = Symbol("l")


eq1 = (vr + vo - 1, 
       ir - ic - il, 
       vr - ir*r,
       # what does 1j mean?  '''1j refers to the imaginary number i=\sqrt{-1}'''
       vo - ic/(1j*omega*c),
       # complete the following line:
       vo - il*(1j*omega*l))
# what does the following line do? 
'''It solves the systems of equations in eq1'''

sol = solve(eq1, (vo, vr, ir, ic, il))
vos = simplify(sol[vo])
vos1 = sol[vo]
# compare the output of the following line if vos = sol[vo]
'''There is no difference becuase sol[vo] is already simplified'''
print vos
print vos1

numvalue = {c: 10**-6, l: 10**-3}
# what does subs()?
'''It substitutes the value to corresponding variable'''

# is vos.subs(c=10**-6, l=10**-3) allowed? Try it. 
'''No since the arguments should be a dictionary'''

vosnum = vos.subs(numvalue)
flist = [vosnum.subs({r: 100.0*3**s}) for s in range(0, 4)]
omega_axis = np.linspace(20000, 43246, 100)
# what does 121 in the following line mean?
'''It creates the first cell in a 1 row by 2 column notional grid ''' 

# what are the other possible parameters of subplot()?
'''The other parameters of subplot() are axisbg, polar and projection'''

plt.figure(figsize=plt.figaspect(0.4))
plt.subplot(121)
# describe (python type, dimensions, etc) of the input parameter/s of zip() below
'''The input should be a tuple, list and array of the same dimension'''

# what does zip(*a) do if a is a 2-D list or numpy array?
'''It returns a list of tuples, where each tuple contains the i-th element from each row of the array a'''

plt.plot(omega_axis, zip(*[[abs(f.subs({omega: o})) for o in omega_axis] 
                                                    for f in flist]))
plt.xlim(20000, 43246)
plt.ylim(0, 1)
plt.xlabel('omega')
plt.ylabel('Abs[vo]')
plt.xticks([20000, 30000, 40000])

plt.subplot(122)
## Replicate Fig. 2.6, right pane following the code for Fig. 2.6, left pane
plt.plot(omega_axis, zip(*[[phase(f.subs({omega: o})) for o in omega_axis] 
                                                    for f in flist]))

plt.xlim(20000, 43246)
plt.ylim(-1.5, 1.5)
plt.xlabel('omega')
plt.ylabel('phase[vo]')
plt.xticks([20000, 30000, 40000])
plt.show()

def vsaw(t, T=1.0): 
    return signal.sawtooth(2*np.pi*t*T)
    '''The function creates a sawtooth voltage'''
    # complete this function   

omegares = 1./np.sqrt(np.prod(numvalue.values()))
alist = (1/np.sqrt(256)) * vsaw(np.arange(256)/256.0)
blist = np.sqrt(256) * np.fft.fft(alist)

def plot3(fac, w):
    # add a docstring for this function
    """
    The function plots the voltage at the output of the filter for various frequency and
    resistance values
    """
    
    omegai = fac * omegares
    # How were the limits of arange() in the following line chosen?
    '''The limits of arange() were chosen so that the size of volist
        is the same with blist. The use of a negative limit in 
        one of the arange() in volist is to shift higher frequencies to low 
        negative frequencies to avoid inaccurate approximation of Vi(t)'''
        
    volist = np.concatenate(([complex(vosnum.subs({omega: omegai*s, r:
                                                   w}).evalf()) 
                                 for s in np.arange(1, 129)],
                             [0.0],
                             [complex(vosnum.subs({omega: omegai*s, r:
                                                   w}).evalf()) 
                                 for s in np.arange(-127, 0)]))
    vtrans = np.fft.ifft(blist * volist)
    plotlist = np.array([[(k+1)/256., vtrans[k%256]] for k in range(768)])
    plt.plot(plotlist[:,0], plotlist[:,1])
    # what does the following line do?
    '''It adds a horizontal line at y=0'''
    
    plt.axhline(0)
    # add labels
    plt.xlabel('$t/T$')
    plt.ylabel('$V_{0}(t)$')
    plt.show()

plot3(1, 2700.0)
plot3(1/3., 200.0)
plot3(3.0, 5.0)

eq2 = (ir * (r + 1/(1j*omega*c) + 1j*omega*l) + vo - 1,
       ir - (1j*omega*c + 1/(1j*omega*l)) * vo)
sol2 = solve(eq2, (vo, ir))    # complete this line
vos2 = simplify(sol2[vo])
irs = simplify(sol2[ir])

# why should irs be passed to sympy.abs() before squaring?
'''Since irs contains imaginary numbers, we need to take the absolute value
first so that Python can undestand the command.'''
power = (r**2) *( sympy.Abs(irs)**2)
flist3 = [sympy.Abs(vos2.subs(numvalue).subs({r: 10.0*3**s})) 
            for s in range(0, 3)]
omega_axis = np.linspace(10000, 70000, 1000)
lines = plt.plot(omega_axis,zip(*[[abs(f.subs({omega: o})) for o in omega_axis]
                                                    for f in flist3])) # ...
# what does plt.setp() do?
'''It set a property on an artist object.'''
plt.setp(lines[0], lw=2)
plt.setp(lines[1], ls='--')
# add labels and ticks
plt.ylabel('$|V_0|$')
plt.xlabel('$\omega$')
plt.xticks([10000,30000,50000,70000])
plt.minorticks_on()
plt.show()

# replicate fig. 2.10
powerlist = sympy.Abs(power.subs(numvalue).subs({r:10.0}))
plt.plot(omega_axis, [powerlist.subs({omega: o}) for o in omega_axis])
plt.ylabel('$P/P_{0}$')
plt.xlabel('$\omega$')
plt.xticks([10000,30000,50000,70000])
plt.minorticks_on()
plt.show()