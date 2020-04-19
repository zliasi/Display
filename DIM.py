# import
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
import os.path
import matplotlib.ticker as mticker
import math

#keyword and orientation
keyword = 'VHF_trans'
orientation = 'ph'
NP = 'Ag'

# Title (if only "" then the name of the first file)
title = keyword+orientation+NP

# Filetype of plot
end = ".pdf"
title = title + end

# Save to file (y/n)
save_file = "n"

# distancer
if NP == 'test':
    dNP = 7.064328
    if orientation == 'cn':
        dNPtoM = 3.21
    else:
        dNPtoM = 2.86
if NP == 'Au':
    dNP = 7.064328
    if orientation == 'cn':
        dNPtoM = 3.21
    else:
        dNPtoM = 2.86
if NP == 'Ag':
    dNP = 7.076576
    if orientation == 'cn':
        dNPtoM = 3.27
    else:
        dNPtoM = 2.92
if NP == 'Cu':
    dNP = 6.260896
    if orientation == 'cn':
        dNPtoM = 2.95
    else:
        dNPtoM = 2.60
if keyword == 'DHA':
    if orientation == 'mem7':
        dM = 5.60040
    if orientation == 'nh2':
        dM = 3.25918
    if orientation == 'ph':
        dM = 5.85543
    if orientation == 'cn':
        dM = 2.70757
if keyword == 'VHF_cis':
    if orientation == 'mem7':
        dM = 5.67266
    if orientation == 'nh2':
        dM = 2.06182
    if orientation == 'ph':
        dM = 5.36829
    if orientation == 'cn':
        dM = 3.34602
if keyword == 'VHF_trans':
    if orientation == 'mem7':
        dM = 5.74505
    if orientation == 'nh2':
        dM = 2.99884
    if orientation == 'ph':
        dM = 4.59709
    if orientation == 'cn':
        dM = 4.43717

d = (dNP+dM+dNPtoM)/0.52917721092
print(d)
Rc = [d, 0, 0]


# konstanter
Rconstant = 100
R = np.sqrt(Rc[0]**2 + Rc[1]**2 + Rc[2]**2)
evac = 0.7957747151e-1
V = 18403.0
gX = 1.0/3.0
gY = 1.0/3.0
gZ = 1.0/3.0
# For guld
if NP == 'test':
    y = 1.0
    wp = 1.0
    wn = [1.0, 1.0, 1.0, 1.0, 1.0]
    Op = np.sqrt(1.0) * wp
    fn = [1.0, 1.0, 1.0, 1.0, 1.0]
    Tn = [1.0, 1.0, 1.0, 1.0, 1.0]
if NP == 'Au':
    y = 0.053
    wp = 9.03
    wn = [0.415, 0.830, 2.969, 4.304, 13.32]
    Op = np.sqrt(0.760)*wp
    fn = [0.024, 0.010, 0.071, 0.601, 4.384]
    Tn = [0.241, 0.345, 0.870, 2.494, 2.214]
# For sølv
if NP == 'Ag':
    y = 0.048
    wp = 9.01
    wn = [0.816, 4.481, 8.185, 9.083, 20.29]
    Op = np.sqrt(0.845)*wp
    fn = [0.065, 0.124, 0.011, 0.840, 5.646]
    Tn = [3.886, 0.452, 0.065, 0.916, 2.419]

# For kobber
if NP == 'Cu':
    y = 0.030
    wp = 10.83
    wn = [0.291, 2.957, 5.300, 11.18]
    Op = np.sqrt(0.575)*wp
    fn = [0.061, 0.104, 0.723, 0.638]
    Tn = [0.378, 1.056, 3.213, 4.305]

print(y)
print(wp)
print(wn)

# start dipoler for forskellige orienteringer

DHAph = np.array([[0.632754], [2.19743], [-0.054577], [0], [0], [0]])
DHAnh2 = np.array([[-2.19743], [0.632754], [-0.054577], [0], [0], [0]])
DHAcn = np.array([[2.19743], [-0.632754], [-0.054577], [0], [0], [0]])
DHA7mem = np.array([[-0.632754], [-2.19743], [-0.054577], [0], [0], [0]])
VHFcph = np.array([[-0.20655], [-2.3245], [0.28496], [0], [0], [0]])
VHFc7mem = np.array([[0.20655], [2.3245], [0.28496], [0], [0], [0]])
VHFcnh2 = np.array([[-2.3245], [0.20655], [0.28496], [0], [0], [0]])
VHFccn = np.array([[2.3245], [-0.20655], [0.28496], [0], [0], [0]])
VHFt7mem = np.array([[-2.3424], [1.88709], [0.38788], [0], [0], [0]])
VHFtcn = np.array([[2.3424], [-1.88709], [0.38788], [0], [0], [0]])
VHFtph = np.array([[-1.88709], [2.3424], [0.38788], [0], [0], [0]])
VHFtnh2 = np.array([[1.88709], [-2.3424], [0.38788], [0], [0], [0]])

odipDHA = np.array([[0.632754], [2.19743], [-0.054577], [0], [0], [0]])
odipVHFc = np.array([[-0.20655], [-2.3245], [0.28496], [0], [0], [0]])
odipVHFt = np.array([[-2.3424], [1.88709], [0.38788], [0], [0], [0]])

# hiv frekvens og pol for molecule

w = []
aMXXI = []
aMYYI = []
aMZZI = []
aMXXR = []
aMYYR = []
aMZZR = []

for root, dirs, files in os.walk("/home/liasi/Documents/research_practice/imagpol"):
    for file in files:
        if file.endswith(keyword+'.out'):
            with open(os.path.join(root, file), 'r') as F:
                for line in F:
                    if 'XDIPLEN   XDIPLEN' in line:
                        XX = line.strip().split()
                        aMXXR.append(float(XX[4]))
                        aMXXI.append(float(XX[5]))
                        w.append(float(XX[3]))
                    if 'YDIPLEN   YDIPLEN' in line:
                        YY = line.strip().split()
                        aMYYR.append(float(YY[4]))
                        aMYYI.append(float(YY[5]))
                    if 'ZDIPLEN   ZDIPLEN' in line:
                        ZZ = line.strip().split()
                        aMZZR.append(float(ZZ[4]))
                        aMZZI.append(float(ZZ[5]))

w = np.array(w)
w * 27.21138602
w = np.array(w).tolist()
aMXX = []
aMYY = []
aMZZ = []
for i in range(len(aMXXI)):
    aMXX.append(complex(float(aMXXR[i]), float(aMXXI[i])))
    aMYY.append(complex(float(aMYYR[i]), float(aMYYI[i])))
    aMZZ.append(complex(float(aMZZR[i]), float(aMZZI[i])))

# rotation så dipolerne passer

print(w)


def rot(orientation, odipDHA, odipVHFc, odipVHFt, keyword, x, y, z):
    dip = []

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    if keyword == 'DHA':
        if orientation == 'mem7':
            dip.append(-(odipDHA[0]))
            dip.append(-(odipDHA[1]))
            dip.append(odipDHA[2])
            aMXX = -1.0*x
            aMYY = -1.0*y
            aMZZ = z
        if orientation == 'ph':
            dip.append(odipDHA[0])
            dip.append(odipDHA[1])
            dip.append(odipDHA[2])
            aMXX = 1.0*x
            aMYY = 1.0*y
            aMZZ = z

        if orientation == 'nh2':
            dip.append(-(odipDHA[1]))
            dip.append(odipDHA[0])
            dip.append(odipDHA[2])
            aMXX = -1.0*y
            aMYY = 1.0*x
            aMZZ = z

        if orientation == 'cn':
            dip.append(odipDHA[1])
            dip.append(-(odipDHA[0]))
            dip.append(odipDHA[2])
            aMXX = 1.0*y
            aMYY = -1.0*x
            aMZZ = z
    if keyword == 'VHF_cis':
        if orientation == 'mem7':
            dip.append(-(odipVHFc[0]))
            dip.append(-(odipVHFc[1]))
            dip.append(odipVHFc[2])
            aMXX = -1.0*x
            aMYY = -1.0*y
            aMZZ = z
        if orientation == 'ph':
            dip.append(odipVHFc[0])
            dip.append(odipVHFc[1])
            dip.append(odipVHFc[2])
            aMXX = 1.0*x
            aMYY = 1.0*y
            aMZZ = z

        if orientation == 'nh2':
            dip.append((odipVHFc[1]))
            dip.append(-(odipVHFc[0]))
            dip.append(odipVHFc[2])
            aMXX = 1.0*y
            aMYY = -1.0*x
            aMZZ = z

        if orientation == 'cn':
            dip.append(-(odipVHFc[1]))
            dip.append((odipVHFc[0]))
            dip.append(odipVHFc[2])
            aMXX = -1.0*y
            aMYY = 1.0*x
            aMZZ = z

    if keyword == 'VHF_trans':
        if orientation == 'mem7':
            dip.append((odipVHFt[0]))
            dip.append((odipVHFt[1]))
            dip.append(odipVHFt[2])
            aMXX = 1.0*x
            aMYY = 1.0*y
            aMZZ = z
        if orientation == 'ph':
            dip.append(-(odipVHFt[1]))
            dip.append((-odipVHFt[0]))
            dip.append((odipVHFt[2]))
            aMXX = -1.0*y
            aMYY = -1.0*x
            aMZZ = z

        if orientation == 'nh2':
            dip.append((odipVHFt[1]))
            dip.append((odipVHFt[0]))
            dip.append(odipVHFt[2])
            aMXX = 1.0*y
            aMYY = 1.0*x
            aMZZ = z

        if orientation == 'cn':
            dip.append(-(odipVHFt[0]))
            dip.append(-(odipVHFt[1]))
            dip.append(odipVHFt[2])
            aMXX = -1.0*x
            aMYY = -1.0*y
            aMZZ = z
    return dip, aMXX, aMYY, aMZZ


# den funktion vi brugte til NP polarizability
def dielec(w, y, wp, Op, fn, Tn, wn):
    C = 1.0 - Op ** 2.0 / (w * (complex(w, y)))
    c = 0

    for i in range(len(fn)):
        c += (fn[i] * wp ** 2.0) / \
            (complex((wn[i] ** 2.0 - w ** 2.0), -(w * Tn[i])))
    e = C + c
    #print("e:", e)
    return e


# lav elektriske felt
def elec(evac, dip, R, Rc):
    c = 1.0 / (4.0 * np.pi * evac)

    EMX = c * (dip[0] * (1 / (R ** 3) - 3.0 * (Rc[0] ** 2) / (R ** 5)) + dip[1] * (-3.0 * Rc[1] * Rc[0] / (R ** 5)) +
               dip[2] * (-3.0 * Rc[2] * Rc[0] / (R ** 5)))
    EMY = c * (dip[0] * (-3.0 * Rc[0] * Rc[1] / (R ** 5)) + dip[1] * (1 / (R ** 3) - 3.0 * (Rc[1] ** 2) / (R ** 5)) +
               dip[2] * (-3.0 * Rc[2] * Rc[1] / (R ** 5)))
    EMZ = c * (dip[0] * (-3.0 * Rc[0] * Rc[2] / (R ** 5)) + dip[1] * (-3.0 * Rc[1] * Rc[2] / (R ** 5)) + dip[2] * (
        1 / (R ** 3) - 3.0 * (Rc[2] ** 2) / (R ** 5)))
    ENPX = c * (dip[3] * (1 / (R ** 3) - 3.0 * (Rc[0] ** 2) / (R ** 5)) + dip[4] * (-3.0 * Rc[1] * Rc[0] / (R ** 5)) +
                dip[5] * (-3.0 * Rc[2] * Rc[0] / (R ** 5)))
    ENPY = c * (dip[3] * (-3.0 * Rc[0] * Rc[1] / (R ** 5)) + dip[4] * (1 / (R ** 3) - 3.0 * (Rc[1] ** 2) / (R ** 5)) +
                dip[5] * (-3.0 * Rc[2] * Rc[1] / (R ** 5)))
    ENPZ = c * (dip[3] * (-3.0 * Rc[0] * Rc[2] / (R ** 5)) + dip[4] * (-3.0 * Rc[1] * Rc[2] / (R ** 5)) + dip[5] * (
                1 / (R ** 3) - 3.0 * (Rc[2] ** 2) / (R ** 5)))
    Evec = np.array([EMX, EMY, EMZ, ENPX, ENPY, ENPZ])
    return Evec


wlist = w
# lav alpha af molecule


def alphaM(w):
    x = wlist.index(w)
    aM = []
    aM.append(aMXX[x])
    aM.append(aMYY[x])
    aM.append(aMZZ[x])
    aM = np.array(aM)
    return aM


# lav alpha af NP
def alphaNP(w, V, e, gX, gY, gZ):
    e = dielec(w, y, wp, Op, fn, Tn, wn)
    aNPXX = V * ((e - 1.0) / (1.0 + gX * (e - 1.0)))
    aNPYY = V * ((e - 1.0) / (1.0 + gY * (e - 1.0)))
    aNPZZ = V * ((e - 1.0) / (1.0 + gZ * (e - 1.0)))

    return aNPXX, aNPYY, aNPZZ

# lav T


def makeT(R, evac, Rconstant):
    T2XX = 1.0/(4.0*np.pi*evac)*(3*Rc[0]*Rc[0]/(R**5)-1/(R**3))*(math.erf(R/(np.sqrt(2)*Rconstant))-np.sqrt(2/np.pi)*R/Rconstant*math.exp(
        (-R**2)/(2*(Rconstant**2)))-np.sqrt(2/np.pi)*1/(Rconstant**3)*Rc[0]*Rc[0]/(R**2)*math.exp(-(R**2)/(2*(Rconstant**2))))
    T2YY = 1.0/(4.0*np.pi*evac)*(3*Rc[1]*Rc[1]/(R**5)-1/(R**3))*(math.erf(R/(np.sqrt(2)*Rconstant))-np.sqrt(2/np.pi)*R/Rconstant*math.exp(
        (-R**2)/(2*(Rconstant**2)))-np.sqrt(2/np.pi)*1/(Rconstant**3)*Rc[1]*Rc[1]/(R**2)*math.exp(-(R**2)/(2*(Rconstant**2))))
    T2ZZ = 1.0/(4.0*np.pi*evac)*(3*Rc[2]*Rc[2]/(R**5)-1/(R**3))*(math.erf(R/(np.sqrt(2)*Rconstant))-np.sqrt(2/np.pi)*R/Rconstant*math.exp(
        (-R**2)/(2*(Rconstant**2)))-np.sqrt(2/np.pi)*1/(Rconstant**3)*Rc[2]*Rc[2]/(R**2)*math.exp(-(R**2)/(2*(Rconstant**2))))

    #T2 = np.array([[T2XX,0,0],[0,T2YY,0],[0,0,T2ZZ]])
    return T2XX, T2YY, T2ZZ

# lav inv alpha


def makealpha(w):
    aM = alphaM(w)
    print("aM:", aM)
    a = np.array([[aM[0], 0, 0], [0, aM[1], 0], [0, 0, aM[2]]])
    aMinv = linalg.inv(a)
    e = dielec(w, y, wp, Op, fn, Tn, wn)
    aNPXX, aNPYY, aNPZZ = alphaNP(w, V, e, gX, gY, gZ)
    b = np.array([[aNPXX, 0, 0], [0, aNPYY, 0], [0, 0, aNPZZ]])
    aNPinv = linalg.inv(b)
    return aMinv, aNPinv

# lav A matricen


def makeA(w, R):
    e = dielec(w, y, wp, Op, fn, Tn, wn)
    aNPXX, aNPYY, aNPZZ = alphaNP(w, V, e, gX, gY, gZ)
    aM = alphaM(w)
    # print(aMXX)
    # print(aNPXX)
    T2XX, T2YY, T2ZZ = makeT(R, evac, Rconstant)
    A = np.array([[1.0/aM[0], 0, 0, -T2XX, 0, 0], [0, 1.0/aM[1], 0, 0, -T2YY, 0], [0, 0, 1.0/aM[2], 0, 0, -T2ZZ],
                  [-T2XX, 0, 0, 1.0/aNPXX, 0, 0], [0, -T2YY, 0, 0, 1.0/aNPYY, 0], [0, 0, -T2ZZ, 0, 0, 1.0/aNPZZ]])
    # T2 = makeT(R,evac,Rcoprint("aM:",aM)nstant)
    # print(aMinv.shape)
    # print(aNPinv.shape)
    # print(T2.shape)
    #aMinv = np.array(aMinv).tolist()
    #negT2 = np.array(-1.0*T2).tolist()
    #aNPinv = np.array(aNPinv).tolist()
    # np.matrix
    #A1 = aMinv+negT2
    #A2 = negT2+aNPinv
    #A = []
    # A.append(A1)
    # A.append(A2)
    # print(A)
    #A = np.array(A)
    # print(A)
    #A = np.array(A).tolist()

    # print(A.shape)
    # linalg.inv(A)
    B = A

    return A, B


dipXr = []
dipYr = []
dipZr = []
dipXi = []
dipYi = []
dipZi = []
dip, aMXX, aMYY, aMZZ = rot(
    orientation, odipDHA, odipVHFc, odipVHFt, keyword, aMXX, aMYY, aMZZ)
dip.append(0)
dip.append(0)
dip.append(0)

print("w:", w)
E = elec(evac, dip, R, Rc)
for i in w:

    A, B = makeA(i, R)
    print("SHAPE:", B.shape)
    B = linalg.inv(A)
    dipo = B.dot(E)
    print("dipo:", dipo)
    dipXr.append(np.real(dipo[0]))
    dipYr.append(np.real(dipo[1]))
    dipZr.append(np.real(dipo[2]))
    dipXi.append(np.imag(dipo[0]))
    dipYi.append(np.imag(dipo[1]))
    dipZi.append(np.imag(dipo[2]))

# print(dipXr)
# print(dipYr)
# print(dipZr)
# print(dipXi)
# print(dipYi)
# print(dipZi)
dipXr = np.array(dipXr)
dipYr = np.array(dipYr)
dipZr = np.array(dipZr)
dipXi = np.array(dipXi)
dipYi = np.array(dipYi)
dipZi = np.array(dipZi)

print(dipXr)


yx = dipXi/dipXr

yy = dipYi/dipYr

yz = dipZi/dipZr
# print(dipXi)

# ---------------plots------------------------

fig, (ax1, ax2, ax3) = plt.subplots(
    3, 1, sharex=True, sharey=False, figsize=(10, 6))

ax1.plot(w, dipXi, 'b', label='Imaginary')
ax1.set_title('X component', fontsize=10)
ax2.plot(w, dipYi, 'b')
ax2.set_title('Y component', fontsize=10)
ax3.plot(w, dipZi, 'b')
ax3.set_title('Z component', fontsize=10)
ax1.plot(w, dipXr, 'r', label='Real')
f = mticker.ScalarFormatter(useOffset=False, useMathText=True)


def g(w, pos): return "${}$".format(f._formatSciNotation('%1.10e' % w))


ax1.yaxis.set_major_formatter(mticker.FuncFormatter(g))
ax2.plot(w, dipYr, 'r')
f = mticker.ScalarFormatter(useOffset=False, useMathText=True)


def g(w, pos): return "${}$".format(f._formatSciNotation('%1.10e' % w))


ax2.yaxis.set_major_formatter(mticker.FuncFormatter(g))
ax3.plot(w, dipZr, 'r')
# ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%1.e'))
# ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%1.e'))
# ax3.yaxis.set_major_formatter(mtick.FormatStrFormatter('%1.e'))
f = mticker.ScalarFormatter(useOffset=False, useMathText=True)


def g(w, pos): return "${}$".format(f._formatSciNotation('%1.10e' % w))


ax3.yaxis.set_major_formatter(mticker.FuncFormatter(g))
ax3.set_xlabel(r'Frequency [a.u.]', fontsize=15)
ax2.set_ylabel(r'Dipole [a.u.]', fontsize=15)
ax1.legend(loc=0, fontsize=10)
if save_file.lower() == "y" or save_file.lower() == "yes":
    plt.savefig(title, bbox_inches='tight')
plt.show()
# print(w)
