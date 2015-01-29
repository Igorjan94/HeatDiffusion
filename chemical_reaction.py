from math import exp

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from tkinter import *
from matplotlib import cm

# equations:    dX/dt -     D d2X/dz2 =      w(X, T)
#               dT/dt - kappa d2T/dz2 = -Q/C w(X, T)
#               w(X, T) = -K X^alpha exp(-E/RT)

# Physical parameters

R = 8.314
E = 8e+4
K = 1.6e+6
alpha = 1.0
Q = 7e+5
C = 1980.0
rho = 830.0
T0 = 300.0
Tm = T0 + Q / C
lam = 0.13
kappa = lam / (rho * C)
D = kappa
X_left = 0.0
X_right = 1.0
# Computational parameters

method = 'explicit'
dt = 1.0
dz = 0.001
tau = 2000
l = 0.05
time_steps = int(tau / dt)
space_steps = int(l / dz)
rate = 1


def w(X, T):
    return -K * X ** alpha * exp(-E / (R * T))


def solve_tridiagonal(a, b, c, d):
    n = len(d)
    for i in range(n - 1):
        d[i + 1] -= d[i] * a[i] / b[i]
        b[i + 1] -= c[i] * a[i] / b[i]
    for i in reversed(range(n - 1)):
        d[i] -= d[i + 1] * c[i] / b[i + 1]
    return [d[i] / b[i] for i in range(n)]


def initial_X(space_steps):
    return [X_left] + [X_right] * (space_steps - 1)


def initial_T(space_steps):
    return [Tm] + [T0] * (space_steps - 1)


def generic_scheme(method_step):
    result_X, result_T = [initial_X(space_steps)], [initial_T(space_steps)]
    for i in range(time_steps):
        if (i % 100 == 0):
            print("%i/%i iterations finished" % (i, time_steps))
        new_X, new_T = method_step(result_X[-1], result_T[-1])
        result_X.append(new_X)
        result_T.append(new_T)
    return result_X, result_T


def explicit_step(prev_X, prev_T):
    new_X, new_T = [X_left], [Tm]
    for i in range(1, space_steps - 1):
        new_X.append(dt * (D * (prev_X[i - 1] - 2 * prev_X[i] + prev_X[i + 1]) / dz ** 2
                           + w(prev_X[i], prev_T[i])) + prev_X[i])
        new_T.append(dt * (kappa * (prev_T[i - 1] - 2 * prev_T[i] + prev_T[i + 1]) / dz ** 2
                           + w(prev_X[i], prev_T[i])) + prev_T[i])
    new_X.append(X_right)
    new_T.append(T0)
    return new_X, new_T


def implicit_step(prev_X, prev_T):
    n = len(prev_X)

    aX = [-D / dz ** 2 for _ in range(n - 2)] + [0.0]
    bX = [1.0] + [1.0 / dt + 2 * D / dz ** 2 for _ in range(n - 2)] + [1.0]
    cX = [0.0] + [-D / dz ** 2 for _ in range(n - 2)]
    dX = [X_left] + [prev_X[i] / dt + w(prev_X[i], prev_T[i]) for i in range(1, n - 1)] + [X_right]

    aT = [- kappa / dz ** 2 for _ in range(n - 2)] + [0.0]
    bT = [1.0] + [1.0 / dt + 2 * kappa / dz ** 2 for _ in range(n - 2)] + [1.0]
    cT = [0.0] + [-kappa / dz ** 2 for _ in range(n - 2)]
    dT = [Tm] + [prev_T[i] / dt - Q / C * w(prev_X[i], prev_T[i]) for i in range(1, n - 1)] + [T0]

    return solve_tridiagonal(aX, bX, cX, dX), solve_tridiagonal(aT, bT, cT, dT)


def euler_backward_diag_step(prev_X, prev_T):
    n = len(prev_X)

    aX = [-D / dz ** 2 for _ in range(n - 2)] + [0.0]
    bX = [1.0] + [1.0 / dt + 2 * D / dz ** 2 for _ in range(n - 2)] + [1.0]
    cX = [0.0] + [-D / dz ** 2 for _ in range(n - 2)]
    dX = [X_left] + [prev_X[i] / dt + w(prev_X[i], prev_T[i]) for i in range(1, n - 1)] + [X_right]
    ans_X = solve_tridiagonal(aX, bX, cX, dX)

    aT = [- kappa / dz ** 2 for _ in range(n - 2)] + [0.0]
    bT = [1.0] + [1.0 / dt + 2 * kappa / dz ** 2 for _ in range(n - 2)] + [1.0]
    cT = [0.0] + [-kappa / dz ** 2 for _ in range(n - 2)]
    dT = [Tm] + [prev_T[i] / dt - Q / C * w(ans_X[i], prev_T[i]) for i in range(1, n - 1)] + [T0]

    return ans_X, solve_tridiagonal(aT, bT, cT, dT)


def backward_iterations_step(prev_X, prev_T, iterations=5):
    n = len(prev_X)
    ans_X, ans_T = prev_X[:], prev_T[:]
    for _ in range(iterations):
        aX = [-D / dz ** 2 for _ in range(n - 2)] + [0.0]
        bX = [1.0] + [1.0 / dt + 2 * D / dz ** 2 for _ in range(n - 2)] + [1.0]
        cX = [0.0] + [-D / dz ** 2 for _ in range(n - 2)]
        dX = [X_left] + [prev_X[i] / dt + w(ans_X[i], ans_T[i]) for i in range(1, n - 1)] + [X_right]

        aT = [- kappa / dz ** 2 for _ in range(n - 2)] + [0.0]
        bT = [1.0] + [1.0 / dt + 2 * kappa / dz ** 2 for _ in range(n - 2)] + [1.0]
        cT = [0.0] + [-kappa / dz ** 2 for _ in range(n - 2)]
        dT = [Tm] + [prev_T[i] / dt - Q / C * w(ans_X[i], ans_T[i]) for i in range(1, n - 1)] + [T0]

        ans_X, ans_T = solve_tridiagonal(aX, bX, cX, dX), solve_tridiagonal(aT, bT, cT, dT)

    return ans_X, ans_T


def refresh_parameters():
    global R, E, K, alpha, Q, rho, T0, C, lam, D, kappa, time_steps, space_steps
    global dt, dz, time_steps, space_steps, rate, method

    e2f, e2i, e2str = lambda x: float(x.get()), lambda x: int(x.get()), lambda x: x.get()
    rate, = list(map(e2i, [erate]))
    dt, dz, R, E, K, alpha, Q, rho, T0, C, lam, D, tau, l = \
        list(map(e2f, [edt, edz, eR, eE, eK, ealpha, eQ, erho, eT0, eC, elam, eD, etau, el]))
    method, = list(map(e2str, [emethod]))

    kappa = lam / (rho * C)
    time_steps = int(tau / dt)
    space_steps = int(l / dz)


numerical_method_steps = {
    'explicit': explicit_step,
    'implicit': implicit_step,
    # euler_backward_step,
    # euler_backward_diag_step,
    # backward_iterations_step
}


def show_plots(event):
    refresh_parameters()

    if drawX.get() == 1:
        fig1 = plt.figure()
        xx = fig1.add_subplot(111, projection='3d')
        plt.subplots_adjust(left=0.03, right=0.97, top=0.97, bottom=0.03)
        xx.set_xlabel('Time')
        xx.set_ylabel('z')
        xx.set_zlabel('X')

    if drawT.get() == 1:
        fig2 = plt.figure()
        tx = fig2.add_subplot(111, projection='3d')
        plt.subplots_adjust(left=0.03, right=0.97, top=0.97, bottom=0.03)
        tx.set_xlabel('Time')
        tx.set_ylabel('z')
        tx.set_zlabel('T')

    if drawW.get() == 1:
        fig3 = plt.figure()
        wx = fig3.add_subplot(111, projection='3d')
        plt.subplots_adjust(left=0.03, right=0.97, top=0.97, bottom=0.03)
        wx.set_xlabel('Time')
        wx.set_ylabel('z')
        wx.set_zlabel('w')

    method_step = numerical_method_steps[method]
    X, T = generic_scheme(method_step)
    W = [[w(X[n][i], T[n][i]) for i in range(space_steps)] for n in range(time_steps)]
    print('drawing plots...')
    for i in range(time_steps):
        if i % rate == 0:
            time = i * dt

            xs = np.asarray([time] * space_steps)
            ys = np.linspace(0, space_steps * dz, num=space_steps)
            zxs = np.asarray(X[i])
            zts = np.asarray(T[i])
            zws = np.asarray(W[i])

            if drawX.get() == 1:
                xx.plot(xs, ys, zxs, color='r')
            if drawT.get() == 1:
                tx.plot(xs, ys, zts, color='r')
            if drawW.get() == 1:
                wx.plot(xs, ys, zws, color='r')

    if drawX.get() == 1:
        xx.invert_xaxis()
    if drawT.get() == 1:
        tx.invert_xaxis()
    if drawW.get() == 1:
        wx.invert_xaxis()

    plt.show()


def main():
    root = Tk()
    root.resizable(False, False)
    root.geometry('+200+200')

    global drawX, drawT, drawW
    global edt, edz, etau, el, erate, emethod
    global eR, eE, eK, ealpha, eQ, erho, eT0, eTm, eC, elam, eD

    comp = Frame(root)
    comp.pack(side='left')
    phys = Frame(root)
    phys.pack(side='left')
    pad = Frame(root)
    pad.pack(side='left', padx=10)

    lmethod = Label(comp, text='Method:')
    lmethod.pack(side='top', padx=5, pady=5)
    emethod = Entry(comp)
    emethod.insert(0, str(method))
    emethod.pack(side='top')

    ldt = Label(comp, text='dt:')
    ldt.pack(side='top', padx=5, pady=5)
    edt = Entry(comp)
    edt.insert(0, str(dt))
    edt.pack(side='top')

    ldz = Label(comp, text='dz:')
    ldz.pack(side='top', padx=5, pady=5)
    edz = Entry(comp)
    edz.insert(0, str(dz))
    edz.pack(side='top')

    ltau = Label(comp, text='τ:')
    ltau.pack(side='top', padx=5, pady=5)
    etau = Entry(comp)
    etau.insert(0, str(tau))
    etau.pack(side='top')

    ll = Label(comp, text='l:')
    ll.pack(side='top', padx=5, pady=5)
    el = Entry(comp)
    el.insert(0, str(l))
    el.pack(side='top')

    lrate = Label(comp, text='отображать каждый i-й шаг, i:')
    lrate.pack(side='top', padx=5, pady=5)
    erate = Entry(comp)
    erate.insert(0, str(rate))
    erate.pack(side='top')

    drawX = IntVar()
    drawX.set(1)
    cdrawX = Checkbutton(comp, text='график для X', variable=drawX)
    cdrawX.pack(sid='top', padx=5, pady=5)

    drawT = IntVar()
    drawT.set(1)
    cdrawT = Checkbutton(comp, text='график для T', variable=drawT)
    cdrawT.pack(sid='top', padx=5, pady=5)

    drawW = IntVar()
    drawW.set(1)
    cdrawW = Checkbutton(comp, text='график для w', variable=drawW)
    cdrawW.pack(sid='top', padx=5, pady=5)

    lR = Label(phys, text='R (менять с осторожностью!):')
    lR.pack(side='top', padx=5, pady=5)
    eR = Entry(phys)
    eR.insert(0, str(R))
    eR.pack(side='top')

    lE = Label(phys, text='E:')
    lE.pack(side='top', padx=5, pady=5)
    eE = Entry(phys)
    eE.insert(0, str(E))
    eE.pack(side='top')

    lK = Label(phys, text='K:')
    lK.pack(side='top', padx=5, pady=5)
    eK = Entry(phys)
    eK.insert(0, str(K))
    eK.pack(side='top')

    lalpha = Label(phys, text='α:')
    lalpha.pack(side='top', padx=5, pady=5)
    ealpha = Entry(phys)
    ealpha.insert(0, str(alpha))
    ealpha.pack(side='top')

    lQ = Label(phys, text='Q:')
    lQ.pack(side='top', padx=5, pady=5)
    eQ = Entry(phys)
    eQ.insert(0, str(Q))
    eQ.pack(side='top')

    lrho = Label(phys, text='ρ:')
    lrho.pack(side='top', padx=5, pady=5)
    erho = Entry(phys)
    erho.insert(0, str(rho))
    erho.pack(side='top')

    lT0 = Label(phys, text='T0:')
    lT0.pack(side='top', padx=5, pady=5)
    eT0 = Entry(phys)
    eT0.insert(0, str(T0))
    eT0.pack(side='top')

    lC = Label(phys, text='C:')
    lC.pack(side='top', padx=5, pady=5)
    eC = Entry(phys)
    eC.insert(0, str(C))
    eC.pack(side='top')

    llam = Label(phys, text='λ:')
    llam.pack(side='top', padx=5, pady=5)
    elam = Entry(phys)
    elam.insert(0, str(lam))
    elam.pack(side='top')

    lD = Label(phys, text='D:')
    lD.pack(side='top', padx=5, pady=5)
    eD = Entry(phys)
    eD.insert(0, str(D))
    eD.pack(side='top')

    plot = Button(phys, text='Draw')
    plot.bind("<Button-1>", show_plots)
    plot.pack(side='top', padx=10, pady=15)

    root.mainloop()


if __name__ == '__main__':
    main()
