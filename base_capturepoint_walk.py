#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np              # Numerical library
from scipy import *             # Load the scipy functions
from control.matlab import *    # Load the controls systems library
from matplotlib import pyplot as plt

Zc = 0.3
g = 9.8
foot_list = [[0.0, 0.06],[0.1, -0.06],[0.2, 0.06],[0.3, -0.06],[0.4, 0.06],[0.5, -0.06], [0.5, 0.0], [0.5, 0.0]]


def drange(begin, end, step):
    n = begin
    end = end + step
    while n < end:
     yield n
     n += step

class CapturePointWalk():
    dt_list = []
    xc_list = []
    yc_list = []
    cpx_list = []
    cpy_list = []
    ref_cpx_list = []
    ref_cpy_list = []
    px_list = []
    py_list = []
    t = []

    count = 0
    def __init__(self, period, dt, foot, k):
        self.__period = period
        self.__count  = 0
        self.__dt = dt
        self.__foot_list = foot
        self.__k = k
        self.__w = sqrt(g/Zc)

        self.__xc, self.__yc = 0.0, 0.0
        self.__dxc, self.__dyc = 0.0, 0.0

        self.__dxci, self.__dyci = 0.0, 0.0
        self.__xci, self.__yci = 0.0, 0.0
        self.__px, self.__py = 0.0, 0.0

        #euler ver.
        self.__cpx = self.__xc + self.__dxc/self.__w
        self.__cpy = self.__yc + self.__dyc/self.__w

        self.__ddxc = ((self.__xc - self.__cpx * self.__k)*g)/Zc
        self.__ddyc = ((self.__yc - self.__cpy * self.__k)*g)/Zc

    def set_footstep(self):
        for step in range(len(self.__foot_list)):
            for ttt in drange(0.0, self.__period-self.__dt, self.__dt):
                dt_n = float(format(round(self.__period - ttt, 3)))
                self.dt_list.append(dt_n)
                self.ref_cpx_list.append(self.__foot_list[step][0])
                self.ref_cpy_list.append(self.__foot_list[step][1])
                self.__count = self.__count + 1

    def calc_cog_trajectory_rungekutta(self):
        count = 0
        for num in drange(0,self.__count-1, 1):
            time = float(format(round( num*self.__dt, 3)))

            b = exp(self.__w*self.dt_list[count])

            # RungeKutta
            k1x = self.__dt * self.__dxci
            l1x = self.__dt *((self.__xc - self.__px)*g/Zc)
            k1y = self.__dt * self.__dyci
            l1y = self.__dt *((self.__yc - self.__py)*g/Zc)

            k2x = self.__dt * (self.__dxci + l1x/2)
            l2x = self.__dt *((self.__xc - self.__px)*g/Zc)
            k2y = self.__dt * (self.__dyci + l1y/2)
            l2y = self.__dt *((self.__yc - self.__py)*g/Zc)

            k3x = self.__dt * (self.__dxci + l2x/2)
            l3x = self.__dt *((self.__xc - self.__px)*g/Zc)
            k3y = self.__dt * (self.__dyci + l2y/2)
            l3y = self.__dt *((self.__yc - self.__py)*g/Zc)

            k4x = self.__dt * (self.__dxci + l3x)
            l4x = self.__dt *((self.__xc - self.__px)*g/Zc)
            k4y = self.__dt * (self.__dyci + l3y)
            l4y = self.__dt *((self.__yc - self.__py)*g/Zc)

            self.__xci = self.__xci + (k1x + 2*k2x + 2*k3x + k4x)/6
            self.__yci = self.__yci + (k1y + 2*k2y + 2*k3y + k4y)/6

            self.__dxci = self.__dxci + (l1x + 2*l2x + 2*l3x + l4x)/6
            self.__dyci = self.__dyci + (l1y + 2*l2y + 2*l3y + l4y)/6

            # calc capture point
            self.__cpx = self.__xci + self.__dxci/self.__w
            self.__cpy = self.__yci + self.__dyci/self.__w

            #calc zmp
            self.__p_xi  = (1/(1-b))*self.ref_cpx_list[count] - ((b/(1-b))*self.__cpx)
            self.__p_yi  = (1/(1-b))*self.ref_cpy_list[count] - ((b/(1-b))*self.__cpy)

            #append element
            self.xc_list.append(self.__xci)
            self.yc_list.append(self.__yci)
            self.cpx_list.append(self.__cpx)
            self.cpy_list.append(self.__cpy)
            self.px_list.append(self.__p_xi)
            self.py_list.append(self.__p_yi)
            self.t.append(time)

            # update
            self.__xc = self.__xci
            self.__yc = self.__yci
            self.__dxc = self.__dxci
            self.__dyc = self.__dyci
            self.__px = self.__p_xi
            self.__py = self.__p_yi
            count = count + 1

    def calc_cog_trajectory_euler(self):
        count = 0
        for num in drange(0,self.__count-1, 1):
            time = float(format(round( num*self.__dt, 3)))

            b = exp(self.__w*self.dt_list[count])
            dxci  = self.__dxc + self.__ddxc * self.__dt
            dyci  = self.__dyc + self.__ddyc * self.__dt

            xci   = self.__xc  + dxci * self.__dt + self.__ddxc * self.__dt * self.__dt /2
            yci   = self.__yc  + dyci * self.__dt + self.__ddyc * self.__dt * self.__dt /2

            if(time != 0):
                self.__cpx = xci + dxci/self.__w
                self.__cpy = yci + dyci/self.__w

            p_xi  = (1/(1-b))*self.ref_cpx_list[count] - ((b/(1-b))*self.__cpx)
            p_yi  = (1/(1-b))*self.ref_cpy_list[count] - ((b/(1-b))*self.__cpy)

            if(time == 0):
                ddxci = self.__ddxc
                ddyci = self.__ddyc
            else:
                ddxci = ((xci - p_xi)*g)/Zc
                ddyci = ((yci - p_yi)*g)/Zc

            self.xc_list.append(xci)
            self.yc_list.append(yci)
            self.cpx_list.append(self.__cpx)
            self.cpy_list.append(self.__cpy)
            self.px_list.append(p_xi)
            self.py_list.append(p_yi)
            self.t.append(time)

            self.__xc = xci
            self.__yc = yci
            self.__dxc = dxci
            self.__dyc = dyci
            self.__ddxc = ddxci
            self.__ddyc = ddyci

            count = count + 1

    def plot_gait_pattern_list(self):

        c = plt.subplot(3,1,1)
        plt.xlim(-0.1,0.6)
        plt.ylim(-0.10,0.10)
        c.plot(self.xc_list,self.yc_list, color = "red", label="$COM $")
        c.plot(self.cpx_list,self.cpy_list, color = "green", label="$CP $")
        c.plot(self.ref_cpx_list, self.ref_cpy_list, color = "lime", label="$ref CP $")
        c.plot(self.px_list,self.py_list,"*", label="$ZMP $")
        c.set_aspect('equal')
        c.legend(bbox_to_anchor=(1.0, 1.0), loc='upper right', borderaxespad=0, fontsize=5)

        tx = plt.subplot(3,1,2)
        tx.plot(self.t,self.xc_list, color = "red", label="$COMX$")
        tx.plot(self.t,self.cpx_list, color = "green", label="$CPX$")
        tx.plot(self.t,self.ref_cpx_list, color = "lime", label="$ref CPX$")
        tx.plot(self.t,self.px_list, color = "blue", label="$ZMPX$")
        tx.legend(bbox_to_anchor=(0, 1), loc='upper left', borderaxespad=0, fontsize=5)

        ty = plt.subplot(3,1,3)
        ty.plot(self.t,self.yc_list, color = "red", label="$COMY $")
        ty.plot(self.t,self.cpy_list, color = "green", label="$CPY$")
        ty.plot(self.t,self.ref_cpy_list, color = "lime", label="$ref CPY$")
        ty.plot(self.t,self.py_list, color = "blue", label="$ZMPY$")
        ty.legend(bbox_to_anchor=(0, 1), loc='upper left', borderaxespad=0, fontsize=5)
        plt.savefig('capturepoint_walk.png')
        plt.show()


def main():
    CP = CapturePointWalk(0.4,0.01,foot_list,1)
    CP.set_footstep()
    CP.calc_cog_trajectory_rungekutta()
    #CP.calc_cog_trajectory_euler()
    CP.plot_gait_pattern_list()

if __name__ == '__main__':
    main()
