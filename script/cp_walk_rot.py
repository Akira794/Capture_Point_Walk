#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np              # Numerical library
from scipy import *             # Load the scipy functions
from control.matlab import *    # Load the controls systems library
from matplotlib import pyplot as plt


def drange(begin, end, step):
    n = begin
    end = end + step
    while n < end:
     yield n
     n += step

Zc = 0.3
g = 9.8

dist = [0.6, 0.05, ( 3.14/2.0) ]
max_step_x = 0.1 #
max_step_y = 0.05 #
max_step_w = 0.2 #
foot_y = 0.06
period = 0.4 #
forward_period = 1.0 #
calculate_period = 4.0
dt = 0.01

prefx = []
prefy = []
foot_list = []

step_num_x = int(abs(dist[0])/max_step_x)
step_num_y = int(abs(dist[1])/max_step_y)
step_num_w = int(abs(dist[2])/max_step_w)
step_num = max([step_num_x, step_num_y, step_num_w])

if (step_num_x == step_num):
    step_x = sign(dist[0])*min(abs(dist[0]), max_step_x)
    step_y = dist[1]/abs(dist[0])*max_step_x
    step_w = dist[2]/abs(dist[0])*max_step_x
elif (step_num_y == step_num):
    step_x = dist[0]/abs(dist[1])*max_step_y
    step_y = sign(dist[1])*min(abs(dist[1]), max_step_y)
    step_w = dist[2]/abs(dist[1])*max_step_y
elif (step_num_w == step_num):
    step_x = dist[0]/abs(dist[2])*max_step_w
    step_y = dist[1]/abs(dist[2])*max_step_w
    step_w = sign(dist[2])*min(abs(dist[2]), max_step_w)

step = [0.0, 0.0]
foot = [step]
foot.append([0.0, foot_y])
rot = 0

for i in drange(1, step_num, 1):
    shift_y = foot_y * ((divmod(i,2)[1])*2-1)*2*(-1)
    foot_rd = [step_x, step_y+shift_y]
    foot_fd = [foot_rd[0]*cos(rot)-foot_rd[1]*sin(rot), foot_rd[0]*sin(rot)+foot_rd[1]*cos(rot)]
    rot = rot + step_w
    foot.append([foot[i][0] + foot_fd[0], foot[i][1] + foot_fd[1]])

foot_rd = [dist[0]-step_x*step_num, dist[1]-step_y*step_num-shift_y]
foot_fd = [foot_rd[0]*cos(rot)-foot_rd[1]*sin(rot), foot_rd[0]*sin(rot)+foot_rd[1]*cos(rot)]
foot.append([foot[i+1][0] + foot_fd[0], foot[i+1][1] + foot_fd[1]])

rot = dist[2]
foot_rd = [0,shift_y/2]
foot_fd = [foot_rd[0]*cos(rot)-foot_rd[1]*sin(rot), foot_rd[0]*sin(rot)+foot_rd[1]*cos(rot)]
foot.append([foot[i+2][0] + foot_fd[0], foot[i+2][1] + foot_fd[1]])
foot.append([foot[i+2][0] + foot_fd[0], foot[i+2][1] + foot_fd[1]])

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
    regionx_list = []
    regiony_list = []

    count = 0
    def __init__(self, period, dt, foot):
        self.__period = period
        self.__count  = 0
        self.__dt = dt
        self.__foot_list = foot
        self.__w = sqrt(g/Zc)

        self.__xc, self.__yc = 0.0, 0.0
        self.__dxc, self.__dyc = 0.0, 0.0

        self.__dxci, self.__dyci = 0.0, 0.0
        self.__xci, self.__yci = 0.0, 0.0
        self.__px, self.__py = 0.0, 0.0

    def set_footstep(self):
        for step in range(len(self.__foot_list)):
            for ttt in drange(0.0, self.__period-self.__dt, self.__dt):
                dt_n = float(format(round(self.__period - ttt, 3)))
                self.dt_list.append(dt_n)
                self.ref_cpx_list.append(self.__foot_list[step][0])
                self.ref_cpy_list.append(self.__foot_list[step][1])
                self.__count = self.__count + 1

    def calc_cog_trajectory(self):
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

    def plot_gait_pattern(self):

        plt.xlim(-0.1,0.5)
        plt.ylim(-0.10,0.5)
        plt.axes().set_aspect('equal')
        plt.plot(self.xc_list,self.yc_list, color = "red", label="$COM $")
        plt.plot(self.cpx_list,self.cpy_list, color = "green", label="$CP $")
        plt.plot(self.ref_cpx_list, self.ref_cpy_list, color = "lime", label="$ref CP $")
        plt.plot(self.px_list,self.py_list,"*", label="$ZMP $")
        plt.legend(bbox_to_anchor=(0.0, 1.0), loc='upper left', borderaxespad=0, fontsize=10)
        plt.savefig('capturepoint_walk_rot.png')
        plt.show()

    def plot_gait_pattern_list(self):

        tx = plt.subplot(2,1,1)
        tx.plot(self.t,self.xc_list, color = "red", label="$COMX$")
        tx.plot(self.t,self.cpx_list, color = "green", label="$CPX$")
        tx.plot(self.t,self.ref_cpx_list, color = "lime", label="$ref CPX$")
        tx.plot(self.t,self.px_list, color = "blue", label="$ZMPX$")
        tx.legend(bbox_to_anchor=(0, 1), loc='upper left', borderaxespad=0, fontsize=10)

        ty = plt.subplot(2,1,2)
        ty.plot(self.t,self.yc_list, color = "red", label="$COMY $")
        ty.plot(self.t,self.cpy_list, color = "green", label="$CPY$")
        ty.plot(self.t,self.ref_cpy_list, color = "lime", label="$ref CPY$")
        ty.plot(self.t,self.py_list, color = "blue", label="$ZMPY$")
        ty.legend(bbox_to_anchor=(0, 1), loc='upper left', borderaxespad=0, fontsize=10)
        plt.savefig('capturepoint_walk_tx_ty.png')
        plt.show()

def main():
    CP = CapturePointWalk(0.4,0.01,foot)
    CP.set_footstep()
    CP.calc_cog_trajectory()
    CP.plot_gait_pattern()
    CP.plot_gait_pattern_list()

if __name__ == '__main__':
    main()
