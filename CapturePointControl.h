#ifndef _CAPTUREPOINTCONTROL_H_
#define _CAPTUREPOINTCONTROL_H_

#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;
using namespace Eigen;

static const float Zc = 0.3;
static const float ACCELALETION_GRAVITY = 9.81;

class CapturePointControl
{
    private:
        float period;
        float dt;
        int count;
        float x,y;
        float xi,yi;
		float dx,dy;
		float dxi,dyi;
		float px,py;
        float cpxi, cpyi;
        float pxi, pyi;
        float b;

        float kx1, kx2, kx3, kx4, ky1, ky2, ky3, ky4;
        float lx1, lx2, lx3, lx4, ly1, ly2, ly3, ly4;


    public:
        float omega;
        vector<Vector3f> foot_step_list;
        vector<float>ref_dt_list;
        vector<float> cog_list_x;
        vector<float> cog_list_y;
        vector<float> cp_list_x;
        vector<float> cp_list_y;
        vector<float> ref_cp_list_x;
        vector<float> ref_cp_list_y;
        vector<float> zmp_list_x;
        vector<float> zmp_list_y;
        vector<float> time_list;
    public:
        CapturePointControl(float _period, float _dt ) : period(_period), dt(_dt)
        {
            omega = sqrt(ACCELALETION_GRAVITY/Zc);
            /* STEP1*/
            x = 0.0;
            y = 0.0;
            dx = 0.0;
            dy = 0.0;
            xi = x; yi = y;
            dxi = dx; dyi = dy;
            px = 0.0; py = 0.0;
            pxi = px; pyi = py;
            count = 0;

        }
        void set_footstep();
        void ref_cp_patterngenerator();
        void RungeKutta();
        void capture_point_control(int time);
        void update();
        void calc_cog_trajectory();
        void plot_gait_pattern_list();
        void plot_gait_pattern_list_x_time();
        void plot_gait_pattern_list_y_time();
};
#endif
