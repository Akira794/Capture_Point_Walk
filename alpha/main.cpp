#include "CapturePointControl.h"
#include <string>
#include <sstream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>

int main(int argc, char *argv[])
{
    CapturePointControl capture_point_control(0.32, 0.01);
    capture_point_control.set_footstep();
    capture_point_control.ref_cp_petterngenerator();
    capture_point_control.calc_cog_trajectory();
    capture_point_control.plot_gait_pattern_list();
    //capture_point_control.plot_gait_pattern_list_x_time();
    //capture_point_control.plot_gait_pattern_list_y_time();
    return 0;
}
