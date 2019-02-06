#include "FootPlanner.h"
#include "CapturePointControl.h"
#include <string>
#include <sstream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>

// Change for your robot and walking cycle
static const double FOOT_WIDTH			= 0.050;
static const double HEIGHT_ZC			= 0.27;
static const double ZMP_RANGE			= 0.060;
static const double WALKING_HALF_CYCLE	= 0.340;
static const double SAMPLING_TIME		= 0.010;
static const double MAX_X_STEP			= 0.060;
static const double MAX_Y_STEP			= 0.060;
static const double MAX_W_STEP          = 0.200;

double rad2deg(double radian){
    return radian * 180/M_PI;
}

double deg2rad(int deg){
    return deg * M_PI / 180;
}

int main(int argc, char *argv[])
{
#if 0
    if (argc != 10 && argc != 11) {
		std::cout << "Generate a walking pattern." << std::endl <<
			"Usage: GenerateWalkPattern <target_x> <target_y> <target_th> <initial_speed_x> <initial_speed_y> " <<
			"<initial_cog_x> <initial_cog_y> <support leg 0:right 1:left> <0:start 1:walking 2:stop> [<filename>]" << std::endl <<
			"[Caution] If you do not input the <filename>, " << std::endl <<
			"1) a walking pattern for a half of the walking cycle is calculated, and" << std::endl <<
			"2) the cofficients for the Pre-Calculated Preview Control is displayed." << std::endl;
		return -1;
	}
	double target_x = atof(argv[1]);
	double target_y = atof(argv[2]);
    double target_th = deg2rad(atof(argv[3]));
	double speed_x = atof(argv[4]);
	double speed_y = atof(argv[5]);
	double cog_x = atof(argv[6]);
	double cog_y = atof(argv[7]);
	FootStatus foot_status = (atoi(argv[8]) == 1) ? RightLeg : LeftLeg;
	WalkingStatus walking_status = (atoi(argv[9]) == 1) ? Walking : Start;
	if (atoi(argv[9]) >= 2) walking_status = Stop;
	int stop_status = atoi(argv[9]) - 2;	//0:first half,1:second half
	//std::cout << "stop status:" << stop_status << std::endl;
	bool save_file = (argc == 11) ? true : false;
	std::ofstream writing_file;
	if (save_file) {
		writing_file.open(argv[10]);
		writing_file << "time(s) ref_x(m) zmp_x(m) cog_x(m) ref_y(m) zmp_y(m) cog_y(m) vel_x(m) vel_y(m) acc_x(m) acc_y(m)" << std::endl;
	}
	Vector2d com_pos, com_vel, com_acc;
	com_pos(0) = cog_x;
	com_pos(1) = cog_y;
	com_vel(0) = speed_x;
	com_vel(1) = speed_y;
	int sampling_num_half_cycle = (WALKING_HALF_CYCLE + SAMPLING_TIME / 2) / SAMPLING_TIME;

    FootPlanner plan_node(SAMPLING_TIME, WALKING_HALF_CYCLE);
    plan_node.SetFootStepParameter(MAX_X_STEP, MAX_Y_STEP, MAX_W_STEP, FOOT_WIDTH, WALKING_HALF_CYCLE);
    plan_node.SetTargetPos(target_x, target_y, target_th, foot_status, walking_status);

    #if 0
    	std::cout << "Results of foot steps" << std::endl;
    	for (int i = 0; i < plan_node.foot_step_list.size(); i ++)
    		std::cout << plan_node.foot_step_list[i][0] << ","
            << plan_node.foot_step_list[i][1] << ", "
            << plan_node.foot_step_list[i][2] << ", "
            << plan_node.foot_step_list[i][3] * 180/3.14 << "(deg) "<< ", "

            << std::endl;
    #endif
#endif
    CapturePointControl capture_point_control(WALKING_HALF_CYCLE, SAMPLING_TIME);
    capture_point_control.set_footstep();
    capture_point_control.ref_cp_patterngenerator();
    capture_point_control.calc_cog_trajectory();
    capture_point_control.plot_gait_pattern_list();
    capture_point_control.plot_gait_pattern_list_x_time();
    capture_point_control.plot_gait_pattern_list_y_time();
    return 0;
}
