#include "CapturePointControl.h"

void CapturePointControl::set_footstep(){
    this->foot_step_list.push_back(Vector3f(0.0, 0.06, 0.0));
	this->foot_step_list.push_back(Vector3f(0.06, -0.06, 0.0));
	this->foot_step_list.push_back(Vector3f(0.12, 0.06, 0.0));
	this->foot_step_list.push_back(Vector3f(0.18, -0.06, 0.0));
	this->foot_step_list.push_back(Vector3f(0.24, 0.06, 0.0));
	this->foot_step_list.push_back(Vector3f(0.30, -0.06, 0.0));
	this->foot_step_list.push_back(Vector3f(0.30, 0.0, 0.0));
    this->foot_step_list.push_back(Vector3f(0.30, 0.0, 0.0));

}

void CapturePointControl::ref_cp_patterngenerator(){
    size_t foot_step_count = 1;
	size_t size = foot_step_list.size();
    for(int step = 0; step < size; step++){
        for(double t = 0.0; t <= period-dt; t+= dt){
            ref_dt_list.push_back(period-t);
            ref_cp_list_x.push_back(foot_step_list[step][0]);
            ref_cp_list_y.push_back(foot_step_list[step][1]);
            count += 1;
        }
    }
}

void CapturePointControl::RungeKutta(){

    //TODO classにすべき
    kx1 = dt * dx;
    lx1 = dt * ( x - px ) * ACCELALETION_GRAVITY/Zc;
    ky1 = dt * dy;
    ly1 = dt * ( y - py ) * ACCELALETION_GRAVITY/Zc;

    kx2 = dt * (dx + (lx1)/2.0f);
    lx2 = dt * (( x - px ) * ACCELALETION_GRAVITY/Zc);
    ky2 = dt * (dy + (ly1)/2.0f);
    ly2 = dt * (( y - py ) * ACCELALETION_GRAVITY/Zc);

    kx3 = dt * (dx + (lx2)/2.0f);
    lx3 = dt * (( x - px ) * ACCELALETION_GRAVITY/Zc);
    ky3 = dt * (dy + (ly2)/2.0f);
    ly3 = dt * (( y - py ) * ACCELALETION_GRAVITY/Zc);

    kx4 = dt * (dx + lx3);
    lx4 = dt * (( x - px ) * ACCELALETION_GRAVITY/Zc);
    ky4 = dt * (dy + ly3);
    ly4 = dt * (( y - py ) * ACCELALETION_GRAVITY/Zc);

    xi = x + (kx1 + 2*kx2 + 2*kx3 + kx4)/6.0f;
    yi = y + (ky1 + 2*ky2 + 2*ky3 + ky4)/6.0f;

    dxi = dx + (lx1 + 2*lx2 + 2*lx3 + lx4)/6.0f;
    dyi = dy + (ly1 + 2*ly2 + 2*ly3 + ly4)/6.0f;

}

void CapturePointControl::capture_point_control(int time){
    //Capture Point 式(7)
    cpxi = xi + (dxi/omega);
    cpyi = yi + (dyi/omega);

    b = exp(omega*ref_dt_list[time]);
    //Capture Point Control 式(12)
    pxi = (1/(1-b))*ref_cp_list_x[time] - ((b/(1-b))*cpxi);
    pyi = (1/(1-b))*ref_cp_list_y[time] - ((b/(1-b))*cpyi);

    time_list.push_back(time);
}

void CapturePointControl::update(){
    cog_list_x.push_back(xi);
    cog_list_y.push_back(yi);
    cp_list_x.push_back(cpxi);
    cp_list_y.push_back(cpyi);
    zmp_list_x.push_back(pxi);
    zmp_list_y.push_back(pyi);

    //update status
    x = xi;
    y = yi;
    dx = dxi;
    dy = dyi;
    px = pxi;
    py = pyi;
}

void CapturePointControl::calc_cog_trajectory(){
    for(int t = 0; t< count; t++){
        RungeKutta();
        capture_point_control(t);
        update();
    }
}

void CapturePointControl::plot_gait_pattern_list(){
    FILE *gp = popen("gnuplot -persist\n", "w");

	fprintf(gp, "set xlabel \"x [m]\"\n");
	fprintf(gp, "set ylabel \"y [m]\"\n");
    fprintf(gp, "set size ratio -1\n");
    fprintf(gp, "plot '-' with lines lw 5 lt 7 title \"COM\", '-' with lines lw 2 lt 2 title \"Ref-CP\", '-' with lines lw 2 lt 4 title \"CP\", '-' with points lw 1 lt 6 title \"ZMP\",\n");
	for(std::size_t i=0;i<cog_list_x.size();i++) fprintf(gp, "%f\t%f\n", cog_list_x[i], cog_list_y[i]); fprintf(gp,"e\n");
    for(std::size_t i=0;i<ref_cp_list_x.size();i++) fprintf(gp, "%f\t%f\n", ref_cp_list_x[i], ref_cp_list_y[i]); fprintf(gp, "e\n");
    for(std::size_t i=0;i<cp_list_x.size();i++) fprintf(gp, "%f\t%f\n", cp_list_x[i], cp_list_y[i]); fprintf(gp, "e\n");
	for(std::size_t i=0;i<zmp_list_x.size();i++) fprintf(gp, "%f\t%f\n", zmp_list_x[i], zmp_list_y[i]); fprintf(gp, "e\n");
}

void CapturePointControl::plot_gait_pattern_list_x_time(){
    FILE *gp = popen("gnuplot -persist\n", "w");

    fprintf(gp, "set xlabel \"x [m]\"\n");
    fprintf(gp, "set ylabel \"y [m]\"\n");
    fprintf(gp, "set size ratio -1\n");
    fprintf(gp, "plot '-' with lines lw 3 lt 7 title \"COM\", '-' with lines lw 3 lt 2 title \"Ref-CP\", '-' with lines lw 3 lt 4 title \"CP\", '-' with lines lw 3 lt 6 title \"ZMP\",\n");
    for(std::size_t i=0;i<time_list.size();i++) fprintf(gp, "%f\t%f\n", time_list[i], cog_list_x[i]); fprintf(gp,"e\n");
    for(std::size_t i=0;i<time_list.size();i++) fprintf(gp, "%f\t%f\n", time_list[i], ref_cp_list_x[i]); fprintf(gp, "e\n");
    for(std::size_t i=0;i<time_list.size();i++) fprintf(gp, "%f\t%f\n", time_list[i], cp_list_x[i]); fprintf(gp, "e\n");
    for(std::size_t i=0;i<time_list.size();i++) fprintf(gp, "%f\t%f\n", time_list[i], zmp_list_x[i]); fprintf(gp, "e\n");
    fprintf(gp,"exit\n");
}

void CapturePointControl::plot_gait_pattern_list_y_time(){
    FILE *gp = popen("gnuplot -persist\n", "w");

    fprintf(gp, "set xlabel \"x [m]\"\n");
    fprintf(gp, "set ylabel \"y [m]\"\n");
    fprintf(gp, "set size ratio -1\n");
    fprintf(gp, "plot '-' with lines lw 3 lt 7 title \"COM\", '-' with lines lw 3 lt 2 title \"Ref-CP\", '-' with lines lw 3 lt 4 title \"CP\", '-' with lines lw 3 lt 6 title \"ZMP\",\n");
    for(std::size_t i=0;i<time_list.size();i++) fprintf(gp, "%f\t%f\n", time_list[i], cog_list_y[i]); fprintf(gp,"e\n");
    for(std::size_t i=0;i<time_list.size();i++) fprintf(gp, "%f\t%f\n", time_list[i], ref_cp_list_y[i]); fprintf(gp, "e\n");
    for(std::size_t i=0;i<time_list.size();i++) fprintf(gp, "%f\t%f\n", time_list[i], cp_list_y[i]); fprintf(gp, "e\n");
    for(std::size_t i=0;i<time_list.size();i++) fprintf(gp, "%f\t%f\n", time_list[i], zmp_list_y[i]); fprintf(gp, "e\n");
    fprintf(gp,"exit\n");
}
