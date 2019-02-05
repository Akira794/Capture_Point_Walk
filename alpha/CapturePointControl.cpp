#include "CapturePointControl.h"

void CapturePointControl::set_footstep(){
    this->foot_step_list.push_back(Vector3f(0.0, 0.06, 0.0));
	this->foot_step_list.push_back(Vector3f(0.1, -0.06, 0.0));
	this->foot_step_list.push_back(Vector3f(0.2, 0.06, 0.0));
	this->foot_step_list.push_back(Vector3f(0.3, -0.06, 0.0));
	this->foot_step_list.push_back(Vector3f(0.4, 0.06, 0.0));
	this->foot_step_list.push_back(Vector3f(0.5, -0.06, 0.0));
	this->foot_step_list.push_back(Vector3f(0.5, 0.0, 0.0));
    this->foot_step_list.push_back(Vector3f(0.5, 0.0, 0.0));

}

void CapturePointControl::ref_cp_petterngenerator(){
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

void CapturePointControl::calc_cog_trajectory(){
    int cnt = 0;
    for(int t = 0; t< count; t++){

        kx1 = dt * dxi;
        lx1 = dt * ( xi - pxi ) * ACCELALETION_GRAVITY/Zc;
        ky1 = dt * dyi;
        ly1 = dt * ( yi - pyi ) * ACCELALETION_GRAVITY/Zc;

        kx2 = dt * (dxi + (lx1)/2.0f);
        lx2 = dt * (( xi - pxi ) * ACCELALETION_GRAVITY/Zc);
        ky2 = dt * (dyi + (ly1)/2.0f);
        ly2 = dt * (( yi - pyi ) * ACCELALETION_GRAVITY/Zc);

        kx3 = dt * (dxi + (lx2)/2.0f);
        lx3 = dt * (( xi - pxi ) * ACCELALETION_GRAVITY/Zc);
        ky3 = dt * (dyi + (ly2)/2.0f);
        ly3 = dt * (( yi - pyi ) * ACCELALETION_GRAVITY/Zc);

        kx4 = dt * (dxi + lx3);
        lx4 = dt * (( xi - pxi ) * ACCELALETION_GRAVITY/Zc);
        ky4 = dt * (dyi + ly3);
        ly4 = dt * (( yi - pyi ) * ACCELALETION_GRAVITY/Zc);

        xi = xi + (kx1 + 2*kx2 + 2*kx3 + kx4)/6.0f;
        yi = yi + (ky1 + 2*ky2 + 2*ky3 + ky4)/6.0f;

        dxi = dxi + (lx1 + 2*lx2 + 2*lx3 + lx4)/6.0f;
        dyi = dyi + (ly1 + 2*ly2 + 2*ly3 + ly4)/6.0f;

        /*calc capture point*/
        cpxi = xi + (dxi/omega);
        cpyi = yi + (dyi/omega);

        b = exp(omega*ref_dt_list[t]);//cnt

        pxi = (1/(1-b))*ref_cp_list_x[cnt] - ((b/(1-b))*cpxi);
        pyi = (1/(1-b))*ref_cp_list_y[cnt] - ((b/(1-b))*cpyi);

        cog_list_x.push_back(xi);
        cog_list_y.push_back(yi);
        cp_list_x.push_back(cpxi);
        cp_list_y.push_back(cpyi);
        zmp_list_x.push_back(pxi);
        zmp_list_y.push_back(pyi);
        time_list.push_back(t);

        cnt += 1;
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
