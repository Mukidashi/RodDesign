#ifndef _DISPLAY_FUNCTIONS_H_
#define _DISPLAY_FUNCTIONS_H_


class curve_on_mesh;
class DiscreteRod;

void choose_color_data(int cflg);


void set_plastic_color(float cval[3]);


void disp_sphere(float org[3], float rad, float zoom);


void norm_for_curve_polygon(float fr1[4][3], float fr2[4][3], float nd[6][3]);


void disp_arrow(float end1[3], float end2[3], float rad, float zoom);


void disp_connection_polygon(float end1[3], float end2[3], float norm1[3], float norm2[3], float rad, float hei, float zoom);


void disp_curve_polygon(class curve_on_mesh *curv, float zoom,float off = 0.0);


void disp_rod_polygon(class DiscreteRod *rod, float zoom);


void calc_quat_for_mouse_event(float mx, float my, float quat0[4]);


void jetColorMap(float val, float rgb[3]);

#endif