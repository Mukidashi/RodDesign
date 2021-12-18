#ifndef _BASIC_GEOMETRIC_CALCULATION_H_
#define _BASIC_GEOMETRIC_CALCULATION_H_


bool solve_three_by_two_system(float A[3][2], float b[3], float sol[2]);


bool solve_three_by_three_system(float A[3][3], float b[3], float[3]);


bool solve_three_by_three_system_use_psude_inverse(float A[3][3], float b[3], float[3]);


void line_and_triangle_intersection(float orig[3], float dire[3], float fp[3][3], float cc[3]);


float closest_point_of_triangle(float pc[3], float fps[3][3], float cord[3]);


bool intersection_of_triangle_with_line(float orig[3], float dire[3], float fps[3][3], float cord[3]);


void derive_quaternion_from_rot_mat(float Rot[3][3],float quat[4]);


void derive_rot_mat_from_quaternion(float quat[4], float Rot[3][3]);


void derive_quaternion_from_rot_mat(double Rot[3][3], double quat[4]);


void derive_rot_mat_from_quaternion(double quat[4], double Rot[3][3]);


bool intersection_of_cuboid_with_line(float orig[3], float dire[3], float rect[3][2], float& depth);


void intersection_of_triangle_with_plane(float orig[3], float norm[3], float fps[3][3], float icc[3], int &state);


float calc_dihedral_angle_of_triangles(float fps[4][3]);


float calc_ratio_of_circumcircle_radius_to_min_edge(float fps[3][3]);


float calc_distance_between_lines(float orig1[3], float dire1[3], float orig2[3], float dire2[3]);


float calc_distance_between_line_and_segment(float orig[3], float dire[3], float end1[3], float end2[3]);


float calc_distance_from_line(float orig[3], float dire[3], float poi[3]);

#endif