#ifndef _LOOP_SUBDIVISION_H_
#define _LOOP_SUBDIVISION_H_

#include "geometric_data_structure.h"

struct eigen_struct{
	Eigen::VectorXcf evals;
	Eigen::MatrixXcf evecs;
	Eigen::MatrixXf Amat;
	int Pmap[3][12];
	int ***Pmat;
};

class loop_subdivision_mesh :public Mesh{
private:
	bool PREPARED,Emap5_EXIST;
	void evaluate_box_spline(float cpoint[12][3], float uvw[3], float bs[3]);
	void evaluate_derivative_of_box_spline(float cpoint[12][3], float uvw[3], float dbs[2][3]);

	void set_control_points_for_box_spline_with_reflection(int fn, int neib[3][6], float cpoint[12][3]);
public:
	loop_subdivision_mesh();
	~loop_subdivision_mesh();
	bool *bcheck;
	std::map<int, struct eigen_struct*> Emap1, Emap2, Emap3;
	std::map<std::pair<int, int>, struct eigen_struct*> Emap4;
	struct eigen_struct* Emap5;
	void subdivision(void);
	void prepare_for_evaluation(void);
	bool evaluate_box_spline_on_triangle(int fn, float uvw[3],float bsp[3]);
	bool evaluate_derivative_of_box_spline_on_triangle(int fn, float uvw[3], float dbs[2][3]);
};


#endif