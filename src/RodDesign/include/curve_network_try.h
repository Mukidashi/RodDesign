#ifndef _CURVE_NETWORK_H_
#define _CURVE_NETWORK_H_

#include "loop_subdivision.h"
#include "Setting.h"

class AABB_struct;
class loop_subdivision_mesh;


struct adaptive_frame{
	float d1[3], d2[3], d3[3];
};


struct curve_node{
	int fn;
	float uv[2], pos[3], pos2[3],pos3[3];
};


struct ConstraintData{
	int pn;
	int cType;
	float pos[3];
};


struct ConnectData{
	int c1, c2;
	int n1, n2;
};


class curve_on_mesh{
	void derivative_of_bend(int in, float leng1,float leng2,float gke[2][2][3],float hkee[2][4][3][3]);
	void move_on_mesh(int &fn, float &Vk, float &Wk, float dv, float dw);
	void move_on_mesh2(int &fn, float &Vk, float &Wk, float &dv, float &dw);
	float bisection_line_search(Eigen::VectorXi fn,Eigen::VectorXd VW,Eigen::VectorXd &dVW);
	float line_search(Eigen::VectorXi fn, Eigen::VectorXd VW, Eigen::VectorXd &dVW);
	void curve_energy(Eigen::VectorXi f0,Eigen::VectorXd VW,Eigen::VectorXd dVW, float alph, float&Eval,float &dE);
	void discent_direction_of_curve_energy(Eigen::VectorXi f0, Eigen::VectorXd VW, Eigen::VectorXd &dVW);
	void calc_derivative_of_twist(int eN, adaptive_frame *mfr, float cent[3], float grm[3], float hesm[3][3]);
public:
	curve_on_mesh();
	~curve_on_mesh();

	class loop_subdivision_mesh *lmesh;
	class AABB_struct *laabb;
	int nodeN;
	struct curve_node *cNode;
	struct adaptive_frame *TPframe, *Mframe;
	//float hei, wid;
	//float *tm, *bk[2];
	struct SimulationParameter *sParam;
	float eav;
	struct ConstraintData *pCons, *dCons;
	int **Clist,*Csize;
	bool CONNECT_LIST_EXIST;

	void prepare_for_minimization();
	bool minimize_discrete_rod_energy(AABB_struct *aabb,int maxIter = 200);
	void update_frame_for_minimization(bool Initialize = false);
	void setPointConstraint(std::vector<struct ConstraintData> pcSet);
	void setDirectionConstraint(std::vector<struct ConstraintData> dcSet);
	void remesh_curve_node(int selectC, int &selectN,int &selectE,class AABB_struct *AABB,std::vector<struct ConnectData> &connect);
};


class CurveNetwork{
private:
	
public:
	CurveNetwork();
	~CurveNetwork();
	void addCurve(class curve_on_mesh *Curv);
	void deleteCurve(int cN);
	void addConnect(int c1, int n1, int c2, int n2);
	void deleteConnect(int cN);
	void setInitConnection();
	bool exportCurveNet(std::string fname, class AABB_struct *AABB);
	bool loadCurveNet(std::string fname,class loop_subdivision_mesh *imesh);
	struct qelem* traverseConnectGraph(int cN, bool *visited, struct qelem *que);


	class loop_subdivision_mesh *lmesh;
	std::vector<class curve_on_mesh *> curSet;
	std::vector<struct ConnectData> Connect;
	struct SimulationParameter *sParam;
	//float wid, hei;
};

#endif