#ifndef _RODNETWORK_H_
#define _RODNETWORK_H_

#include <Eigen/Core>
#include "Setting.h"

struct AdaptiveFrame{
	double d1[3], d2[3], d3[3];
};


struct ConnectedEdge{
	int rN, eN;
	double frame[3][3];
	int relate;
};

struct ConnectedNode{
	int rN, nN;
	double pos[3];
};

struct ConnectionData{
	int edgeN,nodeN;
	struct ConnectedEdge *cEdge;
	struct ConnectedNode *cNode;
	double gp[3], quat[4],norm[3];
	double Mass, Iner[3][3];
	double ogp[3], oquat[4];
};

class CurveNetwork;

class DiscreteRod{

public:
	DiscreteRod();
	~DiscreteRod();
	bool NODE_EXIST, FRAME_EXIST,PREPARED;
	int nodeN;
	//double hei, wid;

	void update_frame();
	void derivative_of_twist_bend(int in, double leng1, double leng2, double gram[2][3], double hesm[4][3][3],
		double grke[2][2][3], double grkt[2][2], double heske[2][4][3][3], double hket[2][2][2][3], double hktt[2][2]);

	double(*Nodes)[3],*oLeng,(*oNodes)[3];
	struct AdaptiveFrame *TPframe, *Mframe,*oTPframe,*oMframe;
	double *tm, (*bk)[2];
	Eigen::VectorXd Mass;
	Eigen::VectorXd twis, vXk, vTk;
	int *nvlist, *evlist;
	struct SimulationParameter *sParam;
};


class RodNetwork{
	void traverseCNlist(int N, int *ccomp,int **CNlist,int *CNsize, int cN);
	double LineSearch(Eigen::VectorXd &pk);
	void evaluateRodNetEnergy(double alph,Eigen::VectorXd pk,double &Eval,double &dE);
	void discentDirection(Eigen::VectorXd &npk);
public:
	RodNetwork();
	~RodNetwork();

	void prepare_for_simulation();
	void simulate_for_animation(float tstep,bool &tstep_is_large, bool &is_converge, bool &is_diverge);
	bool simulateFull();
	bool simulateOneStep(Eigen::VectorXd &disp,Eigen::VectorXd &resF);
	void PointsMatching();
	bool curve_to_rod(class CurveNetwork *cNet,class AABB_struct *iAABB,class loop_subdivision_mesh *lmesh);
	void exportRodNetwork(std::string fname,class Mesh *mesh);

	double old_E, initE;
	class DiscreteRod *Rods;
	struct ConnectionData *Connect;
	int rodN, conN;
	Eigen::VectorXd vTc, Wc,preSol,Lamd;
	int snodeN, slamdN, svarN;
	int *cvlist;
	
	SimulationParameter *sParam;
	double Gpos[3];

};

double twist_angle_calculation(struct AdaptiveFrame fr1, struct AdaptiveFrame fr2);



#endif