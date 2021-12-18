#ifndef _MOMENTREDUCTIONMAP_H_
#define _MOMENTREDUCTIONMAP_H_

#include <Eigen/Core>
#include "Setting.h"

class RodNetwork;
class loop_subdivision_mesh;
class AABB_struct;


struct NormalSample{
	double uv[2], norm[3], pos[3];
};

struct cutPlaneCurve{
	double thet;
	int pN;
	double(*pSet)[3];
};


class WeakRegion{
	void CutSurfaceByPlane(double pos[3], float norm[3], float dir[3], float ***pset, int &pN, float maxMove);

public:
	WeakRegion();
	~WeakRegion();
	void Initialize(SimulationParameter *Param, loop_subdivision_mesh *mesh, AABB_struct *aabb);
	void generateCurvatureMap();
	void generateMapMesh();
	double maximizeNormalCurvature();

	SimulationParameter *sParam;
	loop_subdivision_mesh *lmesh;
	AABB_struct *lAABB;
	Mesh *mMesh;
	bool MmeshExist;
	double *colmap;

	double pos[3], Mom[3];
	double Frame[3][3];
	int theN;
	struct cutPlaneCurve *pCur;
	struct NormalSample *nsamp;
	int nsN;

	double Cmat[2][2], Evec[2][2], Eval[2];
	double kgamm, Tmax;
};



class MomentReductionMap{
	void estimateWeakRegion();
	void resetMapData();

public:
	MomentReductionMap();
	~MomentReductionMap();
	void generateMomentReductionMap();

	class RodNetwork *rNet;
	class loop_subdivision_mesh *lmesh;
	class AABB_struct *lAABB;

	int wrN;
	class WeakRegion *wReg;
};



#endif