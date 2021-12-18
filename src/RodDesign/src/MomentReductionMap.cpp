#define _USE_MATH_DEFINES
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <string>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <Eigen/Core>
#include <vector>
#include <queue>

#include<Eigen/SparseCore>
#include<Eigen/Sparse>
#include<Eigen/Dense>
#include<Eigen/src/SparseCholesky/SimplicialCholesky.h>
#include<Eigen/src/SparseLU/SparseLU.h>
#include<Eigen/src/SparseQR/SparseQR.h>
using namespace Eigen;

using namespace std;

#include "RodNetwork.h"
#include "basic_geometric_calculation.h"
#include "geometric_data_structure.h"
#include "AABB_data_structure.h"
#include "loop_subdivision.h"
#include "MomentReductionMap.h"
#include "Setting.h"


//コンストラクタ
MomentReductionMap::MomentReductionMap(){

	wrN = 0;
}


//デストラクタ
MomentReductionMap::~MomentReductionMap(){

	printf("MR dest start\n");

	resetMapData();

	printf("MR dest end\n");
}



void MomentReductionMap::generateMomentReductionMap(){

	//データの初期化
	resetMapData();


	//Weak Regionの推定
	estimateWeakRegion();


	//曲率を求めてマップの作成
	for (int i = 0; i < wrN; i++){
		wReg[i].Initialize(rNet->sParam, lmesh, lAABB);
		wReg[i].generateCurvatureMap();
	}

}


void MomentReductionMap::estimateWeakRegion(){


	//RodNetのシミュレーションで変位の推定
	Eigen::VectorXd Disp, resF;
	rNet->prepare_for_simulation();
	rNet->simulateOneStep(Disp, resF);


	//変位大のノードを抽出
	int nnum = 0;
	float avMove = 0.0;
	for (int i = 0; i < rNet->rodN; i++){
		for (int j = 0; j < rNet->Rods[i].nodeN; j++){
			int vN = rNet->Rods[i].nvlist[j];
			float move = 0.0;
			for (int k = 0; k < 3; k++){
				move += Disp(vN + k)*Disp(vN + k);
			}
			move = sqrt(move);
			avMove += move;
			nnum++;
		}
	}
	avMove /= (float)nnum;
	float stdV = 0.0;
	float maxDev;
	for (int i = 0; i < rNet->rodN; i++){
		for (int j = 0; j < rNet->Rods[i].nodeN; j++){
			int vN = rNet->Rods[i].nvlist[j];
			float move = 0.0;
			for (int k = 0; k < 3; k++){
				move += Disp(vN + k)*Disp(vN + k);
			}
			move = sqrt(move);
			stdV += (avMove - move)*(avMove - move) / (float)nnum;
			if ((i == 0 && j == 0) || maxDev < move-avMove){
				maxDev = move - avMove;
			}
		}
	}
	stdV = sqrt(stdV);

	int **wlist;
	int **Clist;
	wlist = new int*[rNet->rodN];
	Clist = new int*[rNet->rodN];
	for (int i = 0; i < rNet->rodN; i++){
		wlist[i] = new int[rNet->Rods[i].nodeN];
		Clist[i] = new int[rNet->Rods[i].nodeN];
		for (int j = 0; j < rNet->Rods[i].nodeN; j++){
			wlist[i][j] = -1;
			Clist[i][j] = -1;
		}
	}
	for (int i = 0; i < rNet->conN; i++){
		for (int j = 0; j < rNet->Connect[i].nodeN; j++){
			int rN = rNet->Connect[i].cNode[j].rN;
			int nN = rNet->Connect[i].cNode[j].nN;
			Clist[rN][nN] = i;
		}
	}


	printf("std:%f max%f\n", stdV, maxDev);
	float wRatio = 1.0;
	for (int i = 0; i < rNet->rodN; i++){
		for (int j = 0; j < rNet->Rods[i].nodeN; j++){
			int vN = rNet->Rods[i].nvlist[j];
			float move = 0.0;
			for (int k = 0; k < 3; k++){
				move += Disp(vN + k)*Disp(vN + k);
			}
			move = sqrt(move);
			float dev = move - avMove;
			if (dev > wRatio*stdV) wlist[i][j] = 1;
		}
	}


	//連結成分からWeak Regionの推定
	std::vector<int> wCon;
	std::vector<std::pair<int, int>> wNod;
	for (int i = 0; i < rNet->rodN; i++){
		for (int j = 0; j < rNet->Rods[i].nodeN; j++){
			if (wlist[i][j] == 1){
				int cN = j;
				int wcN = -1;
				int maxN = j;
				int hasiN = -1;
				float maxMove;
				while (1){
					if (wlist[i][cN] != 1) break;
					int vN = rNet->Rods[i].nvlist[cN];
					float move = 0.0;
					for (int k = 0; k < 3; k++){
						move += Disp(vN + k)*Disp(vN + k);
					}
					move = sqrt(move);
					if (cN == j || move > maxMove) {
						maxMove = move;
						maxN = cN;
					}
					if (Clist[i][cN] != -1){
						wcN = Clist[i][cN];
						bool alreadyExist = false;
						for (int k = 0; k < wCon.size(); k++){
							if (wCon[k] == Clist[i][cN])alreadyExist = true;
						}
						if (!alreadyExist)wCon.push_back(Clist[i][cN]);
					}
					if (cN == 0 || cN == rNet->Rods[i].nodeN - 1) hasiN = cN;
					wlist[i][cN] = 0;
					cN++;
					if (cN == rNet->Rods[i].nodeN)break;
				}

				if (wcN == -1 && hasiN == -1){
					wNod.push_back(make_pair(i, maxN));
				}
				if (wcN == -1 && hasiN != -1){
					wNod.push_back(make_pair(i, hasiN));
				}

			}
		}
	}


	//モーメントの分析
	std::vector<int> dCon, dNod;
	for (int i = 0; i < wNod.size(); i++){
		int rN = wNod[i].first;
		int nN = wNod[i].second;
		float pos[3],Mom[3];
		for (int j = 0; j < 3; j++){
			pos[j] = rNet->Rods[rN].Nodes[nN][j];
			Mom[j] = 0.0;
		}
		int fn;
		float cord[3],dbs[2][3];
		bool CPsuccess = lAABB->closest_point_search_with_mesh_info(pos, &fn, cord);
		if (!CPsuccess){
			printf("CP not found\n");
		}
		lmesh->evaluate_derivative_of_box_spline_on_triangle(fn, cord, dbs);
		float nd[3];
		float ndh = 0.0;
		for (int j = 0; j < 3; j++){
			nd[j] = dbs[0][(j + 1) % 3] * dbs[1][(j + 2) % 3] - dbs[0][(j + 2) % 3] * dbs[1][(j + 1) % 3];
			ndh += nd[j] * nd[j];
		}
		float Proj[3][3];
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				if (j != k) Proj[j][k] = 0.0;
				else Proj[j][k] = 1.0;
				Proj[j][k] -= nd[j] * nd[k] / ndh;
			}
		}

		float avMom = 0.0;
		int neiN = 0;
		if (nN != 0){
			int vN = rNet->Rods[rN].nvlist[nN - 1];
			float Mk[3];
			for (int j = 0; j < 3; j++){
				Mom[j] += (rNet->Rods[rN].Nodes[nN - 1][(j + 1) % 3] - pos[(j + 1) % 3])*resF(vN + (j + 2) % 3)
					- (rNet->Rods[rN].Nodes[nN - 1][(j + 2) % 3] - pos[(j + 2) % 3])*resF(vN + (j + 1) % 3);
				//Mom[j] += rNet->Rods[rN].Mframe[nN - 1].d3[j] * resF(vN + 3);
				Mk[j] = (rNet->Rods[rN].Nodes[nN - 1][(j + 1) % 3] - pos[(j + 1) % 3])*resF(vN + (j + 2) % 3)
					- (rNet->Rods[rN].Nodes[nN - 1][(j + 2) % 3] - pos[(j + 2) % 3])*resF(vN + (j + 1) % 3);
			}
			float momV = 0.0;
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 3; k++){
					momV += Mk[j] * Mk[k] * Proj[j][k];
				}
			}
			momV = sqrt(momV);
			avMom += momV;
			neiN++;
		}
		if (nN != rNet->Rods[rN].nodeN - 1){
			int vN = rNet->Rods[rN].nvlist[nN + 1];
			float Mk[3];
			for (int j = 0; j < 3; j++){
				Mom[j] += (rNet->Rods[rN].Nodes[nN + 1][(j + 1) % 3] - pos[(j + 1) % 3])*resF(vN + (j + 2) % 3)
					- (rNet->Rods[rN].Nodes[nN + 1][(j + 2) % 3] - pos[(j + 2) % 3])*resF(vN + (j + 1) % 3);
				//Mom[j] += rNet->Rods[rN].Mframe[nN].d3[j] * resF(vN - 1);
				Mk[j] = (rNet->Rods[rN].Nodes[nN + 1][(j + 1) % 3] - pos[(j + 1) % 3])*resF(vN + (j + 2) % 3)
					- (rNet->Rods[rN].Nodes[nN + 1][(j + 2) % 3] - pos[(j + 2) % 3])*resF(vN + (j + 1) % 3);
			}
			float momV = 0.0;
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 3; k++){
					momV += Mk[j] * Proj[j][k] * Mk[k];
				}
			}
			momV = sqrt(momV);
			avMom += momV;
			neiN++;
		}
		float momV = 0.0;
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				momV += Mom[j] * Proj[j][k] * Mom[k];
			}
		}
		momV = sqrt(momV);
		avMom /= (float)neiN;
		//printf("%d %d %e->%e %d\n", rN, nN, avMom, momV,neiN);
		if (neiN == 1 || avMom > momV) dNod.push_back(i);
		else dNod.push_back(i);
	}

	for (int i = 0; i < wCon.size(); i++){
		float pos[3], Mom[3];
		for (int j = 0; j < 3; j++){
			pos[j] = rNet->Connect[wCon[i]].gp[j];
			Mom[j] = 0.0;
		}
		int fn;
		float cord[3],dbs[2][3];
		bool CPsuccess = lAABB->closest_point_search_with_mesh_info(pos, &fn, cord);
		if (!CPsuccess){
			printf("CP not found\n");
		}
		lmesh->evaluate_derivative_of_box_spline_on_triangle(fn, cord, dbs);
		float nd[3];
		float ndh = 0.0;
		for (int j = 0; j < 3; j++){
			nd[j] = dbs[0][(j + 1) % 3] * dbs[1][(j + 2) % 3] - dbs[0][(j + 2) % 3] * dbs[1][(j + 1) % 3];
			ndh += nd[j] * nd[j];
		}
		float Proj[3][3];
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				if (j != k)Proj[j][k] = 0.0;
				else Proj[j][k] = 1.0;
				Proj[j][k] -= nd[j] * nd[k] / ndh;
			}
		}

		float avMom = 0.0;
		for (int j = 0; j < rNet->Connect[wCon[i]].nodeN; j++){
			int rN = rNet->Connect[wCon[i]].cNode[j].rN;
			int nN = rNet->Connect[wCon[i]].cNode[j].nN;
			int vN = rNet->Rods[rN].nvlist[nN];
			float Mk[3];
			for (int k = 0; k < 3; k++){
				Mom[k] += (rNet->Rods[rN].Nodes[nN][(k + 1) % 3] - pos[(k + 1) % 3])*resF(vN + (k + 2) % 3)
					- (rNet->Rods[rN].Nodes[nN][(k + 2) % 3] - pos[(k + 2) % 3])*resF(vN + (k + 1) % 3);
				Mk[k] += (rNet->Rods[rN].Nodes[nN][(k + 1) % 3] - pos[(k + 1) % 3])*resF(vN + (k + 2) % 3)
					- (rNet->Rods[rN].Nodes[nN][(k + 2) % 3] - pos[(k + 2) % 3])*resF(vN + (k + 1) % 3);
			}
			float momV = 0.0;
			for (int k = 0; k < 3; k++){
				for (int l = 0; l < 3; l++){
					momV += Mk[k] * Mk[l] * Proj[k][l];
				}
			}
			momV = sqrt(momV);
			avMom += momV;
		}
		/*for (int j = 0; j < rNet->Connect[wCon[i]].edgeN; j++){
			int rN = rNet->Connect[wCon[i]].cEdge[j].rN;
			int eN = rNet->Connect[wCon[i]].cEdge[j].eN;
			int vN = rNet->Rods[rN].evlist[eN];
			for (int k = 0; k < 3; k++){
				Mom[k] += rNet->Rods[rN].Mframe[eN].d3[k] * resF(vN);
			}
		}*/
		float momV = 0.0;
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				momV += Mom[j] * Proj[j][k] * Mom[k];
			}
		}
		momV = sqrt(momV);
		avMom /= (float)rNet->Connect[wCon[i]].nodeN;
		//printf("%dcon %e->%e\n", wCon[i], avMom, momV);
		if (avMom >momV)dCon.push_back(i);
	}


	//モーメントの計算
	wrN = wCon.size() + wNod.size() - dCon.size() - dNod.size();
	if(wrN > 0) wReg = new class WeakRegion[wrN];

	wrN = 0;
	for (int i = 0; i < wNod.size(); i++){
		bool needNode = true;
		for (int j = 0; j < dNod.size(); j++){
			if (dNod[j] == i) needNode = false;
		}
		if (!needNode)continue;

		int rN = wNod[i].first;
		int nN = wNod[i].second;
		for (int j = 0; j < 3; j++){
			wReg[wrN].pos[j] = rNet->Rods[rN].Nodes[nN][j];
			wReg[wrN].Mom[j] = 0.0;
		}
		if (nN != 0){
			int vN = rNet->Rods[rN].nvlist[nN - 1];
			for (int j = 0; j < 3; j++){
				wReg[wrN].Mom[j] += (rNet->Rods[rN].Nodes[nN - 1][(j + 1) % 3] - wReg[wrN].pos[(j + 1) % 3])*resF(vN + (j + 2) % 3)
					- (rNet->Rods[rN].Nodes[nN - 1][(j + 2) % 3] - wReg[wrN].pos[(j + 2) % 3])*resF(vN + (j + 1) % 3);
				wReg[wrN].Mom[j] += rNet->Rods[rN].Mframe[nN - 1].d3[j] * resF(vN + 3);
			}
		}
		if (nN != rNet->Rods[rN].nodeN - 1){
			int vN = rNet->Rods[rN].nvlist[nN + 1];
			for (int j = 0; j < 3; j++){
				wReg[wrN].Mom[j] += (rNet->Rods[rN].Nodes[nN + 1][(j + 1) % 3] - wReg[wrN].pos[(j + 1) % 3])*resF(vN + (j + 2) % 3)
					- (rNet->Rods[rN].Nodes[nN + 1][(j + 2) % 3] - wReg[wrN].pos[(j + 2) % 3])*resF(vN + (j + 1) % 3);
				wReg[wrN].Mom[j] += rNet->Rods[rN].Mframe[nN].d3[j] * resF(vN - 1);
			}
		}
		wrN++;
	}
	for (int i = 0; i < wCon.size(); i++){
		bool needCon = true;
		for (int j = 0; j < dCon.size(); j++){
			if (dCon[j] == i) needCon = false;
		}
		if (!needCon) continue;

		for (int j = 0; j < 3; j++){
			wReg[wrN].pos[j] = rNet->Connect[wCon[i]].gp[j];
			wReg[wrN].Mom[j] = 0.0;
		}
		for (int j = 0; j < rNet->Connect[wCon[i]].nodeN; j++){
			int rN = rNet->Connect[wCon[i]].cNode[j].rN;
			int nN = rNet->Connect[wCon[i]].cNode[j].nN;
			int vN = rNet->Rods[rN].nvlist[nN];
			float momV = 0.0;
			for (int k = 0; k < 3; k++){
				wReg[wrN].Mom[k] += (rNet->Rods[rN].Nodes[nN][(k + 1) % 3] - wReg[wrN].pos[(k + 1) % 3])*resF(vN + (k + 2) % 3)
					- (rNet->Rods[rN].Nodes[nN][(k + 2) % 3] - wReg[wrN].pos[(k + 2) % 3])*resF(vN + (k + 1) % 3);
			}
		}
		for (int j = 0; j < rNet->Connect[wCon[i]].edgeN; j++){
			int rN = rNet->Connect[wCon[i]].cEdge[j].rN;
			int eN = rNet->Connect[wCon[i]].cEdge[j].eN;
			int vN = rNet->Rods[rN].evlist[eN];
			for (int k = 0; k < 3; k++){
				wReg[wrN].Mom[k] += rNet->Rods[rN].Mframe[eN].d3[k] * resF(vN);
			}
		}

		wrN++;
	}



	for (int i = 0; i < rNet->rodN; i++){
		delete[] wlist[i];
		delete[] Clist[i];
	}
	delete[] wlist;
	delete[] Clist;
}


void MomentReductionMap::resetMapData(){

	if (wrN != 0){
		delete[] wReg;
	}

}



///////////
////WeakRegion
///////////

WeakRegion::WeakRegion(){

	theN = 0;

	MmeshExist = false;

}


WeakRegion::~WeakRegion(){

	printf("Weak dest start\n");
	if (theN != 0){
		for (int i = 0; i < theN; i++){
			delete[] pCur[i].pSet;
		}
		delete[] pCur;
	}

	if (MmeshExist){
		delete mMesh;
		delete[] colmap;
	}
	printf("Weak dest end\n");
}


void WeakRegion::Initialize(SimulationParameter *Param, loop_subdivision_mesh *mesh, AABB_struct *aabb){

	sParam = Param;
	lmesh = mesh;
	lAABB = aabb;

}


void WeakRegion::generateCurvatureMap(){

	//フレームの計算
	float np[3], scp[3], duv[2][3];
	for (int j = 0; j < 3; j++){
		np[j] = pos[j];
	}
	int fn;
	float cord[3];
	bool success = lAABB->closest_point_search_with_mesh_info(np, &fn, cord);
	lmesh->evaluate_box_spline_on_triangle(fn, cord, scp);
	lmesh->evaluate_derivative_of_box_spline_on_triangle(fn, cord, duv);

	double nd[3];
	float ndh = 0.0;
	for (int j = 0; j < 3; j++){
		nd[j] = duv[0][(j + 1) % 3] * duv[1][(j + 2) % 3] - duv[0][(j + 2) % 3] * duv[1][(j + 1) % 3];
		ndh += nd[j] * nd[j];
	}
	ndh = sqrt(ndh);
	if (ndh < 1.0e-10) ndh = 1.0e-10;
	double Mh = 0.0;
	double seki = 0.0;
	for (int j = 0; j < 3; j++){
		nd[j] /= ndh;
		Frame[2][j] = nd[j];
		seki += -Frame[2][j] * Mom[j];
	}
	for (int j = 0; j < 3; j++){
		Frame[1][j] = -Mom[j] - seki*nd[j];
		Mh += Frame[1][j] * Frame[1][j];
	}
	Mh = sqrt(Mh);
	if (Mh < 1.0e-10) Mh = 1.0e-10;
	for (int j = 0; j < 3; j++){
		Frame[1][j] /= Mh;
	}
	Mh = 0.0;
	for (int j = 0; j < 3; j++){
		Frame[0][j] = Frame[1][(j + 1) % 3] * Frame[2][(j + 2) % 3] - Frame[1][(j + 2) % 3] * Frame[2][(j + 1) % 3];
		Mh += Frame[0][j] * Frame[0][j];
	}
	Mh = sqrt(Mh);
	for (int j = 0; j < 3; j++){
		Frame[0][j] /= Mh;
	}

	//面切断による曲線集合の作成
	theN = 36;
	pCur = new cutPlaneCurve[theN];

	for (int i = 0; i < theN; i++){
		pCur[i].thet = 2.0*M_PI / (double)theN*(double)i;
		double theta = pCur[i].thet + M_PI / 2.0;
		float norm[3], dir[3], opos[3];
		for (int j = 0; j < 3; j++){
			norm[j] = cos(theta)*Frame[0][j] + sin(theta)*Frame[1][j];
			opos[j] = pos[j];
		}
		float nh = 0.0;
		float dh = 0.0;
		for (int j = 0; j < 3; j++){
			dir[j] = norm[(j + 1) % 3] * Frame[2][(j + 2) % 3] - norm[(j + 2) % 3] * Frame[2][(j + 1) % 3];
			nh += norm[j] * norm[j];
			dh += dir[j] * dir[j];
		}
		nh = sqrt(nh);
		dh = sqrt(dh);
		for (int j = 0; j < 3; j++){
			norm[j] /= nh;
			dir[j] /= dh;
		}
		float **pset;
		int pN;
		float leng = 0.0;
		for (int j = 0; j < 4; j++){
			leng += sParam->sampL[j];
		}
		leng += 1.0;
		CutSurfaceByPlane(pos, norm, dir, &pset, pN, leng);

		pCur[i].pN = pN;
		if (pN != 0) pCur[i].pSet = new double[pN][3];
		for (int j = 0; j < pN; j++){
			for (int k = 0; k < 3; k++){
				pCur[i].pSet[j][k] = pset[j][k];
			}
		}

		if (pN > 0){
			for (int j = 0; j < pN; j++){
				delete[] pset[j];
			}
			delete[] pset;
		}
	}


	//Ｍ推定による曲率テンソルの計算
	//サンプル生成
	int tN = theN / 2;
	int lN = 3;
	nsamp = new struct NormalSample[tN*lN];

	float haba = (sParam->sampL[0] + sParam->sampL[1] +sParam->sampL[2]) / (float)lN;

	nsN = 0;
	for (int i = 0; i < tN; i++){
		int pcN = (2 * i) % theN;

		int stN = 0;
		float hleng = 0.0;
		float accL = 0.0;
		for (int j = 0; j < lN; j++){
			float tarL = ((float)j + 1.0)*haba - accL;
			for (int k = stN; k < pCur[pcN].pN - 1; k++){
				float leng = 0.0;
				for (int l = 0; l < 3; l++){
					leng += pow(pCur[pcN].pSet[k][l] - pCur[pcN].pSet[k + 1][l], 2);
				}
				leng = sqrt(leng);
				if (hleng + leng > tarL){
					float coef = (hleng + leng - tarL) / leng;
					for (int l = 0; l < 3; l++){
						nsamp[nsN].pos[l] = coef*pCur[pcN].pSet[k][l] + (1.0 - coef)*pCur[pcN].pSet[k + 1][l];
					}
					nsamp[nsN].uv[0] = cos(pCur[pcN].thet)*tarL;
					nsamp[nsN].uv[1] = sin(pCur[pcN].thet)*tarL;
					nsN++;
					stN = k;
					break;
				}
				else {
					hleng += leng;
					if (k == pCur[pcN].pN - 2) {
						for (int l = 0; l < 3; l++){
							nsamp[nsN].pos[l] = pCur[pcN].pSet[k + 1][l];
						}
						nsamp[nsN].uv[0] = cos(pCur[pcN].thet)*hleng;
						nsamp[nsN].uv[1] = sin(pCur[pcN].thet)*hleng;
						nsN++;
						stN = 0;
						accL = hleng;
						hleng = 0.0;
					}
				}
			}
		}
	}
	for (int i = 0; i < nsN; i++){
		int fn;
		float cord[3], pc[3], dbs[2][3];
		for (int j = 0; j < 3; j++){
			pc[j] = nsamp[i].pos[j];
		}
		lAABB->closest_point_search_with_mesh_info(pc, &fn, cord);
		lmesh->evaluate_derivative_of_box_spline_on_triangle(fn, cord, dbs);

		double ndh = 0.0;
		for (int j = 0; j < 3; j++){
			nsamp[i].norm[j] = dbs[0][(j + 1) % 3] * dbs[1][(j + 2) % 3] - dbs[0][(j + 2) % 3] * dbs[1][(j + 1) % 3];
			ndh += pow(nsamp[i].norm[j], 2);
		}
		ndh = sqrt(ndh);
		for (int j = 0; j < 3; j++){
			nsamp[i].norm[j] /= ndh;
		}
	}

	//曲率計算
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j++){
			Cmat[i][j] = 0.0;
		}
	}
	int iterN = 0;
	int maxIter = 10;
	double sigm;
	double cmag;
	while (1){
		Eigen::Matrix2d Mmat;
		Eigen::Vector2d bv1, bv2;
		Mmat = Eigen::Matrix2d::Zero();
		bv1 = Eigen::Vector2d::Zero();
		bv2 = Eigen::Vector2d::Zero();

		for (int i = 0; i < nsN; i++){
			for (int j = 0; j < i; j++){
				double pu = 0.0;
				double pv = 0.0;
				double nu = 0.0;
				double nv = 0.0;
				for (int k = 0; k < 3; k++){
					pu += Frame[0][k] * (nsamp[i].pos[k] - nsamp[j].pos[k]);
					pv += Frame[1][k] * (nsamp[i].pos[k] - nsamp[j].pos[k]);
					nu += Frame[0][k] * (nsamp[i].norm[k] - nsamp[j].norm[k]);
					nv += Frame[1][k] * (nsamp[i].norm[k] - nsamp[j].norm[k]);
				}
				double wei = 1.0;
				if (iterN != 0){
					double dif = 0.0;
					dif += pow(Cmat[0][0] * pu + Cmat[0][1] * pv - nu, 2);
					dif += pow(Cmat[1][0] * pu + Cmat[1][1] * pv - nv, 2);
					dif = sqrt(dif);
					wei = 2.0 / pow(1 + dif*dif / sigm / sigm, 2);
					if (wei < 0.2) wei = 0.2;
				}
				Mmat(0, 0) += wei*pu*pu;
				Mmat(0, 1) += wei*pu*pv;
				Mmat(1, 0) += wei*pu*pv;
				Mmat(1, 1) += wei*pv*pv;
				bv1(0) += wei*nu*pu;
				bv1(1) += wei*nu*pv;
				bv2(0) += wei*nv*pu;
				bv2(1) += wei*nv*pv;
			}
		}

		//行列の計算
		Eigen::Vector2d sol1 = Mmat.inverse()*bv1;
		Eigen::Vector2d sol2 = Mmat.inverse()*bv2;

		float dif = 0.0;
		for (int i = 0; i < 2; i++){
			dif += pow(Cmat[0][i] - sol1(i), 2);
			dif += pow(Cmat[1][i] - sol2(i), 2);
			Cmat[0][i] = sol1(i);
			Cmat[1][i] = sol2(i);
		}
		dif = sqrt(dif);
		if (iterN == 0)cmag = dif;


		//パラメータの計算
		std::priority_queue<float> dque;
		for (int i = 0; i < nsN; i++){
			for (int j = 0; j < i; j++){
				double pu = 0.0;
				double pv = 0.0;
				double nu = 0.0;
				double nv = 0.0;
				for (int k = 0; k < 3; k++){
					pu += Frame[0][k] * (nsamp[i].pos[k] - nsamp[j].pos[k]);
					pv += Frame[1][k] * (nsamp[i].pos[k] - nsamp[j].pos[k]);
					nu += Frame[0][k] * (nsamp[i].norm[k] - nsamp[j].norm[k]);
					nv += Frame[1][k] * (nsamp[i].norm[k] - nsamp[j].norm[k]);
				}
				double dif = 0.0;
				dif += pow(Cmat[0][0] * pu + Cmat[0][1] * pv - nu, 2);
				dif += pow(Cmat[1][0] * pu + Cmat[1][1] * pv - nv, 2);
				dif = sqrt(dif);
				dque.push(dif);
			}
		}
		float medV;
		for (int i = 0; i < nsN*(nsN - 1) / 4; i++){
			medV = dque.top();
			dque.pop();
		}
		sigm = medV*1.4826;

		if (iterN != 0 && dif < cmag*0.01) break;
		iterN++;
		if (iterN == maxIter) break;
	}


	//固有値計算により主曲率方向をもとめる
	Eigen::Matrix2d Ctens;
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j++){
			Ctens(i, j) = Cmat[i][j];
		}
	}

	SelfAdjointEigenSolver<Matrix2d> es(Ctens);
	for (int i = 0; i < 2; i++){
		Eval[i] = es.eigenvalues()[i];
		for (int j = 0; j < 2; j++){
			Evec[i][j] = es.eigenvectors()(i, j);
		}
	}

	double eleng = sqrt(Evec[0][0] * Evec[0][0] + Evec[1][0] * Evec[1][0]);
	kgamm = Evec[0][0] / eleng;
	if (kgamm > 1.0) kgamm = 1.0;
	if (kgamm < -1.0) kgamm = -1.0;
	kgamm = acos(kgamm);
	if (Evec[1][0] < 0) kgamm *= -1.0;


	Tmax = maximizeNormalCurvature();
	printf("Eval:%f %f\n", Eval[0], Eval[1]);
	printf("max %f\n", Tmax*180.0 / M_PI);


	//マップのメッシュ作成
	generateMapMesh();

}


double WeakRegion::maximizeNormalCurvature(){

	double the0 = 0.0;
	double ks = cos(the0)*(Eval[0] * cos(-kgamm)*cos(-kgamm) + Eval[1] * sin(-kgamm)*sin(-kgamm));
	if (ks < 0.0){
		the0 += M_PI;
		ks = cos(the0)*(Eval[0] * cos(the0 - kgamm)*cos(the0 - kgamm) + Eval[1] * sin(the0 - kgamm)*sin(the0 - kgamm));
	}

	double thk = the0 - kgamm;
	double coef[4];
	coef[0] = Eval[0] * cos(kgamm);
	coef[1] = Eval[0] * sin(kgamm);
	coef[2] = Eval[1] * cos(kgamm);
	coef[3] = Eval[1] * sin(kgamm);

	Eigen::Matrix3d Amat;
	Eigen::Vector3d bvec, sol;
	int iterN = 0;
	int maxIter = 100;
	while (1){
		double xv = cos(thk);
		double yv = sin(thk);

		double Eval = 0.0;
		Eval = -coef[0] * pow(xv, 3) + coef[1] * xv*xv*yv - coef[2] * xv*yv*yv + coef[3] * pow(yv, 3);

		bvec(0) = 3.0*coef[0] * xv*xv - 2.0*coef[1] * xv*yv + coef[2] * yv*yv;
		bvec(1) = -coef[1] * xv*xv + 2.0*coef[2] * xv*yv - 3.0*coef[3] * yv*yv;
		bvec(2) = 1.0 - xv*xv - yv*yv;


		Amat(0, 0) = -6.0*coef[0] * xv + 2.0*coef[1] * yv;
		Amat(0, 1) = 2.0*coef[1] * xv - 2.0*coef[2] * yv;
		Amat(0, 2) = 2.0*xv;
		Amat(1, 0) = 2.0*coef[1] * xv - 2.0*coef[2] * yv;
		Amat(1, 1) = -2.0*coef[2] * xv + 6.0*coef[3] * yv;
		Amat(1, 2) = 2.0*yv;
		Amat(2, 0) = 2.0*xv;
		Amat(2, 1) = 2.0*yv;
		Amat(2, 2) = 0.0;

		sol = Amat.inverse()*bvec;
		double seki = sol(0)*bvec(0) + sol(1)*bvec(1);
		if (seki < 0.0){
			sol(0) = -sol(0);
			sol(1) = -sol(1);
		}
		double alph = 1.0;
		while (1){
			double nx = xv + alph*sol(0);
			double ny = yv + alph*sol(1);
			double ntk = nx / sqrt(nx*nx + ny*ny);
			if (ntk > 1.0) ntk = 1.0;
			if (ntk < -1.0) ntk = -1.0;
			ntk = acos(ntk);
			if (ny < 0.0) ntk *= -1.0;

			nx = cos(ntk);
			ny = sin(ntk);
			double nEval = -coef[0] * pow(nx, 3) + coef[1] * nx*nx*ny - coef[2] * nx*ny*ny + coef[3] * pow(ny, 3);
			if (nEval < Eval) break;
			alph *= 0.5;
			if (alph < 1.0e-6) {
				alph = 0.0;
				break;
			}
		}

		xv += alph*sol(0);
		yv += alph*sol(1);
		double ntk = xv / sqrt(xv*xv + yv*yv);
		if (ntk > 1.0)ntk = 1.0;
		if (ntk < -1.0)ntk = -1.0;
		ntk = acos(ntk);
		if (yv < 0.0) ntk *= -1.0;

		double dif = fabs(ntk - thk);
		double nx = cos(ntk);
		double ny = sin(ntk);
		double nEval = -coef[0] * pow(nx, 3) + coef[1] * nx*nx*ny - coef[2] * nx*ny*ny + coef[3] * pow(ny, 3);
		//printf("%d %f(%f) %f\n", iterN, nEval, Eval, (ntk + kgamm)*180.0 / M_PI);

		if (dif < M_PI*0.01 / 180.0) break;
		thk = ntk;

		iterN++;
		if (iterN == maxIter) {
			break;
		}

	}

	return thk + kgamm;

}


void WeakRegion::generateMapMesh(){

	if (MmeshExist){
		delete mMesh;
	}
	MmeshExist = true;
	mMesh = new Mesh;

	double minL = sParam->sampL[0] + sParam->sampL[1] / 100.0;
	double maxL = sParam->sampL[0] + sParam->sampL[1] / 2.0;
	double offV = sParam->Hei*0.1;
	int divN = 10;

	float(*Nodes)[3];
	Nodes = new float[(divN + 1)*theN][3];
	int nodeN = 0;

	int(*nlist)[2];
	nlist = new int[theN][2];
	for (int i = 0; i < theN; i++){

		int stN = 0;
		double aLeng = 0.0;
		nlist[i][0] = nodeN;
		for (int j = 0; j < divN + 1; j++){
			double tarL = minL + (double)j*(maxL - minL) / (double)divN;
			for (int k = stN; k < pCur[i].pN - 1; k++){
				double leng = 0.0;
				for (int l = 0; l < 3; l++){
					leng += pow(pCur[i].pSet[k + 1][l] - pCur[i].pSet[k][l], 2);
				}
				leng = sqrt(leng);
				if (aLeng + leng >= tarL){
					double coef = (aLeng + leng - tarL) / leng;
					for (int l = 0; l < 3; l++){
						Nodes[nodeN][l] = coef*pCur[i].pSet[k][l] + (1.0 - coef)*pCur[i].pSet[k + 1][l];
					}
					nodeN++;
					stN = k;
					break;
				}
				else {
					aLeng += leng;
					if (k == pCur[i].pN - 2) stN = k;
				}
			}
		}
		nlist[i][1] = nodeN;

	}

	int(*Faces)[3];
	int faceN = 0;
	Faces = new int[2 * divN*theN][3];
	for (int i = 0; i < theN; i++){
		int tn1 = i;
		int tn2 = (i + 1) % theN;
		if (nlist[tn1][1] == nlist[tn1][0] || nlist[tn2][1] == nlist[tn2][0]) continue;

		for (int j = 0; j < divN; j++){
			int fn1 = nlist[tn1][0] + j;
			int fn2 = nlist[tn2][0] + j;
			if (fn1 + 1 >= nlist[tn1][1] || fn2 + 1 >= nlist[tn2][1]) break;

			Faces[faceN][0] = fn1;
			Faces[faceN][1] = fn1 + 1;
			Faces[faceN][2] = fn2;
			faceN++;
			Faces[faceN][0] = fn1 + 1;
			Faces[faceN][1] = fn2 + 1;
			Faces[faceN][2] = fn2;
			faceN++;
		}
	}

	mMesh->nodeN = nodeN;
	mMesh->Nodes = new float[nodeN][3];
	for (int i = 0; i < nodeN; i++){
		float pos[3];
		for (int j = 0; j < 3; j++){
			mMesh->Nodes[i][j] = Nodes[i][j];
			pos[j] = mMesh->Nodes[i][j];
		}
		int fn;
		float cord[3];
		lAABB->closest_point_search_with_mesh_info(pos, &fn, cord);

		float norm[3];
		float ndh = 0.0;
		for (int j = 0; j < 3; j++){
			norm[j] = (lmesh->Nodes[lmesh->Faces[fn][1]][(j + 1) % 3] - lmesh->Nodes[lmesh->Faces[fn][0]][(j + 1) % 3])
				*(lmesh->Nodes[lmesh->Faces[fn][2]][(j + 2) % 3] - lmesh->Nodes[lmesh->Faces[fn][0]][(j + 2) % 3])
				- (lmesh->Nodes[lmesh->Faces[fn][1]][(j + 2) % 3] - lmesh->Nodes[lmesh->Faces[fn][0]][(j + 2) % 3])
				*(lmesh->Nodes[lmesh->Faces[fn][2]][(j + 1) % 3] - lmesh->Nodes[lmesh->Faces[fn][0]][(j + 1) % 3]);
			ndh += norm[j] * norm[j];
		}
		ndh = sqrt(ndh);
		for (int j = 0; j < 3; j++){
			mMesh->Nodes[i][j] += offV*norm[j] / ndh;
		}
	}
	mMesh->faceN = faceN;
	mMesh->Faces = new int[faceN][3];
	for (int i = 0; i < faceN; i++){
		for (int j = 0; j < 3; j++){
			mMesh->Faces[i][j] = Faces[i][j];
		}
	}


	double ksMax = cos(Tmax)*(Eval[0] * cos(Tmax - kgamm)*cos(Tmax - kgamm) + Eval[1] * sin(Tmax - kgamm)*sin(Tmax - kgamm));
	//ksMax = fabs(Eval[1]);
	colmap = new double[nodeN];
	for (int i = 0; i < theN; i++){
		if (nlist[i][0] == nlist[i][1]) continue;

		double thet = pCur[i].thet;
		double ks = cos(thet)*(Eval[0] * cos(thet - kgamm)*cos(thet - kgamm) + Eval[1] * sin(thet - kgamm)*sin(thet - kgamm));
		if (ksMax < ks) ksMax = ks;
		/*for (int j = nlist[i][0]; j < nlist[i][1]; j++){
			colmap[j] = ks / ksMax;
		}*/
	}
	for (int i = 0; i < theN; i++){
		if (nlist[i][0] == nlist[i][1]) continue;

		double thet = pCur[i].thet;
		double ks = cos(thet)*(Eval[0] * cos(thet - kgamm)*cos(thet - kgamm) + Eval[1] * sin(thet - kgamm)*sin(thet - kgamm));
		for (int j = nlist[i][0]; j < nlist[i][1]; j++){
			colmap[j] = ks / ksMax;
		}
	}

	delete[] Nodes;
	delete[]nlist;
	delete[] Faces;
	
}


void WeakRegion::CutSurfaceByPlane(double pos[3], float norm[3], float dir[3], float ***pset, int &pN, float maxMove){


	int cstate;
	float icc[3];
	float np[3], scp[3], duv[2][3];
	for (int j = 0; j < 3; j++){
		np[j] = pos[j];
	}
	int fn;
	float cord[3];
	bool success = lAABB->closest_point_search_with_mesh_info(np, &fn, cord);

	if (!success){
		printf("CP not found@Cut Surface\n");

	}
	float fps[3][3];
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			fps[i][j] = lmesh->Nodes[lmesh->Faces[fn][i]][j];
		}
	}
	intersection_of_triangle_with_plane(np, norm, fps, icc, cstate);

	float pd[3];
	float pdh = 0.0;
	for (int i = 0; i < 3; i++){
		pd[i] = norm[(i + 1) % 3] * dir[(i + 2) % 3] - norm[(i + 2) % 3] * dir[(i + 1) % 3];
		pdh += pd[i] * pd[i];
	}
	pdh = sqrt(pdh);
	if (pdh < 1.0e-10)pdh = 1.0e-10;
	for (int i = 0; i < 3; i++){
		pd[i] /= pdh;
	}


	//ここから
	pN = 0;

	//スタート点での処理
	std::vector<int> sens;
	std::vector<float> scoef;
	if (cstate == -1) return;
	if (cstate == 1){
		int ipn = -1;
		for (int i = 0; i < 3; i++){
			if (icc[i] == 0.0) ipn = i;
		}
		if (ipn == -1) return;
		int en = lmesh->FElist[fn][(ipn + 2) % 3];
		float coef = 1.0;
		if (lmesh->Edges[en][0] != lmesh->Faces[fn][ipn]) coef = 0.0;
		sens.push_back(en);
		scoef.push_back(coef);
		en = lmesh->FElist[fn][ipn];
		coef = 1.0;
		if (lmesh->Edges[en][0] != lmesh->Faces[fn][ipn]) coef = 0.0;
		sens.push_back(en);
		scoef.push_back(coef);
	}
	else if (cstate == 2){
		int opn = -1;
		for (int i = 0; i < 3; i++){
			if (icc[i]>0.9)opn = i;
		}
		if (opn == -1) return;
		int en = lmesh->FElist[fn][opn];
		float coef = 0.0;
		if (lmesh->Edges[en][0] != lmesh->Faces[fn][opn]) coef = 1.0;
		sens.push_back(en);
		scoef.push_back(coef);
		en = lmesh->FElist[fn][(opn + 2) % 3];
		coef = 0.0;
		if (lmesh->Edges[en][0] != lmesh->Faces[fn][opn]) coef = 1.0;
		sens.push_back(en);
		scoef.push_back(coef);
	}
	else if (cstate == 3){
		return;
	}
	else {
		int intN = 0;
		int oen = -1;
		for (int i = 0; i < 3; i++){
			if (icc[i] >= -0.0 && icc[i] <= 1.0) intN++;
			else oen = i;
		}
		if (oen == -1 && intN != 2) return;
		int en = lmesh->FElist[fn][(oen + 1) % 3];
		float coef = icc[(oen + 1) % 3];
		if (lmesh->Edges[en][0] == lmesh->Faces[fn][(oen + 1) % 3]) coef = 1.0 - coef;
		sens.push_back(en);
		scoef.push_back(coef);
		en = lmesh->FElist[fn][(oen + 2) % 3];
		coef = icc[(oen + 2) % 3];
		if (lmesh->Edges[en][0] == lmesh->Faces[fn][(oen + 2) % 3]) coef = 1.0 - coef;
		sens.push_back(en);
		scoef.push_back(coef);
	}


	//境界追跡
	bool start_is_vertex = false;
	int startp = -1;
	if (cord[0] > 0.9){
		startp = lmesh->Faces[fn][0];
		start_is_vertex = true;;
	}
	else if (cord[1] > 0.9){
		startp = lmesh->Faces[fn][1];
		start_is_vertex = true;
	}
	else if (cord[2] > 0.9){
		startp = lmesh->Faces[fn][2];
		start_is_vertex = true;
	}

	std::vector<std::vector<int>> Cfset, Ceset;
	std::vector<std::vector<float>> Ccoefs;
	for (int i = 0; i < sens.size(); i++){
		std::vector<int> eset, fset;
		std::vector<float> coefs;
		int cen = sens[i];
		int cfn = fn;
		float seki = 0.0;
		for (int j = 0; j < 3; j++){
			seki += dir[j] * (scoef[i] * lmesh->Nodes[lmesh->Edges[cen][0]][j] + (1.0 - scoef[i])*lmesh->Nodes[lmesh->Edges[cen][1]][j] - pos[j]);
		}
		if (seki < 0.0) continue;
		fset.push_back(cfn);
		eset.push_back(sens[i]);
		coefs.push_back(scoef[i]);

		float Mleng = 0.0;
		bool MoveSuccess = false;
		while (1){
			int nfn = lmesh->EFlist[cen][0];
			if (nfn == cfn) nfn = lmesh->EFlist[cen][1];
			if (nfn == -1) break;

			int nen;
			float ncoef;
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 3; k++){
					fps[j][k] = lmesh->Nodes[lmesh->Faces[nfn][j]][k];
				}
			}
			intersection_of_triangle_with_plane(np, norm, fps, icc, cstate);
			if (cstate != 0) printf("state:%d\n", cstate);
			if (cstate == -1) break;
			if (cstate == 1){
				int ipn = -1;
				for (int j = 0; j < 3; j++){
					if (fabs(icc[j]) < 1.0e-10) ipn = j;
				}
				if (ipn == -1)break;
				int ens[2];
				ens[0] = lmesh->FElist[nfn][(ipn + 2) % 3];
				ens[1] = lmesh->FElist[nfn][ipn];
				nen = ens[0];
				if (ens[0] == cen)nen = ens[1];
				ncoef = 0.0;
				if (lmesh->Edges[nen][0] == lmesh->Faces[nfn][ipn])ncoef = 1.0;
			}
			else if (cstate == 2){
				int opn = -1;
				for (int j = 0; j < 3; j++){
					if (icc[j]>0.9)opn = j;
				}
				if (opn == -1) break;
				if (cen == lmesh->FElist[nfn][opn]){
					nen = lmesh->FElist[nfn][(opn + 2) % 3];
					ncoef = 0.0;
					if (lmesh->Edges[nen][0] == lmesh->Faces[nfn][(opn + 2) % 3]) ncoef = 1.0;
				}
				else if (cen == lmesh->FElist[nfn][(opn + 2) % 3]){
					nen = lmesh->FElist[nfn][opn];
					ncoef = 0.0;
					if (lmesh->Edges[nen][0] == lmesh->Faces[nfn][(opn + 1) % 3]) ncoef = 1.0;
				}
				else {
					if (coefs.size() > 1 && (coefs[coefs.size() - 1] < 0.1 || coefs[coefs.size() - 1] > 0.9)){
						int ppn = lmesh->Edges[cen][1];
						if (coefs[coefs.size() - 1] > 0.9) ppn = lmesh->Edges[cen][0];
						int ppen = eset[coefs.size() - 2];
						int fpn = 1;
						if (ppn != lmesh->Faces[nfn][(opn + 1) % 3]) fpn = 2;
						if (lmesh->Edges[ppen][0] == ppn || lmesh->Edges[ppen][1] == ppn){
							if (fpn == 1) {
								nen = lmesh->FElist[nfn][(opn + 2) % 3];
								ncoef = 0.0;
								if (lmesh->Edges[nen][0] == lmesh->Faces[nfn][(opn + 2) % 3]) ncoef = 1.0;
							}
							else {
								nen = lmesh->FElist[nfn][opn];
								ncoef = 0.0;
								if (lmesh->Edges[nen][0] == lmesh->Faces[nfn][(opn + 1) % 3]) ncoef = 1.0;
							}
						}
						else {
							if (fpn == 1){
								nen = lmesh->FElist[nfn][opn];
								ncoef = 0.0;
								if (lmesh->Edges[nen][0] == lmesh->Faces[nfn][(opn + 1) % 3]) ncoef = 1.0;
							}
							else {
								nen = lmesh->FElist[nfn][(opn + 2) % 3];
								ncoef = 0.0;
								if (lmesh->Edges[nen][0] == lmesh->Faces[nfn][(opn + 2) % 3]) ncoef = 1.0;
							}
						}
					}
					else {
						float seki1 = 0.0;
						float seki2 = 0.0;
						for (int j = 0; j < 3; j++){
							seki1 += dir[j] * (lmesh->Nodes[lmesh->Faces[nfn][(opn + 1) % 3]][j] - pos[j]);
							seki2 += dir[j] * (lmesh->Nodes[lmesh->Faces[nfn][(opn + 2) % 3]][j] - pos[j]);
						}
						if (seki1 > seki2){
							nen = lmesh->FElist[nfn][opn];
							ncoef = 0.0;
							if (lmesh->Edges[nen][0] == lmesh->Faces[nfn][(opn + 1) % 3]) ncoef = 1.0;
						}
						else {
							nen = lmesh->FElist[nfn][(opn + 2) % 3];
							ncoef = 0.0;
							if (lmesh->Edges[nen][0] == lmesh->Faces[nfn][(opn + 2) % 3]) ncoef = 1.0;
						}
					}
				}
			}
			else if (cstate == 3){
				float sekis[3], lengs[3];
				int maxN = 0;
				float maxseki;
				for (int j = 0; j < 3; j++){
					sekis[j] = 0.0;
					lengs[j] = 0.0;
					for (int k = 0; k < 3; k++){
						sekis[j] += dir[k] * (lmesh->Nodes[lmesh->Faces[nfn][j]][k] - pos[k]);
						lengs[j] += pd[k] * (lmesh->Nodes[lmesh->Faces[nfn][j]][k] - pos[k]);
					}
					if (j == 0 || maxseki < sekis[j]){
						maxN = j;
						maxseki = sekis[j];
					}
				}
				if (lmesh->FElist[nfn][maxN] == cen){
					nen = lmesh->FElist[nfn][(maxN + 2) % 3];
					ncoef = 0.0;
					if (lmesh->Edges[nen][0] == lmesh->Faces[nfn][maxN])ncoef = 1.0;
				}
				else if (lmesh->FElist[nfn][(maxN + 2) % 3] == cen){
					nen = lmesh->FElist[nfn][maxN];
					ncoef = 0.0;
					if (lmesh->Edges[nen][0] == lmesh->Faces[nfn][maxN]) ncoef = 1.0;
				}
				else {
					if (fabs(lengs[(maxN + 1) % 3]) < fabs(lengs[(maxN + 2) % 3])){
						nen = lmesh->FElist[nfn][maxN];
						ncoef = 0.0;
						if (lmesh->Edges[nen][0] == lmesh->Faces[nfn][maxN]) ncoef = 1.0;
					}
					else {
						nen = lmesh->FElist[nfn][(maxN + 2) % 3];
						ncoef = 0.0;
						if (lmesh->Edges[nen][0] == lmesh->Faces[nfn][maxN])ncoef = 1.0;
					}
				}
			}
			else {
				int cpn = -1;
				for (int j = 0; j < 3; j++){
					if (icc[j] < 0.0 || icc[j] > 1.0) cpn = j;
				}
				if (cpn == -1) break;
				cpn = (cpn + 2) % 3;
				if (lmesh->FElist[nfn][(cpn + 2) % 3] == cen){
					nen = lmesh->FElist[nfn][cpn];
					ncoef = icc[cpn];
					if (lmesh->Edges[nen][0] == lmesh->Faces[nfn][cpn]) ncoef = 1.0 - ncoef;
				}
				else if (lmesh->FElist[nfn][cpn] == cen){
					nen = lmesh->FElist[nfn][(cpn + 2) % 3];
					ncoef = icc[(cpn + 2) % 3];
					if (lmesh->Edges[nen][1] == lmesh->Faces[nfn][cpn]) ncoef = 1.0 - ncoef;
				}
				else {
					break;
				}
			}

			//更新・終端チェック


			if (fset.size() > 10){
				if (nfn == fn) break;
				if (start_is_vertex){
					bool return_neib = false;
					for (int j = 0; j < lmesh->PFsize[startp]; j++){
						if (nfn == lmesh->PFlist[startp][j]) return_neib = true;
					}
					if (return_neib) break;
				}
			}
			if (fset.size() > lmesh->faceN / 2.0) break;

			if (cstate != 0) printf("%d %d %f\n", nfn, cen, ncoef);
			cfn = nfn;
			cen = nen;
			fset.push_back(cfn);
			eset.push_back(cen);
			coefs.push_back(ncoef);

			float pos1[3], pos2[3];
			float leng = 0.0;
			for (int j = 0; j < 3; j++){
				pos1[j] = ncoef*lmesh->Nodes[lmesh->Edges[cen][0]][j] + (1.0 - ncoef)*lmesh->Nodes[lmesh->Edges[cen][1]][j];
				if (eset.size() == 1){
					pos2[j] = pos[j];
				}
				else {
					pos2[j] = coefs[coefs.size() - 2] * lmesh->Nodes[lmesh->Edges[eset[eset.size() - 2]][0]][j] +
						(1.0 - coefs[coefs.size() - 2])*lmesh->Nodes[lmesh->Edges[eset[eset.size() - 2]][1]][j];
				}
				leng += pow(pos1[j] - pos2[j], 2);
			}
			leng = sqrt(leng);
			Mleng += leng;
			if (Mleng > maxMove){
				MoveSuccess = true;
				break;
			}

		}

		if (MoveSuccess){
			Cfset.push_back(fset);
			Ceset.push_back(eset);
			Ccoefs.push_back(coefs);
		}
	}
	if (Cfset.size() == 0) return;

	int cNum = -1;
	float maxDist;
	for (int i = 0; i < Cfset.size(); i++){
		float hleng = 0.0;
		float dist = 0.0;;
		float p1[3], p2[3];
		for (int j = 0; j < 3; j++){
			p1[j] = pos[j];
		}
		for (int j = 0; j < Ceset[i].size(); j++){
			int en = Ceset[i][j];
			float coef = Ccoefs[i][j];
			float leng = 0.0;
			float seki = 0.0;
			for (int k = 0; k < 3; k++){
				p2[k] = coef*lmesh->Nodes[lmesh->Edges[en][0]][k] + (1.0 - coef)*lmesh->Nodes[lmesh->Edges[en][1]][k];
				leng += pow(p1[k] - p2[k], 2);
				seki += dir[k] * ((p1[k] + p2[k]) / 2.0 - pos[k]);
			}
			leng = sqrt(leng);
			hleng += leng;
			for (int k = 0; k < 3; k++){
				p1[k] = p2[k];
			}
			if (hleng > maxMove || j == Ceset[i].size() - 1){
				dist = seki;
				break;
			}
		}
		if (i == 0 || dist > maxDist){
			cNum = i;
			maxDist = dist;
		}
	}

	pN = Cfset[cNum].size() + 1;
	float **nset;
	nset = new float*[pN];
	for (int i = 0; i < pN; i++){
		nset[i] = new float[3];
	}
	for (int i = 0; i < 3; i++){
		nset[0][i] = pos[i];
	}
	for (int i = 0; i < Cfset[cNum].size(); i++){
		int en = Ceset[cNum][i];
		for (int j = 0; j < 3; j++){
			nset[i + 1][j] = Ccoefs[cNum][i] * lmesh->Nodes[lmesh->Edges[en][0]][j] +
				(1.0 - Ccoefs[cNum][i]) *lmesh->Nodes[lmesh->Edges[en][1]][j];
		}
	}
	*pset = nset;

}