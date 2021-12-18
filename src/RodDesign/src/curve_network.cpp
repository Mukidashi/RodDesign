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
#include<Eigen/src/Eigenvalues/EigenSolver.h>
#include<Eigen/Eigenvalues>
#include <vector>
#include <map>
using namespace std;

#include<Eigen/SparseCore>
#include<Eigen/Sparse>
#include<Eigen/Dense>
#include<Eigen/src/SparseCholesky/SimplicialCholesky.h>
#include<Eigen/src/SparseLU/SparseLU.h>
using namespace Eigen;

#include "curve_network.h"
#include "basic_geometric_calculation.h"
#include "AABB_data_structure.h"
#include "loop_subdivision.h"
#include "Setting.h"


curve_on_mesh::curve_on_mesh(){

	nodeN = 0;
	CONNECT_LIST_EXIST = false;

}


curve_on_mesh::~curve_on_mesh(){

	printf("curve_on_mesh dest start\n");

	if (nodeN == 0) return;

	delete[] cNode;
	delete[] Mframe;
	delete[] TPframe;
	delete[] pCons;
	delete[] dCons;

	if (CONNECT_LIST_EXIST){
		for (int i = 0; i < nodeN; i++){
			if(Csize[i] > 0) delete[] Clist[i];
		}
		delete[] Clist;
		delete[] Csize;
	}

	nodeN = 0;

	printf("curve_on_mesh dest end\n");
}


void curve_on_mesh::setPointConstraint(std::vector<struct ConstraintData> pcSet){

	for (int i = 0; i < pcSet.size(); i++){
		int pn = pcSet[i].pn;
		pCons[pn].cType = 1;
		for (int j = 0; j < 3; j++){
			pCons[pn].pos[j] = pcSet[i].pos[j];
		}
	}

}


void  curve_on_mesh::setDirectionConstraint(std::vector<struct ConstraintData> dcSet){

	for (int i = 0; i < dcSet.size(); i++){
		int pn = dcSet[i].pn;
		dCons[pn].cType = 1;
		float leng = 0.0;
		for (int j = 0; j < 3; j++){
			dCons[pn].pos[j] = dcSet[i].pos[j];
			leng += dCons[pn].pos[j] * dCons[pn].pos[j];
		}
		leng = sqrt(leng);
		if (leng != 0.0){
			for (int j = 0; j < 3; j++){
				dCons[pn].pos[j] /= leng;
			}
		}
	}

}


//最小化のための準備
void curve_on_mesh::prepare_for_minimization(){
	
	pCons = new struct ConstraintData[nodeN];
	dCons = new struct ConstraintData[nodeN - 1];
	for (int i = 0; i < nodeN; i++){
		pCons[i].cType = -1;
		if (i != nodeN - 1) dCons[i].cType = -1;
	}



	//フレームの設定
	TPframe = new struct adaptive_frame[nodeN - 1];
	Mframe = new struct adaptive_frame[nodeN - 1];

	for (int i = 0; i < nodeN - 1; i++){
		float leng = 0.0;
		for (int j = 0; j < 3; j++){
			TPframe[i].d3[j] = cNode[i + 1].pos[j] - cNode[i].pos[j];
			leng += pow(cNode[i + 1].pos[j] - cNode[i].pos[j], 2);
		}
		leng = sqrt(leng);
		for (int j = 0; j < 3; j++){
			TPframe[i].d3[j] /= leng;
		}
		float duv1[2][3], duv2[2][3];
		float uvw1[3], uvw2[3];
		uvw1[1] = cNode[i + 1].uv[0];
		uvw1[2] = cNode[i + 1].uv[1];
		uvw1[0] = 1.0 - uvw1[1] - uvw1[2];
		uvw2[1] = cNode[i].uv[0];
		uvw2[2] = cNode[i].uv[1];
		uvw2[0] = 1.0 - uvw2[1] - uvw2[2];
		lmesh->evaluate_derivative_of_box_spline_on_triangle(cNode[i + 1].fn, uvw1, duv1);
		lmesh->evaluate_derivative_of_box_spline_on_triangle(cNode[i].fn, uvw2, duv2);

		float nd1[3], nd2[3];
		float seki1 = 0.0;
		float seki2 = 0.0;
		for (int j = 0; j < 3; j++){
			nd1[j] = duv1[0][(j + 1) % 3] * duv1[1][(j + 2) % 3] - duv1[0][(j + 2) % 3] * duv1[1][(j + 1) % 3];
			nd2[j] = duv2[0][(j + 1) % 3] * duv2[1][(j + 2) % 3] - duv2[0][(j + 2) % 3] * duv2[1][(j + 1) % 3];
			seki1 += nd1[j] * TPframe[i].d3[j];
			seki2 += nd2[j] * TPframe[i].d3[j];
		}
		float leng1 = 0.0;
		float leng2 = 0.0;
		for (int j = 0; j < 3; j++){
			nd1[j] -= seki1*TPframe[i].d3[j];
			nd2[j] -= seki2*TPframe[i].d3[j];
			leng1 += nd1[j] * nd1[j];
			leng2 += nd2[j] * nd2[j];
		}
		leng1 = sqrt(leng1);
		leng2 = sqrt(leng2);
		if (leng1 < 1.0e-7) leng1 = 1.0e-7;
		if (leng2 < 1.0e-7) leng2 = 1.0e-7;
		seki1 = 0.0;
		float oseki[3];
		for (int j = 0; j < 3; j++){
			nd1[j] /= leng1;
			nd2[j] /= leng2;
			seki1 += nd1[j] * nd2[j];
			oseki[j] = nd1[(j + 1) % 3] * nd2[(j + 2) % 3] - nd1[(j + 2) % 3] * nd2[(j + 1) % 3];
		}
		seki2 = 0.0;
		for (int j = 0; j < 3; j++){
			seki2 += oseki[j] * TPframe[i].d3[j];
		}

		if (seki1 > 1.0) seki1 = 1.0;
		if (seki1 < -1.0) seki1 = -1.0;
		seki1 = acos(seki1);
		if (seki2 < 0.0) seki1 = -seki1;
		float quat[4], Rot[3][3];
		quat[0] = cos(seki1 / 4.0);
		for (int j = 0; j < 3; j++){
			quat[j + 1] = sin(seki1 / 4.0)*TPframe[i].d3[j];
		}
		derive_rot_mat_from_quaternion(quat, Rot);
		for (int j = 0; j < 3; j++){
			TPframe[i].d2[j] = 0.0;
			for (int k = 0; k < 3; k++){
				TPframe[i].d2[j] += Rot[j][k] * nd1[k];
			}
		}
		leng = 0.0;
		for (int j = 0; j < 3; j++){
			TPframe[i].d1[j] = TPframe[i].d2[(j + 1) % 3] * TPframe[i].d3[(j + 2) % 3]
				- TPframe[i].d2[(j + 2) % 3] * TPframe[i].d3[(j + 1) % 3];
			leng += pow(TPframe[i].d1[j], 2);
		}
		leng = sqrt(leng);
		if (leng < 1.0e-7) leng = 1.0e-7;
		for (int j = 0; j < 3; j++){
			TPframe[i].d1[j] /= leng;
		}

		for (int j = 0; j < 3; j++){
			Mframe[i].d1[j] = TPframe[i].d1[j];
			Mframe[i].d2[j] = TPframe[i].d2[j];
			Mframe[i].d3[j] = TPframe[i].d3[j];
		}
	}


}


//rod energyの最小化
bool curve_on_mesh::minimize_discrete_rod_energy(int maxIter){


	Eigen::VectorXd Fint(2 * nodeN), dVW(2 * nodeN), VW(2 * nodeN),Xk(3*nodeN);
	Eigen::VectorXi fn(nodeN);
	Eigen::VectorXd X0(3 * nodeN), VW0(2 * nodeN);
	Eigen::VectorXi f0(nodeN);

	update_frame_for_minimization();

	float(*duv)[2][3];
	duv = new float[nodeN][2][3];

	for (int i = 0; i < nodeN; i++){
		VW(2 * i) = cNode[i].uv[0];
		VW(2 * i + 1) = cNode[i].uv[1];
		fn(i) = cNode[i].fn;
		VW0(2 * i) = cNode[i].uv[0];
		VW0(2 * i + 1) = cNode[i].uv[1];
		f0(i) = cNode[i].fn;
		for (int j = 0; j < 3; j++){
			Xk(3 * i + j) = cNode[i].pos[j];
			X0(3 * i + j) = cNode[i].pos[j];
		}
	}

	int iterN = 0;
	bool is_converged = true;
	while (1){
		Eigen::SparseMatrix<double> Kmat(2 * nodeN, 2 * nodeN);
		std::vector<Eigen::Triplet<double>> kmlist;
		Fint = Eigen::VectorXd::Zero(2 * nodeN);

#pragma omp parallel for
		for (int i = 0; i < nodeN; i++){
			float uvw[3];
			uvw[1] = cNode[i].uv[0];
			uvw[2] = cNode[i].uv[1];
			uvw[0] = 1.0 - uvw[1] - uvw[2];
			lmesh->evaluate_derivative_of_box_spline_on_triangle(cNode[i].fn, uvw, duv[i]);
		}
		for (int i = 0; i < nodeN; i++){
			if (pCons[i].cType != -1){
				float leng = 0.0;
				for (int j = 0; j < 3; j++){
					leng += pow(cNode[i].pos[j] - pCons[i].pos[j], 2);
				}
				leng = sqrt(leng);
				if (leng < 1.0e-3* sParam->Wid) pCons[i].cType = 2;
			}
		}
		for (int i = 0; i < nodeN - 1; i++){
			if (dCons[i].cType != -1){
				float ed[3];
				float seki1 = 0.0;
				float edh = 0.0;
				for (int j = 0; j < 3; j++){
					seki1 += dCons[i].pos[j] * Mframe[i].d2[j];
					ed[j] = cNode[i + 1].pos[j] - cNode[i].pos[j];
					edh += ed[j] * ed[j];
				}
				edh = sqrt(edh);
				if (edh < 1.0e-10) edh = 1.0e-10;
				float seki2 = 0.0;
				float dh = 0.0;
				for (int j = 0; j < 3; j++){
					dCons[i].pos[j] -= seki1*Mframe[i].d2[j];
					seki2 += dCons[i].pos[j] * ed[j];
					dh += dCons[i].pos[j] * dCons[i].pos[j];
				}
				dh = sqrt(dh);
				for (int j = 0; j < 3; j++){
					if (dh > 1.0e-10 && seki2 > 0) dCons[i].pos[j] /= dh;
					else if (dh > 1.0e-10) dCons[i].pos[j] /= -dh;
					else dCons[i].pos[j] = ed[j] / edh;
				}

			}
		}



		//ベクトル・行列への代入
		double Across = sParam->Hei*sParam->Wid;
		double ks = sParam->Es*Across;
		double kt = sParam->Gt*Across * (sParam->Wid*sParam->Wid + sParam->Hei*sParam->Hei) / 12.0;
		double kben[2][2];
		kben[0][0] = sParam->Eb11*Across*sParam->Wid*sParam->Wid / 12.0;
		kben[1][1] = sParam->Eb22 *Across*sParam->Hei*sParam->Hei / 12.0;
		kben[0][1] = sParam->Eb12*Across*sParam->Hei*sParam->Wid / 12.0;
		kben[1][0] = sParam->Eb12*Across*sParam->Wid*sParam->Hei / 12.0;

		float Ebend = 0.0;
		float Ediff = 0.0;
		float Econs = 0.0;

		eav = 0.0;
		for (int i = 0; i < nodeN - 1; i++){
			float leng = 0.0;
			for (int j = 0; j < 3; j++){
				leng += pow(cNode[i + 1].pos[j] - cNode[i].pos[j], 2);
			}
			leng = sqrt(leng);
			eav += leng / (float)(nodeN - 1);
		}
		//テスト
		ks *= 0.01;
		for (int i = 0; i < nodeN - 1; i++){
			float leng = 0.0;
			float ed[3];
			for (int j = 0; j < 3; j++){
				ed[j] = cNode[i + 1].pos[j] - cNode[i].pos[j];
				leng += pow(cNode[i + 1].pos[j] - cNode[i].pos[j], 2);
			}
			leng = sqrt(leng);
			Ediff += 0.5*ks*pow(leng - eav, 2) / eav;

			float Fi[6], Ki[6][6];
			for (int j = 0; j < 3; j++){
				float fval = ks / eav / leng*(leng - eav)*ed[j];
				Fi[j] = fval;
				Fi[3 + j] = -fval;
				for (int k = 0; k < 3; k++){
					float val = ks / pow(leng, 3)*ed[j] * ed[k];
					if (j == k) val += ks / eav*(1.0 - eav / leng);
					Ki[j][k] = val;
					Ki[3 + j][k] = -val;
					Ki[j][3 + k] = -val;
					Ki[3 + j][k + 3] = val;
				}
			}
			for (int j = 0; j < 2; j++){
				for (int k = 0; k < 2; k++){
					float val = 0.0;
					for (int l = 0; l < 3; l++){
						val += Fi[3 * j + l] * duv[i + j][k][l];
					}
					Fint(2 * (i + j) + k) += val;
					for (int l = 0; l < 2; l++){
						for (int m = 0; m < 2; m++){
							float val = 0.0;
							for (int n = 0; n < 3; n++){
								for (int o = 0; o < 3; o++){
									val += Ki[3 * j + n][3 * k + o] * duv[i + j][l][n] * duv[i + k][m][o];
								}
							}
							kmlist.push_back(Triplet<double>(2 * (i + j) + l, 2 * (i + k) + m, val));
						}
					}
				}
			}
		}
		ks *= 100.0;

		for (int i = 1; i < nodeN - 1; i++){
			float kb[3];
			float kai = 1.0;
			for (int j = 0; j < 3; j++){
				kb[j] = TPframe[i - 1].d3[(j + 1) % 3] * TPframe[i].d3[(j + 2) % 3]
					- TPframe[i - 1].d3[(j + 2) % 3] * TPframe[i].d3[(j + 1) % 3];
				kai += TPframe[i - 1].d3[j] * TPframe[i].d3[j];
			}
			for (int j = 0; j < 3; j++){
				kb[j] = 2.0*kb[j] / kai;
			}
			float ki[2];
			ki[0] = 0.0;
			ki[1] = 0.0;
			float leng1 = 0.0;
			float leng2 = 0.0;
			for (int j = 0; j < 3; j++){
				ki[0] += 0.5*(Mframe[i - 1].d2[j] + Mframe[i].d2[j])*kb[j];
				ki[1] -= 0.5*(Mframe[i - 1].d1[j] + Mframe[i].d1[j])*kb[j];
				leng1 += pow(cNode[i].pos[j] - cNode[i - 1].pos[j], 2);
				leng2 += pow(cNode[i + 1].pos[j] - cNode[i].pos[j], 2);
			}
			leng1 = sqrt(leng1);
			leng2 = sqrt(leng2);
			float leni = 0.5*(leng1 + leng2);

			float gke[2][2][3], hkee[2][4][3][3];
			derivative_of_bend(i, leng1, leng2, gke, hkee);

			for (int j = 0; j < 2; j++){
				for (int k = 0; k < 2; k++){
					Ebend += 0.5*ki[j] * kben[j][k] * ki[k] / leni;
				}
			}
			//Ebend += 0.5*ki[0] * kben[0][0] * ki[0] / leni;
			
			float Fi[9],Ki[9][9];
			float Fi2[9], Ki2[9][9];
			for (int j = 0; j < 9; j++){
				Fi[j] = 0.0;
				Fi2[j] = 0.0;
				for (int k = 0; k < 9; k++){
					Ki[j][k] = 0.0;
					Ki2[j][k] = 0.0;
				}
			}
			for (int j = 0; j < 2; j++){
				for (int k = 0; k < 3; k++){
					float val = 0.0;
					for (int l = 0; l < 2; l++){
						for (int m = 0; m < 2; m++){
							val += gke[l][j][k] * kben[l][m] * ki[m] / leni;
						}
					}
					Fi[3 * j + k] += val;
					Fi[3 * (j + 1) + k] += -val;
					
					Fi2[3 * j + k] += val;
					Fi2[3 * (j + 1) + k] += -val;
				}
				for (int k = 0; k < 2; k++){
					for (int l = 0; l < 3; l++){
						for (int m = 0; m < 3; m++){
							float val = (gke[0][j][l] * kben[0][0] * gke[0][k][m] + hkee[0][2 * j + k][l][m] * kben[0][0] * ki[0]) / leni;
							val = 0.0;
							for (int n = 0; n < 2; n++){
								for (int o = 0; o < 2; o++){
									val += (gke[n][j][l] * kben[n][o] * gke[o][k][m] + hkee[n][2 * j + k][l][m] * kben[n][o] * ki[o]) / leni;
								}
							}
							Ki[3 * j + l][3 * k + m] += val;
							Ki[3 * (j + 1) + l][3 * k + m] -= val;
							Ki[3 * j + l][3 * (k + 1) + m] -= val;
							Ki[3 * (j + 1) + l][3 * (k + 1) + m] += val;

							Ki2[3 * j + l][3 * k + m] += val;
							Ki2[3 * (j + 1) + l][3 * k + m] -= val;
							Ki2[3 * j + l][3 * (k + 1) + m] -= val;
							Ki2[3 * (j + 1) + l][3 * (k + 1) + m] += val;
						}
					}
				}
			}

			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 2; k++){
					float fval = 0.0;
					for (int l = 0; l < 3; l++){
						fval += Fi[3 * j + l] * duv[i - 1 + j][k][l];
					}
					Fint(2 * (i - 1 + j) + k) += fval;
				}
				for (int k = 0; k < 3; k++){
					for (int l = 0; l < 2; l++){
						for (int m = 0; m < 2; m++){
							float kval = 0.0;
							for (int n = 0; n < 3; n++){
								for (int o = 0; o < 3; o++){
									kval += Ki[3 * j + n][3 * k + o] * duv[i - 1 + j][l][n] * duv[i - 1 + k][m][o];
								}
							}
							kmlist.push_back(Eigen::Triplet<double>(2 * (i - 1 + j) + l, 2 * (i - 1 + k) + m, kval));
						}
					}
				}
			}

		}

		//位置拘束
		float mu1 = 0.01*ks;
		for (int i = 0; i < nodeN; i++){
			if (pCons[i].cType == 1){
				float Fi[3],Ki[3][3];
				for (int j = 0; j < 3; j++){
					Fi[j] = -mu1*(cNode[i].pos[j] - pCons[i].pos[j]);
					for (int k = 0; k < 3; k++){
						if(j == k) Ki[j][k] = mu1;
						else Ki[j][k] = 0.0;
					}
					Econs += 0.5*mu1*pow(cNode[i].pos[j] - pCons[i].pos[j], 2);
				}
				for (int j = 0; j < 2; j++){
					float val = 0.0;
					for (int k = 0; k < 3; k++){
						val += Fi[k] * duv[i][j][k];
					}
					Fint(2 * i + j) += val;
					for (int k = 0; k < 2; k++){
						val = 0.0;
						for (int l = 0; l < 3; l++){
							for (int m = 0; m < 3; m++){
								val += Ki[l][m] * duv[i][j][l] * duv[i][k][m];
							}
						}
						kmlist.push_back(Triplet<double>(2 * i + j, 2 * i + k, val));
					}
				}
			}
		}

		//エッジ方向の拘束
		float mu2 = ks*0.01*eav*(float)(nodeN - 1.0);
		float Edcons = 0.0;
		for (int i = 0; i < nodeN - 1; i++){
			if (dCons[i].cType != -1){
				float ed[3];
				float leng = 0.0;
				float seki = 0.0;
				for (int j = 0; j < 3; j++){
					ed[j] = cNode[i + 1].pos[j] - cNode[i].pos[j];
					leng += ed[j] * ed[j];
					seki += ed[j] * dCons[i].pos[j];
				}
				leng = sqrt(leng);
				float theta = seki;
				if (leng != 0) theta /= leng;
				if (theta > 1.0) theta = 1.0;
				if (theta < -1.0) theta = -1.0;
				theta = acos(theta);
				Econs += 0.5*mu2*theta*theta;
				Edcons += 0.5*mu2*theta*theta;
				/*Econs += 0.5*mu2*pow(seki / leng - 1, 2);
				Edcons += 0.5*mu2*pow(seki / leng - 1, 2);*/

				float gT[3],KT[3][3];
				for (int j = 0; j < 3; j++){
					if (fabs(sin(theta)) > 1.0e-6) gT[j] = -(dCons[i].pos[j] - seki / leng / leng*ed[j]) / leng / sin(theta);
					else gT[j] = 0.0;
				}
				for (int j = 0; j < 3; j++){
					for (int k = 0; k < 3; k++){
						if (fabs(sin(theta))>1.0e-6){
							KT[j][k] = -cos(theta) / sin(theta)*gT[j] * gT[k];
							KT[j][k] += (dCons[i].pos[j] * ed[k] + ed[j] * dCons[i].pos[k] - 3.0*seki / leng / leng*ed[j] * ed[k]) / sin(theta) / pow(leng, 3);
							if (j == k) KT[j][k] += seki / sin(theta) / pow(leng, 3);
						}
						else{
							KT[j][k] = 0.0;
						}
					}
				}
				float Fi[6], Ki[6][6];
				for (int j = 0; j < 3; j++){
					Fi[j] = mu2*theta*gT[j];
					Fi[3 + j] = -mu2*theta*gT[j];
					for (int k = 0; k < 3; k++){
						float kval = mu2*gT[j] * gT[k] + mu2*theta*KT[j][k];
						Ki[j][k] = kval;
						Ki[3 + j][k] = -kval;
						Ki[j][3 + k] = -kval;
						Ki[3 + j][3 + k] = kval;
					}
				}
				/*for (int j = 0; j < 3; j++){
					Fi[j] = mu2 / leng*(seki / leng - 1.0)*(dCons[i].pos[j] - seki / leng / leng*ed[j]);
					Fi[3 + j] = -mu2 / leng*(seki / leng - 1.0)*(dCons[i].pos[j] - seki / leng / leng*ed[j]);
					for (int k = 0; k < 3; k++){
						float kval = mu2 / leng / leng*(dCons[i].pos[j] - seki / leng / leng*ed[j])*(dCons[i].pos[k] - seki / leng / leng*ed[k]);
						kval -= mu2 / pow(leng, 3)*(seki / leng - 1.0)*(dCons[i].pos[j] * ed[k] + ed[j] * dCons[i].pos[k] - 3.0*seki / leng / leng*ed[j] * ed[k]);
						if (j == k) kval -= mu2 / pow(leng, 3)*(seki / leng - 1.0)*seki;
					}
				}*/
				
				for (int j = 0; j < 2; j++){
					for (int k = 0; k < 2; k++){
						float fval = 0.0;
						for (int l = 0; l < 3; l++){
							fval += Fi[3 * j + l] * duv[i + j][k][l];
						}
						Fint(2 * (i + j) + k) += fval;

						for (int l = 0; l < 2; l++){
							for (int m = 0; m < 2; m++){
								float kval = 0.0;
								for (int n = 0; n < 3; n++){
									for (int o = 0; o < 3; o++){
										kval += Ki[3 * j + n][3 * k + o] * duv[i + j][l][n] * duv[i + k][m][o];
									}
								}
								kmlist.push_back(Triplet<double>(2 * (i + j) + l, 2 * (i + k) + m, kval));
							}
						}
					}
				}

			}
		}


		Kmat.setFromTriplets(kmlist.begin(), kmlist.end());

		//拘束の反映
		for (int i = 0; i < nodeN; i++){
			if (pCons[i].cType == 2){
				for (int j = 0; j < 2; j++){
					for (SparseMatrix<double>::InnerIterator it(Kmat, 2 * i + j); it; ++it){
						if (it.row() == 2 * i + j) it.valueRef() = 1.0;
						else it.valueRef() = 0.0;
					}
					Fint(2 * i + j) = 0.0;
				}
			}
		}
		for (int i = 0; i < 2 * nodeN; i++){
			for (SparseMatrix<double>::InnerIterator it(Kmat, i); it; ++it){
				int nN = (int)(it.row() / 2.0);
				if (pCons[nN].cType == 2){
					if (it.row() != i) it.valueRef() = 0.0;
				}
			}
		}



		//線形方程式の求解
		/*SparseLU<SparseMatrix<double>> sLU;
		sLU.compute(Kmat);
		if (sLU.info() != Success){
			printf("sparse LU failed\n");

		}
		dVW = sLU.solve(Fint);*/
		SimplicialLDLT<SparseMatrix<double>> LDLT;
		LDLT.compute(Kmat);
		if (LDLT.info() != Success){
			printf("LDLT failed\n");
		}
		dVW = LDLT.solve(Fint);

		//Line searchによるステップサイズの決定
		float alph = bisection_line_search(fn, VW, dVW);
		//float alph = line_search(fn, VW, dVW);

		//更新・収束判定
		float move = 0.0;
		float max_move = 0.0;
		float fval = 0.0;
#pragma omp parallel for
		for (int i = 0; i < nodeN; i++){
			move_on_mesh(cNode[i].fn, cNode[i].uv[0], cNode[i].uv[1], alph*dVW(2 * i), alph*dVW(2 * i + 1));
			float uvw[3], pos[3];
			uvw[1] = cNode[i].uv[0];
			uvw[2] = cNode[i].uv[1];
			uvw[0] = 1.0 - uvw[1] - uvw[2];
			lmesh->evaluate_box_spline_on_triangle(cNode[i].fn, uvw, pos);
			for (int j = 0; j < 3; j++){
				Xk(3 * i + j) = pos[j];
			}
		}

		for (int i = 0; i < nodeN;i++){
			float leng = 0.0;
			for (int j = 0; j < 3; j++){
				leng += pow(cNode[i].pos[j] - Xk(3*i+j), 2);
			}
			move += leng;
			if (max_move < sqrt(leng)) max_move = sqrt(leng);
			fval += Fint(2 * i)*Fint(2 * i) + Fint(2 * i + 1)*Fint(2 * i + 1);

			fn(i) = cNode[i].fn;
			VW(2 * i) = cNode[i].uv[0];
			VW(2 * i + 1) = cNode[i].uv[1];
			for (int j = 0; j < 3; j++){
				cNode[i].pos[j] = Xk(3 * i + j);
			}
		}
		//printf("%d E:%f %f f:%f move:%f\n",iterN,Ebend,Ediff,sqrt(fval),max_move);
		//printf("%d E:%f(%f %f %f,%f) move:%f\n",iterN,Ebend+Ediff+Econs,Ebend,Ediff,Econs,Edcons,max_move);

		if (max_move < 1.0e-3*sParam->Wid) break;

		update_frame_for_minimization();

		iterN++;
		if (iterN == maxIter) {
			is_converged = false;
			break;
		}
	}

	delete[] duv;

	if (!is_converged){
		for (int i = 0; i < nodeN; i++){
			cNode[i].fn = f0(i);
			for (int j = 0; j < 3; j++){
				cNode[i].pos[j] = X0(3 * i + j);
				if (j != 2) cNode[i].uv[j] = VW0(2 * i + j);
			}
		}
		update_frame_for_minimization();

		return false;
	}
	return true;


}


//フレーム・位置の更新（最小化各反復後）
void curve_on_mesh::update_frame_for_minimization(bool Initialize){

#pragma omp parallel for
	for (int i = 0; i < nodeN - 1; i++){
		float leng = 0.0;
		for (int j = 0; j < 3; j++){
			Mframe[i].d3[j] = cNode[i + 1].pos[j] - cNode[i].pos[j];
			leng += pow(Mframe[i].d3[j], 2);
		}
		leng = sqrt(leng);
		if (leng < 1.0e-7) leng = 1.0e-7;
		for (int j = 0; j < 3; j++){
			Mframe[i].d3[j] /= leng;
		}

		float uvw1[3], uvw2[3];
		float duv1[2][3], duv2[2][3];
		uvw1[1] = cNode[i + 1].uv[0];
		uvw1[2] = cNode[i+1].uv[1];
		uvw1[0] = 1.0 - uvw1[1] - uvw1[2];
		uvw2[1] = cNode[i].uv[0];
		uvw2[2] = cNode[i].uv[1];
		uvw2[0] = 1.0 - uvw2[1] - uvw2[2];
		lmesh->evaluate_derivative_of_box_spline_on_triangle(cNode[i+1].fn, uvw1, duv1);
		lmesh->evaluate_derivative_of_box_spline_on_triangle(cNode[i].fn, uvw2, duv2);

		float nd1[3], nd2[3];
		float seki1 = 0.0;
		float seki2 = 0.0;
		for (int j = 0; j < 3; j++){
			nd1[j] = duv1[0][(j + 1) % 3] * duv1[1][(j + 2) % 3]
				- duv1[0][(j + 2) % 3] * duv1[1][(j + 1) % 3];
			nd2[j] = duv2[0][(j + 1) % 3] * duv2[1][(j + 2) % 3]
				- duv2[0][(j + 2) % 3] * duv2[1][(j + 1) % 3];
			seki1 += nd1[j] * Mframe[i].d3[j];
			seki2 += nd2[j] * Mframe[i].d3[j];
		}
		float leng1 = 0.0;
		float leng2 = 0.0;
		for (int j = 0; j < 3; j++){
			nd1[j] -= seki1*Mframe[i].d3[j];
			nd2[j] -= seki2*Mframe[i].d3[j];
			leng1 += nd1[j] * nd1[j];
			leng2 += nd2[j] * nd2[j];
		}
		leng1 = sqrt(leng1);
		leng2 = sqrt(leng2);
		if (leng1 < 1.0e-7) leng1 = 1.0e-7;
		if (leng2 < 1.0e-7) leng2 = 1.0e-7;
		float theta = 0.0;
		for (int j = 0; j < 3; j++){
			nd1[j] /= leng1;
			nd2[j] /= leng2;
			theta += nd1[j] * nd2[j];
		}
		float seki = 0.0;
		for (int j = 0; j<3; j++){
			seki += (nd1[(j + 1) % 3] * nd2[(j + 2) % 3] - nd1[(j + 1) % 3] * nd2[(j + 2) % 3])*Mframe[i].d3[j];
		}
		if (theta > 1.0) theta = 1.0;
		if (theta < -1.0) theta = -1.0;
		theta = acos(theta);
		if (seki < 0.0) theta *= -1.0;

		float quat[4], Rot[3][3];
		quat[0] = cos(theta / 4.0);
		for (int j = 0; j < 3; j++){
			quat[1 + j] = sin(theta / 4.0)*Mframe[i].d3[j];
		}
		derive_rot_mat_from_quaternion(quat, Rot);
		for (int j = 0; j < 3; j++){
			Mframe[i].d2[j] = 0.0;
			for (int k = 0; k < 3; k++){
				Mframe[i].d2[j] += Rot[j][k] * nd1[k];
			}
		}
		leng = 0.0;
		for (int j = 0; j < 3; j++){
			Mframe[i].d1[j] = Mframe[i].d2[(j + 1) % 3] * Mframe[i].d3[(j + 2) % 3]
				- Mframe[i].d2[(j + 2) % 3] * Mframe[i].d3[(j + 1) % 3];
			leng += pow(Mframe[i].d1[j], 2);
		}
		leng = sqrt(leng);
		if (leng < 1.0e-7) leng = 1.0e-7;
		for (int j = 0; j < 3; j++){
			Mframe[i].d1[j] /= leng;
		}

		if (Initialize){
			for (int j = 0; j < 3; j++){
				TPframe[i].d1[j] = Mframe[i].d1[j];
				TPframe[i].d2[j] = Mframe[i].d2[j];
				TPframe[i].d3[j] = Mframe[i].d3[j];
			}
			continue;
		}

		float axis[3];
		float nseki = 0.0;
		float axh = 0.0;
		for (int j = 0; j < 3; j++){
			axis[j] = TPframe[i].d3[(j + 1) % 3] * Mframe[i].d3[(j + 2) % 3]
				- TPframe[i].d3[(j + 2) % 3] * Mframe[i].d3[(j + 1) % 3];
			nseki += TPframe[i].d3[j] * Mframe[i].d3[j];
			axh += axis[j] * axis[j];
		}
		axh = sqrt(axh);

		
		if (axh > sin(1.0e-5)){
			float quat[4];
			if (nseki > 1.0) nseki = 1.0;
			if (nseki < -1.0) nseki = -1.0;
			nseki = acos(nseki);
			quat[0] = cos(nseki / 2.0);
			for (int j = 0; j < 3; j++){
				quat[1 + j] = sin(nseki/2.0)*axis[j] / axh;
			}
			derive_rot_mat_from_quaternion(quat, Rot);
		}
		else {
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 3; k++){
					if (j == k && nseki >0.0)Rot[j][k] = 1.0;
					else if (j == k) Rot[j][k] = -1.0;
					else Rot[j][k] = 0.0;
				}
			}
		}

		float d1[3], d2[3];
		for (int j = 0; j < 3; j++){
			d1[j] = 0.0;
			d2[j] = 0.0;
			for (int k = 0; k < 3; k++){
				d1[j] += Rot[j][k] * TPframe[i].d1[k];
				d2[j] += Rot[j][k] * TPframe[i].d2[k];
			}
		}
		for (int j = 0; j < 3; j++){
			TPframe[i].d1[j] = d1[j];
			TPframe[i].d2[j] = d2[j];
			TPframe[i].d3[j] = Mframe[i].d3[j];
		}

	}


}


//微分の計算
void curve_on_mesh::derivative_of_bend(int in, float leng1, float leng2, float gke[2][2][3], float hkee[2][4][3][3]){

	float kb[3], tit[3];
	float kai = 1.0;
	for (int i = 0; i < 3; i++){
		kb[i] = Mframe[in - 1].d3[(i + 1) % 3] * Mframe[in].d3[(i + 2) % 3]
			- Mframe[in - 1].d3[(i + 2) % 3] * Mframe[in].d3[(i + 1) % 3];
		kai += Mframe[in - 1].d3[i] * Mframe[in].d3[i];
	}
	for (int i = 0; i < 3; i++){
		kb[i] = 2.0*kb[i] / kai;
		tit[i] = (Mframe[in - 1].d3[i] + Mframe[in].d3[i]) / kai;
	}

	float ki[2], tid1[3], tid2[3];
	ki[0] = 0.0;
	ki[1] = 0.0;
	for (int i = 0; i < 3; i++){
		ki[0] += 0.5*(Mframe[in - 1].d2[i] + Mframe[in].d2[i])*kb[i];
		ki[1] -= 0.5*(Mframe[in - 1].d1[i] + Mframe[in].d1[i])*kb[i];
		tid1[i] = (Mframe[in - 1].d1[i] + Mframe[in].d1[i]) / kai;
		tid2[i] = (Mframe[in - 1].d2[i] + Mframe[in].d2[i]) / kai;
	}
	float oseki[2][2][3];
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 3; j++){
			oseki[i][0][j] = Mframe[in - 1 + i].d3[(j + 1) % 3] * tid1[(j + 2) % 3]
				- Mframe[in - 1 + i].d3[(j + 2) % 3] * tid1[(j + 1) % 3];
			oseki[i][1][j] = Mframe[in - 1 + i].d3[(j + 1) % 3] * tid2[(j + 2) % 3]
				- Mframe[in - 1 + i].d3[(j + 2) % 3] * tid2[(j + 1) % 3];
		}
	}

	for (int i = 0; i < 3; i++){
		gke[0][0][i] = (-ki[0] * tit[i] + oseki[1][1][i]) / leng1;
		gke[0][1][i] = (-ki[0] * tit[i] - oseki[0][1][i]) / leng2;
		gke[1][0][i] = (-ki[1] * tit[i] - oseki[1][0][i]) / leng1;
		gke[1][1][i] = (-ki[1] * tit[i] + oseki[0][0][i]) / leng2;
		for (int j = 0; j < 3; j++){
			hkee[0][0][i][j] = (2.0*ki[0] * tit[i] * tit[j] - oseki[1][1][i] * tit[j] - tit[i] * oseki[1][1][j]) / leng1 / leng1
				+ ki[0] * Mframe[in - 1].d3[i] * Mframe[in - 1].d3[j] / kai / leng1 / leng1;
			hkee[0][1][i][j] = (2.0*ki[0] * tit[i] * tit[j] - oseki[1][1][i] * tit[j] + tit[i] * oseki[0][1][j]) / leng1 / leng2
				- ki[0] * Mframe[in - 1].d3[i] * Mframe[in].d3[j] / kai / leng1 / leng2;
			hkee[0][2][i][j] = (2.0*ki[0] * tit[i] * tit[j] - tit[i] * oseki[1][1][j] + oseki[0][1][i] * tit[j]) / leng1 / leng2
				- ki[0] * Mframe[in].d3[i] * Mframe[in - 1].d3[j] / kai / leng1 / leng2;
			hkee[0][3][i][j] = (2.0*ki[0] * tit[i] * tit[j] + oseki[0][1][i] * tit[j] + tit[i] * oseki[0][1][j]) / leng2 / leng2
				+ ki[0] * Mframe[in].d3[i] * Mframe[in].d3[j] / kai / leng2 / leng2;
			hkee[1][0][i][j] = (2.0*ki[1] * tit[i] * tit[j] + oseki[1][0][i] * tit[j] + tit[i] * oseki[1][0][j]) / leng1 / leng1
				+ ki[1] * Mframe[in - 1].d3[i] * Mframe[in - 1].d3[j] / kai / leng1 / leng1;
			hkee[1][1][i][j] = (2.0*ki[1] * tit[i] * tit[j] + oseki[1][0][i] * tit[j] - tit[i] * oseki[0][0][j]) / leng1 / leng2
				- ki[1] * Mframe[in - 1].d3[i] * Mframe[in].d3[j] / kai / leng1 / leng2;
			hkee[1][2][i][j] = (2.0*ki[1] * tit[i] * tit[j] + tit[i] * oseki[1][0][j] - oseki[0][0][i] * tit[j]) / leng1 / leng2
				- ki[1] * Mframe[in].d3[i] * Mframe[in - 1].d3[j] / kai / leng1 / leng2;
			hkee[1][3][i][j] = (2.0*ki[1] * tit[i] * tit[j] - oseki[0][0][i] * tit[j] - tit[i] * oseki[0][0][j]) / leng2 / leng2
				+ ki[1] * Mframe[in].d3[i] * Mframe[in].d3[j] / kai / leng2 / leng2;
			if (i == j){
				hkee[0][0][i][j] -= ki[0] / kai / leng1 / leng1;
				hkee[0][1][i][j] -= ki[0] / kai / leng1 / leng2;
				hkee[0][2][i][j] -= ki[0] / kai / leng1 / leng2;
				hkee[0][3][i][j] -= ki[0] / kai / leng2 / leng2;
				hkee[1][0][i][j] -= ki[1] / kai / leng1 / leng1;
				hkee[1][1][i][j] -= ki[1] / kai / leng1 / leng2;
				hkee[1][2][i][j] -= ki[1] / kai / leng1 / leng2;
				hkee[1][3][i][j] -= ki[1] / kai / leng2 / leng2;
			}
			else {
				float sig = 1.0;
				if (j == (i + 1) % 3) sig = -1.0;
				hkee[0][1][i][j] += -sig*tid2[3 - i - j] / leng1 / leng2;
				hkee[0][2][i][j] += sig*tid2[3 - i - j] / leng1 / leng2;
				hkee[1][1][i][j] += sig*tid1[3 - i - j] / leng1 / leng2;
				hkee[1][2][i][j] += -sig*tid1[3 - i - j] / leng1 / leng2;
			}
			//論文記載のヘッシアン
			hkee[0][0][i][j] += (kb[i] * Mframe[in - 1].d2[j] + Mframe[in - 1].d2[i] * kb[j]) / 4.0 / leng1 / leng1;
			hkee[0][3][i][j] += (kb[i] * Mframe[in].d2[j] + Mframe[in].d2[i] * kb[j]) / 4.0 / leng2 / leng2;
			hkee[1][0][i][j] -= (kb[i] * Mframe[in - 1].d1[j] + Mframe[in - 1].d1[i] * kb[j]) / 4.0 / leng1 / leng1;
			hkee[1][3][i][j] -= (kb[i] * Mframe[in].d1[j] + Mframe[in].d1[i] * kb[j]) / 4.0 / leng2 / leng2;
		}
	}


}


//メッシュ上の移動
void curve_on_mesh::move_on_mesh(int &fn, float &Vk, float &Wk, float dv, float dw){


	if (Vk < 0.0 || Wk < 0.0 || Vk + Wk > 1.0){
		printf("point is outside at start @move on mesh\n");
		return;
	}

	float uvw[3], duvw[3];
	uvw[0] = 1.0 - Vk - Wk;
	uvw[1] = Vk;
	uvw[2] = Wk;
	duvw[0] = -dv - dw;
	duvw[1] = dv;
	duvw[2] = dw;

	Eigen::MatrixXf SS(2, 2);
	Eigen::VectorXf Sv(2),dv2(2);

	int iterN = 0;
	int cfn = fn;
	while (1){

		//侵出点の特定
		int eN1 = -1;
		float intK, edK1;
		if (fabs(duvw[2]) > 1.0e-5){
			float kv = -uvw[2] / duvw[2];
			if (kv > 1.0e-7 && kv <= 1.0){
				float inv = uvw[1] + kv*duvw[1];
				if (inv >= 0.0 && inv <= 1.0){
					eN1 = 0;
					intK = kv;
					edK1 = 1.0 - inv;
				}
			}
		}
		if (fabs(duvw[0]) > 1.0e-5){
			float kv = -uvw[0] / duvw[0];
			if (kv > 1.0e-7 && kv <= 1.0){
				float inv = uvw[1] + kv*duvw[1];
				if (inv >= 0.0 && inv <= 1.0){
					eN1 = 1;
					intK = kv;
					edK1 = inv;
				}
			}
		}
		if (fabs(duvw[1]) > 1.0e-5){
			float kv = -uvw[1] / duvw[1];
			if (kv > 1.0e-7 && kv <= 1.0){
				float inw = uvw[2] + kv*duvw[2];
				if (inw >= 0.0 && inw <= 1.0){
					eN1 = 2;
					intK = kv;
					edK1 = inw;
				}
			}
		}
		if (eN1 == -1){
			fn = cfn;
			Vk = uvw[1] + duvw[1];
			Wk = uvw[2] + duvw[2];
			if (Vk > 1.0) Vk = 1.0;
			else if (Vk < 0.0)Vk = 0.0;
			if (Wk > 1.0) Wk = 1.0;
			else if (Wk < 0.0) Wk = 0.0;
			if (Vk + Wk > 1.0){
				float wa = Vk + Wk;
				Vk = Vk / wa;
				Wk = Wk / wa;
			}
			break;
		}

		//侵出面の特定
		float edK = edK1;
		int cen = lmesh->FElist[cfn][eN1];
		if (lmesh->Edges[cen][0] != lmesh->Faces[cfn][eN1])edK = 1.0 - edK1;
		int nfn = lmesh->EFlist[cen][0];
		if (nfn == cfn) nfn = lmesh->EFlist[cen][1];
		if (nfn == -1){
			printf("boundary is found@move on mesh\n");
			fn = cfn;
			Vk = uvw[1] + intK*duvw[1];
			Wk = uvw[2] + intK*duvw[2];
			break;
		}

		int eN2 = -1;
		for (int i = 0; i < 3; i++){
			if (lmesh->Edges[cen][0] == lmesh->Faces[nfn][i] && lmesh->Edges[cen][1] == lmesh->Faces[nfn][(i + 1) % 3]) eN2 = i;
			if (lmesh->Edges[cen][1] == lmesh->Faces[nfn][i] && lmesh->Edges[cen][0] == lmesh->Faces[nfn][(i + 1) % 3]) eN2 = i;
		}
		float edK2 = edK;
		if (lmesh->Edges[cen][0] != lmesh->Faces[nfn][eN2]) edK2 = 1.0 - edK;

		float uvw1[3], uvw2[3];
		uvw1[eN1] = edK1;
		uvw1[(eN1 + 1) % 3] = 1.0 - edK1;
		uvw1[(eN1 + 2) % 3] = 0.0;
		uvw2[eN2] = edK2;
		uvw2[(eN2 + 1) % 3] = 1.0 - edK2;
		uvw2[(eN2 + 2) % 3] = 0.0;

		float dbs1[2][3], dbs2[2][3];
		lmesh->evaluate_derivative_of_box_spline_on_triangle(cfn, uvw1, dbs1);
		lmesh->evaluate_derivative_of_box_spline_on_triangle(nfn, uvw2, dbs2);

		for (int i = 0; i < 3; i++){
			duvw[i] = (1.0 - intK)*duvw[i];
		}

		for (int i = 0; i < 2; i++){
			Sv(i) = 0.0;
			for (int j = 0; j < 2; j++){
				SS(i, j) = 0.0;
				for (int k = 0; k < 3; k++){
					SS(i, j) += dbs2[i][k] * dbs2[j][k];
					Sv(i) += dbs2[i][k] * dbs1[j][k] * duvw[1 + j];
				}
			}
		}
		dv2 = SS.inverse()*Sv;

		//fn・dvの更新
		cfn = nfn;
		for (int i = 0; i < 3; i++){
			uvw[i] = uvw2[i];
		}
		duvw[0] = - dv2(0) - dv2(1);
		duvw[1] = dv2(0);
		duvw[2] = dv2(1);
		
		//チェック用
		iterN++;
		if (iterN == lmesh->faceN) break;
	}


}


//メッシュ上の移動（移動方向の更新有）
void curve_on_mesh::move_on_mesh2(int &fn, float &Vk, float &Wk, float &dv, float &dw){


	if (Vk < 0.0 || Wk < 0.0 || Vk + Wk > 1.0){
		printf("point is outside at start @move on mesh\n");
		return;
	}

	float uvw[3], duvw[3];
	uvw[0] = 1.0 - Vk - Wk;
	uvw[1] = Vk;
	uvw[2] = Wk;
	duvw[0] = -dv - dw;
	duvw[1] = dv;
	duvw[2] = dw;

	Eigen::MatrixXf SS(2, 2);
	Eigen::VectorXf Sv(2), dv2(2);

	int iterN = 0;
	int cfn = fn;
	while (1){

		//侵出点の特定
		int eN1 = -1;
		float intK, edK1;
		if (fabs(duvw[2]) > 1.0e-5){
			float kv = -uvw[2] / duvw[2];
			if (kv > 1.0e-7 && kv <= 1.0){
				float inv = uvw[1] + kv*duvw[1];
				if (inv >= 0.0 && inv <= 1.0){
					eN1 = 0;
					intK = kv;
					edK1 = 1.0 - inv;
				}
			}
		}
		if (fabs(duvw[0]) > 1.0e-5){
			float kv = -uvw[0] / duvw[0];
			if (kv > 1.0e-7 && kv <= 1.0){
				float inv = uvw[1] + kv*duvw[1];
				if (inv >= 0.0 && inv <= 1.0){
					eN1 = 1;
					intK = kv;
					edK1 = inv;
				}
			}
		}
		if (fabs(duvw[1]) > 1.0e-5){
			float kv = -uvw[1] / duvw[1];
			if (kv > 1.0e-7 && kv <= 1.0){
				float inw = uvw[2] + kv*duvw[2];
				if (inw >= 0.0 && inw <= 1.0){
					eN1 = 2;
					intK = kv;
					edK1 = inw;
				}
			}
		}
		if (eN1 == -1){
			fn = cfn;
			Vk = uvw[1] + duvw[1];
			Wk = uvw[2] + duvw[2];
			if (Vk > 1.0) Vk = 1.0;
			else if (Vk < 0.0)Vk = 0.0;
			if (Wk > 1.0) Wk = 1.0;
			else if (Wk < 0.0) Wk = 0.0;
			if (Vk + Wk > 1.0){
				float wa = Vk + Wk;
				Vk = Vk / wa;
				Wk = Wk / wa;
			}
			dv = duvw[1];
			dw = duvw[2];
			break;
		}

		//侵出面の特定
		float edK = edK1;
		int cen = lmesh->FElist[cfn][eN1];
		if (lmesh->Edges[cen][0] != lmesh->Faces[cfn][eN1])edK = 1.0 - edK1;
		int nfn = lmesh->EFlist[cen][0];
		if (nfn == cfn) nfn = lmesh->EFlist[cen][1];
		if (nfn == -1){
			printf("boundary is found@move on mesh\n");
			fn = cfn;
			Vk = uvw[1] + intK*duvw[1];
			Wk = uvw[2] + intK*duvw[2];
			dv = duvw[1];
			dw = duvw[2];
			break;
		}

		int eN2 = -1;
		for (int i = 0; i < 3; i++){
			if (lmesh->Edges[cen][0] == lmesh->Faces[nfn][i] && lmesh->Edges[cen][1] == lmesh->Faces[nfn][(i + 1) % 3]) eN2 = i;
			if (lmesh->Edges[cen][1] == lmesh->Faces[nfn][i] && lmesh->Edges[cen][0] == lmesh->Faces[nfn][(i + 1) % 3]) eN2 = i;
		}
		float edK2 = edK;
		if (lmesh->Edges[cen][0] != lmesh->Faces[nfn][eN2]) edK2 = 1.0 - edK;

		float uvw1[3], uvw2[3];
		uvw1[eN1] = edK1;
		uvw1[(eN1 + 1) % 3] = 1.0 - edK1;
		uvw1[(eN1 + 2) % 3] = 0.0;
		uvw2[eN2] = edK2;
		uvw2[(eN2 + 1) % 3] = 1.0 - edK2;
		uvw2[(eN2 + 2) % 3] = 0.0;

		float dbs1[2][3], dbs2[2][3];
		lmesh->evaluate_derivative_of_box_spline_on_triangle(cfn, uvw1, dbs1);
		lmesh->evaluate_derivative_of_box_spline_on_triangle(nfn, uvw2, dbs2);

		for (int i = 0; i < 3; i++){
			duvw[i] = (1.0 - intK)*duvw[i];
		}

		for (int i = 0; i < 2; i++){
			Sv(i) = 0.0;
			for (int j = 0; j < 2; j++){
				SS(i, j) = 0.0;
				for (int k = 0; k < 3; k++){
					SS(i, j) += dbs2[i][k] * dbs2[j][k];
					Sv(i) += dbs2[i][k] * dbs1[j][k] * duvw[1 + j];
				}
			}
		}
		dv2 = SS.inverse()*Sv;

		//fn・dvの更新
		cfn = nfn;
		for (int i = 0; i < 3; i++){
			uvw[i] = uvw2[i];
		}
		duvw[0] = -dv2(0) - dv2(1);
		duvw[1] = dv2(0);
		duvw[2] = dv2(1);

		////チェック用
		//iterN++;
		//if (iterN == 10 * lmesh->faceN) {
		//	printf();
		//	break;
		//}
	}


}


//二分法によるline search
float curve_on_mesh::bisection_line_search(Eigen::VectorXi f0, Eigen::VectorXd VW0, Eigen::VectorXd &dVW){

	float Eval,dE;
	float alph = 0.0;
	curve_energy(f0,VW0,dVW, alph, Eval,dE);
	//printf("  initial val: %f\n",Eval);


	//line searchループ
	alph = 1.0;
	int iterN = 0;
	bool is_converge = true;
	while (1){
		float Ek, dE;
		curve_energy(f0, VW0, dVW, alph, Ek, dE);
		//printf("%d val:%f\n",iterN+1,Ek);
		
		if (Ek < Eval) break;
		alph /= 2.0;
		iterN++;
		if (alph < 1.0e-3) {
			//printf("line search not converge\n");
			is_converge = false;
			break;
		}
	}

	if (!is_converge){
		Eigen::VectorXd dVW2;
		dVW2 = Eigen::VectorXd::Zero(2*nodeN);
		discent_direction_of_curve_energy(f0, VW0, dVW2);
		for (int i = 0; i < nodeN; i++){
			float odh = sqrt(dVW(2 * i)*dVW(2 * i) + dVW(2 * i + 1)*dVW(2 * i + 1));
			float dh = sqrt(dVW2(2 * i)*dVW2(2 * i) + dVW2(2 * i + 1)*dVW2(2 * i + 1));
			if (dh > 1.0e-7){
				dVW2(2 * i) = dVW2(2 * i) / dh*odh;
				dVW2(2 * i + 1) = dVW2(2 * i + 1) / dh*odh;
			}
			else {
				dVW2(2 * i) = -dVW(2 * i);
				dVW2(2 * i + 1) = -dVW(2 * i + 1);
			}
		}

		dVW = dVW2;
		alph = 1.0;
	    iterN = 0;
		is_converge = true;
		while (1){
			float Ek, dE;
			curve_energy(f0, VW0, dVW, alph, Ek, dE);
			//printf("%d val:%f\n",iterN+1,Ek);

			if (Ek < Eval) break;
			alph /= 2.0;
			iterN++;
			if (alph < 1.0e-4) {
				//printf("line search not converge for discent direction\n");
				is_converge = false;
				break;
			}
		}
		if (!is_converge) alph = 0.0;
	}

	return alph;
}


//line search
float curve_on_mesh::line_search(Eigen::VectorXi fn, Eigen::VectorXd VW, Eigen::VectorXd &dVW){

	float Phik, Phik1, dPhik, dPhik1;
	float Phi0;
	float dPhi0;
	//curve_energy(fn, VW, dVW, 1.0, Phi0, dPhi0);
	//printf("full %f %f\n", Phi0, dPhi0);
	curve_energy(fn, VW, dVW, 0.0, Phi0, dPhi0);
	if (dPhi0 > 0) {
		printf("delX is not discent direction %f %f\n", Phi0, dPhi0);
		float tP, tdP;
		Eigen::VectorXd dVW2;
		dVW2 = Eigen::VectorXd::Zero(2 * nodeN);
		discent_direction_of_curve_energy(fn, VW, dVW2);
		for (int i = 0; i < nodeN; i++){
			float odh = sqrt(dVW(2 * i)*dVW(2 * i) + dVW(2 * i + 1)*dVW(2 * i + 1));
			float dh = sqrt(dVW2(2 * i)*dVW2(2 * i) + dVW2(2 * i + 1)*dVW2(2 * i + 1));
			if (dh > 1.0e-7){
				dVW2(2 * i) = dVW2(2 * i) / dh*odh;
				dVW2(2 * i + 1) = dVW2(2 * i + 1) / dh*odh;
			}
			else {
				dVW2(2 * i) = -dVW(2 * i);
				dVW2(2 * i + 1) = -dVW(2 * i + 1);
			}
		}
		
		for (int i = 0; i < 2*nodeN; i++){
			dVW(i) = -dVW(i);
			dVW(i) = dVW2(i);
		}
		curve_energy(fn, VW, dVW, 0.0, Phi0, dPhi0);
	}
	//printf("init %f %f\n", Phi0, dPhi0);
	float al_max = 5.0;
	float alph = 1.0;
	float alp1 = 0.0;
	float c1 = 1.0e-4;
	float c2 = 0.9;
	Phik1 = Phi0;
	dPhik1 = dPhi0;
	int iterN = 5;
	while (1){
		curve_energy(fn, VW, dVW, alph, Phik, dPhik);

		bool scheck = true;
		if (Phik > Phi0 + c1*alph*dPhi0 || (iterN > 0 && Phik > Phik1)){
			float al_lo, al_hi, alpj;
			float Ph_lo, Ph_hi, dPh_lo, dPh_hi;
			al_lo = alp1;
			al_hi = alph;
			Ph_lo = Phik1;
			Ph_hi = Phik;
			dPh_lo = dPhik1;
			dPh_hi = dPhik;

			int innerN = 0;
			bool do_success = false;
			float Phij, dPhij;
			while (1){
				float d1 = dPh_lo + dPh_hi - 3.0*(Ph_lo - Ph_hi) / (al_lo - al_hi);
				float d2 = sqrt(d1*d1 - dPh_lo*dPh_hi);
				if (al_lo > al_hi) d2 *= -1.0;
				alpj = al_hi - (al_hi - al_lo)*(dPh_hi + d2 - d1) / (dPh_hi - dPh_lo + 2.0*d2);
				if ((al_lo - alpj)*(al_hi - alpj) > 0 || fabs(al_lo - alpj) < 0.01*fabs(al_lo - al_hi)
					|| fabs(al_hi - alpj) < 0.01*fabs(al_lo - al_hi) || __isnan(alpj)){
					alpj = (al_lo + al_hi) / 2.0;
				}
				//printf(" alph:%f %f %f %f %f\n", alpj, d1, d2, al_lo, al_hi);

				curve_energy(fn, VW, dVW, alpj, Phij, dPhij);
				//printf("%f(%f) %f(%f)-%f(%f)\n", alpj, Phij, al_lo, Ph_lo, al_hi, Ph_hi);

				//printf("  %d %f %f %f\n", innerN, alpj, Phij, dPhij);
				if (Phij > Phi0 + c1*alpj*dPhi0 || Phij >= Ph_lo){
					al_hi = alpj;
					Ph_hi = Phij;
					dPh_hi = dPhij;
				}
				else {
					//if (fabs(dPhij) <= -c2*dPhi0){
					if (dPhij >= c2*dPhi0){
						alph = alpj;
						do_success = true;
						break;
					}
					if (dPhij*(al_hi - al_lo) >= 0){
						al_hi = al_lo;
						Ph_hi = Ph_lo;
						dPh_hi = dPh_lo;
					}
					al_lo = alpj;
					Ph_lo = Phij;
					dPh_lo = dPhij;
				}

				innerN++;
				if (innerN > 10) break;
				if (fabs(al_lo - al_hi) < 1.0e-5*alph){
					//printf("zoom in line search failed\n");
					alph = (al_lo + al_hi) / 2.0;
					break;
				}
			}

			if (do_success) break;
			if (Phij < Phi0){
				alph = alpj;
				break;
			}
			al_max = alpj;
			alph = al_max*0.1;
			scheck = false;
		}

		if (scheck && fabs(dPhik) <= -c2*dPhi0){
			break;
		}
		if (scheck && dPhik >= 0){
			float al_lo, al_hi, alpj;
			float Ph_lo, Ph_hi, dPh_lo, dPh_hi;
			al_hi = alp1;
			al_lo = alph;
			Ph_hi = Phik1;
			Ph_lo = Phik;
			dPh_hi = dPhik1;
			dPh_lo = dPhik;

			int innerN = 0;
			bool do_success = false;
			float Phij, dPhij;
			while (1){
				float d1 = dPh_lo + dPh_hi - 3.0*(Ph_lo - Ph_hi) / (al_lo - al_hi);
				float d2 = sqrt(d1*d1 - dPh_lo*dPh_hi);
				if (al_lo > al_hi) d2 *= -1.0;
				alpj = al_hi - (al_hi - al_lo)*(dPh_hi + d2 - d1) / (dPh_hi - dPh_lo + 2.0*d2);
				if ((al_lo - alpj)*(al_hi - alpj) > 0 || fabs(al_lo - alpj) < 0.01*fabs(al_lo - al_hi)
					|| fabs(al_hi - alpj) < 0.01*fabs(al_lo - al_hi)){
					alpj = (al_lo + al_hi) / 2.0;
				}

				curve_energy(fn, VW, dVW, alpj, Phij, dPhij);
				//printf("  %d %f %f %f\n", innerN, alpj, Phij, dPhij);
				if (Phij > Phi0 + c1*alpj*dPhi0 || Phij >= Ph_lo){
					al_hi = alpj;
					Ph_hi = Phij;
					dPh_hi = dPhij;
				}
				else {
					//if (fabs(dPhij) <= -c2*dPhi0){
					if (dPhij >= c2*dPhi0){
						alph = alpj;
						do_success = true;
						break;
					}
					if (dPhij*(al_hi - al_lo) >= 0){
						al_hi = al_lo;
						Ph_hi = Ph_lo;
						dPh_hi = dPh_lo;
					}
					al_lo = alpj;
					Ph_lo = Phij;
					dPh_lo = dPhij;
				}

				innerN++;
				if (innerN >= 10) break;
				if (fabs(al_lo - al_hi) < 1.0e-5*alph){
					//printf("zoom in line search failed\n");
					alph = (al_lo + al_hi) / 2.0;
					break;
				}
			}
			if (do_success) break;
			if (Phij < Phi0){
				alph = alpj;
				break;
			}
			al_max = alpj;
			alph = 0.5*al_max;
		}

		if (1.2*alph >= al_max || al_max < 1.0e-5) break;
		alp1 = alph;
		alph *= 1.2;

		Phik1 = Phik;
		dPhik1 = dPhik;
		iterN++;
	}


	curve_energy(fn, VW, dVW, alph, Phi0, dPhi0);
	//printf(" final%f %f(%f)\n", Phi0, dPhi0,alph);

	return alph;

}


//エネルギーの計算
void curve_on_mesh::curve_energy(Eigen::VectorXi f0, Eigen::VectorXd VW, Eigen::VectorXd dVW, float alph, float &Eval,float &dE){

	float(*duv)[2][3];
	duv = new float[nodeN][2][3];

	float (*vwdr)[2];
	vwdr = new float[nodeN][2];

	Eigen::VectorXf Xk(3 * nodeN);
#pragma omp parallel for
	for (int i = 0; i < nodeN; i++){
		int fn = f0(i);
		float uvw[3],pos[3];
		uvw[1] = VW(2 * i);
		uvw[2] = VW(2 * i + 1);
		float dire[2];
		dire[0] = alph*dVW(2 * i);
		dire[1] = alph*dVW(2 * i + 1);
		move_on_mesh2(fn, uvw[1], uvw[2], dire[0], dire[1]);
		if (fabs(alph) > 1.0e-7){
			float odh = sqrt(dVW(2 * i)*dVW(2 * i) + dVW(2 * i + 1)*dVW(2 * i + 1));
			float dh = sqrt(dire[0] * dire[0] + dire[1] * dire[1]);
			if (dh > 1.0e-7){
				vwdr[i][0] = dire[0] / dh*odh;
				vwdr[i][1] = dire[1] / dh*odh;
			}
			else {
				vwdr[i][0] = odh / sqrt(2.0);
				vwdr[i][1] = odh / sqrt(2.0);
			}
		}
		else {
			vwdr[i][0] = dVW(2 * i);
			vwdr[i][1] = dVW(2 * i + 1);
		}
		
		uvw[0] = 1.0 - uvw[1] - uvw[2];
		lmesh->evaluate_box_spline_on_triangle(fn, uvw, pos);
		for (int j = 0; j < 3; j++){
			Xk(3 * i + j) = pos[j];
		}
		lmesh->evaluate_derivative_of_box_spline_on_triangle(fn, uvw, duv[i]);
	}

	struct adaptive_frame *mframe;
	mframe = new struct adaptive_frame[nodeN - 1];
#pragma omp parallel for
	for (int i = 0; i < nodeN - 1; i++){
		float leng = 0.0;
		for (int j = 0; j < 3; j++){
			mframe[i].d3[j] = Xk(3 * i + 3 + j) - Xk(3 * i + j);
			leng += pow(Xk(3 * i + 3 + j) - Xk(3 * i + j), 2);
		}
		leng = sqrt(leng);
		if (leng < 1.0e-7) leng = 1.0e-7;
		for (int j = 0; j < 3; j++){
			mframe[i].d3[j] /= leng;
		}

		float nd1[3], nd2[3];
		float seki1 = 0.0;
		float seki2 = 0.0;
		for (int j = 0; j < 3; j++){
			nd1[j] = duv[i + 1][0][(j + 1) % 3] * duv[i + 1][1][(j + 2) % 3]
				- duv[i + 1][0][(j + 2) % 3] * duv[i + 1][1][(j + 1) % 3];
			nd2[j] = duv[i][0][(j + 1) % 3] * duv[i][1][(j + 2) % 3]
				- duv[i][0][(j + 2) % 3] * duv[i][1][(j + 1) % 3];
			seki1 += nd1[j] * mframe[i].d3[j];
			seki2 += nd2[j] * mframe[i].d3[j];
		}
		float leng1 = 0.0;
		float leng2 = 0.0;
		for (int j = 0; j < 3; j++){
			nd1[j] -= seki1*mframe[i].d3[j];
			nd2[j] -= seki2*mframe[i].d3[j];
			leng1 += nd1[j] * nd1[j];
			leng2 += nd2[j] * nd2[j];
		}
		leng1 = sqrt(leng1);
		leng2 = sqrt(leng2);
		if (leng1 < 1.0e-7) leng1 = 1.0e-7;
		if (leng2 < 1.0e-7) leng2 = 1.0e-7;
		float theta = 0.0;
		float oseki[3];
		for (int j = 0; j < 3; j++){
			nd1[j] /= leng1;
			nd2[j] /= leng2;
			theta += nd1[j] * nd2[j];
		}
		float seki = 0.0;
		for (int j = 0; j<3; j++){
			seki += (nd1[(j + 1) % 3] * nd2[(j + 2) % 3] - nd1[(j + 1) % 3] * nd2[(j + 2) % 3])*mframe[i].d3[j];
		}
		if (theta > 1.0) theta = 1.0;
		if (theta < -1.0) theta = -1.0;
		theta = acos(theta);
		if (seki < 0.0) theta *= -1.0;
		
		float quat[4], Rot[3][3];
		quat[0] = cos(theta / 4.0);
		for (int j = 0; j < 3; j++){
			quat[1 + j] = sin(theta / 4.0)*mframe[i].d3[j];
		}
		derive_rot_mat_from_quaternion(quat, Rot);
		for (int j = 0; j < 3; j++){
			mframe[i].d2[j] = 0.0;
			for (int k = 0; k < 3; k++){
				mframe[i].d2[j] += Rot[j][k] * nd1[k];
			}
		}
		leng = 0.0;
		for (int j = 0; j < 3; j++){
			mframe[i].d1[j] = mframe[i].d2[(j + 1) % 3] * mframe[i].d3[(j + 2) % 3]
				- mframe[i].d2[(j + 2) % 3] * mframe[i].d3[(j + 1) % 3];
			leng += pow(mframe[i].d1[j], 2);
		}
		leng = sqrt(leng);
		if (leng < 1.0e-7) leng = 1.0e-7;
		for (int j = 0; j < 3; j++){
			mframe[i].d1[j] /= leng;
		}
	}



	//ベクトル・行列への代入
	double Across = sParam->Hei*sParam->Wid;
	double ks = sParam->Es*Across;
	double kt = sParam->Gt*Across * (sParam->Wid*sParam->Wid + sParam->Hei*sParam->Hei) / 12.0;
	double kben[2][2];
	kben[0][0] = sParam->Eb11*Across*sParam->Wid*sParam->Wid / 12.0;
	kben[1][1] = sParam->Eb22 *Across*sParam->Hei*sParam->Hei / 12.0;
	kben[0][1] = sParam->Eb12*Across*sParam->Hei*sParam->Wid / 12.0;
	kben[1][0] = sParam->Eb12*Across*sParam->Wid*sParam->Hei / 12.0;

	Eval = 0.0;
	dE = 0.0;
	ks *= 0.01;
	for (int i = 0; i < nodeN - 1; i++){
		float leng = 0.0;
		float ed[3];
		for (int j = 0; j < 3; j++){
			ed[j] = Xk(3 * i + 3 + j) - Xk(3 * i + j);
			leng += ed[j] * ed[j];
		}
		leng = sqrt(leng);
		Eval += 0.5*ks*pow(leng - eav, 2) / eav;

		float Fi[6], Ki[6][6];
		for (int j = 0; j < 3; j++){
			float fval = ks / eav / leng*(leng - eav)*ed[j];
			Fi[j] = -fval;
			Fi[3 + j] = fval;
		}
		for (int j = 0; j < 2; j++){
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 3; l++){
					dE += Fi[3 * j + l] * duv[i + j][k][l] * vwdr[i + j][k];
				}
			}
		}
	}
	ks *= 100.0;
	for (int i = 1; i < nodeN - 1; i++){
		float kb[3];
		float kai = 1.0;
		for (int j = 0; j < 3; j++){
			kb[j] = mframe[i - 1].d3[(j + 1) % 3] * mframe[i].d3[(j + 2) % 3]
				- mframe[i - 1].d3[(j + 2) % 3] * mframe[i].d3[(j + 1) % 3];
			kai += mframe[i - 1].d3[j] * mframe[i].d3[j];
		}
		for (int j = 0; j < 3; j++){
			kb[j] = 2.0*kb[j] / kai;
		}
		float ki[2];
		ki[0] = 0.0;
		ki[1] = 0.0;
		float leng1 = 0.0;
		float leng2 = 0.0;
		float tit[3], tid1[3], tid2[3];
		float e1[3], e2[3];
		for (int j = 0; j < 3; j++){
			ki[0] += 0.5*(mframe[i - 1].d2[j] + mframe[i].d2[j])*kb[j];
			ki[1] -= 0.5*(mframe[i - 1].d1[j] + mframe[i].d1[j])*kb[j];
			leng1 += pow(Xk(3 * i + j) - Xk(3 * i - 3 + j), 2);
			leng2 += pow(Xk(3 * i + 3 + j) - Xk(3 * i + j), 2);
			e1[j] = Xk(3 * i + j) - Xk(3 * i - 3 + j);
			e2[j] = Xk(3 * i + 3 + j) - Xk(3 * i + j);
			tit[j] = (mframe[i - 1].d3[j] + mframe[i].d3[j]) / kai;
			tid1[j] = (mframe[i - 1].d1[j] + mframe[i].d1[j]) / kai;
			tid2[j] = (mframe[i - 1].d2[j] + mframe[i].d2[j]) / kai;
		}
		leng1 = sqrt(leng1);
		leng2 = sqrt(leng2);
		float leni = 0.5*(leng1 + leng2);

		for (int j = 0; j < 2; j++){
			for (int k = 0; k < 2; k++){
				Eval += 0.5*ki[j] * kben[j][k] * ki[k] / leni;
			}
		}
		//Eval += 0.5*ki[0] * kben[0][0] * ki[0] / leni;
		//Eval += 0.5*kd*pow(leng1 - leng2, 2) / leni;

		float gke[2][2][3];
		for (int j = 0; j < 3; j++){
			gke[0][0][j] = (-ki[0] * tit[j] + mframe[i].d3[(j + 1) % 3] * tid2[(j + 2) % 3] - mframe[i].d3[(j + 2) % 3] * tid2[(j + 1) % 3]) / leng1;
			gke[0][1][j] = (-ki[0] * tit[j] - mframe[i - 1].d3[(j + 1) % 3] * tid2[(j + 2) % 3] + mframe[i - 1].d3[(j + 2) % 3] * tid2[(j + 1) % 3]) / leng2;
			gke[1][0][j] = (-ki[1] * tit[j] - mframe[i].d3[(j + 1) % 3] * tid1[(j + 2) % 3] + mframe[i].d3[(j + 2) % 3] * tid1[(j + 1) % 3]) / leng1;
			gke[1][1][j] = (-ki[1] * tit[j] + mframe[i - 1].d3[(j + 1) % 3] * tid1[(j + 2) % 3] - mframe[i - 1].d3[(j + 2) % 3] * tid1[(j + 1) % 3]) / leng2;
		}
		float gE[9];
		for (int j = 0; j < 9; j++){
			gE[j] = 0.0;
		}
		for (int j = 0; j < 2; j++){
			for (int k = 0; k < 3; k++){
				float val = 0.0;
				for (int l = 0; l < 2; l++){
					for (int m = 0; m < 2; m++){
						val += gke[l][j][k] * kben[l][m] * ki[m] / leni;
					}
				}
				gE[3 * j + k] += -val;
				gE[3 * (j + 1) + k] += val;
			}
		}
		/*for (int j = 0; j < 3; j++){
			gE[j] += -kd*(leng1 - leng2) / leng1*e1[j] / leni;
			gE[3 + j] += kd*(leng1 - leng2) / leng1*e1[j] / leni
				- kd*(leng2 - leng1) / leng2*e2[j] / leni;
			gE[6 + j] += kd*(leng2 - leng1) / leng2*e2[j] / leni;
		}*/

		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				for (int l = 0; l < 2; l++){
					dE += gE[3 * j + k] * duv[i - 1 + j][l][k] * vwdr[i - 1 + j][l];
				}
			}
		}

	}

	float mu1 = 0.01*ks;
	for (int i = 0; i < nodeN; i++){
		if (pCons[i].cType == 1){
			float gE[3];
			for (int j = 0; j < 3; j++){
				gE[j] = mu1*(Xk(3 * i + j) - pCons[i].pos[j], 2);
				Eval += 0.5*mu1*pow(Xk(3 * i + j) - pCons[i].pos[j], 2);
			}
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 2; k++){
					dE += gE[j] * duv[i][k][j] * vwdr[i][k];
				}
			}
		}
	}

	float mu2 = ks*0.01*eav*(float)(nodeN - 1.0);
	for (int i = 0; i < nodeN - 1; i++){
		if (dCons[i].cType != -1){
			float ed[3];
			float leng = 0.0;
			float seki = 0.0;
			for (int j = 0; j < 3; j++){
				ed[j] = Xk(3 * i + 3 + j) - Xk(3 * i + j);
				leng += ed[j] * ed[j];
				seki += ed[j] * dCons[i].pos[j];
			}
			leng = sqrt(leng);
			float theta = seki;
			if (leng != 0) theta /= leng;
			if (theta > 1.0) theta = 1.0;
			if (theta < -1.0) theta = -1.0;
			theta = acos(theta);
			Eval += 0.5*mu2*theta*theta;
			//Eval += 0.5*mu2*pow(seki / leng - 1.0, 2);

			float gE[6];
			for (int j = 0; j < 3; j++){
				float fval = 0.0;
				if (fabs(sin(theta)) > 1.0e-6) fval = -(dCons[i].pos[j] - seki / leng / leng*ed[j]) / leng / sin(theta);
				gE[j] = -mu2*theta*fval;
				gE[3 + j] = mu2*theta*fval;

				/*gE[j] = -mu2 / leng*(seki / leng - 1.0)*(dCons[i].pos[j] - seki / leng / leng*ed[j]);
				gE[3+j] = mu2 / leng*(seki / leng - 1.0)*(dCons[i].pos[j] - seki / leng / leng*ed[j]);*/
			}
			for (int j = 0; j < 2; j++){
				for (int k = 0; k < 2; k++){
					for (int l = 0; l < 3; l++){
						dE += gE[3 * j + l] * duv[i + j][k][l] * vwdr[i + j][k];
					}
				}
			}

		}
	}



	delete[] duv;
	delete[] mframe;
	delete[] vwdr;

}


void curve_on_mesh::discent_direction_of_curve_energy(Eigen::VectorXi f0, Eigen::VectorXd VW, Eigen::VectorXd &dVW){

	float(*duv)[2][3];
	duv = new float[nodeN][2][3];

	Eigen::VectorXf Xk(3 * nodeN);
#pragma omp parallel for
	for (int i = 0; i < nodeN; i++){
		int fn = f0(i);
		float uvw[3],pos[3];
		uvw[1] = VW(2 * i);
		uvw[2] = VW(2 * i + 1);
		uvw[0] = 1.0 - uvw[1] - uvw[2];
		lmesh->evaluate_box_spline_on_triangle(fn, uvw, pos);
		for (int j = 0; j < 3; j++){
			Xk(3 * i + j) = pos[j];
		}
		lmesh->evaluate_derivative_of_box_spline_on_triangle(fn, uvw, duv[i]);
	}

	struct adaptive_frame *mframe;
	mframe = new struct adaptive_frame[nodeN - 1];
	for (int i = 0; i < nodeN - 1; i++){
		float leng = 0.0;
		for (int j = 0; j < 3; j++){
			mframe[i].d3[j] = Xk(3 * i + 3 + j) - Xk(3 * i + j);
			leng += pow(Xk(3 * i + 3 + j) - Xk(3 * i + j), 2);
		}
		leng = sqrt(leng);
		if (leng < 1.0e-7) leng = 1.0e-7;
		for (int j = 0; j < 3; j++){
			mframe[i].d3[j] /= leng;
		}

		float nd1[3], nd2[3];
		float seki1 = 0.0;
		float seki2 = 0.0;
		for (int j = 0; j < 3; j++){
			nd1[j] = duv[i + 1][0][(j + 1) % 3] * duv[i + 1][1][(j + 2) % 3]
				- duv[i + 1][0][(j + 2) % 3] * duv[i + 1][1][(j + 1) % 3];
			nd2[j] = duv[i][0][(j + 1) % 3] * duv[i][1][(j + 2) % 3]
				- duv[i][0][(j + 2) % 3] * duv[i][1][(j + 1) % 3];
			seki1 += nd1[j] * mframe[i].d3[j];
			seki2 += nd2[j] * mframe[i].d3[j];
		}
		float leng1 = 0.0;
		float leng2 = 0.0;
		for (int j = 0; j < 3; j++){
			nd1[j] -= seki1*mframe[i].d3[j];
			nd2[j] -= seki2*mframe[i].d3[j];
			leng1 += nd1[j] * nd1[j];
			leng2 += nd2[j] * nd2[j];
		}
		leng1 = sqrt(leng1);
		leng2 = sqrt(leng2);
		if (leng1 < 1.0e-7) leng1 = 1.0e-7;
		if (leng2 < 1.0e-7) leng2 = 1.0e-7;
		float theta = 0.0;
		for (int j = 0; j < 3; j++){
			nd1[j] /= leng1;
			nd2[j] /= leng2;
			theta += nd1[j] * nd2[j];
		}
		float seki = 0.0;
		for (int j = 0; j<3; j++){
			seki += (nd1[(j + 1) % 3] * nd2[(j + 2) % 3] - nd1[(j + 2) % 3] * nd2[(j + 1) % 3])*mframe[i].d3[j];
		}
		if (theta > 1.0) theta = 1.0;
		if (theta < -1.0) theta = -1.0;
		theta = acos(theta);
		if (seki < 0.0) theta *= -1.0;

		float quat[4], Rot[3][3];
		quat[0] = cos(theta / 4.0);
		for (int j = 0; j < 3; j++){
			quat[1 + j] = sin(theta / 4.0)*mframe[i].d3[j];
		}
		derive_rot_mat_from_quaternion(quat, Rot);
		for (int j = 0; j < 3; j++){
			mframe[i].d2[j] = 0.0;
			for (int k = 0; k < 3; k++){
				mframe[i].d2[j] += Rot[j][k] * nd1[k];
			}
		}
		leng = 0.0;
		for (int j = 0; j < 3; j++){
			mframe[i].d1[j] = mframe[i].d2[(j + 1) % 3] * mframe[i].d3[(j + 2) % 3]
				- mframe[i].d2[(j + 2) % 3] * mframe[i].d3[(j + 1) % 3];
			leng += pow(mframe[i].d1[j], 2);
		}
		leng = sqrt(leng);
		if (leng < 1.0e-7) leng = 1.0e-7;
		for (int j = 0; j < 3; j++){
			mframe[i].d1[j] /= leng;
		}
	}

	//dVW = Eigen::VectorXd::Zero(3 * nodeN);

	//ベクトル・行列への代入
	double Across = sParam->Hei*sParam->Wid;
	double ks = sParam->Es*Across;
	double kt = sParam->Gt*Across * (sParam->Wid*sParam->Wid + sParam->Hei*sParam->Hei) / 12.0;
	double kben[2][2];
	kben[0][0] = sParam->Eb11*Across*sParam->Wid*sParam->Wid / 12.0;
	kben[1][1] = sParam->Eb22 *Across*sParam->Hei*sParam->Hei / 12.0;
	kben[0][1] = sParam->Eb12*Across*sParam->Hei*sParam->Wid / 12.0;
	kben[1][0] = sParam->Eb12*Across*sParam->Wid*sParam->Hei / 12.0;

	float kd = 1.0;
	ks *= 0.01;
	for (int i = 0; i < nodeN - 1; i++){
		float leng = 0.0;
		float ed[3];
		for (int j = 0; j < 3; j++){
			ed[j] = Xk(3 * i + 3 + j) - Xk(3 * i + j);
			leng += ed[j] * ed[j];
		}
		leng = sqrt(leng);

		float Fi[6], Ki[6][6];
		for (int j = 0; j < 3; j++){
			float fval = ks / eav / leng*(leng - eav)*ed[j];
			Fi[j] = -fval;
			Fi[3 + j] = fval;
		}
		for (int j = 0; j < 2; j++){
			for (int k = 0; k < 2; k++){ 
				for (int l = 0; l < 3; l++){
					dVW(2 * (i + j) + k) += -Fi[3 * j + l] * duv[i + j][k][l];
				}
			}
		}
	}
	ks *= 100.0;
	for (int i = 1; i < nodeN - 1; i++){
		float kb[3];
		float kai = 1.0;
		for (int j = 0; j < 3; j++){
			kb[j] = mframe[i - 1].d3[(j + 1) % 3] * mframe[i].d3[(j + 2) % 3]
				- mframe[i - 1].d3[(j + 2) % 3] * mframe[i].d3[(j + 1) % 3];
			kai += mframe[i - 1].d3[j] * mframe[i].d3[j];
		}
		for (int j = 0; j < 3; j++){
			kb[j] = 2.0*kb[j] / kai;
		}
		float ki[2];
		ki[0] = 0.0;
		ki[1] = 0.0;
		float leng1 = 0.0;
		float leng2 = 0.0;
		float tit[3], tid1[3], tid2[3];
		float e1[3], e2[3];
		for (int j = 0; j < 3; j++){
			ki[0] += 0.5*(mframe[i - 1].d2[j] + mframe[i].d2[j])*kb[j];
			ki[1] -= 0.5*(mframe[i - 1].d1[j] + mframe[i].d1[j])*kb[j];
			leng1 += pow(Xk(3 * i + j) - Xk(3 * i - 3 + j), 2);
			leng2 += pow(Xk(3 * i + 3 + j) - Xk(3 * i + j), 2);
			e1[j] = Xk(3 * i + j) - Xk(3 * i - 3 + j);
			e2[j] = Xk(3 * i + 3 + j) - Xk(3 * i + j);
			tit[j] = (mframe[i - 1].d3[j] + mframe[i].d3[j]) / kai;
			tid1[j] = (mframe[i - 1].d1[j] + mframe[i].d1[j]) / kai;
			tid2[j] = (mframe[i - 1].d2[j] + mframe[i].d2[j]) / kai;
		}
		leng1 = sqrt(leng1);
		leng2 = sqrt(leng2);
		float leni = 0.5*(leng1 + leng2);

		float gke[2][2][3];
		for (int j = 0; j < 3; j++){
			gke[0][0][j] = (-ki[0] * tit[j] + mframe[i].d3[(j + 1) % 3] * tid2[(j + 2) % 3] - mframe[i].d3[(j + 2) % 3] * tid2[(j + 1) % 3]) / leng1;
			gke[0][1][j] = (-ki[0] * tit[j] - mframe[i - 1].d3[(j + 1) % 3] * tid2[(j + 2) % 3] + mframe[i - 1].d3[(j + 2) % 3] * tid2[(j + 1) % 3]) / leng2;
			gke[1][0][j] = (-ki[1] * tit[j] - mframe[i].d3[(j + 1) % 3] * tid1[(j + 2) % 3] + mframe[i].d3[(j + 2) % 3] * tid1[(j + 1) % 3]) / leng1;
			gke[1][1][j] = (-ki[1] * tit[j] + mframe[i - 1].d3[(j + 1) % 3] * tid1[(j + 2) % 3] - mframe[i - 1].d3[(j + 2) % 3] * tid1[(j + 1) % 3]) / leng2;
		}
		float gE[9];
		for (int j = 0; j < 9; j++){
			gE[j] = 0.0;
		}
		for (int j = 0; j < 2; j++){
			for (int k = 0; k < 3; k++){
				float val = 0.0;
				for (int l = 0; l < 2; l++){
					for (int m = 0; m < 2; m++){
						val += gke[l][j][k] * kben[l][m] * ki[m] / leni;
					}
				}
				gE[3 * j + k] += -val;
				gE[3 * (j + 1) + k] += val;
			}
		}

		/*for (int j = 0; j < 3; j++){
			gE[j] += -kd*(leng1 - leng2) / leng1*e1[j] / leni;
			gE[3 + j] += kd*(leng1 - leng2) / leng1*e1[j] / leni
				- kd*(leng2 - leng1) / leng2*e2[j] / leni;
			gE[6 + j] += kd*(leng2 - leng1) / leng2*e2[j] / leni;
		}*/

		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 3; l++){
					dVW(2 * (i - 1 + j) + k) += -gE[3 * j + l] * duv[i - 1 + j][k][l];
				}
			}
		}

	}

	//位置拘束
	float mu1 = 0.01*ks;
	for (int i = 0; i < nodeN; i++){
		if (pCons[i].cType == 2){
			dVW(2 * i) = 0.0;
			dVW(2 * i + 1) = 0.0;
		}
		else if (pCons[i].cType == 1){
			float gE[3];
			for (int j = 0; j < 3; j++){
				gE[j] = mu1*(Xk(3 * i + j) - pCons[i].pos[j]);
			}
			for (int j = 0; j < 2; j++){
				for (int k = 0; k < 3; k++){
					dVW(2 * i + j) += -gE[k] * duv[i][j][k];
				}
			}
		}
	}

	float mu2 = ks*0.01*eav*(float)(nodeN - 1.0);
	for (int i = 0; i < nodeN - 1; i++){
		if (dCons[i].cType != -1){
			float ed[3];
			float leng = 0.0;
			float seki = 0.0;
			for (int j = 0; j < 3; j++){
				ed[j] = Xk(3 * i + 3 + j) - Xk(3 * i + j);
				leng += ed[j] * ed[j];
				seki += ed[j] * dCons[i].pos[j];
			}
			leng = sqrt(leng);
			float theta = seki;
			/*if (leng < 1.0e-7) leng = 1.0e-7;
			theta /= leng;*/
			if (leng != 0) theta /= leng;
			if (theta > 1.0) theta = 1.0;
			if (theta < -1.0) theta = -1.0;
			theta = acos(theta);

			float gE[6];
			for (int j = 0; j < 3; j++){
				float fval = 0.0;
				if (fabs(sin(theta)) > 1.0e-6) fval = -(dCons[i].pos[j] - seki / leng / leng*ed[j]) / leng / sin(theta);
				gE[j] = -mu2*theta*fval;
				gE[3 + j] = mu2*theta*fval;

				/*gE[j] = -mu2 / leng*(seki / leng - 1.0)*(dCons[i].pos[j] - seki / leng / leng*ed[j]);
				gE[3 + j] = mu2 / leng*(seki / leng - 1.0)*(dCons[i].pos[j] - seki / leng / leng*ed[j]);*/
			}
			for (int j = 0; j < 2; j++){
				for (int k = 0; k < 2; k++){
					for (int l = 0; l < 3; l++){
						dVW(2 * (i + j) + k) += -gE[3 * j + l] * duv[i + j][k][l];
					}
				}
			}

		}
	}

	delete[] duv;
	delete[] mframe;

}


void curve_on_mesh::remesh_curve_node(int selectC,int &selectN, int &selectE, class AABB_struct *AABB,std::vector<struct ConnectData> &connect){


	//短いエッジを消す
	while (1){
		float minL;
		int eNum = -1;
		for (int i = 0; i < nodeN - 1; i++){
			float leng = 0.0;
			for (int j = 0; j < 3; j++){
				leng += pow(cNode[i + 1].pos[j] - cNode[i].pos[j], 2);
			}
			leng = sqrt(leng);
			bool dcheck = true;
			if (pCons[i].cType != -1 && pCons[i + 1].cType != -1) dcheck = false;
			if (pCons[i].cType != -1 && i + 1 == selectN) dcheck = false;
			if (pCons[i + 1].cType != -1 && i == selectN) dcheck = false;
			if (i == 0 && (pCons[i + 1].cType != -1 || selectN == 1)) dcheck = false;
			if (i == nodeN - 2 && (pCons[i].cType != -1 || selectN == nodeN-2)) dcheck = false;
			if (i == selectE) dcheck = false;
			if (dcheck){
				if (eNum == -1 || minL > leng){
					minL = leng;
					eNum = i;
				}
			}
		}

		if (eNum == -1 || nodeN == 4) break;
		if (minL > 1.5*sParam->Wid) break;

		int delN = -1;
		if (pCons[eNum].cType != -1 || eNum == selectN) delN = eNum + 1;
		if (pCons[eNum + 1].cType != -1 || eNum + 1 == selectN) delN = eNum;
		if (eNum == 0) delN = 1;
		if (eNum == nodeN - 2) delN = nodeN - 2;
		if (delN == -1){
			float leng1 = 0.0;
			float leng2 = 0.0;
			for (int i = 0; i < 3; i++){
				leng1 += pow(cNode[eNum - 1].pos[i] - cNode[eNum].pos[i], 2);
				leng2 += pow(cNode[eNum + 1].pos[i] - cNode[eNum + 2].pos[i], 2);
			}
			if (leng1 > leng2) delN = eNum + 1;
			else delN = eNum;
		}

		if (selectN > delN) selectN--;
		if (selectE >= delN) selectE--;

		struct curve_node *NCnode;
		NCnode = new struct curve_node[nodeN - 1];
		struct ConstraintData *NPcon, *NDcon;
		NPcon = new struct ConstraintData[nodeN - 1];
		NDcon = new struct ConstraintData[nodeN - 2];
		int **NClist, *NCsize;
		if (CONNECT_LIST_EXIST){
			NClist = new int*[nodeN - 1];
			NCsize = new int[nodeN - 1];
		}
		int nnodeN = 0;
		for (int i = 0; i < nodeN; i++){
			if (i != delN){
				NCnode[nnodeN].fn = cNode[i].fn;
				for (int j = 0; j < 3; j++){
					NCnode[nnodeN].pos[j] = cNode[i].pos[j];
					if (j != 2) NCnode[nnodeN].uv[j] = cNode[i].uv[j];
				}
				NPcon[nnodeN].cType = pCons[i].cType;
				if (pCons[i].cType != -1){
					for (int j = 0; j < 3; j++){
						NPcon[nnodeN].pos[j] = pCons[i].pos[j];
					}
				}
				if (CONNECT_LIST_EXIST){
					NCsize[nnodeN] = Csize[i];
					if (Csize[i] > 0){
						NClist[nnodeN] = new int[Csize[i]];
						for (int j = 0; j < Csize[i]; j++){
							NClist[nnodeN][j] = Clist[i][j];
							if (i != nnodeN){
								if (connect[Clist[i][j]].c1 == selectC && connect[Clist[i][j]].n1 == i) connect[Clist[i][j]].n1 = nnodeN;
								if (connect[Clist[i][j]].c2 == selectC && connect[Clist[i][j]].n2 == i) connect[Clist[i][j]].n2 = nnodeN;
							}
						}
					}
				}
				nnodeN++;
			}
		}

		int nedgeN = 0;
		for (int i = 0; i < nodeN - 1; i++){
			if (i != delN){
				if (i != delN - 1){
					NDcon[nedgeN].cType = dCons[i].cType;
					if (dCons[i].cType != -1){
						for (int j = 0; j < 3; j++){
							NDcon[nedgeN].pos[j] = dCons[i].pos[j];
						}
					}
				}
				else {
					if (dCons[i].cType == -1 && dCons[i + 1].cType == -1){
						NDcon[nedgeN].cType = -1;
					}
					else if (dCons[i].cType == -1 && dCons[i + 1].cType != -1){
						NDcon[nedgeN].cType = dCons[i + 1].cType;
						for (int j = 0; j < 3; j++){
							NDcon[nedgeN].pos[j] = dCons[i + 1].pos[j];
						}
					}
					else if (dCons[i].cType != -1 && dCons[i + 1].cType == -1){
						NDcon[nedgeN].cType = dCons[i].cType;
						for (int j = 0; j < 3; j++){
							NDcon[nedgeN].pos[j] = dCons[i].pos[j];
						}
					}
					else {
						NDcon[nedgeN].cType = 1;
						for (int j = 0; j < 3; j++){
							NDcon[nedgeN].pos[j] = (dCons[i].pos[j] + dCons[i + 1].pos[j]) / 2.0;
						}
					}
				}
				nedgeN++;
			}
		}

		delete[] cNode;
		delete[] pCons;
		delete[] dCons;
		cNode = NCnode;
		pCons = NPcon;
		dCons = NDcon;

		if (CONNECT_LIST_EXIST){
			for (int i = 0; i < nodeN; i++){
				if (Csize[i] > 0) delete[] Clist[i];
			}
			delete[] Csize;
			delete[] Clist;
			Clist = NClist;
			Csize = NCsize;
		}

		nodeN--;
		delete[] Mframe;
		delete[] TPframe;
		Mframe = new adaptive_frame[nodeN - 1];
		TPframe = new adaptive_frame[nodeN - 1];
		update_frame_for_minimization(true);
		if(!minimize_discrete_rod_energy(20)) return;
	}

	//長いエッジを分割する
	while (1){
		float maxL;
		int maxN = -1;
		for (int i = 0; i < nodeN - 1; i++){
			float leng = 0.0;
			for (int j = 0; j < 3; j++){
				leng += pow(cNode[i + 1].pos[j] - cNode[i].pos[j], 2);
			}
			leng = sqrt(leng);
			if ((pCons[i].cType == -1 || pCons[i + 1].cType == -1)&&dCons[i].cType == -1){
				if (maxN == -1 || leng > maxL){
					maxN = i;
					maxL = leng;
				}
			}
		}

		if (maxN == -1 || nodeN == 200) break;
		if (maxL < 5.0*sParam->Wid) break;
		
		if (selectN > maxN) selectN++;
		if (selectE > maxN) selectE++;

		struct curve_node *NCnode;
		NCnode = new struct curve_node[nodeN + 1];
		struct ConstraintData *NPcon, *NDcon;
		NPcon = new struct ConstraintData[nodeN + 1];
		NDcon = new struct ConstraintData[nodeN];

		int **NClist, *NCsize;
		if (CONNECT_LIST_EXIST){
			NClist = new int*[nodeN + 1];
			NCsize = new int[nodeN + 1];
		}

		int NnodeN = 0;
		bool point_is_located = true;
		for (int i = 0; i < nodeN; i++){
			NCnode[NnodeN].fn = cNode[i].fn;
			for (int j = 0; j < 3; j++){
				NCnode[NnodeN].pos[j] = cNode[i].pos[j];
				if (j != 2) NCnode[NnodeN].uv[j] = cNode[i].uv[j];
			}
			NPcon[NnodeN].cType = pCons[i].cType;
			if (pCons[i].cType != -1){
				for (int j = 0; j < 3; j++){
					NPcon[NnodeN].pos[j] = pCons[i].pos[j];
				}
			}
			if (CONNECT_LIST_EXIST){
				NCsize[NnodeN] = Csize[i];
				if (Csize[i] > 0){
					NClist[NnodeN] = new int[Csize[i]];
					for (int j = 0; j < Csize[i]; j++){
						NClist[NnodeN][j] = Clist[i][j];
						if (i != NnodeN){
							if (connect[Clist[i][j]].c1 == selectC && connect[Clist[i][j]].n1 == i) connect[Clist[i][j]].n1 = NnodeN;
							if (connect[Clist[i][j]].c2 == selectC && connect[Clist[i][j]].n2 == i) connect[Clist[i][j]].n2 = NnodeN;
						}
					}
				}
			}
			NnodeN++;
			if (i == maxN){
				float pos[3];
				for (int j = 0; j < 3; j++){
					pos[j] = (cNode[i].pos[j] + cNode[i + 1].pos[j]) / 2.0;
				}
				int fn;
				float cord[3];
				point_is_located = AABB->closest_point_search_with_mesh_info(pos, &fn, cord);
				if (cord[0] < 0.0 || cord[1] < 0.0 || cord[2] < 0.0 || cord[0] + cord[1] + cord[2] > 1.0){
					point_is_located = false;
				}
				NCnode[NnodeN].fn = fn;
				NCnode[NnodeN].uv[0] = cord[1];
				NCnode[NnodeN].uv[1] = cord[2];
				lmesh->evaluate_box_spline_on_triangle(fn, cord, NCnode[NnodeN].pos);

				float leng = 0.0;
				for (int j = 0; j < 3; j++){
					leng += pow(NCnode[NnodeN].pos[j] - pos[j], 2);
				}
				leng = sqrt(leng);
				if (leng > 5*sParam->Wid) {
					point_is_located = false;
				}
				NPcon[NnodeN].cType = -1;

				if(CONNECT_LIST_EXIST) NCsize[NnodeN] = 0;
				NnodeN++;
			}
		}

		if (!point_is_located){
			delete[] NCnode;
			delete[] NPcon;
			delete[] NDcon;
			break;
		}

		int nedgeN = 0;
		for (int i = 0; i < nodeN - 1; i++){
			NDcon[nedgeN].cType = dCons[i].cType;
			if (dCons[nedgeN].cType != -1){
				for (int j = 0; j < 3; j++){
					NDcon[nedgeN].pos[j] = dCons[i].pos[j];
				}
			}
			nedgeN++;
			if (i == maxN){
				NDcon[nedgeN].cType = dCons[i].cType;
				if (dCons[i].cType != -1){
					for (int j = 0; j < 3; j++){
						NDcon[nedgeN].pos[j] = dCons[i].pos[j];
					}
				}
				nedgeN++;
			}
		}


		delete[] cNode;
		delete[] pCons;
		delete[] dCons;
		cNode = NCnode;
		pCons = NPcon;
		dCons = NDcon;

		if (CONNECT_LIST_EXIST){
			for (int i = 0; i < nodeN; i++){
				if(Csize[i] > 0) delete[] Clist[i];
			}
			delete[] Csize;
			delete[] Clist;
			Clist = NClist;
			Csize = NCsize;
		}

		nodeN++;
		delete[] Mframe;
		delete[] TPframe;
		Mframe = new adaptive_frame[nodeN - 1];
		TPframe = new adaptive_frame[nodeN - 1];
		update_frame_for_minimization(true);
		if (!minimize_discrete_rod_energy(20)) return;

	}
}


///////////////
///CurveNetwork用
///////////////

//コンストラクタ
CurveNetwork::CurveNetwork(){

}


//デストラクタ
CurveNetwork::~CurveNetwork(){

	printf("CurveNet dest start\n");

	for (int i = 0; i < curSet.size(); i++){
		delete curSet[i];
	}

	printf("CurveNet dest end\n");

}


//曲線の追加
void CurveNetwork::addCurve(class curve_on_mesh *Curv){

	class curve_on_mesh *nCur;
	nCur = new class curve_on_mesh;

	nCur->nodeN = Curv->nodeN;
	nCur->cNode = new curve_node[nCur->nodeN];
	for (int i = 0; i < nCur->nodeN; i++){
		nCur->cNode[i].fn = Curv->cNode[i].fn;
		float uvw[3];
		uvw[1] = Curv->cNode[i].uv[0];
		uvw[2] = Curv->cNode[i].uv[1];
		uvw[0] = 1.0 - uvw[1] - uvw[2];
		lmesh->evaluate_box_spline_on_triangle(nCur->cNode[i].fn, uvw, nCur->cNode[i].pos);
		nCur->cNode[i].uv[0] = uvw[1];
		nCur->cNode[i].uv[1] = uvw[2];
	}
	nCur->lmesh = lmesh;
	nCur->sParam = sParam;

	curSet.push_back(nCur);

}


//曲線の削除
void CurveNetwork::deleteCurve(int cN){

	if (curSet.empty() || cN >= curSet.size()) return;

	if (!Connect.empty()){
		int *cmap;
		cmap = new int[curSet.size()];
		for (int i = 0; i < curSet.size(); i++){
			if (i < cN) cmap[i] = i;
			else if (i == cN) cmap[i] = -1;
			else cmap[i] = i - 1;
		}

		//接続リストの更新	
		int *cnmap;
		cnmap = new int[Connect.size()];
		int cnum = 0;
		for (int i = 0; i < Connect.size(); i++){
			if (Connect[i].c1 == cN || Connect[i].c2 == cN) cnmap[i] = -1;
			else {
				cnmap[i] = cnum;
				cnum++;
			}
		}
		
		for (int i = 0; i < Connect.size(); i++){
			if (cnmap[i] != -1){
				if (cnmap[i] != i){
					int c1 = Connect[i].c1;
					int c2 = Connect[i].c2;
					for (int j = 0; j < curSet[c1]->Csize[Connect[i].n1]; j++){
						if (curSet[c1]->Clist[Connect[i].n1][j] == i) curSet[c1]->Clist[Connect[i].n1][j] = cnmap[i];
					}
					for (int j = 0; j < curSet[c2]->Csize[Connect[i].n2]; j++){
						if (curSet[c2]->Clist[Connect[i].n2][j] == i)curSet[c2]->Clist[Connect[i].n2][j] = cnmap[i];
					}

					Connect[i].c1 = cmap[c1];
					Connect[i].c2 = cmap[c2];
				}
			}
			else {
				int c1 = Connect[i].c1;
				int n1 = Connect[i].n1;
				int c2 = Connect[i].c2;
				int n2 = Connect[i].n2;
				std::vector<int> nclist1,nclist2;
				for (int j = 0; j < curSet[c1]->Csize[n1]; j++){
					if (curSet[c1]->Clist[n1][j] != i) nclist1.push_back(curSet[c1]->Clist[n1][j]);
				}
				if (curSet[c1]->Csize[n1] != 0) delete[] curSet[c1]->Clist[n1];
				curSet[c1]->Csize[n1] = nclist1.size();
				if (!nclist1.empty()){
					int *NClist;
					NClist = new int[nclist1.size()];
					for (int j = 0; j < nclist1.size(); j++){
						NClist[j] = nclist1[j];
					}
					curSet[c1]->Clist[n1] = NClist;
				}
				for (int j = 0; j < curSet[c2]->Csize[n2]; j++){
					if (curSet[c2]->Clist[n2][j] != i) nclist2.push_back(curSet[c2]->Clist[n2][j]);
				}
				if (curSet[c2]->Csize[n2] != 0) delete[] curSet[c2]->Clist[n2];
				curSet[c2]->Csize[n2] = nclist2.size();
				if (!nclist2.empty()){
					int *NClist;
					NClist = new int[nclist2.size()];
					for (int j = 0; j < nclist2.size(); j++){
						NClist[j] = nclist2[j];
					}
					curSet[c2]->Clist[n2] = NClist;
				}
			}
		}

		int itrN = Connect.size();
		auto itr = Connect.begin();
		for (int i = 0; i < itrN; i++){
			if (cnmap[i] == -1){
				Connect.erase(itr);
			}
			else {
				itr++;
			}
		}

	}


	//曲線オブジェクトの削除
	int itrN = 0;
	auto itr = curSet.begin();
	while (itr != curSet.end()){
		if (itrN == cN){
			curSet.erase(itr);
			break;
		}
		itrN++;
		itr++;
	}
}


//コネクションの追加
void CurveNetwork::addConnect(int c1, int n1, int c2, int n2){

	if (!curSet[c1]->CONNECT_LIST_EXIST){
		curSet[c1]->Csize = new int[curSet[c1]->nodeN];
		curSet[c1]->Clist = new int*[curSet[c1]->nodeN];
		for (int i = 0; i < curSet[c1]->nodeN; i++){
			curSet[c1]->Csize[i] = 0;
		}
		curSet[c1]->CONNECT_LIST_EXIST = true;
	}
	if (!curSet[c2]->CONNECT_LIST_EXIST){
		curSet[c2]->Csize = new int[curSet[c2]->nodeN];
		curSet[c2]->Clist = new int*[curSet[c2]->nodeN];
		for (int i = 0; i < curSet[c2]->nodeN; i++){
			curSet[c2]->Csize[i] = 0;
		}
		curSet[c2]->CONNECT_LIST_EXIST = true;
	}

	int connN = Connect.size();
	
	int *NClist1;
	NClist1 = new int[curSet[c1]->Csize[n1] + 1];
	for (int i = 0; i < curSet[c1]->Csize[n1]; i++){
		NClist1[i] = curSet[c1]->Clist[n1][i];
	}
	NClist1[curSet[c1]->Csize[n1]] = connN;
	if (curSet[c1]->Csize[n1] != 0) delete[] curSet[c1]->Clist[n1];
	curSet[c1]->Csize[n1]++;
	curSet[c1]->Clist[n1] = NClist1;
	
	int *NClist2;
	NClist2 = new int[curSet[c2]->Csize[n2] + 1];
	for (int i = 0; i < curSet[c2]->Csize[n2]; i++){
		NClist2[i] = curSet[c2]->Clist[n2][i];
	}
	NClist2[curSet[c2]->Csize[n2]] = connN;
	if (curSet[c2]->Csize[n2] != 0) delete[] curSet[c2]->Clist[n2];
	curSet[c2]->Csize[n2]++;
	curSet[c2]->Clist[n2] = NClist2;

	struct ConnectData cdata;
	cdata.c1 = c1;
	cdata.c2 = c2;
	cdata.n1 = n1;
	cdata.n2 = n2;

	Connect.push_back(cdata);

	curSet[c1]->pCons[n1].cType = 2;
	curSet[c2]->pCons[n2].cType = 2;
	for (int i = 0; i < 3; i++){
		curSet[c1]->pCons[n1].pos[i] = curSet[c1]->cNode[n1].pos[i];
		curSet[c2]->pCons[n2].pos[i] = curSet[c2]->cNode[n2].pos[i];
	}

}


//コネクションの削除
void CurveNetwork::deleteConnect(int cN){

	for (int i = cN; i < Connect.size(); i++){
		class curve_on_mesh *cur1 = curSet[Connect[i].c1];
		class curve_on_mesh *cur2 = curSet[Connect[i].c2];
		int n1 = Connect[i].n1;
		int n2 = Connect[i].n2;
		if (i > cN){
			for (int j = 0; j < cur1->Csize[n1]; j++){
				if (cur1->Clist[n1][j] == i) cur1->Clist[n1][j] = i - 1;
			}
			for (int j = 0; j < cur2->Csize[n2]; j++){
				if (cur2->Clist[n2][j] == i) cur2->Clist[n2][j] = i - 1;
			}
		}
		else {
			std::vector<int> NClist1,NClist2;
			for (int j = 0; j < cur1->Csize[n1]; j++){
				if (cur1->Clist[n1][j] != i) NClist1.push_back(cur1->Clist[n1][j]);
			}
			if (cur1->Csize[n1] > 0) delete[] cur1->Clist[n1];
			cur1->Csize[n1] = NClist1.size();
			if (!NClist1.empty()){
				int *NCset;
				NCset = new int[NClist1.size()];
				for (int j = 0; j < NClist1.size(); j++){
					NCset[j] = NClist1[j];
				}
				cur1->Clist[n1] = NCset;
			}
			for (int j = 0; j < cur2->Csize[n2]; j++){
				if (cur2->Clist[n2][j] != i) NClist2.push_back(cur2->Clist[n2][j]);
			}
			if (cur2->Csize[n2] > 0)delete[] cur2->Clist[n2];
			cur2->Csize[n2] = NClist2.size();
			if (!NClist2.empty()){
				int *NCset;
				NCset = new int[NClist2.size()];
				for (int j = 0; j < NClist2.size(); j++){
					NCset[j] = NClist2[j];
				}
				cur2->Clist[n2] = NCset;
			}
		}
	}



	auto itr = Connect.begin();
	int itrN = 0;
	while (itr != Connect.end()){
		if (itrN == cN) {
			Connect.erase(itr);
			break;
		}
		itrN++;
		itr++;
	}

}


//初期コネクションの設定
void CurveNetwork::setInitConnection(){

	for (int i = 0; i < curSet.size(); i++){
		class curve_on_mesh *scur = curSet[i];
		for (int j = i; j < curSet.size(); j++){
			class curve_on_mesh *ocur = curSet[j];

			std::vector<std::pair<int, int>> Npair;
			for (int k = 0; k < scur->nodeN; k++){
				std::vector<int> myNei;
				std::vector<float> myLeng;
				int startN = 0;
				if (j == i) startN = k + 3;
				for (int l = startN; l < ocur->nodeN; l++){
					float leng = 0.0;
					for (int m = 0; m < 3; m++){
						leng += pow(scur->cNode[k].pos[m] - ocur->cNode[l].pos[m], 2);
					}
					leng = sqrt(leng);
					if (leng < sParam->Wid*3.0){
						myNei.push_back(l);
						myLeng.push_back(leng);
					}
				}
				if (!myNei.empty()){
					std::vector <std::vector<int>> Ngroup;
					std::vector<int> nhozon;
					for (int l = 0; l < myNei.size(); l++){
						if (nhozon.empty() || myNei[l] < nhozon[nhozon.size() - 1] + 3){
							nhozon.push_back(l);
						}
						else {
							std::vector<int> ncopy;
							for (int m = 0; m < nhozon.size(); m++){
								ncopy.push_back(nhozon[m]);
							}
							Ngroup.push_back(ncopy);
							nhozon.clear();
							nhozon.push_back(l);
						}
					}
					Ngroup.push_back(nhozon);
					for (int l = 0; l < Ngroup.size(); l++){
						if (Ngroup[l].size() == 1){
							Npair.push_back(make_pair(k, myNei[Ngroup[l][0]]));
						}
						else if(!Ngroup[l].empty()){
							int Enode = -1;
							int minN = -1;
							float minL;
							for (int m = 0; m < Ngroup[l].size(); m++){
								if (myNei[Ngroup[l][m]] == 0 || myNei[Ngroup[l][m]] == ocur->nodeN - 1){
									Enode = m;
								}
								if (minN == -1 || minL > myLeng[Ngroup[l][m]]){
									minN = m;
									minL = myLeng[Ngroup[l][m]];
								}
							}
							if (Enode != -1) Npair.push_back(make_pair(k, myNei[Ngroup[l][Enode]]));
							else if (minN != -1) Npair.push_back(make_pair(k, myNei[Ngroup[l][minN]]));
						}
					}
				}
			}

			if (Npair.empty()) continue;

			std::vector<int> glist1,glist2;
			for (int k = 0; k < Npair.size() + 1 ; k++){
				if (k != Npair.size() && (glist1.empty() || Npair[k].first < glist1[glist1.size() - 1] + 2)){
					glist1.push_back(Npair[k].first);
					int psec = Npair[k].second;
					bool pexist = false;
					auto it = glist2.begin();
					for (int l = 0; l < glist2.size(); l++){
						if (psec > glist2[l]) ++it;
						if (glist2[l] == psec) pexist = true;
					}
					if (!pexist) glist2.insert(it,psec);
				}
				else {
					std::vector<std::vector<int>> igroup;
					std::vector<int> ihozon;
					for (int l = 0; l < glist2.size(); l++){
						int gele = glist2[l];
						if (ihozon.empty() || gele < ihozon[ihozon.size() - 1] + 2){
							ihozon.push_back(gele);
						}
						else {
							std::vector<int> ncopy;
							for (int m = 0; m < ihozon.size(); m++){
								ncopy.push_back(ihozon[m]);
							}
							igroup.push_back(ncopy);
							ihozon.clear();
							ihozon.push_back(gele);
						}
					}
					igroup.push_back(ihozon);

					int repN = -1;
					/*for (int l = 0; l < glist1.size(); l++){
						if (glist1[l] == 0 || glist1[l] == scur->nodeN - 1) repN = glist1[l];
					}*/
					for (int l = 0; l < igroup.size(); l++){
						int repN2 = -1;
						/*for (int m = 0; m < igroup[l].size(); m++){
							if (igroup[l][m] == 0 || igroup[l][m] == ocur->nodeN - 1) repN2 = igroup[l][m];
						}*/
						struct ConnectData cdata;
						cdata.c1 = i;
						cdata.c2 = j;
						if (repN != -1 && repN2 != -1){
							cdata.n1 = repN;
							cdata.n2 = repN2;
						}
						else if (repN != -1){
							int minN = -1;
							float minL;
							for (int m = 0; m < igroup[l].size(); m++){
								float leng = 0.0;
								for (int n = 0; n < 3; n++){
									leng += pow(scur->cNode[repN].pos[n] - ocur->cNode[igroup[l][m]].pos[n], 2);
								}
								leng = sqrt(leng);
								if (minN == -1 || minL > leng){
									minN = igroup[l][m];
									minL = leng;
								}
							}
							cdata.n1 = repN;
							cdata.n2 = minN;
						}
						else if (repN2 != -1){
							int minN = -1;
							float minL;
							for (int m = 0; m < glist1.size(); m++){
								float leng = -0.0;
								for (int n = 0; n < 3; n++){
									leng += pow(scur->cNode[glist1[m]].pos[n] - ocur->cNode[repN2].pos[n], 2);
								}
								leng = sqrt(leng);
								if (minN == -1 || minL > leng){
									minN = glist1[m];
									minL = leng;
								}
							}
							cdata.n1 = minN;
							cdata.n2 = repN2;
						}
						else {
							int minN1 = -1;
							int minN2 = -1;
							float minL;
							for (int m = 0; m < glist1.size(); m++){
								for (int n = 0; n < igroup[l].size(); n++){
									float leng = 0.0;
									for (int o = 0; o < 3; o++){
										leng += pow(scur->cNode[glist1[m]].pos[o] - ocur->cNode[igroup[l][n]].pos[o], 2);
									}
									leng = sqrt(leng);
									if (minN1 == -1 || minL > leng){
										minN1 = glist1[m];
										minN2 = igroup[l][n];
										minL = leng;
									}
								}
							}
							cdata.n1 = minN1;
							cdata.n2 = minN2;
						}
						Connect.push_back(cdata);
					}

					glist1.clear();
					glist2.clear();
					if (k != Npair.size()) {
						glist1.push_back(Npair[k].first);
						glist2.push_back(Npair[k].second);
					}
				}
			}

		}
	}

	//各曲線のコネクションリストの作成
	for (int i = 0; i < curSet.size(); i++){
		if (curSet[i]->CONNECT_LIST_EXIST){
			for (int j = 0; j < curSet[i]->nodeN; j++){
				if (curSet[i]->Csize[j] > 0) delete[] curSet[i]->Clist[j];
			}
			delete[] curSet[i]->Clist;
			delete[] curSet[i]->Csize;
		}
		curSet[i]->CONNECT_LIST_EXIST = true;
		curSet[i]->Clist = new int*[curSet[i]->nodeN];
		curSet[i]->Csize = new int[curSet[i]->nodeN];
		for (int j = 0; j < curSet[i]->nodeN; j++){
			curSet[i]->Csize[j] = 0;
		}
	}
	for (int i = 0; i < Connect.size(); i++){
		class curve_on_mesh *cur1 = curSet[Connect[i].c1];
		class curve_on_mesh *cur2 = curSet[Connect[i].c2];
		cur1->Csize[Connect[i].n1]++;
		cur2->Csize[Connect[i].n2]++;
	}
	for (int i = 0; i < curSet.size(); i++){
		for (int j = 0; j < curSet[i]->nodeN; j++){
			if (curSet[i]->Csize[j] > 0) curSet[i]->Clist[j] = new int[curSet[i]->Csize[j]];
			curSet[i]->Csize[j] = 0;
		}
	}
	for (int i = 0; i < Connect.size(); i++){
		class curve_on_mesh *cur1 = curSet[Connect[i].c1];
		class curve_on_mesh *cur2 = curSet[Connect[i].c2];
		cur1->Clist[Connect[i].n1][cur1->Csize[Connect[i].n1]] = i;
		cur2->Clist[Connect[i].n2][cur2->Csize[Connect[i].n2]] = i;
		cur1->Csize[Connect[i].n1]++;
		cur2->Csize[Connect[i].n2]++;
	}
	for (int i = 0; i < curSet.size(); i++){
		for (int j = 0; j < curSet[i]->nodeN; j++){
			if (curSet[i]->Csize[j] > 0){
				curSet[i]->pCons[j].cType = 2;
				for (int k = 0; k < 3; k++){
					curSet[i]->pCons[j].pos[k] = curSet[i]->cNode[j].pos[k];
				}
			}
		}
	}

}


//曲線ネットワークの出力
bool CurveNetwork::exportCurveNet(std::string fname, class AABB_struct *AABB){

	if (curSet.empty() || Connect.empty()) return false;

	//曲線のリメッシュ（必要なら実装）

	//コネクションのジオメトリを考慮した接続トポロジーの整理（必要なら実装）




	//グラフ構造に基づくコネクション関係の整理
	int *cmap;
	bool *visited;
	cmap = new int[Connect.size()];
	visited = new bool[Connect.size()];
	for (int i = 0; i < Connect.size(); i++){
		cmap[i] = -1;
		visited[i] = false;
	}

	int conNum = 0;
	for (int i = 0; i < Connect.size(); i++){
		if (!visited[i]){
			struct qelem *startq;
			startq = new struct qelem;
			startq->ln = i;
			struct qelem *endq = traverseConnectGraph(i, visited, startq);

			struct qelem *cque = endq;
			while (1){
				cmap[cque->ln] = conNum;
				if(cque == startq) break;
				cque = cque->back;
			}
			conNum++;
		}
		if (cmap[i] == -1) printf("%d connection is not traversed\n");
	}

	std::vector<std::vector<std::pair<int,int>>> Cgroup;
	for (int i = 0; i < conNum; i++){
		std::vector<std::pair<int,int>> Cpair;
		for (int j = 0; j < Connect.size(); j++){
			if (cmap[j] == i){
				bool exist1 = false;
				bool exist2 = false;
				for (int k = 0; k < Cpair.size(); k++){
					if (Cpair[k].first == Connect[j].c1 && Cpair[k].second == Connect[j].n1) exist1 = true;
					if (Cpair[k].first == Connect[j].c2 && Cpair[k].second == Connect[j].n2) exist2 = true;
				}
				if (!exist1) Cpair.push_back(make_pair(Connect[j].c1, Connect[j].n1));
				if (!exist2) Cpair.push_back(make_pair(Connect[j].c2, Connect[j].n2));
			}
		}
		Cgroup.push_back(Cpair);
	}

	float(*gps)[3],(*nds)[3];
	gps = new float[conNum][3];
	nds = new float[conNum][3];
	for (int i = 0; i < Cgroup.size(); i++){
		for (int j = 0; j < 3; j++){
			gps[i][j] = 0.0;
		}
		for (int j = 0; j < Cgroup[i].size(); j++){
			int cN = Cgroup[i][j].first;
			int nN = Cgroup[i][j].second;
			for (int k = 0; k < 3; k++){
				gps[i][k] += curSet[cN]->cNode[nN].pos[k] / (float)Cgroup[i].size();
			}
		}
		int fn;
		float cord[3],dbs[2][3];
		AABB->closest_point_search_with_mesh_info(gps[i], &fn, cord);
		lmesh->evaluate_derivative_of_box_spline_on_triangle(fn, cord, dbs);

		float norm[3];
		float nh = 0.0;
		for (int j = 0; j < 3; j++){
			norm[j] = dbs[0][(j + 1) % 3] * dbs[1][(j + 2) % 3] - dbs[0][(j + 2) % 3] * dbs[1][(j + 1) % 3];
			nh += norm[j] * norm[j];
		}
		nh = sqrt(nh);
		if (nh < 1.0e-10) nh = 1.0e-10;
		for (int j = 0; j < 3; j++){
			nds[i][j] = norm[j] / nh;
		}
	}


	delete[] cmap;
	delete[] visited;


	//ファイル出力
	FILE *fp = fopen(fname.c_str(),"w");
	if (fp == NULL){
		return false;
	}

	fprintf(fp,"#Mesh\n");
	fprintf(fp,"%d %d\n",lmesh->nodeN,lmesh->faceN);
	for (int i = 0; i < lmesh->nodeN; i++){
		fprintf(fp, "%d %f %f %f\n",i+1,lmesh->Nodes[i][0],lmesh->Nodes[i][1],lmesh->Nodes[i][2]);
	}
	for (int i = 0; i < lmesh->faceN; i++){
		fprintf(fp,"%d %d %d %d\n",i+1,lmesh->Faces[i][0],lmesh->Faces[i][1],lmesh->Faces[i][2]);
	}

	fprintf(fp,"#Curves\n");
	fprintf(fp,"%d\n",curSet.size());
	for (int i = 0; i < curSet.size(); i++){
		class curve_on_mesh *cur = curSet[i];
		fprintf(fp,"%d %d\n",i,cur->nodeN);
		for (int j = 0; j < cur->nodeN; j++){
			fprintf(fp, "%d %f %f %f %d %f %f\n", j, cur->cNode[j].pos[0], cur->cNode[j].pos[1], cur->cNode[j].pos[2],
				cur->cNode[j].fn,cur->cNode[j].uv[0],cur->cNode[j].uv[1]);
		}
		for (int j = 0; j < cur->nodeN - 1; j++){
			fprintf(fp,"%d %f %f %f",j,cur->TPframe[j].d1[0],cur->TPframe[j].d1[1],cur->TPframe[j].d1[2]);
			fprintf(fp, " %f %f %f", cur->TPframe[j].d2[0], cur->TPframe[j].d2[1], cur->TPframe[j].d2[2]);
			fprintf(fp, " %f %f %f\n", cur->TPframe[j].d3[0], cur->TPframe[j].d3[1], cur->TPframe[j].d3[2]);
			fprintf(fp, " %f %f %f", cur->Mframe[j].d1[0], cur->Mframe[j].d1[1], cur->Mframe[j].d1[2]);
			fprintf(fp, " %f %f %f", cur->Mframe[j].d2[0], cur->Mframe[j].d2[1], cur->Mframe[j].d2[2]);
			fprintf(fp, " %f %f %f\n", cur->Mframe[j].d3[0], cur->Mframe[j].d3[1], cur->Mframe[j].d3[2]);
		}
	}

	fprintf(fp,"#Connection\n");
	fprintf(fp,"%d\n",Cgroup.size());
	for (int i = 0; i < Cgroup.size(); i++){
		fprintf(fp,"%d ",Cgroup[i].size());
		fprintf(fp, "%f %f %f %f %f %f\n", gps[i][0], gps[i][1], gps[i][2], nds[i][0], nds[i][1], nds[i][2]);
		for (int j = 0; j < Cgroup[i].size(); j++){
			fprintf(fp, "%d %d ", Cgroup[i][j].first, Cgroup[i][j].second);
		}
		fprintf(fp,"\n");
	}

	fclose(fp);

	return true;
}


//曲線ネットワークの読み込み
bool CurveNetwork::loadCurveNet(std::string fname,class loop_subdivision_mesh *imesh){

	ifstream file;
	file.open(fname, ios::in);

	if (!file.is_open()){
		std::cout << "File open error!!" << endl;
	}

	std::string line;
	getline(file, line);

	file >> imesh->nodeN >> imesh->faceN;
	if (imesh->nodeN > 0) imesh->Nodes = new float[imesh->nodeN][3];
	if (imesh->faceN > 0) imesh->Faces = new int[imesh->faceN][3];
	for (int i = 0; i < imesh->nodeN; i++){
		int nN;
		file >> nN >> imesh->Nodes[i][0] >> imesh->Nodes[i][1] >> imesh->Nodes[i][2];
	}
	for (int i = 0; i < imesh->faceN; i++){
		int fN;
		file >> fN >> imesh->Faces[i][0] >> imesh->Faces[i][1] >> imesh->Faces[i][2];
	}

	imesh->construct_data_structure_without_edge();
	imesh->prepare_for_evaluation();
	lmesh = imesh;


	getline(file, line);
	getline(file, line);

	int curN;
	file >> curN;

	for (int i = 0; i < curN; i++){
		class curve_on_mesh *cmesh;
		cmesh = new class curve_on_mesh;
		int cN;
		file >> cN >> cmesh->nodeN;
		if (cmesh->nodeN>0) cmesh->cNode = new struct curve_node[cmesh->nodeN];
		for (int j = 0; j < cmesh->nodeN; j++){
			int nN, fN;
			float uvw[3];
			float xv, yv, zv;
			file >> nN >> xv >> yv >> zv >> fN >> uvw[1] >> uvw[2];
			uvw[0] = 1.0 - uvw[1] - uvw[2];
			imesh->evaluate_box_spline_on_triangle(fN, uvw, cmesh->cNode[j].pos);
			cmesh->cNode[j].fn = fN;
			cmesh->cNode[j].uv[0] = uvw[1];
			cmesh->cNode[j].uv[1] = uvw[2];
		}
		for (int j = 0; j < cmesh->nodeN - 1; j++){
			int fN;
			float fxv, fyv, fzv;
			file >> fN >> fxv >> fyv >> fzv >> fxv >> fyv >> fzv >> fxv >> fyv >> fzv;
			file >> fxv >> fyv >> fzv >> fxv >> fyv >> fzv >> fxv >> fyv >> fzv;
		}
		cmesh->lmesh = imesh;
		cmesh->sParam = sParam;
		cmesh->prepare_for_minimization();

		cmesh->pCons[0].cType = 1;
		cmesh->pCons[cmesh->nodeN - 1].cType = 1;
		for (int j = 0; j < 3; j++){
			cmesh->pCons[0].pos[j] = cmesh->cNode[0].pos[j];
			cmesh->pCons[cmesh->nodeN - 1].pos[j] = cmesh->cNode[cmesh->nodeN - 1].pos[j];
		}
		curSet.push_back(cmesh);
	}

	getline(file, line);
	getline(file, line);

	int cN;
	file >> cN;
	for (int i = 0; i < cN; i++){
		int ccN;
		float xv, yv, zv;
		file >> ccN >> xv >> yv >> zv >> xv >> yv >> zv;
		std::vector<std::pair<int, int>> cpairs;
		for (int j = 0; j < ccN; j++){
			int cn, nn;
			file >> cn >> nn;
			cpairs.push_back(make_pair(cn, nn));
		}
		for (int j = 0; j < ccN - 1; j++){
			addConnect(cpairs[0].first, cpairs[0].second, cpairs[1 + j].first, cpairs[1 + j].second);
		}
	}

	return true;
}


//コネクション関係の走査
struct qelem* CurveNetwork::traverseConnectGraph(int cN, bool *visited, struct qelem *que){

	visited[cN] = true;

	class curve_on_mesh *cur1 = curSet[Connect[cN].c1];
	class curve_on_mesh *cur2 = curSet[Connect[cN].c2];
	int n1 = Connect[cN].n1;
	int n2 = Connect[cN].n2;

	struct qelem *cque = que;
	for (int i = 0; i < cur1->Csize[n1]; i++){
		if (!visited[cur1->Clist[n1][i]]){
			struct qelem *nque;
			nque = new struct qelem;
			nque->ln = cur1->Clist[n1][i];
			nque->back = cque;
			cque = traverseConnectGraph(cur1->Clist[n1][i], visited, nque);
		}
	}
	for (int i = 0; i < cur2->Csize[n2]; i++){
		if (!visited[cur2->Clist[n2][i]]){
			struct qelem *nque;
			nque = new struct qelem;
			nque->ln = cur2->Clist[n2][i];
			nque->back = cque;
			cque = traverseConnectGraph(cur2->Clist[n2][i], visited, nque);
		}
	}

	return cque;

}