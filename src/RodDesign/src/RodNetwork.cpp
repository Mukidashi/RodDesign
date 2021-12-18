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
#include "curve_network.h"
#include "Setting.h"

////////////
////discrete rodのメンバ関数
////////////

DiscreteRod::DiscreteRod(){

	NODE_EXIST = false;
	FRAME_EXIST = false;
	PREPARED = false;

	//hei = ROD_HEI;
	//wid = ROD_WID;
}


DiscreteRod::~DiscreteRod(){

	if (NODE_EXIST){
		delete[] Nodes;
	}

	if (FRAME_EXIST){
		delete[] Mframe;
		delete[] TPframe;
	}

	if (PREPARED){
		delete[] bk;
		delete[] tm;
		delete[] oLeng;
		delete[] nvlist;
		delete[] evlist;
	}
}


void DiscreteRod::update_frame(){

	for (int i = 0; i < nodeN - 1; i++){
		float d1[3], d2[3], d3[3];
		float leng = 0.0;
		for (int j = 0; j < 3; j++){
			leng += pow(Nodes[i + 1][j] - Nodes[i][j], 2);
		}
		leng = sqrt(leng);
		for (int j = 0; j < 3; j++){
			d3[j] = (Nodes[i + 1][j] - Nodes[i][j]) / leng;
		}
		float ax[3];
		float nseki = 0.0;
		float axh = 0.0;
		for (int j = 0; j < 3; j++){
			ax[j] = TPframe[i].d3[(j + 1) % 3] * d3[(j + 2) % 3] - TPframe[i].d3[(j + 2) % 3] * d3[(j + 1) % 3];
			nseki += TPframe[i].d3[j] * d3[j];
			axh += ax[j] * ax[j];
		}
		axh = sqrt(axh);

		float Rot[3][3];
		if (axh > sin(1.0e-5)){
			float quat[4];
			if (nseki > 1.0) nseki = 1.0;
			if (nseki < -1.0)nseki = -1.0;
			nseki = acos(nseki);
			quat[0] = cos(nseki / 2.0);
			for (int j = 0; j < 3; j++){
				quat[1 + j] = sin(nseki / 2.0)*ax[j] / axh;
			}
			derive_rot_mat_from_quaternion(quat, Rot);
		}
		else {
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 3; k++){
					if (j == k && nseki >= 0.0) Rot[j][k] = 1.0;
					else if (j == k) Rot[j][k] = -1.0;
					else Rot[j][k] = 0.0;
				}
			}
		}

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
			TPframe[i].d3[j] = d3[j];
			Mframe[i].d3[j] = d3[j];
		}

		float quat[4], Rot2[3][3];
		quat[0] = cos(twis(i) / 2.0);
		for (int j = 0; j < 3; j++){
			quat[j + 1] = sin(twis(i) / 2.0)*Mframe[i].d3[j];
		}
		derive_rot_mat_from_quaternion(quat, Rot2);
		for (int j = 0; j < 3; j++){
			d1[j] = 0.0;
			d2[j] = 0.0;
			for (int k = 0; k < 3; k++){
				d1[j] += Rot2[j][k] * TPframe[i].d1[k];
				d2[j] += Rot2[j][k] * TPframe[i].d2[k];
			}
		}
		for (int j = 0; j < 3; j++){
			Mframe[i].d1[j] = d1[j];
			Mframe[i].d2[j] = d2[j];
		}
	}


}


//ねじり・曲げの位置についての微分計算
void DiscreteRod::derivative_of_twist_bend(int in, double leng1, double leng2, double gram[2][3], double hesm[4][3][3],
	double grke[2][2][3], double grkt[2][2], double heske[2][4][3][3], double hket[2][2][2][3], double hktt[2][2]){

	double kb[3], tit[3];
	double kai = 1.0;
	for (int i = 0; i < 3; i++){
		kb[i] = Mframe[in - 1].d3[(i + 1) % 3] * Mframe[in].d3[(i + 2) % 3]
			- Mframe[in - 1].d3[(i + 2) % 3] * Mframe[in].d3[(i + 1) % 3];
		kai += Mframe[in - 1].d3[i] * Mframe[in].d3[i];
	}
	for (int i = 0; i < 3; i++){
		kb[i] = 2.0*kb[i] / kai;
		tit[i] = (Mframe[in - 1].d3[i] + Mframe[in].d3[i]) / kai;
	}

	for (int j = 0; j < 3; j++){
		gram[0][j] = kb[j] / leng1 / 2.0;
		gram[1][j] = kb[j] / leng2 / 2.0;
		for (int k = 0; k < 3; k++){
			hesm[0][j][k] = -kb[j] * (TPframe[in - 1].d3[k] + tit[k]) / 2.0 / leng1 / leng1;
			hesm[1][j][k] = -kb[j] * tit[k] / 2.0 / leng1 / leng2;
			hesm[2][j][k] = -kb[j] * tit[k] / 2.0 / leng1 / leng2;
			hesm[3][j][k] = -kb[j] * (TPframe[in].d3[k] + tit[k]) / 2.0 / leng2 / leng2;
			if (j != k){
				double sig = 1.0;
				if (k == (j + 1) % 3) sig = -1.0;
				hesm[0][j][k] += -sig*TPframe[in].d3[3 - j - k] / leng1 / leng1 / kai;
				hesm[1][j][k] += sig*TPframe[in - 1].d3[3 - j - k] / leng1 / leng2 / kai;
				hesm[2][j][k] += -sig*TPframe[in].d3[3 - j - k] / leng1 / leng2 / kai;
				hesm[3][j][k] += sig*TPframe[in - 1].d3[3 - j - k] / leng2 / leng2 / kai;
			}
			//論文記載のヘッシアン
			hesm[0][j][k] = -(kb[j] * (TPframe[in - 1].d3[k] + tit[k]) + (TPframe[in - 1].d3[j] + tit[j])*kb[k]) / 4.0 / leng1 / leng1;
			hesm[3][j][k] = -(kb[j] * (TPframe[in].d3[k] + tit[k]) + (TPframe[in].d3[j] + tit[j])*kb[k]) / 4.0 / leng2 / leng2;
		}
	}


	double ki[2], tid1[3], tid2[3];
	ki[0] = 0.0;
	ki[1] = 0.0;
	for (int i = 0; i < 3; i++){
		ki[0] += 0.5*(Mframe[in - 1].d2[i] + Mframe[in].d2[i])*kb[i];
		ki[1] -= 0.5*(Mframe[in - 1].d1[i] + Mframe[in].d1[i])*kb[i];
		tid1[i] = (Mframe[in - 1].d1[i] + Mframe[in].d1[i]) / kai;
		tid2[i] = (Mframe[in - 1].d2[i] + Mframe[in].d2[i]) / kai;
	}
	double oseki[2][2][3];
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 3; j++){
			oseki[i][0][j] = Mframe[in - 1 + i].d3[(j + 1) % 3] * tid1[(j + 2) % 3]
				- Mframe[in - 1 + i].d3[(j + 2) % 3] * tid1[(j + 1) % 3];
			oseki[i][1][j] = Mframe[in - 1 + i].d3[(j + 1) % 3] * tid2[(j + 2) % 3]
				- Mframe[in - 1 + i].d3[(j + 2) % 3] * tid2[(j + 1) % 3];
		}
	}


	grkt[0][0] = 0.0;
	grkt[0][1] = 0.0;
	grkt[1][0] = 0.0;
	grkt[1][1] = 0.0;
	for (int j = 0; j < 3; j++){
		grke[0][0][j] = (-ki[0] * tit[j] + oseki[1][1][j]) / leng1;
		grke[0][1][j] = (-ki[0] * tit[j] - oseki[0][1][j]) / leng2;
		grke[1][0][j] = (-ki[1] * tit[j] - oseki[1][0][j]) / leng1;
		grke[1][1][j] = (-ki[1] * tit[j] + oseki[0][0][j]) / leng2;
		grkt[0][0] += -0.5*kb[j] * Mframe[in - 1].d1[j];
		grkt[0][1] += -0.5*kb[j] * Mframe[in].d1[j];
		grkt[1][0] += -0.5*kb[j] * Mframe[in - 1].d2[j];
		grkt[1][1] += -0.5*kb[j] * Mframe[in].d2[j];
		for (int k = 0; k < 3; k++){
			heske[0][0][j][k] = (2.0 * ki[0] * tit[j] * tit[k] - oseki[1][1][j] * tit[k] - tit[j] * oseki[1][1][k]) / leng1 / leng1
				+ ki[0] * (Mframe[in - 1].d3[j] * Mframe[in - 1].d3[k]) / leng1 / leng1 / kai;
			heske[0][1][j][k] = -ki[0] * Mframe[in - 1].d3[j] * Mframe[in].d3[k] / kai / leng1 / leng2 + (2.0*ki[0] * tit[j] * tit[k]
				- oseki[1][1][j] * tit[k] + tit[j] * oseki[0][1][k]) / leng1 / leng2;
			heske[0][2][j][k] = -ki[0] * Mframe[in].d3[j] * Mframe[in - 1].d3[k] / kai / leng1 / leng2 + (2.0*ki[0] * tit[j] * tit[k]
				- tit[j] * oseki[1][1][k] + oseki[0][1][j] * tit[k]) / leng1 / leng2;
			heske[0][3][j][k] = (2.0*ki[0] * tit[j] * tit[k] + oseki[0][1][j] * tit[k] + tit[j] * oseki[0][1][k]) / leng2 / leng2
				+ ki[0] * (Mframe[in].d3[j] * Mframe[in].d3[k]) / leng2 / leng2 / kai;
			heske[1][0][j][k] = (2.0*ki[1] * tit[j] * tit[k] + oseki[1][0][j] * tit[k] + tit[j] * oseki[1][0][k]) / leng1 / leng1
				+ ki[1] * Mframe[in - 1].d3[j] * Mframe[in - 1].d3[k] / leng1 / leng1 / kai;
			heske[1][1][j][k] = -ki[1] * Mframe[in - 1].d3[j] * Mframe[in].d3[k] / leng1 / leng2 / kai
				+ (2.0*ki[1] * tit[j] * tit[k] + oseki[1][0][j] * tit[k] - tit[j] * oseki[0][0][k]) / leng1 / leng2;
			heske[1][2][j][k] = -ki[1] * Mframe[in].d3[j] * Mframe[in - 1].d3[k] / leng1 / leng2 / kai
				+ (2.0*ki[1] * tit[j] * tit[k] + tit[j] * oseki[1][0][k] - oseki[0][0][j] * tit[k]) / leng1 / leng2;
			heske[1][3][j][k] = (2.0*ki[1] * tit[j] * tit[k] - oseki[0][0][j] * tit[k] - tit[j] * oseki[0][0][k]) / leng2 / leng2
				+ ki[1] * Mframe[in].d3[j] * Mframe[in].d3[k] / leng2 / leng2 / kai;
			if (j == k) {
				heske[0][0][j][k] += -ki[0] / leng1 / leng1 / kai;
				heske[0][1][j][k] += -ki[0] / leng1 / leng2 / kai;
				heske[0][2][j][k] += -ki[0] / leng1 / leng2 / kai;
				heske[0][3][j][k] += -ki[0] / leng2 / leng2 / kai;
				heske[1][0][j][k] += -ki[1] / leng1 / leng1 / kai;
				heske[1][1][j][k] += -ki[1] / leng1 / leng2 / kai;
				heske[1][2][j][k] += -ki[1] / leng1 / leng2 / kai;
				heske[1][3][j][k] += -ki[1] / leng2 / leng2 / kai;
			}
			else {
				double sig = 1.0;
				if (k == (j + 1) % 3) sig = -1.0;
				heske[0][1][j][k] += -sig*tid2[3 - j - k] / leng1 / leng2;
				heske[0][2][j][k] += sig*tid2[3 - j - k] / leng1 / leng2;
				heske[1][1][j][k] += sig*tid1[3 - j - k] / leng1 / leng2;
				heske[1][2][j][k] += -sig*tid1[3 - j - k] / leng1 / leng2;
			}
			//論文記載のヘッシアン
			heske[0][0][j][k] += (kb[j] * Mframe[in - 1].d2[k] + Mframe[in - 1].d2[j] * kb[k]) / 4.0 / leng1 / leng1;
			heske[0][3][j][k] += (kb[j] * Mframe[in].d2[k] + Mframe[in].d2[j] * kb[k]) / 4.0 / leng2 / leng2;
			heske[1][0][j][k] -= (kb[j] * Mframe[in - 1].d1[k] + Mframe[in - 1].d1[j] * kb[k]) / 4.0 / leng1 / leng1;
			heske[1][3][j][k] -= (kb[j] * Mframe[in].d1[k] + Mframe[in].d1[j] * kb[k]) / 4.0 / leng2 / leng2;
		}
		for (int k = 0; k < 2; k++){
			hket[0][0][k][j] = -(Mframe[in].d3[(j + 1) % 3] * Mframe[in - 1 + k].d1[(j + 2) % 3] - Mframe[in].d3[(j + 2) % 3] * Mframe[in - 1 + k].d1[(j + 1) % 3]) / kai / leng1;
			hket[1][0][k][j] = -(Mframe[in].d3[(j + 1) % 3] * Mframe[in - 1 + k].d2[(j + 2) % 3] - Mframe[in].d3[(j + 2) % 3] * Mframe[in - 1 + k].d2[(j + 1) % 3]) / leng1 / kai;
			hket[0][1][k][j] = (Mframe[in - 1].d3[(j + 1) % 3] * Mframe[in - 1 + k].d1[(j + 2) % 3] - Mframe[in - 1].d3[(j + 2) % 3] * Mframe[in - 1 + k].d1[(j + 1) % 3]) / kai / leng2;
			hket[1][1][k][j] = (Mframe[in - 1].d3[(j + 1) % 3] * Mframe[in - 1 + k].d2[(j + 2) % 3] - Mframe[in - 1].d3[(j + 2) % 3] * Mframe[in - 1 + k].d2[(j + 1) % 3]) / kai / leng2;
			for (int l = 0; l < 3; l++){
				hket[0][0][k][j] += 0.5*kb[l] * Mframe[in - 1 + k].d1[l] * tit[j] / leng1;
				hket[1][0][k][j] += 0.5*kb[l] * Mframe[in - 1 + k].d2[l] * tit[j] / leng1;
				hket[0][1][k][j] += 0.5*kb[l] * Mframe[in - 1 + k].d1[l] * tit[j] / leng2;
				hket[1][1][k][j] += 0.5*kb[l] * Mframe[in - 1 + k].d2[l] * tit[j] / leng2;
			}
		}
	}
	for (int j = 0; j < 2; j++){
		hktt[0][j] = 0.0;
		hktt[1][j] = 0.0;
		for (int k = 0; k < 3; k++){
			hktt[0][j] += -0.5*kb[k] * Mframe[in - 1 + j].d2[k];
			hktt[1][j] += 0.5*kb[k] * Mframe[in - 1 + j].d1[k];
		}
	}

}


void derivative_of_twist_bend(int in, AdaptiveFrame *Mf, double leng1, double leng2, double gram[2][3], double hesm[4][3][3], double grke[2][2][3], double grkt[2][2],
	double heske[2][4][3][3], double hket[2][2][2][3], double hktt[2][2]){

	double kb[3], tit[3];
	double kai = 1.0;
	for (int i = 0; i < 3; i++){
		kb[i] = Mf[in - 1].d3[(i + 1) % 3] * Mf[in].d3[(i + 2) % 3]
			- Mf[in - 1].d3[(i + 2) % 3] * Mf[in].d3[(i + 1) % 3];
		kai += Mf[in - 1].d3[i] * Mf[in].d3[i];
	}
	for (int i = 0; i < 3; i++){
		kb[i] = 2.0*kb[i] / kai;
		tit[i] = (Mf[in - 1].d3[i] + Mf[in].d3[i]) / kai;
	}

	for (int j = 0; j < 3; j++){
		gram[0][j] = kb[j] / leng1 / 2.0;
		gram[1][j] = kb[j] / leng2 / 2.0;
		for (int k = 0; k < 3; k++){
			hesm[0][j][k] = -kb[j] * (Mf[in - 1].d3[k] + tit[k]) / 2.0 / leng1 / leng1;
			hesm[1][j][k] = -kb[j] * tit[k] / 2.0 / leng1 / leng2;
			hesm[2][j][k] = -kb[j] * tit[k] / 2.0 / leng1 / leng2;
			hesm[3][j][k] = -kb[j] * (Mf[in].d3[k] + tit[k]) / 2.0 / leng2 / leng2;
			if (j != k){
				double sig = 1.0;
				if (k == (j + 1) % 3) sig = -1.0;
				hesm[0][j][k] += -sig*Mf[in].d3[3 - j - k] / leng1 / leng1 / kai;
				hesm[1][j][k] += sig*Mf[in - 1].d3[3 - j - k] / leng1 / leng2 / kai;
				hesm[2][j][k] += -sig*Mf[in].d3[3 - j - k] / leng1 / leng2 / kai;
				hesm[3][j][k] += sig*Mf[in - 1].d3[3 - j - k] / leng2 / leng2 / kai;
			}
			//論文記載のヘッシアン
			hesm[0][j][k] = -(kb[j] * (Mf[in - 1].d3[k] + tit[k]) + (Mf[in - 1].d3[j] + tit[j])*kb[k]) / 4.0 / leng1 / leng1;
			hesm[3][j][k] = -(kb[j] * (Mf[in].d3[k] + tit[k]) + (Mf[in].d3[j] + tit[j])*kb[k]) / 4.0 / leng2 / leng2;
		}
	}


	double ki[2], tid1[3], tid2[3];
	ki[0] = 0.0;
	ki[1] = 0.0;
	for (int i = 0; i < 3; i++){
		ki[0] += 0.5*(Mf[in - 1].d2[i] + Mf[in].d2[i])*kb[i];
		ki[1] -= 0.5*(Mf[in - 1].d1[i] + Mf[in].d1[i])*kb[i];
		tid1[i] = (Mf[in - 1].d1[i] + Mf[in].d1[i]) / kai;
		tid2[i] = (Mf[in - 1].d2[i] + Mf[in].d2[i]) / kai;
	}
	double oseki[2][2][3];
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 3; j++){
			oseki[i][0][j] = Mf[in - 1 + i].d3[(j + 1) % 3] * tid1[(j + 2) % 3]
				- Mf[in - 1 + i].d3[(j + 2) % 3] * tid1[(j + 1) % 3];
			oseki[i][1][j] = Mf[in - 1 + i].d3[(j + 1) % 3] * tid2[(j + 2) % 3]
				- Mf[in - 1 + i].d3[(j + 2) % 3] * tid2[(j + 1) % 3];
		}
	}


	grkt[0][0] = 0.0;
	grkt[0][1] = 0.0;
	grkt[1][0] = 0.0;
	grkt[1][1] = 0.0;
	for (int j = 0; j < 3; j++){
		grke[0][0][j] = (-ki[0] * tit[j] + oseki[1][1][j]) / leng1;
		grke[0][1][j] = (-ki[0] * tit[j] - oseki[0][1][j]) / leng2;
		grke[1][0][j] = (-ki[1] * tit[j] - oseki[1][0][j]) / leng1;
		grke[1][1][j] = (-ki[1] * tit[j] + oseki[0][0][j]) / leng2;
		grkt[0][0] += -0.5*kb[j] * Mf[in - 1].d1[j];
		grkt[0][1] += -0.5*kb[j] * Mf[in].d1[j];
		grkt[1][0] += -0.5*kb[j] * Mf[in - 1].d2[j];
		grkt[1][1] += -0.5*kb[j] * Mf[in].d2[j];
		for (int k = 0; k < 3; k++){
			heske[0][0][j][k] = (2.0 * ki[0] * tit[j] * tit[k] - oseki[1][1][j] * tit[k] - tit[j] * oseki[1][1][k]) / leng1 / leng1
				+ ki[0] * (Mf[in - 1].d3[j] * Mf[in - 1].d3[k]) / leng1 / leng1 / kai;
			heske[0][1][j][k] = -ki[0] * Mf[in - 1].d3[j] * Mf[in].d3[k] / kai / leng1 / leng2 + (2.0*ki[0] * tit[j] * tit[k]
				- oseki[1][1][j] * tit[k] + tit[j] * oseki[0][1][k]) / leng1 / leng2;
			heske[0][2][j][k] = -ki[0] * Mf[in].d3[j] * Mf[in - 1].d3[k] / kai / leng1 / leng2 + (2.0*ki[0] * tit[j] * tit[k]
				- tit[j] * oseki[1][1][k] + oseki[0][1][j] * tit[k]) / leng1 / leng2;
			heske[0][3][j][k] = (2.0*ki[0] * tit[j] * tit[k] + oseki[0][1][j] * tit[k] + tit[j] * oseki[0][1][k]) / leng2 / leng2
				+ ki[0] * (Mf[in].d3[j] * Mf[in].d3[k]) / leng2 / leng2 / kai;
			heske[1][0][j][k] = (2.0*ki[1] * tit[j] * tit[k] + oseki[1][0][j] * tit[k] + tit[j] * oseki[1][0][k]) / leng1 / leng1
				+ ki[1] * Mf[in - 1].d3[j] * Mf[in - 1].d3[k] / leng1 / leng1 / kai;
			heske[1][1][j][k] = -ki[1] * Mf[in - 1].d3[j] * Mf[in].d3[k] / leng1 / leng2 / kai
				+ (2.0*ki[1] * tit[j] * tit[k] + oseki[1][0][j] * tit[k] - tit[j] * oseki[0][0][k]) / leng1 / leng2;
			heske[1][2][j][k] = -ki[1] * Mf[in].d3[j] * Mf[in - 1].d3[k] / leng1 / leng2 / kai
				+ (2.0*ki[1] * tit[j] * tit[k] + tit[j] * oseki[1][0][k] - oseki[0][0][j] * tit[k]) / leng1 / leng2;
			heske[1][3][j][k] = (2.0*ki[1] * tit[j] * tit[k] - oseki[0][0][j] * tit[k] - tit[j] * oseki[0][0][k]) / leng2 / leng2
				+ ki[1] * Mf[in].d3[j] * Mf[in].d3[k] / leng2 / leng2 / kai;
			if (j == k) {
				heske[0][0][j][k] += -ki[0] / leng1 / leng1 / kai;
				heske[0][1][j][k] += -ki[0] / leng1 / leng2 / kai;
				heske[0][2][j][k] += -ki[0] / leng1 / leng2 / kai;
				heske[0][3][j][k] += -ki[0] / leng2 / leng2 / kai;
				heske[1][0][j][k] += -ki[1] / leng1 / leng1 / kai;
				heske[1][1][j][k] += -ki[1] / leng1 / leng2 / kai;
				heske[1][2][j][k] += -ki[1] / leng1 / leng2 / kai;
				heske[1][3][j][k] += -ki[1] / leng2 / leng2 / kai;
			}
			else {
				double sig = 1.0;
				if (k == (j + 1) % 3) sig = -1.0;
				heske[0][1][j][k] += -sig*tid2[3 - j - k] / leng1 / leng2;
				heske[0][2][j][k] += sig*tid2[3 - j - k] / leng1 / leng2;
				heske[1][1][j][k] += sig*tid1[3 - j - k] / leng1 / leng2;
				heske[1][2][j][k] += -sig*tid1[3 - j - k] / leng1 / leng2;
			}
			//論文記載のヘッシアン
			heske[0][0][j][k] += (kb[j] * Mf[in - 1].d2[k] + Mf[in - 1].d2[j] * kb[k]) / 4.0 / leng1 / leng1;
			heske[0][3][j][k] += (kb[j] * Mf[in].d2[k] + Mf[in].d2[j] * kb[k]) / 4.0 / leng2 / leng2;
			heske[1][0][j][k] -= (kb[j] * Mf[in - 1].d1[k] + Mf[in - 1].d1[j] * kb[k]) / 4.0 / leng1 / leng1;
			heske[1][3][j][k] -= (kb[j] * Mf[in].d1[k] + Mf[in].d1[j] * kb[k]) / 4.0 / leng2 / leng2;
		}
		for (int k = 0; k < 2; k++){
			hket[0][0][k][j] = -(Mf[in].d3[(j + 1) % 3] * Mf[in - 1 + k].d1[(j + 2) % 3] - Mf[in].d3[(j + 2) % 3] * Mf[in - 1 + k].d1[(j + 1) % 3]) / kai / leng1;
			hket[1][0][k][j] = -(Mf[in].d3[(j + 1) % 3] * Mf[in - 1 + k].d2[(j + 2) % 3] - Mf[in].d3[(j + 2) % 3] * Mf[in - 1 + k].d2[(j + 1) % 3]) / leng1 / kai;
			hket[0][1][k][j] = (Mf[in - 1].d3[(j + 1) % 3] * Mf[in - 1 + k].d1[(j + 2) % 3] - Mf[in - 1].d3[(j + 2) % 3] * Mf[in - 1 + k].d1[(j + 1) % 3]) / kai / leng2;
			hket[1][1][k][j] = (Mf[in - 1].d3[(j + 1) % 3] * Mf[in - 1 + k].d2[(j + 2) % 3] - Mf[in - 1].d3[(j + 2) % 3] * Mf[in - 1 + k].d2[(j + 1) % 3]) / kai / leng2;
			for (int l = 0; l < 3; l++){
				hket[0][0][k][j] += 0.5*kb[l] * Mf[in - 1 + k].d1[l] * tit[j] / leng1;
				hket[1][0][k][j] += 0.5*kb[l] * Mf[in - 1 + k].d2[l] * tit[j] / leng1;
				hket[0][1][k][j] += 0.5*kb[l] * Mf[in - 1 + k].d1[l] * tit[j] / leng2;
				hket[1][1][k][j] += 0.5*kb[l] * Mf[in - 1 + k].d2[l] * tit[j] / leng2;
			}
		}
	}
	for (int j = 0; j < 2; j++){
		hktt[0][j] = 0.0;
		hktt[1][j] = 0.0;
		for (int k = 0; k < 3; k++){
			hktt[0][j] += -0.5*kb[k] * Mf[in - 1 + j].d2[k];
			hktt[1][j] += 0.5*kb[k] * Mf[in - 1 + j].d1[k];
		}
	}
}




////////////
////RodNetworkのメンバ関数
////////////


RodNetwork::RodNetwork(){

	rodN = 0;
	conN = 0;

}


RodNetwork::~RodNetwork(){


	printf("RodNet dest start\n");
	if (rodN > 0){
		delete[] Rods;
	}

	if (conN > 0){
		for (int i = 0; i < conN; i++){
			if (Connect[i].edgeN > 0) delete[] Connect[i].cEdge;
			if (Connect[i].nodeN > 0)delete[] Connect[i].cNode;
		}

		delete[] Connect;
	}
	printf("RodNet dest start\n");

}


void RodNetwork::prepare_for_simulation(){

	old_E = -1.0;

	for (int i = 0; i < rodN; i++){
		Rods[i].vXk = Eigen::VectorXd::Zero(3 * Rods[i].nodeN);
		Rods[i].vTk = Eigen::VectorXd::Zero(Rods[i].nodeN - 1);
		Rods[i].twis.resize(Rods[i].nodeN - 1);
		for (int j = 0; j < Rods[i].nodeN - 1; j++){
			double theta = 0.0;
			double nd[3];
			double seki = 0.0;
			for (int k = 0; k < 3; k++){
				theta += Rods[i].Mframe[j].d1[k] * Rods[i].TPframe[j].d1[k];
				nd[k] = Rods[i].TPframe[j].d1[(k + 1) % 3] * Rods[i].Mframe[j].d1[(k + 2) % 3]
					- Rods[i].TPframe[j].d1[(k + 2) % 3] * Rods[i].Mframe[j].d1[(k + 1) % 3];
				seki += nd[k] * Rods[i].TPframe[j].d3[k];
			}
			if (theta > 1.0) theta = 1.0;
			if (theta < -1.0) theta = -1.0;
			theta = acos(theta);
			if (seki < 0.0) theta *= -1.0;
			Rods[i].twis(j) = theta;
		}
		Rods[i].update_frame();
	}

	//各ロッドの設定
	for (int i = 0; i < rodN; i++){
		if (Rods[i].PREPARED){
			delete[] Rods[i].oLeng;
			delete[] Rods[i].bk;
			delete[] Rods[i].tm;
			delete[] Rods[i].oNodes;
			delete[] Rods[i].oTPframe;
			delete[] Rods[i].oMframe;
		}

		Rods[i].oLeng = new double[Rods[i].nodeN - 1];
		Rods[i].bk = new double[Rods[i].nodeN - 2][2];
		Rods[i].tm = new double[Rods[i].nodeN - 2];

		double fleng = 0.0;
		for (int j = 0; j < Rods[i].nodeN - 1; j++){
			double leng = 0.0;
			for (int k = 0; k < 3; k++){
				leng += pow(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k], 2);
			}
			leng = sqrt(leng);
			fleng += leng;
			Rods[i].oLeng[j] = leng;

			if (j == Rods[i].nodeN - 2) continue;
			Rods[i].bk[j][0] = 0.0;
			Rods[i].bk[j][1] = 0.0;
			Rods[i].tm[j] = 0.0;
		}
		Rods[i].oNodes = new double[Rods[i].nodeN][3];
		Rods[i].oTPframe = new AdaptiveFrame[Rods[i].nodeN - 1];
		Rods[i].oMframe = new AdaptiveFrame[Rods[i].nodeN - 1];
		for (int j = 0; j < Rods[i].nodeN; j++){
			for (int k = 0; k < 3; k++){
				Rods[i].oNodes[j][k] = Rods[i].Nodes[j][k];
				if (j != Rods[i].nodeN - 1){
					Rods[i].oTPframe[j].d1[k] = Rods[i].TPframe[j].d1[k];
					Rods[i].oTPframe[j].d2[k] = Rods[i].TPframe[j].d2[k];
					Rods[i].oTPframe[j].d3[k] = Rods[i].TPframe[j].d3[k];
					Rods[i].oMframe[j].d1[k] = Rods[i].Mframe[j].d1[k];
					Rods[i].oMframe[j].d2[k] = Rods[i].Mframe[j].d2[k];
					Rods[i].oMframe[j].d3[k] = Rods[i].Mframe[j].d3[k];
				}
			}
		}

		Rods[i].PREPARED = true;
	}

	//各ロッドの集中質量
	for (int i = 0; i < rodN; i++){
		Rods[i].Mass.resize(4 * Rods[i].nodeN - 1);
		Rods[i].Mass = Eigen::VectorXd::Zero(4 * Rods[i].nodeN - 1);
		for (int j = 0; j < Rods[i].nodeN - 1; j++){
			double leng = Rods[i].oLeng[j];
			double wei = sParam->Density*leng / 2.0*sParam->Wid*sParam->Hei;
			for (int k = 0; k < 3; k++){
				Rods[i].Mass(4 * j + k) += wei;
				Rods[i].Mass(4 * j + 4 + k) += wei;
			}
			Rods[i].Mass(4 * j + 3) += sParam->Density*(pow(sParam->Wid, 3)*sParam->Hei + pow(sParam->Hei, 3)*sParam->Wid)*leng / 24.0;
		}
	}


	//コネクションの設定
	for (int i = 0; i < conN; i++){
		//質量・慣性能率
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				Connect[i].Iner[j][k] = 0.0;
			}
		}
		double fleng = 0.0;
		for (int j = 0; j < Connect[i].edgeN; j++){
			int rN = Connect[i].cEdge[j].rN;
			int eN = Connect[i].cEdge[j].eN;
			double leng1 = 0.0;
			double leng2 = 0.0;
			for (int k = 0; k < 3; k++){
				leng1 += pow(Connect[i].gp[k] - Rods[rN].Nodes[eN][k], 2);
				leng2 += pow(Connect[i].gp[k] - Rods[rN].Nodes[eN + 1][k], 2);
			}
			fleng += sqrt(leng1) + sqrt(leng2);
		}
		fleng /= 2.0*(double)Connect[i].edgeN;
		Connect[i].Mass = rho2*fleng*fleng*M_PI*sParam->Hei;
		Connect[i].Iner[0][0] = (fleng*fleng / 4.0 + sParam->Hei*sParam->Hei / 12.0)*Connect[i].Mass;
		Connect[i].Iner[1][1] = (fleng*fleng / 4.0 + sParam->Hei*sParam->Hei / 12.0)*Connect[i].Mass;
		Connect[i].Iner[2][2] = (fleng*fleng + fleng*fleng) / 4.0*Connect[i].Mass;

		//フレームの設定
		double frame[3][3];
		int maxN;
		double maxV;
		double seki = 0.0;
		for (int j = 0; j < 3; j++){
			frame[2][j] = Connect[i].norm[j];
			frame[0][j] = Connect[i].norm[j];
			if (j == 0 || maxV < frame[2][j]){
				maxN = j;
				maxV = frame[2][j];
			}
		}
		frame[0][maxN] = -(pow(frame[0][(maxN + 1) % 3], 2) + pow(frame[0][(maxN + 2) % 3], 2)) / frame[2][maxN];

		double leng = 0.0;
		for (int j = 0; j < 3; j++){
			leng += pow(frame[0][j], 2);
		}
		leng = sqrt(leng);
		if (leng < 1.0e-10) leng = 1.0e-10;
		for (int j = 0; j < 3; j++){
			frame[0][j] /= leng;
		}
		leng = 0.0;
		for (int j = 0; j < 3; j++){
			frame[1][j] = frame[2][(j + 1) % 3] * frame[0][(j + 2) % 3] - frame[2][(j + 2) % 3] * frame[0][(j + 1) % 3];
			leng += frame[1][j] * frame[1][j];
		}
		leng = sqrt(leng);
		for (int j = 0; j < 3; j++){
			frame[1][j] /= leng;
		}

		//相対座標系
		for (int j = 0; j < Connect[i].edgeN; j++){
			int rN = Connect[i].cEdge[j].rN;
			int eN = Connect[i].cEdge[j].eN;
			for (int k = 0; k < 3; k++){
				Connect[i].cEdge[j].frame[0][k] = 0.0;
				Connect[i].cEdge[j].frame[1][k] = 0.0;
				Connect[i].cEdge[j].frame[2][k] = 0.0;
				for (int l = 0; l < 3; l++){
					Connect[i].cEdge[j].frame[0][k] += frame[k][l] * Rods[rN].Mframe[eN].d1[l];
					Connect[i].cEdge[j].frame[1][k] += frame[k][l] * Rods[rN].Mframe[eN].d2[l];
					Connect[i].cEdge[j].frame[2][k] += frame[k][l] * Rods[rN].Mframe[eN].d3[l];
				}
			}
		}
		for (int j = 0; j < Connect[i].nodeN; j++){
			int rN = Connect[i].cNode[j].rN;
			int nN = Connect[i].cNode[j].nN;
			for (int k = 0; k < 3; k++){
				Connect[i].cNode[j].pos[k] = 0.0;
				for (int l = 0; l < 3; l++){
					Connect[i].cNode[j].pos[k] += frame[k][l] * (Rods[rN].Nodes[nN][l] - Connect[i].gp[l]);
				}
			}
		}

		double Rot[3][3];
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				Rot[j][k] = frame[k][j];
			}
		}
		derive_quaternion_from_rot_mat(Rot, Connect[i].quat);

		for (int j = 0; j < 4; j++){
			Connect[i].oquat[j] = Connect[i].quat[j];
			if (j != 3) Connect[i].ogp[j] = Connect[i].gp[j];
		}
	}


	//変数マップの作成
	svarN = 0;
	for (int i = 0; i < rodN; i++){
		Rods[i].nvlist = new int[Rods[i].nodeN];
		Rods[i].evlist = new int[Rods[i].nodeN - 1];
		for (int j = 0; j < Rods[i].nodeN; j++){
			Rods[i].nvlist[j] = svarN;
			svarN += 3;
			if (j != Rods[i].nodeN - 1) {
				Rods[i].evlist[j] = svarN;
				svarN++;
			}
		}
	}

	snodeN = svarN;
	svarN += 6*conN;

	cvlist = new int[conN];
	slamdN = 0;
	for (int i = 0; i < conN; i++){
		cvlist[i] = svarN + slamdN;
		slamdN += 3 * Connect[i].nodeN + Connect[i].edgeN;
	}
	slamdN += 6;
	svarN += slamdN;


	//剛体速度初期化
	vTc = Eigen::VectorXd::Zero(3 * conN);
	Wc = Eigen::VectorXd::Zero(3 * conN);

	Lamd = Eigen::VectorXd::Zero(slamdN);

	preSol = Eigen::VectorXd::Zero(svarN - snodeN - 6 * conN);
}


//剛体変換拘束非対応
void RodNetwork::simulate_for_animation(float tstep,bool &tstep_is_large, bool &is_converged,bool &is_diverged){

	//剛性マトリックス・残差ベクトルの計算
	Eigen::SparseMatrix<double> Kmat(svarN, svarN);
	Eigen::VectorXd bvec(svarN), sol(svarN), resF(snodeN);
	std::vector<Triplet<double>> kmlist;
	bvec = Eigen::VectorXd::Zero(svarN);
	resF = Eigen::VectorXd::Zero(snodeN);

	for (int i = 0; i < rodN; i++){
		for (int j = 0; j < Rods[i].nodeN; j++){
			int vN = Rods[i].nvlist[j];
			for (int k = 0; k < 3; k++){
				bvec(vN + k) += tstep*Rods[i].vXk(3 * j + k);
				//if (k == 2) bvec(vN + k) -= tstep*tstep * 9800;
			}
			if (j != Rods[i].nodeN - 1){
				bvec(vN + 3) += tstep*Rods[i].vTk(j);
			}
		}
	}

	for (int i = 0; i < conN; i++){
		for (int j = 0; j < 3; j++){
			bvec(snodeN + 6 * i + j) += tstep*vTc(3 * i + j);
			bvec(snodeN + 6 * i + 3 + j) += Wc(3 * i + j);

			kmlist.push_back(Triplet<double>(snodeN + 6 * i + j, snodeN + 6 * i + j, 1.0));
			kmlist.push_back(Triplet<double>(snodeN + 6 * i + 3 + j, snodeN + 6 * i + 3 + j, 1.0));
		}
		Eigen::MatrixXd Iner(3, 3);
		Eigen::VectorXd wk(3), Iw(3), wIw(3);
		for (int j = 0; j < 3; j++){
			wk(j) = Wc(3 * i + j);
			for (int k = 0; k < 3; k++){
				Iner(j, k) = Connect[i].Iner[j][k];
			}
		}
		Iw = Iner*wk;
		for (int j = 0; j < 3; j++){
			wIw(j) = wk((j + 1) % 3) * Iw((j + 2) % 3) - wk((j + 2) % 3)*Iw((j + 1) % 3);
		}
		Iw = Iner.inverse() * wIw;
		for (int j = 0; j < 3; j++){
			bvec(snodeN + 6 * i + 3 + j) -= tstep*Iw(j);
		}
	}


	//剛性パラメータ
	double Across = sParam->Hei*sParam->Wid;
	double ks = sParam->Es*Across;
	double kt = sParam->Gt*Across * (sParam->Wid*sParam->Wid + sParam->Hei*sParam->Hei) / 12.0;
	double kben[2][2];
	kben[0][0] = sParam->Eb11*Across*sParam->Wid*sParam->Wid / 12.0;
	kben[1][1] = sParam->Eb22 *Across*sParam->Hei*sParam->Hei / 12.0;
	kben[0][1] = sParam->Eb12*Across*sParam->Hei*sParam->Wid / 12.0;
	kben[1][0] = sParam->Eb12*Across*sParam->Wid*sParam->Hei / 12.0;

	double Estre = 0.0;
	double Ebend = 0.0;
	double Etwis = 0.0;
	double Eext = 0.0;

	for (int i = 0; i < rodN; i++){
		//ストレッチ成分
		for (int j = 0; j < Rods[i].nodeN - 1; j++){
			int vN = Rods[i].nvlist[j];
			double leng = 0.0;
			double leng0 = Rods[i].oLeng[j];
			for (int k = 0; k < 3; k++){
				leng += pow(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k], 2);
			}
			leng = sqrt(leng);

			Estre += 0.5*ks*pow(leng / leng0 - 1.0, 2)*leng0;
			for (int k = 0; k < 3; k++){
				bvec(vN + k) += ks*(1.0 / leng0 - 1.0 / leng)*(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k]) * tstep*tstep / Rods[i].Mass(4 * j + k);
				bvec(vN + 4 + k) -= ks*(1.0 / leng0 - 1.0 / leng)*(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k]) * tstep*tstep / Rods[i].Mass(4 * j + 4 + k);
				resF(vN + k) += ks*(1.0 / leng0 - 1.0 / leng)*(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k]);
				resF(vN + 4 + k) -= ks*(1.0 / leng0 - 1.0 / leng)*(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k]);
				for (int l = 0; l < 3; l++){
					double val = ks / pow(leng, 3)*(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k])*(Rods[i].Nodes[j + 1][l] - Rods[i].Nodes[j][l]);
					if (k == l) val += ks*(1.0 / leng0 - 1.0 / leng);
					val *= tstep*tstep;
					kmlist.push_back(Triplet<double>(vN + k, vN + l, val / Rods[i].Mass(4 * j + k)));
					kmlist.push_back(Triplet<double>(vN + k, vN + 4 + l, -val / Rods[i].Mass(4 * j + k)));
					kmlist.push_back(Triplet<double>(vN + 4 + k, vN + l, -val / Rods[i].Mass(4 * j + 4 + k)));
					kmlist.push_back(Triplet<double>(vN + 4 + k, vN + 4 + l, val / Rods[i].Mass(4 * j + 4 + k)));
				}
			}
		}


		//ねじり・曲げ成分
		for (int j = 1; j < Rods[i].nodeN - 1; j++){
			int vN = Rods[i].nvlist[j];
			double mbar = twist_angle_calculation(Rods[i].TPframe[j - 1], Rods[i].TPframe[j]);
			double tmk = Rods[i].twis(j) - Rods[i].twis(j - 1) + mbar;
			double tm0 = Rods[i].tm[j - 1];
			double kb[3], ki[2];
			double kai = 1.0;
			for (int k = 0; k < 3; k++){
				kb[k] = Rods[i].Mframe[j - 1].d3[(k + 1) % 3] * Rods[i].Mframe[j].d3[(k + 2) % 3]
					- Rods[i].Mframe[j - 1].d3[(k + 2) % 3] * Rods[i].Mframe[j].d3[(k + 1) % 3];
				kai += Rods[i].Mframe[j - 1].d3[k] * Rods[i].Mframe[j].d3[k];
			}
			ki[0] = 0.0;
			ki[1] = 0.0;
			for (int k = 0; k < 3; k++){
				kb[k] = 2.0*kb[k] / kai;
				ki[0] += 0.5*(Rods[i].Mframe[j - 1].d2[k] + Rods[i].Mframe[j].d2[k])*kb[k];
				ki[1] -= 0.5*(Rods[i].Mframe[j - 1].d1[k] + Rods[i].Mframe[j].d1[k])*kb[k];
			}

			double leng1 = 0.0;
			double leng2 = 0.0;
			for (int k = 0; k < 3; k++){
				leng1 += pow(Rods[i].Nodes[j][k] - Rods[i].Nodes[j - 1][k], 2);
				leng2 += pow(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k], 2);
			}
			leng1 = sqrt(leng1);
			leng2 = sqrt(leng2);
			double leni = (leng1 + leng2) / 2.0;
			
			Etwis += 0.5*kt*(tmk - tm0)*(tmk - tm0) / leni;
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 2; l++){
					Ebend += 0.5*(ki[k] - Rods[i].bk[j - 1][k])*kben[k][l]*(ki[l] - Rods[i].bk[j - 1][l]) / leni;
				}
			}

			double grme[2][3], hesm[4][3][3];
			double grke[2][2][3], heske[2][4][3][3];
			double grkt[2][2], hket[2][2][2][3], hktt[2][2];
			Rods[i].derivative_of_twist_bend(j, leng1,leng2, grme, hesm, grke, grkt, heske, hket, hktt);

			for (int k = 0; k < 3; k++){
				bvec(vN - 4 + k) += kt*(tmk - tm0) / leni*grme[0][k] * tstep*tstep / Rods[i].Mass(4 * j - 4 + k);
				bvec(vN + k) += kt*(tmk - tm0) / leni*(-grme[0][k] + grme[1][k])*tstep*tstep / Rods[i].Mass(4 * j + k);
				bvec(vN + 4 + k) += -kt*(tmk - tm0) / leni*grme[1][k] * tstep*tstep / Rods[i].Mass(4 * j + 4 + k);
				resF(vN - 4 + k) += kt*(tmk - tm0) / leni*grme[0][k];
				resF(vN + k) += kt*(tmk - tm0) / leni*(-grme[0][k] + grme[1][k]);
				resF(vN + 4 + k) += -kt*(tmk - tm0) / leni*grme[1][k];
				for (int l = 0; l < 2; l++){
					for (int m = 0; m < 2; m++){
						bvec(vN - 4 + k) += grke[l][0][k] * kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni*tstep*tstep / Rods[i].Mass(4 * j - 4 + k);
						bvec(vN + k) += (-grke[l][0][k] + grke[l][1][k])*kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni*tstep*tstep / Rods[i].Mass(4 * j + k);
						bvec(vN + 4 + k) += -grke[l][1][k] * kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni*tstep*tstep / Rods[i].Mass(4 * j + 4 + k);
						resF(vN - 4 + k) += grke[l][0][k] * kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
						resF(vN + k) += (-grke[l][0][k] + grke[l][1][k])*kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
						resF(vN + 4 + k) += -grke[l][1][k] * kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
					}
				}
			}
			bvec(vN - 1) += kt*(tmk - tm0) / leni*tstep*tstep / Rods[i].Mass(4 * j - 1);
			bvec(vN + 3) -= kt*(tmk - tm0) / leni*tstep*tstep / Rods[i].Mass(4 * j + 3);
			resF(vN - 1) += kt*(tmk - tm0) / leni;
			resF(vN + 3) -= kt*(tmk - tm0) / leni;
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 2; l++){
					bvec(vN - 1) -= grkt[k][0] * kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni*tstep*tstep / Rods[i].Mass(4 * j - 1);
					bvec(vN + 3) -= grkt[k][1] * kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni*tstep*tstep / Rods[i].Mass(4 * j + 3);
					resF(vN - 1) -= grkt[k][0] * kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni;
					resF(vN + 3) -= grkt[k][1] * kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni;
				}
			}
			for (int k = 0; k < 3; k++){
				for (int l = 0; l < 3; l++){
					for (int m = 0; m < 2; m++){
						for (int n = 0; n < 2; n++){
							double val = kt*(grme[m][k] * grme[n][l] + (tmk - tm0)*hesm[2 * m + n][k][l]) / leni;
							val *= tstep*tstep;
							kmlist.push_back(Triplet<double>(vN + 4 * (m - 1) + k, vN + 4 * (n - 1) + l, val / Rods[i].Mass(4 * (j - 1 + m) + k)));
							kmlist.push_back(Triplet<double>(vN + 4 * (m - 1) + k, vN + 4 * n + l, -val / Rods[i].Mass(4 * (j - 1 + m) + k)));
							kmlist.push_back(Triplet<double>(vN + 4 * m + k, vN + 4 * (n - 1) + l, -val / Rods[i].Mass(4 * (j + m) + k)));
							kmlist.push_back(Triplet<double>(vN + 4 * m + k, vN + 4 * n + l, val / Rods[i].Mass(4 * (j + m) + k)));

							val = 0.0;
							for (int o = 0; o < 2; o++){
								for (int p = 0; p < 2; p++){
									val += (grke[o][m][k] * kben[o][p] * grke[p][n][l] + heske[o][2 * m + n][k][l] * kben[o][p] * (ki[p] - Rods[i].bk[j - 1][p])) / leni;
								}
							}
							val *= tstep*tstep;
							kmlist.push_back(Triplet<double>(vN + 4 * (m - 1) + k, vN + 4 * (n - 1) + l, val / Rods[i].Mass(4 * (j - 1 + m) + k)));
							kmlist.push_back(Triplet<double>(vN + 4 * (m - 1) + k, vN + 4 * n + l, -val / Rods[i].Mass(4 * (j - 1 + m) + k)));
							kmlist.push_back(Triplet<double>(vN + 4 * m + k, vN + 4 * (n - 1) + l, -val / Rods[i].Mass(4 * (j + m) + k)));
							kmlist.push_back(Triplet<double>(vN + 4 * m + k, vN + 4 * n + l, val / Rods[i].Mass(4 * (j + m) + k)));
						}
					}
				}
			}
			kmlist.push_back(Triplet<double>(vN - 1, vN - 1, kt / leni*tstep*tstep / Rods[i].Mass(4 * j - 1)));
			kmlist.push_back(Triplet<double>(vN - 1, vN + 3, -kt / leni*tstep*tstep / Rods[i].Mass(4 * j - 1)));
			kmlist.push_back(Triplet<double>(vN + 3, vN - 1, -kt / leni*tstep*tstep / Rods[i].Mass(4 * j + 3)));
			kmlist.push_back(Triplet<double>(vN + 3, vN + 3, kt / leni*tstep*tstep / Rods[i].Mass(4 * j + 3)));
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 2; l++){
					float val = 0.0;
					for (int m = 0; m < 2; m++){
						for (int n = 0; n < 2; n++){
							val += grkt[m][k] * kben[m][n] * grkt[n][l] / leni;
							if (k == l) val += hktt[m][k] * kben[m][n] * (ki[n] - Rods[i].bk[j - 1][n]) / leni;
						}
					}
					kmlist.push_back(Triplet<double>(vN + 4 * (k - 1) + 3, vN + 4 * (l - 1) + 3, val*tstep*tstep / Rods[i].Mass(4 * (j - 1 + k) + 3)));
				}
			}
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 3; l++){
					double val = kt*grme[k][l] / leni*tstep*tstep;
					kmlist.push_back(Triplet<double>(vN + 4 * (k - 1) + l, vN - 1, val / Rods[i].Mass(4 * (j - 1 + k) + l)));
					kmlist.push_back(Triplet<double>(vN + 4 * k + l, vN - 1, -val / Rods[i].Mass(4 * (j + k) + l)));
					kmlist.push_back(Triplet<double>(vN + 4 * (k - 1) + l, vN + 3, -val / Rods[i].Mass(4 * (j - 1 + k) + l)));
					kmlist.push_back(Triplet<double>(vN + 4 * k + l, vN + 3, val / Rods[i].Mass(4 * (j + k) + l)));

					kmlist.push_back(Triplet<double>(vN - 1, vN + 4 * (k - 1) + l, val / Rods[i].Mass(4 * j - 1)));
					kmlist.push_back(Triplet<double>(vN - 1, vN + 4 * k + l, -val / Rods[i].Mass(4 * j - 1)));
					kmlist.push_back(Triplet<double>(vN + 3, vN + 4 * (k - 1) + l, -val / Rods[i].Mass(4 * j + 3)));
					kmlist.push_back(Triplet<double>(vN + 3, vN + 4 * k + l, val / Rods[i].Mass(4 * j + 3)));
				}
				for (int l = 0; l < 2; l++){
					for (int m = 0; m < 3; m++){
						double val = 0.0;
						for (int n = 0; n < 2; n++){
							for (int o = 0; o < 2; o++){
								val += (grke[n][k][m] * kben[n][o] * grkt[o][l] + hket[n][k][l][m] * kben[n][o] * (ki[o] - Rods[i].bk[j - 1][o])) / leni;
							}
						}
						kmlist.push_back(Triplet<double>(vN + 4 * (k - 1) + m, vN + 4 * (l - 1) + 3, -val*tstep*tstep / Rods[i].Mass(4 * (j - 1 + k) + m)));
						kmlist.push_back(Triplet<double>(vN + 4 * k + m, vN + 4 * (l - 1) + 3, val*tstep*tstep / Rods[i].Mass(4 * (j + k) + m)));
						kmlist.push_back(Triplet<double>(vN + 4 * (l - 1) + 3, vN + 4 * (k - 1) + m, -val*tstep*tstep / Rods[i].Mass(4 * (j - 1 + l) + 3)));
						kmlist.push_back(Triplet<double>(vN + 4 * (l - 1) + 3, vN + 4 * k + m, val*tstep*tstep / Rods[i].Mass(4 * (j - 1 + l) + 3)));
					}
				}
			}

		}
	}

	for (int i = 0; i < snodeN; i++){
		kmlist.push_back(Triplet<double>(i, i, 1.0));
	}


	//拘束関連
	for (int i = 0; i < conN; i++){
		int vN = cvlist[i];
		for (int j = 0; j < Connect[i].nodeN; j++){
			int rN = Connect[i].cNode[j].rN;
			int nN = Connect[i].cNode[j].nN;
			int cvn = Rods[rN].nvlist[nN];
			double Cp[3],crossX[3][3];
			for (int k = 0; k < 3; k++){
				Cp[k] = Rods[rN].Nodes[nN][k] - Connect[i].gp[k];
				for (int l = 0; l < 3; l++){
					if (k != l && l == (k + 1) % 3) crossX[k][l] = -Connect[i].cNode[j].pos[3 - k - l];
					else if (k != l) crossX[k][l] = Connect[i].cNode[j].pos[3 - k - l];
					else crossX[k][l] = 0.0;
				}
			}
			
			Eigen::MatrixXf grC(3, 3), Iner(3, 3), ingC(3, 3);
			double Rot0[3][3];
			derive_rot_mat_from_quaternion(Connect[i].quat, Rot0);
			for (int k = 0; k < 3; k++){
				for (int l = 0; l < 3; l++){
					Cp[k] -= Rot0[k][l] * Connect[i].cNode[j].pos[l];
					Iner(k, l) = Connect[i].Iner[k][l];
					grC(k, l) = 0.0;
					for (int m = 0; m < 3; m++){
						grC(k, l) -= Rot0[k][m] * crossX[l][m]*tstep;
					}
				}
			}
			ingC = Iner.inverse()*grC.transpose();
			for (int k = 0; k < 3; k++){
				bvec(vN + 3 * j + k) -= Cp[k];
			}

			//printf("%d %d %f %f %f\n", i, j, Cp[0], Cp[1], Cp[2]);

			for (int k = 0; k < 3; k++){
				kmlist.push_back(Triplet<double>(cvn + k, vN + 3 * j + k, tstep*tstep / Rods[rN].Mass(4 * nN + k)));
				kmlist.push_back(Triplet<double>(vN + 3 * j + k, cvn + k, 1.0));
				kmlist.push_back(Triplet<double>(snodeN + 6 * i + k, vN + 3 * j + k, -tstep*tstep / Connect[i].Mass));
				kmlist.push_back(Triplet<double>(vN + 3 * j + k, snodeN + 6 * i + k, -1.0));
				for (int l = 0; l < 3; l++){
					kmlist.push_back(Triplet<double>(snodeN + 6 * i + 3 + k, vN + 3 * j + l, tstep*ingC(k, l)));
					kmlist.push_back(Triplet<double>(vN + 3 * j + k, snodeN + 6 * i + 3 + l, grC(k, l)));
				}
			}
		}
		vN += 3*Connect[i].nodeN;
		for (int j = 0; j < Connect[i].edgeN; j++){
			int rN = Connect[i].cEdge[j].rN;
			int eN = Connect[i].cEdge[j].eN;
			int cvn = Rods[rN].evlist[eN];
			double Cf = 0.0;

			double Rot0[3][3], Rd1[3];
			derive_rot_mat_from_quaternion(Connect[i].quat, Rot0);

			Eigen::MatrixXd Iner(3, 3), invI(3, 3);
			for (int k = 0; k < 3; k++){
				Rd1[k] = 0.0;
				for (int l = 0; l < 3; l++){
					Rd1[k] += Rot0[k][l] * Connect[i].cEdge[j].frame[0][l];
					Iner(k, l) = Connect[i].Iner[k][l];
				}
			}
			invI = Iner.inverse();
			double seki = 0.0;
			for (int k = 0; k < 3; k++){
				seki += Rd1[k] * Rods[rN].Mframe[eN].d3[k];
			}
			double leng = 0.0;
			for (int k = 0; k < 3; k++){
				Rd1[k] -= seki*Rods[rN].Mframe[eN].d3[k];
				leng += Rd1[k] * Rd1[k];
			}
			leng = sqrt(leng);
			if (leng >1.0e-7){
				for (int k = 0; k < 3; k++){
					Rd1[k] /= leng;
				}
				double theta = 0.0;
				double nd[3];
				for (int k = 0; k < 3; k++){
					theta += Rd1[k] * Rods[rN].Mframe[eN].d1[k];
					nd[k] = Rd1[(k + 1) % 3] * Rods[rN].Mframe[eN].d1[(k + 2) % 3]
						- Rd1[(k + 2) % 3] * Rods[rN].Mframe[eN].d1[(k + 1) % 3];
				}
				seki = 0.0;
				for (int k = 0; k < 3; k++){
					seki += nd[k] * Rods[rN].Mframe[eN].d3[k];
				}
				if (theta > 1.0) theta = 1.0;
				if (theta < -1.0) theta = -1.0;
				if (seki > 0) Cf = acos(theta);
				else Cf = -acos(theta);
			}
			bvec(vN + j) -= Cf;

			//printf("%d %d %f\n",i,j,Cf);

			double igC[3];
			for (int k = 0; k < 3; k++){
				igC[k] = 0.0;
				for (int l = 0; l < 3; l++){
					igC[k] += -tstep*invI(k, l)*Connect[i].cEdge[j].frame[2][l];
				}
			}

			kmlist.push_back(Triplet<double>(cvn, vN + j, tstep*tstep / Rods[rN].Mass(4 * eN + 3)));
			kmlist.push_back(Triplet<double>(vN + j, cvn, 1.0));
			for (int k = 0; k < 3; k++){
				kmlist.push_back(Triplet<double>(snodeN + 6 * i + 3 + k, vN + j, tstep*igC[k]));
				kmlist.push_back(Triplet<double>(vN + j, snodeN + 6 * i + 3 + k, -tstep*Connect[i].cEdge[j].frame[2][k]));
			}
		}
	}

	double wei = 1.0e-9;
	for (int i = snodeN + 6 * conN; i < svarN; i++){
		kmlist.push_back(Triplet<double>(i, i, wei));
	}

	Kmat.setFromTriplets(kmlist.begin(), kmlist.end());

	
	//端点拘束
	for (int i = 0; i < 1; i++){
		int vN1 = Rods[i].nvlist[0];
		for (int j = 0; j < 4; j++){
			for (SparseMatrix<double>::InnerIterator it(Kmat, vN1 + j); it; ++it){
				if (it.row() != vN1+j) it.valueRef() = 0.0;
				else it.valueRef() = 1.0;
			}
			bvec(vN1 + j) = 0.0;
		}

		vN1 = Rods[i].nvlist[1];
		for (int j = 0; j < 3; j++){
			for (SparseMatrix<double>::InnerIterator it(Kmat, vN1 + j); it; ++it){
				if (it.row() != vN1 + j) it.valueRef() = 0.0;
				else it.valueRef() = 1.0;
			}
			bvec(vN1 + j) = 0.0;
		}
	}
	for (int i = 0; i < svarN; i++){
		for (SparseMatrix<double>::InnerIterator it(Kmat, i); it; ++it){
			for (int j = 0; j < 1; j++){
				int vN1 = Rods[j].nvlist[0];
				if (it.row() >= vN1 && it.row() < vN1 + 4 && it.row() != i) it.valueRef() = 0.0;

				vN1 = Rods[j].nvlist[1];
				if (it.row() >= vN1 && it.row() < vN1 + 3 && it.row() != i) it.valueRef() = 0.0;
			}
		}
	}



	//行列計算
	SparseLU<SparseMatrix<double>> sLU;
	sLU.compute(Kmat);
	bool csuccess = (sLU.info() == Success);
	if (csuccess){
		sol = sLU.solve(bvec);
		Eigen::VectorXd dif = Kmat*sol - bvec;
		if (dif.norm() > 1.0e-4)csuccess = false;
	}
	if (!csuccess){
		bool CSuccess = false;
		double addmass = 1.0;
		printf("sparse LU failed\n");
		while (1){
			for (int i = 0; i < snodeN + 6 * conN; i++){
				for (SparseMatrix<double>::InnerIterator it(Kmat, i); it; ++it){
					if (i == it.row()){
						bool is_constrained = false;
						for (int j = 0; j < 1; j++){
							int vN1 = Rods[j].nvlist[0];
							if (it.row() >= vN1 && it.row() < vN1 + 4) is_constrained = true;
							vN1 = Rods[j].nvlist[1];
							if (it.row() >= vN1 && it.row() < vN1 + 3) is_constrained = true;
						}
						if (!is_constrained){
							it.valueRef() += addmass;
						}
					}
				}
			}
			sLU.compute(Kmat);
			CSuccess = (sLU.info() == Success);
			if (CSuccess) {
				sol = sLU.solve(bvec);
				Eigen::VectorXd dif = Kmat*sol - bvec;
				if (dif.norm() < 1.0e-4) break;
				CSuccess = false;
			}
			addmass *= 2.0;
			if (addmass > 1.0e5) break;
		}
		printf("add:%f\n",addmass);
		//if (addmass > 1.0e3) tstep_is_large = true;
		if (!CSuccess){
			printf("Simulation Diverge\n");
			is_diverged = true;
			tstep_is_large = true;
			return;
		}
	}

	//Eigen::VectorXd dif = Kmat*sol - bvec;
	//cout <<"norm:" << dif.norm() <<"("<<bvec.norm()<<")"<<endl;


	double Epot = 0.0;
	double Nmove_max = 0.0;
	double Erot_max = 0.0;
	for (int i = 0; i < rodN; i++){
		for (int j = 0; j < Rods[i].nodeN; j++){
			int vN = Rods[i].nvlist[j];
			double move = 0.0;
			for (int k = 0; k < 3; k++){
				Rods[i].Nodes[j][k] += sol(vN + k);
				Rods[i].vXk(3 * j + k) = sol(vN + k) / tstep;
				move += sol(vN + k)*sol(vN + k);
			}
			move = sqrt(move);
			if (move > Nmove_max) Nmove_max = move;

			if (j == Rods[i].nodeN - 1) continue;
			Rods[i].twis(j) += sol(vN + 3);
			Rods[i].vTk(j) = sol(vN + 3) / tstep;

			if (Erot_max < fabs(sol(vN + 3))) Erot_max = sol(vN + 3);
		}
	}
	double Rmove_max = 0.0;
	double Rrot_max = 0.0;
	double lamV = 0.0;
	for (int i = 0; i < conN; i++){
		double move = 0.0;
		for (int j = 0; j < 3; j++){
			Connect[i].gp[j] += sol(snodeN + 6 * i + j);
			vTc(3 * i + j) = sol(snodeN + 6 * i + j) / tstep;
			move += sol(snodeN + 6 * i + j)*sol(snodeN + 6 * i + j);
		}
		move = sqrt(move);
		if (move > Rmove_max) Rmove_max = move;
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				Epot += 0.5*Connect[i].Iner[j][k] * sol(snodeN + 6 * i + 3 + j)*sol(snodeN + 6 * i + 3 + k);
			}
		}
		
		double theta = 0.0;
		double avel[3];
		for (int j = 0; j < 3; j++){
			Wc(3 * i + j) = sol(snodeN + 6 * i + 3 + j);
			avel[j] = sol(snodeN + 6 * i + 3 + j);
			theta += avel[j] * avel[j];
		}
		theta = sqrt(theta);
		if (theta > 1.0e-5){
			for (int j = 0; j < 3; j++){
				avel[j] /= theta;
			}
		}
		else {
			avel[0] = 1.0; avel[1] = 0.0; avel[2] = 0.0;
		}
		theta *= tstep;
		if (Rrot_max < fabs(theta)) Rrot_max = fabs(theta);

		double quat[4];
		quat[0] = cos(theta / 2.0);
		for (int j = 0; j < 3; j++){
			quat[1 + j] = sin(theta / 2.0)*avel[j];
		}
		double Rot[3][3], Rot0[3][3], nRot[3][3];
		derive_rot_mat_from_quaternion(quat, Rot);
		derive_rot_mat_from_quaternion(Connect[i].quat, Rot0);
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				nRot[j][k] = 0.0;
				for (int l = 0; l < 3; l++){
					nRot[j][k] += Rot0[j][l] * Rot[l][k];
				}
			}
		}
		derive_quaternion_from_rot_mat(nRot, Connect[i].quat);

		//resF用
		int vN = cvlist[i];
		for (int j = 0; j < Connect[i].nodeN; j++){
			int rN = Connect[i].cNode[j].rN;
			int nN = Connect[i].cNode[j].nN;
			for (int k = 0; k < 3; k++){
				lamV += pow(preSol(vN + 3 * j + k - snodeN - 6 * conN) - sol(vN + 3 * j + k), 2);
				preSol(vN + 3 * j + k - snodeN - 6 * conN) = sol(vN + 3 * j + k);
				resF(Rods[rN].nvlist[nN] + k) -= preSol(vN + 3 * j + k - snodeN - 6 * conN);
				//resF(Rods[rN].nvlist[nN] + k) = 0.0;
			}
		}
		vN += 3 * Connect[i].nodeN;
		for (int j = 0; j < Connect[i].edgeN; j++){
			int rN = Connect[i].cEdge[j].rN;
			int eN = Connect[i].cEdge[j].eN;
			lamV += pow(preSol(vN + j - snodeN - 6 * conN) - sol(vN + j), 2);
			preSol(vN + j - snodeN - 6 * conN) = sol(vN + j);
			resF(Rods[rN].evlist[eN]) -= preSol(vN + j - snodeN - 6 * conN);
			//resF(Rods[rN].evlist[eN]) = 0.0;
		}
	}

	printf("Energy:%f (%f %f %f) \n",Estre+Ebend+Etwis,Estre,Ebend,Etwis);
	printf("Potential:%f (%f %f %f %f)\n", Epot, Nmove_max, Erot_max*180.0 / M_PI, Rmove_max, Rrot_max*180.0 / M_PI);
	cout << "F res:" << resF.norm() <<" "<< preSol.norm()<<" ("<<sqrt(lamV)<<")"<< endl;
	if (Nmove_max < sParam->Wid*0.01 && Erot_max < M_PI*0.1 / 180.0 && Rmove_max < sParam->Wid*0.01 && Rrot_max < M_PI*0.1 / 180.0){
		is_converged = true;
	}
	if (old_E > 0.0 && Estre + Ebend + Etwis > initE * 100){
		is_diverged = true;
	}

	//kinematic dumping
	//old_E = Estre + Ebend + Etwis + 1.0;
	if (old_E < 0.0) initE = Estre + Ebend + Etwis;
	if (old_E >= 0.0 && old_E > Estre + Ebend + Etwis){
		for (int i = 0; i < rodN; i++){
			for (int j = 0; j < Rods[i].nodeN; j++){
				for (int k = 0; k < 3; k++){
					Rods[i].vXk(3 * j + k) = 0.0;
				}
				if (j == Rods[i].nodeN - 1) continue;
				Rods[i].vTk(j) = 0.0;
			}
		}
		for (int i = 0; i < conN; i++){
			for (int j = 0; j < 3; j++){
				vTc(3 * i + j) = 0.0;
				Wc(3 * i + j) = 0.0;
			}
		}

	}
	old_E = Estre + Ebend + Etwis;

	for (int i = 0; i < rodN; i++){
		Rods[i].update_frame();
	}


}


//準静的に解く
bool RodNetwork::simulateFull(){

	Eigen::VectorXd bvec(svarN), sol(svarN), resF(snodeN);

	int RnodeN = 0;
	for (int i = 0; i < 3; i++){
		Gpos[i] = 0.0; 
	}
	for (int i = 0; i < rodN; i++){
		for (int j = 0; j < Rods[i].nodeN; j++){
			for (int k = 0; k < 3; k++){
				Gpos[k] += Rods[i].Nodes[j][k];
			}
			RnodeN++;
		}
	}
	for (int i = 0; i < 3; i++){
		Gpos[i] /= (double)RnodeN;
	}

	double Ldamp = 0.0;
	double LamSize = -1.0;
	bool is_converged = false;
	int iterN = 0;
	int maxIter = 500;
	while (1){
		//剛性マトリックス・残差ベクトルの計算
		Eigen::SparseMatrix<double> Kmat(svarN, svarN);
		std::vector<Triplet<double>> kmlist;
		bvec = Eigen::VectorXd::Zero(svarN);
		resF = Eigen::VectorXd::Zero(snodeN);

		/*for (int i = 0; i < conN; i++){
			bvec(snodeN + 6 * i + 2) -= Connect[i].Mass * 9800;
		}*/

		//剛性パラメータ
		double Across = sParam->Hei*sParam->Wid;
		double ks = sParam->Es*Across;
		double kt = sParam->Gt*Across * (sParam->Wid*sParam->Wid + sParam->Hei*sParam->Hei) / 12.0;
		double kben[2][2];
		kben[0][0] = sParam->Eb11*Across*sParam->Wid*sParam->Wid / 12.0;
		kben[1][1] = sParam->Eb22 *Across*sParam->Hei*sParam->Hei / 12.0;
		kben[0][1] = sParam->Eb12*Across*sParam->Hei*sParam->Wid / 12.0;
		kben[1][0] = sParam->Eb12*Across*sParam->Wid*sParam->Hei / 12.0;

		double Estre = 0.0;
		double Ebend = 0.0;
		double Etwis = 0.0;
		double Eext = 0.0;
		double Econ = 0.0; 

		double varCf = 0.0;
		double varCp = 0.0;

		for (int i = 0; i < rodN; i++){
			//ストレッチ成分
			for (int j = 0; j < Rods[i].nodeN - 1; j++){
				int vN = Rods[i].nvlist[j];
				double leng = 0.0;
				double leng0 = Rods[i].oLeng[j];
				for (int k = 0; k < 3; k++){
					leng += pow(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k], 2);
				}
				leng = sqrt(leng);

				Estre += 0.5*ks*pow(leng / leng0 - 1.0, 2)*leng0;
				for (int k = 0; k < 3; k++){
					bvec(vN + k) += ks*(1.0 / leng0 - 1.0 / leng)*(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k]);
					bvec(vN + 4 + k) -= ks*(1.0 / leng0 - 1.0 / leng)*(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k]);
					resF(vN + k) += ks*(1.0 / leng0 - 1.0 / leng)*(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k]);
					resF(vN + 4 + k) -= ks*(1.0 / leng0 - 1.0 / leng)*(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k]);
					for (int l = 0; l < 3; l++){
						double val = ks / pow(leng, 3)*(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k])*(Rods[i].Nodes[j + 1][l] - Rods[i].Nodes[j][l]);
						if (k == l) val += ks*(1.0 / leng0 - 1.0 / leng);
						kmlist.push_back(Triplet<double>(vN + k, vN + l, val));
						kmlist.push_back(Triplet<double>(vN + k, vN + 4 + l, -val));
						kmlist.push_back(Triplet<double>(vN + 4 + k, vN + l, -val));
						kmlist.push_back(Triplet<double>(vN + 4 + k, vN + 4 + l, val));
					}
				}
			}


			//ねじり・曲げ成分
			for (int j = 1; j < Rods[i].nodeN - 1; j++){
				int vN = Rods[i].nvlist[j];
				double mbar = twist_angle_calculation(Rods[i].TPframe[j - 1], Rods[i].TPframe[j]);
				double tmk = Rods[i].twis(j) - Rods[i].twis(j - 1) + mbar;
				double tm0 = Rods[i].tm[j - 1];
				double kb[3], ki[2];
				double kai = 1.0;
				for (int k = 0; k < 3; k++){
					kb[k] = Rods[i].Mframe[j - 1].d3[(k + 1) % 3] * Rods[i].Mframe[j].d3[(k + 2) % 3]
						- Rods[i].Mframe[j - 1].d3[(k + 2) % 3] * Rods[i].Mframe[j].d3[(k + 1) % 3];
					kai += Rods[i].Mframe[j - 1].d3[k] * Rods[i].Mframe[j].d3[k];
				}
				ki[0] = 0.0;
				ki[1] = 0.0;
				for (int k = 0; k < 3; k++){
					kb[k] = 2.0*kb[k] / kai;
					ki[0] += 0.5*(Rods[i].Mframe[j - 1].d2[k] + Rods[i].Mframe[j].d2[k])*kb[k];
					ki[1] -= 0.5*(Rods[i].Mframe[j - 1].d1[k] + Rods[i].Mframe[j].d1[k])*kb[k];
				}

				double leng1 = 0.0;
				double leng2 = 0.0;
				for (int k = 0; k < 3; k++){
					leng1 += pow(Rods[i].Nodes[j][k] - Rods[i].Nodes[j - 1][k], 2);
					leng2 += pow(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k], 2);
				}
				leng1 = sqrt(leng1);
				leng2 = sqrt(leng2);
				double leni = (leng1 + leng2) / 2.0;

				Etwis += 0.5*kt*(tmk - tm0)*(tmk - tm0) / leni;
				for (int k = 0; k < 2; k++){
					for (int l = 0; l < 2; l++){
						Ebend += 0.5*(ki[k] - Rods[i].bk[j - 1][k])*kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni;
					}
				}

				double grme[2][3], hesm[4][3][3];
				double grke[2][2][3], heske[2][4][3][3];
				double grkt[2][2], hket[2][2][2][3], hktt[2][2];
				Rods[i].derivative_of_twist_bend(j, leng1, leng2, grme, hesm, grke, grkt, heske, hket, hktt);

				for (int k = 0; k < 3; k++){
					bvec(vN - 4 + k) += kt*(tmk - tm0) / leni*grme[0][k];
					bvec(vN + k) += kt*(tmk - tm0) / leni*(-grme[0][k] + grme[1][k]);
					bvec(vN + 4 + k) += -kt*(tmk - tm0) / leni*grme[1][k];
					resF(vN - 4 + k) += kt*(tmk - tm0) / leni*grme[0][k];
					resF(vN + k) += kt*(tmk - tm0) / leni*(-grme[0][k] + grme[1][k]);
					resF(vN + 4 + k) += -kt*(tmk - tm0) / leni*grme[1][k];
					for (int l = 0; l < 2; l++){
						for (int m = 0; m < 2; m++){
							bvec(vN - 4 + k) += grke[l][0][k] * kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
							bvec(vN + k) += (-grke[l][0][k] + grke[l][1][k])*kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
							bvec(vN + 4 + k) += -grke[l][1][k] * kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
							resF(vN - 4 + k) += grke[l][0][k] * kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
							resF(vN + k) += (-grke[l][0][k] + grke[l][1][k])*kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
							resF(vN + 4 + k) += -grke[l][1][k] * kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
						}
					}
				}
				bvec(vN - 1) += kt*(tmk - tm0) / leni;
				bvec(vN + 3) -= kt*(tmk - tm0) / leni;
				resF(vN - 1) += kt*(tmk - tm0) / leni;
				resF(vN + 3) -= kt*(tmk - tm0) / leni;
				for (int k = 0; k < 2; k++){
					for (int l = 0; l < 2; l++){
						bvec(vN - 1) -= grkt[k][0] * kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni;
						bvec(vN + 3) -= grkt[k][1] * kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni;
						resF(vN - 1) -= grkt[k][0] * kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni;
						resF(vN + 3) -= grkt[k][1] * kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni;
					}
				}
				for (int k = 0; k < 3; k++){
					for (int l = 0; l < 3; l++){
						for (int m = 0; m < 2; m++){
							for (int n = 0; n < 2; n++){
								double val = kt*(grme[m][k] * grme[n][l] + (tmk - tm0)*hesm[2 * m + n][k][l]) / leni;
								kmlist.push_back(Triplet<double>(vN + 4 * (m - 1) + k, vN + 4 * (n - 1) + l, val));
								kmlist.push_back(Triplet<double>(vN + 4 * (m - 1) + k, vN + 4 * n + l, -val));
								kmlist.push_back(Triplet<double>(vN + 4 * m + k, vN + 4 * (n - 1) + l, -val));
								kmlist.push_back(Triplet<double>(vN + 4 * m + k, vN + 4 * n + l, val));

								val = 0.0;
								for (int o = 0; o < 2; o++){
									for (int p = 0; p < 2; p++){
										val += (grke[o][m][k] * kben[o][p] * grke[p][n][l] + heske[o][2 * m + n][k][l] * kben[o][p] * (ki[p] - Rods[i].bk[j - 1][p])) / leni;
									}
								}
								kmlist.push_back(Triplet<double>(vN + 4 * (m - 1) + k, vN + 4 * (n - 1) + l, val));
								kmlist.push_back(Triplet<double>(vN + 4 * (m - 1) + k, vN + 4 * n + l, -val));
								kmlist.push_back(Triplet<double>(vN + 4 * m + k, vN + 4 * (n - 1) + l, -val));
								kmlist.push_back(Triplet<double>(vN + 4 * m + k, vN + 4 * n + l, val));
							}
						}
					}
				}
				kmlist.push_back(Triplet<double>(vN - 1, vN - 1, kt / leni));
				kmlist.push_back(Triplet<double>(vN - 1, vN + 3, -kt / leni));
				kmlist.push_back(Triplet<double>(vN + 3, vN - 1, -kt / leni));
				kmlist.push_back(Triplet<double>(vN + 3, vN + 3, kt / leni));
				for (int k = 0; k < 2; k++){
					for (int l = 0; l < 2; l++){
						float val = 0.0;
						for (int m = 0; m < 2; m++){
							for (int n = 0; n < 2; n++){
								val += grkt[m][k] * kben[m][n] * grkt[n][l] / leni;
								if (k == l) val += hktt[m][k] * kben[m][n] * (ki[n] - Rods[i].bk[j - 1][n]) / leni;
							}
						}
						kmlist.push_back(Triplet<double>(vN + 4 * (k - 1) + 3, vN + 4 * (l - 1) + 3, val));
					}
				}
				for (int k = 0; k < 2; k++){
					for (int l = 0; l < 3; l++){
						double val = kt*grme[k][l] / leni;
						kmlist.push_back(Triplet<double>(vN + 4 * (k - 1) + l, vN - 1, val ));
						kmlist.push_back(Triplet<double>(vN + 4 * k + l, vN - 1, -val));
						kmlist.push_back(Triplet<double>(vN + 4 * (k - 1) + l, vN + 3, -val));
						kmlist.push_back(Triplet<double>(vN + 4 * k + l, vN + 3, val));

						kmlist.push_back(Triplet<double>(vN - 1, vN + 4 * (k - 1) + l, val));
						kmlist.push_back(Triplet<double>(vN - 1, vN + 4 * k + l, -val));
						kmlist.push_back(Triplet<double>(vN + 3, vN + 4 * (k - 1) + l, -val));
						kmlist.push_back(Triplet<double>(vN + 3, vN + 4 * k + l, val));
					}
					for (int l = 0; l < 2; l++){
						for (int m = 0; m < 3; m++){
							double val = 0.0;
							for (int n = 0; n < 2; n++){
								for (int o = 0; o < 2; o++){
									val += (grke[n][k][m] * kben[n][o] * grkt[o][l] + hket[n][k][l][m] * kben[n][o] * (ki[o] - Rods[i].bk[j - 1][o])) / leni;
								}
							}
							kmlist.push_back(Triplet<double>(vN + 4 * (k - 1) + m, vN + 4 * (l - 1) + 3, -val));
							kmlist.push_back(Triplet<double>(vN + 4 * k + m, vN + 4 * (l - 1) + 3, val));
							kmlist.push_back(Triplet<double>(vN + 4 * (l - 1) + 3, vN + 4 * (k - 1) + m, -val));
							kmlist.push_back(Triplet<double>(vN + 4 * (l - 1) + 3, vN + 4 * k + m, val));
						}
					}
				}

			}
		}

		//拘束関連
		int cpN = 0.0;
		int cfN = 0.0;
		for (int i = 0; i < conN; i++){
			int vN = cvlist[i];
			cpN += Connect[i].nodeN;
			for (int j = 0; j < Connect[i].nodeN; j++){
				int rN = Connect[i].cNode[j].rN;
				int nN = Connect[i].cNode[j].nN;
				int cvn = Rods[rN].nvlist[nN];
				double Cp[3], crossX[3][3];
				for (int k = 0; k < 3; k++){
					Cp[k] = Rods[rN].Nodes[nN][k] - Connect[i].gp[k];
					for (int l = 0; l < 3; l++){
						if (k != l && l == (k + 1) % 3) crossX[k][l] = -Connect[i].cNode[j].pos[3 - k - l];
						else if (k != l) crossX[k][l] = Connect[i].cNode[j].pos[3 - k - l];
						else crossX[k][l] = 0.0;
					}
				}

				Eigen::MatrixXf grC(3, 3);
				double Rot0[3][3];
				derive_rot_mat_from_quaternion(Connect[i].quat, Rot0);
				for (int k = 0; k < 3; k++){
					for (int l = 0; l < 3; l++){
						Cp[k] -= Rot0[k][l] * Connect[i].cNode[j].pos[l];
						grC(k, l) = 0.0;
						for (int m = 0; m < 3; m++){
							grC(k, l) -= Rot0[k][m] * crossX[l][m];
						}
					}
				}
				varCp += sqrt(Cp[0] * Cp[0] + Cp[1] * Cp[1] + Cp[2] * Cp[2]);

				for (int k = 0; k < 3; k++){
					Econ += Cp[k] * Lamd(vN - cvlist[0] + 3 * j + k);
					bvec(vN + 3 * j + k) -= Cp[k];
					bvec(cvn + k) -= 1.0*Lamd(vN - cvlist[0] + 3 * j + k);
					bvec(snodeN + 6 * i + k) -= -1.0*Lamd(vN - cvlist[0] + 3 * j + k);
					for (int l = 0; l < 3; l++){
						bvec(snodeN + 6 * i + 3 + k) -= grC(l, k)*Lamd(vN - cvlist[0] + 3 * j + l);
					}
				}

				for (int k = 0; k < 3; k++){
					kmlist.push_back(Triplet<double>(cvn + k, vN + 3 * j + k, 1.0));
					kmlist.push_back(Triplet<double>(vN + 3 * j + k, cvn + k, 1.0));
					kmlist.push_back(Triplet<double>(snodeN + 6 * i + k, vN + 3 * j + k, -1.0));
					kmlist.push_back(Triplet<double>(vN + 3 * j + k, snodeN + 6 * i + k, -1.0));
					for (int l = 0; l < 3; l++){
						kmlist.push_back(Triplet<double>(snodeN + 6 * i + 3 + k, vN + 3 * j + l, grC(l, k)));
						kmlist.push_back(Triplet<double>(vN + 3 * j + k, snodeN + 6 * i + 3 + l, grC(k, l)));
					}
				}
			}
			vN += 3 * Connect[i].nodeN;
			cfN += Connect[i].edgeN;
			for (int j = 0; j < Connect[i].edgeN; j++){
				int rN = Connect[i].cEdge[j].rN;
				int eN = Connect[i].cEdge[j].eN;
				int cvn = Rods[rN].evlist[eN];
				double Cf = 0.0;

				double Rot0[3][3], Rd1[3];
				derive_rot_mat_from_quaternion(Connect[i].quat, Rot0);

				for (int k = 0; k < 3; k++){
					Rd1[k] = 0.0;
					for (int l = 0; l < 3; l++){
						Rd1[k] += Rot0[k][l] * Connect[i].cEdge[j].frame[0][l];
					}
				}
				double seki = 0.0;
				for (int k = 0; k < 3; k++){
					seki += Rd1[k] * Rods[rN].Mframe[eN].d3[k];
				}
				double leng = 0.0;
				for (int k = 0; k < 3; k++){
					Rd1[k] -= seki*Rods[rN].Mframe[eN].d3[k];
					leng += Rd1[k] * Rd1[k];
				}
				leng = sqrt(leng);
				if (leng >1.0e-7){
					for (int k = 0; k < 3; k++){
						Rd1[k] /= leng;
					}
					double theta = 0.0;
					double nd[3];
					for (int k = 0; k < 3; k++){
						theta += Rd1[k] * Rods[rN].Mframe[eN].d1[k];
						nd[k] = Rd1[(k + 1) % 3] * Rods[rN].Mframe[eN].d1[(k + 2) % 3]
							- Rd1[(k + 2) % 3] * Rods[rN].Mframe[eN].d1[(k + 1) % 3];
					}
					seki = 0.0;
					for (int k = 0; k < 3; k++){
						seki += nd[k] * Rods[rN].Mframe[eN].d3[k];
					}
					if (theta > 1.0) theta = 1.0;
					if (theta < -1.0) theta = -1.0;
					if (seki > 0) Cf = acos(theta);
					else Cf = -acos(theta);
				}

				varCf += Cf;
				Econ += Cf*Lamd(vN - cvlist[0] + j);
				bvec(vN + j) -= Cf;
				bvec(cvn) -= 1.0*Lamd(vN - cvlist[0] + j);
				for (int k = 0; k < 3; k++){
					bvec(snodeN + 6 * i + 3 + k) -= -Connect[i].cEdge[j].frame[2][k] * Lamd(vN - cvlist[0] + j);
				}

				kmlist.push_back(Triplet<double>(cvn, vN + j, 1.0));
				kmlist.push_back(Triplet<double>(vN + j, cvn, 1.0));
				for (int k = 0; k < 3; k++){
					kmlist.push_back(Triplet<double>(snodeN + 6 * i + 3 + k, vN + j, -Connect[i].cEdge[j].frame[2][k]));
					kmlist.push_back(Triplet<double>(vN + j, snodeN + 6 * i + 3 + k, -Connect[i].cEdge[j].frame[2][k]));
				}
			}
		}

		double c1[3];
		for (int i = 0; i < 3; i++){
			c1[i] = -Gpos[i] * (double)RnodeN;
		}
		for (int i = 0; i < rodN; i++){
			for (int j = 0; j < Rods[i].nodeN; j++){
				int vN = Rods[i].nvlist[j];
				for (int k = 0; k < 3; k++){
					kmlist.push_back(Triplet<double>(vN + k, svarN - 6 + k, 1.0));
					kmlist.push_back(Triplet<double>(svarN - 6 + k, vN + k, 1.0));
					c1[k] += Rods[i].Nodes[j][k];
					for (int l = 0; l < 3; l++){
						double val = Rods[i].Nodes[j][3 - k - l] - Gpos[3 - k - l];
						if (k != l && l == (k + 1) % 3){
							kmlist.push_back(Triplet<double>(vN + k, svarN - 3 + l, -val));
							kmlist.push_back(Triplet<double>(svarN - 3 + k, vN + l, val));
						}
						else if (k != l){
							kmlist.push_back(Triplet<double>(vN + k, svarN - 3 + l, val));
							kmlist.push_back(Triplet<double>(svarN - 3 + k, vN + l, -val));
						}
					}
				}
			}
		}
		for (int i = 0; i < 3; i++){
			bvec(svarN - 6 + i) -= c1[i];
		}


		varCp /= (float)cpN;
		varCf /= (float)cfN;

		//double wei = 1.0e-7;
		double wei = Ldamp;
		for (int i = snodeN + 6 * conN; i < svarN - 6; i++){
			kmlist.push_back(Triplet<double>(i, i, wei));
			bvec(i) -= wei*Lamd(i - cvlist[0]);
		}

		Kmat.setFromTriplets(kmlist.begin(), kmlist.end());


		//端点拘束
		/*for (int i = 0; i < 1; i++){
			int vN1 = Rods[i].nvlist[0];
			for (int j = 0; j < 4; j++){
				for (SparseMatrix<double>::InnerIterator it(Kmat, vN1 + j); it; ++it){
					if (it.row() != vN1 + j) it.valueRef() = 0.0;
					else it.valueRef() = 1.0;
				}
				bvec(vN1 + j) = 0.0;
			}

			vN1 = Rods[i].nvlist[1];
			for (int j = 0; j < 3; j++){
				for (SparseMatrix<double>::InnerIterator it(Kmat, vN1 + j); it; ++it){
					if (it.row() != vN1 + j) it.valueRef() = 0.0;
					else it.valueRef() = 1.0;
				}
				bvec(vN1 + j) = 0.0;
			}
		}
		for (int i = 0; i < svarN; i++){
			for (SparseMatrix<double>::InnerIterator it(Kmat, i); it; ++it){
				for (int j = 0; j < 1; j++){
					int vN1 = Rods[j].nvlist[0];
					if (it.row() >= vN1 && it.row() < vN1 + 4 && it.row() != i) it.valueRef() = 0.0;

					vN1 = Rods[j].nvlist[1];
					if (it.row() >= vN1 && it.row() < vN1 + 3 && it.row() != i) it.valueRef() = 0.0;
				}
			}
		}*/


		//行列計算
		SparseLU<SparseMatrix<double>> sLU;
		sLU.compute(Kmat);
		bool csuccess = (sLU.info() == Success);
		if (csuccess){
			sol = sLU.solve(bvec);
			Eigen::VectorXd dif = Kmat*sol - bvec;
			if (dif.norm() > 1.0e-4)csuccess = false;
		}
		if (!csuccess){
			bool CSuccess = false;
			double addmass = 1.0;
			printf("sparse LU failed\n");
			while (1){
				double addDamp = Ldamp*9.0;
				if (Ldamp < 1.0e-12) {
					Ldamp = 1.0e-12;
					addDamp = 1.0e-11;
				}
				Ldamp *= 10.0;
				for (int i = snodeN + 6 * conN; i < svarN - 6; i++){
					bool diag_exist = false;
					for (SparseMatrix<double>::InnerIterator it(Kmat, i); it; ++it){
						if (it.row() == i){
							it.valueRef() += addDamp;
							diag_exist = true;
						}
					}
					if (!diag_exist)Kmat.insert(i, i) = Ldamp;
					bvec(i) -= addDamp*Lamd(i - cvlist[0]);
				}
				sLU.compute(Kmat);
				CSuccess = (sLU.info() == Success);
				if (CSuccess) {
					sol = sLU.solve(bvec);
					Eigen::VectorXd dif = Kmat*sol - bvec;
					if (dif.norm() < 1.0e-4) break;
					CSuccess = false;
				}
				if (LamSize > 0 && Ldamp*LamSize > sParam->Wid*0.001) break;
				else if (LamSize < 0 && Ldamp > 1.0e-6) break;
			}

			while (1){
				for (int i = 0; i < snodeN + 6 * conN; i++){
					for (SparseMatrix<double>::InnerIterator it(Kmat, i); it; ++it){
						if (i == it.row()){
							bool is_constrained = false;
							/*for (int j = 0; j < 1; j++){
								int vN1 = Rods[j].nvlist[0];
								if (it.row() >= vN1 && it.row() < vN1 + 4) is_constrained = true;
								vN1 = Rods[j].nvlist[1];
								if (it.row() >= vN1 && it.row() < vN1 + 3) is_constrained = true;
							}*/
							if (!is_constrained){
								it.valueRef() += addmass;
							}
						}
					}
				}
				sLU.compute(Kmat);
				CSuccess = (sLU.info() == Success);
				if (CSuccess) {
					sol = sLU.solve(bvec);
					Eigen::VectorXd dif = Kmat*sol - bvec;
					if (dif.norm() < 1.0e-4) break;
					CSuccess = false;
				}
				addmass *= 2.0;
				if (addmass > 1.0e5) break;
			}
			printf("Damp%e add:%f\n",Ldamp, addmass);
			//if (addmass > 1.0e3) tstep_is_large = true;
			if (!CSuccess){
				printf("Simulation Diverge\n");
				is_converged = false;
				break;
			}
		}


		//ステップサイズの決定
		double alph = LineSearch(sol);
		if (alph < 1.0e-5){
			Eigen::VectorXd npk;
			npk = Eigen::VectorXd::Zero(svarN);
			discentDirection(npk);
			for (int i = 0; i < rodN; i++){
				for (int j = 0; j < Rods[i].nodeN; j++){
					int vN = Rods[i].nvlist[j];
					for (int k = 0; k < 3; k++){
						npk(vN + k) -= sol(svarN - 6 + k);
						for (int l = 0; l < 3; l++){
							if (k != l && l == (k + 1) % 3){
								npk(vN + k) -= -(Rods[i].Nodes[j][3 - k - l] - Gpos[3 - k - l])*sol(svarN - 3 + l);
							}
							else if (k != l){
								npk(vN + k) -= (Rods[i].Nodes[j][3 - k - l] - Gpos[3 - k - l])*sol(svarN - 3 + l);
							}
						}
					}
				}
			}

			double pkh = 0.0;
			double npkh = 0.0;
			for (int i = 0; i < svarN; i++){
				pkh += sol(i)*sol(i);
				npkh += npk(i)*npk(i);
			}
			pkh = sqrt(pkh);
			npkh = sqrt(npkh);
			if (npkh < 1.0e-8) npkh = 1.0e-8;
			for (int i = 0; i < svarN - 6; i++){
				sol(i) = npk(i) / npkh*pkh;
			}
			alph = LineSearch(sol);
		}



		//更新
		double Nmove_max = 0.0;
		double Erot_max = 0.0;
		for (int i = 0; i < rodN; i++){
			for (int j = 0; j < Rods[i].nodeN; j++){
				int vN = Rods[i].nvlist[j];
				double dif[3];
				double move = 0.0;
				for (int k = 0; k < 3; k++){
					dif[k] = Rods[i].Nodes[j][k] - Gpos[k];
					resF(vN + k) -= sol(svarN - 6 + k);
					
					Rods[i].Nodes[j][k] += alph*sol(vN + k);
					move += alph*alph*sol(vN + k)*sol(vN + k);
				}
				move = sqrt(move);
				if (move > Nmove_max) Nmove_max = move;

				if (j == Rods[i].nodeN - 1) continue;
				Rods[i].twis(j) += alph*sol(vN + 3);
				if (Erot_max < fabs(alph*sol(vN + 3))) Erot_max = fabs(alph*sol(vN + 3));

				for (int k = 0; k < 3; k++){
					for (int l = 0; l < 3; l++){
						if (k != l && l == (k + 1) % 3){
							resF(vN + k) -= -sol(svarN - 3 + l)*dif[3 - k - l];
						}
						else if (k != l){
							resF(vN + k) -= sol(svarN - 3 + l)*dif[3 - k - l];
						}
					}
				}
			}
		}
		double Rmove_max = 0.0;
		double Rrot_max = 0.0;
		double lamV = 0.0;
		for (int i = 0; i < conN; i++){
			double move = 0.0;
			for (int j = 0; j < 3; j++){
				Connect[i].gp[j] += alph*sol(snodeN + 6 * i + j);
				move += alph*alph*sol(snodeN + 6 * i + j)*sol(snodeN + 6 * i + j);
			}
			move = sqrt(move);
			if (move > Rmove_max) Rmove_max = move;

			double theta = 0.0;
			double avel[3];
			for (int j = 0; j < 3; j++){
				Wc(3 * i + j) = alph*sol(snodeN + 6 * i + 3 + j);
				avel[j] = alph*sol(snodeN + 6 * i + 3 + j);
				theta += avel[j] * avel[j];
			}
			theta = sqrt(theta);
			if (theta > 1.0e-5){
				for (int j = 0; j < 3; j++){
					avel[j] /= theta;
				}
			}
			else {
				avel[0] = 1.0; avel[1] = 0.0; avel[2] = 0.0;
			}
			if (Rrot_max < fabs(theta)) Rrot_max = fabs(theta);

			double quat[4];
			quat[0] = cos(theta / 2.0);
			for (int j = 0; j < 3; j++){
				quat[1 + j] = sin(theta / 2.0)*avel[j];
			}
			double Rot[3][3], Rot0[3][3], nRot[3][3];
			derive_rot_mat_from_quaternion(quat, Rot);
			derive_rot_mat_from_quaternion(Connect[i].quat, Rot0);
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 3; k++){
					nRot[j][k] = 0.0;
					for (int l = 0; l < 3; l++){
						nRot[j][k] += Rot0[j][l] * Rot[l][k];
					}
				}
			}
			derive_quaternion_from_rot_mat(nRot, Connect[i].quat);

			//resF用
			int vN = cvlist[i];
			for (int j = 0; j < Connect[i].nodeN; j++){
				for (int k = 0; k < 3; k++){
					//Lamd(vN - cvlist[0] + 3 * j + k) += alph*sol(vN + 3 * j + k);
					Lamd(vN - cvlist[0] + 3 * j + k) += sol(vN + 3 * j + k);
				}

				int rN = Connect[i].cNode[j].rN;
				int nN = Connect[i].cNode[j].nN;
				for (int k = 0; k < 3; k++){
					//lamV += pow(alph*sol(vN + 3 * j + k), 2);
					lamV += pow(sol(vN + 3 * j + k), 2);
					resF(Rods[rN].nvlist[nN] + k) -= preSol(vN + 3 * j + k - snodeN - 6 * conN);
					preSol(vN + 3 * j + k - snodeN - 6 * conN) = Lamd(vN - cvlist[0] + 3 * j + k);
				}
			}
			vN += 3 * Connect[i].nodeN;
			for (int j = 0; j < Connect[i].edgeN; j++){
				//Lamd(vN - cvlist[0] + j) += alph*sol(vN + j);
				Lamd(vN - cvlist[0] + j) += sol(vN + j);

				int rN = Connect[i].cEdge[j].rN;
				int eN = Connect[i].cEdge[j].eN;
				//lamV += pow(alph*sol(vN + j), 2);
				lamV += pow(sol(vN + j), 2);
				resF(Rods[rN].evlist[eN]) -= preSol(vN + j - snodeN - 6 * conN);
				preSol(vN + j - cvlist[0]) = Lamd(vN - cvlist[0] + j);
				//resF(Rods[rN].evlist[eN]) = 0.0;
			}
		}

		printf("iter%d Energy:%f (%f %f %f %f) \n", iterN,Estre + Ebend + Etwis + Econ, Estre, Ebend, Etwis,Econ);
		printf("Move: %f %f %f %f\n", Nmove_max, Erot_max*180.0 / M_PI, Rmove_max, Rrot_max*180.0 / M_PI);
		cout << "F res:" << resF.norm() << " " << preSol.norm() << " (" << sqrt(lamV) << "," << LamSize << ")" << endl;
		printf("cons %f %f damp%e\n\n", varCp, varCf,Ldamp);
		
		if (Nmove_max < sParam->Wid*0.01 && Erot_max < M_PI*0.1 / 180.0 && Rmove_max < sParam->Wid*0.01 && Rrot_max < M_PI*0.1 / 180.0){
			is_converged = true;
			break;
		}
		if (old_E > 0.0 && Estre + Ebend + Etwis > initE * 100){
			is_converged = false;
			break;
		}

		for (int i = 0; i < rodN; i++){
			Rods[i].update_frame();
		}

		LamSize = -1.0;
		for (int i = 0; i < slamdN - 6; i++){
			if (fabs(Lamd(i)) > LamSize) LamSize = fabs(Lamd(i));
		}
		if (Ldamp *LamSize > sParam->Wid * 0.001) Ldamp = sParam->Wid*0.001 / LamSize;

		iterN++;
		if (iterN == maxIter) break;
	}

	return is_converged;

}


double RodNetwork::LineSearch(Eigen::VectorXd &pk){

	double Einit;
	double Phik, Phik1, dPhik, dPhik1;
	double Phi0;
	double dPhi0;
	evaluateRodNetEnergy(0.0, pk, Phi0, dPhi0);
	double tE, tdE;
	if (dPhi0 > 0.0) {
		//printf("delX is not discent direction %f %f\n", Phi0, dPhi0);
		Eigen::VectorXd npk;
		npk = Eigen::VectorXd::Zero(svarN);
		discentDirection(npk);
		double pkh = 0.0;
		double npkh = 0.0;
		for (int i = 0; i < svarN; i++){
			pkh += pk(i)*pk(i);
			npkh += npk(i)*npk(i);
		}
		pkh = sqrt(pkh);
		npkh = sqrt(npkh);
		if (npkh < 1.0e-8) npkh = 1.0e-8;
		for (int i = 0; i < svarN; i++){
			pk(i) = npk(i) / npkh*pkh;
		}
		int vN1 = Rods[0].nvlist[0];
		int vN2 = Rods[0].nvlist[1];
		for (int i = 0; i < 4; i++){
			pk(vN1 + i) = 0.0;
			if (i != 3)pk(vN2 + i) = 0.0;
		}
		evaluateRodNetEnergy(0.0, pk, Phi0, dPhi0);
	}
	//printf("  init %f %f\n", Phi0, dPhi0);

	double al_max = 5.0;
	double alph = 1.0;
	double alp1 = 0.0;
	double c1 = 1.0e-4;
	double c2 = 0.9;
	Phik1 = Phi0;
	dPhik1 = dPhi0;
	int iterN = 0;
	while (1){
		evaluateRodNetEnergy(alph, pk, Phik, dPhik);

		bool scheck = true;
		if (Phik > Phi0 + c1*alph*dPhi0 || Phik > Phik1){
			double al_lo, al_hi, alpj;
			double Ph_lo, Ph_hi, dPh_lo, dPh_hi;
			al_lo = alp1;
			al_hi = alph;
			Ph_lo = Phik1;
			Ph_hi = Phik;
			dPh_lo = dPhik1;
			dPh_hi = dPhik;

			int innerN = 0;
			bool do_success = false;
			double Phij, dPhij;
			while (1){
				double d1 = dPh_lo + dPh_hi - 3.0*(Ph_lo - Ph_hi) / (al_lo - al_hi);
				double d2 = sqrt(d1*d1 - dPh_lo*dPh_hi);
				if (al_lo > al_hi) d2 *= -1.0;
				alpj = al_hi - (al_hi - al_lo)*(dPh_hi + d2 - d1) / (dPh_hi - dPh_lo + 2.0*d2);
				if ((al_lo - alpj)*(al_hi - alpj) > 0 || fabs(al_lo - alpj) < 0.01*fabs(al_lo - al_hi)
					|| fabs(al_hi - alpj) < 0.01*fabs(al_lo - al_hi) || __isnan(alpj)){
					alpj = (al_lo + al_hi) / 2.0;
				}
				//printf(" alph:%f %f %f %f %f\n", alpj, d1, d2, al_lo, al_hi);

				evaluateRodNetEnergy(alpj, pk, Phij, dPhij);
				//printf(" zoom1 %f(%f) %f(%f)-%f(%f)\n", alpj, Phij, al_lo, Ph_lo, al_hi, Ph_hi);
				//printf("  %d %f %f %f\n", innerN, alpj, Phij, dPhij);
				if (Phij > Phi0 + c1*alpj*dPhi0 || Phij >= Ph_lo){
					al_hi = alpj;
					Ph_hi = Phij;
					dPh_hi = dPhij;
				}
				else {
					if (fabs(dPhij) <= -c2*dPhi0){
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
					alpj = (al_lo + al_hi) / 2.0;
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
			double al_lo, al_hi, alpj;
			double Ph_lo, Ph_hi, dPh_lo, dPh_hi;
			al_hi = alp1;
			al_lo = alph;
			Ph_hi = Phik1;
			Ph_lo = Phik;
			dPh_hi = dPhik1;
			dPh_lo = dPhik;

			int innerN = 0;
			bool do_success = false;
			double Phij, dPhij;
			while (1){
				double d1 = dPh_lo + dPh_hi - 3.0*(Ph_lo - Ph_hi) / (al_lo - al_hi);
				double d2 = sqrt(d1*d1 - dPh_lo*dPh_hi);
				if (al_lo > al_hi) d2 *= -1.0;
				alpj = al_hi - (al_hi - al_lo)*(dPh_hi + d2 - d1) / (dPh_hi - dPh_lo + 2.0*d2);
				if ((al_lo - alpj)*(al_hi - alpj) > 0 || fabs(al_lo - alpj) < 0.01*fabs(al_lo - al_hi)
					|| fabs(al_hi - alpj) < 0.01*fabs(al_lo - al_hi)){
					alpj = (al_lo + al_hi) / 2.0;
				}


				evaluateRodNetEnergy(alpj, pk, Phij, dPhij);
				//printf(" zoom2 %f(%f) %f(%f)-%f(%f)\n", alpj, Phij, al_lo, Ph_lo, al_hi, Ph_hi);
				if (Phij > Phi0 + c1*alpj*dPhi0 || Phij >= Ph_lo){
					al_hi = alpj;
					Ph_hi = Phij;
					dPh_hi = dPhij;
				}
				else {
					if (fabs(dPhij) <= -c2*dPhi0){
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
					alpj = (al_lo + al_hi) / 2.0;
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

		if (1.2*alph >= al_max || al_max < 1.0e-5){
			//printf("not Found\n");
			break;
		}
		alp1 = alph;
		alph *= 1.2;

		Phik1 = Phik;
		dPhik1 = dPhik;
		iterN++;
	}

	double Phif, dPhif;
	evaluateRodNetEnergy(alph, pk, Phif, dPhif);
	//printf("  final%f %f(%f)\n", Phif, dPhif, alph);

	if (Phif > Phi0) alph = 0.0;
	return alph;

}


void RodNetwork::evaluateRodNetEnergy(double alph, Eigen::VectorXd pk, double &Eval, double &dE){

	Eigen::VectorXd lamk(slamdN);
	double **Xk, **tk, (*cXk)[3], (*cQk)[4];
	struct AdaptiveFrame **mfk, **tpk;
	Xk = new double*[rodN];
	tk = new double*[rodN];
	mfk = new struct AdaptiveFrame*[rodN];
	tpk = new struct AdaptiveFrame*[rodN];
	for (int i = 0; i < rodN; i++){
		Xk[i] = new double[3*Rods[i].nodeN];
		tk[i] = new double[Rods[i].nodeN - 1];
		mfk[i] = new struct AdaptiveFrame[Rods[i].nodeN - 1];
		tpk[i] = new struct AdaptiveFrame[Rods[i].nodeN - 1];
		for (int j = 0; j < Rods[i].nodeN; j++){
			int vN = Rods[i].nvlist[j];
			for (int k = 0; k < 3; k++){
				Xk[i][3 * j + k] = Rods[i].Nodes[j][k] + alph*pk(vN + k);
			}
			if (j == Rods[i].nodeN - 1) continue;
			tk[i][j] = Rods[i].twis(j) + alph*pk(vN + 3);
		}
	}
	cXk = new double[conN][3];
	cQk = new double[conN][4];
	for (int i = 0; i < conN; i++){
		double axis[3];
		double theta = 0.0;
		for (int j = 0; j < 3; j++){
			cXk[i][j] = Connect[i].gp[j] + alph*pk(snodeN + 6 * i + j);
			axis[j] = alph*pk(snodeN + 6 * i + 3 + j);
			theta += axis[j] * axis[j];
		}
		theta = sqrt(theta);
		if (theta > 1.0e-10){
			for (int j = 0; j < 3; j++){
				axis[j] /= theta;
			}
		}
		else {
			axis[0] = 1.0; axis[1] = 0.0; axis[2] = 0.0;
		}
		double quat[4], Rot[3][3], Rot0[3][3], nRot[3][3];
		quat[0] = cos(theta / 2.0);
		for (int j = 0; j < 3; j++){
			quat[j + 1] = sin(theta / 2.0)*axis[j];
		}
		derive_rot_mat_from_quaternion(quat, Rot);
		derive_rot_mat_from_quaternion(Connect[i].quat, Rot0);

		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				nRot[j][k] = 0.0;
				for (int l = 0; l < 3; l++){
					nRot[j][k] += Rot0[j][l] * Rot[l][k];
				}
			}
		}
		derive_quaternion_from_rot_mat(nRot, quat);
		for (int j = 0; j < 4; j++){
			cQk[i][j] = quat[j];
		}
	}

	for (int i = snodeN + 6 * conN; i < svarN; i++){
		//lamk(i - cvlist[0]) = Lamd(i - cvlist[0]) + alph*pk(i);
		lamk(i - cvlist[0]) = Lamd(i - cvlist[0]) + pk(i);
	}


	//剛性パラメータ
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

	for (int i = 0; i < rodN; i++){
		//フレームの更新
		for (int j = 0; j < Rods[i].nodeN - 1; j++){
			double d3[3];
			double dh3 = 0.0;
			for (int k = 0; k < 3; k++){
				d3[k] = Xk[i][3 * j + 3 + k] - Xk[i][3 * j + k];
				dh3 += d3[k] * d3[k];
			}
			dh3 = sqrt(dh3);
			double axis[3];
			double seki = 0.0;
			double axh = 0.0;
			for (int k = 0; k < 3; k++){
				axis[k] = Rods[i].TPframe[j].d3[(k + 1) % 3] * d3[(k + 2) % 3]
					- Rods[i].TPframe[j].d3[(k + 2) % 3] * d3[(k + 1) % 3];
				axh += axis[k] * axis[k];
				seki += Rods[i].TPframe[j].d3[k] * d3[k];
			}
			axh = sqrt(axh);
			seki /= dh3;
			if (seki > 1.0)seki = 1.0;
			if (seki < -1.0)seki = -1.0;
			seki = acos(seki);

			double quat[4], Rot[3][3];
			if (axh > 1.0e-10){
				quat[0] = cos(seki / 2.0);
				for (int k = 0; k < 3; k++){
					quat[1 + k] = sin(seki / 2.0)*axis[k] / axh;
				}
			}
			else {
				quat[0] = 1.0;
				for (int k = 0; k < 3; k++){
					quat[1 + k] = 0.0;
				}
			}
			derive_rot_mat_from_quaternion(quat, Rot);

			for (int k = 0; k < 3; k++){
				tpk[i][j].d1[k] = 0.0;
				tpk[i][j].d2[k] = 0.0;
				tpk[i][j].d3[k] = d3[k] / dh3;
				for (int l = 0; l < 3; l++){
					tpk[i][j].d1[k] += Rot[k][l] * Rods[i].TPframe[j].d1[l];
					tpk[i][j].d2[k] += Rot[k][l] * Rods[i].TPframe[j].d2[l];
				}
			}
			quat[0] = cos(tk[i][j] / 2.0);
			for (int k = 0; k < 3; k++){
				quat[1 + k] = sin(tk[i][j] / 2.0)*tpk[i][j].d3[k];
			}
			derive_rot_mat_from_quaternion(quat, Rot);
			for (int k = 0; k < 3; k++){
				mfk[i][j].d1[k] = 0.0;
				mfk[i][j].d2[k] = 0.0;
				mfk[i][j].d3[k] = tpk[i][j].d3[k];
				for (int l = 0; l < 3; l++){
					mfk[i][j].d1[k] += Rot[k][l] * tpk[i][j].d1[l];
					mfk[i][j].d2[k] += Rot[k][l] * tpk[i][j].d2[l];
				}
			}
		}

		//ストレッチ成分
		for (int j = 0; j < Rods[i].nodeN - 1; j++){
			int vN = Rods[i].nvlist[j];
			double leng = 0.0;
			double leng0 = Rods[i].oLeng[j];
			for (int k = 0; k < 3; k++){
				leng += pow(Xk[i][3 * j + 3 + k] - Xk[i][3 * j + k], 2);
			}
			leng = sqrt(leng);

			Eval += 0.5*ks*pow(leng / leng0 - 1.0, 2)*leng0;
			for (int k = 0; k < 3; k++){
				dE += -ks*(1.0 / leng0 - 1.0 / leng)*(Xk[i][3 * j + 3 + k] - Xk[i][3 * j + k])*pk(vN + k);
				dE += ks*(1.0 / leng0 - 1.0 / leng)*(Xk[i][3 * j + 3 + k] - Xk[i][3 * j + k])*pk(vN + 4 + k);
			}
		}

		//ねじり・曲げ成分
		for (int j = 1; j < Rods[i].nodeN - 1; j++){
			int vN = Rods[i].nvlist[j];
			double mbar = twist_angle_calculation(tpk[i][j - 1], tpk[i][j]);
			double tmk = tk[i][j] - tk[i][j - 1] + mbar;
			double tm0 = Rods[i].tm[j - 1];
			double kb[3], ki[2];
			double kai = 1.0;
			for (int k = 0; k < 3; k++){
				kb[k] = mfk[i][j - 1].d3[(k + 1) % 3] * mfk[i][j].d3[(k + 2) % 3] - mfk[i][j - 1].d3[(k + 2) % 3] * mfk[i][j].d3[(k + 1) % 3];
				kai += mfk[i][j - 1].d3[k] * mfk[i][j].d3[k];
			}
			ki[0] = 0.0;
			ki[1] = 0.0;
			for (int k = 0; k < 3; k++){
				kb[k] = 2.0*kb[k] / kai;
				ki[0] += 0.5*(mfk[i][j - 1].d2[k] + mfk[i][j].d2[k])*kb[k];
				ki[1] -= 0.5*(mfk[i][j - 1].d1[k] + mfk[i][j].d1[k])*kb[k];
			}

			double leng1 = 0.0;
			double leng2 = 0.0;
			for (int k = 0; k < 3; k++){
				leng1 += pow(Xk[i][3 * j + k] - Xk[i][3 * j - 3 + k], 2);
				leng2 += pow(Xk[i][3 * j + 3 + k] - Xk[i][3 * j + k], 2);
			}
			leng1 = sqrt(leng1);
			leng2 = sqrt(leng2);
			double leni = (leng1 + leng2) / 2.0;

			Eval += 0.5*kt*(tmk - tm0)*(tmk - tm0) / leni;
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 2; l++){
					Eval += 0.5*(ki[k] - Rods[i].bk[j - 1][k])*kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni;
				}
			}

			double grme[2][3], hesm[4][3][3];
			double grke[2][2][3], heske[2][4][3][3];
			double grkt[2][2], hket[2][2][2][3], hktt[2][2];
			Rods[i].derivative_of_twist_bend(j, leng1, leng2, grme, hesm, grke, grkt, heske, hket, hktt);

			for (int k = 0; k < 3; k++){
				dE += -kt*(tmk - tm0) / leni*grme[0][k] * pk(vN - 4 + k);
				dE += -kt*(tmk - tm0) / leni*(-grme[0][k] + grme[1][k])*pk(vN + k);
				dE += kt*(tmk - tm0) / leni*grme[1][k] * pk(vN + 4 + k);
				for (int l = 0; l < 2; l++){
					for (int m = 0; m < 2; m++){
						dE += -grke[l][0][k] * kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni*pk(vN - 4 + k);
						dE += -(-grke[l][0][k] + grke[l][1][k])*kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni*pk(vN + k);
						dE += grke[l][1][k] * kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni*pk(vN + 4 + k);
					}
				}
			}
			dE += -kt*(tmk - tm0) / leni*pk(vN - 1);
			dE += kt*(tmk - tm0) / leni*pk(vN + 3);
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 2; l++){
					dE += grkt[k][0] * kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni *pk(vN - 1);
					dE += grkt[k][1] * kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni*pk(vN + 3);
				}
			}

		}
	}

	//拘束関連
	for (int i = 0; i < conN; i++){
		int vN = cvlist[i];
		for (int j = 0; j < Connect[i].nodeN; j++){
			int rN = Connect[i].cNode[j].rN;
			int nN = Connect[i].cNode[j].nN;
			int cvn = Rods[rN].nvlist[nN];
			double Cp[3], crossX[3][3];
			for (int k = 0; k < 3; k++){
				Cp[k] = Xk[rN][3 * nN + k] - cXk[i][k];
				for (int l = 0; l < 3; l++){
					if (k != l && l == (k + 1) % 3) crossX[k][l] = -Connect[i].cNode[j].pos[3 - k - l];
					else if (k != l) crossX[k][l] = Connect[i].cNode[j].pos[3 - k - l];
					else crossX[k][l] = 0.0;
				}
			}

			Eigen::MatrixXf grC(3, 3);
			double Rot0[3][3];
			derive_rot_mat_from_quaternion(cQk[i], Rot0);
			for (int k = 0; k < 3; k++){
				for (int l = 0; l < 3; l++){
					Cp[k] -= Rot0[k][l] * Connect[i].cNode[j].pos[l];
					grC(k, l) = 0.0;
					for (int m = 0; m < 3; m++){
						grC(k, l) -= Rot0[k][m] * crossX[l][m];
					}
				}
			}
			for (int k = 0; k < 3; k++){
				Eval += Cp[k] * lamk(vN - cvlist[0] + 3 * j + k);
				dE += lamk(vN - cvlist[0] + 3 * j + k)*1.0*pk(cvn + k);
				dE += -lamk(vN - cvlist[0] + 3 * j + k)*1.0*pk(snodeN + 6 * i + k);
				//dE += Cp[k] * pk(vN + 3 * j + k);
				for (int l = 0; l < 3; l++){
					dE += lamk(vN - cvlist[0] + 3 * j + l)*grC(l, k)*pk(snodeN + 6 * i + 3 + k);
				}
			}
		}
		vN += 3 * Connect[i].nodeN;
		for (int j = 0; j < Connect[i].edgeN; j++){
			int rN = Connect[i].cEdge[j].rN;
			int eN = Connect[i].cEdge[j].eN;
			int cvn = Rods[rN].evlist[eN];
			double Cf = 0.0;

			double Rot0[3][3], Rd1[3];
			derive_rot_mat_from_quaternion(cQk[i], Rot0);

			for (int k = 0; k < 3; k++){
				Rd1[k] = 0.0;
				for (int l = 0; l < 3; l++){
					Rd1[k] += Rot0[k][l] * Connect[i].cEdge[j].frame[0][l];
				}
			}
			double seki = 0.0;
			for (int k = 0; k < 3; k++){
				seki += Rd1[k] * mfk[rN][eN].d3[k];
			}
			double leng = 0.0;
			for (int k = 0; k < 3; k++){
				Rd1[k] -= seki*mfk[rN][eN].d3[k];
				leng += Rd1[k] * Rd1[k];
			}
			leng = sqrt(leng);
			if (leng >1.0e-7){
				for (int k = 0; k < 3; k++){
					Rd1[k] /= leng;
				}
				double theta = 0.0;
				double nd[3];
				for (int k = 0; k < 3; k++){
					theta += Rd1[k] * mfk[rN][eN].d1[k];
					nd[k] = Rd1[(k + 1) % 3] * mfk[rN][eN].d1[(k + 2) % 3]
						- Rd1[(k + 2) % 3] * mfk[rN][eN].d1[(k + 1) % 3];
				}
				seki = 0.0;
				for (int k = 0; k < 3; k++){
					seki += nd[k] * mfk[rN][eN].d3[k];
				}
				if (theta > 1.0) theta = 1.0;
				if (theta < -1.0) theta = -1.0;
				if (seki > 0) Cf = acos(theta);
				else Cf = -acos(theta);
			}
			Eval += Cf*lamk(vN - cvlist[0] + j);
			//dE += Cf*pk(vN + j);
			dE += lamk(vN - cvlist[0] + j)*1.0*pk(cvn);
			for (int k = 0; k < 3; k++){
				dE += -Connect[i].cEdge[j].frame[2][k] * lamk(vN - cvlist[0] + j)*pk(snodeN + 6 * i + 3 + k);
			}
		}
	}

	double c1[3], c2[3];
	for (int i = 0; i < 3; i++){
		c1[i] = 0.0;
		c2[i] = 0.0;
	}
	for (int i = 0; i < rodN; i++){
		for (int j = 0; j < Rods[i].nodeN; j++){
			int vN = Rods[i].nvlist[j];
			for (int k = 0; k < 3; k++){
				c1[k] += Xk[i][3 * j + k] - Gpos[k];
				dE += lamk(slamdN - 6 + k)*pk(vN + k);
				for (int l = 0; l < 3; l++){
					if (k != l && l == (k + 1) % 3){
						c2[k] += (Rods[i].Nodes[j][3 - k - l] - Gpos[3 - k - l])*alph*pk(vN + l);
						dE += lamk(slamdN - 3 + k)*(Rods[i].Nodes[j][3 - k - l] - Gpos[3 - k - l])*pk(vN + l);
					}
					else if (k != l){
						c2[k] += -(Rods[i].Nodes[j][3 - k - l] - Gpos[3 - k - l])*alph*pk(vN + l);
						dE += -lamk(slamdN - 3 + k)*(Rods[i].Nodes[j][3 - k - l] - Gpos[3 - k - l])*pk(vN + l);
					}
				}
			}
		}
	}
	for (int i = 0; i < 3; i++){
		Eval += lamk(slamdN - 6 + i)*c1[i];
		Eval += lamk(slamdN - 3 + i)*c2[i];
	}


	for (int i = 0; i < rodN; i++){
		delete[] Xk[i];
		delete[] tk[i];
		delete[] mfk[i];
		delete[] tpk[i];
	}
	delete[] Xk;
	delete[] tk;
	delete[] mfk;
	delete[] tpk;
	delete[] cXk;
	delete[] cQk;

}


void RodNetwork::discentDirection(Eigen::VectorXd &npk){


	//剛性パラメータ
	double Across = sParam->Hei*sParam->Wid;
	double ks = sParam->Es*Across;
	double kt = sParam->Gt*Across * (sParam->Wid*sParam->Wid + sParam->Hei*sParam->Hei) / 12.0;
	double kben[2][2];
	kben[0][0] = sParam->Eb11*Across*sParam->Wid*sParam->Wid / 12.0;
	kben[1][1] = sParam->Eb22 *Across*sParam->Hei*sParam->Hei / 12.0;
	kben[0][1] = sParam->Eb12*Across*sParam->Hei*sParam->Wid / 12.0;
	kben[1][0] = sParam->Eb12*Across*sParam->Wid*sParam->Hei / 12.0;

	for (int i = 0; i < svarN; i++){
		npk(i) = 0.0;
	}

	for (int i = 0; i < rodN; i++){

		//ストレッチ成分
		for (int j = 0; j < Rods[i].nodeN - 1; j++){
			int vN = Rods[i].nvlist[j];
			double leng = 0.0;
			double leng0 = Rods[i].oLeng[j];
			for (int k = 0; k < 3; k++){
				leng += pow(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k], 2);
			}
			leng = sqrt(leng);

			for (int k = 0; k < 3; k++){
				npk(vN + k) += ks*(1.0 / leng0 - 1.0 / leng)*(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k]);
				npk(vN + 4 + k) += -ks*(1.0 / leng0 - 1.0 / leng)*(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k]);
			}
		}

		//ねじり・曲げ成分
		for (int j = 1; j < Rods[i].nodeN - 1; j++){
			int vN = Rods[i].nvlist[j];
			double mbar = twist_angle_calculation(Rods[i].TPframe[j - 1], Rods[i].TPframe[j]);
			double tmk = Rods[i].twis(j) - Rods[i].twis(j - 1) + mbar;
			double tm0 = Rods[i].tm[j - 1];
			double kb[3], ki[2];
			double kai = 1.0;
			for (int k = 0; k < 3; k++){
				kb[k] = Rods[i].Mframe[j - 1].d3[(k + 1) % 3] * Rods[i].Mframe[j].d3[(k + 2) % 3]
					- Rods[i].Mframe[j - 1].d3[(k + 2) % 3] * Rods[i].Mframe[j].d3[(k + 1) % 3];
				kai += Rods[i].Mframe[j - 1].d3[k] * Rods[i].Mframe[j].d3[k];
			}
			ki[0] = 0.0;
			ki[1] = 0.0;
			for (int k = 0; k < 3; k++){
				kb[k] = 2.0*kb[k] / kai;
				ki[0] += 0.5*(Rods[i].Mframe[j - 1].d2[k] + Rods[i].Mframe[j].d2[k])*kb[k];
				ki[1] -= 0.5*(Rods[i].Mframe[j - 1].d1[k] + Rods[i].Mframe[j].d1[k])*kb[k];
			}

			double leng1 = 0.0;
			double leng2 = 0.0;
			for (int k = 0; k < 3; k++){
				leng1 += pow(Rods[i].Nodes[j][k] - Rods[i].Nodes[j - 1][k], 2);
				leng2 += pow(Rods[i].Nodes[j + 1][k] - Rods[i].Nodes[j][k], 2);
			}
			leng1 = sqrt(leng1);
			leng2 = sqrt(leng2);
			double leni = (leng1 + leng2) / 2.0;

			double grme[2][3], hesm[4][3][3];
			double grke[2][2][3], heske[2][4][3][3];
			double grkt[2][2], hket[2][2][2][3], hktt[2][2];
			Rods[i].derivative_of_twist_bend(j, leng1, leng2, grme, hesm, grke, grkt, heske, hket, hktt);

			for (int k = 0; k < 3; k++){
				npk(vN - 4 + k) += kt*(tmk - tm0) / leni*grme[0][k];
				npk(vN + k) += kt*(tmk - tm0) / leni*(-grme[0][k] + grme[1][k]);
				npk(vN + 4 + k) += -kt*(tmk - tm0) / leni*grme[1][k];
				for (int l = 0; l < 2; l++){
					for (int m = 0; m < 2; m++){
						npk(vN - 4 + k) += grke[l][0][k] * kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
						npk(vN + k) += (-grke[l][0][k] + grke[l][1][k])*kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
						npk(vN + 4 + k) += -grke[l][1][k] * kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
					}
				}
			}
			npk(vN - 1) += kt*(tmk - tm0) / leni;
			npk(vN + 3) += -kt*(tmk - tm0) / leni;
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 2; l++){
					npk(vN - 1) += -grkt[k][0] * kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni;
					npk(vN + 3) += -grkt[k][1] * kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni;
				}
			}

		}
	}

	//拘束関連
	for (int i = 0; i < conN; i++){
		int vN = cvlist[i];
		for (int j = 0; j < Connect[i].nodeN; j++){
			int rN = Connect[i].cNode[j].rN;
			int nN = Connect[i].cNode[j].nN;
			int cvn = Rods[rN].nvlist[nN];
			double Cp[3], crossX[3][3];
			for (int k = 0; k < 3; k++){
				Cp[k] = Rods[rN].Nodes[nN][k] - Connect[i].gp[k];
				for (int l = 0; l < 3; l++){
					if (k != l && l == (k + 1) % 3) crossX[k][l] = -Connect[i].cNode[j].pos[3 - k - l];
					else if (k != l) crossX[k][l] = Connect[i].cNode[j].pos[3 - k - l];
					else crossX[k][l] = 0.0;
				}
			}

			Eigen::MatrixXf grC(3, 3);
			double Rot0[3][3];
			derive_rot_mat_from_quaternion(Connect[i].quat, Rot0);
			for (int k = 0; k < 3; k++){
				for (int l = 0; l < 3; l++){
					Cp[k] -= Rot0[k][l] * Connect[i].cNode[j].pos[l];
					grC(k, l) = 0.0;
					for (int m = 0; m < 3; m++){
						grC(k, l) -= Rot0[k][m] * crossX[l][m];
					}
				}
			}

			for (int k = 0; k < 3; k++){
				npk(cvn + k) += -Lamd(vN - cvlist[0] + 3 * j + k)*1.0;
				npk(snodeN + 6 * i + k) += Lamd(vN - cvlist[0] + 3 * j + k)*1.0;
				//npk(vN + 3 * j + k) += -Cp[k];
				for (int l = 0; l < 3; l++){
					npk(snodeN + 6 * i + 3 + k) += -Lamd(vN - cvlist[0] + 3 * j + l)*grC(l, k);
				}
			}
		}
		vN += 3 * Connect[i].nodeN;
		for (int j = 0; j < Connect[i].edgeN; j++){
			int rN = Connect[i].cEdge[j].rN;
			int eN = Connect[i].cEdge[j].eN;
			int cvn = Rods[rN].evlist[eN];
			double Cf = 0.0;

			double Rot0[3][3], Rd1[3];
			derive_rot_mat_from_quaternion(Connect[i].quat, Rot0);

			for (int k = 0; k < 3; k++){
				Rd1[k] = 0.0;
				for (int l = 0; l < 3; l++){
					Rd1[k] += Rot0[k][l] * Connect[i].cEdge[j].frame[0][l];
				}
			}
			double seki = 0.0;
			for (int k = 0; k < 3; k++){
				seki += Rd1[k] * Rods[rN].Mframe[eN].d3[k];
			}
			double leng = 0.0;
			for (int k = 0; k < 3; k++){
				Rd1[k] -= seki*Rods[rN].Mframe[eN].d3[k];
				leng += Rd1[k] * Rd1[k];
			}
			leng = sqrt(leng);
			if (leng >1.0e-7){
				for (int k = 0; k < 3; k++){
					Rd1[k] /= leng;
				}
				double theta = 0.0;
				double nd[3];
				for (int k = 0; k < 3; k++){
					theta += Rd1[k] * Rods[rN].Mframe[eN].d1[k];
					nd[k] = Rd1[(k + 1) % 3] * Rods[rN].Mframe[eN].d1[(k + 2) % 3]
						- Rd1[(k + 2) % 3] * Rods[rN].Mframe[eN].d1[(k + 1) % 3];
				}
				seki = 0.0;
				for (int k = 0; k < 3; k++){
					seki += nd[k] * Rods[rN].Mframe[eN].d3[k];
				}
				if (theta > 1.0) theta = 1.0;
				if (theta < -1.0) theta = -1.0;
				if (seki > 0) Cf = acos(theta);
				else Cf = -acos(theta);
			}

			//npk(vN + j) += -Cf;
			npk(cvn) += -Lamd(vN - cvlist[0] + j);
			for (int k = 0; k < 3; k++){
				npk(snodeN + 6 * i + 3 + k) += Connect[i].cEdge[j].frame[2][k] * Lamd(vN - cvlist[0] + j);
			}
		}
	}

}


bool RodNetwork::curve_to_rod(class CurveNetwork *cNet,class AABB_struct *iaabb,class loop_subdivision_mesh *lmesh){


	if (cNet->curSet.empty() || cNet->Connect.empty()) return false;

	//元の配列の解放
	if (rodN > 0){
		delete[] Rods;
		rodN = 0;
	}
	if (conN > 0){
		delete[] Connect;
		delete[] cvlist;
	}
	

	//曲線の変換
	int **CNlist, *CNsize;
	CNlist = new int*[cNet->curSet.size()];
	CNsize = new int[cNet->curSet.size()];
	for (int i = 0; i < cNet->curSet.size(); i++){
		CNsize[i] = 0;
	}
	for (int i = 0; i < cNet->Connect.size(); i++){
		CNsize[cNet->Connect[i].c1]++;
		CNsize[cNet->Connect[i].c2]++;
	}
	for (int i = 0; i < cNet->curSet.size(); i++){
		if (CNsize[i] > 0) CNlist[i] = new int[CNsize[i]];
		CNsize[i] = 0;
	}
	for (int i = 0; i < cNet->Connect.size(); i++){
		int c1 = cNet->Connect[i].c1;
		int c2 = cNet->Connect[i].c2;
		bool exist = false;
		for (int j = 0; j < CNsize[c1]; j++){
			if (CNlist[c1][j] == c2) exist = true;
		}
		if (!exist){
			CNlist[c1][CNsize[c1]] = c2;
			CNsize[c1]++;
		}
		exist = false;
		for (int j = 0; j < CNsize[c2]; j++){
			if (CNlist[c2][j] == c1) exist = true;
		}
		if (!exist){
			CNlist[c2][CNsize[c2]] = c1;
			CNsize[c2]++;
		}
	}
	
	int *ccomp;
	ccomp = new int[cNet->curSet.size()];
	for (int i = 0; i < cNet->curSet.size(); i++){
		ccomp[i] = -1;
	}
	int cmN = 0;
	vector<int> compS;
	for (int i = 0; i < cNet->curSet.size(); i++){
		if (ccomp[i] == -1){
			traverseCNlist(i,ccomp,CNlist,CNsize,cmN);
			cmN++;
			compS.push_back(0);
		}
	}
	for (int i = 0; i < cNet->curSet.size(); i++){
		compS[ccomp[i]]++;
	}
	int maxN,maxS;
	for (int i = 0; i < cmN; i++){
		if (i == 0 || maxS < compS[i]){
			maxS = compS[i];
			maxN = i;
		}
	}

	rodN = 0;
	for (int i = 0; i < cNet->curSet.size(); i++){
		if (ccomp[i] == maxN){
			rodN++;
		}
	}
	Rods = new DiscreteRod[rodN];
	rodN = 0;
	for (int i = 0; i < cNet->curSet.size(); i++){
		if (ccomp[i] == maxN){
			Rods[rodN].nodeN = cNet->curSet[i]->nodeN;
			Rods[rodN].Nodes = new double[Rods[rodN].nodeN][3];
			Rods[rodN].TPframe = new AdaptiveFrame[Rods[rodN].nodeN - 1];
			Rods[rodN].Mframe = new AdaptiveFrame[Rods[rodN].nodeN - 1];
			for (int j = 0; j < Rods[rodN].nodeN; j++){
				for (int k = 0; k < 3; k++){
					Rods[rodN].Nodes[j][k] = (double)cNet->curSet[i]->cNode[j].pos[k];
					if (j != Rods[rodN].nodeN - 1){
						Rods[rodN].TPframe[j].d1[k] = (double)cNet->curSet[i]->TPframe[j].d1[k];
						Rods[rodN].TPframe[j].d2[k] = (double)cNet->curSet[i]->TPframe[j].d2[k];
						Rods[rodN].TPframe[j].d3[k] = (double)cNet->curSet[i]->TPframe[j].d3[k];
						Rods[rodN].Mframe[j].d1[k] = (double)cNet->curSet[i]->Mframe[j].d1[k];
						Rods[rodN].Mframe[j].d2[k] = (double)cNet->curSet[i]->Mframe[j].d2[k];
						Rods[rodN].Mframe[j].d3[k] = (double)cNet->curSet[i]->Mframe[j].d3[k];
					}
				}
			}
			Rods[rodN].sParam = sParam;
			Rods[rodN].FRAME_EXIST = true;
			Rods[rodN].NODE_EXIST = true;
			ccomp[i] = rodN;
			rodN++;
		}
		else {
			ccomp[i] = -1;
		}
	}

	for (int i = 0; i < cNet->curSet.size(); i++){
		if (CNsize[i] > 0) delete[] CNlist[i];
	}

	delete[] CNsize;
	delete[] CNlist;


	//コネクションの整理
	int *cmap;
	bool *visited;
	cmap = new int[cNet->Connect.size()];
	visited = new bool[cNet->Connect.size()];
	for (int i = 0; i < cNet->Connect.size(); i++){
		cmap[i] = -1;
		visited[i] = false;
	}

	int conNum = 0;
	for (int i = 0; i < cNet->Connect.size(); i++){
		if (!visited[i]){
			struct qelem *startq;
			startq = new struct qelem;
			startq->ln = i;
			struct qelem *endq = cNet->traverseConnectGraph(i, visited, startq);

			struct qelem *cque = endq;
			while (1){
				cmap[cque->ln] = conNum;
				if (cque == startq) break;
				cque = cque->back;
			}
			conNum++;
		}
		if (cmap[i] == -1) printf("%d connection is not traversed\n",i);
	}

	std::vector<std::vector<std::pair<int, int>>> Cgroup;
	for (int i = 0; i < conNum; i++){
		std::vector<std::pair<int, int>> Cpair;
		for (int j = 0; j < cNet->Connect.size(); j++){
			if (cmap[j] == i && ccomp[cNet->Connect[i].c1] != -1 && ccomp[cNet->Connect[i].c2] != -1){
				bool exist1 = false;
				bool exist2 = false;
				for (int k = 0; k < Cpair.size(); k++){
					if (Cpair[k].first == ccomp[cNet->Connect[j].c1] && Cpair[k].second == cNet->Connect[j].n1) exist1 = true;
					if (Cpair[k].first == ccomp[cNet->Connect[j].c2] && Cpair[k].second == cNet->Connect[j].n2) exist2 = true;
				}
				if (!exist1) Cpair.push_back(make_pair(ccomp[cNet->Connect[j].c1], cNet->Connect[j].n1));
				if (!exist2) Cpair.push_back(make_pair(ccomp[cNet->Connect[j].c2], cNet->Connect[j].n2));
			}
		}
		if(!Cpair.empty()) Cgroup.push_back(std::vector<std::pair<int,int>>(Cpair));
	}

	conN = Cgroup.size();
	Connect = new ConnectionData[conN];
	for (int i = 0; i < conN; i++){
		float gp[3];
		for (int j = 0; j < 3; j++){
			gp[j] = 0.0;
		}
		for (int j = 0; j < Cgroup[i].size(); j++){
			int rN = Cgroup[i][j].first;
			int nN = Cgroup[i][j].second;
			for (int k = 0; k < 3; k++){
				gp[k] += Rods[rN].Nodes[nN][k] / (double)(Cgroup[i].size());
			}
		}
		int fn;
		float cord[3], norm[3];
		bool success = iaabb->closest_point_search_with_mesh_info(gp, &fn, cord);
		if (success){
			float dbs[2][3];
			lmesh->evaluate_derivative_of_box_spline_on_triangle(fn, cord, dbs);
			float nh = 0.0;
			for (int k = 0; k < 3; k++){
				norm[k] = dbs[0][(k + 1) % 3] * dbs[1][(k + 2) % 3] - dbs[0][(k + 2) % 3] * dbs[1][(k + 1) % 3];
				nh += norm[k] * norm[k];
			}
			nh = sqrt(nh);
			if (nh < 1.0e-10) success = false;
			if (success){
				for (int k = 0; k < 3; k++){
					Connect[i].norm[k] = (double)norm[k] / nh;
				}
			}
		}
		if (!success){
			for (int k = 0; k < 3; k++){
				Connect[i].norm[k] = 0.0;
			}
			Connect[i].norm[0] = 1.0;
		}
		for (int k = 0; k < 3; k++){
			Connect[i].gp[k] = (double)gp[k];
		}

		std::vector<std::pair<int, int>> cnode;
		for (int j = 0; j < Cgroup[i].size(); j++){
			int cN = Cgroup[i][j].first;
			int nN = Cgroup[i][j].second;
			std::vector<std::pair<int, int>> ehozon;
			if (nN != 0) ehozon.push_back(make_pair(cN, nN - 1));
			if (nN != Rods[cN].nodeN - 1) ehozon.push_back(make_pair(cN, nN));
			for (int k = 0; k < ehozon.size(); k++){
				int cN = ehozon[k].first;
				int eN = ehozon[k].second;
				bool is_duplicate = false;
				for (int l = 0; l < cnode.size(); l++){
					if (cnode[l].first == cN && cnode[l].second == eN) is_duplicate = true;
				}
				if (!is_duplicate) cnode.push_back(ehozon[k]);
			}
		}
		Connect[i].edgeN = cnode.size();
		Connect[i].cEdge = new ConnectedEdge[cnode.size()];
		for (int j = 0; j < cnode.size(); j++){
			Connect[i].cEdge[j].rN = cnode[j].first;
			Connect[i].cEdge[j].eN = cnode[j].second;
		}
	}

	for (int i = 0; i < conN; i++){
		std::vector<std::pair<int, int>> cnode;
		for (int j = 0; j < Connect[i].edgeN; j++){
			bool N1exist = false;
			bool N2exist = false;
			for (int k = 0; k < cnode.size(); k++){
				if (cnode[k].first == Connect[i].cEdge[j].rN && cnode[k].second == Connect[i].cEdge[j].eN) N1exist = true;
				if (cnode[k].first == Connect[i].cEdge[j].rN && cnode[k].second == Connect[i].cEdge[j].eN + 1) N2exist = true;
			}
			if (!N1exist){
				cnode.push_back(make_pair(Connect[i].cEdge[j].rN, Connect[i].cEdge[j].eN));
			}
			if (!N2exist){
				cnode.push_back(make_pair(Connect[i].cEdge[j].rN, Connect[i].cEdge[j].eN + 1));
			}
			N1exist = false;
			N2exist = false;
			for (int k = 0; k < Cgroup[i].size(); k++){
				if (Connect[i].cEdge[j].rN == Cgroup[i][k].first && Connect[i].cEdge[j].eN == Cgroup[i][k].second) N1exist = true;
				if (Connect[i].cEdge[j].rN == Cgroup[i][k].first && Connect[i].cEdge[j].eN + 1 == Cgroup[i][k].second) N2exist = true;
			} 
			if (N1exist && !N2exist){
				Connect[i].cEdge[j].relate = 0;
			}
			else if (N2exist && !N1exist){
				Connect[i].cEdge[j].relate = 1;
			}
			else {
				Connect[i].cEdge[j].relate = -1;
			}

		}
		Connect[i].nodeN = cnode.size();
		Connect[i].cNode = new ConnectedNode[cnode.size()];
		for (int j = 0; j < cnode.size(); j++){
			Connect[i].cNode[j].rN = cnode[j].first;
			Connect[i].cNode[j].nN = cnode[j].second;
		}
	}

	delete[] cmap;
	delete[] visited;
	delete[] ccomp;

	return true;

}


void RodNetwork::exportRodNetwork(std::string fname,class Mesh *mesh){

	//ファイル出力
	FILE *fp = fopen(fname.c_str(), "w");
	if (fp == NULL){
		printf("File Open Error(%s)\n",fname);
		return;
	}

	fprintf(fp, "%d %d\n",rodN,conN);
	for (int i = 0; i < rodN; i++){
		fprintf(fp,"%d %d\n",i,Rods[i].nodeN);
		for (int j = 0; j < Rods[i].nodeN; j++){
			fprintf(fp,"%d %f %f %f\n",j,Rods[i].Nodes[j][0],Rods[i].Nodes[j][1],Rods[i].Nodes[j][2]);
		}
		for (int j = 0; j < Rods[i].nodeN - 1; j++){
			fprintf(fp, "%d %f %f %f", j, Rods[i].TPframe[j].d1[0], Rods[i].TPframe[j].d1[1], Rods[i].TPframe[j].d1[2]);
			fprintf(fp," %f %f %f %f %f %f\n",Rods[i].TPframe[j].d2[0],Rods[i].TPframe[j].d2[1],Rods[i].TPframe[j].d2[2]
				,Rods[i].TPframe[j].d3[0],Rods[i].TPframe[j].d3[1],Rods[i].TPframe[j].d3[2]);
			fprintf(fp, " %f %f %f", Rods[i].Mframe[j].d1[0], Rods[i].Mframe[j].d1[1], Rods[i].Mframe[j].d1[2]);
			fprintf(fp, " %f %f %f %f %f %f\n", Rods[i].Mframe[j].d2[0], Rods[i].Mframe[j].d2[1], Rods[i].Mframe[j].d2[2]
				, Rods[i].Mframe[j].d3[0], Rods[i].Mframe[j].d3[1], Rods[i].Mframe[j].d3[2]);
		}
	}
	for (int i = 0; i < conN; i++){
		fprintf(fp, "%d %d %d\n", i, Connect[i].nodeN, Connect[i].edgeN);
		fprintf(fp, " %f %f %f %f %f %f\n", Connect[i].gp[0], Connect[i].gp[1], Connect[i].gp[2], Connect[i].norm[0], Connect[i].norm[1], Connect[i].norm[2]);
		for (int j = 0; j < Connect[i].nodeN; j++){
			fprintf(fp, " %d %d", Connect[i].cNode[j].rN, Connect[i].cNode[j].nN);
		}
		fprintf(fp, "\n");
		for (int j = 0; j < Connect[i].edgeN; j++){
			fprintf(fp, " %d %d %d", Connect[i].cEdge[j].rN, Connect[i].cEdge[j].eN,Connect[i].cEdge[j].relate);
		}
		fprintf(fp,"\n");
	}

	fprintf(fp,"%d %d\n",mesh->nodeN,mesh->faceN);
	for (int i = 0; i < mesh->nodeN; i++){
		fprintf(fp,"%d %f %f %f\n",i+1,mesh->Nodes[i][0],mesh->Nodes[i][1],mesh->Nodes[i][2]);
	}
	for (int i = 0; i < mesh->faceN; i++){
		fprintf(fp,"%d %d %d %d\n",i+1,mesh->Faces[i][0],mesh->Faces[i][1],mesh->Faces[i][2]);
	}

	fclose(fp);


	//パラメータの出力
	std::string fpass = fname;
	size_t pos1;
	pos1 = fname.rfind('\\');
	if (pos1 != std::string::npos){
		fpass = fname.substr(0,pos1 + 1);
	}
	else {
		pos1 = fname.rfind('/');
		if (pos1 != std::string::npos){
			fpass = fname.substr(0, pos1 + 1);
		}
	}
	std::string pfname = fpass + "Parameter.spf";
	
	fp = fopen(pfname.c_str(),"w");
	if (fp == NULL){
		printf("File Open Error(%s)\n",pfname);
		return;
	}

	fprintf(fp, "%f %f %f\n", sParam->Wid, sParam->Hei,sParam->Density*1000.0);
	fprintf(fp, "%f %f %f %f %f\n", sParam->Es, sParam->Gt, sParam->Eb11, sParam->Eb22, sParam->Eb12);

	fclose(fp);

}
 

////////////
////その他の関数
////////////

//点群位置合わせ
void RodNetwork::PointsMatching(){

	//剛体変換の計算
	double gpl[3], gpr[3];
	for (int i = 0; i < 3; i++){
		gpl[i] = 0.0;
		gpr[i] = 0.0;
	}

	int pnum = 0;
	for (int i = 0; i < rodN; i++){
		for (int j = 0; j < Rods[i].nodeN; j++){
			for (int k = 0; k < 3; k++){
				gpl[k] += Rods[i].Nodes[j][k];
				gpr[k] += Rods[i].oNodes[j][k];
			}
			pnum++;
		}
	}

	for (int i = 0; i < 3; i++){
		gpl[i] /= (double)pnum;
		gpr[i] /= (double)pnum;
	}

	double Mmat[3][3];
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			Mmat[i][j] = 0.0;
		}
	}
	for (int i = 0; i < rodN; i++){
		for (int j = 0; j < Rods[i].nodeN; j++){
			for (int k = 0; k < 3; k++){
				for (int l = 0; l < 3; l++){
					Mmat[k][l] += (Rods[i].Nodes[j][k] - gpl[k])*(Rods[i].oNodes[j][l] - gpr[l]);
				}
			}
		}
	}

	Eigen::Matrix4d Nmat;
	Eigen::Vector4d qvec;
	Nmat(0, 0) = Mmat[0][0] + Mmat[1][1] + Mmat[2][2];
	for (int i = 0; i < 3; i++){
		Nmat(0, i + 1) = Mmat[(i + 1) % 3][(i + 2) % 3] - Mmat[(i + 2) % 3][(i + 1) % 3];
		Nmat(i + 1, 0) = Mmat[(i + 1) % 3][(i + 2) % 3] - Mmat[(i + 2) % 3][(i + 1) % 3];
	}
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			if (i == j) Nmat(i + 1, j + 1) = Mmat[i][i] - Mmat[(i + 1) % 3][(i + 2) % 3] - Mmat[(i + 2) % 3][(i + 2) % 3];
			else Nmat(i + 1, j + 1) = Mmat[i][j] + Mmat[j][i];
		}
	}

	SelfAdjointEigenSolver<Matrix4d> eigensolver(Nmat);
	if (eigensolver.info() != Success) printf("Eigen Decomp failed\n");
	qvec = eigensolver.eigenvectors().col(3);

	double quat[4], Rot[3][3];
	double qh = 0.0;
	for (int i = 0; i < 4; i++){
		quat[i] = qvec(i);
		qh += quat[i] * quat[i];
	}
	qh = sqrt(qh);
	if (qh < 1.0e-10) qh = 1.0e-10;
	for (int i = 0; i < 4; i++){
		quat[i] /= qh;
	}
	derive_rot_mat_from_quaternion(quat, Rot);

	double Tran[3];
	for (int i = 0; i < 3; i++){
		Tran[i] = gpr[i];
		for (int j = 0; j < 3; j++){
			Tran[i] -= Rot[i][j] * gpl[j];
		}
	}

	//剛体変換の適用
	for (int i = 0; i < rodN; i++){
		for (int j = 0; j < Rods[i].nodeN; j++){
			double rpos[3];
			for (int k = 0; k < 3; k++){
				rpos[k] = Tran[k];
				for (int l = 0; l < 3; l++){
					rpos[k] += Rot[k][l] * Rods[i].Nodes[j][l];
				}
			}
			for (int k = 0; k < 3; k++){
				Rods[i].Nodes[j][k] = rpos[k];
			}
			if (j != Rods[i].nodeN - 1){
				double Tframe[3][3], Mframe[3][3];
				for (int k = 0; k < 3; k++){
					Tframe[0][k] = 0.0;
					Tframe[1][k] = 0.0;
					Tframe[2][k] = 0.0;
					Mframe[0][k] = 0.0;
					Mframe[1][k] = 0.0;
					Mframe[2][k] = 0.0;
					for (int l = 0; l < 3; l++){
						Tframe[0][k] += Rot[k][l] * Rods[i].TPframe[j].d1[l];
						Tframe[1][k] += Rot[k][l] * Rods[i].TPframe[j].d2[l];
						Tframe[2][k] += Rot[k][l] * Rods[i].TPframe[j].d3[l];
						Mframe[0][k] += Rot[k][l] * Rods[i].Mframe[j].d1[l];
						Mframe[1][k] += Rot[k][l] * Rods[i].Mframe[j].d2[l];
						Mframe[2][k] += Rot[k][l] * Rods[i].Mframe[j].d3[l];
					}
				}
				for (int k = 0; k < 3; k++){
					Rods[i].TPframe[j].d1[k] = Tframe[0][k];
					Rods[i].TPframe[j].d2[k] = Tframe[1][k];
					Rods[i].TPframe[j].d3[k] = Tframe[2][k];
					Rods[i].Mframe[j].d1[k] = Mframe[0][k];
					Rods[i].Mframe[j].d2[k] = Mframe[1][k];
					Rods[i].Mframe[j].d3[k] = Mframe[2][k];
				}
			}
		}
	}

	for (int i = 0; i < conN; i++){
		double Rot0[3][3], nRot[3][3];
		double ngp[3];
		derive_rot_mat_from_quaternion(Connect[i].quat, Rot0);
		for (int j = 0; j < 3; j++){
			ngp[j] = Tran[j];
			for (int k = 0; k < 3; k++){
				nRot[j][k] = 0.0;
				ngp[j] += Rot[j][k] * Connect[i].gp[k];
				for (int l = 0; l < 3; l++){
					nRot[j][k] += Rot[j][l] * Rot0[l][k];
				}
			}
		}
		derive_quaternion_from_rot_mat(nRot, Connect[i].quat);
		for (int j = 0; j < 3; j++){
			Connect[i].gp[j] = ngp[j];
		}
	}


}


void RodNetwork::traverseCNlist(int N, int *ccomp, int**CNlist, int*CNsize, int cmN){

	ccomp[N] = cmN;

	for (int i = 0; i < CNsize[N]; i++){
		if (ccomp[CNlist[N][i]] == -1){
			traverseCNlist(CNlist[N][i], ccomp, CNlist, CNsize, cmN);
		}
	}
}


bool RodNetwork::simulateOneStep(Eigen::VectorXd &Disp,Eigen::VectorXd &resF){

	Eigen::VectorXd bvec(svarN), sol(svarN);

	int RnodeN = 0;
	for (int i = 0; i < 3; i++){
		Gpos[i] = 0.0;
	}
	for (int i = 0; i < rodN; i++){
		for (int j = 0; j < Rods[i].nodeN; j++){
			int vN = Rods[i].nvlist[j];
			for (int k = 0; k < 3; k++){
				Gpos[k] += Rods[i].oNodes[j][k];
			}
			RnodeN++;
		}
	}
	for (int i = 0; i < 3; i++){
		Gpos[i] /= (double)RnodeN;
	}

	double Ldamp = 0.0;
	double LamSize = -1.0;
	double Mu = 1.0e-3;

	//剛性マトリックス・残差ベクトルの計算
	Eigen::SparseMatrix<double> Kmat(svarN, svarN);
	std::vector<Triplet<double>> kmlist;
	bvec = Eigen::VectorXd::Zero(svarN);
	resF = Eigen::VectorXd::Zero(snodeN);


	//剛性パラメータ
	double Across = sParam->Hei*sParam->Wid;
	double ks = sParam->Es*Across;
	double kt = sParam->Gt*Across * (sParam->Wid*sParam->Wid + sParam->Hei*sParam->Hei) / 12.0;
	double kben[2][2];
	kben[0][0] = sParam->Eb11*Across*sParam->Wid*sParam->Wid / 12.0;
	kben[1][1] = sParam->Eb22 *Across*sParam->Hei*sParam->Hei / 12.0;
	kben[0][1] = sParam->Eb12*Across*sParam->Hei*sParam->Wid / 12.0;
	kben[1][0] = sParam->Eb12*Across*sParam->Wid*sParam->Hei / 12.0;

	double Estre = 0.0;
	double Ebend = 0.0;
	double Etwis = 0.0;
	double Eext = 0.0;
	double Econ = 0.0;

	double varCf = 0.0;
	double varCp = 0.0;

	for (int i = 0; i < rodN; i++){
		//ねじり計算
		Eigen::VectorXd otwis(Rods[i].nodeN - 1);
		for (int j = 0; j < Rods[i].nodeN - 1; j++){
			double thet = 0.0;
			double seki = 0.0;
			for (int k = 0; k < 3; k++){
				thet += Rods[i].oMframe[j].d1[k] * Rods[i].oTPframe[j].d1[k];
				seki += Rods[i].oMframe[j].d3[k] * (Rods[i].oTPframe[j].d1[(k + 1) % 3] * Rods[i].oMframe[j].d1[(k + 2) % 3]
					- Rods[i].oTPframe[j].d1[(k + 2) % 3] * Rods[i].oMframe[j].d1[(k + 1) % 3]);
			}
			if (thet > 1.0)thet = 1.0;
			if (thet < -1.0)thet = -1.0;
			thet = acos(thet);
			if (seki < 0.0) thet *= -1.0;
			otwis(j) = thet;
		}

		//ストレッチ成分
		for (int j = 0; j < Rods[i].nodeN - 1; j++){
			int vN = Rods[i].nvlist[j];
			double leng = 0.0;
			double leng0 = Rods[i].oLeng[j];
			for (int k = 0; k < 3; k++){
				leng += pow(Rods[i].oNodes[j + 1][k] - Rods[i].oNodes[j][k], 2);
			}
			leng = sqrt(leng);

			Estre += 0.5*ks*pow(leng / leng0 - 1.0, 2)*leng0;
			for (int k = 0; k < 3; k++){
				bvec(vN + k) += ks*(1.0 / leng0 - 1.0 / leng)*(Rods[i].oNodes[j + 1][k] - Rods[i].oNodes[j][k]);
				bvec(vN + 4 + k) -= ks*(1.0 / leng0 - 1.0 / leng)*(Rods[i].oNodes[j + 1][k] - Rods[i].oNodes[j][k]);
				resF(vN + k) += ks*(1.0 / leng0 - 1.0 / leng)*(Rods[i].oNodes[j + 1][k] - Rods[i].oNodes[j][k]);
				resF(vN + 4 + k) -= ks*(1.0 / leng0 - 1.0 / leng)*(Rods[i].oNodes[j + 1][k] - Rods[i].oNodes[j][k]);
				for (int l = 0; l < 3; l++){
					double val = ks / pow(leng, 3)*(Rods[i].oNodes[j + 1][k] - Rods[i].oNodes[j][k])*(Rods[i].oNodes[j + 1][l] - Rods[i].oNodes[j][l]);
					if (k == l) val += ks*(1.0 / leng0 - 1.0 / leng);
					kmlist.push_back(Triplet<double>(vN + k, vN + l, val));
					kmlist.push_back(Triplet<double>(vN + k, vN + 4 + l, -val));
					kmlist.push_back(Triplet<double>(vN + 4 + k, vN + l, -val));
					kmlist.push_back(Triplet<double>(vN + 4 + k, vN + 4 + l, val));
				}
			}
		}


		//ねじり・曲げ成分
		for (int j = 1; j < Rods[i].nodeN - 1; j++){
			int vN = Rods[i].nvlist[j];
			double mbar = twist_angle_calculation(Rods[i].oTPframe[j - 1], Rods[i].oTPframe[j]);
			double tmk = otwis(j) - otwis(j - 1) + mbar;
			double tm0 = Rods[i].tm[j - 1];
			double kb[3], ki[2];
			double kai = 1.0;
			for (int k = 0; k < 3; k++){
				kb[k] = Rods[i].oMframe[j - 1].d3[(k + 1) % 3] * Rods[i].oMframe[j].d3[(k + 2) % 3]
					- Rods[i].oMframe[j - 1].d3[(k + 2) % 3] * Rods[i].oMframe[j].d3[(k + 1) % 3];
				kai += Rods[i].oMframe[j - 1].d3[k] * Rods[i].oMframe[j].d3[k];
			}
			ki[0] = 0.0;
			ki[1] = 0.0;
			for (int k = 0; k < 3; k++){
				kb[k] = 2.0*kb[k] / kai;
				ki[0] += 0.5*(Rods[i].oMframe[j - 1].d2[k] + Rods[i].oMframe[j].d2[k])*kb[k];
				ki[1] -= 0.5*(Rods[i].oMframe[j - 1].d1[k] + Rods[i].oMframe[j].d1[k])*kb[k];
			}

			double leng1 = 0.0;
			double leng2 = 0.0;
			for (int k = 0; k < 3; k++){
				leng1 += pow(Rods[i].oNodes[j][k] - Rods[i].oNodes[j - 1][k], 2);
				leng2 += pow(Rods[i].oNodes[j + 1][k] - Rods[i].oNodes[j][k], 2);
			}
			leng1 = sqrt(leng1);
			leng2 = sqrt(leng2);
			double leni = 0.5*(Rods[i].oLeng[j - 1] + Rods[i].oLeng[j]);

			Etwis += 0.5*kt*(tmk - tm0)*(tmk - tm0) / leni;
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 2; l++){
					Ebend += 0.5*(ki[k] - Rods[i].bk[j - 1][k])*kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni;
				}
			}

			double grme[2][3], hesm[4][3][3];
			double grke[2][2][3], heske[2][4][3][3];
			double grkt[2][2], hket[2][2][2][3], hktt[2][2];
			derivative_of_twist_bend(j, Rods[i].oMframe, leng1, leng2, grme, hesm, grke, grkt, heske, hket, hktt);

			for (int k = 0; k < 3; k++){
				bvec(vN - 4 + k) += kt*(tmk - tm0) / leni*grme[0][k];
				bvec(vN + k) += kt*(tmk - tm0) / leni*(-grme[0][k] + grme[1][k]);
				bvec(vN + 4 + k) += -kt*(tmk - tm0) / leni*grme[1][k];
				resF(vN - 4 + k) += kt*(tmk - tm0) / leni*grme[0][k];
				resF(vN + k) += kt*(tmk - tm0) / leni*(-grme[0][k] + grme[1][k]);
				resF(vN + 4 + k) += -kt*(tmk - tm0) / leni*grme[1][k];
				for (int l = 0; l < 2; l++){
					for (int m = 0; m < 2; m++){
						bvec(vN - 4 + k) += grke[l][0][k] * kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
						bvec(vN + k) += (-grke[l][0][k] + grke[l][1][k])*kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
						bvec(vN + 4 + k) += -grke[l][1][k] * kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
						resF(vN - 4 + k) += grke[l][0][k] * kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
						resF(vN + k) += (-grke[l][0][k] + grke[l][1][k])*kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
						resF(vN + 4 + k) += -grke[l][1][k] * kben[l][m] * (ki[m] - Rods[i].bk[j - 1][m]) / leni;
					}
				}
			}
			bvec(vN - 1) += kt*(tmk - tm0) / leni;
			bvec(vN + 3) -= kt*(tmk - tm0) / leni;
			resF(vN - 1) += kt*(tmk - tm0) / leni;
			resF(vN + 3) -= kt*(tmk - tm0) / leni;
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 2; l++){
					bvec(vN - 1) -= grkt[k][0] * kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni;
					bvec(vN + 3) -= grkt[k][1] * kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni;
					resF(vN - 1) -= grkt[k][0] * kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni;
					resF(vN + 3) -= grkt[k][1] * kben[k][l] * (ki[l] - Rods[i].bk[j - 1][l]) / leni;
				}
			}
			for (int k = 0; k < 3; k++){
				for (int l = 0; l < 3; l++){
					for (int m = 0; m < 2; m++){
						for (int n = 0; n < 2; n++){
							double val = kt*(grme[m][k] * grme[n][l] + (tmk - tm0)*hesm[2 * m + n][k][l]) / leni;
							kmlist.push_back(Triplet<double>(vN + 4 * (m - 1) + k, vN + 4 * (n - 1) + l, val));
							kmlist.push_back(Triplet<double>(vN + 4 * (m - 1) + k, vN + 4 * n + l, -val));
							kmlist.push_back(Triplet<double>(vN + 4 * m + k, vN + 4 * (n - 1) + l, -val));
							kmlist.push_back(Triplet<double>(vN + 4 * m + k, vN + 4 * n + l, val));

							val = 0.0;
							for (int o = 0; o < 2; o++){
								for (int p = 0; p < 2; p++){
									val += (grke[o][m][k] * kben[o][p] * grke[p][n][l] + heske[o][2 * m + n][k][l] * kben[o][p] * (ki[p] - Rods[i].bk[j - 1][p])) / leni;
								}
							}
							kmlist.push_back(Triplet<double>(vN + 4 * (m - 1) + k, vN + 4 * (n - 1) + l, val));
							kmlist.push_back(Triplet<double>(vN + 4 * (m - 1) + k, vN + 4 * n + l, -val));
							kmlist.push_back(Triplet<double>(vN + 4 * m + k, vN + 4 * (n - 1) + l, -val));
							kmlist.push_back(Triplet<double>(vN + 4 * m + k, vN + 4 * n + l, val));
						}
					}
				}
			}
			kmlist.push_back(Triplet<double>(vN - 1, vN - 1, kt / leni));
			kmlist.push_back(Triplet<double>(vN - 1, vN + 3, -kt / leni));
			kmlist.push_back(Triplet<double>(vN + 3, vN - 1, -kt / leni));
			kmlist.push_back(Triplet<double>(vN + 3, vN + 3, kt / leni));
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 2; l++){
					float val = 0.0;
					for (int m = 0; m < 2; m++){
						for (int n = 0; n < 2; n++){
							val += grkt[m][k] * kben[m][n] * grkt[n][l] / leni;
							if (k == l) val += hktt[m][k] * kben[m][n] * (ki[n] - Rods[i].bk[j - 1][n]) / leni;
						}
					}
					kmlist.push_back(Triplet<double>(vN + 4 * (k - 1) + 3, vN + 4 * (l - 1) + 3, val));
				}
			}
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 3; l++){
					double val = kt*grme[k][l] / leni;
					kmlist.push_back(Triplet<double>(vN + 4 * (k - 1) + l, vN - 1, val));
					kmlist.push_back(Triplet<double>(vN + 4 * k + l, vN - 1, -val));
					kmlist.push_back(Triplet<double>(vN + 4 * (k - 1) + l, vN + 3, -val));
					kmlist.push_back(Triplet<double>(vN + 4 * k + l, vN + 3, val));

					kmlist.push_back(Triplet<double>(vN - 1, vN + 4 * (k - 1) + l, val));
					kmlist.push_back(Triplet<double>(vN - 1, vN + 4 * k + l, -val));
					kmlist.push_back(Triplet<double>(vN + 3, vN + 4 * (k - 1) + l, -val));
					kmlist.push_back(Triplet<double>(vN + 3, vN + 4 * k + l, val));
				}
				for (int l = 0; l < 2; l++){
					for (int m = 0; m < 3; m++){
						double val = 0.0;
						for (int n = 0; n < 2; n++){
							for (int o = 0; o < 2; o++){
								val += (grke[n][k][m] * kben[n][o] * grkt[o][l] + hket[n][k][l][m] * kben[n][o] * (ki[o] - Rods[i].bk[j - 1][o])) / leni;
							}
						}
						kmlist.push_back(Triplet<double>(vN + 4 * (k - 1) + m, vN + 4 * (l - 1) + 3, -val));
						kmlist.push_back(Triplet<double>(vN + 4 * k + m, vN + 4 * (l - 1) + 3, val));
						kmlist.push_back(Triplet<double>(vN + 4 * (l - 1) + 3, vN + 4 * (k - 1) + m, -val));
						kmlist.push_back(Triplet<double>(vN + 4 * (l - 1) + 3, vN + 4 * k + m, val));
					}
				}
			}

		}
	}

	//拘束関連
	int cpN = 0.0;
	int cfN = 0.0;
	for (int i = 0; i < conN; i++){
		int vN = cvlist[i];
		cpN += Connect[i].nodeN;
		for (int j = 0; j < Connect[i].nodeN; j++){
			int rN = Connect[i].cNode[j].rN;
			int nN = Connect[i].cNode[j].nN;
			int cvn = Rods[rN].nvlist[nN];
			double Cp[3], crossX[3][3];
			for (int k = 0; k < 3; k++){
				Cp[k] = Rods[rN].oNodes[nN][k] - Connect[i].ogp[k];
				for (int l = 0; l < 3; l++){
					if (k != l && l == (k + 1) % 3) crossX[k][l] = -Connect[i].cNode[j].pos[3 - k - l];
					else if (k != l) crossX[k][l] = Connect[i].cNode[j].pos[3 - k - l];
					else crossX[k][l] = 0.0;
				}
			}

			Eigen::MatrixXf grC(3, 3);
			double Rot0[3][3];
			derive_rot_mat_from_quaternion(Connect[i].oquat, Rot0);
			for (int k = 0; k < 3; k++){
				for (int l = 0; l < 3; l++){
					Cp[k] -= Rot0[k][l] * Connect[i].cNode[j].pos[l];
					grC(k, l) = 0.0;
					for (int m = 0; m < 3; m++){
						grC(k, l) -= Rot0[k][m] * crossX[l][m];
					}
				}
			}
			varCp += sqrt(Cp[0] * Cp[0] + Cp[1] * Cp[1] + Cp[2] * Cp[2]);

			for (int k = 0; k < 3; k++){
				//Econ += Cp[k] * oLamd(vN - cvlist[0] + 3 * j + k);
				bvec(vN + 3 * j + k) -= Cp[k];
				//bvec(cvn + k) -= 1.0*oLamd(vN - cvlist[0] + 3 * j + k);
				//bvec(snodeN + 6 * i + k) -= -1.0*oLamd(vN - cvlist[0] + 3 * j + k);
				for (int l = 0; l < 3; l++){
					//bvec(snodeN + 6 * i + 3 + k) -= grC(l, k)*oLamd(vN - cvlist[0] + 3 * j + l);
				}
			}

			for (int k = 0; k < 3; k++){
				kmlist.push_back(Triplet<double>(cvn + k, vN + 3 * j + k, 1.0));
				kmlist.push_back(Triplet<double>(vN + 3 * j + k, cvn + k, 1.0));
				kmlist.push_back(Triplet<double>(snodeN + 6 * i + k, vN + 3 * j + k, -1.0));
				kmlist.push_back(Triplet<double>(vN + 3 * j + k, snodeN + 6 * i + k, -1.0));
				for (int l = 0; l < 3; l++){
					kmlist.push_back(Triplet<double>(snodeN + 6 * i + 3 + k, vN + 3 * j + l, grC(l, k)));
					kmlist.push_back(Triplet<double>(vN + 3 * j + k, snodeN + 6 * i + 3 + l, grC(k, l)));
				}
			}
		}
		vN += 3 * Connect[i].nodeN;
		cfN += Connect[i].edgeN;
		for (int j = 0; j < Connect[i].edgeN; j++){
			int rN = Connect[i].cEdge[j].rN;
			int eN = Connect[i].cEdge[j].eN;
			int cvn = Rods[rN].evlist[eN];
			double Cf = 0.0;

			double Rot0[3][3], Rd[3][3], rd1[3], rd3[3], rd1h[3];
			derive_rot_mat_from_quaternion(Connect[i].oquat, Rot0);
			for (int k = 0; k < 3; k++){
				rd3[k] = 0.0;
				rd1h[k] = 0.0;
				for (int l = 0; l < 3; l++){
					rd3[k] += Rot0[k][l] * Connect[i].cEdge[j].frame[2][l];
					rd1h[k] += Rot0[k][l] * Connect[i].cEdge[j].frame[0][l];
				}
			}
			double axis[3];
			double seki = 0.0;
			double axh = 0.0;
			for (int k = 0; k < 3; k++){
				axis[k] = rd3[(k + 1) % 3] * Rods[rN].oMframe[eN].d3[(k + 2) % 3]
					- rd3[(k + 2) % 3] * Rods[rN].oMframe[eN].d3[(k + 1) % 3];
				seki += rd3[k] * Rods[rN].oMframe[eN].d3[k];
				axh += axis[k] * axis[k];
			}
			axh = sqrt(axh);
			double quat[4];
			if (axh > 1.0e-10){
				if (seki > 1.0) seki = 1.0;
				if (seki < -1.0) seki = -1.0;
				seki = acos(seki);
				quat[0] = cos(seki / 2.0);
				for (int k = 0; k < 3; k++){
					quat[1 + k] = sin(seki / 2.0)*axis[k] / axh;
				}
			}
			else {
				quat[0] = 1.0; quat[1] = 0.0; quat[2] = 0.0; quat[3] = 0.0;
			}
			derive_rot_mat_from_quaternion(quat, Rd);

			for (int k = 0; k < 3; k++){
				rd1[k] = 0.0;
				for (int l = 0; l < 3; l++){
					rd1[k] += Rd[k][l] * rd1h[l];
				}
			}

			double theta = 0.0;
			double nd[3];
			for (int k = 0; k < 3; k++){
				theta += rd1[k] * Rods[rN].oMframe[eN].d1[k];
				nd[k] = rd1[(k + 1) % 3] * Rods[rN].oMframe[eN].d1[(k + 2) % 3]
					- rd1[(k + 2) % 3] * Rods[rN].oMframe[eN].d1[(k + 1) % 3];
			}
			seki = 0.0;
			for (int k = 0; k < 3; k++){
				seki += nd[k] * Rods[rN].oMframe[eN].d3[k];
			}
			if (theta > 1.0) theta = 1.0;
			if (theta < -1.0) theta = -1.0;
			if (seki > 0) Cf = acos(theta);
			else Cf = -acos(theta);

			varCf += Cf;
			//Econ += Cf*oLamd(vN - cvlist[0] + j);
			bvec(vN + j) -= Cf;
			//bvec(cvn) -= 1.0*oLamd(vN - cvlist[0] + j);
			for (int k = 0; k < 3; k++){
				//bvec(snodeN + 6 * i + 3 + k) -= -Connect[i].cEdge[j].frame[2][k] * oLamd(vN - cvlist[0] + j);
			}

			kmlist.push_back(Triplet<double>(cvn, vN + j, 1.0));
			kmlist.push_back(Triplet<double>(vN + j, cvn, 1.0));
			for (int k = 0; k < 3; k++){
				kmlist.push_back(Triplet<double>(snodeN + 6 * i + 3 + k, vN + j, -Connect[i].cEdge[j].frame[2][k]));
				kmlist.push_back(Triplet<double>(vN + j, snodeN + 6 * i + 3 + k, -Connect[i].cEdge[j].frame[2][k]));
			}
		}
	}


	//剛体変換の拘束
	double c1[3];
	for (int i = 0; i < 3; i++){
		c1[i] = -Gpos[i] * (double)RnodeN;
	}
	for (int i = 0; i < rodN; i++){
		for (int j = 0; j < Rods[i].nodeN; j++){
			int vN = Rods[i].nvlist[j];
			for (int k = 0; k < 3; k++){
				kmlist.push_back(Triplet<double>(vN + k, svarN - 6 + k, 1.0));
				kmlist.push_back(Triplet<double>(svarN - 6 + k, vN + k, 1.0));
				c1[k] += Rods[i].oNodes[j][k];
				for (int l = 0; l < 3; l++){
					double val = Rods[i].oNodes[j][3 - k - l] - Gpos[3 - k - l];
					if (k != l&& l == (k + 1) % 3){
						kmlist.push_back(Triplet<double>(vN + k, svarN - 3 + l, -val));
						kmlist.push_back(Triplet<double>(svarN - 3 + k, vN + l, val));
					}
					else if (k != l){
						kmlist.push_back(Triplet<double>(vN + k, svarN - 3 + l, val));
						kmlist.push_back(Triplet<double>(svarN - 3 + k, vN + l, -val));
					}
				}
			}
		}
	}
	for (int i = 0; i < 3; i++){
		bvec(svarN - 6 + i) -= c1[i];
	}

	varCp /= (float)cpN;
	varCf /= (float)cfN;

	//double wei = 1.0e-7;
	double wei = Ldamp;
	for (int i = snodeN + 6 * conN; i < svarN - 6; i++){
		kmlist.push_back(Triplet<double>(i, i, wei));
		//bvec(i) -= wei*oLamd(i - cvlist[0]);
	}

	Kmat.setFromTriplets(kmlist.begin(), kmlist.end());


	//行列計算
	SparseLU<SparseMatrix<double>> sLU;
	sLU.compute(Kmat);
	bool csuccess = (sLU.info() == Success);
	if (csuccess){
		sol = sLU.solve(bvec);
		Eigen::VectorXd dif = Kmat*sol - bvec;
		if (dif.norm() > 1.0e-4*bvec.norm())csuccess = false;
	}
	if (!csuccess){
		bool CSuccess = false;
		double addmass = 1.0;
		printf("sparse LU failed\n");
		while (1){
			double addDamp = Ldamp*9.0;
			if (Ldamp < 1.0e-4) {
				Ldamp = 1.0e-4;
				addDamp = 1.0e-3;
			}
			Ldamp *= 10.0;
			for (int i = snodeN + 6 * conN; i < svarN - 6; i++){
				bool diag_exist = false;
				for (SparseMatrix<double>::InnerIterator it(Kmat, i); it; ++it){
					if (it.row() == i){
						it.valueRef() += addDamp;
						diag_exist = true;
					}
				}
				if (!diag_exist)Kmat.insert(i, i) = Ldamp;
				//bvec(i) -= addDamp*oLamd(i - cvlist[0]);
			}
			sLU.compute(Kmat);
			CSuccess = (sLU.info() == Success);
			if (CSuccess) {
				sol = sLU.solve(bvec);
				Eigen::VectorXd dif = Kmat*sol - bvec;
				if (dif.norm() < 1.0e-4) break;
				CSuccess = false;
			}
			if (LamSize > 0 && Ldamp*LamSize > sParam->Wid*0.001) break;
			else if (LamSize < 0 && Ldamp > 1.0e-6) break;
		}

		while (1){
			if (CSuccess) break;
			for (int i = 0; i < snodeN + 6 * conN; i++){
				for (SparseMatrix<double>::InnerIterator it(Kmat, i); it; ++it){
					if (i == it.row()){
						bool is_constrained = false;
						if (!is_constrained){
							it.valueRef() += addmass;
						}
					}
				}
			}
			sLU.compute(Kmat);
			CSuccess = (sLU.info() == Success);
			if (CSuccess) {
				sol = sLU.solve(bvec);
				Eigen::VectorXd dif = Kmat*sol - bvec;
				if (dif.norm() < 1.0e-4) break;
				CSuccess = false;
			}
			addmass *= 2.0;
			if (addmass > 1.0e5) break;
		}
		//if (addmass > 1.0e3) tstep_is_large = true;
		if (!CSuccess){
			//printf("Simulation Diverge\n");
			return false;
		}
	}

	double alph = LineSearch(sol);

	Disp = alph*sol;

}


//d3方向を合わせこんだ時のd3周りの回転量
double twist_angle_calculation(struct AdaptiveFrame fr1, struct AdaptiveFrame fr2){

	double ax[3];
	double axh = 0.0;
	double nseki = 0.0;
	for (int i = 0; i < 3; i++){
		ax[i] = fr1.d3[(i + 1) % 3] * fr2.d3[(i + 2) % 3] - fr1.d3[(i + 2) % 3] * fr2.d3[(i + 1) % 3];
		axh += ax[i] * ax[i];
		nseki += fr1.d3[i] * fr2.d3[i];
	}
	axh = sqrt(axh);

	double Rot[3][3];
	if (axh > sin(1.0e-5)){
		double quat[4];
		if (nseki > 1.0) nseki = 1.0;
		if (nseki < -1.0) nseki = -1.0;
		nseki = acos(nseki);
		quat[0] = cos(nseki / 2.0);
		for (int i = 0; i < 3; i++){
			quat[i + 1] = sin(nseki / 2.0)*ax[i] / axh;
		}
		derive_rot_mat_from_quaternion(quat, Rot);
	}
	else {
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				if (i == j) Rot[i][j] = 1.0;
				else Rot[i][j] = 0.0;
			}
		}
	}

	double Rd1[3];
	for (int i = 0; i < 3; i++){
		Rd1[i] = 0.0;
		for (int j = 0; j < 3; j++){
			Rd1[i] += Rot[i][j] * fr1.d1[j];
		}
	}

	nseki = 0.0;
	double oseki[3];
	for (int i = 0; i < 3; i++){
		nseki += Rd1[i] * fr2.d1[i];
		oseki[i] = Rd1[(i + 1) % 3] * fr2.d1[(i + 2) % 3] - Rd1[(i + 2) % 3] * fr2.d1[(i + 1) % 3];
	}
	if (nseki > 1.0) nseki = 1.0;
	if (nseki < -1.0) nseki = -1.0;
	double sig = 0.0;
	for (int i = 0; i < 3; i++){
		sig += oseki[i] * fr2.d3[i];
	}
	if (sig >= 0)return acos(nseki);
	else return -acos(nseki);

}
