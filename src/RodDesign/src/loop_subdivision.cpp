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

#include "geometric_data_structure.h"
#include "loop_subdivision.h"


//コンストラクタ
loop_subdivision_mesh::loop_subdivision_mesh(){

	PREPARED = false;
	Emap5_EXIST = false;
}

//デストラクタ
loop_subdivision_mesh::~loop_subdivision_mesh(){

	for (auto itr = Emap1.begin(); itr != Emap1.end(); ++itr){
		delete itr->second;
	}

	for (auto itr = Emap2.begin(); itr != Emap2.end(); ++itr){
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 12; j++){
				delete[] itr->second->Pmat[i][j];
			}
			delete[] itr->second->Pmat[i];
		}
		delete[] itr->second->Pmat;
		delete itr->second;
	}

	for (auto itr = Emap3.begin(); itr != Emap3.end(); ++itr){
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 12; j++){
				delete[] itr->second->Pmat[i][j];
			}
			delete[] itr->second->Pmat[i];
		}
		delete[] itr->second->Pmat;
		delete itr->second;
	}

	for (auto itr = Emap4.begin(); itr != Emap4.end(); ++itr){
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 12; j++){
				delete[] itr->second->Pmat[i][j];
			}
			delete[] itr->second->Pmat[i];
		}
		delete[] itr->second->Pmat;
		delete itr->second;
	}

	if (Emap5_EXIST){
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 12; j++){
				delete[] Emap5->Pmat[i][j];
			}
			delete[] Emap5->Pmat[i];
		}
		delete[] Emap5->Pmat;
		delete Emap5;
	}

	if(PREPARED) delete[] bcheck;
}


//subdivisionの実行
void loop_subdivision_mesh::subdivision(void){

	bool *nbcheck;
	nbcheck = new bool[nodeN + edgeN];

	//点の追加と位置の更新
	float (*newN)[3];
	newN = new float[nodeN + edgeN][3];

	for (int i = 0; i < edgeN; i++){
		if (EFlist[i][1] != -1){
			int fp1 = -1;
			int fp2 = -1;
			for (int j = 0; j < 3; j++){
				if (Faces[EFlist[i][0]][j] != Edges[i][0] && Faces[EFlist[i][0]][j] != Edges[i][1]) fp1 = Faces[EFlist[i][0]][j];
				if (Faces[EFlist[i][1]][j] != Edges[i][0] && Faces[EFlist[i][1]][j] != Edges[i][1]) fp2 = Faces[EFlist[i][1]][j];
			}
			if (fp1 == -1 || fp2 == -1) printf("opposite face point not found\n");
			for (int j = 0; j < 3; j++){
				newN[nodeN + i][j] = (3.0*Nodes[Edges[i][0]][j] + 3.0*Nodes[Edges[i][1]][j] + Nodes[fp1][j] + Nodes[fp2][j]) / 8.0;
			}
			nbcheck[nodeN + i] = false;
		}
		else {
			for (int j = 0; j < 3; j++){
				newN[nodeN + i][j] = (Nodes[Edges[i][0]][j] + Nodes[Edges[i][1]][j]) / 2.0;
			}
			nbcheck[nodeN + i] = true;
		}
	}

	for (int i = 0; i < nodeN; i++){
		if (!bcheck[i]){
			float beta = (0.625 - pow(3.0 + 2.0*cos(2.0*M_PI / (float)PEsize[i]), 2) / 64.0) / (float)PEsize[i];
			for (int j = 0; j < 3; j++){
				newN[i][j] = (1.0 - beta*(float)PEsize[i])*Nodes[i][j];
			}
			for (int j = 0; j < PEsize[i]; j++){
				int opn = Edges[PElist[i][j]][0];
				if (i == opn) opn = Edges[PElist[i][j]][1];
				for (int j = 0; j < 3; j++){
					newN[i][j] += beta*Nodes[opn][j];
				}
 			}
			nbcheck[i] = false;
		}
		else {
			for (int j = 0; j < 3; j++){
				newN[i][j] = 0.75*Nodes[i][j];
			}
			int bp1 = -1;
			int bp2 = -1;
			for (int j = 0; j < PEsize[i]; j++){
				if (EFlist[PElist[i][j]][1] == -1){
					int oen = Edges[PElist[i][j]][0];
					if (oen == i) oen = Edges[PElist[i][j]][1];
					if (bp1 == -1) bp1 = oen;
					else bp2 = oen;
				}
			}
			if (bp2 == -1) printf("boundary node not found\n");
			for (int j = 0; j < 3; j++){
				newN[i][j] += 0.125*(Nodes[bp1][j] + Nodes[bp2][j]);
			}
			nbcheck[i] = true;
		}
	}

	delete[] Nodes;
	delete[] bcheck;
	Nodes = newN;
	bcheck = nbcheck;
	//nodeN += edgeN;


	//エッジ・面の更新
	int(*newE)[2], (*newF)[3];
	newE = new int[2 * edgeN + 3 * faceN][2];
	newF = new int[4 * faceN][3];

	for (int i = 0; i < edgeN; i++){
		newE[i][0] = Edges[i][0];
		newE[i][1] = nodeN + i;
		newE[edgeN + i][0] = nodeN + i;
		newE[edgeN + i][1] = Edges[i][1];
	}
	for (int i = 0; i < faceN; i++){
		for (int j = 0; j < 3; j++){
			newE[2 * edgeN + 3 * i + j][0] = nodeN + FElist[i][j];
			newE[2 * edgeN + 3 * i + j][1] = nodeN + FElist[i][(j + 1) % 3];
			newF[4 * i + j][0] = Faces[i][j];
			newF[4 * i + j][1] = nodeN + FElist[i][j];
			newF[4 * i + j][2] = nodeN + FElist[i][(j + 2) % 3];
			newF[4 * i + 3][j] = nodeN + FElist[i][j];
		}
	}
	delete[] Edges;
	delete[] Faces;
	Edges = newE;
	Faces = newF;
	//edgeN = 2 * edgeN + 3 * faceN;
	//faceN = 4 * faceN;

	
	//接続関係の更新
	for (int i = 0; i < nodeN; i++){
		delete[] PElist[i];
		delete[] PFlist[i];
	}
	delete[] PEsize;
	delete[] PElist;
	delete[] PFsize;
	delete[] PFlist;
	delete[] FElist, EFlist;

	nodeN += edgeN;
	edgeN = 2 * edgeN + 3 * faceN;
	faceN = 4 * faceN;
	CONSTRUCT_FLG = false;
	construct_data_structure();

}


//spline曲面の評価のための準備
void loop_subdivision_mesh::prepare_for_evaluation(void){

	PREPARED = true;

	bcheck = new bool[nodeN];
	for (int i = 0; i < nodeN; i++){
		bcheck[i] = false;
	}
	bool bound_exist = false;
	for (int i = 0; i < edgeN; i++){
		if (EFlist[i][1] == -1){
			bcheck[Edges[i][0]] = true;
			bcheck[Edges[i][1]] = true;
			bound_exist = true;
		}
	}

	if (bound_exist) subdivision();

	std::vector<int> Vlist1, Vlist2, Vlist3;
	std::vector<std::pair<int, int>> Vlist4;
	bool type3_exist = false;
	while (1){
		bool subdivision_is_necessary = false;
		
		//三角形のカテゴリーチェック
		for (int i = 0; i < faceN; i++){
			bool extp[3];
			int extpN = 0;
			int exN = -1;
			for (int j = 0; j < 3; j++){
				if (!bcheck[Faces[i][j]]){
					if (PFsize[Faces[i][j]] == 6) extp[j] = false;
					else {
						extp[j] = true;
						extpN++;
						exN = j;
					}
				}
				else {
					if (PFsize[Faces[i][j]] == 3) extp[j] = false;
					else {
						extp[j] = true;
						extpN++;
						exN = j;
					}
				}
			}
			if (extpN > 1){
				subdivision_is_necessary = true;
				break;
			}
			if (exN != -1 && !bcheck[Faces[i][exN]]){
				int PEN = PEsize[Faces[i][exN]];
				for (int j = 0; j < Vlist1.size(); j++){
					if (PEN == Vlist1[j]) PEN = -1;
				}
				if (PEN != -1) Vlist1.push_back(PEN);
			}
			else if (exN != -1 && bcheck[Faces[i][exN]]){
				if (!bcheck[Faces[i][(exN + 1) % 3]] && !bcheck[Faces[i][(exN + 2) % 3]]){
					int cfn = i;
					int cpn = Faces[i][(exN + 1) % 3];
					int cen = FElist[i][exN];
					int lvn = 1;
					while (1){
						int nfn = EFlist[cen][0];
						if (nfn == cfn) nfn = EFlist[cen][1];
						if (nfn == -1) break;
						int opn = -1;
						cen = -1;
						for (int j = 0; j < 3; j++){
							if (Faces[nfn][j] != Faces[i][exN] && Faces[nfn][j] != cpn){
								opn = Faces[nfn][j];
								cen = FElist[nfn][(j + 2) % 3];
							}
						}
						if (cen == -1)printf("cen not found@preparation\n");
						cfn = nfn;
						cpn = opn;
						lvn++;
					}
					cfn = i;
					cpn = Faces[i][(exN + 2) % 3];
					cen = FElist[i][(exN + 2) % 3];
					int rvn = 1;
					while (1){
						int nfn = EFlist[cen][0];
						if (nfn == cfn) nfn = EFlist[cen][1];
						if (nfn == -1) break;
						int opn = -1;
						cen = -1;
						for (int j = 0; j < 3; j++){
							if (Faces[nfn][j] != Faces[i][exN] && Faces[nfn][j] != cpn){
								opn = Faces[nfn][j];
								cen = FElist[nfn][j];
							}
						}
						if (cen == -1) printf("cen not found@preparation\n");
						cfn = nfn;
						cpn = opn;
						rvn++;
					}
					Vlist4.push_back(make_pair(lvn, lvn + rvn));
				}
				else if (bcheck[Faces[i][(exN + 1) % 3]] && !bcheck[Faces[i][(exN + 2) % 3]]){
					if (EFlist[FElist[i][exN]][1] != -1) printf("Assumption not satisfied\n");
					int PEN = PEsize[Faces[i][exN]];
					for (int j = 0; j < Vlist2.size(); j++){
						if (PEN == Vlist2[j]) PEN = -1;
					}
					if (PEN != -1) Vlist2.push_back(PEN);
				}
				else if (!bcheck[Faces[i][(exN + 1) % 3]] && bcheck[Faces[i][(exN + 2) % 3]]){
					if (EFlist[FElist[i][(exN + 2) % 3]][1] != -1) printf("Assumption not satisfied\n");
					int PEN = PEsize[Faces[i][exN]];
					for (int j = 0; j < Vlist3.size(); j++){
						if (PEN == Vlist3[j]) PEN = -1;
					}
					if (PEN != -1) Vlist3.push_back(PEN);
				}
				else {
					if (EFlist[FElist[i][exN]][1] != -1 || EFlist[FElist[i][(exN + 2) % 3]][1] != -1) printf("Assumption not satisfied\n");
					type3_exist = true;
				}
			}
		}

		if (subdivision_is_necessary){
			//printf("subdivide for regularize\n");
			subdivision();
		}
		else break;
	}


	//境界点を持たない三角形のEigen structure構築
	for (int i = 0; i < Vlist1.size(); i++){
		Emap1[Vlist1[i]] = new struct eigen_struct;

		int valen = Vlist1[i];
		Eigen::MatrixXf Smat(valen + 6, valen + 6);
		Smat = Eigen::MatrixXf::Zero(valen + 6, valen + 6);
		Emap1[valen]->Amat = Eigen::MatrixXf::Zero(valen + 12, valen + 6);

		float beta = (0.625 - pow(3 + 2.0*cos(2.0*M_PI / (float)valen), 2) / 64.0) / (float)valen;
		Smat(0, 0) = 1.0 - (float)valen * beta;
		for (int j = 0; j < valen; j++){
			Smat(0, 1 + j) = beta;
			Smat(1 + j, 0) = 0.375;
			Smat(1 + j, 1 + j) = 0.375;
			Smat(1 + j, 1 + (j + 1) % valen) = 0.125;
			Smat(1 + j, 1 + (j + valen - 1) % valen) = 0.125;
		}
		for (int j = 0; j < 3; j++){
			int off = 0;
			if (j == 1) off = 1;
			else if (j == 2) off = -1;
			Smat(valen + 1 + 2 * j, 0) = 0.125;
			Smat(valen + 1 + 2 * j, valen + 1 + 2 * j) = 0.125;
			Smat(valen + 1 + 2 * j, 1 + (valen - 1 + off) % valen) = 0.375;
			Smat(valen + 1 + 2 * j, 1 + (valen + off) % valen) = 0.375;

			if (j != 2){
				Smat(valen + 2 + 2 * j, 1 + (valen - j) % valen) = 0.625;
				Smat(valen + 2 + 2 * j, 0) = 0.0625;
				Smat(valen + 2 + 2 * j, 2 - j) = 0.0625;
				Smat(valen + 2 + 2 * j, valen - j) = 0.0625;
				Smat(valen + 2 + 2 * j, valen + 1) = 0.0625;
				Smat(valen + 2 + 2 * j, valen + 2 + 2 * j) = 0.0625;
				Smat(valen + 2 + 2 * j, valen + 3 + 2 * j) = 0.0625;
			}
		}

		for (int j = 0; j < valen + 6; j++){
			for (int k = 0; k < valen + 6; k++){
				Emap1[valen]->Amat(j, k) = Smat(j, k);
			}
		}
		for (int j = 0; j < 3; j++){
			Emap1[valen]->Amat(valen + 6 + j, 1) = 0.375;
			Emap1[valen]->Amat(valen + 6 + j, valen + 1 + j) = 0.375;
			Emap1[valen]->Amat(valen + 6 + j, valen + j) = 0.125;
			if(j != 2) Emap1[valen]->Amat(valen + 6 + j, valen + 2 + j) = 0.125;
			else Emap1[valen]->Amat(valen + 6 + j, 2) = 0.125;

			Emap1[valen]->Amat(valen + 9 + j, valen) = 0.375;
			if (j != 0) Emap1[valen]->Amat(valen + 9 + j, valen + 3 + j) = 0.375;
			else Emap1[valen]->Amat(valen + 9 + j, valen + 1) = 0.375;
			if (j != 2) Emap1[valen]->Amat(valen + 9 + j, valen + 4 + j) = 0.125;
			else Emap1[valen]->Amat(valen + 9 + j, valen - 1) = 0.125;
			if (j != 0) Emap1[valen]->Amat(valen + 9 + j, valen - 2 + 3 * j) = 0.125;
			else Emap1[valen]->Amat(valen + 9 + j, 1) = 0.125;
		}

		//固有値解析
		Eigen::EigenSolver<Eigen::MatrixXf> ES(Smat);
		Emap1[valen]->evals = ES.eigenvalues();
		Emap1[valen]->evecs = ES.eigenvectors();

		if (ES.info() != Eigen::Success){
			printf("Eigen decomposion failed (v=%d no bound)\n", valen);
		}

		//Permutation matrixの設定
		Emap1[valen]->Pmap[0][0] = 2;  Emap1[valen]->Pmap[0][1] = 0;  Emap1[valen]->Pmap[0][2] = valen + 3;
		Emap1[valen]->Pmap[0][3] = 1;  Emap1[valen]->Pmap[0][4] = valen;  Emap1[valen]->Pmap[0][5] = valen + 8;
		Emap1[valen]->Pmap[0][6] = valen + 2;  Emap1[valen]->Pmap[0][7] = valen + 1;  Emap1[valen]->Pmap[0][8] = valen + 4;
		Emap1[valen]->Pmap[0][9] = valen + 7;  Emap1[valen]->Pmap[0][10] = valen + 6;  Emap1[valen]->Pmap[0][11] = valen + 9;
		Emap1[valen]->Pmap[1][0] = valen + 9;  Emap1[valen]->Pmap[1][1] = valen + 6;  Emap1[valen]->Pmap[1][2] = valen + 4;
		Emap1[valen]->Pmap[1][3] = valen + 1;  Emap1[valen]->Pmap[1][4] = valen + 2;  Emap1[valen]->Pmap[1][5] = valen + 5;
		Emap1[valen]->Pmap[1][6] = valen;  Emap1[valen]->Pmap[1][7] = 1;   Emap1[valen]->Pmap[1][8] = valen + 3;
		Emap1[valen]->Pmap[1][9] = valen - 1;  Emap1[valen]->Pmap[1][10] = 0;  Emap1[valen]->Pmap[1][11] = 2;
		Emap1[valen]->Pmap[2][0] = 0;  Emap1[valen]->Pmap[2][1] = valen - 1;  Emap1[valen]->Pmap[2][2] = 1;
		Emap1[valen]->Pmap[2][3] = valen;  Emap1[valen]->Pmap[2][4] = valen + 5;  Emap1[valen]->Pmap[2][5] = valen + 2;
		Emap1[valen]->Pmap[2][6] = valen + 1;  Emap1[valen]->Pmap[2][7] = valen + 4;  Emap1[valen]->Pmap[2][8] = valen + 11;
		Emap1[valen]->Pmap[2][9] = valen + 6;  Emap1[valen]->Pmap[2][10] = valen + 9;  Emap1[valen]->Pmap[2][11] = valen + 10;
	}


	//extraordinary点のみ境界
	for (int i = 0; i < Vlist4.size(); i++){
		Emap4[Vlist4[i]] = new struct eigen_struct;
		int vn = Vlist4[i].first;
		int valen = Vlist4[i].second;
		Eigen::MatrixXf Smat(valen + 6, valen + 6);
		Smat = Eigen::MatrixXf::Zero(valen + 6, valen + 6);
		Emap4[Vlist4[i]]->Amat = Eigen::MatrixXf::Zero(valen + 12, valen + 6);

		Smat(0, 0) = 0.75;
		Smat(0, vn) = 0.125;
		Smat(0, vn + 1) = 0.125;
		for (int j = 0; j < 2; j++){
			Smat(vn + j, 0) = 0.5;
			Smat(vn + j, vn + j) = 0.5;
		}
		for (int j = 0; j < vn - 1; j++){
			Smat(1 + j, 0) = 0.375;
			Smat(1 + j, 1 + j) = 0.375;
			Smat(1 + j, 2 + j) = 0.125;
			Smat(1 + j, 1 + (valen - 1 + j) % valen) = 0.125;
		}
		for (int j = vn + 2; j < valen + 1; j++){
			Smat(j, 0) = 0.375;
			Smat(j, j) = 0.375;
			Smat(j, j - 1) = 0.125;
			Smat(j, 1 + j % valen) = 0.125;
		}
		for (int j = 0; j < 3; j++){
			int off = 0;
			if (j == 1) off = 1;
			else if (j == 2) off = -1;
			Smat(valen + 1 + 2 * j, 0) = 0.125;
			Smat(valen + 1 + 2 * j, valen + 1 + 2 * j) = 0.125;
			Smat(valen + 1 + 2 * j, 1 + (valen + off - 1) % valen) = 0.375;
			Smat(valen + 1 + 2 * j, 1 + (valen + off) % valen) = 0.375;
		}
		Smat(valen + 2, 1) = 0.625;
		Smat(valen + 4, valen) = 0.625;
		for (int j = 0; j < 3; j++){
			Smat(valen + 2, valen + 1 + j) = 0.0625;
			if (j != 2) Smat(valen + 2, 2 * j) = 0.0625;
			else Smat(valen + 2, valen) = 0.0625;
			if (j != 2) Smat(valen + 4, j) = 0.0625;
			else Smat(valen + 4, valen - 1) = 0.0625;
			if (j != 2) Smat(valen + 4, valen + 4 + j) = 0.0625;
			else Smat(valen + 4, valen + 1) = 0.0625;
		}

		for (int j = 0; j < valen + 6; j++){
			for (int k = 0; k < valen + 6; k++){
				Emap4[Vlist4[i]]->Amat(j, k) = Smat(j, k);
			}
		}
		for (int j = 0; j < 3; j++){
			Emap4[Vlist4[i]]->Amat(valen + 6 + j, 1) = 0.375;
			Emap4[Vlist4[i]]->Amat(valen + 6 + j, valen + 1 + j) = 0.375;
			Emap4[Vlist4[i]]->Amat(valen + 6 + j, valen + j) = 0.125;
			if (j != 2) Emap4[Vlist4[i]]->Amat(valen + 6 + j, valen + 2 + j) = 0.125;
			else Emap4[Vlist4[i]]->Amat(valen + 6 + j, 2) = 0.125;
			Emap4[Vlist4[i]]->Amat(valen + 9 + j, valen) = 0.375;
			if (j != 0) Emap4[Vlist4[i]]->Amat(valen + 9 + j, valen + 3 + j) = 0.375;
			else Emap4[Vlist4[i]]->Amat(valen + 9 + j, valen + 1) = 0.375;
			if (j != 2) Emap4[Vlist4[i]]->Amat(valen + 9 + j, valen + 4 + j) = 0.125;
			else Emap4[Vlist4[i]]->Amat(valen + 9 + j, valen - 1) = 0.125;
			if (j != 0) Emap4[Vlist4[i]]->Amat(valen + 9 + j, valen - 2 + 3 * j) = 0.125;
			else Emap4[Vlist4[i]]->Amat(valen + 9 + j, 1) = 0.125;
		}

		//固有値解析
		Eigen::EigenSolver<Eigen::MatrixXf> ES(Smat);
		Emap4[Vlist4[i]]->evals = ES.eigenvalues();
		Emap4[Vlist4[i]]->evecs = ES.eigenvectors();

		if (ES.info() != Eigen::Success){
			printf("Eigen decomposion failed (v=%d bound type2-A)\n", valen);
		}

		//Permutation matrixの設定
		Emap4[Vlist4[i]]->Pmat = new int**[3];
		for (int j = 0; j < 3; j++){
			Emap4[Vlist4[i]]->Pmat[j] = new int*[12];
			for (int k = 0; k < 12; k++){
				Emap4[Vlist4[i]]->Pmat[j][k] = new int[valen + 12];
				for (int l = 0; l < valen + 12; l++){
					Emap4[Vlist4[i]]->Pmat[j][k][l] = 0;
				}
			}
		}
		int ***perm;
		perm = Emap4[Vlist4[i]]->Pmat;
		perm[0][0][2] = 1;  perm[0][1][0] = 1;  perm[0][2][valen + 3] = 1;  perm[0][3][1] = 1;  perm[0][4][valen] = 1;
		perm[0][5][valen + 8] = 1;  perm[0][6][valen + 2] = 1;  perm[0][7][valen + 1] = 1;  perm[0][8][valen + 4] = 1;
		perm[0][9][valen + 7] = 1;  perm[0][10][valen + 6] = 1;  perm[0][11][valen + 9] = 1; 
		perm[1][0][valen + 9] = 1;  perm[1][1][valen + 6] = 1;  perm[1][2][valen + 4] = 1;  perm[1][3][valen + 1] = 1;
		perm[1][4][valen + 2] = 1;  perm[1][5][valen + 5] = 1;  perm[1][6][valen] = 1;  perm[1][7][1] = 1;
		perm[1][8][valen + 3] = 1;  perm[1][9][valen - 1] = 1;  perm[1][10][0] = 1;  perm[1][11][2] = 1;
		perm[2][0][0] = 1;  perm[2][1][valen - 1] = 1;  perm[2][2][1] = 1;  perm[2][3][valen] = 1; 
		perm[2][4][valen + 5] = 1;  perm[2][5][valen + 2] = 1;  perm[2][6][valen + 1] = 1;  perm[2][7][valen + 4] = 1;
		perm[2][8][valen + 11] = 1;  perm[2][9][valen + 6] = 1;  perm[2][10][valen + 9] = 1;  perm[2][11][valen + 10] = 1;

	}


	//境界エッジ1つのみ（その１）
	for (int i = 0; i < Vlist2.size(); i++){
		Emap2[Vlist2[i]] = new struct eigen_struct;
		int valen = Vlist2[i];
		Eigen::MatrixXf Smat(valen + 5, valen + 5);
		Smat = Eigen::MatrixXf::Zero(valen + 5, valen + 5);
		Emap2[valen]->Amat = Eigen::MatrixXf::Zero(valen + 10, valen + 5);

		Smat(0, 0) = 0.75;
		Smat(0, 1) = 0.125;
		Smat(0, 2) = 0.125;
		Smat(1, 0) = 0.5;
		Smat(1, 1) = 0.5;
		Smat(2, 0) = 0.5;
		Smat(2, 2) = 0.5;
		for (int j = 0; j < valen - 2; j++){
			Smat(j + 3, 0) = 0.375;
			Smat(j + 3, j + 3) = 0.375;
			Smat(j + 3, 2 + j) = 0.125;
			Smat(j + 3, 1 + (3 + j) % valen) = 0.125;
		}
		Smat(valen + 1, 0) = 0.125;
		Smat(valen + 1, valen + 1) = 0.125;
		Smat(valen + 1, 1) = 0.375;
		Smat(valen + 1, valen) = 0.375;
		Smat(valen + 2, 1) = 0.75;
		Smat(valen + 2, 0) = 0.125;
		Smat(valen + 2, valen + 2) = 0.125;
		Smat(valen + 4, valen - 1) = 0.375;
		Smat(valen + 4, valen) = 0.375;
		Smat(valen + 4, 0) = 0.125;
		Smat(valen + 4, valen + 4) = 0.125;
		Smat(valen + 3, valen) = 0.625;
		for (int j = 0; j < 3; j++){
			if (j != 2)Smat(valen + 3, j) = 0.0625;
			else Smat(valen + 3, valen - 1) = 0.0625;
			if (j != 0) Smat(valen + 3, valen + 2 + j) = 0.0625;
			else Smat(valen + 3, valen + 1) = 0.0625;
		}


		for (int j = 0; j < valen + 5; j++){
			for (int k = 0; k < valen + 5; k++){
				Emap2[valen]->Amat(j, k) = Smat(j, k);
			}
		}

		Emap2[valen]->Amat(valen + 6, 1) = 0.5;
		Emap2[valen]->Amat(valen + 6, valen + 2) = 0.5;
		Emap2[valen]->Amat(valen + 5, 1) = 0.375;
		Emap2[valen]->Amat(valen + 5, valen + 1) = 0.375;
		Emap2[valen]->Amat(valen + 5, valen + 2) = 0.125;
		Emap2[valen]->Amat(valen + 5, valen) = 0.125;
		for (int j = 0; j < 3; j++){
			Emap2[valen]->Amat(valen + 7 + j, valen) = 0.375;
			if(j != 0) Emap2[valen]->Amat(valen + 7 + j, valen + 2 + j) = 0.375;
			else Emap2[valen]->Amat(valen + 7 + j, valen + 1) = 0.375;
			if (j != 2) Emap2[valen]->Amat(valen + 7 + j, valen + 3 + j) = 0.125;
			else Emap2[valen]->Amat(valen + 7 + j, valen - 1) = 0.125;
			if (j != 0) Emap2[valen]->Amat(valen + 7 + j, valen + 2 * j - 1) = 0.125;
			else Emap2[valen]->Amat(valen + 7 + j, 1) = 0.125;
		}

		//固有値解析
		Eigen::EigenSolver<Eigen::MatrixXf> ES(Smat);
		Emap2[valen]->evals = ES.eigenvalues();
		Emap2[valen]->evecs = ES.eigenvectors();

		if (ES.info() != Eigen::Success){
			printf("Eigen decomposion failed (v=%d bound type2-A)\n", valen);
		}

		//Permutation matrixの設定
		Emap2[valen]->Pmat = new int**[3];
		for (int j = 0; j < 3; j++){
			Emap2[valen]->Pmat[j] = new int*[12];
			for (int k = 0; k < 12; k++){
				Emap2[valen]->Pmat[j][k] = new int[valen + 10];
				for (int l = 0; l < valen + 10; l++){
					Emap2[valen]->Pmat[j][k][l] = 0;
				}
			}
		}
		int ***perm;
		perm = Emap2[valen]->Pmat;
		perm[0][0][0] = 1;  perm[0][0][1] = 1;  perm[0][0][valen] = -1;  perm[0][1][0] = 1;
		perm[0][2][1] = 1;  perm[0][2][valen + 2] = 1;  perm[0][2][valen + 1] = -1;  perm[0][3][1] = 1;
		perm[0][4][valen] = 1;  perm[0][5][valen + 2] = 1;  perm[0][5][valen + 6] = 1;  perm[0][5][valen + 5] = -1;
		perm[0][6][valen + 2] = 1;  perm[0][7][valen + 1] = 1;  perm[0][8][valen + 3] = 1;  perm[0][9][valen + 6] = 1;
		perm[0][10][valen + 5] = 1;  perm[0][11][valen + 7] = 1;
		perm[1][0][valen + 7] = 1;  perm[1][1][valen + 5] = 1;  perm[1][2][valen + 3] = 1;  perm[1][3][valen + 1] = 1;
		perm[1][4][valen + 2] = 1;  perm[1][5][valen + 4] = 1;  perm[1][6][valen] = 1;  perm[1][7][1] = 1;
		perm[1][8][valen + 2] = 1;  perm[1][8][1] = 1;  perm[1][8][valen + 1] = -1;  perm[1][9][valen - 1] = 1;
		perm[1][10][0] = 1; perm[1][11][0] = 1;  perm[1][11][1] = 1;  perm[1][11][valen] = -1;
		perm[2][0][0] = 1;  perm[2][1][valen - 1] = 1;  perm[2][2][1] = 1; perm[2][3][valen] = 1;
		perm[2][4][valen + 4] = 1;  perm[2][5][valen + 2] = 1;  perm[2][6][valen + 1] = 1; perm[2][7][valen + 3] = 1;
		perm[2][8][valen + 9] = 1;  perm[2][9][valen + 5] = 1; perm[2][10][valen + 7] = 1;  perm[2][11][valen + 8] = 1;
	}

	//境界エッジ1つのみ（その２）
	for (int i = 0; i < Vlist3.size(); i++){
		Emap3[Vlist3[i]] = new struct eigen_struct;
		int valen = Vlist3[i];
		Eigen::MatrixXf Smat(valen + 5, valen + 5);
		Smat = Eigen::MatrixXf::Zero(valen + 5, valen + 5);
		Emap3[valen]->Amat = Eigen::MatrixXf::Zero(valen + 10, valen + 5);

		Smat(0, 0) = 0.75;
		Smat(0, valen - 1) = 0.125;
		Smat(0, valen) = 0.125;
		for (int j = 0; j < valen - 2; j++){
			Smat(1 + j, 0) = 0.375;
			Smat(1 + j, 1 + j) = 0.375;
			Smat(1 + j, 1 + (valen - 1 + j) % valen) = 0.125;
			Smat(1 + j, 2 + j) = 0.125;
		}
		for (int j = 0; j < 2; j++){
			Smat(valen - 1 + j, 0) = 0.5;
			Smat(valen - 1 + j, valen - 1 + j) = 0.5;
			Smat(valen + 1 + 2 * j, 0) = 0.125;
			Smat(valen + 1 + 2 * j, 1 + (valen - 1 + j) % valen) = 0.375;
			Smat(valen + 1 + 2 * j, 1 + j) = 0.375;
			Smat(valen + 1 + 2 * j, valen + 1 + 2 * j) = 0.125;
		}
		Smat(valen + 2, 1) = 0.625;
		for (int j = 0; j < 3; j++){
			Smat(valen + 2, valen + 1 + j) = 0.0625;
			if (j != 2)Smat(valen + 2, 2 * j) = 0.0625;
			else Smat(valen + 2, valen) = 0.0625;
		}
		Smat(valen + 4, valen) = 0.75;
		Smat(valen + 4, 0) = 0.125;
		Smat(valen + 4, valen + 4) = 0.125;

		for (int j = 0; j < valen + 5; j++){
			for (int k = 0; k < valen + 5; k++){
				Emap3[valen]->Amat(j, k) = Smat(j, k);
			}
		}
		for (int j = 0; j < 3; j++){
			Emap3[valen]->Amat(valen + 5 + j, 1) = 0.375;
			Emap3[valen]->Amat(valen + 5 + j, valen + 1 + j) = 0.375;
			Emap3[valen]->Amat(valen + 5 + j, valen + j) = 0.125;
			if (j != 2) Emap3[valen]->Amat(valen + 5 + j, valen + 2 + j) = 0.125;
			else Emap3[valen]->Amat(valen + 5 + j, 2) = 0.125;
		}
		Emap3[valen]->Amat(valen + 8, valen) = 0.375;
		Emap3[valen]->Amat(valen + 8, valen + 1) = 0.375;
		Emap3[valen]->Amat(valen + 8, 1) = 0.125;
		Emap3[valen]->Amat(valen + 8, valen + 4) = 0.125;
		Emap3[valen]->Amat(valen + 9, valen) = 0.5;
		Emap3[valen]->Amat(valen + 9, valen + 4) = 0.5;

		//固有値解析
		Eigen::EigenSolver<Eigen::MatrixXf> ES(Smat);
		Emap3[valen]->evals = ES.eigenvalues();
		Emap3[valen]->evecs = ES.eigenvectors();

		if (ES.info() != Eigen::Success){
			printf("Eigen decomposion failed (v=%d bound type2-B)\n", valen);
		}

		//Permutation matrixの設定
		Emap3[valen]->Pmat = new int**[3];
		for (int j = 0; j < 3; j++){
			Emap3[valen]->Pmat[j] = new int*[12];
			for (int k = 0; k < 12; k++){
				Emap3[valen]->Pmat[j][k] = new int[valen + 10];
				for (int l = 0; l < valen + 10; l++){
					Emap3[valen]->Pmat[j][k][l] = 0;
				}
			}
		}
		int ***perm;
		perm = Emap3[valen]->Pmat;
		perm[0][0][2] = 1;  perm[0][1][0] = 1;  perm[0][2][valen + 3] = 1;  perm[0][3][1] = 1;  perm[0][4][valen] = 1;
		perm[0][5][valen + 7] = 1;  perm[0][6][valen + 2] = 1;  perm[0][7][valen + 1] = 1;  perm[0][8][valen + 4] = 1;
		perm[0][9][valen + 6] = 1;  perm[0][10][valen + 5] = 1;  perm[0][11][valen + 8] = 1;
		perm[1][0][valen + 8] = 1;  perm[1][1][valen + 5] = 1;  perm[1][2][valen + 4] = 1;  perm[1][3][valen + 1] = 1;
		perm[1][4][valen + 2] = 1;  perm[1][5][valen] = 1;  perm[1][5][valen + 4] = 1;  perm[1][5][valen + 1] = -1;
		perm[1][6][valen] = 1;  perm[1][7][1] = 1;  perm[1][8][valen + 3] = 1;  perm[1][9][0] = 1;  perm[1][9][valen] = 1;
		perm[1][9][1] = -1;  perm[1][10][0] = 1;  perm[1][11][2] = 1;
		perm[2][0][0] = 1;  perm[2][1][0] = 1;  perm[2][1][valen] = 1;  perm[2][1][1] = -1;  perm[2][2][1] = 1;
		perm[2][3][valen] = 1;  perm[2][4][valen] = 1;  perm[2][4][valen + 4] = 1;  perm[2][4][valen + 1] = -1;
		perm[2][5][valen + 2] = 1;  perm[2][6][valen + 1] = 1;  perm[2][7][valen + 4] = 1;  perm[2][8][valen + 4] = 1;
		perm[2][8][valen + 9] = 1;  perm[2][8][valen + 8] = -1;  perm[2][9][valen + 5] = 1;  perm[2][10][valen + 8] = 1;
		perm[2][11][valen + 9] = 1;
	}


	//2辺共に境界
	if (type3_exist){
		Emap5_EXIST = true;
		Emap5 = new struct eigen_struct;
		Eigen::MatrixXf Smat(6, 6);
		Smat = Eigen::MatrixXf::Zero(6, 6);
		Emap5->Amat = Eigen::MatrixXf::Zero(10, 10);

		Smat(0, 0) = 0.75;
		Smat(0, 1) = 0.125;
		Smat(0, 2) = 0.125;
		for (int i = 0; i < 2; i++){
			Smat(1 + i, 0) = 0.5;
			Smat(1 + i, 1 + i) = 0.5;
			Smat(4 + i, 0) = 0.125;
			Smat(4 + i, 1 + i) = 0.75;
			Smat(4 + i, 4 + i) = 0.125;
		}
		Smat(3, 0) = 0.125;
		Smat(3, 1) = 0.375;
		Smat(3, 2) = 0.375;
		Smat(3, 3) = 0.125;

		for (int i = 0; i < 6; i++){
			for (int j = 0; j < 6; j++){
				Emap5->Amat(i, j) = Smat(i, j);
			}
		}
		for (int i = 0; i < 2; i++){
			Emap5->Amat(6 + 2 * i, 3) = 0.375;
			Emap5->Amat(6 + 2 * i, 1 + i) = 0.375;
			Emap5->Amat(6 + 2 * i, 2 - i) = 0.125;
			Emap5->Amat(6 + 2 * i, 4 + i) = 0.125;
			Emap5->Amat(7 + 2 * i, 1 + i) = 0.5;
			Emap5->Amat(7 + 2 * i, 4 + i) = 0.5;
		}

		//固有値解析
		Eigen::EigenSolver<Eigen::MatrixXf> ES(Smat);
		Emap5->evals = ES.eigenvalues();
		Emap5->evecs = ES.eigenvectors();
		
		if (ES.info() != Eigen::Success){
			printf("Eigen decomposion failed (bound type3)\n");
		}

		//Permutation matrixの設定
		Emap5->Pmat = new int**[3];
		for (int i = 0; i < 3; i++){
			Emap5->Pmat[i] = new int*[12];
			for (int j = 0; j < 12; j++){
				Emap5->Pmat[i][j] = new int[10];
				for (int k = 0; k < 10; k++){
					Emap5->Pmat[i][j][k] = 0;
				}
			}
		}
		int ***perm;
		perm = Emap5->Pmat;
		perm[0][0][0] = 1;  perm[0][0][1] = 1;  perm[0][0][2] = -1;  perm[0][1][0] = 1;  perm[0][2][1] = 1; perm[0][2][4] = 1;   perm[0][2][3] = -1;
		perm[0][3][1] = 1;  perm[0][4][2] = 1;  perm[0][5][4] = 1;  perm[0][5][7] = 1;  perm[0][5][6] = -1;  perm[0][6][4] = 1;  perm[0][7][3] = 1;
		perm[0][8][5] = 1;  perm[0][9][7] = 1;  perm[0][10][6] = 1;  perm[0][11][8] = 1;
		perm[1][0][8] = 1;  perm[1][1][6] = 1;  perm[1][2][5] = 1;  perm[1][3][3] = 1;  perm[1][4][4] = 1;  perm[1][5][2] = 1;  perm[1][5][5] = 1;
		perm[1][5][3] = -1;  perm[1][6][2] = 1;  perm[1][7][1] = 1;  perm[1][8][1] = 1;  perm[1][8][4] = 1;  perm[1][8][3] = -1;  perm[1][9][0] = 1;
		perm[1][9][2] = 1;  perm[1][9][1] = -1;  perm[1][10][0] = 1;  perm[1][11][0] = 1;  perm[1][11][1] = 1;  perm[1][11][2] = -1;
		perm[2][0][0] = 1;  perm[2][1][0] = 1;  perm[2][1][2] = 1;  perm[2][1][1] = -1;  perm[2][2][1] = 1;  perm[2][3][2] = 1;  perm[2][4][2] = 1;
		perm[2][4][5] = 1;  perm[2][4][3] = -1;  perm[2][5][4] = 1;  perm[2][6][3] = 1;  perm[2][7][5] = 1;  perm[2][8][5] = 1;  perm[2][8][9] = 1;
		perm[2][8][8] = -1;  perm[2][9][6] = 1;  perm[2][10][8] = 1;  perm[2][11][9] = 1;
	}
}


//spline曲面の評価
bool loop_subdivision_mesh::evaluate_box_spline_on_triangle(int fn, float uvw[3], float bsp[3]){

	//三角形のカテゴリーチェック
	bool extp[3],ebound[3];
	int exN = -1;
	int extpN = 0;
	for (int i = 0; i < 3; i++){
		if (!bcheck[Faces[fn][i]]){
			if (PFsize[Faces[fn][i]] == 6) extp[i] = false;
			else {
				extp[i] = true;
				exN = i;
				extpN++;
			}
		}
		else {
			if (PFsize[Faces[fn][i]] == 3) extp[i] = false;
			else {
				extp[i] = true;
				exN = i;
				extpN++;
			}
		}
		if (EFlist[FElist[fn][i]][1] == -1) ebound[i] = true;
		else ebound[i] = false;
	}

	if (extpN > 1) {
		printf("more subdivision is necessary\n");
		return false;
	}

	//すべてregular
	if (extpN == 0){

		int neib[3][6];
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 6; j++){
				neib[i][j] = -1;
			}
			neib[i][0] = Faces[fn][(i + 1) % 3];
			neib[i][5] = Faces[fn][(i + 2) % 3];
			int cpn = neib[i][0];
			int cen = FElist[fn][i];
			int cfn = fn;
			int neiN = 1;
			while (1){
				int nfn = EFlist[cen][0];
				if (nfn == cfn) nfn = EFlist[cen][1];
				if (nfn == -1 || nfn == fn) break;
				
				int opn = -1;
				cen = -1;
				for (int j = 0; j < 3; j++){
					if (Faces[nfn][j] != Faces[fn][i] && Faces[nfn][j] != cpn) {
						opn = Faces[nfn][j];
						cen = FElist[nfn][(j + 2) % 3];
					}
				}
				if (opn == -1 || cen == -1) printf("cen or ofn not found@regular case\n");

				if (neiN == 6) printf("valence is over 6@regular case\n");
				neib[i][neiN] = opn;
				neiN++;
				cpn = opn;
				cfn = nfn;
				if (bcheck[Faces[fn][i]] && neiN == 4) printf("extraordinary boundary is found@regular boundary case\n");
			}
			if (bcheck[Faces[fn][i]]){
				cpn = neib[i][5];
				cen = FElist[fn][(i + 2) % 3];
				cfn = fn;
				neiN = 4;
				while (1){
					int nfn = EFlist[cen][0];
					if (nfn == cfn) nfn = EFlist[cen][1];
					if (nfn == -1 || nfn == fn) break;

					int ofn = -1;
					cen = -1;
					for (int j = 0; j < 3; j++){
						if (Faces[nfn][j] != Faces[fn][i] && Faces[nfn][j] != cpn) {
							ofn = Faces[nfn][j];
							cen = FElist[nfn][j];
						}
					}
					if (ofn == -1) printf("cen and ofn not found@regular boundary case\n");

					neib[i][neiN] = ofn;
					neiN--;
					cpn = ofn;
					cfn = nfn;
					if (neiN == 1){
						printf("valence is over 6@regular boundary case\n");
					}
				}
			}
		}

		//box-spline評価用のcontrol point設定
		float cpset[12][3];
		set_control_points_for_box_spline_with_reflection(fn,neib, cpset);
		evaluate_box_spline(cpset, uvw, bsp);
	}
	//Stamの方法
	else {
		//隣接リストの決定
		vector<vector<int>> neib;
		for (int i = 0; i < 3; i++){
			vector<int> nbi;
			nbi.push_back(Faces[fn][(exN + i + 1) % 3]);
			int cpn = nbi[0];
			int cen = FElist[fn][(exN + i) % 3];
			int cfn = fn;
			while (1){
				int nfn = EFlist[cen][0];
				if (nfn == cfn) nfn = EFlist[cen][1];
				if (nfn == -1 || nfn == fn) break;

				int opn = -1;
				cen = -1;
				for (int j = 0; j < 3; j++){
					if (Faces[nfn][j] != Faces[fn][(exN + i) % 3] && Faces[nfn][j] != cpn){
						opn = Faces[nfn][j];
						cen = FElist[nfn][(j + 2) % 3];
					}
				}
				if (opn == -1) printf("cen & opn not found@extraordinary case\n");

				nbi.push_back(opn);
				cpn = opn;
				cfn = nfn;
			}
			if (bcheck[Faces[fn][(exN + i) % 3]]){
				vector<int> inbi;
				inbi.push_back(Faces[fn][(exN + i + 2) % 3]);
				cpn = inbi[0];
				cen = FElist[fn][(exN + i + 2) % 3];
				cfn = fn;
				while (1){
					int nfn = EFlist[cen][0];
					if (nfn == cfn) nfn = EFlist[cen][1];
					if (nfn == -1 || nfn == fn) break;

					int opn = -1;
					cen = -1;
					for (int j = 0; j < 3; j++){
						if (Faces[nfn][j] != Faces[fn][(exN + i) % 3] && Faces[nfn][j] != cpn){
							opn = Faces[nfn][j];
							cen = FElist[nfn][j];
						}
					}
					if (opn == -1) printf("opn & cen not found@extraordinary case\n");

					inbi.push_back(opn);
					cpn = opn;
					cfn = nfn;
				}
				nbi.push_back(-1);
				for (int j = 0; j < inbi.size(); j++){
					nbi.push_back(inbi[inbi.size() - 1 - j]);
				}
			}
			neib.push_back(nbi);
		}

		//パラメータの変換とlayerの決定
		float suvw[3];
		for (int j = 0; j < 3; j++){
			suvw[j] = uvw[(exN + j) % 3];
		}
		int lvm;
		if (fabs(suvw[1] + suvw[2]) > 1.0e-6)lvm = floor(1.0 - log(suvw[1] + suvw[2]) / log(2));
		else {
			lvm = floor(1.0 - log(1.0e-6) / log(2));
			suvw[1] = 5.0e-7;
			suvw[2] = 5.0e-7;
		}
		int Pk;
		suvw[1] *= pow(2.0, lvm - 1);
		suvw[2] *= pow(2.0, lvm - 1);
		if (suvw[1] > 0.5){
			Pk = 0;
			suvw[1] = 2.0*suvw[1] - 1.0;
			suvw[2] = 2.0*suvw[2];
		}
		else if (suvw[2] > 0.5){
			Pk = 2;
			suvw[1] = 2.0*suvw[1];
			suvw[2] = 2.0*suvw[2] - 1.0;
		}
		else {
			Pk = 1;
			suvw[1] = 1.0 - 2.0*suvw[1];
			suvw[2] = 1.0 - 2.0*suvw[2];
		}
		suvw[0] = 1.0 - suvw[1] - suvw[2];


		if (!bcheck[Faces[fn][0]] && !bcheck[Faces[fn][1]] && !bcheck[Faces[fn][2]]){

			//すべて内点の場合
			float (*C0)[3];
			C0 = new float[neib[0].size() + 6][3];
			int valen = neib[0].size();
			for (int j = 0; j < 3; j++){
				C0[0][j] = Nodes[Faces[fn][exN]][j];
				for (int k = 0; k < valen; k++){
					C0[1 + k][j] = Nodes[neib[0][k]][j];
				}
				for (int k = 0; k < 3; k++){
					C0[valen + 1 + k][j] = Nodes[neib[1][k + 1]][j];
				}
				C0[valen + 4][j] = Nodes[neib[2][3]][j];
				C0[valen + 5][j] = Nodes[neib[2][2]][j];
			}

			//サンプル点の生成
			Eigen::MatrixXcf Sm(valen + 6, valen + 6),Dmat;
			Eigen::VectorXcf mD(valen + 6);
			for (int j = 0; j < valen + 6; j++){
				mD(j) = pow(Emap1[valen]->evals(j), lvm - 1);
			}
			Dmat = mD.asDiagonal();
			Sm = Emap1[valen]->evecs * Dmat*Emap1[valen]->evecs.inverse();

			Eigen::MatrixXf Cm(valen+12,3);
			for (int j = 0; j < valen + 12; j++){
				for (int k = 0; k < 3; k++){
					Cm(j, k) = 0.0;
					for (int l = 0; l < valen + 6; l++){
						for (int m = 0; m < valen + 6; m++){
							Cm(j, k) += Emap1[valen]->Amat(j, l)*Sm(l, m).real()*C0[m][k];
						}
					}
				}
			}

			//Permutation matrixとスプライン基底の乗算
			float cpset[12][3];
			for (int j = 0; j < 12; j++){
				for (int k = 0; k < 3; k++){
					cpset[j][k] = Cm(Emap1[valen]->Pmap[Pk][j], k);
				}
			}
			evaluate_box_spline(cpset, suvw, bsp);

			delete[] C0;
		}
		else if (bcheck[Faces[fn][exN]] && !bcheck[Faces[fn][(exN+1)%3]] && !bcheck[Faces[fn][(exN+2)%3]]){
			//境界点を含む場合 type1
			int valen = neib[0].size() - 1; 
			int vn = -1;
			for (int j = 0; j < neib[0].size(); j++){
				if (neib[0][j] == -1) vn = j;
			}
			std::pair<int, int> Vpair;
			Vpair = make_pair(vn, valen);
			
			float(*C0)[3];
			C0 = new float[valen + 6][3];
			int nN = 0;
			for (int j = 0; j < 3; j++){
				C0[0][j] = Nodes[Faces[fn][exN]][j];
				nN = 0;
				for (int k = 0; k < neib[0].size(); k++){
					if (neib[0][k] != -1) {
						C0[1 + nN][j] = Nodes[neib[0][k]][j];
						nN++;
					}
				}
				for (int k = 0; k < 3; k++){
					C0[valen + 1 + k][j] = Nodes[neib[1][1 + k]][j];
					if (k != 2) C0[valen + 4 + k][j] = Nodes[neib[2][3 - k]][j];
				}
			}
			if (nN != valen) printf("neighbor size is wrong@type1\n");

			//サンプル点の生成
			Eigen::MatrixXcf Sm(valen + 6, valen + 6), Dmat;
			Eigen::VectorXcf mD(valen + 6);
			for (int j = 0; j < valen + 6; j++){
				mD(j) = pow(Emap4[Vpair]->evals(j), lvm - 1);
			}
			Dmat = mD.asDiagonal();
			Sm = Emap4[Vpair]->evecs * Dmat*Emap4[Vpair]->evecs.inverse();

			Eigen::MatrixXf Cm(valen + 12, 3);
			for (int j = 0; j < valen + 12; j++){
				for (int k = 0; k < 3; k++){
					Cm(j, k) = 0.0;
					for (int l = 0; l < valen + 6; l++){
						for (int m = 0; m < valen + 6; m++){
							Cm(j, k) += Emap4[Vpair]->Amat(j, l)*Sm(l, m).real()*C0[m][k];
						}
					}
				}
			}

			//Permutation matrixとスプライン基底の乗算
			float cpset[12][3];
			for (int j = 0; j < 12; j++){
				for (int k = 0; k < 3; k++){
					cpset[j][k] = 0.0;
					for (int l = 0; l < valen + 12; l++){
						cpset[j][k] += (float)Emap4[Vpair]->Pmat[Pk][j][l] * Cm(l, k);
					}
				}
			}
			evaluate_box_spline(cpset, suvw, bsp);
			
			delete[] C0;
			
		}
		else if (bcheck[Faces[fn][exN]] && bcheck[Faces[fn][(exN + 1) % 3]] && !bcheck[Faces[fn][(exN + 2) % 3]]){
			//境界エッジを含む場合 type2-A
			float(*C0)[3];
			C0 = new float[neib[0].size() + 4][3];
			int valen = neib[0].size() - 1;
			for (int j = 0; j < 3; j++){
				C0[0][j] = Nodes[Faces[fn][exN]][j];
				C0[1][j] = Nodes[neib[0][0]][j];
				for (int k = 2; k < neib[0].size(); k++){
					C0[k][j] = Nodes[neib[0][k]][j];
				}
				for (int k = 0; k < 2; k++){
					C0[valen + 1 + k][j] = Nodes[neib[1][1 + k]][j];
					C0[valen + 3 + k][j] = Nodes[neib[2][3 - k]][j];
				}
			}
			if (neib[0][1] != -1 || neib[1][3] != -1 )printf("boundary position is wrong @type2-A\n");
			if (neib[1].size() != 5 || neib[2].size() != 6) printf("size okashii\n");

			//サンプル点の生成
			Eigen::MatrixXcf Sm(valen + 5, valen + 5), Dmat;
			Eigen::VectorXcf mD(valen + 5);
			for (int j = 0; j < valen + 5; j++){
				mD(j) = pow(Emap2[valen]->evals(j), lvm - 1);
			}
			Dmat = mD.asDiagonal();
			Sm = Emap2[valen]->evecs * Dmat*Emap2[valen]->evecs.inverse();

			Eigen::MatrixXf Cm(valen + 10, 3);
			for (int j = 0; j < valen + 10; j++){
				for (int k = 0; k < 3; k++){
					Cm(j, k) = 0.0;
					for (int l = 0; l < valen + 5; l++){
						for (int m = 0; m < valen + 5; m++){
							Cm(j, k) += Emap2[valen]->Amat(j, l)*Sm(l, m).real()*C0[m][k];
						}
					}
				}
			}

			//Permutation matrixとスプライン基底の乗算
			float cpset[12][3];
			for (int j = 0; j < 12; j++){
				for (int k = 0; k < 3; k++){
					cpset[j][k] = 0.0;
					for (int l = 0; l < valen + 10; l++){
						cpset[j][k] += (float)Emap2[valen]->Pmat[Pk][j][l] * Cm(l, k);
					}
				}
			}
			evaluate_box_spline(cpset, suvw, bsp);

			delete[] C0;

		}
		else if (bcheck[Faces[fn][exN]] && !bcheck[Faces[fn][(exN + 1) % 3]] && bcheck[Faces[fn][(exN + 2) % 3]]){
			//境界エッジを含む場合 type2-B
			float(*C0)[3];
			C0 = new float[neib[0].size() + 4][3];
			int valen = neib[0].size() - 1;
			for (int j = 0; j < 3; j++){
				C0[0][j] = Nodes[Faces[fn][exN]][j];
				for (int k = 0; k < neib[0].size() - 2; k++){
					C0[1 + k][j] = Nodes[neib[0][k]][j];
				}
				C0[valen][j] = Nodes[neib[0][valen]][j];
				for (int k = 0; k < 3; k++){
					C0[valen + 1 + k][j] = Nodes[neib[1][1 + k]][j];
				}
				C0[valen + 4][j] = Nodes[neib[2][2]][j];
			}
			if (neib[0][valen-1] != -1 || neib[2][1] != -1)printf("boundary position is wrong @type2-B\n");

			//サンプル点の生成
			Eigen::MatrixXcf Sm(valen + 5, valen + 5), Dmat;
			Eigen::VectorXcf mD(valen + 5);
			for (int j = 0; j < valen + 5; j++){
				mD(j) = pow(Emap3[valen]->evals(j), lvm - 1);
			}
			Dmat = mD.asDiagonal();
			Sm = Emap3[valen]->evecs * Dmat*Emap3[valen]->evecs.inverse();

			Eigen::MatrixXf Cm(valen + 10, 3);
			for (int j = 0; j < valen + 10; j++){
				for (int k = 0; k < 3; k++){
					Cm(j, k) = 0.0;
					for (int l = 0; l < valen + 5; l++){
						for (int m = 0; m < valen + 5; m++){
							Cm(j, k) += Emap3[valen]->Amat(j, l)*Sm(l, m).real()*C0[m][k];
						}
					}
				}
			}

			//Permutation matrixとスプライン基底の乗算
			float cpset[12][3];
			for (int j = 0; j < 12; j++){
				for (int k = 0; k < 3; k++){
					cpset[j][k] = 0.0;
					for (int l = 0; l < valen + 10; l++){
						cpset[j][k] += (float)Emap3[valen]->Pmat[Pk][j][l] * Cm(l, k);
					}
				}
			}
			evaluate_box_spline(cpset, suvw, bsp);

			delete[] C0;
		}
		else {
			//境界エッジを含む場合 type3
			float(*C0)[3];
			C0 = new float[6][3];
			int valen = 2;
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 3; k++){
					C0[k][j] = Nodes[Faces[fn][(exN + k) % 3]][j];
				}
				for (int k = 0; k < 2; k++){
					C0[3 + k][j] = Nodes[neib[1][1 + k]][j];
				}
				C0[5][j] = Nodes[neib[2][2]][j];
			}
			if (neib[0][1] != -1 || neib[1][3] != -1 || neib[2][1] != -1)printf("boundary position is wrong @type3\n");

			//サンプル点の生成
			Eigen::MatrixXcf Sm(6, 6), Dmat;
			Eigen::VectorXcf mD(6);
			for (int j = 0; j < 6; j++){
				mD(j) = pow(Emap5->evals(j), lvm - 1);
			}
			Dmat = mD.asDiagonal();
			Sm = Emap5->evecs * Dmat*Emap5->evecs.inverse();

			Eigen::MatrixXf Cm(10, 3);
			for (int j = 0; j < 10; j++){
				for (int k = 0; k < 3; k++){
					Cm(j, k) = 0.0;
					for (int l = 0; l < 6; l++){
						for (int m = 0; m < 6; m++){
							Cm(j, k) += Emap5->Amat(j, l)*Sm(l, m).real()*C0[m][k];
						}
					}
				}
			}

			//Permutation matrixとスプライン基底の乗算
			float cpset[12][3];
			for (int j = 0; j < 12; j++){
				for (int k = 0; k < 3; k++){
					cpset[j][k] = 0.0;
					for (int l = 0; l < 10; l++){
						cpset[j][k] += (float)Emap5->Pmat[Pk][j][l] * Cm(l, k);
					}
				}
			}
			evaluate_box_spline(cpset, suvw, bsp);

			delete[] C0;
		}
	}


	return true;

}


//splineの微分の評価
bool loop_subdivision_mesh::evaluate_derivative_of_box_spline_on_triangle(int fn, float uvw[3], float dbs[2][3]){


	//三角形のカテゴリーチェック
	bool extp[3], ebound[3];
	int exN = -1;
	int extpN = 0;
	for (int i = 0; i < 3; i++){
		if (!bcheck[Faces[fn][i]]){
			if (PFsize[Faces[fn][i]] == 6) extp[i] = false;
			else {
				extp[i] = true;
				exN = i;
				extpN++;
			}
		}
		else {
			if (PFsize[Faces[fn][i]] == 3) extp[i] = false;
			else {
				extp[i] = true;
				exN = i;
				extpN++;
			}
		}
		if (EFlist[FElist[fn][i]][1] == -1) ebound[i] = true;
		else ebound[i] = false;
	}

	if (extpN > 1) {
		printf("more subdivision is necessary\n");
		return false;
	}

	//すべてregular
	if (extpN == 0){

		int neib[3][6];
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 6; j++){
				neib[i][j] = -1;
			}
			neib[i][0] = Faces[fn][(i + 1) % 3];
			neib[i][5] = Faces[fn][(i + 2) % 3];
			int cpn = neib[i][0];
			int cen = FElist[fn][i];
			int cfn = fn;
			int neiN = 1;
			while (1){
				int nfn = EFlist[cen][0];
				if (nfn == cfn) nfn = EFlist[cen][1];
				if (nfn == -1 || nfn == fn) break;

				int opn = -1;
				cen = -1;
				for (int j = 0; j < 3; j++){
					if (Faces[nfn][j] != Faces[fn][i] && Faces[nfn][j] != cpn) {
						opn = Faces[nfn][j];
						cen = FElist[nfn][(j + 2) % 3];
					}
				}
				if (opn == -1 || cen == -1) printf("cen or ofn not found@regular case\n");

				if (neiN == 6) printf("valence is over 6@regular case\n");
				neib[i][neiN] = opn;
				neiN++;
				cpn = opn;
				cfn = nfn;
				if (bcheck[Faces[fn][i]] && neiN == 4) printf("extraordinary boundary is found@regular boundary case\n");
			}
			if (bcheck[Faces[fn][i]]){
				cpn = neib[i][5];
				cen = FElist[fn][(i + 2) % 3];
				cfn = fn;
				neiN = 4;
				while (1){
					int nfn = EFlist[cen][0];
					if (nfn == cfn) nfn = EFlist[cen][1];
					if (nfn == -1 || nfn == fn) break;

					int ofn = -1;
					cen = -1;
					for (int j = 0; j < 3; j++){
						if (Faces[nfn][j] != Faces[fn][i] && Faces[nfn][j] != cpn) {
							ofn = Faces[nfn][j];
							cen = FElist[nfn][j];
						}
					}
					if (ofn == -1) printf("cen and ofn not found@regular boundary case\n");

					neib[i][neiN] = ofn;
					neiN--;
					cpn = ofn;
					cfn = nfn;
					if (neiN == 1){
						printf("valence is over 6@regular boundary case\n");
					}
				}
			}
		}

		//box-spline評価用のcontrol point設定
		float cpset[12][3];
		set_control_points_for_box_spline_with_reflection(fn, neib, cpset);
		evaluate_derivative_of_box_spline(cpset, uvw, dbs);
	}
	//Stamの方法
	else {
		//隣接リストの決定
		vector<vector<int>> neib;
		for (int i = 0; i < 3; i++){
			vector<int> nbi;
			nbi.push_back(Faces[fn][(exN + i + 1) % 3]);
			int cpn = nbi[0];
			int cen = FElist[fn][(exN + i) % 3];
			int cfn = fn;
			while (1){
				int nfn = EFlist[cen][0];
				if (nfn == cfn) nfn = EFlist[cen][1];
				if (nfn == -1 || nfn == fn) break;

				int opn = -1;
				cen = -1;
				for (int j = 0; j < 3; j++){
					if (Faces[nfn][j] != Faces[fn][(exN + i) % 3] && Faces[nfn][j] != cpn){
						opn = Faces[nfn][j];
						cen = FElist[nfn][(j + 2) % 3];
					}
				}
				if (opn == -1) printf("cen & opn not found@extraordinary case\n");

				nbi.push_back(opn);
				cpn = opn;
				cfn = nfn;
			}
			if (bcheck[Faces[fn][(exN + i) % 3]]){
				vector<int> inbi;
				inbi.push_back(Faces[fn][(exN + i + 2) % 3]);
				cpn = inbi[0];
				cen = FElist[fn][(exN + i + 2) % 3];
				cfn = fn;
				while (1){
					int nfn = EFlist[cen][0];
					if (nfn == cfn) nfn = EFlist[cen][1];
					if (nfn == -1 || nfn == fn) break;

					int opn = -1;
					cen = -1;
					for (int j = 0; j < 3; j++){
						if (Faces[nfn][j] != Faces[fn][(exN + i) % 3] && Faces[nfn][j] != cpn){
							opn = Faces[nfn][j];
							cen = FElist[nfn][j];
						}
					}
					if (opn == -1) printf("opn & cen not found@extraordinary case\n");

					inbi.push_back(opn);
					cpn = opn;
					cfn = nfn;
				}
				nbi.push_back(-1);
				for (int j = 0; j < inbi.size(); j++){
					nbi.push_back(inbi[inbi.size() - 1 - j]);
				}
			}
			neib.push_back(nbi);
		}

		//パラメータの変換とlayerの決定
		float suvw[3];
		for (int j = 0; j < 3; j++){
			suvw[j] = uvw[(exN + j) % 3];
		}
		int lvm;
		if (fabs(suvw[1] + suvw[2]) > 1.0e-6)lvm = floor(1.0 - log(suvw[1] + suvw[2]) / log(2));
		else {
			lvm = floor(1.0 - log(1.0e-6) / log(2));
			suvw[1] = 5.0e-7;
			suvw[2] = 5.0e-7;
		}
		int Pk;
		suvw[1] *= pow(2.0, lvm - 1);
		suvw[2] *= pow(2.0, lvm - 1);
		if (suvw[1] > 0.5){
			Pk = 0;
			suvw[1] = 2.0*suvw[1] - 1.0;
			suvw[2] = 2.0*suvw[2];
		}
		else if (suvw[2] > 0.5){
			Pk = 2;
			suvw[1] = 2.0*suvw[1];
			suvw[2] = 2.0*suvw[2] - 1.0;
		}
		else {
			Pk = 1;
			suvw[1] = 1.0 - 2.0*suvw[1];
			suvw[2] = 1.0 - 2.0*suvw[2];
		}
		suvw[0] = 1.0 - suvw[1] - suvw[2];


		if (!bcheck[Faces[fn][0]] && !bcheck[Faces[fn][1]] && !bcheck[Faces[fn][2]]){

			//すべて内点の場合
			float(*C0)[3];
			C0 = new float[neib[0].size() + 6][3];
			int valen = neib[0].size();
			for (int j = 0; j < 3; j++){
				C0[0][j] = Nodes[Faces[fn][exN]][j];
				for (int k = 0; k < valen; k++){
					C0[1 + k][j] = Nodes[neib[0][k]][j];
				}
				for (int k = 0; k < 3; k++){
					C0[valen + 1 + k][j] = Nodes[neib[1][k + 1]][j];
				}
				C0[valen + 4][j] = Nodes[neib[2][3]][j];
				C0[valen + 5][j] = Nodes[neib[2][2]][j];
			}

			//サンプル点の生成
			Eigen::MatrixXcf Sm(valen + 6, valen + 6), Dmat;
			Eigen::VectorXcf mD(valen + 6);
			for (int j = 0; j < valen + 6; j++){
				mD(j) = pow(Emap1[valen]->evals(j), lvm - 1);
			}
			Dmat = mD.asDiagonal();
			Sm = Emap1[valen]->evecs * Dmat*Emap1[valen]->evecs.inverse();

			Eigen::MatrixXf Cm(valen + 12, 3);
			for (int j = 0; j < valen + 12; j++){
				for (int k = 0; k < 3; k++){
					Cm(j, k) = 0.0;
					for (int l = 0; l < valen + 6; l++){
						for (int m = 0; m < valen + 6; m++){
							Cm(j, k) += Emap1[valen]->Amat(j, l)*Sm(l, m).real()*C0[m][k];
						}
					}
				}
			}

			//Permutation matrixとスプライン基底の乗算
			float cpset[12][3];
			for (int j = 0; j < 12; j++){
				for (int k = 0; k < 3; k++){
					cpset[j][k] = Cm(Emap1[valen]->Pmap[Pk][j], k);
				}
			}
			float mdbs[2][3];
			evaluate_derivative_of_box_spline(cpset, suvw, mdbs);
			float nin = pow(2.0, lvm);
			for (int j = 0; j < 2; j++){
				for (int k = 0; k < 3; k++){
					if(Pk != 1) mdbs[j][k] = nin*mdbs[j][k];
					else mdbs[j][k] = -nin*mdbs[j][k];
				}
			}
			if (exN == 0){
				for (int i = 0; i < 2; i++){
					for (int j = 0; j < 3; j++){
						dbs[i][j] = mdbs[i][j];
					}
				}
			}
			else if (exN == 1){
				for (int i = 0; i < 3; i++){
					dbs[0][i] = -mdbs[1][i];
					dbs[1][i] = mdbs[0][i] - mdbs[1][i];
				}
			}
			else {
				for (int i = 0; i < 3; i++){
					dbs[0][i] = mdbs[1][i] - mdbs[0][i];
					dbs[1][i] = -mdbs[0][i];
				}
			}
			delete[] C0;
		}
		else if (bcheck[Faces[fn][exN]] && !bcheck[Faces[fn][(exN + 1) % 3]] && !bcheck[Faces[fn][(exN + 2) % 3]]){
			//境界点を含む場合 type1
			int valen = neib[0].size() - 1;
			int vn = -1;
			for (int j = 0; j < neib[0].size(); j++){
				if (neib[0][j] == -1) vn = j;
			}
			std::pair<int, int> Vpair;
			Vpair = make_pair(vn, valen);

			float(*C0)[3];
			C0 = new float[valen + 6][3];
			int nN = 0;
			for (int j = 0; j < 3; j++){
				C0[0][j] = Nodes[Faces[fn][exN]][j];
				nN = 0;
				for (int k = 0; k < neib[0].size(); k++){
					if (neib[0][k] != -1) {
						C0[1 + nN][j] = Nodes[neib[0][k]][j];
						nN++;
					}
				}
				for (int k = 0; k < 3; k++){
					C0[valen + 1 + k][j] = Nodes[neib[1][1 + k]][j];
					if (k != 2) C0[valen + 4 + k][j] = Nodes[neib[2][3 - k]][j];
				}
			}
			if (nN != valen) printf("neighbor size is wrong@type1\n");

			//サンプル点の生成
			Eigen::MatrixXcf Sm(valen + 6, valen + 6), Dmat;
			Eigen::VectorXcf mD(valen + 6);
			for (int j = 0; j < valen + 6; j++){
				mD(j) = pow(Emap4[Vpair]->evals(j), lvm - 1);
			}
			Dmat = mD.asDiagonal();
			Sm = Emap4[Vpair]->evecs * Dmat*Emap4[Vpair]->evecs.inverse();

			Eigen::MatrixXf Cm(valen + 12, 3);
			for (int j = 0; j < valen + 12; j++){
				for (int k = 0; k < 3; k++){
					Cm(j, k) = 0.0;
					for (int l = 0; l < valen + 6; l++){
						for (int m = 0; m < valen + 6; m++){
							Cm(j, k) += Emap4[Vpair]->Amat(j, l)*Sm(l, m).real()*C0[m][k];
						}
					}
				}
			}

			//Permutation matrixとスプライン基底の乗算
			float cpset[12][3];
			for (int j = 0; j < 12; j++){
				for (int k = 0; k < 3; k++){
					cpset[j][k] = 0.0;
					for (int l = 0; l < valen + 12; l++){
						cpset[j][k] += (float)Emap4[Vpair]->Pmat[Pk][j][l] * Cm(l, k);
					}
				}
			}
			float mdbs[2][3];
			evaluate_derivative_of_box_spline(cpset, suvw, mdbs);
			float nin = pow(2.0, lvm);
			for (int j = 0; j < 2; j++){
				for (int k = 0; k < 3; k++){
					if (Pk != 1) mdbs[j][k] = nin*mdbs[j][k];
					else mdbs[j][k] = -nin*mdbs[j][k];
				}
			}
			if (exN == 0){
				for (int i = 0; i < 2; i++){
					for (int j = 0; j < 3; j++){
						dbs[i][j] = mdbs[i][j];
					}
				}
			}
			else if (exN == 1){
				for (int i = 0; i < 3; i++){
					dbs[0][i] = -mdbs[1][i];
					dbs[1][i] = mdbs[0][i] - mdbs[1][i];
				}
			}
			else {
				for (int i = 0; i < 3; i++){
					dbs[0][i] = mdbs[1][i] - mdbs[0][i];
					dbs[1][i] = -mdbs[0][i];
				}
			}
			delete[] C0;
		}
		else if (bcheck[Faces[fn][exN]] && bcheck[Faces[fn][(exN + 1) % 3]] && !bcheck[Faces[fn][(exN + 2) % 3]]){
			//境界エッジを含む場合 type2-A
			float(*C0)[3];
			C0 = new float[neib[0].size() + 4][3];
			int valen = neib[0].size() - 1;
			for (int j = 0; j < 3; j++){
				C0[0][j] = Nodes[Faces[fn][exN]][j];
				C0[1][j] = Nodes[neib[0][0]][j];
				for (int k = 2; k < neib[0].size(); k++){
					C0[k][j] = Nodes[neib[0][k]][j];
				}
				for (int k = 0; k < 2; k++){
					C0[valen + 1 + k][j] = Nodes[neib[1][1 + k]][j];
					C0[valen + 3 + k][j] = Nodes[neib[2][3 - k]][j];
				}
			}
			if (neib[0][1] != -1 || neib[1][3] != -1)printf("boundary position is wrong @type2-A\n");
			if (neib[1].size() != 5 || neib[2].size() != 6) printf("size okashii\n");

			//サンプル点の生成
			Eigen::MatrixXcf Sm(valen + 5, valen + 5), Dmat;
			Eigen::VectorXcf mD(valen + 5);
			for (int j = 0; j < valen + 5; j++){
				mD(j) = pow(Emap2[valen]->evals(j), lvm - 1);
			}
			Dmat = mD.asDiagonal();
			Sm = Emap2[valen]->evecs * Dmat*Emap2[valen]->evecs.inverse();

			Eigen::MatrixXf Cm(valen + 10, 3);
			for (int j = 0; j < valen + 10; j++){
				for (int k = 0; k < 3; k++){
					Cm(j, k) = 0.0;
					for (int l = 0; l < valen + 5; l++){
						for (int m = 0; m < valen + 5; m++){
							Cm(j, k) += Emap2[valen]->Amat(j, l)*Sm(l, m).real()*C0[m][k];
						}
					}
				}
			}

			//Permutation matrixとスプライン基底の乗算
			float cpset[12][3];
			for (int j = 0; j < 12; j++){
				for (int k = 0; k < 3; k++){
					cpset[j][k] = 0.0;
					for (int l = 0; l < valen + 10; l++){
						cpset[j][k] += (float)Emap2[valen]->Pmat[Pk][j][l] * Cm(l, k);
					}
				}
			}
			float mdbs[2][3];
			evaluate_derivative_of_box_spline(cpset, suvw, mdbs);
			float nin = pow(2.0, lvm);
			for (int j = 0; j < 2; j++){
				for (int k = 0; k < 3; k++){
					if (Pk != 1) mdbs[j][k] = nin*mdbs[j][k];
					else mdbs[j][k] = -nin*mdbs[j][k];
				}
			}
			if (exN == 0){
				for (int i = 0; i < 2; i++){
					for (int j = 0; j < 3; j++){
						dbs[i][j] = mdbs[i][j];
					}
				}
			}
			else if (exN == 1){
				for (int i = 0; i < 3; i++){
					dbs[0][i] = -mdbs[1][i];
					dbs[1][i] = mdbs[0][i] - mdbs[1][i];
				}
			}
			else {
				for (int i = 0; i < 3; i++){
					dbs[0][i] = mdbs[1][i] - mdbs[0][i];
					dbs[1][i] = -mdbs[0][i];
				}
			}
			delete[] C0;
		}
		else if (bcheck[Faces[fn][exN]] && !bcheck[Faces[fn][(exN + 1) % 3]] && bcheck[Faces[fn][(exN + 2) % 3]]){
			//境界エッジを含む場合 type2-B
			float(*C0)[3];
			C0 = new float[neib[0].size() + 4][3];
			int valen = neib[0].size() - 1;
			for (int j = 0; j < 3; j++){
				C0[0][j] = Nodes[Faces[fn][exN]][j];
				for (int k = 0; k < neib[0].size() - 2; k++){
					C0[1 + k][j] = Nodes[neib[0][k]][j];
				}
				C0[valen][j] = Nodes[neib[0][valen]][j];
				for (int k = 0; k < 3; k++){
					C0[valen + 1 + k][j] = Nodes[neib[1][1 + k]][j];
				}
				C0[valen + 4][j] = Nodes[neib[2][2]][j];
			}
			if (neib[0][valen - 1] != -1 || neib[2][1] != -1)printf("boundary position is wrong @type2-B\n");

			//サンプル点の生成
			Eigen::MatrixXcf Sm(valen + 5, valen + 5), Dmat;
			Eigen::VectorXcf mD(valen + 5);
			for (int j = 0; j < valen + 5; j++){
				mD(j) = pow(Emap3[valen]->evals(j), lvm - 1);
			}
			Dmat = mD.asDiagonal();
			Sm = Emap3[valen]->evecs * Dmat*Emap3[valen]->evecs.inverse();

			Eigen::MatrixXf Cm(valen + 10, 3);
			for (int j = 0; j < valen + 10; j++){
				for (int k = 0; k < 3; k++){
					Cm(j, k) = 0.0;
					for (int l = 0; l < valen + 5; l++){
						for (int m = 0; m < valen + 5; m++){
							Cm(j, k) += Emap3[valen]->Amat(j, l)*Sm(l, m).real()*C0[m][k];
						}
					}
				}
			}

			//Permutation matrixとスプライン基底の乗算
			float cpset[12][3];
			for (int j = 0; j < 12; j++){
				for (int k = 0; k < 3; k++){
					cpset[j][k] = 0.0;
					for (int l = 0; l < valen + 10; l++){
						cpset[j][k] += (float)Emap3[valen]->Pmat[Pk][j][l] * Cm(l, k);
					}
				}
			}
			float mdbs[2][3];
			evaluate_derivative_of_box_spline(cpset, suvw, mdbs);
			float nin = pow(2.0, lvm);
			for (int j = 0; j < 2; j++){
				for (int k = 0; k < 3; k++){
					if (Pk != 1) mdbs[j][k] = nin*mdbs[j][k];
					else mdbs[j][k] = -nin*mdbs[j][k];
				}
			}
			if (exN == 0){
				for (int i = 0; i < 2; i++){
					for (int j = 0; j < 3; j++){
						dbs[i][j] = mdbs[i][j];
					}
				}
			}
			else if (exN == 1){
				for (int i = 0; i < 3; i++){
					dbs[0][i] = -mdbs[1][i];
					dbs[1][i] = mdbs[0][i] - mdbs[1][i];
				}
			}
			else {
				for (int i = 0; i < 3; i++){
					dbs[0][i] = mdbs[1][i] - mdbs[0][i];
					dbs[1][i] = -mdbs[0][i];
				}
			}
			delete[] C0;
		}
		else {
			//境界エッジを含む場合 type3
			float(*C0)[3];
			C0 = new float[6][3];
			int valen = 2;
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 3; k++){
					C0[k][j] = Nodes[Faces[fn][(exN + k) % 3]][j];
				}
				for (int k = 0; k < 2; k++){
					C0[3 + k][j] = Nodes[neib[1][1 + k]][j];
				}
				C0[5][j] = Nodes[neib[2][2]][j];
			}
			if (neib[0][1] != -1 || neib[1][3] != -1 || neib[2][1] != -1)printf("boundary position is wrong @type3\n");

			//サンプル点の生成
			Eigen::MatrixXcf Sm(6, 6), Dmat;
			Eigen::VectorXcf mD(6);
			for (int j = 0; j < 6; j++){
				mD(j) = pow(Emap5->evals(j), lvm - 1);
			}
			Dmat = mD.asDiagonal();
			Sm = Emap5->evecs * Dmat*Emap5->evecs.inverse();

			Eigen::MatrixXf Cm(10, 3);
			for (int j = 0; j < 10; j++){
				for (int k = 0; k < 3; k++){
					Cm(j, k) = 0.0;
					for (int l = 0; l < 6; l++){
						for (int m = 0; m < 6; m++){
							Cm(j, k) += Emap5->Amat(j, l)*Sm(l, m).real()*C0[m][k];
						}
					}
				}
			}

			//Permutation matrixとスプライン基底の乗算
			float cpset[12][3];
			for (int j = 0; j < 12; j++){
				for (int k = 0; k < 3; k++){
					cpset[j][k] = 0.0;
					for (int l = 0; l < 10; l++){
						cpset[j][k] += (float)Emap5->Pmat[Pk][j][l] * Cm(l, k);
					}
				}
			}
			float mdbs[2][3];
			evaluate_derivative_of_box_spline(cpset, suvw, mdbs);
			float nin = pow(2.0, lvm);
			for (int j = 0; j < 2; j++){
				for (int k = 0; k < 3; k++){
					if (Pk != 1) mdbs[j][k] = nin*mdbs[j][k];
					else mdbs[j][k] = -nin*mdbs[j][k];
				}
			}
			if (exN == 0){
				for (int i = 0; i < 2; i++){
					for (int j = 0; j < 3; j++){
						dbs[i][j] = mdbs[i][j];
					}
				}
			}
			else if (exN == 1){
				for (int i = 0; i < 3; i++){
					dbs[0][i] = -mdbs[1][i];
					dbs[1][i] = mdbs[0][i] - mdbs[1][i];
				}
			}
			else {
				for (int i = 0; i < 3; i++){
					dbs[0][i] = mdbs[1][i] - mdbs[0][i];
					dbs[1][i] = -mdbs[0][i];
				}
			}
			delete[] C0;
		}
	}


	return true;

}


//制御点からbox-splineの計算
void loop_subdivision_mesh::evaluate_box_spline(float cpset[12][3], float uvw[3], float bs[3]){

	float basis[12];
	basis[0] = (pow(uvw[0], 4) + 2.0*pow(uvw[0], 3)*uvw[1]) / 12.0;
	basis[1] = (pow(uvw[0], 4) + 2.0*pow(uvw[0], 3)*uvw[2]) / 12.0;
	basis[5] = (pow(uvw[1], 4) + 2.0*pow(uvw[1], 3)*uvw[0]) / 12.0;
	basis[9] = (pow(uvw[1], 4) + 2.0*pow(uvw[1], 3)*uvw[2]) / 12.0;
	basis[8] = (pow(uvw[2], 4) + 2.0*pow(uvw[2], 3)*uvw[0]) / 12.0;
	basis[11] = (pow(uvw[2], 4) + 2.0*pow(uvw[2], 3)*uvw[1]) / 12.0;

	for (int i = 0; i < 3; i++){
		int bN = 2;
		if (i == 1) bN = 10;
		else if (i == 2) bN = 4;
		float p1 = uvw[i];
		float p2 = uvw[(i + 1) % 3];
		float p3 = uvw[(i + 2) % 3];
		basis[bN] = (pow(p1, 4) + 2.0*pow(p1, 3)*p3 + 6.0*pow(p1, 3)*p2 + 6.0*p1*p1*p2*p3 + 12.0*p1*p1*p2*p2 + 6.0*p1*p2*p2*p3 
			+ 6.0*p1*pow(p2, 3) + 2.0*pow(p2, 3)*p3 + pow(p2, 4)) / 12.0;

		bN = 3;
		if (i == 1) bN = 6;
		else if (i == 2) bN = 7;
		float pn0 = uvw[i];
		float opn[2];
		for (int j = 0; j < 2; j++){
			opn[j] = uvw[(i + 1 + j) % 3];
		}
		basis[bN] = (6.0*pow(pn0, 4) + 60.0*pn0 * pn0 * opn[0] * opn[1] + 12.0*opn[0] * opn[0] * opn[1] * opn[1]) / 12.0;
		for (int j = 0; j < 2; j++){
			basis[bN] += (36.0*pn0*opn[j] * opn[(j + 1) % 2] * opn[(j + 1) % 2] + 24.0*pn0*pn0*opn[j] * opn[j] + 24.0*pow(pn0, 3)*opn[j]
				+ 8.0*pn0 * pow(opn[j], 3) + 6.0*opn[(j + 1) % 2] * pow(opn[j], 3) + pow(opn[j], 4)) / 12.0;
		}
	}
	for (int i = 0; i < 3; i++){
		bs[i] = 0.0;
		for (int j = 0; j < 12; j++){
			bs[i] += cpset[j][i] * basis[j];
		}
	}


}


//制御点からbox-splineの微分の計算
void loop_subdivision_mesh::evaluate_derivative_of_box_spline(float cpset[12][3], float uvw[3], float dbs[2][3]){

	float basis[3][12];
	basis[0][0] = (4.0*pow(uvw[0], 3) + 6.0*uvw[0] * uvw[0] * uvw[1]) / 12.0;
	basis[0][1] = (4.0*pow(uvw[0], 3) + 6.0*uvw[0] * uvw[0] * uvw[2]) / 12.0;
	basis[0][5] = pow(uvw[1], 3) / 6.0;
	basis[0][9] = 0.0;
	basis[0][8] = pow(uvw[2], 3) / 6.0;
	basis[0][11] = 0.0;

	basis[1][0] = pow(uvw[0], 3) / 6.0;
	basis[1][1] = 0.0;
	basis[1][5] = (4.0*pow(uvw[1], 3) + 6.0*uvw[1] * uvw[1] * uvw[0]) / 12.0;
	basis[1][9] = (4.0*pow(uvw[1], 3) + 6.0*uvw[1] * uvw[1] * uvw[2]) / 12.0;
	basis[1][8] = 0.0;
	basis[1][11] = pow(uvw[2], 3) / 6.0;

	basis[2][0] = 0.0;
	basis[2][1] = pow(uvw[0], 3) / 6.0;
	basis[2][5] = 0.0;
	basis[2][9] = pow(uvw[1], 3) / 6.0;
	basis[2][8] = (4.0*pow(uvw[2], 3) + 6.0*uvw[2] * uvw[2] * uvw[0]) / 12.0;
	basis[2][11] = (4.0*pow(uvw[2], 3) + 6.0*uvw[2] * uvw[2] * uvw[1]) / 12.0;

	for (int i = 0; i < 3; i++){
		int bN = 2;
		if (i == 1) bN = 10;
		else if (i == 2) bN = 4;
		float p1 = uvw[i];
		float p2 = uvw[(i + 1) % 3];
		float p3 = uvw[(i + 2) % 3];
		basis[i][bN] = (4.0*pow(p1, 3) + 6.0*p1*p1*p3 + 18.0*p1*p1*p2 + 12.0*p1*p2*p3 + 24.0*p1*p2*p2 + 6.0*p2*p2*p3 + 6.0*pow(p2, 3)) / 12.0;
		basis[(i + 1) % 3][bN] = (6.0*pow(p1, 3) + 6.0*p1*p1*p3 + 24.0*p1*p1*p2 + 12.0*p1*p2*p3 + 18.0*p1*p2*p2 + 6.0*p2*p2*p3 + 4.0*pow(p2, 3)) / 12.0;
		basis[(i + 2) % 3][bN] = (2.0*pow(p1, 3) + 6.0*p1*p1*p2 + 6.0*p1*p2*p2 + 2.0*pow(p2, 3)) / 12.0;

		bN = 3;
		if (i == 1) bN = 6;
		else if (i == 2) bN = 7;
		float pn0 = uvw[i];
		float opn[2];
		for (int j = 0; j < 2; j++){
			opn[j] = uvw[(i + 1 + j) % 3];
		}
		basis[i][bN] = (24.0*pow(pn0, 3) + 120.0*pn0*opn[0] * opn[1]) / 12.0;
		basis[(i + 1) % 3][bN] = (60.0*pn0*pn0*opn[1] + 24.0*opn[0] * opn[1] * opn[1]) / 12.0;
		basis[(i + 2) % 3][bN] = (60.0*pn0*pn0*opn[0] + 24.0*opn[0] * opn[0] * opn[1]) / 12.0;
		for (int j = 0; j < 2; j++){
			basis[i][bN] += (36.0*opn[j] * opn[(j + 1) % 2] * opn[(j + 1) % 2] + 48.0*pn0*opn[j] * opn[j] + 72.0*pn0*pn0*opn[j]
				+ 8.0*pow(opn[j], 3)) / 12.0;
			basis[(i + 1 + j) % 3][bN] += (36.0*pn0*opn[(j + 1) % 2] * opn[(j + 1) % 2] + 48.0*pn0*pn0*opn[j] + 24.0*pow(pn0, 3)
				+ 24.0*pn0*opn[j] * opn[j] + 18.0*opn[(j + 1) % 2] * opn[j] * opn[j] + 4.0*pow(opn[j], 3)) / 12.0;
			basis[(i + 1 + (j + 1) % 2) % 3][bN] += (72.0*pn0*opn[j] * opn[(j + 1) % 2] + 6.0*pow(opn[j], 3)) / 12.0;
		}
	}


	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 3; j++){
			dbs[i][j] = 0.0;
			for (int k = 0; k < 12; k++){
				dbs[i][j] += cpset[k][j] * (basis[1 + i][k] - basis[0][k]);
			}
		}
	}


}


//近傍のリストからスプライン用の制御点を求める
void loop_subdivision_mesh::set_control_points_for_box_spline_with_reflection(int fn,int neib[3][6], float cpoint[12][3]){

	float neibC[3][6][3];
	for (int i = 0; i < 3; i++){
		int first_boundN = -1;
		for (int j = 0; j < 6; j++){
			if (neib[i][j] != -1){
				for (int k = 0; k < 3; k++){
					neibC[i][j][k] = Nodes[neib[i][j]][k];
				}
			}
			else if (first_boundN == -1) first_boundN = j;
		}
		if (first_boundN != -1){
			int fbN = first_boundN;
			for (int j = 0; j < 3; j++){
				neibC[i][fbN][j] = Nodes[Faces[fn][i]][j] + neibC[i][(fbN + 5) % 6][j] - neibC[i][(fbN + 4) % 6][j];
				neibC[i][(fbN + 1) % 6][j] = Nodes[Faces[fn][i]][j] + neibC[i][(fbN + 2) % 6][j] - neibC[i][(fbN + 3) % 6][j];
			}
		}
	}

	for (int i = 0; i < 3; i++){
		cpoint[0][i] = neibC[0][2][i];
		cpoint[1][i] = neibC[0][3][i];
		cpoint[2][i] = neibC[0][1][i];
		cpoint[3][i] = Nodes[Faces[fn][0]][i];
		cpoint[4][i] = neibC[0][4][i];
		cpoint[5][i] = neibC[1][3][i];
		cpoint[6][i] = Nodes[Faces[fn][1]][i];
		cpoint[7][i] = Nodes[Faces[fn][2]][i];
		cpoint[8][i] = neibC[2][2][i];
		cpoint[9][i] = neibC[1][2][i];
		cpoint[10][i] = neibC[1][1][i];
		cpoint[11][i] = neibC[2][3][i];
	}

}