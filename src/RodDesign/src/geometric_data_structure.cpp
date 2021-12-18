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

using namespace std;

#include "geometric_data_structure.h"

//////
//Meshメンバ関数
//////


Mesh::Mesh(){
	nodeN = 0;
	edgeN = 0;
	faceN = 0;
	CONSTRUCT_FLG = false;
}


Mesh::~Mesh(){


	if (CONSTRUCT_FLG == true){
		for (int i = 0; i < nodeN; i++){
			delete[] PFlist[i];
			delete[] PElist[i];
		}
		delete[] PFlist;
		delete[] PFsize;
		delete[] PElist;
		delete[] PEsize;

		delete[] EFlist;
		delete[] FElist;
	}
	if (nodeN != 0) delete[] Nodes;
	if (edgeN != 0) delete[] Edges;
	if(faceN != 0) delete[] Faces;


	nodeN = 0;
	edgeN = 0;
	faceN = 0;
	CONSTRUCT_FLG = false;
}


void Mesh::construct_data_structure(void){

	if (CONSTRUCT_FLG == true){
		return;
	}

	CONSTRUCT_FLG = true;

	//頂点->面リスト
	PFlist = new int*[nodeN];
	PFsize = new int[nodeN];
	for (int i = 0; i < nodeN; i++){
		PFsize[i] = 0;
	}
	for (int i = 0; i < faceN; i++){
		for (int j = 0; j < 3; j++){
			PFsize[Faces[i][j]]++;
		}
	}
	for (int i = 0; i < nodeN; i++){
		PFlist[i] = new int[PFsize[i]];
		PFsize[i] = 0;
	}
	for (int i = 0; i < faceN; i++){
		for (int j = 0; j < 3; j++){
			PFlist[Faces[i][j]][PFsize[Faces[i][j]]] = i;
			PFsize[Faces[i][j]]++;
		}
	}

	//頂点->エッジリスト
	PElist = new int*[nodeN];
	PEsize = new int[nodeN];
	for (int i = 0; i < nodeN; i++){
		PEsize[i] = 0;
	}
	for (int i = 0; i < edgeN; i++){
		PEsize[Edges[i][0]]++;
		PEsize[Edges[i][1]]++;
	}
	for (int i = 0; i < nodeN; i++){
		PElist[i] = new int[PEsize[i]];
		PEsize[i] = 0;
	}
	for (int i = 0; i < edgeN; i++){
		PElist[Edges[i][0]][PEsize[Edges[i][0]]] = i;
		PEsize[Edges[i][0]]++;
		PElist[Edges[i][1]][PEsize[Edges[i][1]]] = i;
		PEsize[Edges[i][1]]++;
	}

	
	//面->エッジ データ
	int damN = 0;
	FElist = new int[faceN][3];
	for (int i = 0; i < faceN; i++){
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < j; k++){
				int jn = Faces[i][j];
				int kn = Faces[i][k];
				int en = -1;
				for (int l = 0; l < PEsize[jn]; l++){
					if (Edges[PElist[jn][l]][0] == jn && Edges[PElist[jn][l]][1] == kn) en = PElist[jn][l];
					else if (Edges[PElist[jn][l]][1] == jn && Edges[PElist[jn][l]][0] == kn) en = PElist[jn][l];
					if (en != -1) break;
				}
				FElist[i][j + k - 1] = en;
				if (en == -1) printf("Okashii!!\n");
			}
		}
		int ehozon = FElist[i][1];
		FElist[i][1] = FElist[i][2];
		FElist[i][2] = ehozon;
	}


	//エッジ->面 データ
	EFlist = new int[edgeN][2];
	for (int i = 0; i < edgeN; i++){
		EFlist[i][0] = -1;
		EFlist[i][1] = -1;
	}
	for (int i = 0; i < faceN; i++){
		for (int j = 0; j < 3; j++){
			int p1 = Faces[i][j];
			int p2 = Faces[i][(j + 1) % 3];
			int en = FElist[i][j];
			if (p1 == Edges[en][0] && p2 == Edges[en][1]) {
				if (EFlist[en][0] == -1) EFlist[en][0] = i;
				else EFlist[en][1] = i;
			}
			else {
				if (EFlist[en][1] == -1) EFlist[en][1] = i;
				else EFlist[en][0] = i;
			}
		}
	}
	for (int i = 0; i < edgeN; i++){
		if (EFlist[i][0] == -1){
			EFlist[i][0] = EFlist[i][1];
			EFlist[i][1] = -1;

			int phozon = Edges[i][0];
			Edges[i][0] = Edges[i][1];
			Edges[i][1] = phozon;
		}
	}


	//PElist,PFlistの整理
	for (int i = 0; i < nodeN; i++){
		if (PEsize[i] == PFsize[i] || PEsize[i] == (PFsize[i] + 1)){
			int *PEhozon, *PFhozon;
			PEhozon = new int[PEsize[i]];
			PFhozon = new int[PFsize[i]];
			for (int j = 0; j < PEsize[i]; j++){
				PEhozon[j] = -1;
			}
			for (int j = 0; j < PFsize[i]; j++){
				PFhozon[j] = -1;
			}

			int be1 = -1;
			int be2 = -1;
			for (int j = 0; j < PEsize[i]; j++){
				int en = PElist[i][j];
				if (EFlist[en][1] == -1){
					if (Edges[en][0] == i && be1 == -1) be1 = en;
					else if (be2 == -1) be2 = en;
					else be1 = en;
				}
			}
			int oe = be1;
			if (be1 == -1) oe = PElist[i][0];
			int of = EFlist[oe][0];
			int fpn = -1;
			for (int j = 0; j < 3; j++){
				if (Faces[of][j] == i)fpn = j;
			}
			if (oe != FElist[of][fpn]){
				of = EFlist[oe][1];
			}

			int curN = 0;
			int cue = oe;
			int cuf = of;
			while (1){
				PEhozon[curN] = cue;
				PFhozon[curN] = cuf;

				fpn = -1;
				for (int j = 0; j < 3; j++){
					if (Faces[cuf][j] == i) fpn = j;
				}
				if (fpn != -1) {
					if (cue != FElist[cuf][(fpn + 2) % 3]) cue = FElist[cuf][(fpn + 2) % 3];
					else cue = FElist[cuf][fpn];
				}
				else {
					break;
				}
				if (EFlist[cue][0] != cuf) cuf = EFlist[cue][0];
				else cuf = EFlist[cue][1];
				if (cuf == -1){
					if (curN < (PEsize[i] - 1)) PEhozon[curN + 1] = cue;
					break;
				}
				if (cuf == of){
					break;
				}
				curN++;

				if (curN == PEsize[i] || curN == PFsize[i]) break;
			}

			if (PEhozon[PEsize[i] - 1] != -1 && PFhozon[PFsize[i] - 1] != -1){
				for (int j = 0; j < PEsize[i]; j++){
					PElist[i][j] = PEhozon[j];
				}
				for (int j = 0; j < PFsize[i]; j++){
					PFlist[i][j] = PFhozon[j];
				}
			}

			delete[] PEhozon;
			delete[] PFhozon;

		}
	}

}


void Mesh::construct_data_structure_without_edge(void){

	if (CONSTRUCT_FLG == true) return;

	CONSTRUCT_FLG = true;

	//頂点->面リスト
	PFlist = new int*[nodeN];
	PFsize = new int[nodeN];
	for (int i = 0; i < nodeN; i++){
		PFsize[i] = 0;
	}

	for (int i = 0; i < faceN; i++){
		for (int j = 0; j < 3; j++){
			PFsize[Faces[i][j]]++;
		}
	}

	for (int i = 0; i < nodeN; i++){
		PFlist[i] = new int[PFsize[i]];
		PFsize[i] = 0;
	}

	for (int i = 0; i < faceN; i++){
		for (int j = 0; j < 3; j++){
			int pn = Faces[i][j];
			PFlist[pn][PFsize[pn]] = i;
			PFsize[pn]++;
		}
	}

	//頂点->頂点 データ
	bool *bound;
	bound = new bool[nodeN];
	PElist = new int*[nodeN];
	PEsize = new int[nodeN];
	for (int i = 0; i < nodeN; i++){
		PEsize[i] = 0;
		PElist[i] = new int[PFsize[i]];
		bound[i] = false;
	}

	for (int i = 0; i < faceN; i++){
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < j; k++){
				int jn = Faces[i][j];
				int kn = Faces[i][k];
				bool ex_check = false;
				for (int l = 0; l < PEsize[jn]; l++){
					if (PElist[jn][l] == kn){
						ex_check = true;
					}
				}
				if (ex_check == false){
					if (PEsize[jn] == PFsize[jn]){
						int *PEhozon;
						PEhozon = new int[PFsize[jn]*2];
						for (int l = 0; l < PFsize[jn]; l++){
							PEhozon[l] = PElist[jn][l];
						}
						delete[] PElist[jn];
						PElist[jn] = PEhozon;
						bound[jn] = true;
					}
					PElist[jn][PEsize[jn]] = kn;
					PEsize[jn]++;
				}

				ex_check = false;
				for (int l = 0; l < PEsize[kn]; l++){
					if (PElist[kn][l] == jn){
						ex_check = true;
					}
				}
				if (ex_check == false){
					if (PEsize[kn] == PFsize[kn]){
						int *PEhozon;
						PEhozon = new int[PFsize[kn] * 2];
						for (int l = 0; l < PFsize[kn]; l++){
							PEhozon[l] = PElist[kn][l];
						}
						delete[] PElist[kn];
						PElist[kn] = PEhozon;
						bound[kn] = true;
					}
					PElist[kn][PEsize[kn]] = jn;
					PEsize[kn]++;
				}
			}
		}
	}

	for (int i = 0; i < nodeN; i++){
		if (bound[i] == true){
			int *PEhozon;
			PEhozon = new int[PEsize[i]];
			for (int j = 0; j < PEsize[i]; j++){
				PEhozon[j] = PElist[i][j];
			}
			delete[] PElist[i];
			PElist[i] = PEhozon;
		}
	}

	//頂点->エッジ、エッジ　データ
	edgeN = 0;
	for (int i = 0; i < nodeN; i++){
		edgeN += PEsize[i];
	}
	edgeN = (int)(edgeN / 2.0);

	Edges = new int[edgeN][2];
	edgeN = 0;
	for (int i = 0; i < nodeN; i++){
		for (int j = 0; j < PEsize[i]; j++){
			int jn = PElist[i][j];
			if (i < jn){
				Edges[edgeN][0] = i;
				Edges[edgeN][1] = jn;
				edgeN++;
			}
		}
	}

	for (int i = 0; i < nodeN; i++){
		PEsize[i] = 0;
	}
	for (int i = 0; i < edgeN; i++){
		int p1 = Edges[i][0];
		int p2 = Edges[i][1];
		PElist[p1][PEsize[p1]] = i;
		PElist[p2][PEsize[p2]] = i;
		PEsize[p1]++;
		PEsize[p2]++;
	}


	//面->エッジ データ
	FElist = new int[faceN][3];
	for (int i = 0; i < faceN; i++){
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < j; k++){
				int jn = Faces[i][j];
				int kn = Faces[i][k];
				int en = -1;
				for (int l = 0; l < PEsize[jn]; l++){
					if (Edges[PElist[jn][l]][0] == jn && Edges[PElist[jn][l]][1] == kn) en = PElist[jn][l];
					else if (Edges[PElist[jn][l]][1] == jn && Edges[PElist[jn][l]][0] == kn) en = PElist[jn][l];
					if (en != -1) break;
				}
				FElist[i][j + k - 1] = en;
				if (en == -1) printf("Okashii!!\n");
			}
		}
		int ehozon = FElist[i][1];
		FElist[i][1] = FElist[i][2];
		FElist[i][2] = ehozon;
	}


	//エッジ->面 データ
	EFlist = new int[edgeN][2];
	for (int i = 0; i < edgeN; i++){
		EFlist[i][0] = -1;
		EFlist[i][1] = -1;
	}
	for (int i = 0; i < faceN; i++){
		for (int j = 0; j < 3; j++){
			int p1 = Faces[i][j];
			int p2 = Faces[i][(j + 1) % 3];
			int en = FElist[i][j];
			if (p1 == Edges[en][0] && p2 == Edges[en][1]) {
				if (EFlist[en][0] == -1) EFlist[en][0] = i;
				else EFlist[en][1] = i;
			}
			else {
				if (EFlist[en][1] == -1) EFlist[en][1] = i;
				else EFlist[en][0] = i;
			}
		}
	}
	for (int i = 0; i < edgeN; i++){
		if (EFlist[i][0] == -1){
			EFlist[i][0] = EFlist[i][1];
			EFlist[i][1] = -1;

			int phozon = Edges[i][0];
			Edges[i][0] = Edges[i][1];
			Edges[i][1] = phozon;
		}
	}


	//PElist,PFlistの整理
	for (int i = 0; i < nodeN; i++){
		if (PEsize[i] == PFsize[i] || PEsize[i] == (PFsize[i] + 1)){
			int *PEhozon, *PFhozon;
			PEhozon = new int[PEsize[i]];
			PFhozon = new int[PFsize[i]];
			for (int j = 0; j < PEsize[i]; j++){
				PEhozon[j] = -1;
			}
			for (int j = 0; j < PFsize[i]; j++){
				PFhozon[j] = -1;
			}

			int be1 = -1;
			int be2 = -1;
			for (int j = 0; j < PEsize[i]; j++){
				int en = PElist[i][j];
				if (EFlist[en][1] == -1){
					if (Edges[en][0] == i && be1 == -1) be1 = en;
					else if (be2 == -1) be2 = en;
					else be1 = en;
				}
			}
			int oe = be1;
			if (be1 == -1) oe = PElist[i][0];
			int of = EFlist[oe][0];
			int fpn = -1;
			for (int j = 0; j < 3; j++){
				if (Faces[of][j] == i)fpn = j;
			}
			if (oe != FElist[of][fpn]){
				of = EFlist[oe][1];
			}

			int curN = 0;
			int cue = oe;
			int cuf = of;
			while (1){
				PEhozon[curN] = cue;
				PFhozon[curN] = cuf;

				fpn = -1;
				for (int j = 0; j < 3; j++){
					if (Faces[cuf][j] == i) fpn = j;
				}
				if (fpn != -1) {
					if (cue != FElist[cuf][(fpn + 2) % 3]) cue = FElist[cuf][(fpn + 2) % 3];
					else cue = FElist[cuf][fpn];
				}
				else {
					break;
				}
				if (EFlist[cue][0] != cuf) cuf = EFlist[cue][0];
				else cuf = EFlist[cue][1];
				if (cuf == -1){
					if (curN < (PEsize[i] - 1)) PEhozon[curN + 1] = cue;
					break;
				}
				if (cuf == of){
					break;
				}
				curN++;

				if (curN == PEsize[i] || curN == PFsize[i]) break;
			}

			if (PEhozon[PEsize[i] - 1] != -1 && PFhozon[PFsize[i] - 1] != -1){
				for (int j = 0; j < PEsize[i]; j++){
					PElist[i][j] = PEhozon[j];
				}
				for (int j = 0; j < PFsize[i]; j++){
					PFlist[i][j] = PFhozon[j];
				}
			}

			delete[] PEhozon;
			delete[] PFhozon;

		}
	}
}


void Mesh::offset_surface_mesh(float offv){


	//境界エッジ列の取得
	bool *esearch;
	esearch = new bool[edgeN];
	for (int i = 0; i < edgeN; i++){
		if (EFlist[i][1] == -1) esearch[i] = false;
		else esearch[i] = true;
	}


	std::vector<struct qelem *> eseq;

	std::vector<int > eseqN;
	struct qelem *sz;
	sz = (struct qelem*)malloc(sizeof(struct qelem));
	sz->ln = -1;
	sz->next = sz;
	while (1){

		int root_edge = -1;
		for (int i = 0; i < edgeN; i++){
			if (esearch[i] == false){
				root_edge = i;
				break;
			}
		}
		if (root_edge == -1) break;

		struct qelem *felem;
		felem = (struct qelem*) malloc(sizeof(struct qelem));
		felem->ln = root_edge;
		felem->next = sz;
		esearch[root_edge] = true;

		struct qelem *ceq;
		ceq = felem;
		int eleN = 1;
		while (1){
			int oen = ceq->ln;
			int ofn = EFlist[oen][0];

			int cen = oen;
			int cfn = ofn;
			while (1){
				int echeck;
				for (int i = 0; i < 3; i++){
					if (FElist[cfn][i] == cen) echeck = i;
				}
				cen = FElist[cfn][(echeck + 1) % 3];
				if (cfn != EFlist[cen][0]) cfn = EFlist[cen][0];
				else cfn = EFlist[cen][1];

				if (cfn == -1) break;
			}
			if (cfn == ofn) break;
			if (esearch[cen]) break;
			if (cen == root_edge){
				ceq->next = felem;
				break;
			}

			struct qelem *nelem;
			nelem = (struct qelem*)malloc(sizeof(struct qelem));
			nelem->ln = cen;
			nelem->next = sz;
			ceq->next = nelem;
			ceq = nelem;
			esearch[cen] = true;
			eleN++;
		}

		//逆方向の探索
		if (ceq->next != felem){
			ceq = felem;
			while (1){
				int oen = ceq->ln;
				int ofn = EFlist[oen][0];

				int cen = oen;
				int cfn = ofn;
				while (1){
					int echeck;
					for (int i = 0; i < 3; i++){
						if (FElist[cfn][i] == cen) echeck = i;
					}
					cen = FElist[cfn][(echeck + 2) % 3];
					if (cfn != EFlist[cen][0]) cfn = EFlist[cen][0];
					else cfn = EFlist[cen][1];

					if (cfn == -1) break;
				}
				if (cfn == ofn) break;
				if (esearch[cen]) break;

				struct qelem *nelem;
				nelem = (struct qelem*)malloc(sizeof(struct qelem));
				nelem->ln = cen;
				nelem->next = ceq;
				ceq = nelem;
				esearch[cen] = true;
				eleN++;
			}
			felem = ceq;
		}
		eseq.push_back(felem);
		eseqN.push_back(eleN);
	}
	delete[] esearch;


	
	int addFN = faceN;
	int addNN = nodeN;
	for (int i = 0; i < eseqN.size(); i++){
		addFN += 2 * eseqN[i];
	}


	float(*nhozon)[3];
	int (*fhozon)[3];
	nhozon = new float[nodeN + addNN][3];
	fhozon = new int[faceN + addFN][3];

	//頂点の追加
	for (int i = 0; i < nodeN; i++){
		for (int j = 0; j < 3; j++){
			nhozon[i][j] = Nodes[i][j];
		}
	}
	for (int i = 0; i < nodeN; i++){
		float odire[3];
		for (int j = 0; j < 3; j++){
			odire[j] = 0.0;
		}
		for (int j = 0; j < PFsize[i]; j++){
			int fn = PFlist[i][j];
			float ab[3], ac[3], nd[3];
			for (int k = 0; k < 3; k++){
				ab[k] = Nodes[Faces[fn][1]][k] - Nodes[Faces[fn][0]][k];
				ac[k] = Nodes[Faces[fn][2]][k] - Nodes[Faces[fn][0]][k];
			}
			for (int k = 0; k < 3; k++){
				nd[k] = ab[(k + 1) % 3] * ac[(k + 2) % 3] - ab[(k + 2) % 3] * ac[(k + 1) % 3];
			}
			float ndh = sqrt(nd[0] * nd[0] + nd[1] * nd[1] + nd[2] * nd[2]);
			if (ndh < 1.0e-5) ndh = 1.0e-5;
			for (int k = 0; k < 3; k++){
				odire[k] += nd[k] / ndh;
			}
		}
		float h = sqrt(odire[0] * odire[0] + odire[1] * odire[1] + odire[2] * odire[2]);
		if (h < 1.0e-5) h = 1.0e-5;
		for (int j = 0; j < 3; j++){
			nhozon[i + nodeN][j] = Nodes[i][j] - offv*odire[j] / h;
		}
	}

	//三角形の追加
	for (int i = 0; i < faceN; i++){
		for (int j =0 ; j < 3; j++){
			fhozon[i][j] = Faces[i][j];
		}
	}
	for (int i = 0; i < faceN; i++){
		fhozon[i + faceN][0] = Faces[i][0] + nodeN;
		fhozon[i + faceN][1] = Faces[i][2] + nodeN;
		fhozon[i + faceN][2] = Faces[i][1] + nodeN;
	}

	int addN = 2 * faceN;
	for (int i = 0; i < eseq.size(); i++){
		struct qelem *ceq;
		ceq = eseq[i];
		while (1){
			int cen = ceq->ln;
			int cfn = EFlist[cen][0];
			int opn;
			for (int j = 0; j < 3; j++){
				if (Faces[cfn][j] != Edges[cen][0] && Faces[cfn][j] != Edges[cen][1]) opn = j;
			}
			int p1 = Faces[cfn][(opn + 1) % 3];
			int p2 = Faces[cfn][(opn + 2) % 3];

			float ab[3], ac[3], ad[3], nd1[3], nd2[3];
			for (int j = 0; j < 3; j++){
				ab[j] = nhozon[p1 + nodeN][j] - nhozon[p1][j];
				ac[j] = nhozon[p2 + nodeN][j] - nhozon[p1][j];
				ad[j] = nhozon[p2][j] - nhozon[p1][j];
			}
			for (int j = 0; j < 3; j++){
				nd1[j] = ab[(j + 1) % 3] * ac[(j + 2) % 3] - ab[(j + 2) % 3] * ac[(j + 1) % 3];
				nd2[j] = ac[(j + 1) % 3] * ad[(j + 2) % 3] - ac[(j + 1) % 3] * ad[(j + 1) % 3];
			}
			float nd1h = sqrt(nd1[0] * nd1[0] + nd1[1] * nd1[1] + nd1[2] * nd1[2]);
			float nd2h = sqrt(nd2[0] * nd2[0] + nd2[1] * nd2[1] + nd2[2] * nd2[2]);
			if (nd1h < 1.0e-5) nd1h = 1.0e-5;
			if (nd2h < 1.0e-5) nd2h = 1.0e-5;
			for (int j = 0; j < 3; j++){
				nd1[j] /= nd1h;
				nd2[j] /= nd2h;
			}
			float seki1 = nd1[0] * nd2[0] + nd1[1] * nd2[1] + nd1[2] * nd2[2];

			for (int j = 0; j < 3; j++){
				ab[j] = nhozon[p1][j] - nhozon[p2][j];
				ac[j] = nhozon[p1 + nodeN][j] - nhozon[p2][j];
				ad[j] = nhozon[p2+nodeN][j] - nhozon[p1][j];
			}
			for (int j = 0; j < 3; j++){
				nd1[j] = ab[(j + 1) % 3] * ac[(j + 2) % 3] - ab[(j + 2) % 3] * ac[(j + 1) % 3];
				nd2[j] = ac[(j + 1) % 3] * ad[(j + 2) % 3] - ac[(j + 1) % 3] * ad[(j + 1) % 3];
			}
			nd1h = sqrt(nd1[0] * nd1[0] + nd1[1] * nd1[1] + nd1[2] * nd1[2]);
			nd2h = sqrt(nd2[0] * nd2[0] + nd2[1] * nd2[1] + nd2[2] * nd2[2]);
			if (nd1h < 1.0e-5) nd1h = 1.0e-5;
			if (nd2h < 1.0e-5) nd2h = 1.0e-5;
			for (int j = 0; j < 3; j++){
				nd1[j] /= nd1h;
				nd2[j] /= nd2h;
			}
			float seki2 = nd1[0] * nd2[0] + nd1[1] * nd2[1] + nd1[2] * nd2[2];

			if (seki1 >= seki2){
				fhozon[addN][0] = p1;
				fhozon[addN][1] = p1 + nodeN;
				fhozon[addN][2] = p2 + nodeN;
				fhozon[addN + 1][0] = p1;
				fhozon[addN + 1][1] = p2 + nodeN;
				fhozon[addN + 1][2] = p2;
			}
			else{
				fhozon[addN][0] = p1;
				fhozon[addN][1] = p1 + nodeN;
				fhozon[addN][2] = p2;
				fhozon[addN + 1][0] = p2;
				fhozon[addN + 1][1] = p1 + nodeN;
				fhozon[addN + 1][2] = p2 + nodeN;
			}
			addN += 2;

			if (ceq->next == eseq[i] || ceq->next == sz) break;
			ceq = ceq->next;
		}
	}

	delete[] Nodes,Faces,Edges;
	Nodes = nhozon;
	Faces = fhozon;

	if (CONSTRUCT_FLG){
		for (int i = 0; i < nodeN; i++){
			delete[] PFlist[i];
			delete[] PElist[i];
		}
		delete[] PFlist;
		delete[] PElist;
		delete[] PFsize;
		delete[] PEsize;
		delete[] EFlist;
		delete[] FElist;
		CONSTRUCT_FLG = false;
	}

	nodeN += addNN;
	faceN += addFN;


	//メモリの解放
	for (int i = 0; i < eseq.size(); i++){
		struct qelem *ceq;
		ceq = eseq[i]->next;
		while (1){
			if (ceq == sz || ceq == eseq[i]) break;
			struct qelem *ehozon;
			ehozon = ceq->next;
			free(ceq);
			ceq = ehozon;
		}
		free(eseq[i]);
	}
	free(sz);



}


void Mesh::mesh_copy(class Mesh *cmesh){

	nodeN = cmesh->nodeN;
	if (nodeN > 0){
		Nodes = new float[nodeN][3];
		for (int i = 0; i < nodeN; i++){
			for (int j = 0; j < 3; j++){
				Nodes[i][j] = cmesh->Nodes[i][j];
			}
		}
	}
	edgeN = cmesh->edgeN;
	if (edgeN > 0){
		Edges = new int[edgeN][2];
		for (int i = 0; i < edgeN; i++){
			Edges[i][0] = cmesh->Edges[i][0];
			Edges[i][1] = cmesh->Edges[i][1];
		}
	}
	faceN = cmesh->faceN;
	if (faceN > 0){
		Faces = new int[faceN][3];
		for (int i = 0; i < faceN; i++){
			for (int j = 0; j < 3; j++){
				Faces[i][j] = cmesh->Faces[i][j];
			}
		}
	}

	if (cmesh->CONSTRUCT_FLG){
		PFlist = new int*[nodeN];
		PElist = new int*[nodeN];
		PFsize = new int[nodeN];
		PEsize = new int[nodeN];
		for (int i = 0; i < nodeN; i++){
			PFsize[i] = cmesh->PFsize[i];
			PEsize[i] = cmesh->PEsize[i];
			if (PFsize[i] > 0) PFlist[i] = new int[PFsize[i]];
			if (PEsize[i] > 0) PElist[i] = new int[PEsize[i]];
			for (int j = 0; j < PFsize[i]; j++){
				PFlist[i][j] = cmesh->PFlist[i][j];
			}
			for (int j = 0; j < PEsize[i]; j++){
				PElist[i][j] = cmesh->PElist[i][j];
			}
		}
		EFlist = new int[edgeN][2];
		for (int i = 0; i < edgeN; i++){
			EFlist[i][0] = cmesh->EFlist[i][0];
			EFlist[i][1] = cmesh->EFlist[i][1];
		}
		FElist = new int[faceN][3];
		for (int i = 0; i < faceN; i++){
			for (int j = 0; j < 3; j++){
				FElist[i][j] = cmesh->FElist[i][j];
			}
		}
	}

}


//////
//その他基本データ処理関数
/////


//alias wavefront objファイルの読み込み
void read_obj_file(std::string fname, class Mesh *mesh){

	ifstream file;
	file.open(fname, ios::in);

	if (!file.is_open()){
		cout << "File open error!!" << endl;
	}

	float(*Nodes)[3];
	int(*Faces)[3];
	int nlimit = 500;
	int flimit = 500;
	Nodes = new float[nlimit][3];
	Faces = new int[flimit][3];

	mesh->nodeN = 0;
	mesh->faceN = 0;

	string dtype;
	while (file >> dtype && !file.eof()){
		if (dtype == "v"){
			float px, py, pz;
			file >> px >> py >> pz;
			Nodes[mesh->nodeN][0] = px;
			Nodes[mesh->nodeN][1] = py;
			Nodes[mesh->nodeN][2] = pz;
			mesh->nodeN++;

			if (mesh->nodeN == nlimit){
				float(*nhozon)[3];
				nhozon = new float[2 * nlimit][3];
				for (int i = 0; i < nlimit; i++){
					for (int j = 0; j < 3; j++){
						nhozon[i][j] = Nodes[i][j];
					}
				}
				delete[] Nodes;
				Nodes = nhozon;
				nlimit = 2 * nlimit;
			}
		}
		else if (dtype == "vn"){
			float vnx, vny, vnz;
			file >> vnx >> vny >> vnz;
		}
		else if (dtype == "vt"){
			float u, v;
			file >> u >> v;
		}
		else if (dtype == "f"){
			char cline[1024];
			file.getline(cline, 1024);
			istringstream istr(cline);
			int i0, i1, i2;

			istr >> i0;
			if (istr.get() == '/') {
				if (istr.peek() != '/') {
					int ti;
					istr >> ti;
				}
				if (istr.get() == '/') {
					int ni;
					istr >> ni;
				}
			}

			istr >> i1;
			if (istr.get() == '/') {
				if (istr.peek() != '/') {
					int ti;
					istr >> ti;
				}
				if (istr.get() == '/') {
					int ni;
					istr >> ni;
				}
			}

			while (!istr.eof()){
				istr >> i2;
				if (istr.get() == '/') {
					if (istr.peek() != '/') {
						int ti;
						istr >> ti;
					}
					if (istr.get() == '/') {
						int ni;
						istr >> ni;
					}
				}
			}
			Faces[mesh->faceN][0] = i0 - 1;
			Faces[mesh->faceN][1] = i1 - 1;
			Faces[mesh->faceN][2] = i2 - 1;
			mesh->faceN++;

			if (mesh->faceN == flimit){
				int(*fhozon)[3];
				fhozon = new int[flimit * 2][3];
				for (int i = 0; i < flimit; i++){
					for (int j = 0; j < 3; j++){
						fhozon[i][j] = Faces[i][j];
					}
				}
				delete[] Faces;
				Faces = fhozon;
				flimit = flimit * 2;
			}

		}
		else {
			string line;
			getline(file, line);
		}
	}
	file.close();

	mesh->Nodes = new float[mesh->nodeN][3];
	for (int i = 0; i < mesh->nodeN; i++){
		for (int j = 0; j < 3; j++){
			mesh->Nodes[i][j] = Nodes[i][j];
		}
	}

	mesh->Faces = new int[mesh->faceN][3];
	for (int i = 0; i < mesh->faceN; i++){
		for (int j = 0; j < 3; j++){
			mesh->Faces[i][j] = Faces[i][j];
		}
	}

	delete[] Nodes;
	delete[] Faces;
}


//libiglメッシュフォーマットの読み込み
void read_libigl_mesh_format(Eigen::MatrixXf V, Eigen::MatrixXi F, class Mesh *mesh){

	mesh->nodeN = V.rows();
	mesh->Nodes = new float[mesh->nodeN][3];
	for (int i = 0; i < mesh->nodeN; i++){
		for (int j = 0; j < 3; j++){
			mesh->Nodes[i][j] = V(i, j);
		}
	}

	mesh->faceN = F.rows();
	mesh->Faces = new int[mesh->faceN][3];
	for (int i = 0; i < mesh->faceN; i++){
		for (int j = 0; j < 3; j++){
			mesh->Faces[i][j] = F(i, j);
		}
	}

}


//alias wavefront objファイルの書き出し
void output_obj_file(std::string fname, class Mesh *mesh){

	FILE *fp = fopen(fname.c_str(), "w");
	if (fp == NULL){
		printf("File open error\n");
	}

	//fprintf(fp, "mtllib ./");
	//fprintf(fp,"%s.mtl\n",fname.c_str());

	fprintf(fp,"\n");
	for (int i = 0; i < mesh->nodeN; i++){
		fprintf(fp,"v ");
		fprintf(fp,"%f %f %f\n",mesh->Nodes[i][0],mesh->Nodes[i][1],mesh->Nodes[i][2]);
	}
	fprintf(fp,"# %d vertices, 0 vertices normals\n",mesh->nodeN);
	fprintf(fp,"\n");

	for (int i = 0; i < mesh->faceN; i++){
		fprintf(fp,"f %d//%d ",mesh->Faces[i][0]+1,mesh->Faces[i][0]+1);
		fprintf(fp,"%d//%d ",mesh->Faces[i][1]+1,mesh->Faces[i][1]+1);
		fprintf(fp, "%d//%d\n",mesh->Faces[i][2] + 1, mesh->Faces[i][2] + 1);
	}
	fprintf(fp,"# %d faces, 0 coords texture\n",mesh->faceN);
	fprintf(fp,"\n #End of File");
}


//二面角によるラベリング
int labeling_by_dihedral_angle(class Mesh *mesh, int *label, float ang_thres){

	struct qelem *selem, *zelem;

	zelem = (struct qelem*) malloc(sizeof(struct qelem));
	if (zelem == NULL){
		printf("Memory allocation error\n");
	}
	zelem->ln = -1;

	bool *fcheck;
	fcheck = new bool[mesh->faceN];
	for (int i = 0; i < mesh->faceN; i++){
		fcheck[i] = false;
	}

	int labelN = 0;
	while (1){
		int sfn = -1;
		for (int i = 0; i < mesh->faceN; i++){
			if (fcheck[i] == false){
				sfn = i;
				break;
			}
		}
		if (sfn == -1) break;

		selem = (struct qelem*) malloc(sizeof(struct qelem));
		if (selem == NULL){
			printf("Memory allocation error\n");
		}
		selem->next = zelem;
		selem->ln = sfn;
		label[sfn] = labelN;
		fcheck[sfn] = true;

		while (1){
			int fn = selem->ln;
			float ab[3], ac[3], fnd[3];
			for (int i = 0; i < 3; i++){
				ab[i] = mesh->Nodes[mesh->Faces[fn][1]][i] - mesh->Nodes[mesh->Faces[fn][0]][i];
				ac[i] = mesh->Nodes[mesh->Faces[fn][2]][i] - mesh->Nodes[mesh->Faces[fn][0]][i];
			}
			for (int i = 0; i < 3; i++){
				fnd[i] = ab[(i + 1) % 3] * ac[(i + 2) % 3] - ab[(i + 2) % 3] * ac[(i + 1) % 3];
			}
			float fndh = sqrt(fnd[0] * fnd[0] + fnd[1] * fnd[1] + fnd[2] * fnd[2]);
			if (fndh < 1.0e-5) fndh = 1.0e-5;
			for (int i = 0; i < 3; i++){
				fnd[i] /= fndh;
			}

			for (int i = 0; i < 3; i++){
				int en = mesh->FElist[fn][i];
				int fn2 = mesh->EFlist[en][0];
				if (fn2 == fn) fn2 = mesh->EFlist[en][1];
				if (fn2 != -1 && !fcheck[fn2]){
					for (int j = 0; j < 3; j++){
						ab[j] = mesh->Nodes[mesh->Faces[fn2][1]][j] - mesh->Nodes[mesh->Faces[fn2][0]][j];
						ac[j] = mesh->Nodes[mesh->Faces[fn2][2]][j] - mesh->Nodes[mesh->Faces[fn2][0]][j];
					}
					float fnd2[3];
					for (int j = 0; j < 3; j++){
						fnd2[j] = ab[(j + 1) % 3] * ac[(j + 2) % 3] - ab[(j + 2) % 3] * ac[(j + 1) % 3];
					}
					float fndh = sqrt(fnd2[0] * fnd2[0] + fnd2[1] * fnd2[1] + fnd2[2] * fnd2[2]);
					if (fndh < 1.0e-5) fndh = 1.0e-5;
					for (int j = 0; j < 3; j++){
						fnd2[j] /= fndh;
					}
					float seki = fnd[0] * fnd2[0] + fnd[1] * fnd2[1] + fnd[2] * fnd2[2];
					if (seki > cos(ang_thres)){
						struct qelem *nelem;
						nelem = (struct qelem *) malloc(sizeof(struct qelem));
						if (nelem == NULL){
							printf("Memory allocation error\n");
						}
						nelem->ln = fn2;
						nelem->next = selem->next;
						selem->next = nelem;
						label[fn2] = labelN;
						fcheck[fn2] = true;
					}
				}
			}

			struct qelem *qhozon;
			qhozon = selem->next;
			free(selem);

			if (qhozon == zelem) break;
			selem = qhozon;
		}


		labelN++;
	}
	free(zelem);

	return labelN;

}


//指定ラベル領域をalias wavefront obj形式で出力
void output_labeled_mesh_by_obj(std::string fname, class Mesh *mesh, int *label, int ln){
	
	int *nnode;
	nnode = new int[mesh->nodeN];
	for (int i = 0; i < mesh->nodeN; i++){
		nnode[i] = -1;
	}
	int nnodeN = 0;
	int nfaceN = 0;
	for (int i = 0; i < mesh->faceN; i++){
		if (label[i] == ln){
			nfaceN++;
			for (int j = 0; j < 3; j++){
				if (nnode[mesh->Faces[i][j]] == -1){
					nnode[mesh->Faces[i][j]] = 1;
					nnodeN++;
				}
			}
		}
	}
	nnodeN = 0;
	for (int i = 0; i < mesh->nodeN; i++){
		if (nnode[i] != -1){
			nnode[i] = nnodeN;
			nnodeN++;
		}
	}

	FILE *fp = fopen(fname.c_str(),"w");
	if (fp == NULL){
		printf("File open error\n");
	}

	fprintf(fp, "\n");
	for (int i = 0; i < mesh->nodeN; i++){
		if (nnode[i] != -1){
			fprintf(fp,"v ");
			fprintf(fp,"%f %f %f\n",mesh->Nodes[i][0],mesh->Nodes[i][1],mesh->Nodes[i][2]);
		}
	}
	fprintf(fp, "# %d vertices, 0 vertices normals\n", nnodeN);
	fprintf(fp, "\n");

	for (int i = 0; i < mesh->faceN; i++){
		if (label[i] == ln){
			fprintf(fp,"f %d//%d ",nnode[mesh->Faces[i][0]]+1,nnode[mesh->Faces[i][0]]+1);
			fprintf(fp,"%d//%d ",nnode[mesh->Faces[i][1]]+1,nnode[mesh->Faces[i][1]]+1);
			fprintf(fp,"%d//%d\n",nnode[mesh->Faces[i][2]]+1,nnode[mesh->Faces[i][2]]+1);
		}
	}
	fprintf(fp, "# %d faces, 0 coords texture\n", nfaceN);
	fprintf(fp, "\n #End of File");

}