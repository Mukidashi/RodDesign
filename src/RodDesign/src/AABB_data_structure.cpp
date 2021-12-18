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
#include <vector>
#include <cfloat>
#include <Eigen/Core>

using namespace std;

#include "geometric_data_structure.h"
#include "AABB_data_structure.h"
#include "basic_geometric_calculation.h"


////
//AABB_structクラスのメンバ関数
////

AABB_struct::AABB_struct(){
	STRUCT_EXIST = false;
}


AABB_struct::~AABB_struct(){
	
	if (STRUCT_EXIST){
		free_aabb_node(root);
		free(znode);
		STRUCT_EXIST = false;
	}
}


void AABB_struct::free_aabb_node(struct aabb_node *cnode){

	if (cnode->leaf != -1) {
		free(cnode);
		return;
	}

	for (int i = 0; i < 8; i++){
		if (cnode->snode[i] != znode){
			free_aabb_node(cnode->snode[i]);
		}
	}
	free(cnode);

}


void AABB_struct::build_AABB_struct(Mesh *mesh){

	if (STRUCT_EXIST){
		free_aabb_node(root);
		free(znode);
	}
	
	mymesh = mesh;

	znode = (struct aabb_node*)malloc(sizeof(struct aabb_node));
	if (znode == NULL){
		printf("Memory allocation error!\n");
	}

	root = (struct aabb_node*)malloc(sizeof(struct aabb_node));
	if (root == NULL){
		printf("Memory allocation error!\n");
	}
	for (int i = 0; i < 8; i++){
		root->snode[i] = znode;
	}
	root->leaf = -1;

	int(*dwlist)[2];
	dwlist = new int[mymesh->faceN][2];
	for (int i = 0; i < mymesh->faceN; i++){
		dwlist[i][0] = 0;
		dwlist[i][1] = 0;
	}

	recursive_mesh_divide(root, 0, 0, dwlist);

	delete[] dwlist;

	STRUCT_EXIST = true;
}


void AABB_struct::update_AABB_range(){

	if (!STRUCT_EXIST){
		printf("AABB structure does not exist\n");
		return;
	}

	
	float box[3][2];
	recursive_range_specification(root, box);
}


//最近点探索（最近点を含む三角形の座標と重心座標を返す）
bool AABB_struct::closest_point_search(float pc[3], float fps[3][3], float cord[3]){

	float eave = 0.0;
	for (int i = 0; i < mymesh->edgeN; i++){
		float ab[3];
		for (int j = 0; j < 3; j++){
			ab[j] = mymesh->Nodes[mymesh->Edges[i][0]][j] - mymesh->Nodes[mymesh->Edges[i][1]][j];
		}
		eave += sqrt(ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2]) /(float) mymesh->edgeN;
	}

	float crange = eave;
	float fullrange = 0.0;
	for (int i = 0; i < 3; i++){
		fullrange += pow(root->aabb[i][1] - root->aabb[i][0], 2);
	}
	fullrange = sqrt(fullrange);
	
	struct qelem *bottom;
	bottom = (struct qelem*)malloc(sizeof(struct qelem));
	if (bottom == NULL){
		printf("Memory allocation error\n");
	}
	bottom->ln = -1;
	struct qelem *ctop;
	int iterN = 0;
	while (1){
		float cbox[3][2];
		for (int i = 0; i<3; i++){
			cbox[i][0] = pc[i] - crange;
			cbox[i][1] = pc[i] + crange;
		}
		ctop = recursive_intersect_triangle_search(root, cbox, bottom);

		if (ctop != bottom) break;
		crange = 2 * crange;
		if (crange > fullrange) break;
		iterN++;
	}

	if (ctop == bottom){
		free(bottom);
		return false;
	}

	float mink;
	int minN = -1;
	struct qelem *celem;
	celem = bottom->next;
	while (1){
		float tfp[3][3], mycord[3];
		int fn = celem->ln;
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				tfp[j][k] = mymesh->Nodes[mymesh->Faces[fn][j]][k];
			}
		}
		float kyori = closest_point_of_triangle(pc, tfp, mycord);
		if (minN == -1){
			minN = fn;
			mink = kyori;
			for (int j = 0; j<3; j++){
				cord[j] = mycord[j];
			}
		}
		else if (mink > kyori){
			minN = fn;
			mink = kyori;
			for (int j = 0; j<3; j++){
				cord[j] = mycord[j];
			}
		}

		if (celem == ctop) break;
		celem = celem->next;
	}

	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			fps[i][j] = mymesh->Nodes[mymesh->Faces[minN][i]][j];
		}
	}

	free_qelem(bottom,ctop);

	//エッジ上に最近点がある場合の処理
	bool on_edge = false;
	bool on_vertex = false;
	int opn, ovn;
	for (int i = 0; i < 3; i++){
		if (fabs(cord[i]) < 1.0e-5 && cord[(i + 1) % 3] >= 0 && cord[(i + 2) % 3] >= 0){
			on_edge = true;
			opn = i;
		}
		if (fabs(cord[i]-1.0) < 1.0e-4 && fabs(cord[(i + 1) % 3]) < 1.0e-5 && fabs(cord[(i + 2) % 3]) < 1.0e-5){
			on_vertex = true;
			ovn = i;
		}
	}

	if (!on_vertex && on_edge){
		int fn = minN;
		int en = mymesh->FElist[fn][(opn+1)%3];
		int ofn = mymesh->EFlist[en][0];
		if (ofn == fn) ofn = mymesh->EFlist[en][1];

		if (ofn != -1){
			float edir[3],norm[3],fnp[3];
			float op1[3], op2[3], tp[3];
			int opn2;
			for (int i = 0; i < 3; i++){
				if (mymesh->Faces[ofn][i] != mymesh->Edges[en][0] && mymesh->Faces[ofn][i] != mymesh->Edges[en][1]) opn2 = i;
			}
			for (int i = 0; i < 3; i++){
				edir[i] = mymesh->Nodes[mymesh->Edges[en][1]][i] - mymesh->Nodes[mymesh->Edges[en][0]][i];
				op1[i] = mymesh->Nodes[mymesh->Faces[fn][opn]][i] - mymesh->Nodes[mymesh->Edges[en][0]][i];
				op2[i] = mymesh->Nodes[mymesh->Faces[ofn][opn2]][i] - mymesh->Nodes[mymesh->Edges[en][0]][i];
				tp[i] = cord[0] * fps[0][i] + cord[1] * fps[1][i] + cord[2] * fps[2][i];
				tp[i] = pc[i] - mymesh->Nodes[mymesh->Edges[en][0]][i];
				fnp[i] = (fps[1][(i + 1) % 3] - fps[0][(i + 1) % 3])*(fps[2][(i + 2) % 3] - fps[0][(i + 2) % 3]);
				fnp[i] -= (fps[1][(i + 2) % 3] - fps[0][(i + 2) % 3])*(fps[2][(i + 1) % 3] - fps[0][(i + 1) % 3]);
			}
			float edirh = sqrt(edir[0] * edir[0] + edir[1] * edir[1] + edir[2] * edir[2]);
			if (edirh < 1.0e-5) edirh = 1.0e-5;
			for (int i = 0; i < 3; i++){
				edir[i] /= edirh;
			}
			float seki1 = edir[0] * op1[0] + edir[1] * op1[1] + edir[2] * op1[2];
			for (int i = 0; i < 3; i++){
				op1[i] -= seki1 * edir[i];
			}
			seki1 = sqrt(op1[0] * op1[0] + op1[1] * op1[1] + op1[2] * op1[2]);
			if (seki1 < 1.0e-5) seki1 = 1.0e-5;
			for (int i = 0; i < 3; i++){
				op1[i] /= seki1;
			}
			for (int i = 0; i < 3; i++){
				norm[i] = edir[(i + 1) % 3] * op1[(i + 2) % 3] - edir[(i + 2) % 3] * op1[(i + 1) % 3];
			}
			float normh = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
			if (normh < 1.0e-5) normh = 1.0e-5;
			for (int i = 0; i < 3; i++){
				norm[i] /= normh;
			}
			seki1 = fnp[0] * norm[0] + fnp[1] * norm[1] + fnp[2] * norm[2];
			if (seki1 < 0.0){
				for (int i = 0; i < 3; i++){
					norm[i] = -norm[i];
				}
			}

			float ecod1[2], ecod2[2];
			float tph = sqrt(tp[0] * tp[0] + tp[1] * tp[1] + tp[2] * tp[2]);
			float op2h = sqrt(op2[0] * op2[0] + op2[1] * op2[1] + op2[2] * op2[2]);
			if (tph < 1.0e-5) tph = 1.0e-5;
			if (op2h < 1.0e-5) op2h = 1.0e-5;
			ecod1[0] = (op1[0] * tp[0] + op1[1] * tp[1] + op1[2] * tp[2]) / tph;
			ecod1[1] = (norm[0] * tp[0] + norm[1] * tp[1] + norm[2] * tp[2])/tph;
			ecod2[0] = (op1[0] * op2[0] + op1[1] * op2[1] + op1[2] * op2[2])/op2h;
			ecod2[1] = (norm[0] * op2[0] + norm[1] * op2[1] + norm[2] * op2[2])/op2h;
			if (ecod1[0] > 1.0) ecod1[0] = 1.0;
			else if (ecod1[0] < -1.0) ecod1[0] = -1.0;
			if (ecod2[0] > 1.0) ecod2[0] = 1.0;
			else if (ecod2[0] < -1.0) ecod2[0] = -1.0;
			ecod1[0] = acos(ecod1[0]);
			if (ecod1[1] < 0.0) ecod1[0] = -ecod1[0];
			ecod2[0] = acos(ecod2[0]);
			if (ecod2[1] < 0.0) ecod2[0] = -ecod2[0];

			bool pc_is_inside;
			if (ecod2[0] >= 0.0){
				if (ecod1[0] >= 0.0 && ecod1[0] <= ecod2[0]) pc_is_inside = true;
				else pc_is_inside = false;
			}
			else {
				if (ecod1[0] < 0.0 && ecod1[0] >= ecod2[0]) pc_is_inside = false;
				else pc_is_inside = true;
			}

			bool need_to_interchange = false;
			if (pc_is_inside && ecod1[0] < 0.0) need_to_interchange = true;
			else if (!pc_is_inside && ecod1[0] > 0.0) need_to_interchange = true;

			if (need_to_interchange){
				for (int i = 0; i < 3; i++){
					for (int j = 0; j < 3; j++){
						fps[i][j] = mymesh->Nodes[mymesh->Faces[ofn][i]][j];
					}
				}

				float chozon[3];
				chozon[opn2] = cord[opn];
				chozon[(opn2 + 2) % 3] = cord[(opn + 1) % 3];
				chozon[(opn2 + 1) % 3] = cord[(opn + 2) % 3];
				for (int i = 0; i < 3; i++){
					cord[i] = chozon[i];
				}
			}
		}

	}


	//最近点がメッシュ頂点と一致しているときの処理
	if (on_vertex){
		int fn = minN;
		int pn = mymesh->Faces[minN][ovn];
		float scord[2], scord2[2];
		for (int i = 0; i < 2; i++){
			scord[i] = 0.0;
			scord2[i] = 0.0;
		}
		for (int i = 0; i < mymesh->PFsize[pn]; i++){
			fn = mymesh->PFlist[pn][i];
			float fnd[3];
			float ab[3], ac[3];
			float mycord[2],mycord2[2];
			for (int j = 0; j < 3; j++){
				ab[j] = mymesh->Nodes[mymesh->Faces[fn][1]][j] - mymesh->Nodes[mymesh->Faces[fn][0]][j];
				ac[j] = mymesh->Nodes[mymesh->Faces[fn][2]][j] - mymesh->Nodes[mymesh->Faces[fn][0]][j];
			}
			for (int j = 0; j < 3; j++){
				fnd[j] = ab[(j + 1) % 3] * ac[(j + 2) % 3] - ab[(j + 2) % 3] * ab[(j + 1) % 3];
			}
			float fndh = sqrt(fnd[0] * fnd[0] + fnd[1] * fnd[1] + fnd[2] * fnd[2]);
			if (fndh < 1.0e-5) fndh = 1.0e-5;
			for (int j = 0; j < 3; j++){
				fnd[j] /= fndh;
			}
			if (fnd[0] > 1.0) fnd[0] = 1.0;
			if (fnd[0] < -1.0) fnd[0] = -1.0;
			mycord[0] = acos(fnd[0]);
			mycord2[0] = M_PI - acos(fnd[0]);
			if (fnd[1] < 0.0) mycord[0] = -mycord[0];
			if (fnd[1] >= 0.0) mycord2[0] = -mycord2[0];

			if (fnd[2] > 1.0) fnd[2] = 1.0;
			if (fnd[2] < -1.0) fnd[2] = -1.0;
			mycord[1] = acos(fnd[2]);
			mycord2[1] = M_PI - acos(fnd[2]);

			if (i == 0){
				for (int j = 0; j<2; j++){
					scord[j] = mycord[j];
					scord2[j] = mycord2[j];
				}
			}
			else {
				if (scord[0] > mycord[0] && (scord[0] - mycord[0]) > M_PI){
					scord[0] = (i*scord[0] + mycord[0] + 2 * M_PI) / ((float)(i + 1));
				}
				else if (scord[0] < mycord[0] && (mycord[0] - scord[0]) > M_PI){
					scord[0] = (i*scord[0] + mycord[0] - 2 * M_PI) / ((float)(i + 1));
				}
				else {
					scord[0] = (i*scord[0] + mycord[0]) / ((float)(i + 1));
				}
				if (scord[0] > M_PI) scord[0] -= 2 * M_PI;
				else if (scord[0] < -M_PI) scord[0] += 2 * M_PI;
				scord[1] = (i*scord[1] + mycord[1]) / ((float)(i + 1));

				if (scord2[0] > mycord2[0] && (scord2[0] - mycord2[0]) > M_PI){
					scord2[0] = (i*scord2[0] + mycord2[0] + 2 * M_PI) / ((float)(i + 1));
				}
				else if (scord2[0] < mycord2[0] && (mycord2[0] - scord2[0]) > M_PI){
					scord2[0] = (i*scord2[0] + mycord2[0] - 2 * M_PI) / ((float)i + 1);
				}
				else {
					scord2[0] = (i*scord2[0] + mycord2[0]) / ((float)(i + 1));
				}
				if (scord2[0] > M_PI) scord2[0] -= 2 * M_PI;
				else if (scord2[0] < -M_PI) scord2[0] += 2 * M_PI;
				scord2[1] = (i*scord2[1] + mycord2[1]) / ((float)(i + 1));
			}
		}
		float sloc[3], sloc2[3];
		sloc[0] = sin(scord[1])*cos(scord[0]);
		sloc[1] = sin(scord[1])*sin(scord[0]);
		sloc[2] = cos(scord[1]);
		sloc2[0] = sin(scord2[1])*cos(scord2[0]);
		sloc2[1] = sin(scord2[1])*sin(scord2[0]);
		sloc2[2] = cos(scord2[1]);
		
		float tpn[3];
		for (int i = 0; i < 3; i++){
			tpn[i] = pc[i] - mymesh->Nodes[pn][i];
		}
		float tpnh = sqrt(tpn[0] * tpn[0] + tpn[1] * tpn[1] + tpn[2] * tpn[2]);
		if (tpnh < 1.0e-5) tpnh = 1.0e-5;
		for (int i = 0; i < 3; i++){
			tpn[i] /= tpnh;
		}
		float kyori1 = sqrt(pow(tpn[0] - sloc[0], 2) + pow(tpn[1] - sloc[1], 2) + pow(tpn[2] - sloc[2], 2));
		float kyori2 = sqrt(pow(tpn[0] - sloc2[0], 2) + pow(tpn[1] - sloc2[1], 2) + pow(tpn[2] - sloc2[2], 2));

		bool is_outside = true;
		if (kyori1 > kyori2){
			is_outside = false;
		} 

		int nfn = -1;
		for (int i = 0; i < mymesh->PFsize[pn]; i++){
			fn = mymesh->PFlist[pn][i];
			float fnd[3];
			float ab[3], ac[3];
			for (int j = 0; j < 3; j++){
				ab[j] = mymesh->Nodes[mymesh->Faces[fn][1]][j] - mymesh->Nodes[mymesh->Faces[fn][0]][j];
				ac[j] = mymesh->Nodes[mymesh->Faces[fn][2]][j] - mymesh->Nodes[mymesh->Faces[fn][0]][j];
			}
			for (int j = 0; j < 3; j++){
				fnd[j] = ab[(j + 1) % 3] * ac[(j + 2) % 3] - ab[(j + 2) % 3] * ab[(j + 1) % 3];
			}
			float fndh = sqrt(fnd[0] * fnd[0] + fnd[1] * fnd[1] + fnd[2] * fnd[2]);
			if (fndh < 1.0e-5) fndh = 1.0e-5;
			for (int j = 0; j < 3; j++){
				fnd[j] /= fndh;
			}
			float seki = fnd[0] * tpn[0] + fnd[1] * tpn[1] + fnd[2] * tpn[2];
			if (is_outside && seki >= 0.0){
				nfn = fn;
				break;
			}
			else if (!is_outside && seki <= 0.0){
				nfn = fn;
				break;
			}
		}
		if (nfn == -1){
			printf("suitable face does not exist\n");
			return false;
		}
		else {
			int pnn = 0;
			for (int i = 0; i < 3; i++){
				if (pn == mymesh->Faces[nfn][i]) pnn = i;
				cord[i] = 0.0;
				for (int j = 0; j < 3; j++){
					fps[i][j] = mymesh->Nodes[mymesh->Faces[nfn][i]][j];
				}
			}
			cord[pnn] = 1.0;
		}

	}


	return true;
}


//最近点探索（最近点を含む三角形のインデックスと重心座標を返す）
bool AABB_struct::closest_point_search_with_mesh_info(float pc[3], int *msfn, float cord[3]){

	float eave = 0.0;
	for (int i = 0; i < mymesh->edgeN; i++){
		float ab[3];
		for (int j = 0; j < 3; j++){
			ab[j] = mymesh->Nodes[mymesh->Edges[i][0]][j] - mymesh->Nodes[mymesh->Edges[i][1]][j];
		}
		eave += sqrt(ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2]) / (float)mymesh->edgeN;
	}

	float crange = eave;
	float fullrange = 0.0;
	for (int i = 0; i < 3; i++){
		fullrange += pow(root->aabb[i][1] - root->aabb[i][0], 2);
	}
	fullrange = sqrt(fullrange);

	struct qelem *bottom;
	bottom = (struct qelem*)malloc(sizeof(struct qelem));
	if (bottom == NULL){
		printf("Memory allocation error\n");
	}
	bottom->ln = -1;
	struct qelem *ctop;
	int iterN = 0;
	while (1){
		float cbox[3][2];
		for (int i = 0; i<3; i++){
			cbox[i][0] = pc[i] - crange;
			cbox[i][1] = pc[i] + crange;
		}
		ctop = recursive_intersect_triangle_search(root, cbox, bottom);

		if (ctop != bottom) break;
		crange = 2 * crange;
		if (crange > fullrange) break;
		iterN++;
	}

	if (ctop == bottom){
		free(bottom);
		return false;
	}

	float mink;
	int minN = -1;
	struct qelem *celem;
	celem = bottom->next;
	while (1){
		float tfp[3][3], mycord[3];
		int fn = celem->ln;
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				tfp[j][k] = mymesh->Nodes[mymesh->Faces[fn][j]][k];
			}
		}
		float kyori = closest_point_of_triangle(pc, tfp, mycord);
		if (minN == -1){
			minN = fn;
			mink = kyori;
			for (int j = 0; j<3; j++){
				cord[j] = mycord[j];
			}
		}
		else if (mink > kyori){
			minN = fn;
			mink = kyori;
			for (int j = 0; j<3; j++){
				cord[j] = mycord[j];
			}
		}

		if (celem == ctop) break;
		celem = celem->next;
	}

	*msfn = minN;

	free_qelem(bottom, ctop);

	//エッジ上に最近点がある場合の処理
	bool on_edge = false;
	bool on_vertex = false;
	int opn, ovn;
	for (int i = 0; i < 3; i++){
		if (fabs(cord[i]) < 1.0e-5 && cord[(i + 1) % 3] >= 0 && cord[(i + 2) % 3] >= 0){
			on_edge = true;
			opn = i;
		}
		if (fabs(cord[i] - 1.0) < 1.0e-4 && fabs(cord[(i + 1) % 3]) < 1.0e-5 && fabs(cord[(i + 2) % 3]) < 1.0e-5){
			on_vertex = true;
			ovn = i;
		}
	}

	if (!on_vertex && on_edge){
		int fn = minN;
		int en = mymesh->FElist[fn][(opn + 1) % 3];
		int ofn = mymesh->EFlist[en][0];
		if (ofn == fn) ofn = mymesh->EFlist[en][1];

		if (ofn != -1){
			float edir[3], norm[3], fnp[3],mfp[3][3];
			float op1[3], op2[3], tp[3];
			int opn2;
			for (int i = 0; i < 3; i++){
				if (mymesh->Faces[ofn][i] != mymesh->Edges[en][0] && mymesh->Faces[ofn][i] != mymesh->Edges[en][1]) opn2 = i;
				for (int j = 0; j < 3; j++){
					mfp[i][j] = mymesh->Nodes[mymesh->Faces[minN][i]][j];
				}
			}
			for (int i = 0; i < 3; i++){
				edir[i] = mymesh->Nodes[mymesh->Edges[en][1]][i] - mymesh->Nodes[mymesh->Edges[en][0]][i];
				op1[i] = mymesh->Nodes[mymesh->Faces[fn][opn]][i] - mymesh->Nodes[mymesh->Edges[en][0]][i];
				op2[i] = mymesh->Nodes[mymesh->Faces[ofn][opn2]][i] - mymesh->Nodes[mymesh->Edges[en][0]][i];
				tp[i] = pc[i] - mymesh->Nodes[mymesh->Edges[en][0]][i];
				fnp[i] = (mfp[1][(i + 1) % 3] - mfp[0][(i + 1) % 3])*(mfp[2][(i+2)%3] - mfp[0][(i+2)%3]);
				fnp[i] -= (mfp[1][(i + 2) % 3] - mfp[0][(i + 2) % 3])*(mfp[2][(i + 1) % 3] - mfp[0][(i + 1) % 3]);
			}
			float edirh = sqrt(edir[0] * edir[0] + edir[1] * edir[1] + edir[2] * edir[2]);
			if (edirh < 1.0e-5) edirh = 1.0e-5;
			for (int i = 0; i < 3; i++){
				edir[i] /= edirh;
			}
			float seki1 = edir[0] * op1[0] + edir[1] * op1[1] + edir[2] * op1[2];
			for (int i = 0; i < 3; i++){
				op1[i] -= seki1 * edir[i];
			}
			seki1 = sqrt(op1[0] * op1[0] + op1[1] * op1[1] + op1[2] * op1[2]);
			if (seki1 < 1.0e-5) seki1 = 1.0e-5;
			for (int i = 0; i < 3; i++){
				op1[i] /= seki1;
			}
			for (int i = 0; i < 3; i++){
				norm[i] = edir[(i + 1) % 3] * op1[(i + 2) % 3] - edir[(i + 2) % 3] * op1[(i + 1) % 3];
			}
			float normh = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
			if (normh < 1.0e-5) normh = 1.0e-5;
			for (int i = 0; i < 3; i++){
				norm[i] /= normh;
			}
			seki1 = fnp[0] * norm[0] + fnp[1] * norm[1] + fnp[2] * norm[2];
			if (seki1 < 0.0){
				for (int i = 0; i < 3; i++){
					norm[i] = -norm[i];
				}
			}

			float ecod1[2], ecod2[2];
			float tph = sqrt(tp[0] * tp[0] + tp[1] * tp[1] + tp[2] * tp[2]);
			float op2h = sqrt(op2[0] * op2[0] + op2[1] * op2[1] + op2[2] * op2[2]);
			if (tph < 1.0e-5) tph = 1.0e-5;
			if (op2h < 1.0e-5) op2h = 1.0e-5;
			ecod1[0] = (op1[0] * tp[0] + op1[1] * tp[1] + op1[2] * tp[2]) / tph;
			ecod1[1] = (norm[0] * tp[0] + norm[1] * tp[1] + norm[2] * tp[2]) / tph;
			ecod2[0] = (op1[0] * op2[0] + op1[1] * op2[1] + op1[2] * op2[2]) / op2h;
			ecod2[1] = (norm[0] * op2[0] + norm[1] * op2[1] + norm[2] * op2[2]) / op2h;
			if (ecod1[0] > 1.0) ecod1[0] = 1.0;
			else if (ecod1[0] < -1.0) ecod1[0] = -1.0;
			if (ecod2[0] > 1.0) ecod2[0] = 1.0;
			else if (ecod2[0] < -1.0) ecod2[0] = -1.0;
			ecod1[0] = acos(ecod1[0]);
			if (ecod1[1] < 0.0) ecod1[0] = -ecod1[0];
			ecod2[0] = acos(ecod2[0]);
			if (ecod2[1] < 0.0) ecod2[0] = -ecod2[0];

			bool pc_is_inside;
			if (ecod2[0] >= 0.0){
				if (ecod1[0] >= 0.0 && ecod1[0] <= ecod2[0]) pc_is_inside = true;
				else pc_is_inside = false;
			}
			else {
				if (ecod1[0] < 0.0 && ecod1[0] >= ecod2[0]) pc_is_inside = false;
				else pc_is_inside = true;
			}

			bool need_to_interchange = false;
			if (pc_is_inside && ecod1[0] < 0.0) need_to_interchange = true;
			else if (!pc_is_inside && ecod1[0] > 0.0) need_to_interchange = true;

			if (need_to_interchange){
				*msfn = ofn;

				float chozon[3];
				chozon[opn2] = cord[opn];
				chozon[(opn2 + 2) % 3] = cord[(opn + 1) % 3];
				chozon[(opn2 + 1) % 3] = cord[(opn + 2) % 3];
				for (int i = 0; i < 3; i++){
					cord[i] = chozon[i];
				}
			}
		}

	}


	//最近点がメッシュ頂点と一致しているときの処理
	if (on_vertex){
		int fn = minN;
		int pn = mymesh->Faces[minN][ovn];
		float scord[2], scord2[2];
		for (int i = 0; i < 2; i++){
			scord[i] = 0.0;
			scord2[i] = 0.0;
		}
		for (int i = 0; i < mymesh->PFsize[pn]; i++){
			fn = mymesh->PFlist[pn][i];
			float fnd[3];
			float ab[3], ac[3];
			float mycord[2], mycord2[2];
			for (int j = 0; j < 3; j++){
				ab[j] = mymesh->Nodes[mymesh->Faces[fn][1]][j] - mymesh->Nodes[mymesh->Faces[fn][0]][j];
				ac[j] = mymesh->Nodes[mymesh->Faces[fn][2]][j] - mymesh->Nodes[mymesh->Faces[fn][0]][j];
			}
			for (int j = 0; j < 3; j++){
				fnd[j] = ab[(j + 1) % 3] * ac[(j + 2) % 3] - ab[(j + 2) % 3] * ab[(j + 1) % 3];
			}
			float fndh = sqrt(fnd[0] * fnd[0] + fnd[1] * fnd[1] + fnd[2] * fnd[2]);
			if (fndh < 1.0e-5) fndh = 1.0e-5;
			for (int j = 0; j < 3; j++){
				fnd[j] /= fndh;
			}
			if (fnd[0] > 1.0) fnd[0] = 1.0;
			if (fnd[0] < -1.0) fnd[0] = -1.0;
			mycord[0] = acos(fnd[0]);
			mycord2[0] = M_PI - acos(fnd[0]);
			if (fnd[1] < 0.0) mycord[0] = -mycord[0];
			if (fnd[1] >= 0.0) mycord2[0] = -mycord2[0];

			if (fnd[2] > 1.0) fnd[2] = 1.0;
			if (fnd[2] < -1.0) fnd[2] = -1.0;
			mycord[1] = acos(fnd[2]);
			mycord2[1] = M_PI - acos(fnd[2]);

			if (i == 0){
				for (int j = 0; j<2; j++){
					scord[j] = mycord[j];
					scord2[j] = mycord2[j];
				}
			}
			else {
				if (scord[0] > mycord[0] && (scord[0] - mycord[0]) > M_PI){
					scord[0] = (i*scord[0] + mycord[0] + 2 * M_PI) / ((float)(i + 1));
				}
				else if (scord[0] < mycord[0] && (mycord[0] - scord[0]) > M_PI){
					scord[0] = (i*scord[0] + mycord[0] - 2 * M_PI) / ((float)(i + 1));
				}
				else {
					scord[0] = (i*scord[0] + mycord[0]) / ((float)(i + 1));
				}
				if (scord[0] > M_PI) scord[0] -= 2 * M_PI;
				else if (scord[0] < -M_PI) scord[0] += 2 * M_PI;
				scord[1] = (i*scord[1] + mycord[1]) / ((float)(i + 1));

				if (scord2[0] > mycord2[0] && (scord2[0] - mycord2[0]) > M_PI){
					scord2[0] = (i*scord2[0] + mycord2[0] + 2 * M_PI) / ((float)(i + 1));
				}
				else if (scord2[0] < mycord2[0] && (mycord2[0] - scord2[0]) > M_PI){
					scord2[0] = (i*scord2[0] + mycord2[0] - 2 * M_PI) / ((float)i + 1);
				}
				else {
					scord2[0] = (i*scord2[0] + mycord2[0]) / ((float)(i + 1));
				}
				if (scord2[0] > M_PI) scord2[0] -= 2 * M_PI;
				else if (scord2[0] < -M_PI) scord2[0] += 2 * M_PI;
				scord2[1] = (i*scord2[1] + mycord2[1]) / ((float)(i + 1));
			}
		}
		float sloc[3], sloc2[3];
		sloc[0] = sin(scord[1])*cos(scord[0]);
		sloc[1] = sin(scord[1])*sin(scord[0]);
		sloc[2] = cos(scord[1]);
		sloc2[0] = sin(scord2[1])*cos(scord2[0]);
		sloc2[1] = sin(scord2[1])*sin(scord2[0]);
		sloc2[2] = cos(scord2[1]);

		float tpn[3];
		for (int i = 0; i < 3; i++){
			tpn[i] = pc[i] - mymesh->Nodes[pn][i];
		}
		float tpnh = sqrt(tpn[0] * tpn[0] + tpn[1] * tpn[1] + tpn[2] * tpn[2]);
		if (tpnh < 1.0e-5) tpnh = 1.0e-5;
		for (int i = 0; i < 3; i++){
			tpn[i] /= tpnh;
		}
		float kyori1 = sqrt(pow(tpn[0] - sloc[0], 2) + pow(tpn[1] - sloc[1], 2) + pow(tpn[2] - sloc[2], 2));
		float kyori2 = sqrt(pow(tpn[0] - sloc2[0], 2) + pow(tpn[1] - sloc2[1], 2) + pow(tpn[2] - sloc2[2], 2));

		bool is_outside = true;
		if (kyori1 > kyori2){
			is_outside = false;
		}

		int nfn = -1;
		for (int i = 0; i < mymesh->PFsize[pn]; i++){
			fn = mymesh->PFlist[pn][i];
			float fnd[3];
			float ab[3], ac[3];
			for (int j = 0; j < 3; j++){
				ab[j] = mymesh->Nodes[mymesh->Faces[fn][1]][j] - mymesh->Nodes[mymesh->Faces[fn][0]][j];
				ac[j] = mymesh->Nodes[mymesh->Faces[fn][2]][j] - mymesh->Nodes[mymesh->Faces[fn][0]][j];
			}
			for (int j = 0; j < 3; j++){
				fnd[j] = ab[(j + 1) % 3] * ac[(j + 2) % 3] - ab[(j + 2) % 3] * ab[(j + 1) % 3];
			}
			float fndh = sqrt(fnd[0] * fnd[0] + fnd[1] * fnd[1] + fnd[2] * fnd[2]);
			if (fndh < 1.0e-5) fndh = 1.0e-5;
			for (int j = 0; j < 3; j++){
				fnd[j] /= fndh;
			}
			float seki = fnd[0] * tpn[0] + fnd[1] * tpn[1] + fnd[2] * tpn[2];
			if (is_outside && seki >= 0.0){
				nfn = fn;
				break;
			}
			else if (!is_outside && seki <= 0.0){
				nfn = fn;
				break;
			}
		}
		if (nfn == -1){
			printf("suitable face does not exist\n");
		}
		else {
			int pnn = 0;
			for (int i = 0; i < 3; i++){
				if (pn == mymesh->Faces[nfn][i]) pnn = i;
				cord[i] = 0.0;
				*msfn = nfn;
			}
			cord[pnn] = 1.0;
		}

	}


	return true;





}


//直線と交差する三角形の探索
bool AABB_struct::search_intersect_triangle_with_line(float orig[3], float dire[3], std::vector<int> *intF){


	struct qelem *bottom, *top;
	bottom = (struct qelem*) malloc(sizeof(struct qelem));
	if (bottom == NULL){
		printf("Memory allocation error\n");
	}
	bottom->ln = -1;

	top = recursive_search_intersect_triangle_with_line(root, orig, dire, bottom);

	if (top == bottom){
		free(bottom);
		return true;
	}
	struct qelem *cele;
	cele = bottom->next;
	while (1){
		int fn = cele->ln;
		float fps[3][3], cord[3];
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				fps[i][j] = mymesh->Nodes[mymesh->Faces[fn][i]][j];
			}
		}
		bool inter_check = intersection_of_triangle_with_line(orig, dire, fps, cord);

		if (inter_check && cord[0] >= 0.0 && cord[1] >= 0.0 && cord[2] >= 0.0) intF->push_back(fn);

		if (cele == top) break;
		struct qelem *hozon;
		hozon = cele->next;
		free(cele);
		cele = hozon;
	}
	free(top);
	free(bottom);

	return true;
}


//AABB構築の際に再帰的に呼び出す関数
void AABB_struct::recursive_mesh_divide(aabb_node *cnode, int depth, int wnum, int(*dwlist)[2]){

	float gp[3];
	for (int i = 0; i < 3; i++){
		gp[i] = 0.0;
	}

	int setN = 0;
	for (int i = 0; i < mymesh->faceN; i++){
		if (dwlist[i][0] == depth && dwlist[i][1] == wnum){
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 3; k++){
					gp[k] += mymesh->Nodes[mymesh->Faces[i][j]][k] / 3.0;
				}
			}
			setN++;
		}
	}

	if (setN > 0){
		for (int i = 0; i < 3; i++){
			gp[i] /= (float)setN;
		}
	}
	else {
		printf("Okashii@recursive_mesh_divide\n");
		return;
	}

	if (setN == 1){
		for (int i = 0; i < mymesh->faceN; i++){
			if (dwlist[i][0] == depth && dwlist[i][1] == wnum){
				cnode->leaf = i;
				dwlist[i][0] = -1;
				dwlist[i][1] = -1;
				break;
			}
		}
		return;
	}

	int cset[8];
	for (int i = 0; i < 8; i++){
		cset[i] = 0;
	}
	for (int i = 0; i < mymesh->faceN; i++){
		if (dwlist[i][0] == depth && dwlist[i][1] == wnum){
			float fp[3];
			for (int j = 0; j < 3; j++){
				fp[j] = 0.0;
				for (int k = 0; k < 3; k++){
					fp[j] += mymesh->Nodes[mymesh->Faces[i][k]][j] / 3.0;
				}
			}

			dwlist[i][0] = depth + 1;
			if (fp[2] <= gp[2]){
				if (fp[1] <= gp[1]){
					if (fp[0] <= gp[0]){
						dwlist[i][1] = 0;
						cset[0]++;
					}
					else {
						dwlist[i][1] = 1;
						cset[1]++;
					}
				}
				else {
					if (fp[0] <= gp[0]){
						dwlist[i][1] = 2;
						cset[2]++;
					}
					else {
						dwlist[i][1] = 3;
						cset[3]++;
					}
				}
			}
			else{
				if (fp[1] <= gp[1]){
					if (fp[0] <= gp[0]){
						dwlist[i][1] = 4;
						cset[4]++;
					}
					else {
						dwlist[i][1] = 5;
						cset[5]++;
					}
				}
				else {
					if (fp[0] <= gp[0]){
						dwlist[i][1] = 6;
						cset[6]++;
					}
					else {
						dwlist[i][1] = 7;
						cset[7]++;
					}
				}
			}
		}
	}

	//三角形重心が一致した場合の対策
	int all_in_one = false;
	for (int i = 0; i < 8; i++){
		if (cset[i] == setN) all_in_one = true;
	}
	if (all_in_one){
		for (int i = 0; i < 8; i++){
			cset[i] = 0;
		}
		int snum = 0;
		for (int i = 0; i < mymesh->faceN; i++){
			if (dwlist[i][0] == (depth+1)){
				dwlist[i][1] = snum;
				cset[snum]++;
				snum = (snum + 1) % 8;
			}
		}
	}


	for (int i = 0; i < 8; i++){
		if (cset[i] > 0){
			struct aabb_node *nnode;
			nnode = (struct aabb_node*)malloc(sizeof(struct aabb_node));
			if (nnode == NULL){
				printf("Memory allocation error!\n");
			}
			for (int j = 0; j < 8; j++){
				nnode->snode[j] = znode;
			}
			nnode->leaf = -1;
			cnode->snode[i] = nnode;
			recursive_mesh_divide(nnode, depth + 1, i, dwlist);
		}
	}

}


//AABBの値の更新の際に再帰的に呼び出す関数
void AABB_struct::recursive_range_specification(struct aabb_node *cnode, float cbox[3][2]){

	if (cnode->leaf != -1){
		int fn = cnode->leaf;
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				if (j == 0){
					cbox[i][0] = mymesh->Nodes[mymesh->Faces[fn][j]][i];
					cbox[i][1] = mymesh->Nodes[mymesh->Faces[fn][j]][i];
				}
				else {
					if (cbox[i][0] > mymesh->Nodes[mymesh->Faces[fn][j]][i]) cbox[i][0] = mymesh->Nodes[mymesh->Faces[fn][j]][i];
					if (cbox[i][1] < mymesh->Nodes[mymesh->Faces[fn][j]][i]) cbox[i][1] = mymesh->Nodes[mymesh->Faces[fn][j]][i];
				}
			}
		}
		for (int i = 0; i < 3; i++){
			cnode->aabb[i][0] = cbox[i][0];
			cnode->aabb[i][1] = cbox[i][1];
		}
		return;
	}

	bool fupdate = false;
	for (int i = 0; i < 8; i++){
		if (cnode->snode[i] != znode){
			float nbox[3][2];
			recursive_range_specification(cnode->snode[i], nbox);
			if (!fupdate){
				for (int j = 0; j < 3; j++){
					cbox[j][0] = nbox[j][0];
					cbox[j][1] = nbox[j][1];
				}
				fupdate = true;
			}
			else {
				for (int j = 0; j < 3; j++){
					if (cbox[j][0] > nbox[j][0]) cbox[j][0] = nbox[j][0];
					if (cbox[j][1] < nbox[j][1]) cbox[j][1] = nbox[j][1];
				}
			}
		}
	}
	if (fupdate){
		for (int i = 0; i < 3; i++){
			cnode->aabb[i][0] = cbox[i][0];
			cnode->aabb[i][1] = cbox[i][1];
		}
	}

}


//最近点探索の際に再帰的に呼び出す関数
struct qelem* AABB_struct::recursive_intersect_triangle_search(struct aabb_node *cnode, float cbox[3][2], struct qelem *ctop){
	bool intersect_check = true;
	for (int i = 0; i < 3; i++){
		if (cbox[i][0] > cnode->aabb[i][1] || cbox[i][1] < cnode->aabb[i][0]) intersect_check = false;
	}
	if (!intersect_check) return ctop;

	if (cnode->leaf != -1){
		struct qelem *ntop;
		ntop = (struct qelem*)malloc(sizeof(struct qelem));
		if (ntop == NULL){
			printf("Memory allocation error\n");
		}
		ntop->ln = cnode->leaf;
		ctop->next = ntop;
		//printf("%d\n",ntop->ln);
		return ntop;
	}

	struct qelem *ntop,*thozon;
	ntop = ctop;
	for (int i = 0; i < 8; i++){
		if (cnode->snode[i] != znode){
			thozon = recursive_intersect_triangle_search(cnode->snode[i], cbox, ntop);
			//if (thozon != ntop) printf("->%d\n",thozon->ln);
			ntop = thozon;
		}
	}
	return ntop;
}


//直線と交差する三角形探索のために再帰的に呼び出す関数
struct qelem* AABB_struct::recursive_search_intersect_triangle_with_line(struct aabb_node *cnode, float orig[3], float dire[3], struct qelem *ctop){

	//AABBと直線の交差の判定
	bool intersect_flg = false;
	float k1, k2;
	float pip1[2], pip2[2];
	if (fabs(dire[0]) >= FLT_MIN){
		k1 = (cnode->aabb[0][0] - orig[0]) / dire[0];
		k2 = (cnode->aabb[0][1] - orig[0]) / dire[0];
		pip1[0] = orig[1] + k1*dire[1];
		pip1[1] = orig[2] + k1*dire[2];
		pip2[0] = orig[1] + k2*dire[1];
		pip2[1] = orig[2] + k2*dire[2];
		if (pip1[0] >= cnode->aabb[1][0] && pip1[0] <= cnode->aabb[1][1] && pip1[1] >= cnode->aabb[2][0] && pip1[1] <= cnode->aabb[2][1]) intersect_flg = true;
		if (pip2[0] >= cnode->aabb[1][0] && pip2[0] <= cnode->aabb[1][1] && pip2[1] >= cnode->aabb[2][0] && pip2[1] <= cnode->aabb[2][1]) intersect_flg = true;
	}
	if (!intersect_flg && fabs(dire[1]) >= FLT_MIN){
		k1 = (cnode->aabb[1][0] - orig[1]) / dire[1];
		k2 = (cnode->aabb[1][1] - orig[1]) / dire[1];
		pip1[0] = orig[0] + k1*dire[0];
		pip1[1] = orig[2] + k1*dire[2];
		pip2[0] = orig[0] + k2*dire[0];
		pip2[1] = orig[2] + k2*dire[2];
		if (pip1[0] >= cnode->aabb[0][0] && pip1[0] <= cnode->aabb[0][1] && pip1[1] >= cnode->aabb[2][0] && pip1[1] <= cnode->aabb[2][1]) intersect_flg = true;
		if (pip2[0] >= cnode->aabb[0][0] && pip2[0] <= cnode->aabb[0][1] && pip2[1] >= cnode->aabb[2][0] && pip2[1] <= cnode->aabb[2][1]) intersect_flg = true;
	}
	if (!intersect_flg && fabs(dire[2]) >= FLT_MIN){
		k1 = (cnode->aabb[2][0] - orig[2]) / dire[2];
		k2 = (cnode->aabb[2][1] - orig[2]) / dire[2];
		pip1[0] = orig[0] + k1 * dire[0];
		pip1[1] = orig[1] + k1*dire[1];
		pip2[0] = orig[0] + k2*dire[0];
		pip2[1] = orig[1] + k2*dire[1];
		if (pip1[0] >= cnode->aabb[0][0] && pip1[0] <= cnode->aabb[0][1] && pip1[1] >= cnode->aabb[1][0] && pip1[1] <= cnode->aabb[1][1]) intersect_flg = true;
		if (pip2[0] >= cnode->aabb[0][0] && pip2[0] <= cnode->aabb[0][1] && pip2[1] >= cnode->aabb[1][0] && pip2[1] <= cnode->aabb[1][1]) intersect_flg = true;
	}

	if (!intersect_flg) return ctop;

	if (cnode->leaf != -1){
		struct qelem *ntop;
		ntop = (struct qelem*)malloc(sizeof(struct qelem));
		if (ntop == NULL){
			printf("Memory allocation error\n");
		}
		ntop->ln = cnode->leaf;
		ctop->next = ntop;
		return ntop;
	}

	struct qelem *ntop, *thozon;
	ntop = ctop;
	for (int i = 0; i < 8; i++){
		if (cnode->snode[i] != znode){
			thozon = recursive_search_intersect_triangle_with_line(cnode->snode[i], orig, dire, ntop);
			ntop = thozon;
		}
	}
	return ntop;

}



////
//oriented_basic_AABB_structクラスのメンバ関数
////


oriented_basic_AABB_struct::oriented_basic_AABB_struct(){
	

	STRUCT_EXIST = false;
	MAX_DEPTH = 16;

}


oriented_basic_AABB_struct::~oriented_basic_AABB_struct(){

	if (STRUCT_EXIST){
		free_basic_aabb_node(root);
		free(znode);
		STRUCT_EXIST = false;
	}
}


void oriented_basic_AABB_struct::free_basic_aabb_node(struct basic_aabb_node *cnode){

	if (cnode->LEAF_FLG){
		free(cnode);
		return;
	}

	for (int i = 0; i < 8; i++){
		free_basic_aabb_node(cnode->dnode[i]);
	}
	free(znode);

}


//AABB構造の構築
void oriented_basic_AABB_struct::build_AABB_struct(float size, float orient[3][3],float orig[3],int depth){

	if (STRUCT_EXIST){
		free_basic_aabb_node(root);
		free(znode);
	}

	if (depth > MAX_DEPTH) depth = MAX_DEPTH;

	//初期化
	FULL_SIZE = size;
	for (int i = 0; i < 3; i++){
		Origin[i] = orig[i];
	}
	float axh = sqrt(pow(orient[0][0], 2) + pow(orient[0][1], 2) + pow(orient[0][2], 2));
	if (axh < FLT_MIN) {
		printf("orientation1 is too short\n");
		axh = FLT_MIN;
	}
	for (int i = 0; i < 3; i++){
		Axis[0][i] = orient[0][i] / axh;
	}
	float seki = Axis[0][0] * orient[1][0] + Axis[0][1] * orient[1][1] + Axis[0][2] * orient[1][2];
	for (int i = 0; i < 3; i++){
		Axis[1][i] = orient[1][i] - seki*Axis[0][i];
	}
	axh = sqrt(pow(Axis[1][0], 2) + pow(Axis[1][1], 2) + pow(Axis[1][2], 2));
	if (axh < FLT_MIN){
		printf("orientation2 is parallel to orientation1\n");
		axh = FLT_MIN;
	}
	for (int i = 0; i < 3; i++){
		Axis[1][i] /= axh;
	}
	for (int i = 0; i < 3; i++){
		Axis[2][i] = Axis[0][(i + 1) % 3] * Axis[1][(i + 2) % 3] - Axis[0][(i + 2) % 3] * Axis[1][(i + 1) % 3];
	}
	axh = sqrt(pow(Axis[2][0], 2) + pow(Axis[2][1], 2) + pow(Axis[2][2], 2));
	for (int i = 0; i < 3; i++){
		Axis[2][i] /= axh;
	}


	znode = (struct basic_aabb_node*) malloc(sizeof(struct basic_aabb_node));
	if (znode == NULL){
		printf("Memory allocation error\n");
	}
	znode->LEAF_FLG = false;

	root = (struct basic_aabb_node*)malloc(sizeof(struct basic_aabb_node));
	if (root == NULL){
		printf("Memory allocation error\n");
	}
	root->LEAF_FLG = false;
	for (int i = 0; i < 8; i++){
		root->dnode[i] = znode;
	}
	for (int i = 0; i < 3; i++){
		root->aabb[i][0] = -size / 2.0;
		root->aabb[i][1] = size / 2.0;
	}
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j++){
			for (int k = 0; k < 2; k++){
				int snn = 4 * i + 2 * j + k;
				for (int l = 0; l < 3; l++){
					root->sumip[snn][l] = Origin[l];
					if (k == 0) root->sumip[snn][l] += root->aabb[0][0] * Axis[0][l];
					else root->sumip[snn][l] += root->aabb[0][1] * Axis[0][l];
					if (j == 0) root->sumip[snn][l] += root->aabb[1][0] * Axis[1][l];
					else root->sumip[snn][l] += root->aabb[1][1] * Axis[1][l];
					if (i == 0) root->sumip[snn][l] += root->aabb[2][0] * Axis[2][l];
					else root->sumip[snn][l] += root->aabb[2][1] * Axis[2][l];
				}
			}
		}
	}

	recursive_space_division(root, depth, 0);


}


//メッシュの挿入
void oriented_basic_AABB_struct::insert_mesh_elemets(class Mesh *mesh){

	mymesh = mesh;

	struct qelem *ftop,*cbottom;
	ftop = (struct qelem*)malloc(sizeof(struct qelem));
	if (ftop == NULL){
		printf("Memory allocation error\n");
	}
	ftop->ln = -1;
	cbottom = ftop;
	for (int i = 0; i < mymesh->faceN; i++){
		float fps[3][3];
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				fps[j][k] = 0.0;
				for (int l = 0; l < 3; l++){
					fps[j][k] += (mymesh->Nodes[mymesh->Faces[i][j]][l] - Origin[l])*Axis[k][l];
				}
			}
		}
		bool is_inside = intersection_check_of_aabb_and_triagle(root->aabb, fps);

		if (is_inside){
			struct qelem *nef;
			nef = (struct qelem*)malloc(sizeof(struct qelem));
			if (nef == NULL){
				printf("Memory allocation error\n");
			}
			nef->ln = i;
			nef->next = cbottom;
			cbottom = nef;
		}
	}

	recursive_assign_elements(root, cbottom);

	while (1){
		if (cbottom->ln == -1) break;
		struct qelem *chozon;
		chozon = cbottom->next;
		free(cbottom);
		cbottom = chozon;
	}
	free(ftop);
}


//直線上にあるメッシュ内部ノードの特定
void oriented_basic_AABB_struct::fill_intermediate_node_along_line(float orig[3], float dire[3]){

	struct qelem_basic_aabb *top,*bottom;
	bottom = (struct qelem_basic_aabb*)malloc(sizeof(struct qelem_basic_aabb));
	if (bottom == NULL){
		printf("Memory allocation error\n");
	}
	bottom->mynode = znode;

	top = recursive_search_intersect_node_with_line(root, orig, dire, bottom);


	//セルの並び換え
	std::vector<int> intF;
	std::vector<float> intPos;
	struct qelem_basic_aabb *ctop,*nlist;
	nlist = (struct qelem_basic_aabb*)malloc(sizeof(struct qelem_basic_aabb));
	if (nlist == NULL){
		printf("Memory allocation error\n");
	}
	nlist->next = bottom;
	ctop = top;
	while (1){
		if (ctop == bottom) break;
		
		//ノード内部三角形の登録
		for (int i = 0; i < ctop->mynode->fmemN; i++){
			int fn = ctop->mynode->fmemb[i];
			intF.push_back(fn);
		}

		//AABBノードの並べ替え
		struct qelem_basic_aabb *nhozon,*nelem;
		nhozon = nlist;
		while (1){
			if (nhozon->next->mynode == znode) break;
			if (nhozon->next->mynode->aabb[0][0] >= ctop->mynode->aabb[0][0]) break;
			nhozon = nhozon->next;
		}
		nelem = (struct qelem_basic_aabb*)malloc(sizeof(struct qelem_basic_aabb));
		if (nelem == NULL){
			printf("Memory allocation errorn\n");
		}
		nelem->mynode = ctop->mynode;
		nelem->next = nhozon->next;
		nhozon->next = nelem;


		ctop = ctop->next;
	}

	if (intF.size() > 0){
		//重複三角形の除去と交点の計算
		int *fset;
		int fnum = 0;
		fset = new int[intF.size()];
		for (int i = 0; i < intF.size(); i++){
			bool dupli_check = false;
			for (int j = 0; j < fnum; j++){
				if (fset[j] == intF[i]) dupli_check = true;
			}
			if (!dupli_check){
				fset[fnum] = intF[i];
				fnum++;
			}
		}
		for (int i = 0; i < fnum; i++){
			int fn = fset[i];
			float fps[3][3], cord[3];
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 3; k++){
					fps[j][k] = 0.0;
					for (int l = 0; l < 3; l++){
						fps[j][k] += (mymesh->Nodes[mymesh->Faces[fn][j]][l] - Origin[l])*Axis[k][l];
					}
				}
			}
			bool is_intersect = intersection_of_triangle_with_line(orig, dire, fps, cord);

			if (is_intersect && cord[0] >= 0.0 && cord[1] >= 0.0 && cord[2] >= 0.0){
				float tp = cord[0] * fps[0][0] + cord[1] * fps[1][0] + cord[2] * fps[2][0];
				intPos.push_back(tp);
			}
		}
		delete[] fset;

		//内部にあるAABBの特定
		if (intPos.size() > 0 && intPos.size() % 2 == 0){
			int int_num = intPos.size();
			float *Pos;
			Pos = new float[int_num];
			for (int i = 0; i < int_num; i++){
				float phozon = intPos[i];
				for (int j = 0; j < i; j++){
					if (phozon < Pos[j]){
						float thozon = phozon;
						phozon = Pos[j];
						Pos[j] = thozon;
					}
				}
				Pos[i] = phozon;
			}

			int addN = 0;
			for (int i = 0; i < int_num / 2; i++){
				int inn1, inn2;
				inn1 = -1;
				inn2 = -1;
				struct qelem_basic_aabb *nhozon;
				nhozon = nlist->next;
				int iterN = 0;
				while (1){
					if (nhozon == bottom) break;
					if (nhozon->mynode->aabb[0][0] <= Pos[2 * i] && nhozon->mynode->aabb[0][1] >= Pos[2 * i]) inn1 = iterN;
					if (nhozon->mynode->aabb[0][0] <= Pos[2 * i + 1] && nhozon->mynode->aabb[0][1] >= Pos[2 * i + 1]) inn2 = iterN;
					iterN++;
					nhozon = nhozon->next;
				}
				nhozon = nlist->next;
				iterN = 0;
				while (1){
					if (nhozon == bottom) break;
					if (iterN >= inn1 && iterN <= inn2) nhozon->mynode->is_full = true;
					iterN++;
					nhozon = nhozon->next;
				}
			}
		}
		else if(intPos.size() > 0){
			printf("odd intersect\n");
		}
	}

	


	//メモリの解放
	struct qelem_basic_aabb *qhozon;
	qhozon = top;
	while (1){
		if (qhozon == bottom) break;
		struct qelem_basic_aabb *nele;
		nele = qhozon->next;
		free(qhozon);
		qhozon = nele;
	}
	qhozon = nlist;
	while (1){
		if (qhozon == bottom) break;
		struct qelem_basic_aabb *nele;
		nele = qhozon->next;
		free(qhozon);
		qhozon = nele;
	}
	free(bottom);


}


//AABB構築の際に再帰的に呼ばれる関数
void oriented_basic_AABB_struct::recursive_space_division(struct basic_aabb_node *cnode, int max_depth, int depth){

	if (depth == max_depth){
		cnode->LEAF_FLG = true;
		cnode->is_full = false;
		cnode->fmemN = 0;
		return;
	}

	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j++){
			for (int k = 0; k < 2; k++){
				int dnn = i * 4 + j * 2 + k;
				struct basic_aabb_node *nnode;
				nnode = (struct basic_aabb_node*)malloc(sizeof(struct basic_aabb_node));
				if (nnode == NULL){
					printf("Memory allocation error\n");
				}
				nnode->LEAF_FLG = false;
				for (int l = 0; l < 8; l++){
					nnode->dnode[l] = znode;
				}
				if (k == 0){
					nnode->aabb[0][0] = cnode->aabb[0][0];
					nnode->aabb[0][1] = (cnode->aabb[0][0] + cnode->aabb[0][1])/ 2.0;
				}
				else {
					nnode->aabb[0][0] = (cnode->aabb[0][0] + cnode->aabb[0][1]) / 2.0;
					nnode->aabb[0][1] = cnode->aabb[0][1];
				}
				if (j == 0){
					nnode->aabb[1][0] = cnode->aabb[1][0];
					nnode->aabb[1][1] = (cnode->aabb[1][0] + cnode->aabb[1][1]) / 2.0;
				}
				else {
					nnode->aabb[1][0] = (cnode->aabb[1][0] + cnode->aabb[1][1]) / 2.0;
					nnode->aabb[1][1] = cnode->aabb[1][1];
				}
				if (i == 0){
					nnode->aabb[2][0] = cnode->aabb[2][0];
					nnode->aabb[2][1] = (cnode->aabb[2][0] + cnode->aabb[2][1]) / 2.0;
				}
				else {
					nnode->aabb[2][0] = (cnode->aabb[2][0] + cnode->aabb[2][1]) / 2.0;
					nnode->aabb[2][1] = cnode->aabb[2][1];
				}
				for (int l = 0; l < 2; l++){
					for (int m = 0; m < 2; m++){
						for (int n = 0; n < 2; n++){
							int snn = 4 * l + 2 * m + n;
							for (int o = 0; o < 3; o++){
								nnode->sumip[snn][o] = Origin[o];
								if (n == 0) nnode->sumip[snn][o] += nnode->aabb[0][0] * Axis[0][o];
								else nnode->sumip[snn][o] += nnode->aabb[0][1] * Axis[0][o];
								if (m == 0) nnode->sumip[snn][o] += nnode->aabb[1][0] * Axis[1][o];
								else nnode->sumip[snn][o] += nnode->aabb[1][1] * Axis[1][o];
								if (l == 0) nnode->sumip[snn][o] += nnode->aabb[2][0] * Axis[2][o];
								else nnode->sumip[snn][o] += nnode->aabb[2][1] * Axis[2][o];
							}
						}
					}
				}

				cnode->dnode[dnn] = nnode;
				recursive_space_division(nnode, max_depth, depth + 1);
			}
		}
	}

}


//メッシュ挿入の際に再帰的に呼ばれる関数
void oriented_basic_AABB_struct::recursive_assign_elements(struct basic_aabb_node *cnode, struct qelem *celem){

	if (cnode->LEAF_FLG){
		struct qelem *cf;
		cf = celem;
		int fmemN = 0;
		while (1){
			if (cf->ln == -1) break;
			fmemN++;
			cf = cf->next;
		}
		if (fmemN > 0) {
			cnode->fmemN = fmemN;
			cnode->fmemb = new int[fmemN];
		}
		cf = celem;
		fmemN = 0;
		while (1){
			if (cf->ln == -1) break;
			cnode->fmemb[fmemN] = cf->ln;
			fmemN++;
			cf = cf->next;
		}

		if (fmemN != 0) cnode->is_full = true;
		return;
	}

	struct qelem **ftops,**fbottoms;
	ftops = new struct qelem*[8];
	fbottoms = new struct qelem*[8];
	for (int i = 0; i < 8; i++){
		ftops[i] = (struct qelem*)malloc(sizeof(struct qelem));
		if (ftops[i] == NULL){
			printf("Memory allocation error\n");
		}
		ftops[i]->ln = -1;
		fbottoms[i] = ftops[i];
	}

	struct qelem *cfele;
	cfele = celem;
	while (1){
		if (cfele->ln == -1) break;
		int fn = cfele->ln;
		float fps[3][3];
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				fps[j][k] = 0.0;
				for (int l = 0; l < 3; l++){
					fps[j][k] += (mymesh->Nodes[mymesh->Faces[fn][j]][l] - Origin[l])*Axis[k][l];
				}
			}
		}
		for (int j = 0; j < 8; j++){
			bool is_inside = intersection_check_of_aabb_and_triagle(cnode->dnode[j]->aabb, fps);
			if (is_inside){
				struct qelem *nfe;
				nfe = (struct qelem*)malloc(sizeof(struct qelem));
				if (nfe == NULL){
					printf("Memory allocation error\n");
				}
				nfe->ln = fn;
				nfe->next = fbottoms[j];
				fbottoms[j] = nfe;
			}
		}
		cfele = cfele->next;
	}
	for (int i = 0; i < 8; i++){
		if (fbottoms[i]->ln != -1) recursive_assign_elements(cnode->dnode[i], fbottoms[i]);

		while (1){
			if (fbottoms[i]->ln == -1) break;
			struct qelem *fhozon;
			fhozon = fbottoms[i]->next;
			free(fbottoms[i]);
			fbottoms[i] = fhozon;
		}
		free(ftops[i]);
	}
}


//直線と交差するノード探索のために再帰的に呼ばれる関数
struct qelem_basic_aabb* oriented_basic_AABB_struct::recursive_search_intersect_node_with_line(struct basic_aabb_node *cnode, float orig[3], float dire[3], struct qelem_basic_aabb *caabb){

	//AABBと直線の交差の判定
	bool intersect_flg = false;
	float k1, k2;
	float pip1[2], pip2[2];
	if (fabs(dire[0]) >= FLT_MIN){
		k1 = (cnode->aabb[0][0] - orig[0]) / dire[0];
		k2 = (cnode->aabb[0][1] - orig[0]) / dire[0];
		pip1[0] = orig[1] + k1*dire[1];
		pip1[1] = orig[2] + k1*dire[2];
		pip2[0] = orig[1] + k2*dire[1];
		pip2[1] = orig[2] + k2*dire[2];
		if (pip1[0] >= cnode->aabb[1][0] && pip1[0] <= cnode->aabb[1][1] && pip1[1] >= cnode->aabb[2][0] && pip1[1] <= cnode->aabb[2][1]) intersect_flg = true;
		if (pip2[0] >= cnode->aabb[1][0] && pip2[0] <= cnode->aabb[1][1] && pip2[1] >= cnode->aabb[2][0] && pip2[1] <= cnode->aabb[2][1]) intersect_flg = true;
	}
	if (!intersect_flg && fabs(dire[1]) >= FLT_MIN){
		k1 = (cnode->aabb[1][0] - orig[1]) / dire[1];
		k2 = (cnode->aabb[1][1] - orig[1]) / dire[1];
		pip1[0] = orig[0] + k1*dire[0];
		pip1[1] = orig[2] + k1*dire[2];
		pip2[0] = orig[0] + k2*dire[0];
		pip2[1] = orig[2] + k2*dire[2];
		if (pip1[0] >= cnode->aabb[0][0] && pip1[0] <= cnode->aabb[0][1] && pip1[1] >= cnode->aabb[2][0] && pip1[1] <= cnode->aabb[2][1]) intersect_flg = true;
		if (pip2[0] >= cnode->aabb[0][0] && pip2[0] <= cnode->aabb[0][1] && pip2[1] >= cnode->aabb[2][0] && pip2[1] <= cnode->aabb[2][1]) intersect_flg = true;
	}
	if (!intersect_flg && fabs(dire[2]) >= FLT_MIN){
		k1 = (cnode->aabb[2][0] - orig[2]) / dire[2];
		k2 = (cnode->aabb[2][1] - orig[2]) / dire[2];
		pip1[0] = orig[0] + k1 * dire[0];
		pip1[1] = orig[1] + k1*dire[1];
		pip2[0] = orig[0] + k2*dire[0];
		pip2[1] = orig[1] + k2*dire[1];
		if (pip1[0] >= cnode->aabb[0][0] && pip1[0] <= cnode->aabb[0][1] && pip1[1] >= cnode->aabb[1][0] && pip1[1] <= cnode->aabb[1][1]) intersect_flg = true;
		if (pip2[0] >= cnode->aabb[0][0] && pip2[0] <= cnode->aabb[0][1] && pip2[1] >= cnode->aabb[1][0] && pip2[1] <= cnode->aabb[1][1]) intersect_flg = true;
	}

	if (!intersect_flg) return caabb;

	if (cnode->LEAF_FLG){
		struct qelem_basic_aabb *naabb;
		naabb = (struct qelem_basic_aabb*)malloc(sizeof(struct qelem_basic_aabb));
		if (naabb == NULL){
			printf("Memory allocation error\n");
		}
		naabb->mynode = cnode;
		naabb->next = caabb;
		return naabb;
	}

	struct qelem_basic_aabb *naabb, *nhozon;
	naabb = caabb;
	for (int i = 0; i < 8; i++){
		if (cnode->dnode[i] != znode){
			nhozon = recursive_search_intersect_node_with_line(cnode->dnode[i], orig, dire, naabb);
			naabb = nhozon;
		}
	}
	return naabb;

}



////
//共通に使う関数
////
void free_qelem(struct qelem *bottom, struct qelem *ctop){

	struct qelem *celem;
	celem = bottom;
	while (1){
		struct qelem *next;
		next = celem->next;
		free(celem);
		if (next == ctop) break;
		celem = next;
	}
	free(ctop);

}


bool intersection_check_of_aabb_and_triagle(float aabb[3][2], float fps[3][3]){

	float limit[3][2];
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			if (i == 0){
				limit[j][0] = fps[i][j];
				limit[j][1] = fps[i][j];
			}
			else {
				if (limit[j][0] > fps[i][j]) limit[j][0] = fps[i][j];
				if (limit[j][1] < fps[i][j]) limit[j][1] = fps[i][j];
			}
		}
	}

	bool is_outside = false;
	bool is_inside = true;
	for (int j = 0; j < 3; j++){
		if (aabb[j][0] > limit[j][1] || aabb[j][1] < limit[j][0]) is_outside = true;
		if (aabb[j][0] > limit[j][0] || aabb[j][1] < limit[j][1]) is_inside = false;
	}
	if (is_outside) return false;
	if (is_inside) return true;

	is_inside = false;
	for (int i = 0; i < 3; i++){
		bool fp_is_outside = false;
		for (int j = 0; j < 3; j++){
			if (aabb[j][0] > fps[i][j] || aabb[j][1] < fps[i][j]) fp_is_outside = true;
		}
		if (!fp_is_outside) is_inside = true;
	}
	if (is_inside) return true;

	
	//エッジと三角形の交差判定
	float orig[3], dire[3], cc[3];
	int cn[3];
	cn[0] = 0;
	cn[1] = 0;
	cn[2] = 0;
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			orig[j] = aabb[j][cn[j]];
			if (j == i) dire[j] = aabb[j][(cn[j] + 1) % 2] - orig[j];
			else dire[j] = 0.0;
		}
		line_and_triangle_intersection(orig, dire, fps, cc);

		if (cc[0] >= 0.0 && cc[1] >= 0.0 && cc[2] >= 0.0){
			float tp[3];
			for (int j = 0; j < 3; j++){
				tp[j] = cc[0] * fps[0][j] + cc[1] * fps[1][j] + cc[2] * fps[2][j] - orig[j];
			}
			float dirh = dire[0] * dire[0] + dire[1] * dire[1] + dire[2] * dire[2];
			float seki = tp[0] * dire[0] + tp[1] * dire[1] + tp[2] * dire[2];
			if (seki <= dirh && seki >= 0.0) return true;
		}
	}

	cn[0] = 1;
	cn[1] = 1;
	cn[2] = 0;
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			orig[j] = aabb[j][cn[j]];
			if (j == i) dire[j] = aabb[j][(cn[j] + 1) % 2] - orig[j];
			else dire[j] = 0.0;
		}
		line_and_triangle_intersection(orig, dire, fps, cc);

		if (cc[0] >= 0.0 && cc[1] >= 0.0 && cc[2] >= 0.0){
			float tp[3];
			for (int j = 0; j < 3; j++){
				tp[j] = cc[0] * fps[0][j] + cc[1] * fps[1][j] + cc[2] * fps[2][j] - orig[j];
			}
			float dirh = dire[0] * dire[0] + dire[1] * dire[1] + dire[2] * dire[2];
			float seki = tp[0] * dire[0] + tp[1] * dire[1] + tp[2] * dire[2];
			if (seki <= dirh && seki >= 0.0) return true;
		}
	}
	cn[0] = 0;
	cn[1] = 1;
	cn[2] = 1;
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			orig[j] = aabb[j][cn[j]];
			if (j == i) dire[j] = aabb[j][(cn[j] + 1) % 2] - orig[j];
			else dire[j] = 0.0;
		}
		line_and_triangle_intersection(orig, dire, fps, cc);

		if (cc[0] >= 0.0 && cc[1] >= 0.0 && cc[2] >= 0.0){
			float tp[3];
			for (int j = 0; j < 3; j++){
				tp[j] = cc[0] * fps[0][j] + cc[1] * fps[1][j] + cc[2] * fps[2][j] - orig[j];
			}
			float dirh = dire[0] * dire[0] + dire[1] * dire[1] + dire[2] * dire[2];
			float seki = tp[0] * dire[0] + tp[1] * dire[1] + tp[2] * dire[2];
			if (seki <= dirh && seki >= 0.0) return true;
		}
	}

	cn[0] = 1;
	cn[1] = 0;
	cn[2] = 1;
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			orig[j] = aabb[j][cn[j]];
			if (j == i) dire[j] = aabb[j][(cn[j] + 1) % 2] - orig[j];
			else dire[j] = 0.0;
		}
		line_and_triangle_intersection(orig, dire, fps, cc);

		if (cc[0] >= 0.0 && cc[1] >= 0.0 && cc[2] >= 0.0){
			float tp[3];
			for (int j = 0; j < 3; j++){
				tp[j] = cc[0] * fps[0][j] + cc[1] * fps[1][j] + cc[2] * fps[2][j] - orig[j];
			}
			float dirh = dire[0] * dire[0] + dire[1] * dire[1] + dire[2] * dire[2];
			float seki = tp[0] * dire[0] + tp[1] * dire[1] + tp[2] * dire[2];
			if (seki <= dirh && seki >= 0.0) return true;
		}
	}

	return false;
}