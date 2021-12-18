#define _USE_MATH_DEFINES
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <cfloat>
#include <Eigen/Core>
#include<Eigen/IterativeLinearSolvers>
#include<Eigen/Dense>


using namespace std;


//3*2行列の求解
bool solve_three_by_two_system(float A[3][2], float b[3], float sol[2]){

	//行列操作
	//１行目の操作
	int max_val = A[0][0];
	int maxN = 0;
	for (int i = 0; i < 2; i++){
		if (max_val < A[i + 1][0]){
			max_val = A[i + 1][0];
			maxN = i + 1;
		}
	}
	if (maxN != 0){
		float ahozon = A[0][0];
		A[0][0] = A[maxN][0];
		A[maxN][0] = ahozon;
		ahozon = b[0];
		b[0] = b[maxN];
		b[maxN] = ahozon;
	}
	if (fabs(A[0][0]) < FLT_MIN) return false;

	A[0][1] = A[0][1] / A[0][0];
	b[0] = b[0] / A[0][0];
	A[0][0] = 1.0;
	for (int i = 0; i < 2; i++){
		A[i + 1][1] -= A[i + 1][0] * A[0][1];
		b[i + 1] -= A[i + 1][0] * b[0];
		A[i + 1][0] = 0.0;
	}

	//２行目の操作
	if (A[1][1] < A[2][1]){
		float ahozon = A[1][1];
		A[1][1] = A[2][1];
		A[2][1] = ahozon;
		ahozon = b[1];
		b[1] = b[2];
		b[2] = ahozon;
	}
	if (fabs(A[1][1]) < FLT_MIN) return false;

	b[1] /= A[1][1];
	A[1][1] = 1.0;
	b[2] -= A[2][1] * b[1];
	A[2][1] = 0.0;


	//逆代入
	if (fabs(b[2]) < FLT_MIN) return false;
	sol[1] = b[1];
	sol[0] = b[0] - A[0][1] * sol[1];
	return true;
}


//3*3行列の求解
bool solve_three_by_three_system(float A[3][3], float b[3], float sol[3]){

	float invA[3][3], detA;
	detA = 0.0;
	for (int i = 0; i < 3; i++){
		detA += A[i][0] * (A[(i + 1) % 3][1] * A[(i + 2) % 3][2] - A[(i + 1) % 3][2] * A[(i + 2) % 3][1]);
		for (int j = 0; j < 3; j++){
			invA[i][j] = A[(j + 1) % 3][(i + 1) % 3] * A[(j + 2) % 3][(i + 2) % 3] - A[(j + 2) % 3][(i + 1) % 3] * A[(j + 1) % 3][(i + 2) % 3];
		}
	}


	Eigen::MatrixXf Amat(3,3);
	Eigen::VectorXf svec(3),bvec(3);
	for (int i = 0; i < 3; i++){
		bvec(i) = b[i];
		for (int j = 0; j < 3; j++){
			Amat(i, j) = A[i][j];
		}
	}
	//sol = HesM.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(dLvec);
	Eigen::LDLT<Eigen::MatrixXf> decomp(Amat);
	svec = decomp.solve(bvec);


	if (fabs(detA) < FLT_MIN*100) return false;

	for (int i = 0; i < 3; i++){
		sol[i] = svec(i);
	}

	return true;

}

bool solve_three_by_three_system_use_psude_inverse(float A[3][3], float b[3], float sol[3]){

	float AtA[3][3],Atb[3];
	for (int i = 0; i < 3; i++){
		Atb[i] = 0.0;
		for (int j = 0; j < 3; j++){
			Atb[i] += A[j][i] * b[j];
			AtA[i][j] = 0.0;
			for (int k = 0; k < 3; k++){
				AtA[i][j] += A[k][i] * A[k][j];
			}
		}
	}

	Eigen::MatrixXf Amat(3, 3);
	Eigen::VectorXf svec(3), bvec(3);
	for (int i = 0; i < 3; i++){
		bvec(i) = Atb[i];
		for (int j = 0; j < 3; j++){
			Amat(i, j) = AtA[i][j];
		}
	}
	//sol = HesM.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(dLvec);
	Eigen::LDLT<Eigen::MatrixXf> decomp(Amat);
	svec = decomp.solve(bvec);

	for (int i = 0; i < 3; i++){
		sol[i] = svec(i);
	}

	return true;
}


//直線と面の交点計算（重心座標を求める）
void line_and_triangle_intersection(float orig[3],float dire[3],float fp[3][3],float cc[3]){

	float ab[3], ac[3],nd[3];
	float x, y, ck;
	for (int i = 0; i < 3; i++){
		ab[i] = fp[1][i] - fp[0][i];
		ac[i] = fp[2][i] - fp[0][i];
	}
	for (int i = 0; i < 3; i++){
		nd[i] = ab[(i + 1) % 3] * ac[(i + 2) % 3] - ab[(i + 2) % 3] * ac[(i + 1) % 3];
	}
	float ndh = sqrt(nd[0] * nd[0] + nd[1] * nd[1] + nd[2] * nd[2]);
	for (int i = 0; i < 3; i++){
		nd[i] /= ndh;
	}
	float seki1 = 0.0;
	float seki2 = 0.0;
	for (int i = 0; i < 3; i++){
		seki1 += nd[i] * (fp[0][i]-orig[i]);
		seki2 += nd[i] * dire[i];
	}
	if (fabs(seki2) < 1.0e-6){
		if (fabs(seki1) < 1.0e-6){
			//面と直線は同一平面上
			float pnd[3];
			for (int i = 0; i < 3; i++){
				pnd[i] = fp[0][i] - orig[i];
			}
			float seki = pnd[0] * dire[0] + pnd[1] * dire[1] + pnd[2] * dire[2];
			for (int i = 0; i < 3; i++){
				pnd[i] -= seki*dire[i]/(dire[0]*dire[0]+dire[1]*dire[1]+dire[2]*dire[2]);
			}
			float pndh = sqrt(pnd[0] * pnd[0] + pnd[1] * pnd[1] + pnd[2] * pnd[2]);
			if (pndh > 1.0e-5){
				for (int i = 0; i < 3; i++){
					pnd[i] /= pndh;
				}
			}
			else {
				for (int i = 0; i < 3; i++){
					pnd[i] = fp[1][i] - orig[i];
				}
				seki = pnd[0] * dire[0] + pnd[1] * dire[1] + pnd[2] * dire[2];
				for (int i = 0; i < 3; i++){
					pnd[i] -= seki*dire[i] / (dire[0] * dire[0] + dire[1] * dire[1] + dire[2] * dire[2]);
				}
				pndh = sqrt(pnd[0] * pnd[0] + pnd[1] * pnd[1] + pnd[2] * pnd[2]);
				for (int i = 0; i < 3; i++){
					pnd[i] /= pndh;
				}
			}
			seki1 = pnd[0] * (fp[0][0] - orig[0]) + pnd[1] * (fp[0][1] - orig[1]) + pnd[2] * (fp[0][2] - orig[2]);
			seki2 = pnd[0] * (fp[1][0] - orig[0]) + pnd[1] * (fp[1][1] - orig[1]) + pnd[2] * (fp[1][2] - orig[2]);
			float seki3 = pnd[0] * (fp[2][0] - orig[0]) + pnd[1] * (fp[2][1] - orig[1]) + pnd[2] * (fp[2][2] - orig[2]);
			
			if (seki1 > 0.0 && seki2 > 0.0 && seki3 > 0.0){
				//面と直線は交わらない
				cc[0] = -1.0;
				cc[1] = -1.0;
				cc[2] = -1.0;
			}
			else if(seki1 < 0.0 && seki2 < 0.0 && seki3 < 0.0){
				//面と直線は交わらない
				cc[0] = -1.0;
				cc[1] = -1.0;
				cc[2] = -1.0;
			}
			else {
				cc[0] = 1.0 / 3.0;
				cc[1] = 1.0 / 3.0;
				cc[2] = 1.0 / 3.0;
			}
		}
		else {
			//面と直線は交わらない
			cc[0] = -1.0;
			cc[1] = -1.0;
			cc[2] = -1.0;
		}
		return;
	}
	ck = seki1 / seki2;

	float prop[3];
	for (int i = 0; i < 3; i++){
		prop[i] = orig[i] + ck*dire[i];
	}
	float abap = ab[0] * (prop[0]-fp[0][0]) + ab[1] * (prop[1]-fp[0][1]) + ab[2] * (prop[2]-fp[0][2]);
	float acap = ac[0] * (prop[0]-fp[0][0]) + ac[1] * (prop[1]-fp[0][1]) + ac[2] * (prop[2]-fp[0][2]);
	float abab = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
	float acac = ac[0] * ac[0] + ac[1] * ac[1] + ac[2] * ac[2];
	float abac = ab[0] * ac[0] + ab[1] * ac[1] + ab[2] * ac[2];
	x = (abap*acac - acap*abac) / (abab*acac-abac*abac);
	y = (abap*abac - acap*abab) / (abac*abac-abab*acac);

	cc[0] = 1.0 - x - y;
	cc[1] = x;
	cc[2] = y;
}


//与えられた点に三角形上で最も近い点を計算（重心座標系を求める）
float closest_point_of_triangle(float pc[3], float fps[3][3], float cord[3]){

	
	float ab[3], ac[3],nd[3],ah[3];
	for (int i = 0; i < 3; i++){
		ab[i] = fps[1][i] - fps[0][i];
		ac[i] = fps[2][i] - fps[0][i];
	}
	for (int i = 0; i < 3; i++){
		nd[i] = ab[(i + 1) % 3] * ac[(i + 2) % 3] - ab[(i + 2) % 3] * ac[(i + 1) % 3];
	}
	float ndh = sqrt(nd[0] * nd[0] + nd[1] * nd[1] + nd[2] * nd[2]);
	if (ndh > 1.0e-6){
		for (int i = 0; i < 3; i++){
			nd[i] /= ndh;
		}
		float seki = (pc[0] - fps[0][0]) * nd[0] + (pc[1] - fps[0][1]) * nd[1] + (pc[2] - fps[0][2]) * nd[2];
		for (int i = 0; i < 3; i++){
			ah[i] = pc[i] - fps[0][i] - seki*nd[i];
		}

		float abab = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
		float acac = ac[0] * ac[0] + ac[1] * ac[1] + ac[2] * ac[2];
		float abac = ab[0] * ac[0] + ab[1] * ac[1] + ab[2] * ac[2];
		float abah = ab[0] * ah[0] + ab[1] * ah[1] + ab[2] * ah[2];
		float acah = ac[0] * ah[0] + ac[1] * ah[1] + ac[2] * ah[2];
		float k1 = (abah*acac - acah*abac) / (abab*acac - abac*abac);
		float k2 = (abah*abac - acah*abab) / (abac*abac - acac*abab);
		if (k1 > 0.0 && k2 > 0.0 && (1.0 - k1 - k2) > 0.0){
			cord[0] = 1.0 - k1 - k2;
			cord[1] = k1;
			cord[2] = k2;
			float fhp[3];
			for (int i = 0; i < 3; i++){
				fhp[i] = cord[0] * fps[0][i] + cord[1] * fps[1][i] + cord[2] * fps[2][i];
			}
			float kyori = sqrt(pow(fhp[0] - pc[0], 2) + pow(fhp[1] - pc[1], 2) + pow(fhp[2] - pc[2], 2));
			return kyori;
		}
		else {
			int mini = -1;
			float minh;
			for (int i = 0; i < 3; i++){
				float lm[3],pm[3];
				for (int j = 0; j < 3; j++){
					lm[j] = fps[(i + 1) % 3][j] - fps[i][j];
					pm[j] = pc[j] - fps[i][j];
				}
				float seki = (lm[0] * pm[0] + lm[1] * pm[1] + lm[2] * pm[2]);
				float lmh = lm[0] * lm[0] + lm[1] * lm[1] + lm[2] * lm[2];
				if (lmh > 1.0e-6){
					seki /= lmh;
					float kyori;
					float mycord[3];
					if (seki < 0.0){
						mycord[i] = 1.0;
						mycord[(i + 1) % 3] = 0.0;
						mycord[(i + 2) % 3] = 0.0;
						kyori = sqrt(pow(pc[0] - fps[i][0], 2) + pow(pc[1] - fps[i][1], 2) + pow(pc[2] - fps[i][2], 2));
					}
					else if (seki > 1.0){
						mycord[i] = 0.0;
						mycord[(i + 1) % 3] = 1.0;
						mycord[(i + 2) % 3] = 0.0;
						kyori = sqrt(pow(pc[0] - fps[(i + 1) % 3][0], 2) + pow(pc[1] - fps[(i + 1) % 3][1], 2) + pow(pc[2] - fps[(i + 1) % 3][2], 2));
					}
					else {
						mycord[i] = 1.0-seki;
						mycord[(i + 1) % 3] = seki;
						mycord[(i + 2) % 3] = 0.0;
						kyori = pow(pc[0] - (1.0 - seki)*fps[i][0] - seki*fps[(i + 1) % 3][0], 2) + pow(pc[1] - (1.0 - seki)*fps[i][1] - seki*fps[(i + 1) % 3][1], 2);
						kyori += pow(pc[2] - (1.0 - seki)*fps[i][2] - seki*fps[(i + 1) % 3][2], 2);
						kyori = sqrt(kyori);
					}
					if (mini == -1){
						mini = i;
						minh = kyori;
						for (int j = 0; j<3; j++){
							cord[j] = mycord[j];
						}
					}
					else if (minh > kyori){
						mini = i;
						minh = kyori;
						for (int j = 0; j<3; j++){
							cord[j] = mycord[j];
						}
					}
				}
			}
			if (mini != -1){
				return minh;
			}
			else{
				printf("all edge degenerate\n");
				float  gp[3];
				for (int i = 0; i < 3; i++){
					cord[i] = 1.0 / 3.0;
					gp[i] = (fps[0][i] + fps[1][i] + fps[2][i]) / 3.0;
				}
				float kyori = sqrt(pow(gp[0] - pc[0], 2) + pow(gp[1] - pc[1], 2) + pow(gp[2] - pc[2], 2));
				return kyori;
			}
		}

	}
	else {
		printf("triangle degenerate\n");
		float abab = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
		float acac = ac[0] * ac[0] + ac[1] * ac[1] + ac[2] * ac[2];
		if (abab > acac){
			float seki = (ab[0] * (pc[0] - fps[0][0]) + ab[1] * (pc[1] - fps[0][1]) + ab[2] * (pc[2] - fps[0][2]))/abab;
			if (seki < 0){
				cord[0] = 1.0;
				cord[1] = 0.0;
				cord[2] = 0.0;
				float kyori = sqrt(pow(pc[0] - fps[0][0], 2) + pow(pc[1] - fps[0][1], 2) + pow(pc[2] - fps[0][2], 2));
				return kyori;
			}
			else if (seki > 1.0){
				cord[0] = 0.0;
				cord[1] = 1.0;
				cord[2] = 0.0;
				float kyori = sqrt(pow(pc[0] - fps[1][0], 2) + pow(pc[1] - fps[1][1], 2) + pow(pc[2] - fps[1][2], 2));
				return kyori;
			}
			else {
				cord[0] = 1.0 - seki;
				cord[1] = seki;
				cord[2] = 0.0;
				float kyori = pow(pc[0] - cord[0] * fps[0][0] - cord[1] * fps[1][0], 2) + pow(pc[1] - cord[0] * fps[0][1] - cord[1] * fps[1][1], 2);
				kyori += pow(pc[2] - cord[0] * fps[0][2] - cord[1] * fps[1][2], 2);
				return sqrt(kyori);
			}
		}
		else {
			float seki = (ac[0] * (pc[0] - fps[0][0]) + ac[1] * (pc[1] - fps[0][1]) + ac[2] * (pc[2] - fps[0][2]))/acac;
			if (seki < 0){
				cord[0] = 1.0;
				cord[1] = 0.0;
				cord[2] = 0.0;
				float kyori = sqrt(pow(pc[0] - fps[0][0], 2) + pow(pc[1] - fps[0][1], 2) + pow(pc[2] - fps[0][2], 2));
				return kyori;
			}
			else if (seki > 1.0){
				cord[0] = 0.0;
				cord[1] = 0.0;
				cord[2] = 1.0;
				float kyori = sqrt(pow(pc[0] - fps[2][0], 2) + pow(pc[1] - fps[2][1], 2) + pow(pc[2] - fps[2][2], 2));
				return kyori;
			}
			else {
				cord[0] = 1.0 - seki;
				cord[1] = 0.0;
				cord[2] = seki;
				float kyori = pow(pc[0] - cord[0] * fps[0][0] - cord[2] * fps[2][0], 2) + pow(pc[1] - cord[0] * fps[0][1] - cord[2] * fps[2][1], 2);
				kyori += pow(pc[2] - cord[0] * fps[0][2] - cord[2] * fps[2][2], 2);
				return sqrt(kyori);
			}
		}
	}

}


//三角形と直線の交点の計算（重心座標系を求める）
bool intersection_of_triangle_with_line(float orig[3], float dire[3], float fps[3][3], float cord[3]){

	float ab[3], ac[3], nd[3], ap[3];
	for (int i = 0; i < 3; i++){
		ab[i] = fps[1][i] - fps[0][i];
		ac[i] = fps[2][i] - fps[0][i];
		ap[i] = orig[i] - fps[0][i];
	}
	for (int i = 0; i < 3; i++){
		nd[i] = ab[(i + 1) % 3] * ac[(i + 2) % 3] - ab[(i + 2) % 3] * ac[(i + 1) % 3];
	}
	float ndh = sqrt(nd[0] * nd[0] + nd[1] * nd[1] + nd[2] * nd[2]);
	if (ndh >= FLT_MIN){
		for (int i = 0; i < 3; i++){
			nd[i] /= ndh;
		}
		float seki1 = nd[0] * ap[0] + nd[1] * ap[1] + nd[2] * ap[2];
		float seki2 = nd[0] * dire[0] + nd[1] * dire[1] + nd[2] * dire[2];

		if (fabs(seki2) > FLT_MIN){
			float k = -seki1 / seki2;
			float ah[3];
			for (int i = 0; i < 3; i++){
				ah[i] = orig[i] + k*dire[i] - fps[0][i];
			}
			float abab = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
			float acac = ac[0] * ac[0] + ac[1] * ac[1] + ac[2] * ac[2];
			float abac = ab[0] * ac[0] + ab[1] * ac[1] + ab[2] * ac[2];
			float abah = ab[0] * ah[0] + ab[1] * ah[1] + ab[2] * ah[2];
			float acah = ac[0] * ah[0] + ac[1] * ah[1] + ac[2] * ah[2];
			cord[1] = (abah*acac - acah*abac) / (abab*acac - abac*abac);
			cord[2] = (abah*abac - acah*abab) / (abac*abac - abab*acac);
			cord[0] = 1.0 - cord[1] - cord[2];
			if (cord[0] >= 0.0 && cord[1] >= 0.0 && cord[2] >= 0.0) return true;
			else return false;
		}
		else {
			if (fabs(seki1) <= FLT_MIN){
				//三角形の平面上に直線が乗っている
				bool intersect_flg = false;
				for (int i = 0; i < 3; i++){
					float A[3][2], b[3], sol[2];
					for (int j = 0; j < 3; j++){
						A[j][0] = fps[(i + 1) % 3][j] - fps[i][j];
						A[j][1] = -dire[j];
						b[j] = orig[j] - fps[i][j];
					}
					bool sol_flg = solve_three_by_two_system(A, b, sol);
					if (sol_flg && sol[0] >= 0.0 && sol[1] <= 1.0) intersect_flg = true;
				}
				cord[0] = -1.0;
				cord[1] = -1.0;
				cord[2] = -1.0;
				if (intersect_flg) return true;
				else return false;
			}
			else {
				//三角形の平面と直線が平行で交わらない
				cord[0] = -1.0;
				cord[1] = -1.0;
				cord[2] = -1.0;
				return false;
			}
		}

	}
	else {
		//三角形の形状がつぶれている
		float abab = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
		float acac = ac[0] * ac[0] + ac[1] * ac[1] + ac[2] * ac[2];
		if (abab > acac){
			float A[3][2], b[3], sol[2];
			for (int i = 0; i < 3; i++){
				A[i][0] = fps[1][i] - fps[0][i];
				A[i][1] = -dire[i];
				b[i] = orig[i] - fps[0][i];
			}
			bool sol_flg = solve_three_by_two_system(A, b, sol);
			if (sol_flg && sol[0] >= 0.0 && sol[0] <= 1.0){
				cord[0] = 1.0 - sol[0];
				cord[1] = sol[0];
				cord[2] = 0.0;
				return true;
			}
			else {
				cord[0] = -1;
				cord[1] = -1;
				cord[2] = -1;
				return false;
			}
		}
		else{
			float A[3][2], b[3], sol[2];
			for (int i = 0; i < 3; i++){
				A[i][0] = fps[2][i] - fps[0][i];
				A[i][1] = -dire[i];
				b[i] = orig[i] - fps[0][i];
			}
			bool sol_flg = solve_three_by_two_system(A, b, sol);
			if (sol_flg && sol[0] >= 0.0 && sol[0] <= 1.0){
				cord[0] = 1.0 - sol[0];
				cord[1] = 0.0;
				cord[2] = sol[0];
				return true;
			}
			else {
				cord[0] = -1;
				cord[1] = -1;
				cord[2] = -1;
				return false;
			}
		}
	}

}


//回転行列から四元数を求める
void derive_quaternion_from_rot_mat(float Rot[3][3], float quat[4]){

	
	float maxv;
	int maxn;
	for (int j = 0; j < 4; j++){
		if (j == 0){
			quat[j] = 1 + Rot[0][0] + Rot[1][1] + Rot[2][2];
			maxn = 0;
			maxv = quat[j];
		}
		else quat[j] = 1 + Rot[j - 1][j - 1] - Rot[j % 3][j % 3] - Rot[(j + 1) % 3][(j + 1) % 3];
		if (maxv < quat[j]){
			maxn = j;
			maxv = quat[j];
		}
	}
	if (maxn == 0){
		quat[0] = sqrt(quat[0] / 4.0);
		quat[1] = (Rot[2][1] - Rot[1][2]) / quat[0] / 4.0;
		quat[2] = (Rot[0][2] - Rot[2][0]) / quat[0] / 4.0;
		quat[3] = (Rot[1][0] - Rot[0][1]) / quat[0] / 4.0;
	}
	else if (maxn == 1){
		quat[1] = sqrt(quat[1] / 4.0);
		quat[0] = (Rot[2][1] - Rot[1][2]) / quat[1] / 4.0;
		quat[2] = (Rot[0][1] + Rot[1][0]) / quat[1] / 4.0;
		quat[3] = (Rot[0][2] + Rot[2][0]) / quat[1] / 4.0;
	}
	else if (maxn == 2){
		quat[2] = sqrt(quat[2] / 4.0);
		quat[0] = (Rot[0][2] - Rot[2][0]) / quat[2] / 4.0;
		quat[1] = (Rot[0][1] + Rot[1][0]) / quat[2] / 4.0;
		quat[3] = (Rot[1][2] + Rot[2][1]) / quat[2] / 4.0;
	}
	else {
		quat[3] = sqrt(quat[3] / 4.0);
		quat[0] = (Rot[1][0] - Rot[0][1]) / quat[3] / 4.0;
		quat[1] = (Rot[0][2] + Rot[2][0]) / quat[3] / 4.0;
		quat[2] = (Rot[1][2] + Rot[2][1]) / quat[3] / 4.0;
	}

}



//回転行列から四元数を求める
void derive_quaternion_from_rot_mat(double Rot[3][3], double quat[4]){


	double maxv;
	int maxn;
	for (int j = 0; j < 4; j++){
		if (j == 0){
			quat[j] = 1 + Rot[0][0] + Rot[1][1] + Rot[2][2];
			maxn = 0;
			maxv = quat[j];
		}
		else quat[j] = 1 + Rot[j - 1][j - 1] - Rot[j % 3][j % 3] - Rot[(j + 1) % 3][(j + 1) % 3];
		if (maxv < quat[j]){
			maxn = j;
			maxv = quat[j];
		}
	}
	if (maxn == 0){
		quat[0] = sqrt(quat[0] / 4.0);
		quat[1] = (Rot[2][1] - Rot[1][2]) / quat[0] / 4.0;
		quat[2] = (Rot[0][2] - Rot[2][0]) / quat[0] / 4.0;
		quat[3] = (Rot[1][0] - Rot[0][1]) / quat[0] / 4.0;
	}
	else if (maxn == 1){
		quat[1] = sqrt(quat[1] / 4.0);
		quat[0] = (Rot[2][1] - Rot[1][2]) / quat[1] / 4.0;
		quat[2] = (Rot[0][1] + Rot[1][0]) / quat[1] / 4.0;
		quat[3] = (Rot[0][2] + Rot[2][0]) / quat[1] / 4.0;
	}
	else if (maxn == 2){
		quat[2] = sqrt(quat[2] / 4.0);
		quat[0] = (Rot[0][2] - Rot[2][0]) / quat[2] / 4.0;
		quat[1] = (Rot[0][1] + Rot[1][0]) / quat[2] / 4.0;
		quat[3] = (Rot[1][2] + Rot[2][1]) / quat[2] / 4.0;
	}
	else {
		quat[3] = sqrt(quat[3] / 4.0);
		quat[0] = (Rot[1][0] - Rot[0][1]) / quat[3] / 4.0;
		quat[1] = (Rot[0][2] + Rot[2][0]) / quat[3] / 4.0;
		quat[2] = (Rot[1][2] + Rot[2][1]) / quat[3] / 4.0;
	}

}


//四元数から回転行列を求める
void derive_rot_mat_from_quaternion(double quat[4], double Rot[3][3]){

	double quah = sqrt(quat[0] * quat[0] + quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);
	for (int i = 0; i < 4; i++){
		quat[i] /= quah;
	}

	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			if (i == j) Rot[i][j] = quat[0] * quat[0] + quat[i + 1] * quat[j + 1]
				- quat[(i + 1) % 3 + 1] * quat[(j + 1) % 3 + 1] - quat[(i + 2) % 3 + 1] * quat[(j + 2) % 3 + 1];
			else if ((i + 1) % 3 == j) Rot[i][j] = 2.0*(quat[i + 1] * quat[j + 1] - quat[0] * quat[4 - i - j]);
			else Rot[i][j] = 2.0*(quat[i + 1] * quat[j + 1] + quat[0] * quat[4 - i - j]);
		}

	}

}




//四元数から回転行列を求める
void derive_rot_mat_from_quaternion(float quat[4], float Rot[3][3]){

	float quah = sqrt(quat[0] * quat[0] + quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);
	for (int i = 0; i < 4; i++){
		quat[i] /= quah;
	}

	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			if (i == j) Rot[i][j] = quat[0] * quat[0] + quat[i + 1] * quat[j + 1]
				- quat[(i + 1) % 3 + 1] * quat[(j + 1) % 3 + 1] - quat[(i + 2) % 3 + 1] * quat[(j + 2) % 3 + 1];
			else if ((i + 1) % 3 == j) Rot[i][j] = 2.0*(quat[i + 1] * quat[j + 1] - quat[0] * quat[4 - i - j]);
			else Rot[i][j] = 2.0*(quat[i + 1] * quat[j + 1] + quat[0] * quat[4 - i - j]);
		}

	}

}


//直方体と直線の交わりの判定
bool intersection_of_cuboid_with_line(float orig[3], float dire[3], float rect[3][2], float& depth){

	bool intersect = false;
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 2; j++){
			float coeff = rect[i][j] - orig[i];
			if (fabs(dire[i]) > 1.0e-6){
				coeff /= dire[i];
				float ipn[3];
				for (int k = 0; k < 3; k++){
					ipn[k] = orig[k] + coeff*dire[k];
				}
				if (ipn[(i + 1) % 3] >= rect[(i + 1) % 3][0] && ipn[(i + 1) % 3] <= rect[(i + 1) % 3][1]
					&& ipn[(i + 2) % 3] >= rect[(i + 2) % 3][0] && ipn[(i + 2) % 3] <= rect[(i + 2) % 3][1]){
					if (!intersect) {
						intersect = true;
						depth = coeff;
					}
					else if (depth > coeff) depth = coeff;
				}
			}
		}
	}

	return intersect;

}


//平面と三角形の交わりの判定(重心座標系で返す),state( -1:交差なし  0:通常の交差 1:頂点で交差　2:エッジが含まれる　3:全体が含まれる)
void intersection_of_triangle_with_plane(float orig[3], float norm[3], float fps[3][3], float icc[3], int &state){

	state = -1;
	float nh = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);


	float ncord[3];
	for (int i = 0; i < 3; i++){
		float seki = 0.0;
		for (int j = 0; j < 3; j++){
			seki += norm[j] * (fps[i][j] - orig[j]) / nh;
		}
		ncord[i] = seki;
	}
	float ecord[3];
	int intpN = 0;
	int inteN = 0;
	/*for (int i = 0; i < 3; i++){
		float dif = ncord[(i + 1) % 3] - ncord[i];
		if (dif > 0.0 && dif < 1.0e-7) dif = 1.0e-7;
		if (dif < 0.0 && dif > -1.0e-7) dif = -1.0e-7;
		ecord[i] = -ncord[i] / dif;
		icc[i] = ecord[i];
		if (fabs(ncord[i]) < 1.0e-7) intpN++;
		if (ecord[i] >= 0 && ecord[i] <= 1.0) inteN++;
	}*/
	for (int i = 0; i < 3; i++){
		float dif = ncord[(i + 1) % 3] - ncord[i];
		if (dif == 0) ecord[i] = ncord[i] * 1000;
		else ecord[i] = -ncord[i] / dif;
		icc[i] = ecord[i];
		if (ncord[i] == 0) intpN++;
		if (ecord[i] >= -0.0 && ecord[i] <= 1.0) inteN++;
	}
	if (inteN > 1) state = 0;
	if (intpN == 1) state = 1;
	if (intpN == 2) state = 2;
	if (intpN == 3) state = 3;


}


//二面角の計算（凸なら正、凹なら負でrad値を返す//1-2-3,2-1-4が向き付き三角形をなす）
float calc_dihedral_angle_of_triangles(float fps[4][3]){

	float nd1[3], nd2[3];
	float td[3];
	for (int i = 0; i < 3; i++){
		nd1[i] = (fps[1][(i + 1) % 3] - fps[0][(i + 1) % 3])*(fps[2][(i + 2) % 3] - fps[0][(i + 2) % 3])
			-(fps[1][(i + 2) % 3] - fps[0][(i + 2) % 3])*(fps[2][(i + 1) % 3] - fps[0][(i + 1) % 3]);
		nd2[i] = (fps[3][(i + 1) % 3] - fps[0][(i + 1) % 3])*(fps[1][(i + 2) % 3] - fps[0][(i + 2) % 3])
			- (fps[3][(i + 2) % 3] - fps[0][(i + 2) % 3])*(fps[1][(i + 1) % 3] - fps[0][(i + 1) % 3]);
		td[i] = fps[3][i] - fps[0][i];
	}
	float ndh1 = sqrt(nd1[0] * nd1[0] + nd1[1] * nd1[1] + nd1[2] * nd1[2]);
	float ndh2 = sqrt(nd2[0] * nd2[0] + nd2[1] * nd2[1] + nd2[2] * nd2[2]);
	if (ndh1 < 1.0e-6) ndh1 = 1.0e-6;
	if (ndh2 < 1.0e-6) ndh2 = 1.0e-6;
	float cosv = (nd1[0] * nd2[0] + nd1[1] * nd2[1] + nd1[2] * nd2[2]) / ndh1 / ndh2;
	if (cosv > 1.0) cosv = 1.0;
	if (cosv < -1.0) cosv = -1.0;
	float theta = acos(cosv);

	float seki = td[0] * nd1[0] + td[1] * nd1[1] + td[2] * nd1[2];
	if (seki > 0.0)theta = -theta;

	return theta;

}


//三角形の外接円の半径と最小エッジの比（形状評価）の計算
float calc_ratio_of_circumcircle_radius_to_min_edge(float fps[3][3]){

	float rad,minh;
	float fph[3];
	for (int i = 0; i < 3; i++){
		fph[i] = 0.0;
		for (int j = 0; j < 3; j++){
			fph[i] += pow(fps[(i + 1) % 3][j] - fps[i][j], 2);
		}
		fph[i] = sqrt(fph[i]);

		if (i == 0) minh = fph[i];
		else if (minh > fph[i]) minh = fph[i];
	}
	float denom;
	denom = (fph[0] + fph[1] + fph[2])*(-fph[0] + fph[1] + fph[2])*(fph[0] - fph[1] + fph[2])*(fph[0] + fph[1] - fph[2]);
	rad = (fph[0] * fph[1]*fph[2]) / sqrt(denom);

	return minh / rad;
}


//２直線間の距離の計算
float calc_distance_between_lines(float orig1[3], float dire1[3], float orig2[3], float dire2[3]){

	//一応正規化
	float dh1 = sqrt(dire1[0] * dire1[0] + dire1[1] * dire1[1] + dire1[2] * dire1[2]);
	float dh2 = sqrt(dire2[0] * dire2[0] + dire2[1] * dire2[1] + dire2[2] * dire2[2]);
	if (dh1 < 1.0e-10) dh1 = 1.0e-10;
	if (dh2 < 1.0e-10) dh2 = 1.0e-10;
	for (int i = 0; i < 3; i++){
		dire1[i] /= dh1;
		dire2[i] /= dh2;
	}


	//直交軸の計算
	float Xd[3], Yd[3];
	float seki = 0.0;
	for (int i = 0; i < 3; i++){
		Xd[i] = orig2[i] - orig1[i];
		seki += Xd[i] * dire1[i];
	}
	for (int i = 0; i < 3; i++){
		Xd[i] -= seki*dire1[i];
	}
	float Xdh = sqrt(Xd[0] * Xd[0] + Xd[1] * Xd[1] + Xd[2] * Xd[2]);
	if (Xdh < 1.0e-10) return Xdh;
	for (int i = 0; i < 3; i++){
		Xd[i] /= Xdh;
	}
	for (int i = 0; i < 3; i++){
		Yd[i] = dire1[(i + 1) % 3] * Xd[(i + 2) % 3] - dire1[(i + 2) % 3] * Xd[(i + 1) % 3];
	}
	float Ydh = sqrt(Yd[0] * Yd[0] + Yd[1] * Yd[1] + Yd[2] * Yd[2]);
	if (Ydh < 1.0e-10) return Ydh;
	for (int i = 0; i < 3; i++){
		Yd[i] /= Ydh;
	}


	//最短距離の計算
	float ldir[2];
	ldir[0] = Xd[0] * dire2[0] + Xd[1] * dire2[1] + Xd[2] * dire2[2];
	ldir[1] = Yd[0] * dire2[0] + Yd[1] * dire2[1] + Yd[2] * dire2[2];
	float ldh = sqrt(ldir[0] * ldir[0] + ldir[1] * ldir[1]);
	seki = -ldir[0] * Xdh;
	if (ldh < 1.0e-10) return Xdh;
	float ph[2];
	ph[0] = seki / ldh / ldh*ldir[0] + Xdh;
	ph[1] = seki / ldh / ldh*ldir[1];
	float kyori = sqrt(ph[0] * ph[0] + ph[1] * ph[1]);
	return kyori;

}


//直線と線分の距離の計算
float calc_distance_between_line_and_segment(float orig[3], float dire[3], float end1[3], float end2[3]){

	float dh = sqrt(dire[0] * dire[0] + dire[1] * dire[1] + dire[2] * dire[2]);
	if (dh < 1.0e-10) dh = 1.0e-10;
	for (int i = 0; i < 3; i++){
		dire[i] /= dh;
	}


	//直交軸の計算
	float Xd[3], Yd[3];
	float seki = 0.0;
	for (int i = 0; i < 3; i++){
		Xd[i] = end1[i] - orig[i];
		seki += Xd[i] * dire[i];
	}
	for (int i = 0; i < 3; i++){
		Xd[i] -= seki*dire[i];
	}
	float Xdh = sqrt(Xd[0] * Xd[0] + Xd[1] * Xd[1] + Xd[2] * Xd[2]);
	if (Xdh < 1.0e-10) return Xdh;
	for (int i = 0; i < 3; i++){
		Xd[i] /= Xdh;
	}
	for (int i = 0; i < 3; i++){
		Yd[i] = dire[(i + 1) % 3] * Xd[(i + 2) % 3] - dire[(i + 2) % 3] * Xd[(i + 1) % 3];
	}
	float Ydh = sqrt(Yd[0] * Yd[0] + Yd[1] * Yd[1] + Yd[2] * Yd[2]);
	if (Ydh < 1.0e-10) return Ydh;
	for (int i = 0; i < 3; i++){
		Yd[i] /= Ydh;
	}


	//最短距離の計算
	float ldir[2];
	ldir[0] = 0.0;
	ldir[1] = 0.0;
	for (int i = 0; i < 3; i++){
		ldir[0] += Xd[i] * (end2[i] - end1[i]);
		ldir[1] += Yd[i] * (end2[i] - end1[i]);
	}
	float ldh = sqrt(ldir[0] * ldir[0] + ldir[1] * ldir[1]);
	if (ldh < 1.0e-10) return Xdh;
	float k = -Xdh / ldh / ldh*ldir[0];
	if (k < 0.0 || k > 1.0) k = 1.0;
	float kyori;
	kyori = sqrt(pow(k*ldir[0] + Xdh, 2) + k*k*ldir[1] * ldir[1]);
	if (kyori > Xdh) kyori = Xdh;

	return kyori;

}


//点と直線の距離
float calc_distance_from_line(float orig[3], float dire[3], float poi[3]){

	//正規化
	float dh = sqrt(dire[0] * dire[0] + dire[1] * dire[1] + dire[2] * dire[2]);
	if (dh < 1.0e-10) dh = 1.0e-10;
	for (int i = 0; i < 3; i++){
		dire[i] /= dh;
	}


	float seki = 0.0;
	for (int i = 0; i < 3; i++){
		seki += dire[i] * (poi[i] - orig[i]);
	}
	float hp[3];
	for (int i = 0; i < 3; i++){
		hp[i] = poi[i] - orig[i] - seki*dire[i];
	}
	float kyori = sqrt(hp[0] * hp[0] + hp[1] * hp[1] + hp[2] * hp[2]);
	return kyori;

}