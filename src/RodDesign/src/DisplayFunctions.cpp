#include <QtOpenGL>

#include "DisplayFunctions.h"
#include "curve_network.h"
#include "RodNetwork.h"
#include "basic_geometric_calculation.h"

//色の選択
void choose_color_data(int cflg){

	if (cflg == 0){
		//金
		GLfloat mat_specular[] = { 0.0, 0.0, 0.0, 1.0 };
		GLfloat mat_diffuse[] = { 0.5, 0.5, 0.0, 1.0 };
		GLfloat mat_ambient[] = { 0.60, 0.60, 0.50, 1.0 };
		GLfloat mat_shininess = 32.0;
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);
	}
	else if (cflg == 1){
		//グレー
		GLfloat mat_specular[] = { 0.0, 0.0, 0.0, 1.0 };
		GLfloat mat_diffuse[] = { 0.55, 0.55, 0.55, 1.0 };
		GLfloat mat_ambient[] = { 0.70, 0.70, 0.70, 1.0 };
		GLfloat mat_shininess = 32.0;
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);
	}
	else if (cflg == 2){
		//青
		GLfloat mat_specular[] = { 0.0, 0.1, 0.06, 1.0 };
		GLfloat mat_diffuse[] = { 0.0, 0.50980392, 0.50980392, 1.0 };
		GLfloat mat_ambient[] = { 0.50196078, 0.50196078, 0.50196078, 1.0 };
		GLfloat mat_shininess = 32.0;
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);
	}
	else if (cflg == 3){
		//プラスチック 白
		GLfloat mat_specular[] = { 0.70, 0.70, 0.70, 1.0 };
		GLfloat mat_diffuse[] = { 0.55, 0.55, 0.55, 1.0 };
		GLfloat mat_ambient[] = { 0.0, 0.0, 0.0, 1.0 };
		GLfloat mat_shininess = 32.0;
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);
	}
	else if (cflg == 4){
		//クローム
		GLfloat mat_specular[] = { 0.774597, 0.774597, 0.774597, 1.0 };
		GLfloat mat_diffuse[] = { 0.4, 0.4, 0.4, 1.0 };
		GLfloat mat_ambient[] = { 0.25, 0.25, 0.25, 1.0 };
		GLfloat mat_shininess = 76.8;
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);
	}
	else if (cflg == 5){
		//青銅
		//GLfloat mat_specular[] = { 0.393548, 0.271906, 0.166721, 1.0 };
		GLfloat mat_specular[] = { 0.28, 0.21, 0.11, 1.0 };
		GLfloat mat_diffuse[] = { 0.714, 0.4384, 0.18144, 1.0 };
		GLfloat mat_ambient[] = { 0.2125, 0.1275, 0.054, 1.0 };
		//GLfloat mat_shininess = 25.6;
		GLfloat mat_shininess = 32.0;
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);
	}
	else if (cflg == 6){
		//トルコ石
		GLfloat mat_specular[] = { 0.297254, 0.30829, 0.306678, 1.0 };
		GLfloat mat_diffuse[] = { 0.396, 0.74151, 0.69102, 1.0 };
		GLfloat mat_ambient[] = { 0.1, 0.18725, 0.1745, 1.0 };
		GLfloat mat_shininess = 12.8;
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);
	}
	else if (cflg == 7){
		//真珠 （若干変えた）
		GLfloat mat_specular[] = { 0.31, 0.296648, 0.296648, 1.0 };
		GLfloat mat_diffuse[] = { 0.5, 0.429, 0.429, 1.0 };
		GLfloat mat_ambient[] = { 0.21, 0.20725, 0.20725, 1.0 };
		GLfloat mat_shininess = 10.24;
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);
	}
	else if (cflg == 8){
		//緑
		float cval[3];
		cval[0] = 60.0 / 255.0;
		cval[1] = 214.0 / 255.0;
		cval[2] = 29.0 / 255.0;

		GLfloat mat_specular[] = { 0.10, 0.10, 0.10, 1.0 };
		GLfloat mat_diffuse[] = { cval[0], cval[1], cval[2], 1.0 };
		GLfloat mat_ambient[] = { 0.0, 0.0, 0.0, 1.0 };
		GLfloat mat_shininess = 32.0;

		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);
	}


}


//プラスチックの色を定める
void set_plastic_color(float cval[3]){

	GLfloat mat_specular[] = { 0.10, 0.10, 0.10, 1.0 };
	GLfloat mat_diffuse[] = { cval[0], cval[1], cval[2], 1.0 };
	GLfloat mat_ambient[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat mat_shininess = 12.0;

	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);


}


//球の表示
void disp_sphere(float org[3], float rad, float zoom){

	float npole[3], spole[3];
	for (int i = 0; i < 3; i++){
		npole[i] = org[i];
		spole[i] = org[i];
	}
	npole[2] += rad;
	spole[2] -= rad;

	int divN1 = 15;
	int divN2 = 10;
	for (int i = 0; i < divN1; i++){
		float theta1 = 2.0*M_PI / (float)divN1*(float)i;
		float theta2 = 2.0*M_PI / (float)divN1*(float)(i + 1);

		float pset[2][3], cset[2][3];
		for (int j = 0; j < divN2; j++){
			if (j == 0){
				for (int k = 0; k < 3; k++){
					pset[0][k] = npole[k];
				}
			}
			else {
				for (int k = 0; k < 3; k++){
					pset[0][k] = cset[0][k];
					pset[1][k] = cset[1][k];
				}
			}
			if (j != divN2 - 1){
				float rho = M_PI / (float)divN2*(float)(j + 1);
				cset[0][0] = rad*cos(theta1)*sin(rho);
				cset[0][1] = rad*sin(theta1)*sin(rho);
				cset[0][2] = rad*cos(rho);
				cset[1][0] = rad*cos(theta2)*sin(rho);
				cset[1][1] = rad*sin(theta2)*sin(rho);
				cset[1][2] = rad*cos(rho);
				for (int k = 0; k < 3; k++){
					cset[0][k] += org[k];
					cset[1][k] += org[k];
				}
			}
			else {
				for (int k = 0; k < 3; k++){
					cset[0][k] = spole[k];
				}
			}

			float norm[3];
			float nthet = (theta1 + theta2) / 2.0;
			float nrho = M_PI / (float)divN2*(float)(j + 0.5);
			norm[0] = cos(nthet)*sin(nrho);
			norm[1] = sin(nthet)*sin(nrho);
			norm[2] = cos(nrho);

			glBegin(GL_POLYGON);
			glNormal3f(zoom*norm[0], zoom*norm[1], zoom*norm[2]);
			glVertex3fv(cset[0]);
			if (j != divN2 - 1) glVertex3fv(cset[1]);
			if (j != 0) glVertex3fv(pset[1]);
			glVertex3fv(pset[0]);
			glEnd();
		}
	}

}


//曲線ポリゴンの法線計算
void norm_for_curve_polygon(float fr1[4][3], float fr2[4][3], float nd[6][3]){

	for (int i = 0; i < 3; i++){
		nd[0][i] = (fr1[2][(i + 1) % 3] - fr1[0][(i + 1) % 3])*(fr1[3][(i + 2) % 3] - fr1[0][(i + 2) % 3])
			- (fr1[2][(i + 2) % 3] - fr1[0][(i + 2) % 3])*(fr1[3][(i + 1) % 3] - fr1[0][(i + 1) % 3]);
		nd[0][i] += (fr1[3][(i + 1) % 3] - fr1[0][(i + 1) % 3])*(fr1[1][(i + 2) % 3] - fr1[0][(i + 2) % 3])
			- (fr1[3][(i + 2) % 3] - fr1[0][(i + 2) % 3])*(fr1[1][(i + 1) % 3] - fr1[0][(i + 1) % 3]);
		nd[1][i] = (fr1[1][(i + 1) % 3] - fr1[0][(i + 1) % 3])*(fr2[1][(i + 2) % 3] - fr1[0][(i + 2) % 3])
			- (fr1[1][(i + 2) % 3] - fr1[0][(i + 2) % 3])*(fr2[1][(i + 1) % 3] - fr1[0][(i + 1) % 3]);
		nd[1][i] += (fr2[1][(i + 1) % 3] - fr1[0][(i + 1) % 3])*(fr2[0][(i + 2) % 3] - fr1[0][(i + 2) % 3])
			- (fr2[1][(i + 2) % 3] - fr1[0][(i + 2) % 3])*(fr2[0][(i + 1) % 3] - fr1[0][(i + 1) % 3]);
		nd[2][i] = (fr1[3][(i + 1) % 3] - fr1[1][(i + 1) % 3])*(fr2[3][(i + 2) % 3] - fr1[1][(i + 2) % 3])
			- (fr1[3][(i + 2) % 3] - fr1[1][(i + 2) % 3])*(fr2[3][(i + 1) % 3] - fr1[1][(i + 1) % 3]);
		nd[2][i] += (fr2[3][(i + 1) % 3] - fr1[1][(i + 1) % 3])*(fr2[1][(i + 2) % 3] - fr1[1][(i + 2) % 3])
			- (fr2[3][(i + 2) % 3] - fr1[1][(i + 2) % 3])*(fr2[1][(i + 1) % 3] - fr1[1][(i + 1) % 3]);
		nd[3][i] = (fr1[2][(i + 1) % 3] - fr1[3][(i + 1) % 3])*(fr2[2][(i + 2) % 3] - fr1[3][(i + 2) % 3])
			- (fr1[2][(i + 2) % 3] - fr1[3][(i + 2) % 3])*(fr2[2][(i + 1) % 3] - fr1[3][(i + 1) % 3]);
		nd[3][i] += (fr2[2][(i + 1) % 3] - fr1[3][(i + 1) % 3])*(fr2[3][(i + 2) % 3] - fr1[3][(i + 2) % 3])
			- (fr2[2][(i + 2) % 3] - fr1[3][(i + 2) % 3])*(fr2[3][(i + 1) % 3] - fr1[3][(i + 1) % 3]);
		nd[4][i] = (fr1[0][(i + 1) % 3] - fr1[2][(i + 1) % 3])*(fr2[0][(i + 2) % 3] - fr1[2][(i + 2) % 3])
			- (fr1[0][(i + 2) % 3] - fr1[2][(i + 2) % 3])*(fr2[0][(i + 1) % 3] - fr1[2][(i + 1) % 3]);
		nd[4][i] += (fr2[0][(i + 1) % 3] - fr1[2][(i + 1) % 3])*(fr2[2][(i + 2) % 3] - fr1[2][(i + 2) % 3])
			- (fr2[0][(i + 2) % 3] - fr1[2][(i + 2) % 3])*(fr2[2][(i + 1) % 3] - fr1[2][(i + 1) % 3]);
		nd[5][i] = (fr2[1][(i + 1) % 3] - fr2[0][(i + 1) % 3])*(fr2[3][(i + 2) % 3] - fr2[0][(i + 2) % 3])
			- (fr2[1][(i + 2) % 3] - fr2[0][(i + 2) % 3])*(fr2[3][(i + 1) % 3] - fr2[0][(i + 1) % 3]);
		nd[5][i] += (fr2[3][(i + 1) % 3] - fr2[0][(i + 1) % 3])*(fr2[2][(i + 2) % 3] - fr2[0][(i + 2) % 3])
			- (fr2[3][(i + 2) % 3] - fr2[0][(i + 2) % 3])*(fr2[2][(i + 1) % 3] - fr2[0][(i + 1) % 3]);
	}

	for (int i = 0; i < 6; i++){
		float ndh = 0.0;
		for (int j = 0; j < 3; j++){
			ndh += nd[i][j] * nd[i][j];
		}
		ndh = sqrt(ndh);
		for (int j = 0; j < 3; j++){
			nd[i][j] /= ndh;
		}
	}


}


//曲線のポリゴン表示
void disp_curve_polygon(class curve_on_mesh *curv, float zoom,float off){


	/*glPolygonMode(GL_FRONT, GL_FILL);
	glShadeModel(GL_SMOOTH);*/

	float frame1[4][3], frame2[4][3], frame3[4][3];
	for (int i = 0; i < curv->nodeN - 1; i++){
		if (i == 0){
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 2; k++){
					for (int l = 0; l < 2; l++){
						frame1[2 * k + l][j] = curv->cNode[i].pos[j] + off*curv->Mframe[i].d2[j];
						if (l == 0) frame1[2 * k + l][j] -= curv->sParam->Wid*curv->Mframe[i].d1[j] / 2.0;
						else frame1[2 * k + l][j] += curv->sParam->Wid*curv->Mframe[i].d1[j] / 2.0;
						if (k == 0) frame1[2 * k + l][j] -= curv->sParam->Hei*curv->Mframe[i].d2[j] / 2.0;
						else frame1[2 * k + l][j] += curv->sParam->Hei*curv->Mframe[i].d2[j] / 2.0;
					}
				}
			}

		}
		else {
			for (int j = 0; j < 4; j++){
				for (int k = 0; k < 3; k++){
					frame1[j][k] = frame3[j][k];
				}
			}
		}

		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 2; l++){
					frame2[2 * k + l][j] = 0.5*(curv->cNode[i].pos[j] + curv->cNode[i + 1].pos[j]) + off*curv->Mframe[i].d2[j];
					if (l == 0) frame2[2 * k + l][j] -= curv->sParam->Wid*curv->Mframe[i].d1[j] / 2.0;
					else frame2[2 * k + l][j] += curv->sParam->Wid*curv->Mframe[i].d1[j] / 2.0;
					if (k == 0) frame2[2 * k + l][j] -= curv->sParam->Hei*curv->Mframe[i].d2[j] / 2.0;
					else frame2[2 * k + l][j] += curv->sParam->Hei*curv->Mframe[i].d2[j] / 2.0;
				}
			}
		}

		if (i != curv->nodeN - 2){
			float axis[3];
			float axh = 0.0;
			float nseki = 0.0;
			for (int j = 0; j < 3; j++){
				axis[j] = curv->Mframe[i].d3[(j + 1) % 3] * curv->Mframe[i + 1].d3[(j + 2) % 3]
					- curv->Mframe[i].d3[(j + 2) % 3] * curv->Mframe[i + 1].d3[(j + 1) % 3];
				axh += axis[j] * axis[j];
				nseki += curv->Mframe[i].d3[j] * curv->Mframe[i + 1].d3[j];
			}
			axh = sqrt(axh);

			float quat1[4], Rot1[3][3];
			float quat2[4], Rot2[3][3];
			float r1d[3][3], r2d[3];
			if (axh > sin(1.0e-5)){
				if (nseki > 1.0) nseki = 1.0;
				if (nseki < -1.0) nseki = -1.0;
				nseki = acos(nseki);
				quat1[0] = cos(nseki / 4.0);
				quat2[0] = cos(nseki / 2.0);
				for (int j = 0; j < 3; j++){
					quat1[1 + j] = sin(nseki / 4.0)*axis[j] / axh;
					quat2[1 + j] = sin(nseki / 2.0)*axis[j] / axh;
				}
				derive_rot_mat_from_quaternion(quat1, Rot1);
				derive_rot_mat_from_quaternion(quat2, Rot2);

				for (int j = 0; j < 3; j++){
					r1d[0][j] = 0.0;
					r1d[1][j] = 0.0;
					r1d[2][j] = 0.0;
					r2d[j] = 0.0;
					for (int k = 0; k < 3; k++){
						r1d[0][j] += Rot1[j][k] * curv->Mframe[i].d1[k];
						r1d[1][j] += Rot1[j][k] * curv->Mframe[i].d2[k];
						r1d[2][j] += Rot1[j][k] * curv->Mframe[i].d3[k];
						r2d[j] += Rot2[j][k] * curv->Mframe[i].d1[k];
					}
				}
			}
			else {
				for (int j = 0; j < 3; j++){
					r1d[0][j] = curv->Mframe[i].d1[j];
					r1d[1][j] = curv->Mframe[i].d2[j];
					r1d[2][j] = curv->Mframe[i].d3[j];
					r2d[j] = curv->Mframe[i].d1[j];
				}
			}
			nseki = 0.0;
			float oseki[3];
			for (int j = 0; j < 3; j++){
				oseki[j] = r2d[(j + 1) % 3] * curv->Mframe[i + 1].d1[(j + 2) % 3]
					- r2d[(j + 2) % 3] * curv->Mframe[i + 1].d2[(j + 1) % 3];
				nseki += r2d[j] * curv->Mframe[i + 1].d1[j];
			}
			float sig = 0.0;
			for (int j = 0; j < 3; j++){
				sig += oseki[j] * curv->Mframe[i + 1].d3[j];
			}
			if (nseki > 1.0) nseki = 1.0;
			if (nseki < -1.0) nseki = -1.0;
			nseki = acos(nseki);
			if (sig < 0.0)nseki = -nseki;
			quat1[0] = cos(nseki / 4.0);
			for (int j = 0; j < 3; j++){
				quat1[1 + j] = sin(nseki / 4.0) * r1d[2][j];
			}
			derive_rot_mat_from_quaternion(quat1, Rot1);

			float dire[2][3];
			for (int j = 0; j < 3; j++){
				dire[0][j] = 0.0;
				dire[1][j] = 0.0;
				for (int k = 0; k < 3; k++){
					dire[0][j] += Rot1[j][k] * r1d[0][k];
					dire[1][j] += Rot1[j][k] * r1d[1][k];
				}
			}

			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 2; k++){
					for (int l = 0; l < 2; l++){
						frame3[2 * k + l][j] = curv->cNode[i + 1].pos[j]
							+ off*0.5*(curv->Mframe[i].d2[j] + curv->Mframe[i + 1].d2[j]);
						if (l == 0) frame3[2 * k + l][j] -= curv->sParam->Wid*dire[0][j] / 2.0;
						else frame3[2 * k + l][j] += curv->sParam->Wid*dire[0][j] / 2.0;
						if (k == 0) frame3[2 * k + l][j] -= curv->sParam->Hei*dire[1][j] / 2.0;
						else frame3[2 * k + l][j] += curv->sParam->Hei*dire[1][j] / 2.0;
					}
				}
			}

		}
		else {
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 2; k++){
					for (int l = 0; l < 2; l++){
						frame3[2 * k + l][j] = curv->cNode[i + 1].pos[j] + off*curv->Mframe[i].d2[j];
						if (l == 0) frame3[2 * k + l][j] -= curv->sParam->Wid*curv->Mframe[i].d1[j] / 2.0;
						else frame3[2 * k + l][j] += curv->sParam->Wid*curv->Mframe[i].d1[j] / 2.0;
						if (k == 0) frame3[2 * k + l][j] -= curv->sParam->Hei*curv->Mframe[i].d2[j] / 2.0;
						else frame3[2 * k + l][j] += curv->sParam->Hei*curv->Mframe[i].d2[j] / 2.0;
					}
				}
			}
		}

		/*for (int j = 0; j < 4; j++){
		if (i == 0) node_scaling(frame1[j]);
		node_scaling(frame2[j]);
		node_scaling(frame3[j]);
		}*/

		float norm[6][3];
		norm_for_curve_polygon(frame1, frame2, norm);
		for (int j = 0; j < 6; j++){
			for (int k = 0; k < 3; k++){
				norm[j][k] *= zoom;
			}
		}

		glBegin(GL_POLYGON);
		glNormal3fv(norm[0]);
		glVertex3fv(frame1[0]);
		glVertex3fv(frame1[2]);
		glVertex3fv(frame1[3]);
		glVertex3fv(frame1[1]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[1]);
		glVertex3fv(frame1[0]);
		glVertex3fv(frame1[1]);
		glVertex3fv(frame2[1]);
		glVertex3fv(frame2[0]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[2]);
		glVertex3fv(frame1[1]);
		glVertex3fv(frame1[3]);
		glVertex3fv(frame2[3]);
		glVertex3fv(frame2[1]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[3]);
		glVertex3fv(frame1[3]);
		glVertex3fv(frame1[2]);
		glVertex3fv(frame2[2]);
		glVertex3fv(frame2[3]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[4]);
		glVertex3fv(frame1[2]);
		glVertex3fv(frame1[0]);
		glVertex3fv(frame2[0]);
		glVertex3fv(frame2[2]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[5]);
		glVertex3fv(frame2[0]);
		glVertex3fv(frame2[1]);
		glVertex3fv(frame2[3]);
		glVertex3fv(frame2[2]);
		glEnd();


		norm_for_curve_polygon(frame2, frame3, norm);
		for (int j = 0; j < 6; j++){
			for (int k = 0; k < 3; k++){
				norm[j][k] *= zoom;
			}
		}

		glBegin(GL_POLYGON);
		glNormal3fv(norm[0]);
		glVertex3fv(frame2[0]);
		glVertex3fv(frame2[2]);
		glVertex3fv(frame2[3]);
		glVertex3fv(frame2[1]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[1]);
		glVertex3fv(frame2[0]);
		glVertex3fv(frame2[1]);
		glVertex3fv(frame3[1]);
		glVertex3fv(frame3[0]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[2]);
		glVertex3fv(frame2[1]);
		glVertex3fv(frame2[3]);
		glVertex3fv(frame3[3]);
		glVertex3fv(frame3[1]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[3]);
		glVertex3fv(frame2[3]);
		glVertex3fv(frame2[2]);
		glVertex3fv(frame3[2]);
		glVertex3fv(frame3[3]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[4]);
		glVertex3fv(frame2[2]);
		glVertex3fv(frame2[0]);
		glVertex3fv(frame3[0]);
		glVertex3fv(frame3[2]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[5]);
		glVertex3fv(frame3[0]);
		glVertex3fv(frame3[1]);
		glVertex3fv(frame3[3]);
		glVertex3fv(frame3[2]);
		glEnd();
	}

}


//ロッドのポリゴン表示
void disp_rod_polygon(class DiscreteRod *rod, float zoom){

	float frame1[4][3], frame2[4][3], frame3[4][3];
	for (int i = 0; i < rod->nodeN - 1; i++){
		if (i == 0){
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 2; k++){
					for (int l = 0; l < 2; l++){
						frame1[2 * k + l][j] = rod->Nodes[i][j];
						if (l == 0) frame1[2 * k + l][j] -= rod->sParam->Wid*rod->Mframe[i].d1[j] / 2.0;
						else frame1[2 * k + l][j] += rod->sParam->Wid*rod->Mframe[i].d1[j] / 2.0;
						if (k == 0) frame1[2 * k + l][j] -= rod->sParam->Hei*rod->Mframe[i].d2[j] / 2.0;
						else frame1[2 * k + l][j] += rod->sParam->Hei*rod->Mframe[i].d2[j] / 2.0;
					}
				}
			}

		}
		else {
			for (int j = 0; j < 4; j++){
				for (int k = 0; k < 3; k++){
					frame1[j][k] = frame3[j][k];
				}
			}
		}

		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 2; l++){
					frame2[2 * k + l][j] = 0.5*(rod->Nodes[i][j] + rod->Nodes[i + 1][j]);
					if (l == 0) frame2[2 * k + l][j] -= rod->sParam->Wid*rod->Mframe[i].d1[j] / 2.0;
					else frame2[2 * k + l][j] += rod->sParam->Wid*rod->Mframe[i].d1[j] / 2.0;
					if (k == 0) frame2[2 * k + l][j] -= rod->sParam->Hei*rod->Mframe[i].d2[j] / 2.0;
					else frame2[2 * k + l][j] += rod->sParam->Hei*rod->Mframe[i].d2[j] / 2.0;
				}
			}
		}

		if (i != rod->nodeN - 2){
			float axis[3];
			float axh = 0.0;
			float nseki = 0.0;
			for (int j = 0; j < 3; j++){
				axis[j] = rod->Mframe[i].d3[(j + 1) % 3] * rod->Mframe[i + 1].d3[(j + 2) % 3]
					- rod->Mframe[i].d3[(j + 2) % 3] * rod->Mframe[i + 1].d3[(j + 1) % 3];
				axh += axis[j] * axis[j];
				nseki += rod->Mframe[i].d3[j] * rod->Mframe[i + 1].d3[j];
			}
			axh = sqrt(axh);

			float quat1[4], Rot1[3][3];
			float quat2[4], Rot2[3][3];
			float r1d[3][3], r2d[3];
			if (axh > sin(1.0e-5)){
				if (nseki > 1.0) nseki = 1.0;
				if (nseki < -1.0) nseki = -1.0;
				nseki = acos(nseki);
				quat1[0] = cos(nseki / 4.0);
				quat2[0] = cos(nseki / 2.0);
				for (int j = 0; j < 3; j++){
					quat1[1 + j] = sin(nseki / 4.0)*axis[j] / axh;
					quat2[1 + j] = sin(nseki / 2.0)*axis[j] / axh;
				}
				derive_rot_mat_from_quaternion(quat1, Rot1);
				derive_rot_mat_from_quaternion(quat2, Rot2);

				for (int j = 0; j < 3; j++){
					r1d[0][j] = 0.0;
					r1d[1][j] = 0.0;
					r1d[2][j] = 0.0;
					r2d[j] = 0.0;
					for (int k = 0; k < 3; k++){
						r1d[0][j] += Rot1[j][k] * rod->Mframe[i].d1[k];
						r1d[1][j] += Rot1[j][k] * rod->Mframe[i].d2[k];
						r1d[2][j] += Rot1[j][k] * rod->Mframe[i].d3[k];
						r2d[j] += Rot2[j][k] * rod->Mframe[i].d1[k];
					}
				}
			}
			else {
				for (int j = 0; j < 3; j++){
					r1d[0][j] = rod->Mframe[i].d1[j];
					r1d[1][j] = rod->Mframe[i].d2[j];
					r1d[2][j] = rod->Mframe[i].d3[j];
					r2d[j] = rod->Mframe[i].d1[j];
				}
			}
			nseki = 0.0;
			float oseki[3];
			for (int j = 0; j < 3; j++){
				oseki[j] = r2d[(j + 1) % 3] * rod->Mframe[i + 1].d1[(j + 2) % 3]
					- r2d[(j + 2) % 3] * rod->Mframe[i + 1].d2[(j + 1) % 3];
				nseki += r2d[j] * rod->Mframe[i + 1].d1[j];
			}
			float sig = 0.0;
			for (int j = 0; j < 3; j++){
				sig += oseki[j] * rod->Mframe[i + 1].d3[j];
			}
			if (nseki > 1.0) nseki = 1.0;
			if (nseki < -1.0) nseki = -1.0;
			nseki = acos(nseki);
			if (sig < 0.0)nseki = -nseki;
			quat1[0] = cos(nseki / 4.0);
			for (int j = 0; j < 3; j++){
				quat1[1 + j] = sin(nseki / 4.0) * r1d[2][j];
			}
			derive_rot_mat_from_quaternion(quat1, Rot1);

			float dire[2][3];
			for (int j = 0; j < 3; j++){
				dire[0][j] = 0.0;
				dire[1][j] = 0.0;
				for (int k = 0; k < 3; k++){
					dire[0][j] += Rot1[j][k] * r1d[0][k];
					dire[1][j] += Rot1[j][k] * r1d[1][k];
				}
			}

			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 2; k++){
					for (int l = 0; l < 2; l++){
						frame3[2 * k + l][j] = rod->Nodes[i + 1][j];
						if (l == 0) frame3[2 * k + l][j] -= rod->sParam->Wid*dire[0][j] / 2.0;
						else frame3[2 * k + l][j] += rod->sParam->Wid*dire[0][j] / 2.0;
						if (k == 0) frame3[2 * k + l][j] -= rod->sParam->Hei*dire[1][j] / 2.0;
						else frame3[2 * k + l][j] += rod->sParam->Hei*dire[1][j] / 2.0;
					}
				}
			}

		}
		else {
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 2; k++){
					for (int l = 0; l < 2; l++){
						frame3[2 * k + l][j] = rod->Nodes[i + 1][j];
						if (l == 0) frame3[2 * k + l][j] -= rod->sParam->Wid*rod->Mframe[i].d1[j] / 2.0;
						else frame3[2 * k + l][j] += rod->sParam->Wid*rod->Mframe[i].d1[j] / 2.0;
						if (k == 0) frame3[2 * k + l][j] -= rod->sParam->Hei*rod->Mframe[i].d2[j] / 2.0;
						else frame3[2 * k + l][j] += rod->sParam->Hei*rod->Mframe[i].d2[j] / 2.0;
					}
				}
			}
		}

		/*for (int j = 0; j < 4; j++){
		if (i == 0) node_scaling(frame1[j]);
		node_scaling(frame2[j]);
		node_scaling(frame3[j]);
		}*/

		float norm[6][3];
		norm_for_curve_polygon(frame1, frame2, norm);
		for (int j = 0; j < 6; j++){
			for (int k = 0; k < 3; k++){
				norm[j][k] *= zoom;
			}
		}

		glBegin(GL_POLYGON);
		glNormal3fv(norm[0]);
		glVertex3fv(frame1[0]);
		glVertex3fv(frame1[2]);
		glVertex3fv(frame1[3]);
		glVertex3fv(frame1[1]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[1]);
		glVertex3fv(frame1[0]);
		glVertex3fv(frame1[1]);
		glVertex3fv(frame2[1]);
		glVertex3fv(frame2[0]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[2]);
		glVertex3fv(frame1[1]);
		glVertex3fv(frame1[3]);
		glVertex3fv(frame2[3]);
		glVertex3fv(frame2[1]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[3]);
		glVertex3fv(frame1[3]);
		glVertex3fv(frame1[2]);
		glVertex3fv(frame2[2]);
		glVertex3fv(frame2[3]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[4]);
		glVertex3fv(frame1[2]);
		glVertex3fv(frame1[0]);
		glVertex3fv(frame2[0]);
		glVertex3fv(frame2[2]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[5]);
		glVertex3fv(frame2[0]);
		glVertex3fv(frame2[1]);
		glVertex3fv(frame2[3]);
		glVertex3fv(frame2[2]);
		glEnd();


		norm_for_curve_polygon(frame2, frame3, norm);
		for (int j = 0; j < 6; j++){
			for (int k = 0; k < 3; k++){
				norm[j][k] *= zoom;
			}
		}

		glBegin(GL_POLYGON);
		glNormal3fv(norm[0]);
		glVertex3fv(frame2[0]);
		glVertex3fv(frame2[2]);
		glVertex3fv(frame2[3]);
		glVertex3fv(frame2[1]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[1]);
		glVertex3fv(frame2[0]);
		glVertex3fv(frame2[1]);
		glVertex3fv(frame3[1]);
		glVertex3fv(frame3[0]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[2]);
		glVertex3fv(frame2[1]);
		glVertex3fv(frame2[3]);
		glVertex3fv(frame3[3]);
		glVertex3fv(frame3[1]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[3]);
		glVertex3fv(frame2[3]);
		glVertex3fv(frame2[2]);
		glVertex3fv(frame3[2]);
		glVertex3fv(frame3[3]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[4]);
		glVertex3fv(frame2[2]);
		glVertex3fv(frame2[0]);
		glVertex3fv(frame3[0]);
		glVertex3fv(frame3[2]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3fv(norm[5]);
		glVertex3fv(frame3[0]);
		glVertex3fv(frame3[1]);
		glVertex3fv(frame3[3]);
		glVertex3fv(frame3[2]);
		glEnd();
	}


}


//矢印の表示
void disp_arrow(float end1[3], float end2[3], float rad, float zoom){

	//サイズ設定
	float cent[3], ed[3];
	float edh = 0.0;
	for (int i = 0; i < 3; i++){
		cent[i] = (end1[i] + end2[i]) / 2.0;
		ed[i] = end2[i] - end1[i];
		edh += ed[i] * ed[i];
	}
	edh = sqrt(edh);
	for (int i = 0; i < 3; i++){
		ed[i] /= edh;
	}

	float Cwid = rad*1.2;
	float Chei = Cwid*0.85;
	if (edh < Chei*3.0){
		edh = Chei*3.0;
		for (int i = 0; i < 3; i++){
			end1[i] = cent[i] - edh*ed[i] / 2.0;
			end2[i] = cent[i] + edh*ed[i] / 2.0;
		}
	}
	float cyW = rad / 2.0;
	float cyH = edh - 2.0*Chei;

	//基底ベクトル設定
	float bas1[3], bas2[2];
	int maxN = 0;
	for (int i = 1; i < 3; i++){
		if (fabs(ed[i]) > fabs(ed[maxN])) maxN = i;
	}
	float bv = 0.0;
	for (int i = 0; i < 3; i++){
		if (i != maxN) {
			bas1[i] = ed[i];
			bv -= ed[i] * ed[i];
		}
	}
	bas1[maxN] = bv / ed[maxN];
	float seki = ed[0] * bas1[0] + ed[1] * bas1[1] + ed[2] * bas1[2];
	float bh = 0.0;
	for (int i = 0; i < 3; i++){
		bas1[i] -= seki*ed[i];
		bh += bas1[i] * bas1[i];
	}
	bh = sqrt(bh);
	for (int i = 0; i < 3; i++){
		bas1[i] /= bh;
	}
	seki = 0.0;
	for (int i = 0; i < 3; i++){
		bas2[i] = ed[(i + 1) % 3] * bas1[(i + 2) % 3] - ed[(i + 2) % 3] * bas1[(i + 1) % 3];
		seki += bas2[i] * bas2[i];
	}
	seki = sqrt(seki);
	for (int i = 0; i < 3; i++){
		bas2[i] /= seki;
	}

	//描画
	int divN = 10.0;
	for (int i = 0; i < divN; i++){
		float theta1 = 2.0*M_PI / (float)divN*(float)i;
		float theta2 = 2.0*M_PI / (float)divN*(float)(i + 1.0);
		float theta3 = (theta1 + theta2) / 2.0;
		float rect[4][3], norm[3];
		for (int j = 0; j < 2; j++){
			for (int k = 0; k < 3; k++){
				rect[j][k] = cent[k] + cyH / 2.0*ed[k];
				rect[j + 2][k] = cent[k] - cyH / 2.0*ed[k];
			}
		}
		float nh = 0.0;
		for (int j = 0; j < 3; j++){
			rect[0][j] += cyW*cos(theta1)*bas1[j] + cyW*sin(theta1)*bas2[j];
			rect[1][j] += cyW*cos(theta2)*bas1[j] + cyW*sin(theta2)*bas2[j];
			rect[2][j] += cyW*cos(theta1)*bas1[j] + cyW*sin(theta1)*bas2[j];
			rect[3][j] += cyW*cos(theta2)*bas1[j] + cyW*sin(theta2)*bas2[j];
			norm[j] = cos(theta3)*bas1[j] + sin(theta3)*bas2[j];
			nh += norm[j] * norm[j];
		}
		nh = sqrt(nh);
		for (int j = 0; j < 3; j++){
			norm[j] /= nh / zoom;
		}

		glNormal3fv(norm);
		glBegin(GL_POLYGON);
		glVertex3fv(rect[0]);
		glVertex3fv(rect[1]);
		glVertex3fv(rect[3]);
		glVertex3fv(rect[2]);
		glEnd();

		float cirp[4][3];
		for (int j = 0; j < 2; j++){
			for (int k = 0; k < 3; k++){
				rect[j][k] = cent[k] + cyH / 2.0*ed[k];
				rect[j + 2][k] = cent[k] - cyH / 2.0*ed[k];
			}
		}
		float norm1[3], norm2[3];
		float nh1 = 0.0;
		float nh2 = 0.0;
		for (int j = 0; j < 3; j++){
			rect[0][j] += Cwid*cos(theta1)*bas1[j] + Cwid*sin(theta1)*bas2[j];
			rect[1][j] += Cwid*cos(theta2)*bas1[j] + Cwid*sin(theta2)*bas2[j];
			rect[2][j] += Cwid*cos(theta1)*bas1[j] + Cwid*sin(theta1)*bas2[j];
			rect[3][j] += Cwid*cos(theta2)*bas1[j] + Cwid*sin(theta2)*bas2[j];
			cirp[0][j] = cent[j] + cyH / 2.0*ed[j];
			cirp[1][j] = cent[j] - cyH / 2.0*ed[j];
			cirp[2][j] = cirp[0][j] + Chei*ed[j];
			cirp[3][j] = cirp[1][j] - Chei*ed[j];
		}
		for (int j = 0; j < 3; j++){
			norm1[j] = (rect[1][(j + 1) % 3] - rect[0][(j + 1) % 3])*(cirp[2][(j + 2) % 3] - rect[0][(j + 2) % 3])
				- (rect[1][(j + 2) % 3] - rect[0][(j + 2) % 3])*(cirp[2][(j + 1) % 3] - rect[0][(j + 1) % 3]);
			norm2[j] = (rect[2][(j + 1) % 3] - rect[3][(j + 1) % 3])*(cirp[3][(j + 2) % 3] - rect[3][(j + 2) % 3])
				- (rect[2][(j + 2) % 3] - rect[3][(j + 2) % 3])*(cirp[3][(j + 1) % 3] - rect[3][(j + 1) % 3]);
			nh1 += norm1[j] * norm1[j];
			nh2 += norm2[j] * norm2[j];
		}
		nh1 = sqrt(nh1);
		nh2 = sqrt(nh2);
		for (int j = 0; j < 3; j++){
			norm1[j] /= nh1;
			norm2[j] /= nh2;
		}

		glNormal3f(-ed[0] * zoom, -ed[1] * zoom, -ed[2] * zoom);
		glBegin(GL_POLYGON);
		glVertex3fv(cirp[0]);
		glVertex3fv(rect[1]);
		glVertex3fv(rect[0]);
		glEnd();

		glNormal3f(ed[0] * zoom, ed[1] * zoom, ed[2] * zoom);
		glBegin(GL_POLYGON);
		glVertex3fv(cirp[1]);
		glVertex3fv(rect[2]);
		glVertex3fv(rect[3]);
		glEnd();

		glNormal3f(norm1[0] * zoom, norm1[1] * zoom, norm1[2] * zoom);
		glBegin(GL_POLYGON);
		glVertex3fv(rect[1]);
		glVertex3fv(rect[0]);
		glVertex3fv(cirp[2]);
		glEnd();

		glNormal3f(norm2[0] * zoom, norm2[1] * zoom, norm2[2] * zoom);
		glBegin(GL_POLYGON);
		glVertex3fv(rect[2]);
		glVertex3fv(rect[3]);
		glVertex3fv(cirp[3]);
		glEnd();

	}

}


//コネクションの表示
void disp_connection_polygon(float end1[3], float end2[3], float norm1[3], float norm2[3], float rad, float hei, float zoom){

	float bas1[3][3], bas2[3][3];
	float seki1 = 0.0;
	float seki2 = 0.0;
	float edh = 0.0;
	for (int i = 0; i < 3; i++){
		bas1[2][i] = norm1[i];
		bas2[2][i] = norm2[i];
		seki1 += (end2[i] - end1[i])*norm1[i];
		seki2 += (end2[i] - end1[i])*norm2[i];
		edh += pow(end2[i] - end1[i], 2);
	}
	edh = sqrt(edh);
	float bh1 = 0.0;
	float bh2 = 0.0;
	for (int i = 0; i < 3; i++){
		bas1[1][i] = end2[i] - end1[i] - seki1*norm1[i];
		bas2[1][i] = end2[i] - end1[i] - seki2*norm2[i];
		bh1 += bas1[1][i] * bas1[1][i];
		bh2 += bas2[1][i] * bas2[1][i];
	}
	bh1 = sqrt(bh1);
	bh2 = sqrt(bh2);
	if (bh1 < 1.0e-10) bh1 = 1.0e-10;
	if (bh2 < 1.0e-10) bh2 = 1.0e-10;
	for (int i = 0; i < 3; i++){
		bas1[1][i] /= bh1;
		bas2[1][i] /= bh2;
	}
	bh1 = 0.0;
	bh2 = 0.0;
	for (int i = 0; i < 3; i++){
		bas1[0][i] = bas1[1][(i + 1) % 3] * bas1[2][(i + 2) % 3] - bas1[1][(i + 2) % 3] * bas1[2][(i + 1) % 3];
		bas2[0][i] = bas2[1][(i + 1) % 3] * bas2[2][(i + 2) % 3] - bas2[1][(i + 2) % 3] * bas2[2][(i + 1) % 3];
		bh1 += bas1[0][i] * bas1[0][i];
		bh2 += bas2[0][i] * bas2[0][i];
	}
	bh1 = sqrt(bh1);
	bh2 = sqrt(bh2);
	if (bh1 < 1.0e-10) bh1 = 1.0e-10;
	if (bh2 < 1.0e-10) bh2 = 1.0e-10;
	for (int i = 0; i < 3; i++){
		bas1[0][i] /= bh1;
		bas2[0][i] /= bh2;
	}


	//半円の作成
	int divN = 10;
	for (int i = 0; i < divN; i++){
		float theta1 = M_PI / (float)divN *(float)i;
		float theta2 = M_PI / (float)divN*(float)(i + 1.0);
		float theta3 = (theta1 + theta2) / 2.0;
		float cir1[4][3], cir2[4][3];
		float cent1[2][3], cent2[2][3];
		float Rnorm[2][3];
		for (int j = 0; j < 3; j++){
			cir1[0][j] = end1[j] + rad*cos(theta1)*bas1[0][j] - rad*sin(theta1)*bas1[1][j] + hei*bas1[2][j];
			cir1[1][j] = end1[j] + rad*cos(theta2)*bas1[0][j] - rad*sin(theta2)*bas1[1][j] + hei*bas1[2][j];
			cir1[2][j] = end1[j] + rad*cos(theta1)*bas1[0][j] - rad*sin(theta1)*bas1[1][j] + 0.0*bas1[2][j];
			cir1[3][j] = end1[j] + rad*cos(theta2)*bas1[0][j] - rad*sin(theta2)*bas1[1][j] + 0.0*bas1[2][j];
			cir2[0][j] = end2[j] + rad*cos(theta1)*bas2[0][j] + rad*sin(theta1)*bas2[1][j] + hei*bas2[2][j];
			cir2[1][j] = end2[j] + rad*cos(theta2)*bas2[0][j] + rad*sin(theta2)*bas2[1][j] + hei*bas2[2][j];
			cir2[2][j] = end2[j] + rad*cos(theta1)*bas2[0][j] + rad*sin(theta1)*bas2[1][j] + 0.0*bas2[2][j];
			cir2[3][j] = end2[j] + rad*cos(theta2)*bas2[0][j] + rad*sin(theta2)*bas2[1][j] + 0.0*bas2[2][j];
			cent1[0][j] = end1[j] + hei*bas1[2][j];
			cent1[1][j] = end1[j];
			cent2[0][j] = end2[j] + hei*bas2[2][j];
			cent2[1][j] = end2[j];
			Rnorm[0][j] = cos(theta3)*bas1[0][j] - sin(theta3)*bas1[1][j];
			Rnorm[1][j] = cos(theta3)*bas2[0][j] + sin(theta3)*bas2[1][j];
		}

		glBegin(GL_POLYGON);
		glNormal3f(bas1[2][0] * zoom, bas1[2][1] * zoom, bas1[2][2] * zoom);
		glVertex3fv(cir1[1]);
		glVertex3fv(cir1[0]);
		glVertex3fv(cent1[0]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3f(-bas1[2][0] * zoom, -bas1[2][1] * zoom, -bas1[2][2] * zoom);
		glVertex3fv(cir1[2]);
		glVertex3fv(cir1[3]);
		glVertex3fv(cent1[1]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3f(Rnorm[0][0] * zoom, Rnorm[0][1] * zoom, Rnorm[0][2] * zoom);
		glVertex3fv(cir1[0]);
		glVertex3fv(cir1[1]);
		glVertex3fv(cir1[3]);
		glVertex3fv(cir1[2]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3f(bas2[2][0] * zoom, bas2[2][1] * zoom, bas2[2][2] * zoom);
		glVertex3fv(cir2[0]);
		glVertex3fv(cir2[1]);
		glVertex3fv(cent2[0]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3f(-bas2[2][0] * zoom, -bas2[2][1] * zoom, -bas2[2][2] * zoom);
		glVertex3fv(cir2[3]);
		glVertex3fv(cir2[2]);
		glVertex3fv(cent2[1]);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3f(Rnorm[1][0] * zoom, Rnorm[1][1] * zoom, Rnorm[1][2] * zoom);
		glVertex3fv(cir2[0]);
		glVertex3fv(cir2[2]);
		glVertex3fv(cir2[3]);
		glVertex3fv(cir2[1]);
		glEnd();

	}

	//接続部の作成
	float sumi1[4][3], sumi2[4][3];
	for (int i = 0; i < 3; i++){
		sumi1[0][i] = end1[i] - rad*bas1[0][i] + 0.0*bas1[2][i];
		sumi1[1][i] = end1[i] - rad*bas1[0][i] + hei*bas1[2][i];
		sumi1[2][i] = end1[i] + rad*bas1[0][i] + 0.0*bas1[2][i];
		sumi1[3][i] = end1[i] + rad*bas1[0][i] + hei*bas1[2][i];
		sumi2[0][i] = end2[i] - rad*bas2[0][i] + 0.0*bas2[2][i];
		sumi2[1][i] = end2[i] - rad*bas2[0][i] + hei*bas2[2][i];
		sumi2[2][i] = end2[i] + rad*bas2[0][i] + 0.0*bas2[2][i];
		sumi2[3][i] = end2[i] + rad*bas2[0][i] + hei*bas2[2][i];
	}
	float norm[6][3];
	norm_for_curve_polygon(sumi1, sumi2, norm);

	glBegin(GL_POLYGON);
	glNormal3f(norm[0][0] * zoom, norm[0][1] * zoom, norm[0][2] * zoom);
	glVertex3fv(sumi1[0]);
	glVertex3fv(sumi1[2]);
	glVertex3fv(sumi1[3]);
	glVertex3fv(sumi1[1]);
	glEnd();

	glBegin(GL_POLYGON);
	glNormal3f(norm[1][0] * zoom, norm[1][1] * zoom, norm[1][2] * zoom);
	glVertex3fv(sumi1[0]);
	glVertex3fv(sumi1[1]);
	glVertex3fv(sumi2[1]);
	glVertex3fv(sumi2[0]);
	glEnd();

	glBegin(GL_POLYGON);
	glNormal3f(norm[2][0] * zoom, norm[2][1] * zoom, norm[2][2] * zoom);
	glVertex3fv(sumi1[1]);
	glVertex3fv(sumi1[3]);
	glVertex3fv(sumi2[3]);
	glVertex3fv(sumi2[1]);
	glEnd();

	glBegin(GL_POLYGON);
	glNormal3f(norm[3][0] * zoom, norm[3][1] * zoom, norm[3][2] * zoom);
	glVertex3fv(sumi1[3]);
	glVertex3fv(sumi1[2]);
	glVertex3fv(sumi2[2]);
	glVertex3fv(sumi2[3]);
	glEnd();

	glBegin(GL_POLYGON);
	glNormal3f(norm[4][0] * zoom, norm[4][1] * zoom, norm[4][2] * zoom);
	glVertex3fv(sumi1[2]);
	glVertex3fv(sumi1[0]);
	glVertex3fv(sumi2[0]);
	glVertex3fv(sumi2[2]);
	glEnd();

	glBegin(GL_POLYGON);
	glNormal3f(norm[5][0] * zoom, norm[5][1] * zoom, norm[5][2] * zoom);
	glVertex3fv(sumi2[0]);
	glVertex3fv(sumi2[1]);
	glVertex3fv(sumi2[3]);
	glVertex3fv(sumi2[2]);
	glEnd();


}


void calc_quat_for_mouse_event(float mx, float my, float quat0[4]){

	float c = (float)cos(my);
	float s = (float)sin(my);

	float rw = quat0[0];
	float rx = quat0[1];
	float ry = quat0[2];
	float rz = quat0[3];

	float qw = c*rw - s*rx;
	float qx = c*rx + s*rw;
	float qy = c*ry - s*rz;
	float qz = c*rz + s*ry;

	c = (float)cos(mx);
	s = (float)sin(mx);

	rw = c*qw - s*qy;
	rx = c*qx + s*qz;
	ry = c*qy + s*qw;
	rz = c*qz - s*qx;


	float n = (float)sqrt(rw*rw + rx*rx + ry*ry + rz*rz);
	if (n != 0){
		rw /= n;
		rx /= n;
		ry /= n;
		rz /= n;
	}
	else{
		rw = 1.0f;
		rx = ry = rz = 0.0f;
	}

	quat0[0] = rw; quat0[1] = rx; quat0[2] = ry; quat0[3] = rz;
}


void jetColorMap(float val, float rgb[3]){

	if (val > 1.0) val = 1.0;
	if (val < -1.0)val = -1.0;

	float rgb1[3];
	float rgb2[3];
	float coef = 0.0;
	if (val < -0.75){
		//0-1
		coef = fabs(val + 0.75) / 0.25;
		rgb1[0] = 0.0 / 255.0; rgb1[1] = 0.0 / 255.0; rgb1[2] = 144.0 / 255.0;
		rgb2[0] = 0.0 / 255.0; rgb2[1] = 15.0 / 255.0; rgb2[2] = 255.0 / 255.0;
	}
	else if (val < -0.5){
		//1-2
		coef = fabs(val + 0.5) / 0.25;
		rgb1[0] = 0.0 / 255.0; rgb1[1] = 15.0 / 255.0; rgb1[2] = 255.0 / 255.0;
		rgb2[0] = 0.0 / 255.0; rgb2[1] = 144.0 / 255.0; rgb2[2] = 255.0 / 255.0;
	}
	else if (val < -0.25){
		//2-3
		coef = fabs(val + 0.25) / 0.25;
		rgb1[0] = 0.0 / 255.0; rgb1[1] = 144.0 / 255.0; rgb1[2] = 255.0 / 255.0;
		rgb2[0] = 15.0 / 255.0; rgb2[1] = 255.0 / 255.0; rgb2[2] = 238.0 / 255.0;
	}
	else if (val < 0.0){
		//3-4
		coef = fabs(val + 0.0) / 0.25;
		rgb1[0] = 15.0 / 255.0; rgb1[1] = 255.0 / 255.0; rgb1[2] = 238.0 / 255.0;
		rgb2[0] = 144.0 / 255.0; rgb2[1] = 255.0 / 255.0; rgb2[2] = 112.0 / 255.0;
	}
	else if (val < 0.25){
		//4-5
		coef = fabs(val - 0.25) / 0.25;
		rgb1[0] = 144.0 / 255.0; rgb1[1] = 255.0 / 255.0; rgb1[2] = 112.0 / 255.0;
		rgb2[0] = 255.0 / 255.0; rgb2[1] = 238.0 / 255.0; rgb2[2] = 0.0 / 255.0;
	}
	else if (val < 0.5){
		//5-6
		coef = fabs(val - 0.5) / 0.25;
		rgb1[0] = 255.0 / 255.0; rgb1[1] = 238.0 / 255.0; rgb1[2] = 0.0 / 255.0;
		rgb2[0] = 255.0 / 255.0; rgb2[1] = 112.0 / 255.0; rgb2[2] = 0.0 / 255.0;
	}
	else if (val < 0.75){
		//6-7
		coef = fabs(val - 0.75) / 0.25;
		rgb1[0] = 255.0 / 255.0; rgb1[1] = 112.0 / 255.0; rgb1[2] = 0.0 / 255.0;
		rgb2[0] = 238.0 / 255.0; rgb2[1] = 0.0 / 255.0; rgb2[2] = 0.0 / 255.0;
	}
	else {
		//7-8
		coef = fabs(val - 1.0) / 0.25;
		rgb1[0] = 238.0 / 255.0; rgb1[1] = 0.0 / 255.0; rgb1[2] = 0.0 / 255.0;
		rgb2[0] = 127.0 / 255.0; rgb2[1] = 0.0 / 255.0; rgb2[2] = 0.0 / 255.0;
	}

	for (int i = 0; i < 3; i++){
		rgb[i] = coef*rgb1[i] + (1.0 - coef)*rgb2[i];
	}

}



