#include <QtGui>
#include <QtWidgets>
#include <QtOpenGL>

#include "curve_design_display.h"
#include "EditTools.h"
#include "geometric_data_structure.h"
#include "curve_network.h"
#include "AABB_data_structure.h"
#include "basic_geometric_calculation.h"
#include "DisplayFunctions.h"
#include "MomentReductionMap.h"


//コンストラクタ
curve_design_display::curve_design_display(QWidget *parent) : QGLWidget(parent)
{
	setFormat(QGLFormat(QGL::DoubleBuffer | QGL::DepthBuffer));

	shiftX = 0.0;
	shiftY = 0.0;
	shiftZ = 0.0;
	Gpos[0] = 0.0;
	Gpos[1] = 0.0;
	Gpos[2] = 0.0;
	Zoom = 1.0;
	Zdefault = Zoom;
	windowRatio = 1.0;
	windowWid = 100;
	windowHei = 100;

	rw = 0.0;
	rx = 1.0;
	ry = 0.0;
	rz = 0.0;

	Cends = new struct curve_node[2];
	Curv = new class curve_on_mesh;
	Curv->nodeN = 0;
	moveP = new struct curve_node;

	connect(this, SIGNAL(sync_state_call(struct state_param)), parent, SLOT(sync_state_called(struct state_param)));

	myState.eState = PREPSTATE;
	Pressed_Key = Qt::Key_F35;
}


curve_design_display::~curve_design_display(){

	printf("Curve design dest start\n");

	delete[] Cends;
	//delete Curv;
	delete[] moveP;

	printf("Curve design dest end\n");
}


/////////
//設定・準備
/////////


void curve_design_display::initializeGL(){

	glClearColor(0.9, 0.9, 0.9, 1.0);

	GLfloat light_position[] = { 50.0, 25.0, 100.0, 1.0 };
	GLfloat light_ambient[] = { 0.25, 0.25, 0.25, 1.0 };
	GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };

	/* set up ambient, diffuse, and specular components for light 0 */
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);

	glEnable(GL_LIGHT0);   /* enable light 0 */
	glDepthFunc(GL_LEQUAL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_DEPTH_TEST);

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	reset_param_for_mesh();

}


void curve_design_display::resizeGL(int width, int height){
	
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	GLfloat x = GLfloat(width) / height;
	windowRatio = x;
	windowWid = width;
	windowHei = height;
	//glFrustum(-x, x, -1.0, 1.0, 4.0, 15.0);
	glOrtho(-x, x, -1.0, 1.0, -100, 100);
	glMatrixMode(GL_MODELVIEW);

}


void curve_design_display::paintGL(){
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	draw();
}


void curve_design_display::reset_param_for_mesh(){

	if (imesh->nodeN == 0) return;

	float limit[3][2];
	for (int i = 0; i < imesh->nodeN; i++){
		for (int j = 0; j < 3; j++){
			if (i == 0 || limit[j][0] > imesh->Nodes[i][j]) limit[j][0] = imesh->Nodes[i][j];
			if (i == 0 || limit[j][1] < imesh->Nodes[i][j]) limit[j][1] = imesh->Nodes[i][j];
		}
	}

	Gpos[0] = (limit[0][0] + limit[0][1]) / 2.0;
	Gpos[1] = (limit[1][0] + limit[1][1]) / 2.0;
	Gpos[2] = (limit[2][0] + limit[2][1]) / 2.0;

	shiftX = 0.0;
	shiftY = 0.0;
	shiftZ = 0.0;

	float maxV = 0.0;
	for (int i = 0; i < 3; i++){
		if (maxV < limit[i][1] - limit[i][0]){
			maxV = limit[i][1] - limit[i][0];
		}
	}
	Zoom = 1.0 / maxV*2.0;
	Zdefault = Zoom;
}


QSize curve_design_display::sizeHint() const{

	return QSize(1000, 1000);

}



////////////
//イベント関連
////////////

void curve_design_display::mousePressEvent(QMouseEvent *event){
	
	if (myState.eState == ADDCSTATE){
		if ((event->buttons() & Qt::LeftButton) && Pressed_Key == Qt::Key_F35) ADDCSEARCH = true;
		else ADDCSEARCH = false;
	}

	if (myState.eState == EDITCSTATE){
		selectC = -1;
		selectN = -1;
		selectE = -1;
		if ((event->buttons() & Qt::LeftButton) && Pressed_Key == Qt::Key_F35){
			search_selected_element(event->x(), event->y());
			if (selectE != -1){
				class curve_on_mesh *scur = cNet->curSet[selectC];
				float pos[3];
				scur->dCons[selectE].cType = 1;
				for (int i = 0; i < 3; i++){
					scur->dCons[selectE].pos[i] = scur->cNode[selectE + 1].pos[i] - scur->cNode[selectE].pos[i];
					pos[i] = (scur->cNode[selectE + 1].pos[i] + scur->cNode[selectE].pos[i]) / 2.0;
				}
				int fn;
				float cord[3], duv[2][3];
				iAABB->closest_point_search_with_mesh_info(pos, &fn, cord);
				imesh->evaluate_derivative_of_box_spline_on_triangle(fn, cord, duv);
				float ndh = 0.0;
				for (int i = 0; i < 3; i++){
					eNorm[i] = duv[0][(i + 1) % 3] * duv[1][(i + 2) % 3] - duv[0][(i + 2) % 3] * duv[1][(i + 1) % 3];
					ndh += eNorm[i] * eNorm[i];
				}
				ndh = sqrt(ndh);
				if (ndh < 1.0e-10) ndh = 1.0e-10;
 				for (int i = 0; i < 3; i++){
					eNorm[i] /= ndh;
				}
			}
		}
	}

	if (myState.eState == DELCSTATE){
		selectC = -1;
		selectN = -1;
		selectE = -1;
		if ((event->buttons() & Qt::LeftButton) && Pressed_Key == Qt::Key_F35){
			search_selected_element(event->x(), event->y());
		}
	}

	if (myState.eState == DEFCSTATE){
		if ((event->buttons() & Qt::LeftButton) && Pressed_Key == Qt::Key_F35){
			CONNECTSEARCH = true;
			int s1 = selem1[0];
			search_selected_element(event->x(), event->y());
			if (s1 == -1 && selem1[0] == -1) CONNECTSEARCH = false;
			if (s1 != -1 && selem1[1] == -1) CONNECTSEARCH = false;
			if (sconn != -1) CONNECTSEARCH = true;
		}
	}

	lastPos = event->pos();

}


void curve_design_display::mouseMoveEvent(QMouseEvent *event){

	//曲線の端点設定用
	if (myState.eState == ADDCSTATE && (myState.acState != PLANESEARCH && myState.acState != PLANESETTED)){
		if ((event->buttons() & Qt::LeftButton) && Pressed_Key == Qt::Key_F35){
			float leng = sqrt(pow(event->x() - lastPos.x(), 2) + pow(event->x() - lastPos.y(), 2));
			if (leng > 5.0){
				ADDCSEARCH = false;
			} 
			return;
		}
	}
	//曲線の切断面設定用
	else if (myState.eState == ADDCSTATE){
		if ((event->buttons() & Qt::LeftButton) && Pressed_Key == Qt::Key_F35){
			update_cutting_plane(event->x() - lastPos.x(), event->y() - lastPos.y());
			lastPos = event->pos();
			updateGL();
			return;
		}
		else {
			ADDCSEARCH = false;
		}
	}
	//曲線のエディット
	if (myState.eState == EDITCSTATE){
		if ((event->buttons() & Qt::LeftButton) && Pressed_Key == Qt::Key_F35){
			if (myState.ecState == POINTMOVE && selectN != -1){
				struct curve_node *cnode;
				cnode = new struct curve_node;
				cnode->fn = -1;
				search_clicked_point(event->x(), event->y(), cnode);
				if (cnode->fn != -1){
					float leng = 0.0;
					for (int i = 0; i < 3; i++){
						leng += pow(cNet->curSet[selectC]->cNode[selectN].pos[i] - cnode->pos[i], 2);
					}
					leng = sqrt(leng);
					if (leng > cNet->sParam->Wid*5.0){
						//selectC = -1;
						//selectN = -1;
						//moveP->fn = -1;
					}
					else {
						moveP->fn = cnode->fn;
						for (int i = 0; i < 3; i++){
							moveP->pos[i] = cnode->pos[i];
							cNet->curSet[selectC]->pCons[selectN].pos[i] = cnode->pos[i];
							if (i != 2) moveP->uv[i] = cnode->uv[i];
						}
						cNet->curSet[selectC]->pCons[selectN].cType = 1;
						bool is_converged = cNet->curSet[selectC]->minimize_discrete_rod_energy(100);
						cNet->curSet[selectC]->remesh_curve_node(selectC, selectN, selectE, iAABB, cNet->Connect);
						if (!is_converged){
							moveP->fn = cNet->curSet[selectC]->cNode[selectN].fn;
							for (int i = 0; i < 3; i++){
								moveP->pos[i] = cNet->curSet[selectC]->cNode[selectN].pos[i];
								cNet->curSet[selectC]->pCons[selectN].pos[i] = cNet->curSet[selectC]->cNode[selectN].pos[i];
								if (i != 2) moveP->uv[i] = cNet->curSet[selectC]->cNode[selectN].uv[i];
							}
						}
					}
				}
				delete cnode;
				lastPos = event->pos();
				updateGL();
				return;
			}
			if (myState.ecState == EDGEDIRE && selectE != -1){
				update_edge_direction(event->x(),event->y(),event->x()-lastPos.x(), event->y()-lastPos.y());
				lastPos = event->pos();
				updateGL();
				return;
			}
			if (myState.ecState == DELETECONS && selectC != -1){
				float leng = sqrt(pow(lastPos.x() - event->x(), 2) + pow(lastPos.y() - event->y(), 2)) / Zoom;
				if (leng > cNet->sParam->Wid*3.0){
					selectC = -1;
					selectN = -1;
					selectE = -1;
				}
				updateGL();
				return;
			}
		}
	}

	//曲線の消去
	if (myState.eState == DELCSTATE && selectC != -1 && (event->buttons() & Qt::LeftButton) && Pressed_Key == Qt::Key_F35){
		float leng = sqrt(pow(lastPos.x() - event->x(), 2) + pow(lastPos.y() - event->y(), 2)) / Zoom;
		if (leng > cNet->sParam->Wid*3.0){
			selectC = -1;
			selectN = -1;
			selectE = -1;
		}
		updateGL();
		return;
	}
	
	//コネクションの編集
	if (myState.eState == DEFCSTATE && CONNECTSEARCH && (event->buttons() & Qt::LeftButton) && Pressed_Key == Qt::Key_F35){
		float leng = sqrt(pow(lastPos.x() - event->x(), 2) + pow(lastPos.y() - event->y(), 2)) / Zoom;
		if (leng > cNet->sParam->Wid*3.0){
			if (selem2[0] != -1) selem2[0] = -1;
			else selem1[0] = -1;
			CONNECTSEARCH = true;
			sconn = -1;
		}
		updateGL();
		return;
	}

	//通常の操作
	if (event->buttons() & Qt::RightButton){
		float mx = -0.0025f*(event->x() - lastPos.x());
		float my = 0.0025f*(event->y() - lastPos.y());

		float quat[4];
		quat[0] = rw; quat[1] = rx; quat[2] = ry; quat[3] = rz;
		calc_quat_for_mouse_event(mx, my, quat);

		rw = quat[0]; rx = quat[1]; ry = quat[2]; rz = quat[3];
		updateGL();
	}
	else if (event->buttons() && Pressed_Key == Qt::Key_Control){
		Zoom -= 0.0025f*(event->y() - lastPos.y())*Zdefault;
		if (Zoom > 20.0f*Zdefault) Zoom = 20.0f*Zdefault;
		else if (Zoom < 0.05f*Zdefault) Zoom = 0.05f*Zdefault;
		resizeGL(width(), height());
		updateGL();
	}
	else if (event->buttons() && Pressed_Key == Qt::Key_Shift){
		shiftX += 0.0025f*(event->x() - lastPos.x()) / Zoom;
		shiftY -= 0.0025f*(event->y() - lastPos.y()) / Zoom;
		updateGL();
	}
	
	lastPos = event->pos();

}


void curve_design_display::mouseReleaseEvent(QMouseEvent *event){
	

	if (myState.eState == ADDCSTATE && ADDCSEARCH && myState.acState != PLANESEARCH && myState.acState != PLANESETTED){
		struct curve_node *cnode;
		cnode = new struct curve_node;
		cnode->fn = -1;
		search_clicked_point(event->x(), event->y(), cnode);
		bool acstate_changed = false;
		if (cnode->fn != -1){
			int epcount = 0;
			int minN = -1;
			float minL;
			for (int i = 0; i < 2; i++){
				if (Cends[i].fn != -1){
					float leng = 0.0;
					for (int j = 0; j < 3; j++){
						leng += pow(cnode->pos[j] - Cends[i].pos[j], 2);
					}
					leng = sqrt(leng);
					if (minN == -1 || leng < minL){
						minN = i;
						minL = leng;
					}
					epcount++;
				}
			}
			if (minN != -1 && minL < cNet->sParam->Wid){
				Cends[minN].fn = -1;
				epcount--;
			}
			else if (epcount == 0){
				Cends[0].fn = cnode->fn;
				for (int j = 0; j < 3; j++){
					Cends[0].pos[j] = cnode->pos[j];
					if (j != 2) Cends[0].uv[j] = cnode->uv[j];
				}
				epcount++;
			}
			else{
				if (epcount == 1) minN = (minN + 1) % 2;
				Cends[minN].fn = cnode->fn;
				for (int j = 0; j < 3; j++){
					Cends[minN].pos[j] = cnode->pos[j];
					if (j != 2) Cends[minN].uv[j] = cnode->uv[j];
				}
				epcount++;
			}
			if (epcount == 0 && myState.acState != POINTZERO){
				myState.acState = POINTZERO;
				acstate_changed = true;
			}
			if (epcount == 1 && myState.acState != POINTONE){
				myState.acState = POINTONE;
				acstate_changed = true;
			}
			if (epcount == 2 && myState.acState != POINTTWO){
				myState.acState = POINTTWO;
				acstate_changed = true;
			}
		}
		if (acstate_changed) emit sync_state_call(myState);
		delete cnode;
	}

	if (myState.eState == ADDCSTATE && ADDCSEARCH && (myState.acState == PLANESEARCH || myState.acState == PLANESETTED)){
		if (Curv->nodeN > 0) delete[] Curv->cNode;
		Curv->nodeN = 0;
		bool CCEXIST = cross_section_curve_search(Curv);
		bool STATEUP = false;
		if (CCEXIST){
			if (myState.acState == PLANESEARCH) STATEUP = true;
			myState.acState = PLANESETTED;
		}
		else {
			if (myState.acState == PLANESETTED) STATEUP = true;
			myState.acState = PLANESEARCH;
		}
		if(STATEUP) emit sync_state_call(myState);
	}

	if (myState.eState == EDITCSTATE){
		if (myState.ecState == EDGEDIRE && selectE != -1){
			bool is_converged = cNet->curSet[selectC]->minimize_discrete_rod_energy(200);
			cNet->curSet[selectC]->remesh_curve_node(selectC, selectN, selectE, iAABB, cNet->Connect);
			if (!is_converged){
				float ed[3];
				float edh = 0.0;
				for (int i = 0; i < 3; i++){
					ed[i] = cNet->curSet[selectC]->cNode[selectE + 1].pos[i] - cNet->curSet[selectC]->cNode[selectE].pos[i];
					edh += ed[i] * ed[i];
				}
				edh = sqrt(edh);
				if (edh < 1.0e-10) edh = 1.0e-10;
				for (int i = 0; i < 3; i++){
					cNet->curSet[selectC]->dCons[selectE].pos[i] = ed[i] / edh;
				}
			}
		}
		if (myState.ecState == DELETECONS && selectC != -1){
			class curve_on_mesh *scur = cNet->curSet[selectC];
			if (selectN != -1){
				scur->pCons[selectN].cType = -1;
			}
			if (selectE != -1){
				scur->dCons[selectE].cType = -1;
			}
			cNet->curSet[selectC]->minimize_discrete_rod_energy(200);
			cNet->curSet[selectC]->remesh_curve_node(selectC,selectN, selectE, iAABB,cNet->Connect);
		}

		selectC = -1;
		selectN = -1;
		selectE = -1;
	}

	if (myState.eState == DELCSTATE && selectC != -1){
		cNet->deleteCurve(selectC);
		selectC = -1;
		selectN = -1;
		selectE = -1;
	}

	if (myState.eState == DEFCSTATE && CONNECTSEARCH){
		if (selem1[0] != -1 && selem1[0] == selem2[0] && selem1[1] == selem2[1]){
			selem1[0] = -1;
			selem1[1] = -1;
			selem2[0] = -1;
			selem2[1] = -1;
		}
		if (sconn != -1){
			cNet->deleteConnect(sconn);
			sconn = -1;
			selem1[0] = -1;
			selem1[1] = -1;
			selem2[0] = -1;
			selem2[1] = -1;
		}
		else if (selem2[0] != -1){
			cNet->addConnect(selem1[0], selem1[1], selem2[0], selem2[1]);
			selem1[0] = -1;
			selem1[1] = -1;
			selem2[0] = -1;
			selem2[1] = -1;
		}
	}



	updateGL();
}


void curve_design_display::addCurveToNet(){

	int cN = cNet->curSet.size();
	cNet->addCurve(Curv);

	cNet->curSet[cN]->prepare_for_minimization();
	
	std::vector<struct ConstraintData> pcSet, dcSet;
	struct ConstraintData pcons1,pcons2, dcons1,dcons2;
	pcons1.pn = 0;
	pcons2.pn = cNet->curSet[cN]->nodeN - 1;
	for (int i = 0; i < 3; i++){
		pcons1.pos[i] = cNet->curSet[cN]->cNode[0].pos[i];
		pcons2.pos[i] = cNet->curSet[cN]->cNode[cNet->curSet[cN]->nodeN - 1].pos[i];
	}
	pcSet.push_back(pcons1);
	pcSet.push_back(pcons2);
	dcons1.pn = 0;
	dcons2.pn = cNet->curSet[cN]->nodeN - 2;
	for (int i = 0; i < 3; i++){
		dcons1.pos[i] = cNet->curSet[cN]->cNode[1].pos[i] - cNet->curSet[cN]->cNode[0].pos[i];
		dcons2.pos[i] = cNet->curSet[cN]->cNode[dcons2.pn + 1].pos[i] - cNet->curSet[cN]->cNode[dcons2.pn].pos[i];
	}
	dcSet.push_back(dcons1);
	dcSet.push_back(dcons2);
	cNet->curSet[cN]->setPointConstraint(pcSet);
	cNet->curSet[cN]->setDirectionConstraint(dcSet);

	cNet->curSet[cN]->minimize_discrete_rod_energy(200);
	int gomi;
	cNet->curSet[cN]->remesh_curve_node(gomi,gomi,gomi,iAABB,cNet->Connect);

	if (cNet->curSet[cN]->CONNECT_LIST_EXIST){
		for (int i = 0; i < cNet->curSet[cN]->nodeN; i++){
			if (cNet->curSet[cN]->Csize[i] > 0) delete[] cNet->curSet[cN]->Clist[i];
		}
		delete[] cNet->curSet[cN]->Csize;
		delete[] cNet->curSet[cN]->Clist;
	}
	cNet->curSet[cN]->CONNECT_LIST_EXIST;

	cNet->curSet[cN]->Csize = new int[cNet->curSet[cN]->nodeN];
	cNet->curSet[cN]->Clist = new int*[cNet->curSet[cN]->nodeN];
	for (int i = 0; i < cNet->curSet[cN]->nodeN; i++){
		cNet->curSet[cN]->Csize[i] = 0;
	}


}


////////////
//表示関連
////////////


void curve_design_display::draw(){

	float Rot[3][3], quat[4];
	quat[0] = rw; quat[1] = rx; quat[2] = ry; quat[3] = rz;
	derive_rot_mat_from_quaternion(quat, Rot);
	float Trans[3];
	for (int i = 0; i < 3; i++){
		Trans[i] = 0.0;
		for (int j = 0; j < 3; j++){
			Trans[i] -= Rot[i][j] * Gpos[j];
		}
	}
	Trans[0] = (Trans[0] + shiftX)*Zoom;
	Trans[1] = (Trans[1] + shiftY)*Zoom;
	Trans[2] = (Trans[1] + shiftZ)*Zoom;

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(Trans[0], Trans[1], Trans[2]);
	//glTranslatef(shiftX*Zoom, shiftY*Zoom, shiftZ*Zoom);
	glScalef(Zoom, Zoom, Zoom);
	float rh = sqrt(rx*rx + ry*ry + rz*rz);
	if (rh != 0.0) glRotatef(2.0*acos(rw)*180.0/M_PI, rx / rh, ry / rh, rz / rh);

	//printf("wid%f hei%f\n",cNet->hei,cNet->wid);
	display_mesh();

	if (myState.eState == ADDCSTATE) display_curve_ends();
	display_curve_network();
	display_connection();

	display_moment_reduction_map();


	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	displayEditState();
}


void curve_design_display::display_mesh(){


	//glPolygonMode(GL_FRONT, GL_FILL);
	glEnable(GL_LIGHTING);
	for (int i = 0; i < imesh->faceN; i++){
		float fps[3][3];
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				fps[j][k] = imesh->Nodes[imesh->Faces[i][j]][k];
			}
		}

		float norm[3];
		float nh = 0.0;
		for (int j = 0; j < 3; j++){
			norm[j] = (fps[1][(j + 1) % 3] - fps[0][(j + 1) % 3]) *(fps[2][(j + 2) % 3] - fps[0][(j + 2) % 3])
				- (fps[1][(j + 2) % 3] - fps[0][(j + 2) % 3]) *(fps[2][(j + 1) % 3] - fps[0][(j + 1) % 3]);
			nh += norm[j] * norm[j];
		}
		nh = sqrt(nh);
		for (int j = 0; j < 3; j++){
			norm[j] /= nh/Zoom;
		}
		float cval[3];
		cval[0] = 0.98; cval[1] = 0.98; cval[2] = 0.98;
		cval[0] = 180.0 / 256.0; cval[1] = 200 / 256.0; cval[2] = 220.0 / 256.0;
		choose_color_data(8);
		set_plastic_color(cval);
		glBegin(GL_POLYGON);
		glNormal3fv(norm);
		glVertex3fv(fps[0]);
		glVertex3fv(fps[1]);
		glVertex3fv(fps[2]);
		glEnd();
	}

}


void curve_design_display::display_curve_ends(){


	glEnable(GL_LIGHTING);
	for (int i = 0; i < 2; i++){
		if (Cends[i].fn != -1){
			float cval[3];
			cval[0] = 0.9; cval[1] = 0.0; cval[2] = 0.2;
			set_plastic_color(cval);
			
			disp_sphere(Cends[i].pos, cNet->sParam->Wid / 2.0, Zoom);
		}
	}

	if (myState.acState == PLANESEARCH || myState.acState == PLANESETTED){
		glDisable(GL_LIGHTING);
		float pd1[3], pd2[3];
		float pdh1 = 0.0;
		float pdh2 = 0.0;
		for (int i = 0; i < 3; i++){
			pd1[i] = Cends[1].pos[i] - Cends[0].pos[i];
			pdh1 += pd1[i] * pd1[i];
		}
		for (int i = 0; i < 3; i++){
			pd2[i] = pd1[(i + 1) % 3] * cutNorm[(i + 2) % 3] - pd1[(i + 2) % 3] * cutNorm[(i + 1) % 3];
			pdh2 += pd2[i] * pd2[i];
		}
		pdh1 = sqrt(pdh1);
		pdh2 = sqrt(pdh2);

		float sumi[4][3];
		float W = pdh1 * 0.65;
		float L = pdh1;
		for (int i = 0; i < 3; i++){
			sumi[0][i] = Cends[0].pos[i];
			sumi[1][i] = Cends[0].pos[i] + W*pd2[i] /pdh2;
			sumi[2][i] = Cends[1].pos[i] + W*pd2[i] / pdh2;
			sumi[3][i] = Cends[1].pos[i];
		}
		glColor3f(0.6, 0.6, 0.6);
		glColor3f(250.0 / 256.0, 200.0 / 256.0, 20.0 / 256.0);
		//glColor3f(230.0 / 256.0,215.0 / 256.0, 225.0 / 256.0);
		glColor3f(170.0 / 256.0, 240.0 / 256.0, 150.0 / 256.0);
		glBegin(GL_POLYGON);
		glVertex3fv(sumi[0]);
		glVertex3fv(sumi[1]);
		glVertex3fv(sumi[2]);
		glVertex3fv(sumi[3]);
		glEnd();

		glEnable(GL_LIGHTING);
	}


	if (myState.acState == PLANESETTED){
		for (int i = 0; i < Curv->nodeN; i++){
			float cval[3];
			cval[0] = 0.9; cval[1] = 0.0; cval[2] = 0.2;
			set_plastic_color(cval);
			disp_sphere(Curv->cNode[i].pos, cNet->sParam->Wid/2.0, Zoom);
		}
		//for (int i = 0; i < cnNum; i++){
		//	/*glLineWidth(5);
		//	glColor3f(1, 0, 0);
		//	glBegin(GL_LINES);
		//	glVertex3fv(cnSet[i].pos);
		//	glVertex3fv(cnSet[i + 1].pos);
		//	glEnd();*/

		//	float cval[3];
		//	cval[0] = 0.9; cval[1] = 0.0; cval[2] = 0.2;
		//	set_plastic_color(cval);
		//	/*disp_shpere(CNset[i].pos, rod_width, Zoom);*/
		//}
	}
}


void curve_design_display::display_curve_network(){


	glEnable(GL_LIGHTING);
	for (int i = 0; i < cNet->curSet.size(); i++){
		class curve_on_mesh *celem = cNet->curSet[i];
		choose_color_data(3);
		disp_curve_polygon(celem, Zoom, celem->sParam->Hei*0.5);

		float(*sNorm)[3];
		sNorm = new float[celem->nodeN][3];
		for (int j = 0; j < celem->nodeN; j++){
			float uvw[3],duv[2][3];
			uvw[1] = celem->cNode[j].uv[0];
			uvw[2] = celem->cNode[j].uv[1];
			uvw[0] = 1.0 - uvw[1] - uvw[2];
			imesh->evaluate_derivative_of_box_spline_on_triangle(celem->cNode[j].fn, uvw, duv);
			float nh = 0.0;
			for (int k = 0; k < 3; k++){
				sNorm[j][k] = duv[0][(k + 1) % 3] * duv[1][(k + 2) % 3] - duv[0][(k + 2) % 3] * duv[1][(k + 1) % 3];
				nh += sNorm[j][k] * sNorm[j][k];
			}
			nh = sqrt(nh);
			for (int k = 0; k < 3; k++){
				sNorm[j][k] /= nh;
			}
		}
		//頂点の表示
		for (int j = 0; j < celem->nodeN; j++){
			float cval[3];
			cval[0] = 0.0; cval[1] = 0.9; cval[2] = 0.2;
			if (celem->pCons[j].cType != -1){
				cval[0] = 0.2; cval[1] = 0.0; cval[2] = 0.9;
			}
			if (myState.eState == EDITCSTATE && selectC == i && selectN == j){
				cval[0] = 0.9; cval[1] = 0.0; cval[2] = 0.2;
			}
			if (myState.eState == DEFCSTATE){
				if ((selem1[0] == i && selem1[1] == j) || (selem2[0] == i && selem2[1] == j)){
					cval[0] = 0.9; cval[1] = 0.0; cval[2] = 0.2;
				}
			}

			float spos[3];
			float off = 0.5*celem->sParam->Hei;
			for (int k = 0; k < 3; k++){
				spos[k] = celem->cNode[j].pos[k] + sNorm[j][k] * off;
			}

			set_plastic_color(cval);
			disp_sphere(spos, cNet->sParam->Wid/ 2.0, Zoom);
		}
		//拘束エッジの表示
		for (int j = 0; j < celem->nodeN - 1; j++){
			if (celem->dCons[j].cType != -1 || (myState.eState == EDITCSTATE && selectC == i && selectE == j)){
				float cent[3], ed[3];
				float edh = 0.0;
				for (int k = 0; k < 3; k++){
					cent[k] = (celem->cNode[j].pos[k] + celem->cNode[j + 1].pos[k]) / 2.0;
					ed[k] = celem->cNode[j + 1].pos[k] - celem->cNode[j].pos[k];
					ed[k] = celem->dCons[j].pos[k];
					edh += ed[k] * ed[k];
				}
				edh = sqrt(edh);
				float pos1[3],pos2[3];
				for (int k = 0; k < 3; k++){
					pos1[k] = cent[k] + ed[k] / edh*celem->sParam->Wid*1.0 + (sNorm[j][k] + sNorm[j + 1][k]) / 2.0*celem->sParam->Wid*1.0;
					pos2[k] = cent[k] - ed[k] / edh*celem->sParam->Wid*1.0 + (sNorm[j][k] + sNorm[j + 1][k]) / 2.0*celem->sParam->Wid*1.0;
				}
				float cval[3];
				cval[0] = 0.2; cval[1] = 0.0; cval[2] = 0.9;
				if (myState.eState == EDITCSTATE && selectC == i && selectE == j){
					cval[0] = 0.9; cval[1] = 0.0; cval[2] = 0.2;
				}
				set_plastic_color(cval);
				disp_arrow(pos1, pos2, celem->sParam->Wid*0.3, Zoom);
			}
		}
		delete[] sNorm;

	}

	if (myState.eState == EDITCSTATE && selectN != -1 && selectC != -1){
		if (moveP->fn != -1){
			float spos[3], snorm[3];
			float duv[2][3], uvw[3];
			uvw[0] = 1.0 - moveP->uv[0] - moveP->uv[1];
			imesh->evaluate_derivative_of_box_spline_on_triangle(moveP->fn, uvw, duv);
			float nh = 0.0;
			for (int i = 0; i < 3; i++){
				snorm[i] = duv[0][(i + 1) % 3] * duv[1][(i + 2) % 3] - duv[0][(i + 2) % 3] * duv[1][(i + 1) % 3];
				nh += snorm[i] * snorm[i];
			}
			nh = sqrt(nh);
			float off = 0.5*cNet->sParam->Hei;
			for (int i = 0; i < 3; i++){
				spos[i] = moveP->pos[i] + snorm[i] / nh;
			}

			float cval[3];
			cval[0] = 0.9; cval[1] = 0.0; cval[2] = 0.2;
			set_plastic_color(cval);
			disp_sphere(spos, cNet->sParam->Wid / 2.0, Zoom);
		}
	}

}


void curve_design_display::display_connection(){


	for (int i = 0; i < cNet->Connect.size(); i++){
		struct curve_on_mesh *cur1 = cNet->curSet[cNet->Connect[i].c1];
		struct curve_on_mesh *cur2 = cNet->curSet[cNet->Connect[i].c2];
		int n1 = cNet->Connect[i].n1;
		int n2 = cNet->Connect[i].n2;

		float pos1[3],pos2[3];
		for (int j = 0; j < 3; j++){
			pos1[j] = cur1->cNode[n1].pos[j];
			pos2[j] = cur2->cNode[n2].pos[j];
		}
		float uvw1[3], uvw2[3],duv1[2][3],duv2[2][3];
		uvw1[1] = cur1->cNode[n1].uv[0];
		uvw1[2] = cur1->cNode[n1].uv[1];
		uvw1[0] = 1.0 - uvw1[1] - uvw1[2];
		uvw2[1] = cur2->cNode[n2].uv[0];
		uvw2[2] = cur2->cNode[n2].uv[1];
		uvw2[0] = 1.0 - uvw2[1] - uvw2[2];
		imesh->evaluate_derivative_of_box_spline_on_triangle(cur1->cNode[n1].fn, uvw1, duv1);
		imesh->evaluate_derivative_of_box_spline_on_triangle(cur2->cNode[n2].fn, uvw2, duv2);

		float norm1[3], norm2[3];
		float nh1 = 0.0;
		float nh2 = 0.0;
		for (int j = 0; j < 3; j++){
			norm1[j] = duv1[0][(j + 1) % 3] * duv1[1][(j + 2) % 3] - duv1[0][(j + 2) % 3] * duv1[1][(j + 1) % 3];
			norm2[j] = duv2[0][(j + 1) % 3] * duv2[1][(j + 2) % 3] - duv2[0][(j + 2) % 3] * duv2[1][(j + 1) % 3];
			nh1 += norm1[j] * norm1[j];
			nh2 += norm2[j] * norm2[j];
		}
		nh1 = sqrt(nh1);
		nh2 = sqrt(nh2);
		if (nh1 < 1.0e-7) nh1 = 1.0e-7;
		if (nh2 < 1.0e-7) nh2 = 1.0e-7;
		for (int j = 0; j < 3; j++){
			norm1[j] /= nh1;
			norm2[j] /= nh2;
		}

		float cval[3];
		cval[0] = 0.9; cval[1] = 0.0; cval[2] = 0.2;
		if (myState.eState == DEFCSTATE && sconn == i){
			cval[0] = 0.2; cval[1] = 0.0; cval[2] = 0.9;
		}
		set_plastic_color(cval);
		disp_connection_polygon(pos1, pos2, norm1, norm2, cNet->sParam->Wid*1.1, cNet->sParam->Hei*1.05, Zoom);

	}

}


void curve_design_display::display_moment_reduction_map(){

	if (mrMap->wrN == 0)return;

	glDisable(GL_LIGHTING);
	for (int i = 0; i < mrMap->wrN; i++){
		float pos[3];
		for (int j = 0; j < 3; j++){
			pos[j] = mrMap->wReg[i].pos[j];
		}
		float cval[3];
		cval[0] = 1.0; cval[1] = 0.0; cval[2] = 0.0;
		/*disp_sphere(pos, cNet->sParam->Wid / 1.0, Zoom);*/

		glDisable(GL_LIGHTING);
		Mesh *mMesh = mrMap->wReg[i].mMesh;
		for (int j = 0; j < mMesh->faceN; j++){
			float fps[3][3];
			//float cval = 0.0;
			for (int k = 0; k < 3; k++){
				for (int l = 0; l < 3; l++){
					fps[k][l] = mMesh->Nodes[mMesh->Faces[j][k]][l];
				}
				//cval += mrMap->wReg[i].colmap[mMesh->Faces[j][k]] / 3.0;
			}

			/*if (cval > 0.0)glColor3f(cval, 1.0 - cval, 0);
			else glColor3f(0.0, 1.0 + cval, -cval);*/

			float rgb[3];
			glBegin(GL_TRIANGLES);
			jetColorMap(mrMap->wReg[i].colmap[mMesh->Faces[j][0]], rgb);
			glColor3fv(rgb);
			glVertex3fv(fps[0]);
			jetColorMap(mrMap->wReg[i].colmap[mMesh->Faces[j][1]], rgb);
			glColor3fv(rgb);
			glVertex3fv(fps[1]);
			jetColorMap(mrMap->wReg[i].colmap[mMesh->Faces[j][2]], rgb);
			glColor3fv(rgb);
			glVertex3fv(fps[2]);
			glEnd();


		}

		glEnable(GL_LIGHTING);
		float frp[3][3];
		float leng = cNet->sParam->Wid*5.0;
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				frp[j][k] = mrMap->wReg[i].pos[k] + leng*mrMap->wReg[i].Frame[j][k];
			}
		}
		cval[0] = 1.0; cval[1] = 0.0; cval[2] = 0.0;
		set_plastic_color(cval);
		disp_arrow(pos, frp[0], cNet->sParam->Hei / 2.0, Zoom);
		cval[0] = 0.0; cval[1] = 0.0; cval[2] = 1.0;
		set_plastic_color(cval);
		disp_arrow(pos, frp[1], cNet->sParam->Hei / 2.0, Zoom);
		
	}
	glDisable(GL_LIGHTING);

}


void curve_design_display::displayEditState(){

	QFont myFont("Helvetica [Cronyx]", 15);
	QFont myFont2("Helvetica [Cronyx]", 12);

	float psize = myFont.pointSizeF();
	float offH = 15.0*0.75 / (float)windowHei*5.5;
	//printf("%d %f %f \n", windowHei, offH, psize);
	glColor3f(1, 0, 0);
	if (myState.eState == ADDCSTATE){
		QString anno("Add Curves");
		renderText(-0.95*windowRatio, 0.90, 0.0, anno, myFont);
		glColor3f(0, 0, 0);
		if (myState.acState == POINTZERO || myState.acState == POINTONE || myState.acState == POINTTWO){
			QString subt(" Set End Points");
			renderText(-0.95*windowRatio, 0.90-offH, 0.0, subt, myFont2);
		}
		else if (myState.acState == PLANESEARCH || myState.acState == PLANESETTED){
			QString subt(" Set Cutting Plane");
			renderText(-0.95*windowRatio, 0.90-offH, 0.0, subt, myFont2);
		}
	}
	else if (myState.eState == EDITCSTATE){
		QString anno("Edit Curves");
		renderText(-0.95*windowRatio, 0.90, 0.0, anno, myFont);
		glColor3f(0, 0, 0);
		if (myState.ecState == POINTMOVE){
			QString subt(" Node Position");
			renderText(-0.95*windowRatio, 0.90 - offH, 0.0, subt, myFont2);
		}
		else if (myState.ecState == EDGEDIRE){
			QString subt(" Edge Direction");
			renderText(-0.95*windowRatio, 0.90 - offH, 0.0, subt, myFont2);
		}
		else if (myState.ecState == DELETECONS){
			QString subt(" Delete Constraints");
			renderText(-0.95*windowRatio, 0.90 - offH, 0.0, subt, myFont2);
		}
	}
	else if (myState.eState == DEFCSTATE){
		QString anno("Edit Connection");
		renderText(-0.95*windowRatio, 0.90, 0.0, anno, myFont);
	}
	else if (myState.eState == DELCSTATE){
		QString anno("Delete Curves");
		renderText(-0.95*windowRatio, 0.90, 0.0, anno, myFont);
	}
} 


//////////
///その他
/////////

void curve_design_display::state_changed(struct state_param state){

	if (myState.eState != ADDCSTATE && state.eState == ADDCSTATE){
		Cends[0].fn = -1;
		Cends[1].fn = -1;
	}
	if (myState.eState == ADDCSTATE && (myState.acState != PLANESEARCH  && myState.acState != PLANESETTED) && state.acState == PLANESEARCH){
		float np[3], quat[4], Rot[3][3];
		np[0] = 0.0; np[1] = sin(M_PI*0.15); np[2] = cos(M_PI*0.15);
		quat[0] = rw; quat[1] = rx; quat[2] = ry; quat[3] = rz;
		derive_rot_mat_from_quaternion(quat, Rot);
		
		float ed[3],pd[3];
		float seki = 0.0;
		for (int i = 0; i < 3; i++){
			pd[i] = 0.0;
			for (int j = 0; j < 3; j++){
				pd[i] += Rot[j][i] * np[j];
			}
			ed[i] = Cends[1].pos[i] - Cends[0].pos[i];
		}
		for (int i = 0; i < 3; i++){
			cutNorm[i] = pd[(i + 1) % 3] * ed[(i + 2) % 3] - pd[(i + 2) % 3] * ed[(i + 1) % 3];
			seki += cutNorm[i] * cutNorm[i];
		}
		seki = sqrt(seki);
		if (seki != 0){
			for (int i = 0; i < 3; i++){
				cutNorm[i] /= seki;
			}
		}
		else {
			cutNorm[0] = 1.0; cutNorm[1] = 0.0; cutNorm[2] = 0.0;
		}
	}
	if (myState.eState != EDITCSTATE && state.eState == EDITCSTATE){
		selectC = -1;
		selectN = -1;
		selectE = -1;
		moveP->fn = -1;
	}

	if (myState.eState != DELCSTATE && state.eState == DELCSTATE){
		selectC = -1;
		selectN = -1;
		selectE = -1;
	}

	if (myState.eState != DEFCSTATE && state.eState == DEFCSTATE){
		if (cNet->Connect.size() == 0) cNet->setInitConnection();
		selem1[0] = -1;
		selem2[0] = -1;
		sconn = -1;
	}

	myState = state;

	updateGL();
}


void curve_design_display::search_clicked_point(int cx,int cy,struct curve_node *cnode){

	float ppos[3],pdire[3];
	ppos[0] = 2.0*cx / (float)height() - (float)width()/(float)height() - shiftX*Zoom;
	ppos[1] = 1.0 - 2.0*cy / (float)height() - shiftY*Zoom;
	ppos[2] = -100.0 - shiftZ*Zoom;
	pdire[0] = 0.0;
	pdire[1] = 0.0;
	pdire[2] = -1.0;

	float Rot[3][3];
	float quat[4];
	quat[0] = rw;
	quat[1] = rx;
	quat[2] = ry;
	quat[3] = rz;
	derive_rot_mat_from_quaternion(quat,Rot);

	float borg[3], bdire[3];
	for (int i = 0; i < 3; i++){
		borg[i] = 0.0;
		bdire[i] = 0.0;
		for (int j = 0; j < 3; j++){
			borg[i] += Rot[j][i] * ppos[j];
			bdire[i] += Rot[j][i] * pdire[j];
		}
		borg[i] = borg[i]/Zoom + Gpos[i];
	}

	std::vector<int> intF;
	iAABB->search_intersect_triangle_with_line(borg, bdire, &intF);

	int minN = -1;
	float minL;
	for (int i = 0; i < intF.size(); i++){
		float fps[3][3], cord[3];
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				fps[j][k] = imesh->Nodes[imesh->Faces[intF[i]][j]][k]; 
			}
		}
		intersection_of_triangle_with_line(borg, bdire, fps, cord);
		if (cord[0] >= 0.0 && cord[1] >= 0.0 && cord[2] >= 0.0 && cord[0] + cord[1] <= 1.0){
			float intp[3];
			float seki = 0.0;
			for (int j = 0; j < 3; j++){
				intp[j] = cord[0] * fps[0][j] + cord[1] * fps[1][j] + cord[2] * fps[2][j];
				seki += (intp[j] - borg[j])*bdire[j];
			}
			if (minN == -1 || minL > seki){
				minN = i; 
				minL = seki;
				cnode->fn = intF[i];
				cnode->uv[0] = cord[1];
				cnode->uv[1] = cord[2];
				for (int j = 0; j < 3; j++){
					cnode->pos[j] = intp[j];
				}
			}
		}
	}


}


void curve_design_display::search_selected_element(int cx, int cy){

	struct curve_node cnode;
	search_clicked_point(cx, cy, &cnode);

	float minPD,minED;
	int minPN = -1;
	int minEN = -1;
	int minPCN,minECN;
	for (int i = 0; i < cNet->curSet.size(); i++){
		class curve_on_mesh *scur = cNet->curSet[i];
		if (myState.ecState != EDGEDIRE || myState.eState == DELCSTATE || myState.eState == DEFCSTATE){
			for (int j = 0; j < scur->nodeN; j++){
				if (myState.eState == DELCSTATE || myState.eState == DEFCSTATE 
					|| (myState.ecState == POINTMOVE || (scur->pCons[j].cType != -1 && j != 0 && j != scur->nodeN-1))){
					float dist = 0.0;
					for (int k = 0; k < 3; k++){
						dist += pow(scur->cNode[j].pos[k] - cnode.pos[k], 2);
					}
					dist = sqrt(dist);

					if (minPN == -1 || dist < minPD){
						minPN = j;
						minPCN = i;
						minPD = dist;
					}
				}
			}
		} 
		if (myState.ecState != POINTMOVE || myState.eState == DELCSTATE){
			for (int j = 0; j < scur->nodeN - 1; j++){
				bool ccheck = true;
				if (scur->pCons[i].cType != -1 && scur->pCons[i + 1].cType != -1) ccheck = false;
				if (myState.eState == DELCSTATE || ((myState.ecState == EDGEDIRE || scur->dCons[j].cType != -1) && ccheck)){
					float ed[3];
					float edh = 0.0;
					for (int k = 0; k < 3; k++){
						ed[k] = scur->cNode[j + 1].pos[k] - scur->cNode[j].pos[k];
						edh += ed[k] * ed[k];
					}
					edh = sqrt(edh);
					if (edh < 1.0e-10) edh = 1.0e-10;
					float seki = 0.0;
					for (int k = 0; k < 3; k++){
						ed[k] /= edh;
						seki += ed[k] * (cnode.pos[k] - scur->cNode[j].pos[k]);
					}
					float dist = 0.0;
					if (seki < 0){
						for (int k = 0; k<3; k++){
							dist += pow(scur->cNode[j].pos[k] - cnode.pos[k], 2);
						}
					}
					else if (seki > edh){
						for (int k = 0; k < 3; k++){
							dist += pow(scur->cNode[j + 1].pos[k] - cnode.pos[k], 2);
						}
					}
					else {
						for (int k = 0; k < 3; k++){
							dist += pow(cnode.pos[k] - scur->cNode[j].pos[k] - seki*ed[k], 2);
						}
					}
					dist = sqrt(dist);

					if (minEN == -1 || dist < minED){
						minEN = j;
						minECN = i;
						minED = dist;
					}
				}
			}
		}
	}

	if (myState.eState == DEFCSTATE){
		minEN = -1;
		for (int i = 0; i < cNet->Connect.size(); i++){
			class curve_on_mesh *cur1 = cNet->curSet[cNet->Connect[i].c1];
			class curve_on_mesh *cur2 = cNet->curSet[cNet->Connect[i].c2];
			int n1 = cNet->Connect[i].n1;
			int n2 = cNet->Connect[i].n2;
			float ed[3];
			float edh = 0.0;
			for (int j = 0; j < 3; j++){
				ed[j] = cur2->cNode[n2].pos[j] - cur1->cNode[n1].pos[j];
				edh += ed[j] * ed[j];
			}
			edh = sqrt(edh);
			if (edh < 1.0e-10) edh = 1.0e-10;
			float seki = 0.0;
			for (int j = 0; j < 3; j++){
				ed[j] /= edh;
				seki += ed[j] * (cnode.pos[j] - cur1->cNode[n1].pos[j]);
			}
			float dist = 0.0;
			if (seki < 0){
				for (int j = 0; j < 3; j++){
					dist += pow(cnode.pos[j] - cur1->cNode[n1].pos[j], 2);
				}
			}
			else if (seki > edh){
				for (int j = 0; j < 3; j++){
					dist += pow(cnode.pos[j] - cur2->cNode[n2].pos[j], 2);
				}
			}
			else {
				for (int j = 0; j < 3; j++){
					dist += pow(cnode.pos[j] - cur1->cNode[n1].pos[j] - seki*ed[j], 2);
				}
			}
			dist = sqrt(dist);
			if (minEN == -1 || dist < minED){
				minEN = i;
				minED = dist;
			}
			
		}
	}

	if (myState.eState == EDITCSTATE){
		if (myState.ecState == POINTMOVE && minPD < cNet->sParam->Wid*1.5){
			selectC = minPCN;
			selectN = minPN;
		}
		if (myState.ecState == EDGEDIRE && minED < cNet->sParam->Wid*1.5){
			selectC = minECN;
			selectE = minEN;
		}
		if (myState.ecState == DELETECONS){
			if (minPN != -1 && minPD < cNet->sParam->Wid*1.5){
				selectC = minPCN;
				selectN = minPN;
				selectE = -1;
			}
			else if (minEN != -1 && minED < cNet->sParam->Wid*1.5){
				selectC = minECN;
				selectE = minEN;
				selectN = -1;
			}
		}
	}

	if (myState.eState == DELCSTATE){
		if (minPN != -1 && minPD < cNet->sParam->Wid*1.5){
			selectC = minPCN;
			selectN = minPN;
			selectE = -1;
		}
		else if (minEN != -1 && minED < cNet->sParam->Wid*1.5){
			selectC = minECN;
			selectE = minEN;
			selectN = -1;
		}
	}

	if (myState.eState == DEFCSTATE){
		if (selem1[0] == -1 && minEN != -1 && minED < cNet->sParam->Wid*1.5){
			sconn = minEN;
		}
		if (minPD < cNet->sParam->Wid*1.5){
			if (selem1[0] == -1){
				selem1[0] = minPCN;
				selem1[1] = minPN;
			}
			else {
				selem2[0] = minPCN;
				selem2[1] = minPN;
			}
		}
	}

}


void curve_design_display::update_cutting_plane(int dx, int dy){

	float pd[3];
	pd[0] = 2.0*dx / (float)height();
	pd[1] = -2.0*dy / (float)height();
	pd[2] = 0.0;

	float Rot[3][3], Trans[3];
	float quat[4];
	quat[0] = rw;
	quat[1] = rx;
	quat[2] = ry;
	quat[3] = rz;
	derive_rot_mat_from_quaternion(quat, Rot);

	float rpd[3];
	float snd[3];
	for (int i = 0; i < 3; i++){
		rpd[i] = 0.0;
		for (int j = 0; j < 3; j++){
			rpd[i] += Rot[j][i] * pd[j];
		}
		snd[i] = Rot[2][i];
		rpd[i] /= Zoom;
	}
	float ed[3];
	float edh = 0.0;
	for (int i = 0; i < 3; i++){
		ed[i] = Cends[1].pos[i] - Cends[0].pos[i];
		edh += ed[i] * ed[i];
	}
	edh = sqrt(edh);
	for (int i = 0; i < 3; i++){
		ed[i] /= edh;
	}

	float wxr[3];
	float seki = 0.0;
	float wh = 0.0;
	for (int i = 0; i < 3; i++){
		wxr[i] = ed[(i + 1) % 3] * rpd[(i + 2) % 3] - ed[(i + 2) % 3] * rpd[(i + 1) % 3];
		wh += wxr[i] * wxr[i];
		seki += -snd[i] * wxr[i];
	}
	float theta = wh / edh;
	if (theta > M_PI*0.3) theta = M_PI*0.3;
	if (seki < 0.0) theta *= -1.0;

	float nd[3];
	quat[0] = cos(theta/2.0);
	for (int i = 0; i < 3; i++){
		quat[1 + i] = sin(theta / 2.0)*ed[i];
		nd[i] = cutNorm[i];
	}
	derive_rot_mat_from_quaternion(quat, Rot);
	for (int i = 0; i < 3; i++){
		cutNorm[i] = 0.0;
		for (int j = 0; j < 3; j++){
			cutNorm[i] += Rot[i][j] * nd[j];
		}
	}

}


void curve_design_display::update_edge_direction(int cx, int cy, int dx, int dy){

	float pdire[3],ndire[3],ppos[3];
	pdire[0] = (float)dx / (float)height();
	pdire[1] = -(float)dy / (float)height();
	pdire[2] = 0.0;
	ndire[0] = 0.0;
	ndire[1] = 0.0;
	ndire[2] = 1.0;
	ppos[0] = 2.0*cx / (float)height() - (float)width() / (float)height() - shiftX*Zoom;
	ppos[1] = 1.0 - 2.0*cy / (float)height() - shiftY*Zoom;
	ppos[2] = -100.0 - shiftZ*Zoom;

	float Rot[3][3];
	float quat[4];
	quat[0] = rw;
	quat[1] = rx;
	quat[2] = ry;
	quat[3] = rz;
	derive_rot_mat_from_quaternion(quat, Rot);

	float mdire[3],bdire[3];
	float ecent[3];
	for (int i = 0; i < 3; i++){
		bdire[i] = 0.0;
		mdire[i] = 0.0;
		for (int j = 0; j < 3; j++){
			bdire[i] += Rot[j][i] * ndire[j];
			mdire[i] += Rot[j][i] * pdire[j];
		}
		mdire[i] /= Zoom;
		ecent[i] = (cNet->curSet[selectC]->cNode[selectE + 1].pos[i] + cNet->curSet[selectC]->cNode[selectE].pos[i]) / 2.0;
		ecent[i] -= Gpos[i];
	}
	float pcent[3];
	for (int i = 0; i < 3; i++){
		pcent[i] = 0.0;
		for (int j = 0; j < 3; j++){
			pcent[i] += Rot[i][j] * ecent[j];
		}
		pcent[i] *= Zoom;
	}
	pcent[0] += shiftX*Zoom;
	pcent[1] += shiftY*Zoom;
	pcent[2] += shiftZ*Zoom;
	float zval = (ppos[0] - pcent[0])*pdire[1] - (ppos[1] - pcent[1])*pdire[0];

	float theta;
	float seki = 0.0;
	float mleng = 0.0;
	float eleng = 0.0;
	for (int i = 0; i < 3; i++){
		eleng += pow(cNet->curSet[selectC]->cNode[selectE].pos[i] - cNet->curSet[selectC]->cNode[selectE].pos[i], 2);
		seki += bdire[i] * eNorm[i];
		mleng += mdire[i] * mdire[i];
	}
	eleng = sqrt(eleng) / 2.0;
	mleng = sqrt(mleng);

	if (eleng > cNet->sParam->Wid*0.1) theta = mleng / eleng / 2.0;
	else theta = mleng / cNet->sParam->Wid;
	if (theta > 1.0)theta = 1.0;
	if (theta < -1.0) theta = -1.0;
	theta = asin(theta);
	if (theta > M_PI / 6.0) theta = M_PI / 6.0;
	if (zval < 0.0) theta *= -1.0;
	if (seki < 0.0) theta *= -1.0;

	quat[0] = cos(theta / 2.0);
	for (int i = 0; i < 3; i++){
		quat[i + 1] = sin(theta / 2.0)*eNorm[i];
	}
	derive_rot_mat_from_quaternion(quat, Rot);

	float Rdir[3];
	for (int i = 0; i < 3; i++){
		Rdir[i] = 0.0;
		for (int j = 0; j < 3; j++){
			Rdir[i] += Rot[i][j] * cNet->curSet[selectC]->dCons[selectE].pos[j];
		}
	}
	for (int i = 0; i < 3; i++){
		cNet->curSet[selectC]->dCons[selectE].pos[i] = Rdir[i];
	}

}


bool curve_design_display::cross_section_curve_search(class curve_on_mesh *curve){

	float ed[3],pd[3];
	for (int i = 0; i < 3; i++){
		ed[i] = Cends[1].pos[i] - Cends[0].pos[i];
	}
	float pdh = 0.0;
	for (int i = 0; i < 3; i++){
		pd[i] = ed[(i + 1) % 3] * cutNorm[(i + 2) % 3] - ed[(i + 2) % 3] * cutNorm[(i + 1) % 3];
		pdh += pd[i] * pd[i];
	}
	pdh = sqrt(pdh);
	for (int i = 0; i < 3; i++){
		pd[i] /= pdh;
	}

	int cstate;
	int cfn = Cends[0].fn;
	float fps[3][3], icc[3];
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			fps[i][j] = imesh->Nodes[imesh->Faces[cfn][i]][j];
		}
	}
	intersection_of_triangle_with_plane(Cends[0].pos, cutNorm, fps, icc, cstate);

	//スタート点での処理
	std::vector<int> sens;
	std::vector<float> scoef;
	if (cstate == -1) return false;
	if (cstate == 1){
		int ipn = -1;
		for (int i = 0; i < 3; i++){
			if (icc[i] == 0.0) ipn = i;
		}
		if (ipn == -1) return false;
		int en = imesh->FElist[cfn][(ipn + 2) % 3];
		float coef = 1.0;
		if (imesh->Edges[en][0] != imesh->Faces[cfn][ipn]) coef = 0.0;
		sens.push_back(en);
		scoef.push_back(coef);
		en = imesh->FElist[cfn][ipn];
		coef = 1.0;
		if (imesh->Edges[en][0] != imesh->Faces[cfn][ipn]) coef = 0.0;
		sens.push_back(en);
		scoef.push_back(coef);
	}
	else if (cstate == 2){
		int opn = -1;
		for (int i = 0; i < 3; i++){
			if (icc[i]>0.9)opn = i;
		}
		if (opn == -1) return false;
		int en = imesh->FElist[cfn][opn];
		float coef = 0.0;
		if (imesh->Edges[en][0] != imesh->Faces[cfn][opn]) coef = 1.0;
		sens.push_back(en);
		scoef.push_back(coef);
		en = imesh->FElist[cfn][(opn + 2) % 3];
		coef = 0.0;
		if (imesh->Edges[en][0] != imesh->Faces[cfn][opn]) coef = 1.0;
		sens.push_back(en);
		scoef.push_back(coef);
	}
	else if (cstate == 3){
		float ppv[3];
		for (int i = 0; i < 3; i++){
			ppv[i] = 0.0;
			for (int j = 0; j < 3; j++){
				ppv[i] += (fps[i][j] - Cends[0].pos[j])*pd[j];
			}
		}
		int oen = -1;
		int intN = 0;
		for (int i = 0; i < 3; i++){
			float dif = ppv[(i + 1) % 3] - ppv[i];
			if (dif == 0)icc[i] = -ppv[i] * 1000;
			else icc[i] = -ppv[i] / dif;
			if (icc[i] >= -0.0 && icc[i] <= 1.0) intN++;
			else oen = i;
		}
		if (oen == -1 && intN != 2) return false;
		int en = imesh->FElist[cfn][(oen + 1) % 3];
		float coef = icc[(oen + 1) % 3];
		if (imesh->Edges[en][0] == imesh->Faces[cfn][(oen + 1) % 3]) coef = 1.0 - coef;
		sens.push_back(en);
		scoef.push_back(coef);
		en = imesh->FElist[cfn][(oen + 2) % 3];
		coef = icc[(oen + 2) % 3];
		if (imesh->Edges[en][0] == imesh->Faces[cfn][(oen + 2) % 3]) coef = 1.0 - coef;
		sens.push_back(en);
		scoef.push_back(coef);
	}
	else {
		int intN = 0;
		int oen = -1;
		for (int i = 0; i < 3; i++){
			if (icc[i] >= -0.0 && icc[i] <= 1.0) intN++;
			else oen = i;
		}
		if (oen == -1 && intN != 2) return false;
		int en = imesh->FElist[cfn][(oen + 1) % 3];
		float coef = icc[(oen + 1) % 3];
		if (imesh->Edges[en][0] == imesh->Faces[cfn][(oen + 1) % 3]) coef = 1.0 - coef;
		sens.push_back(en);
		scoef.push_back(coef);
		en = imesh->FElist[cfn][(oen + 2) % 3];
		coef = icc[(oen + 2) % 3];
		if (imesh->Edges[en][0] == imesh->Faces[cfn][(oen + 2) % 3]) coef = 1.0 - coef;
		sens.push_back(en);
		scoef.push_back(coef);
	}

	
	//境界追跡
	bool end_is_vertex = false;
	int endp = -1;
	if (Cends[1].uv[0] > 0.9){
		endp = imesh->Faces[Cends[1].fn][1];
		end_is_vertex = true;
	}
	else if (Cends[1].uv[1] > 0.9){
		endp = imesh->Faces[Cends[1].fn][2];
		end_is_vertex = true;
	}
	else if (Cends[1].uv[0] + Cends[1].uv[1] < 0.1){
		endp = imesh->Faces[Cends[1].fn][0];
		end_is_vertex = true;
	}
	bool start_is_vertex = false;
	int startp = -1;
	if (Cends[0].uv[0] > 0.9){
		startp = imesh->Faces[Cends[0].fn][1];
		start_is_vertex = true;;
	}
	else if (Cends[0].uv[1] > 0.9){
		startp = imesh->Faces[Cends[0].fn][2];
		start_is_vertex = true;
	}
	else if (Cends[0].uv[0] + Cends[0].uv[1] < 0.1){
		startp = imesh->Faces[Cends[0].fn][0];
		start_is_vertex = true;
	}


	std::vector<std::vector<int>> Cfset,Ceset;
	std::vector<std::vector<float>> Ccoefs;
	for (int i = 0; i < sens.size(); i++){
		std::vector<int> eset,fset;
		std::vector<float> coefs;
		int cen = sens[i];
		cfn = Cends[0].fn;
		fset.push_back(cfn);
		eset.push_back(sens[i]);
		coefs.push_back(scoef[i]);

		bool reach_other_end = false;
		while (1){
			int nfn = imesh->EFlist[cen][0];
			if (nfn == cfn) nfn = imesh->EFlist[cen][1];
			if (nfn == -1) break;

			int nen;
			float ncoef;
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 3; k++){
					fps[j][k] = imesh->Nodes[imesh->Faces[nfn][j]][k];
				}
			}
			intersection_of_triangle_with_plane(Cends[0].pos, cutNorm, fps, icc, cstate);
			if(cstate != 0) printf("state:%d\n",cstate);
			if (cstate == -1) break; 
			if (cstate == 1){
				int ipn = -1;
				for (int j = 0; j < 3; j++){
					if (icc[j] == 0.0) ipn = j;
				}
				if (ipn == -1)break;
				int ens[2];
				ens[0] = imesh->FElist[nfn][(ipn + 2) % 3];
				ens[1] = imesh->FElist[nfn][ipn];
				nen = ens[0];
				if (ens[0] == cen)nen = ens[1];
				ncoef = 0.0;
				if (imesh->Edges[nen][0] == imesh->Faces[nfn][ipn])ncoef = 1.0;
			}
			else if (cstate == 2){
				int opn = -1;
				for (int j = 0; j < 3; j++){
					if (icc[j]>0.9)opn = j;
				}
				if (opn == -1) break;
				if (cen == imesh->FElist[nfn][opn]){
					nen = imesh->FElist[nfn][(opn + 2) % 3];
					ncoef = 0.0;
					if (imesh->Edges[nen][0] == imesh->Faces[nfn][(opn + 2) % 3]) ncoef = 1.0;
				}
				else if (cen == imesh->FElist[nfn][(opn + 2) % 3]){
					nen = imesh->FElist[nfn][opn];
					ncoef = 0.0;
					if (imesh->Edges[nen][0] == imesh->Faces[nfn][(opn + 1) % 3]) ncoef = 1.0;
				}
				else {
					if (coefs.size() > 1 && (coefs[coefs.size() - 1] < 0.1 || coefs[coefs.size()-1] > 0.9 )){
						int ppn = imesh->Edges[cen][1];
						if (coefs[coefs.size() - 1] > 0.9) ppn = imesh->Edges[cen][0];
						int ppen = eset[coefs.size() - 2];
						int fpn = 1;
						if (ppn != imesh->Faces[nfn][(opn+1)%3]) fpn = 2;
						if (imesh->Edges[ppen][0] == ppn || imesh->Edges[ppen][1] == ppn){
							if (fpn == 1) {
								nen = imesh->FElist[nfn][(opn + 2) % 3];
								ncoef = 0.0;
								if (imesh->Edges[nen][0] == imesh->Faces[nfn][(opn + 2) % 3]) ncoef = 1.0;
							} 
							else {
								nen = imesh->FElist[nfn][opn];
								ncoef = 0.0;
								if (imesh->Edges[nen][0] == imesh->Faces[nfn][(opn + 1) % 3]) ncoef = 1.0;
							}
						}
						else {
							if (fpn == 1){
								nen = imesh->FElist[nfn][opn];
								ncoef = 0.0;
								if (imesh->Edges[nen][0] == imesh->Faces[nfn][(opn + 1) % 3]) ncoef = 1.0;
							}
							else {
								nen = imesh->FElist[nfn][(opn + 2) % 3];
								ncoef = 0.0;
								if (imesh->Edges[nen][0] == imesh->Faces[nfn][(opn + 2) % 3]) ncoef = 1.0;
							}
						}
					}
					else {
						float seki1 = 0.0;
						float seki2 = 0.0;
						for (int j = 0; j < 3; j++){
							seki1 += ed[j] * (imesh->Nodes[imesh->Faces[nfn][(opn + 1) % 3]][j] - Cends[0].pos[j]);
							seki2 += ed[j] * (imesh->Nodes[imesh->Faces[nfn][(opn + 2) % 3]][j] - Cends[0].pos[j]);
						}
						if (seki1 > seki2){
							nen = imesh->FElist[nfn][opn];
							ncoef = 0.0;
							if (imesh->Edges[nen][0] == imesh->Faces[nfn][(opn + 1) % 3]) ncoef = 1.0;
						}
						else {
							nen = imesh->FElist[nfn][(opn + 2) % 3];
							ncoef = 0.0;
							if (imesh->Edges[nen][0] == imesh->Faces[nfn][(opn + 2) % 3]) ncoef = 1.0;
						}
					}
				}
			}
			else if (cstate == 3){
				float sekis[3],lengs[3];
				int maxN = 0;
				float maxseki;
				for (int j = 0; j < 3; j++){
					sekis[j] = 0.0;
					lengs[j] = 0.0;
					for (int k = 0; k < 3; k++){
						sekis[j] += ed[k] * (imesh->Nodes[imesh->Faces[nfn][j]][k] - Cends[0].pos[k]);
						lengs[j] += pd[k] * (imesh->Nodes[imesh->Faces[nfn][j]][k] - Cends[0].pos[k]);
					}
					if (j == 0 || maxseki < sekis[j]){
						maxN = j;
						maxseki = sekis[j];
					}
				}
				if (imesh->FElist[nfn][maxN] == cen){
					nen = imesh->FElist[nfn][(maxN + 2) % 3];
					ncoef = 0.0;
					if (imesh->Edges[nen][0] == imesh->Faces[nfn][maxN])ncoef = 1.0;
				}
				else if (imesh->FElist[nfn][(maxN + 2) % 3] == cen){
					nen = imesh->FElist[nfn][maxN];
					ncoef = 0.0;
					if (imesh->Edges[nen][0] == imesh->Faces[nfn][maxN]) ncoef = 1.0;
				}
				else {
					if (fabs(lengs[(maxN + 1) % 3]) < fabs(lengs[(maxN + 2) % 3])){
						nen = imesh->FElist[nfn][maxN];
						ncoef = 0.0;
						if (imesh->Edges[nen][0] == imesh->Faces[nfn][maxN]) ncoef = 1.0;
					}
					else {
						nen = imesh->FElist[nfn][(maxN + 2) % 3];
						ncoef = 0.0;
						if (imesh->Edges[nen][0] == imesh->Faces[nfn][maxN])ncoef = 1.0;
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
				if (imesh->FElist[nfn][(cpn + 2) % 3] == cen){
					nen = imesh->FElist[nfn][cpn];
					ncoef = icc[cpn];
					if (imesh->Edges[nen][0] == imesh->Faces[nfn][cpn]) ncoef = 1.0 - ncoef;
				}
				else if (imesh->FElist[nfn][cpn] == cen){
					nen = imesh->FElist[nfn][(cpn + 2) % 3];
					ncoef = icc[(cpn + 2) % 3];
					if (imesh->Edges[nen][1] == imesh->Faces[nfn][cpn]) ncoef = 1.0 - ncoef;
				}
				else {
					break;
				}
			}

			//更新・終端チェック
			if (Cends[1].fn == nfn){
				reach_other_end = true;
				break;
			}
			else if (end_is_vertex){
				bool reach_neib = false;
				for (int j = 0; j < imesh->PFsize[endp]; j++){
					if (nfn == imesh->PFlist[endp][j]) reach_neib = true;
				}
				if (reach_neib){
					reach_other_end = true;
					break;
				}
			}
			if (fset.size() > imesh->faceN/4){
				break;
			}
			else if (fset.size() > 10){
				if (nfn == Cends[0].fn) break;
				if (start_is_vertex){
					bool return_neib = false;
					for (int j = 0; j < imesh->PFsize[startp]; j++){
						if (nfn == imesh->PFlist[startp][j]) return_neib = true;
					}
					if (return_neib) break;
				}
			}

			if(cstate != 0) printf("%d %d %f\n",nfn,cen,ncoef);
			cfn = nfn;
			cen = nen;
			fset.push_back(cfn);
			eset.push_back(cen);
			coefs.push_back(ncoef);
		}

		if (reach_other_end){
			Cfset.push_back(fset);
			Ceset.push_back(eset);
			Ccoefs.push_back(coefs);
		}
	}
	if (Cfset.size() == 0) return false;

	int cNum = -1;
	float maxRatio;
	for (int i = 0; i < Cfset.size(); i++){
		float hleng = 0.0;
		float pleng = 0.0;
		float p1[3], p2[3];
		for (int j = 0; j < 3; j++){
			p1[j] = Cends[0].pos[j];
		}
		for (int j = 0; j < Ceset[i].size(); j++){
			int en = Ceset[i][j];
			float coef = Ccoefs[i][j];
			float leng = 0.0;
			float seki = 0.0;
			for (int k = 0; k < 3; k++){
				p2[k] = coef*imesh->Nodes[imesh->Edges[en][0]][k] + (1.0 - coef)*imesh->Nodes[imesh->Edges[en][1]][k];
				leng += pow(p1[k] - p2[k], 2);
				seki += pd[k] * ((p1[k] + p2[k]) / 2.0 - Cends[0].pos[k]);
			}
			leng = sqrt(leng);
			hleng += leng;
			if (seki >= 0.0) pleng += leng;
			for (int k = 0; k < 3; k++){
				p1[k] = p2[k];
			}
		}
		float ratio = pleng / hleng;
		if (i == 0 || ratio > maxRatio){
			cNum = i;
			maxRatio = ratio;
		}
	}

	curve->nodeN = Cfset[cNum].size() + 2;
	curve->cNode = new struct curve_node[curve->nodeN];
	int cnNum = curve->nodeN;
	struct curve_node *cnSet = curve->cNode;
	cnSet[0].fn = Cends[0].fn;
	for (int j = 0; j < 3; j++){
		cnSet[0].pos[j] = Cends[0].pos[j];
		if (j != 2) cnSet[0].uv[j] = Cends[0].uv[j];
	}
	for (int j = 0; j < Cfset[cNum].size(); j++){
		cnSet[1 + j].fn = Cfset[cNum][j];
		int fp1 = -1;
		int fp2 = -1;
		int en = Ceset[cNum][j];
		for (int k = 0; k < 3; k++){
			if (imesh->Faces[cnSet[1 + j].fn][k] == imesh->Edges[en][0]) fp1 = k;
			if (imesh->Faces[cnSet[1 + j].fn][k] == imesh->Edges[en][1]) fp2 = k;
		}
		if (fp1 == -1 || fp2 == -1) printf("okashii@cn search\n");
		for (int k = 0; k < 3; k++){
			cnSet[1 + j].pos[k] = Ccoefs[cNum][j] * imesh->Nodes[imesh->Edges[en][0]][k] +
				(1.0 - Ccoefs[cNum][j])*imesh->Nodes[imesh->Edges[en][1]][k];
		}
		float uvw[3];
		uvw[fp1] = Ccoefs[cNum][j];
		uvw[fp2] = 1.0 - uvw[fp1];
		uvw[3 - fp1 - fp2] = 0.0;
		cnSet[1 + j].uv[0] = uvw[1];
		cnSet[1 + j].uv[1] = uvw[2];
	}
	cnSet[cnNum - 1].fn = Cends[1].fn;
	for (int j = 0; j < 3; j++){
		cnSet[cnNum - 1].pos[j] = Cends[1].pos[j];
		if (j != 2) cnSet[cnNum - 1].uv[j] = Cends[1].uv[j];
	}


	//リサンプリング
	float hleng = 0.0;
	for (int i = 0; i < curve->nodeN - 1; i++){
		float leng = 0.0;
		for (int j = 0; j < 3; j++){
			leng += pow(curve->cNode[i].pos[j] - curve->cNode[i + 1].pos[j], 2);
		}
		leng = sqrt(leng);
		hleng += leng;
	}
	int divN = (int)((float)hleng / cNet->sParam->Wid / 3.0);
	if (divN < 5) divN = 5;
	if (divN > 200) divN = 200;
	float segL = hleng / (float)divN;

	std::vector<struct curve_node> rcSet;
	rcSet.push_back(curve->cNode[0]);
	hleng = 0.0;
	int segN = 0;
	for (int i = 0; i < curve->nodeN - 1; i++){
		float leng = 0.0;
		for (int j = 0; j < 3; j++){
			leng += pow(curve->cNode[i].pos[j] - curve->cNode[i + 1].pos[j], 2);
		}
		leng = sqrt(leng);
		hleng += leng;
		int sN = (int)((float)hleng / (float)segL);
		if (sN != segN){
			for (int j = segN; j < sN; j++){
				float kv = (float)(j + 1)*segL - hleng + leng;
				kv = kv / leng;
				float pos[3];
				for (int k = 0; k < 3; k++){
					pos[k] = (1.0 - kv)*curve->cNode[i].pos[k] + kv*curve->cNode[i + 1].pos[k];
				}
				float cord[3];
				int fn;
				iAABB->closest_point_search_with_mesh_info(pos, &fn, cord);
				struct curve_node cnode;
				cnode.fn = fn;
				cnode.uv[0] = cord[1];
				cnode.uv[1] = cord[2];
				for (int k = 0; k < 3; k++){
					cnode.pos[k] = 0.0;
					for (int l = 0; l < 3; l++){
						cnode.pos[k] += cord[l] * imesh->Nodes[imesh->Faces[fn][l]][k];
					}
				}
				rcSet.push_back(cnode);
			}
		}
		segN = sN;
	}
	delete[] curve->cNode;

	curve->nodeN = rcSet.size();
	curve->cNode = new struct curve_node[rcSet.size()];

	for (int i = 0; i < curve->nodeN; i++){
		curve->cNode[i].fn = rcSet[i].fn;
		for (int j = 0; j < 3; j++){
			curve->cNode[i].pos[j] = rcSet[i].pos[j];
			if (j != 2) curve->cNode[i].uv[j] = rcSet[i].uv[j];
		}
	}

	return true;
}