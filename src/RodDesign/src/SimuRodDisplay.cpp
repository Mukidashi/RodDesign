#include <QtGui>
#include <QtWidgets>
#include <QtOpenGL>
#include <time.h>

#include "SimuRodDisplay.h"
#include "RodNetwork.h"
#include "geometric_data_structure.h"
#include "basic_geometric_calculation.h"
#include "DisplayFunctions.h"


//コンストラクタ
SimuRodDisplay::SimuRodDisplay(QWidget *parent) :QGLWidget(parent)
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

	rw = 0.0;
	rx = 1.0;
	ry = 0.0;
	rz = 0.0;

	simuTimer = new QTimer(this);
	connect(simuTimer, SIGNAL(timeout()), this, SLOT(simuIdleFunc()));

	connect(this, SIGNAL(sync_state_call(struct state_param)), parent, SLOT(sync_state_called(struct state_param)));


}


/////////
//設定・準備
/////////

void SimuRodDisplay::initializeGL(){

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

}


void SimuRodDisplay::resizeGL(int width, int height){

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	GLfloat x = GLfloat(width) / height;
	//glFrustum(-x, x, -1.0, 1.0, 4.0, 15.0);
	glOrtho(-x, x, -1.0, 1.0, -100, 100);
	glMatrixMode(GL_MODELVIEW);

}


void SimuRodDisplay::paintGL(){

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	draw();

}


void SimuRodDisplay::reset_param_for_mesh(){
	
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


QSize SimuRodDisplay::sizeHint() const{
	
	return QSize(500, 500);

}


////////////
//イベント関連
////////////


void SimuRodDisplay::mousePressEvent(QMouseEvent *event){
	
	lastPos = event->pos();

}


void SimuRodDisplay::mouseMoveEvent(QMouseEvent *event){

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


void SimuRodDisplay::mouseReleaseEvent(QMouseEvent *event){

	updateGL();

}


void SimuRodDisplay::simuIdleFunc(){

	clock_t start, end;
	start = clock();

	bool is_converged = false;
	bool is_diverged = false;
	bool tstep_is_large = false;
	rNet->simulate_for_animation(TIMESTEP,tstep_is_large, is_converged, is_diverged);

	end = clock();
	double simu_t = (double)(end - start) / CLOCKS_PER_SEC;

	if (is_diverged){
		TIMESTEP *= 0.1;
	}
	else {
		ANITIME += TIMESTEP;
		FRAMECOUNT++;
	}
	CALTIME += simu_t;
	
	if (tstep_is_large && TIMESTEP >1.0e-3){
		TIMESTEP *= 0.1;
	}
	if (!tstep_is_large) SHORTCOUNT++;
	if (SHORTCOUNT == 10){
		if (TIMESTEP < 1.0) TIMESTEP *= 2.0;
		SHORTCOUNT = 0;
	}

	printf("%d %f %f\n\n", FRAMECOUNT, TIMESTEP, ANITIME);

	//収束または発散時の処理
	//if (FRAMECOUNT == 50){
	//	simuTimer->stop();
	//	//simuTimer->start(1000);
	//	//TIMESTEP = 0.01;
	//}
	if (FRAMECOUNT == 500 || is_converged) {
	//if (FRAMECOUNT == 20) {
		simuTimer->stop();
		myState.sState = SIDLESTATE;
		emit sync_state_call(myState);
	}
	if (is_diverged && TIMESTEP <1.0e-4){
		simuTimer->stop();
		myState.sState = SIDLESTATE;
		emit sync_state_call(myState);
	}

	if (CALTIME > 1.0 / 30.0) {
		updateGL();
		CALTIME = 0.0;
	}
}


////////////
//表示関連
////////////


void SimuRodDisplay::draw(){

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
	if (rh != 0.0) glRotatef(2.0*acos(rw)*180.0 / M_PI, rx / rh, ry / rh, rz / rh);

	display_mesh();
	display_rod_network();

}


void SimuRodDisplay::display_mesh(){

	//glPolygonMode(GL_FRONT, GL_FILL);
	glDisable(GL_DEPTH_TEST);
	//glEnable(GL_LIGHTING);
	glDisable(GL_LIGHTING);
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
			norm[j] /= nh / Zoom;
		}
		float cval[3];
		cval[0] = 0.98; cval[1] = 0.98; cval[2] = 0.98;
		//cval[0] = 180.0 / 256.0; cval[1] = 200 / 256.0; cval[2] = 220.0 / 256.0;
		choose_color_data(8);
		set_plastic_color(cval);
		glColor3f(180.0 / 256.0, 200 / 256.0, 220.0 / 256.0);
		glBegin(GL_POLYGON);
		glNormal3fv(norm);
		glVertex3fv(fps[0]);
		glVertex3fv(fps[1]);
		glVertex3fv(fps[2]);
		glEnd();
	}
	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);

}


void SimuRodDisplay::display_rod_network(){

	glEnable(GL_LIGHTING);
	for (int i = 0; i < rNet->rodN; i++){
		class DiscreteRod *rod = &rNet->Rods[i];
		choose_color_data(3);
		disp_rod_polygon(rod, Zoom);
	}

	for (int i = 0; i < rNet->conN; i++){
		for (int j = 0; j < rNet->Connect[i].nodeN; j++){
			int rN = rNet->Connect[i].cNode[j].rN;
			int nN = rNet->Connect[i].cNode[j].nN;
			float pos[3];
			for (int k = 0; k < 3; k++){
				pos[k] = rNet->Rods[rN].Nodes[nN][k];
			}
			float cval[3];
			cval[0] = 0.9; cval[1] = 0.0; cval[2] = 0.2;
			set_plastic_color(cval);
			disp_sphere(pos, rNet->sParam->Wid*0.45, Zoom);
		}
	}

}


//////////
///その他
/////////


void SimuRodDisplay::state_changed(struct state_param state){

	if (state.sState == SIMUSTATE && myState.sState != SIMUSTATE){
		FRAMECOUNT = 0;
		TIMESTEP = 0.1;
		ANITIME = 0.0;
		CALTIME = 0.0;
		SHORTCOUNT = 0;
		//simuTimer->start(0);
		//rNet->simulateFull();
		myState.sState = SIDLESTATE;
		emit sync_state_call(myState);
	}

	myState = state;

	updateGL();

}
