#include <QtGui>
#include <qdialog.h>
#include <qfiledialog.h>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <qmessagebox.h>
#include <QSplitter>

#include "mainwindow.h"
#include "curve_design_display.h"
#include "SimuRodDisplay.h"
#include "geometric_data_structure.h"
#include "EditTools.h"
#include "AABB_data_structure.h"
#include "curve_network.h"
#include "RodNetwork.h"
#include "MaterialSettingDialog.h"
#include "MomentReductionMap.h"
#include "Setting.h"

//コンストラクタ
MainWindow::MainWindow(QWidget *parent): QMainWindow(parent){

	sParam.Wid = 5.0;
	sParam.Hei = 2.0;
	sParam.Density = 8.12e-7;
	sParam.Es = 1.0e7;
	sParam.Gt = 1.0e5;
	sParam.Eb11 = 1.0e7;
	sParam.Eb22 = 1.0e7;
	sParam.Eb12 = 0.0;

	sParam.sampL[0] = sParam.Wid*2.0;
	sParam.sampL[1] = sParam.Wid*2.0;
	sParam.sampL[2] = sParam.Wid*3.0;
	sParam.sampL[3] = sParam.Wid*3.0;

	imesh = new loop_subdivision_mesh;
	iAABB = new AABB_struct;
	cNet = new CurveNetwork;
	rNet = new RodNetwork;
	mrMap = new MomentReductionMap;

	cNet->sParam = &sParam;
	rNet->sParam = &sParam;
	

	installEventFilter(this);

	Cedit = new curve_design_display(this);
	Cedit->imesh = imesh;
	Cedit->cNet = cNet;
	Cedit->mrMap = mrMap;
	Cedit->resize(300, 300);
	Cedit->installEventFilter(this);

	sDisp = new SimuRodDisplay(this);
	sDisp->imesh = imesh;
	sDisp->rNet = rNet;
	sDisp->resize(100, 100);
	sDisp->installEventFilter(this);

	Etools = new EditTools(this);

	split2 = new QSplitter(Qt::Vertical, this);
	split2->addWidget(Etools);
	split2->addWidget(sDisp);
	split = new QSplitter(Qt::Horizontal, this);
	split->addWidget(Cedit);
	//split->addWidget(Etools);
	split->addWidget(split2);
	setCentralWidget(split);

	createActions();
	createMenus();

	myState.MESH_EXIST = false;
	myState.eState = PREPSTATE;
	myState.acState = POINTZERO;
	myState.sState = SIDLESTATE;
	sync_state(myState);
}


//デストラクタ
MainWindow::~MainWindow(){

	//imesh->~loop_subdivision_mesh();
	//iAABB->~AABB_struct();
	//cNet->~CurveNetwork();

	delete imesh, iAABB, cNet;
}


////////
//設定・準備
////////


bool MainWindow::eventFilter(QObject *target, QEvent *event){
	
	if (event->type() == QEvent::KeyPress){
		QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
		Cedit->Pressed_Key = (Qt::Key)keyEvent->key();
		sDisp->Pressed_Key = (Qt::Key)keyEvent->key();
	}
	else if (event->type() == QEvent::KeyRelease){
		QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
		Cedit->Pressed_Key = Qt::Key_F35;
		sDisp->Pressed_Key = Qt::Key_F35;
	}

	return QMainWindow::eventFilter(target, event);
}


void MainWindow::createActions(){

	openAct = new QAction(tr("&Open"), this);
	openAct->setStatusTip(tr("Open Manifold mesh file"));
	connect(openAct, SIGNAL(triggered()), this, SLOT(open()));

	loadCurvAct = new QAction(tr("&Curves"), this);
	loadCurvAct->setStatusTip(tr("Load Curve Network"));
	connect(loadCurvAct, SIGNAL(triggered()), this, SLOT(loadCurveNet()));

	expCurvAct = new QAction(tr("&Curves"), this);
	expCurvAct->setStatusTip(tr("Export Curve Network"));
	connect(expCurvAct, SIGNAL(triggered()), this, SLOT(expCurveNet()));

	expRodAct = new QAction(tr("&Rods"), this);
	expRodAct->setStatusTip(tr("Export Rod NetWork"));
	connect(expRodAct, SIGNAL(triggered()), this, SLOT(expRodNet()));


	materialAct = new QAction(tr("&Material"), this);
	materialAct->setStatusTip(tr("Set Material Parameter"));
	connect(materialAct, SIGNAL(triggered()), this, SLOT(setMaterialParam()));
}


void MainWindow::createMenus(){
	
	fileMenu = menuBar()->addMenu(tr("&File"));
	fileMenu->addAction(openAct);

	expSubMenu = fileMenu->addMenu(tr("&Export"));
	expSubMenu->addAction(expCurvAct);
	expSubMenu->addAction(expRodAct);

	loadSubMenu = fileMenu->addMenu(tr("&Load"));
	loadSubMenu->addAction(loadCurvAct);

	setMenu = menuBar()->addMenu(tr("&Setting"));
	setMenu->addAction(materialAct);

}


void MainWindow::open(){

	QString fileName = QFileDialog::getOpenFileName(this, tr("Open manifold mesh"), "./Input", tr("Mesh files (*.obj)"));
	if (!fileName.isEmpty()){
		if (loadFile(fileName)) {
			myState.MESH_EXIST = true;
			myState.eState = PREPSTATE;
			sync_state(myState);
		}
	}
}


void MainWindow::expCurveNet(){

	QString fileName = QFileDialog::getSaveFileName(this, tr("Export Curve Network"), "./Output", tr("Curve Network file (*.cn)"));

	if (!fileName.isEmpty()){
		QFile file(fileName);
		bool success = cNet->exportCurveNet((std::string)file.fileName().toLocal8Bit(),iAABB);
		if (!success){
			QMessageBox::warning(this, tr("CurveNetwork"),
				tr("No Curve or No connection"));
		}
	}
}


void MainWindow::expRodNet(){

	if (cNet->curSet.empty() || cNet->Connect.empty()){
		QMessageBox::warning(this, tr("CurveNetwork"),
			tr("No Curve or No connection"));
		return;
	}

	QString fileName = QFileDialog::getSaveFileName(this, tr("Export Rod Network"), "./Output", tr("Rod Network file (*.rn)"));

	if (!fileName.isEmpty()){
		class RodNetwork ORnet;
		ORnet.sParam = &sParam;
		ORnet.curve_to_rod(cNet, iAABB, imesh);
		ORnet.exportRodNetwork((std::string)fileName.toLocal8Bit(),imesh);
	}


}


void MainWindow::loadCurveNet(){

	QString fileName = QFileDialog::getOpenFileName(this, tr("Load Curve Network"), "./Input", tr("Curve Network files (*.cn)"));
	if (!fileName.isEmpty()){
		QFile file(fileName);
		if (!file.open(QIODevice::ReadOnly)){
			QMessageBox::warning(this, tr("Curve Network"), tr("Cannot read file %1:\n%2.")
				.arg(file.fileName()).arg(file.errorString()));
			return;
		}

		imesh->~loop_subdivision_mesh();
		imesh = new loop_subdivision_mesh;
		Cedit->imesh = imesh;
		sDisp->imesh = imesh;
		cNet->~CurveNetwork();
		cNet = new CurveNetwork;
		cNet->sParam = &sParam;
		Cedit->cNet = cNet;
		bool success = cNet->loadCurveNet((std::string)file.fileName().toLocal8Bit(),imesh);
		Cedit->reset_param_for_mesh();
		sDisp->reset_param_for_mesh();

		iAABB->~AABB_struct();
		iAABB = new AABB_struct;
		iAABB->build_AABB_struct(imesh);
		iAABB->update_AABB_range();
		Cedit->iAABB = iAABB;

		if (success){
			myState.MESH_EXIST = true;
			if(!cNet->curSet.empty())myState.eState = IDLESTATE;
			else myState.eState = PREPSTATE;
			sync_state(myState);
		}
	}


}


bool MainWindow::loadFile(const QString &fileName){

	QFile file(fileName);
	if (!file.open(QIODevice::ReadOnly)){
		QMessageBox::warning(this, tr("Mesh"), tr("Cannot read file %1:\n%2.")
			.arg(file.fileName()).arg(file.errorString()));
		return false;
	}

	imesh->~loop_subdivision_mesh();
	imesh = new loop_subdivision_mesh;
	Cedit->imesh = imesh;
	sDisp->imesh = imesh;
	cNet->~CurveNetwork();
	cNet = new CurveNetwork;
	Cedit->cNet = cNet;
	cNet->sParam = &sParam;

	read_obj_file((std::string)file.fileName().toLocal8Bit(), imesh);
	imesh->construct_data_structure_without_edge();
	Cedit->reset_param_for_mesh();
	sDisp->reset_param_for_mesh();
	imesh->prepare_for_evaluation();

	iAABB->~AABB_struct();
	iAABB = new AABB_struct;
	iAABB->build_AABB_struct(imesh);
	iAABB->update_AABB_range();
	Cedit->iAABB = iAABB;
	cNet->lmesh = imesh;

	return true;
}


///////
///ウィジェット間の連絡用
///////


void MainWindow::addCurve(){

	if (myState.eState != ADDCSTATE){
		myState.eState = ADDCSTATE;
		myState.acState = POINTZERO;
	}
	else {
		if(cNet->curSet.size() == 0) myState.eState = PREPSTATE;
		else myState.eState = IDLESTATE;
		myState.acState = POINTZERO;
	}
	
	sync_state(myState);

	//Etools->toggle_for_add_curve();

}


void MainWindow::setCurveEnds(){

	if (myState.acState == POINTTWO){
		myState.acState = PLANESEARCH;
	}

	sync_state(myState);

}


void MainWindow::setCuttingPlane(){

	Cedit->addCurveToNet();
	myState.eState = PREPSTATE;
	if (cNet->curSet.size() > 0) myState.eState = IDLESTATE;
	myState.acState = POINTZERO;
	sync_state(myState);

	myState.eState = ADDCSTATE;
	sync_state(myState);
}


void MainWindow::editCurve(){

	if (myState.eState != EDITCSTATE){
		myState.eState = EDITCSTATE;
		myState.ecState = POINTMOVE;
	}
	else {
		if (cNet->curSet.size() == 0) myState.eState = PREPSTATE;
		else myState.eState = IDLESTATE;
		myState.ecState = POINTMOVE;
	}

	sync_state(myState);
}


void MainWindow::editConnect(){

	if (myState.eState != DEFCSTATE){
		myState.eState = DEFCSTATE;
	}
	else {
		if (cNet->curSet.size() == 0) myState.eState = PREPSTATE;
		else myState.eState = IDLESTATE;
	}
	sync_state(myState);

}


void MainWindow::deleteCurve(){

	if (myState.eState != DELCSTATE){
		myState.eState = DELCSTATE;
	}
	else {
		if (cNet->curSet.size() == 0) myState.eState = PREPSTATE;
		else myState.eState = IDLESTATE;
	}
	sync_state(myState);

}


void MainWindow::simulateRodNet(){

	bool Tsucceed = rNet->curve_to_rod(cNet, iAABB, imesh);
	if (!Tsucceed) return;

	rNet->prepare_for_simulation();
	rNet->simulateFull();
	rNet->PointsMatching();

	myState.sState = SIMUSTATE;
	sync_state(myState);

}


void MainWindow::stabilize(){

	RodNetwork srNet;
	bool Tsucceed = srNet.curve_to_rod(cNet, iAABB, imesh);
	if (!Tsucceed) return;

	srNet.sParam = &sParam;

	mrMap->lAABB = iAABB;
	mrMap->lmesh = imesh;
	mrMap->rNet = &srNet;
	mrMap->generateMomentReductionMap();


}


void MainWindow::optimizeRodParam(){


}


void MainWindow::setMaterialParam(){

	MaterialSettingDialog MatSetDiag(&sParam, this);
	if (!MatSetDiag.exec()){
		return;
	}

}


void MainWindow::editStateChanged(int estate){

	if (estate == 0) myState.ecState = POINTMOVE;
	else if (estate == 1) myState.ecState = EDGEDIRE;
	else myState.ecState = DELETECONS;

	sync_state(myState);
}


void MainWindow::sync_state(struct state_param state){
	
	myState = state;

	Etools->state_changed(state);
	Cedit->state_changed(state);
	sDisp->state_changed(state);
}


void MainWindow::sync_state_called(struct state_param state){

	sync_state(state);

}