#include<QPushButton>
#include<QLabel>
#include<qboxlayout.h>
#include<QGroupBox>
#include<qradiobutton.h>
#include "EditTools.h"

EditTools::EditTools(QWidget *parent) : QFrame(parent){
	
	//曲線追加ツール用
	addCButton = new QPushButton(tr("Add Curve"));
	addCButton->setEnabled(false);
	connect(addCButton, SIGNAL(clicked()), parent, SLOT(addCurve()));

	addText1 = new QLabel(tr("1. Set End Points"));
	//addText1->hide();
	agrButton1 = new QPushButton(tr("Confirm"));
	agrButton1->setEnabled(false);
	//agrButton1->hide();
	connect(agrButton1, SIGNAL(clicked()), parent, SLOT(setCurveEnds()));
	addText2 = new QLabel(tr("2. Set Cutting Plane"));
	//addText2->hide();
	agrButton2 = new QPushButton(tr("Confirm"));
	//agrButton2->hide();
	agrButton2->setEnabled(false);
	connect(agrButton2, SIGNAL(clicked()), parent, SLOT(setCuttingPlane()));

	QVBoxLayout *layout1 = new QVBoxLayout;
	layout1->addWidget(addCButton);
	layout1->addWidget(addText1);
	layout1->addWidget(agrButton1);
	layout1->addWidget(addText2);
	layout1->addWidget(agrButton2);
	addCBox = new QGroupBox;
	addCBox->setLayout(layout1);

	//曲線編集用
	editCButton = new QPushButton(tr("Edit Curve"));
	editCButton->setEnabled(false);
	connect(editCButton, SIGNAL(clicked()), parent, SLOT(editCurve()));

	posRadio = new QRadioButton(tr("Edit Vertex Posion"));
	//posRadio->hide();
	connect(posRadio, SIGNAL(clicked()), this, SLOT(toggled_edit_radiobox()));
	direRadio = new QRadioButton(tr("Edit Edge Direction"));
	//direRadio->hide();
	connect(direRadio, SIGNAL(clicked()), this, SLOT(toggled_edit_radiobox()));
	delRadio = new QRadioButton(tr("Delete Constraint"));
	//delRadio->hide();
	connect(delRadio, SIGNAL(clicked()), this, SLOT(toggled_edit_radiobox()));

	QVBoxLayout *layout2 = new QVBoxLayout;
	layout2->addWidget(editCButton);
	layout2->addWidget(posRadio);
	layout2->addWidget(direRadio);
	layout2->addWidget(delRadio);
	editCBox = new QGroupBox;
	editCBox->setLayout(layout2);
	connect(this, SIGNAL(editStateChanged(int)), parent, SLOT(editStateChanged(int)));

	//曲線削除ツール用
	delCButton = new QPushButton(tr("Delete Curve"));
	delCButton->setEnabled(false);
	connect(delCButton, SIGNAL(clicked()), parent, SLOT(deleteCurve()));


	//コネクションの設定ツール用
	defCButton = new QPushButton(tr("Edit Connection"));
	defCButton->setEnabled(false);
	connect(defCButton, SIGNAL(clicked()), parent, SLOT(editConnect()));
	defText = new QLabel(tr("Select Two Nodes"));
	defText->hide();
	QVBoxLayout *layout3 = new QVBoxLayout;
	layout3->addWidget(defCButton);
	layout3->addWidget(defText);
	defCBox = new QGroupBox;
	defCBox->setLayout(layout3);


	//安定化ツール用
	stbButton = new QPushButton(tr("Stabilize"));
	stbButton->setEnabled(false);
	connect(stbButton, SIGNAL(clicked()), parent, SLOT(stabilize()));

	//シミュレートツール用
	simuButton = new QPushButton(tr("Simulate"));
	simuButton->setEnabled(false);
	connect(simuButton, SIGNAL(clicked()), parent, SLOT(simulateRodNet()));


	//最適化ツール用
	optButton = new QPushButton(tr("Optimize"));
	optButton->setEnabled(false);
	connect(optButton, SIGNAL(clicked()), parent, SLOT(optimizeRodParam()));


	QVBoxLayout *Hlayout = new QVBoxLayout;
	Hlayout->addWidget(addCBox);
	//Hlayout->addWidget(editCBox);
	//Hlayout->addWidget(optPButton);
	Hlayout->addWidget(delCButton);
	//Hlayout->addWidget(defCBox);
	//Hlayout->addWidget(simuButton);
	//Hlayout->addStretch();
	Hlayout->addWidget(defCButton);
	Hlayout->addWidget(simuButton);

	
	QVBoxLayout *Hlayout2 = new QVBoxLayout;
	//Hlayout2->addStretch();
	//Hlayout2->insertLayout(1, Hlayout);
	//Hlayout2->addStretch();
	Hlayout2->addWidget(editCBox);
	Hlayout2->addWidget(stbButton);
	Hlayout2->addWidget(optButton);
	//setLayout(Hlayout2);

	QHBoxLayout *Hlayout3 = new QHBoxLayout;
	Hlayout3->insertLayout(1, Hlayout);
	Hlayout3->insertLayout(2, Hlayout2);

	QVBoxLayout *Hlayout4 = new QVBoxLayout;
	Hlayout4->addStretch();
	Hlayout4->insertLayout(1, Hlayout3);
	Hlayout4->addStretch();
	setLayout(Hlayout4);

	myState.eState = PREPSTATE;
}


void EditTools::toggle_for_add_curve(int sNum){
	
	if(sNum == 1) addCButton->setText("Quit Adding");
	else addCButton->setText("Add Curve");
	/*if(sNum == 1) addText1->show();
	else addText1->hide();
	if (sNum == 1) agrButton1->show();
	else agrButton1->hide();
	if (sNum == 1) addText2->show();
	else addText2->hide();
	if (sNum == 1) agrButton2->show();
	else agrButton2->hide();*/
}


void EditTools::toggle_for_edit_curve(int sNum){

	if (sNum == 1) editCButton->setText("Quit Editting");
	else editCButton->setText("Edit Curve");
	/*if (posRadio->isHidden()) editCButton->setText("Quit Editting");
	if (posRadio->isHidden()) posRadio->show();
	else posRadio->hide();
	if (direRadio->isHidden()) direRadio->show();
	else direRadio->hide();
	if (delRadio->isHidden()) delRadio->show();
	else delRadio->hide();*/
	
}


void EditTools::toggled_edit_radiobox(void){

	int state = 0;
	if (direRadio->isChecked()) state = 1;
	if (delRadio->isChecked()) state = 2;

	emit editStateChanged(state);

}


void EditTools::state_changed(struct state_param state){

	if (state.MESH_EXIST) {
		addCButton->setEnabled(true);
	}
	else {
		addCButton->setEnabled(false);
	}

	if (state.eState == ADDCSTATE){
		if (myState.eState != ADDCSTATE) toggle_for_add_curve(1);
		if (state.acState == POINTTWO){
			agrButton1->setEnabled(true);
		}
		else {
			agrButton1->setEnabled(false);
		}
		if (state.acState == PLANESETTED){
			agrButton2->setEnabled(true);
		}
		else {
			agrButton2->setEnabled(false);
		}
	}
	if (myState.eState == ADDCSTATE && state.eState != ADDCSTATE) toggle_for_add_curve(0);
	

	if (state.eState == EDITCSTATE){
		if (myState.eState != EDITCSTATE){
			toggle_for_edit_curve(1);
			posRadio->setChecked(true);
		}
	}
	if (myState.eState == EDITCSTATE && state.eState != EDITCSTATE) toggle_for_edit_curve(0);

	if (state.eState == DEFCSTATE && myState.eState != DEFCSTATE){
		defCButton->setText("Quit Connecting");
		//defText->show();
	}
	if (state.eState != DEFCSTATE && myState.eState == DEFCSTATE){
		defCButton->setText("Edit Connection");
		//defText->hide();
	}

	if (state.eState == DELCSTATE && myState.eState != DELCSTATE){
		delCButton->setText("Quit delete");
	}
	if (state.eState != DELCSTATE && myState.eState == DELCSTATE){
		delCButton->setText("Delete Curve");
	}


	if (!editCButton->isEnabled() && state.eState == IDLESTATE){
		editCButton->setEnabled(true);
	}
	if (editCButton->isEnabled() && state.eState == PREPSTATE){
		editCButton->setEnabled(false);
	}
	if (!defCButton->isEnabled() && state.eState == IDLESTATE){
		defCButton->setEnabled(true);
	}
	if (defCButton->isEnabled() && state.eState == PREPSTATE){
		defCButton->setEnabled(false);
	}
	if (!delCButton->isEnabled() && state.eState == IDLESTATE){
		delCButton->setEnabled(true);
	}
	if (delCButton->isEnabled() && state.eState == PREPSTATE){
		delCButton->setEnabled(false);
	}

	if (!simuButton->isEnabled() && state.eState == IDLESTATE){
		simuButton->setEnabled(true);
	}
	if (simuButton->isEnabled() && state.eState == PREPSTATE){
		simuButton->setEnabled(false);
	}

	if (!stbButton->isEnabled() && state.eState == IDLESTATE){
		stbButton->setEnabled(true);
	}
	if (stbButton->isEnabled() && state.eState == PREPSTATE){
		stbButton->setEnabled(false);
	}

	if (!optButton->isEnabled() && state.eState == IDLESTATE){
		optButton->setEnabled(true);
	}
	if (optButton->isEnabled() && state.eState == PREPSTATE){
		optButton->setEnabled(false);
	}

	myState = state;

}