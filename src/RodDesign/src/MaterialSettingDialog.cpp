#include<qlayout.h>
#include<qpushbutton.h>
#include<qlineedit.h>
#include<qvalidator.h>
#include<qlabel.h>
#include<qmessagebox.h>
#include<qgroupbox.h>
#include "MaterialSettingDialog.h"
#include "Setting.h"

MaterialSettingDialog::MaterialSettingDialog(struct SimulationParameter *param, QWidget *parent) :QDialog(parent){

	sParam = param;

	//パラメータ設定
	densText = new QLabel(tr("Density (g/mm^3)"));
	densEdit = new QLineEdit(this);
	QDoubleValidator *densValid = new QDoubleValidator(this);
	densEdit->setValidator(densValid);
	QString dnum = QString::number(sParam->Density);
	densEdit->setText(dnum);

	widText = new QLabel(tr("Width (mm)"));
	widEdit = new QLineEdit(this);
	QDoubleValidator *widValid = new QDoubleValidator(this);
	widEdit->setValidator(widValid);
	QString wnum = QString::number(sParam->Wid);
	widEdit->setText(wnum);
	
	heiText = new QLabel(tr("Height (mm)"));
	heiEdit = new QLineEdit(this);
	QDoubleValidator *heiValid = new QDoubleValidator(this);
	heiEdit->setValidator(heiValid);
	QString hnum = QString::number(sParam->Hei);
	heiEdit->setText(hnum);

	stifText = new QLabel(tr("Stiffness"));
	streText = new QLabel(tr("Stretch (kg/mm/s^2)"));
	streEdit = new QLineEdit(this);
	QDoubleValidator *streValid = new QDoubleValidator(this);
	streEdit->setValidator(streValid);
	QString stnum = QString::number(sParam->Es);
	streEdit->setText(stnum);

	twisText = new QLabel(tr("Twist (kg/mm/s^2)"));
	twisEdit = new QLineEdit(this);
	QDoubleValidator *twisValid = new QDoubleValidator(this);
	twisEdit->setValidator(twisValid);
	QString twnum = QString::number(sParam->Gt);
	twisEdit->setText(twnum);

	bendText = new QLabel(tr("Bending (kg/mm/s^2)"));
	bendText1 = new QLabel(tr("B11"));
	bendEdit1 = new QLineEdit(this);
	QDoubleValidator *benValid1 = new QDoubleValidator(this);
	bendEdit1->setValidator(benValid1);
	QString bnum1 = QString::number(sParam->Eb11);
	bendEdit1->setText(bnum1);

	bendText2 = new QLabel(tr("B22"));
	bendEdit2 = new QLineEdit(this);
	QDoubleValidator *benValid2 = new QDoubleValidator(this);
	bendEdit2->setValidator(benValid2);
	QString bnum2 = QString::number(sParam->Eb22);
	bendEdit2->setText(bnum1);

	bendText3 = new QLabel(tr("B12"));
	bendEdit3 = new QLineEdit(this);
	QDoubleValidator *benValid3 = new QDoubleValidator(this);
	bendEdit3->setValidator(benValid3);
	QString bnum3 = QString::number(sParam->Eb12);
	bendEdit3->setText(bnum3);

	stifBox = new QGroupBox(this);

	saveButton = new QPushButton(tr("Save Setting"));
	connect(saveButton, SIGNAL(clicked()), this, SLOT(saveSetting()));


	//レイアウト
	QHBoxLayout *densLayout = new QHBoxLayout;
	densLayout->addWidget(densText);
	densLayout->addWidget(densEdit);

	QHBoxLayout *crossLayout = new QHBoxLayout;
	crossLayout->addWidget(widText);
	crossLayout->addWidget(widEdit);
	crossLayout->addWidget(heiText);
	crossLayout->addWidget(heiEdit);

	QVBoxLayout *stifLayout = new QVBoxLayout;
	
	QHBoxLayout *strLayout = new QHBoxLayout;
	strLayout->addWidget(streText);
	strLayout->addWidget(streEdit);
	strLayout->addWidget(twisText);
	strLayout->addWidget(twisEdit);
	stifLayout->insertLayout(1, strLayout);

	stifLayout->addWidget(bendText);

	QHBoxLayout *bendLayout = new QHBoxLayout;
	bendLayout->addWidget(bendText1);
	bendLayout->addWidget(bendEdit1);
	bendLayout->addWidget(bendText2);
	bendLayout->addWidget(bendEdit2);
	bendLayout->addWidget(bendText3);
	bendLayout->addWidget(bendEdit3);
	stifLayout->insertLayout(3, bendLayout);

	stifBox->setLayout(stifLayout);

	QVBoxLayout *mainLayout = new QVBoxLayout;
	mainLayout->insertLayout(2, densLayout);
	mainLayout->insertLayout(3, crossLayout);
	mainLayout->addWidget(stifBox);
	mainLayout->addWidget(saveButton);
	setLayout(mainLayout);

	setWindowTitle(tr("Material Setting"));

}


void MaterialSettingDialog::saveSetting(){


	//代入
	if (!densEdit->text().isEmpty()){
		sParam->Density = densEdit->text().toDouble();
	}
	if (!widEdit->text().isEmpty()){
		sParam->Wid = widEdit->text().toDouble();
	}
	if (!heiEdit->text().isEmpty()){
		sParam->Hei = heiEdit->text().toDouble();
	}

	if (!streEdit->text().isEmpty()){
		sParam->Es = streEdit->text().toDouble();
	}
	if (!twisEdit->text().isEmpty()){
		sParam->Gt = twisEdit->text().toDouble();
	}
	if (!bendEdit1->text().isEmpty()){
		sParam->Eb11 = bendEdit1->text().toDouble();
	}
	if (!bendEdit2->text().isEmpty()){
		sParam->Eb22 = bendEdit2->text().toDouble();
	}
	if (!bendEdit3->text().isEmpty()){
		sParam->Eb12 = bendEdit3->text().toDouble();
	}
	

	printf("Dens:%e Wid:%e Hei:%e\n", sParam->Density, sParam->Wid, sParam->Hei);

	printf("Str:%e Twis:%e\n", sParam->Es, sParam->Gt);
	printf("Bend:%e %e %e\n", sParam->Eb11, sParam->Eb22, sParam->Eb12);

	accept();
}