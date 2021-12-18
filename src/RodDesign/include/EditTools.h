#ifndef EDITTOOLS_H_
#define EDITTOOLS_H_

#include<QFrame>
#include "Setting.h"

class QPushButton;
class QLabel;
class QGroupBox;
class QRadioButton;

class EditTools :public QFrame{
	Q_OBJECT
public:
	EditTools(QWidget *parent = 0);
	void toggle_for_add_curve(int sNum); 
	void toggle_for_edit_curve(int sNum);
	void state_changed(struct state_param state);
signals:
	void editStateChanged(int estate);
private slots:
	void toggled_edit_radiobox();
private:

	QPushButton *addCButton, *delCButton, *editCButton, *defCButton, *simuButton,*optButton,*stbButton;
	QPushButton *spButton, *sdButton;
	QPushButton *agrButton1, *agrButton2;
	QLabel *addText1,*addText2,*defText; 
	QGroupBox *addCBox,*editCBox,*defCBox;
	QRadioButton *posRadio, *direRadio,*delRadio;
	struct state_param myState;
};

#endif