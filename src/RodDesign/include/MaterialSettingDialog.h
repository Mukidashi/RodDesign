#ifndef _MATERIALSETTINGDIALOG_H_
#define _MATERIALSETTINGDIALOG_H_

#include<qdialog.h>

class QPushButton;
class QLineEdit;
class QLabel;
class QGroupBox;
struct SimulationParameter;

class MaterialSettingDialog : public QDialog
{
	Q_OBJECT
public:
	MaterialSettingDialog(struct SimulationParameter *sParam, QWidget *parent = 0);

	private slots:
	void saveSetting();

private:
	QPushButton *saveButton;
	QLineEdit *densEdit, *widEdit, *heiEdit;
	QLineEdit *streEdit, *twisEdit, *bendEdit1, *bendEdit2, *bendEdit3;
	QLabel *densText, *widText, *heiText;
	QLabel *stifText, *streText, *twisText, *bendText,*bendText1, *bendText2, *bendText3;
	QGroupBox *stifBox;
	struct SimulationParameter *sParam;
};



#endif