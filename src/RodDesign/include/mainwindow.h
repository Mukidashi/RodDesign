#ifndef MAINWINDOW_H_
#define MAINWINDOW_H_

#include<QMainWindow>
#include "Setting.h"

class QAction;
class curve_design_display;
class SimuRodDisplay;
class Mesh;
class loop_subdivision_mesh;
class EditTools;
class QSplitter;
class AABB_struct;
class CurveNetwork;
class RodNetwork;
class MomentReductionMap;

class MainWindow : public QMainWindow{
	Q_OBJECT
public:
	MainWindow(QWidget *parent = Q_NULLPTR);
	~MainWindow();

	void sync_state(struct state_param state);
protected:
	void closeEvent(QCloseEvent *event){};

private slots:
	void open();
	void addCurve();
	void setCurveEnds();
	void setCuttingPlane();
	void sync_state_called(struct state_param state);
	void editCurve();
	void editConnect();
	void deleteCurve();
	void expCurveNet();
	void expRodNet();
	void loadCurveNet();
	void simulateRodNet();
	void stabilize();
	void optimizeRodParam();
	void setMaterialParam();
	void editStateChanged(int estate);
private:
	void createActions();
	void createMenus();
	bool loadFile(const QString &filename);
	bool eventFilter(QObject *target, QEvent *event);

	QMenu *fileMenu,*expSubMenu,*loadSubMenu,*setMenu;
	QAction *openAct, *expCurvAct, *loadCurvAct, *expRodAct, *materialAct;
	QSplitter *split,*split2;
	curve_design_display *Cedit;
	SimuRodDisplay *sDisp;

	EditTools *Etools;
	struct state_param myState;

	loop_subdivision_mesh *imesh;
	AABB_struct *iAABB;
	CurveNetwork *cNet;
	RodNetwork *rNet;
	MomentReductionMap *mrMap;
	SimulationParameter sParam;
};

#endif