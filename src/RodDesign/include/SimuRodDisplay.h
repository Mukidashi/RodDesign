#ifndef _SIMURODDISPLAY_H_
#define _SIMURODDISPLAY_H_


#include <QGLWidget>
#include "Setting.h"

class RodNetwork;
class Mesh;

class SimuRodDisplay :public QGLWidget{
	Q_OBJECT

public:
	SimuRodDisplay(QWidget *parent = 0);
	void reset_param_for_mesh();
	void state_changed(struct state_param state);

	class RodNetwork *rNet;
	Mesh *imesh;
	Qt::Key Pressed_Key;

protected:
	QSize sizeHint() const;
	void initializeGL();
	void resizeGL(int width, int height);
	void paintGL();
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);

signals:
	void sync_state_call(struct state_param state);
private slots :
	void simuIdleFunc();
private:
	void draw();
	void display_mesh();
	void display_rod_network();

	GLfloat shiftX, shiftY, shiftZ;
	GLfloat Gpos[3];
	GLfloat Zoom, Zdefault;
	GLfloat rw, rx, ry, rz;
	QPoint lastPos;
	QTimer *simuTimer;
	int FRAMECOUNT, SHORTCOUNT;
	double TIMESTEP, ANITIME,CALTIME;

	struct state_param myState;
};

#endif
