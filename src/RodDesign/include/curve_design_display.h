#ifndef CURVE_DESIGN_DISPLAY_H_
#define CURVE_DESIGN_DISPLAY_H_

#include <QGLWidget>
#include "Setting.h"

struct curve_node;
class loop_subdivision_mesh;
class EditTools;
class AABB_struct;
class curve_on_mesh;
class CurveNetwork;
class MomentReductionMap;

class curve_design_display :public QGLWidget{
	Q_OBJECT

public:
	curve_design_display(QWidget *parent = 0);
	~curve_design_display();
	void reset_param_for_mesh();
	void state_changed(struct state_param state);
	void addCurveToNet();

	loop_subdivision_mesh *imesh;
	AABB_struct *iAABB;
	class CurveNetwork *cNet;
	class MomentReductionMap *mrMap;
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

private:
	void draw();
	void display_mesh();
	void display_curve_ends();
	void display_curve_network();
	void display_connection();
	void display_moment_reduction_map();
	void displayEditState();
	void search_clicked_point(int cx,int cy, struct curve_node *cnode);
	void search_selected_element(int cx, int cy);
	void update_cutting_plane(int dx, int dy);
	void update_edge_direction(int cx,int cy,int dx, int dy);
	bool cross_section_curve_search(class curve_on_mesh *curve);

	GLfloat shiftX, shiftY, shiftZ;
	GLfloat Gpos[3];
	GLfloat Zoom,Zdefault;
	GLfloat rw, rx, ry, rz;
	QPoint lastPos;
	int windowWid, windowHei;
	float windowRatio;
	bool ADDCSEARCH,CONNECTSEARCH;
	struct state_param myState;
	
	struct curve_node *Cends;
	class curve_on_mesh *Curv;

	int selectC, selectN,selectE;
	struct curve_node *moveP;
	float cutNorm[3],eNorm[3];
	int selem1[2], selem2[2],sconn;
};

#endif
