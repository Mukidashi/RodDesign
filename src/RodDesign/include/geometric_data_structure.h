#ifndef _GEOMETRIC_DATA_STRUCTURE_H_
#define _GEOMETRIC_DATA_STRUCTURE_H_

#include <Eigen/Core>

//queue‚Ì—v‘f
struct qelem{
	struct qelem *next,*back;
	int ln;
};


class Mesh {
protected:
public:
	Mesh();
	~Mesh();
	bool CONSTRUCT_FLG;
	int nodeN, edgeN, faceN;
	float(*Nodes)[3];
	int(*Faces)[3], (*Edges)[2];
	int **PFlist, *PFsize;
	int **PElist, *PEsize;
	int(*FElist)[3], (*EFlist)[2];
	void construct_data_structure(void);
	void construct_data_structure_without_edge(void);
	void offset_surface_mesh(float offv);
	void mesh_copy(class Mesh *cmesh);
};




void read_obj_file(std::string fname, class Mesh *mesh);

void read_libigl_mesh_format(Eigen::MatrixXf V, Eigen::MatrixXi F, class Mesh *mesh);

void output_obj_file(std::string fname, class Mesh *mesh);

int labeling_by_dihedral_angle(class Mesh *mesh, int *label, float ang_thres);

void output_labeled_mesh_by_obj(std::string fname, class Mesh *mesh, int *label, int ln);


#endif