#ifndef _AABB_DATA_STRUCTURE_H_
#define _AABB_DATA_STRUCTURE_H_

#include "geometric_data_structure.h"

//AABBのノード
struct aabb_node{
	struct aabb_node *snode[8];
	float aabb[3][2];
	int leaf;
};


//basic AABB のノード
struct basic_aabb_node{
	struct basic_aabb_node *dnode[8];
	float aabb[3][2];
	float sumip[8][3];
	int *fmemb;
	int fmemN,label;
	bool LEAF_FLG, is_full;
};


//basic AABBのノードのリスト要素
struct qelem_basic_aabb{
	struct basic_aabb_node *mynode;
	struct qelem_basic_aabb *next;
};


//面を１ユニットとしたAABB
class AABB_struct{
	bool STRUCT_EXIST;
	void free_aabb_node(struct aabb_node *cnode);
	void recursive_range_specification(struct aabb_node *cnode, float rbox[3][2]);
	void recursive_mesh_divide(struct aabb_node *cnode, int dept, int wnum, int(*dwlist)[2]);
	struct qelem* recursive_intersect_triangle_search(struct aabb_node *cnode, float cbox[3][2], struct qelem *ctop);
	struct qelem* recursive_search_intersect_triangle_with_line(struct aabb_node *cnode, float orig[3], float dire[3], struct qelem *ctop);
public:
	AABB_struct();
	~AABB_struct();
	Mesh *mymesh;
	struct aabb_node *root, *znode;
	void build_AABB_struct(Mesh *mesh);
	void update_AABB_range();
	bool closest_point_search(float pc[3],float fps[3][3],float cord[3]);
	bool closest_point_search_with_mesh_info(float pc[3], int *fn, float cord[3]);
	bool search_intersect_triangle_with_line(float orig[3], float dire[3], std::vector<int> *intF);
};


//方向付きのベーシックなAABB
class oriented_basic_AABB_struct{
	bool STRUCT_EXIST;
	void recursive_space_division(struct basic_aabb_node *cnode,int max_depth,int depth);
	void recursive_assign_elements(struct basic_aabb_node *cnode, struct qelem *celem);
	void free_basic_aabb_node(struct basic_aabb_node *cnode);
	struct qelem_basic_aabb* recursive_search_intersect_node_with_line(struct basic_aabb_node *cnode,float orig[3],float dire[3],struct qelem_basic_aabb *caabb);
	struct qelem_basic_aabb* recursive_labeling_of_node(struct basic_aabb_node *cnode, struct qelem_basic_aabb *caabb);
public:
	oriented_basic_AABB_struct();
	~oriented_basic_AABB_struct();
	int MAX_DEPTH;
	float FULL_SIZE;
	float Axis[3][3],Origin[3];
	basic_aabb_node *znode, *root;
	Mesh *mymesh;
	void build_AABB_struct(float size, float orient[3][3],float orig[3],int depth);
	void insert_mesh_elemets(class Mesh *mesh);
	void fill_intermediate_node_along_line(float orig[3], float dire[3]);
};




void free_qelem(struct qelem *bottom, struct qelem *ctop);


bool intersection_check_of_aabb_and_triagle(float aabb[3][2], float fps[3][3]);

#endif