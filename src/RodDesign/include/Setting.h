#ifndef SETTING_H_
#define SETTING_H_

enum EntireState{
	PREPSTATE = 0,
	IDLESTATE = 1,
	ADDCSTATE = 2,
	EDITCSTATE = 3,
	DEFCSTATE = 4,
	DELCSTATE = 5,
};

enum AddCState{
	POINTZERO = 0,
	POINTONE = 1,
	POINTTWO = 2,
	PLANESEARCH = 3,
	PLANESETTED = 4,
};

enum EditCState{
	POINTMOVE = 0,
	EDGEDIRE = 1,
	DELETECONS = 2,
};

enum SimuState{
	SIDLESTATE = 0,
	SIMUSTATE = 1,
};


struct state_param{
	bool MESH_EXIST;
	EntireState eState;
	AddCState acState;
	EditCState ecState;
	SimuState sState;
};

//çﬁóøíËêî
//extern const float rho;
//extern const float Young;
//extern const float Poisson;
//extern const float Gshear;
//extern const float rho2;
//
//extern const double ROD_WID;
//extern const double ROD_HEI;


struct SimulationParameter{
	double Density;
	double Es, Gt, Eb11, Eb12, Eb22;
	double Wid, Hei;
	double sampL[2];
};
extern const float rho2;

#endif