// in this class, define all the data, variables, araarys used in the simualtion
// xue 2007
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "malloc.h"

#define MAX_ARRAY 101 // the the distance can be divided into..
#define MAX_ELEMENT 8 // the element number in the sampels
#define MAX_LAYER 4	  // the layers (compound are not same)number of the sample
#define PI 3.1415926

// The follwoing definisiton is for the file recording
#define TRANSMITTEDION 0
#define RBSION 1
#define NORMALION 2	   // stayed inside the sample
#define RECOILEDATOM 3 // this is speicall for the recoild target atoms

#define NORMALTYPE double
#define INTYPE float

#define MONIORING1 // displaying simulation process information or not

// this structure is special for the atom/ion information recodrding
// such as the transmitted ions

typedef struct
{
	int Z,	  // atomic no
		flag; // transmitted? RBS  or ......

	float mass,
		alpha, // direction
		x,	   // location
		y,
		energy;
} AtomType;

class Data
{

private:								   // not used yet
	NORMALTYPE IO[MAX_LAYER][MAX_ELEMENT], //
		K[MAX_LAYER][MAX_ELEMENT],
		KL[MAX_LAYER][MAX_ELEMENT],
		EPSBK[MAX_ELEMENT],
		SBK[MAX_ELEMENT],
		IOUT[101],
		XOUT[101],
		C[MAX_ELEMENT],
		LM[MAX_ELEMENT], // distance between atoms
		LF[MAX_ELEMENT];

public:
	//**************     ion paramters     *************************************

	long ionnumber;
	NORMALTYPE IonEnergyKeV, StoppingEnergy, IonM, IonZ, ALPHA, DA, IZ;

	// cosy=cos(alpha), where alpha is the angle to the surface noramal direction
	// PHI:��λ��
	// PSI��ɢ���
	NORMALTYPE COSIN, SINY, SINE, COSY, PHI, PSI, DELTA, ZETA, DistanceNuclear;
	NORMALTYPE VFERMI,
		LFCTR,
		XLAMB,
		PCOEF[9],
		DisEnergy,
		channelwidth,
		IonEnergyEV,
		ReducedCollEnergy; // instead of ED ReducedCollEnergy: originally epsilong

	// *******************       samples        ********************************

	int MaxLayer, CurrentLayer, num_layers, LayerElementNo[MAX_LAYER];

	NORMALTYPE TargetAtomConc[MAX_LAYER][MAX_ELEMENT], // ���ָ�Ԫ���ڸ����еı���
		TargetAtomMass[MAX_LAYER][MAX_ELEMENT],		   // ��ԭ������
		TargetAtomNumber[MAX_LAYER][MAX_ELEMENT],	   // ����ԭ������
		M1_to_M2[MAX_LAYER][MAX_ELEMENT],			   // M1/M2
		ScreenLen[MAX_LAYER][MAX_ELEMENT],			   // ���γ��ȣ�screen lengths used to estimate sinlge scattering
		FItemp[MAX_LAYER][MAX_ELEMENT],				   // epsilong/energy
		recoilfactor[MAX_LAYER][MAX_ELEMENT],		   // recoil factor

		layerdensity[MAX_LAYER], // layer density
		VF[MAX_LAYER][MAX_ELEMENT],
		MU[MAX_ELEMENT],	// M2/M1
		IONIZ[MAX_ELEMENT], // mean of IO
		H[MAX_ELEMENT],

		LayerThickness[MAX_LAYER], // layer thickness
		LayertoSurface[MAX_LAYER], // depth of the bottom of every layer
		M2[MAX_ELEMENT],
		Z_2[MAX_ELEMENT],
		A[MAX_ELEMENT], // average of AI
		F[MAX_ELEMENT], //

		ARHO[MAX_ELEMENT], // atoms/Ang3
		PMAX[MAX_ELEMENT],
		DamageEFract[MAX_ELEMENT],	 // fract of reconil energy into damage
		ElectronEFract[MAX_ELEMENT]; // fraction of recoil energy into electrctionic loss

	NORMALTYPE SE[4][1001], // se list in 3 layers
		SEO[1001],
		EPSDG[4];

	//******************       recording arrays    ********************************

	// in the following definition, I_ means energy deposted by incident ins;
	// R_ means by the recoiled atoms
	// 1D
	NORMALTYPE ionrange[MAX_ARRAY], // projected in range
		iondistance[MAX_ARRAY],		// ion travel distance before it stops
		IonizatinEDis[MAX_ARRAY],	// depth where energy loss to target ionization occures
		RPHON[MAX_ARRAY],			// depth where recoil atoms create Phonons
		PhoneEDis[MAX_ARRAY],		// depth where energy loss to target phonones occures

		RION[MAX_ARRAY],
		TotalIonization[MAX_ARRAY],
		RVAC[MAX_ARRAY],
		IVAC[MAX_ARRAY];
	// 2D
	NORMALTYPE M[102][31], // ͸�����ӵ������ͽǶ�

		MTOT[101][51],	// �������ֲ�����ά����
		MION[101][51],	// ���ӵ��������ֲ�
		MPART[101][51], // �������ӷֲ�
		MVAC[101][51];	// ��Ѩ�ֲ�

	NORMALTYPE XSUM,
		X2SUM,
		X3SUM,
		X4SUM,
		PLSUN,
		PL2SUM,
		Y4SUM,
		X2YSUM,
		Y2SUMM,
		ICSUM,
		Y2SUM,
		XY2SUM,
		X2Y2SUM,
		PLSUM;

	long RBSIonNo,
		TransIonNo,
		ElapsedIonNo,
		IonCollisionNo; // ������ײ������¼

	//***************      here is some functions  **************************
public:
	NORMALTYPE AMIN1(NORMALTYPE, NORMALTYPE);
	NORMALTYPE AMAX1(NORMALTYPE, NORMALTYPE);
	NORMALTYPE PKURandom();
	void TestRand(); // this will test the random values for 1M and to test its distribution
};
