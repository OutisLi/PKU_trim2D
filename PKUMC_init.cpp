#include "montecarlo.h"
#include "math.h"
#include <stdlib.h> // define exit() function
#include <time.h>

// the following will init the ion information at the begin
// all the message needed to be input from the files have been aready!
void PKU_MC::InitialIon()
{

	E = IonEnergyEV;	// Energy
	COSIN = cos(ALPHA); // direction to the normal
	SINY = COSIN;

	SINE = sin(ALPHA);
	COSY = SINE;

	CurrentLayer = 1;
	IonWayLength = 0;
	IonCollisionNo = 0;
	X = 0;
	Y = 0;
}

void PKU_MC::CurrentCondition(int ioncounter)
{

	if (ElapsedIonNo % 400 == 0)
	{
		printf("Elapsed ion number is: %5d\n", ElapsedIonNo);
	}
}
// all the init are finished here!
void PKU_MC::SimulationInit()
{
	// init the talbles and the parameters,set all their values to 0 if no special data are asigend
	InitialPara();

	// inputting parameters
	init();
	// Se list used in the simulatin for electronic energy deposition calculatin
	InitialSE();
}

// all the init are finished here!
void PKU_MC::InitialPara()
{
	int LayerNo, ElementNo, i, j;
	XSUM = 0;
	X2SUM = 0;
	X3SUM = 0;
	X4SUM = 0;
	Y2SUM = 0;
	XY2SUM = 0;
	X2Y2SUM = 0;
	Y4SUM = 0;
	PLSUM = 0;
	PL2SUM = 0;
	ICSUM = 0;
	RBSIonNo = 0;
	TransIonNo = 0;
	RBSIonE = 0;
	TransIonE = 0;
	for (LayerNo = 0; LayerNo < MAX_LAYER; LayerNo++)
	{
		layerdensity[LayerNo] = 0;
		LayerElementNo[LayerNo] = 0;
		LayerThickness[LayerNo] = 0;
		LayertoSurface[LayerNo] = 0;
		for (ElementNo = 0; ElementNo < MAX_ELEMENT; ElementNo++)
		{
			TargetAtomConc[LayerNo][ElementNo] = 0;	  // ���ָ�Ԫ���ڸ����еı���
			TargetAtomMass[LayerNo][ElementNo] = 0;	  // ��ԭ������
			TargetAtomNumber[LayerNo][ElementNo] = 0; // ����ԭ������
			M1_to_M2[LayerNo][ElementNo] = 0;		  // M1/M2
			ScreenLen[LayerNo][ElementNo] = 0;		  // ���γ��ȣ�screen lengths used to estimate sinalg escattering
			FItemp[LayerNo][ElementNo] = 0;			  // epsilong/energy
			recoilfactor[LayerNo][ElementNo] = 0;	  // recoil factor
			//  IO[LayerNo][ElementNo]=0; //

			//  K[LayerNo][ElementNo]=0;
			//  KL[LayerNo][ElementNo]=0;
			VF[LayerNo][ElementNo] = 0;
		}
	}
	for (ElementNo = 0; ElementNo < MAX_ELEMENT; ElementNo++)
	{
		MU[ElementNo] = 0;
		IONIZ[ElementNo] = 0; // mean of IO
		H[ElementNo] = 0;
		M2[ElementNo] = 0;
		Z_2[ElementNo] = 0;
		ARHO[ElementNo] = 0; // atoms/Ang3
		F[ElementNo] = 0;	 //
		PMAX[ElementNo] = 0;
		DamageEFract[ElementNo] = 0;   // fract of reconil energy int damage
		ElectronEFract[ElementNo] = 0; // fraction of recoil energy into electrctionic loss
	}
	for (i = 0; i < MAX_ARRAY; i++)
	{
		ionrange[i] = 0;	  // projected in range
		iondistance[i] = 0;	  // ion travel distance before it stops
		IonizatinEDis[i] = 0; // depth where energy loss to target ionization occures
		RPHON[i] = 0;		  // depth where recoil atoms create Phonons
		PhoneEDis[i] = 0;	  // depth where energy loss to target phonones occures
	}

	for (i = 0; i < 101; i++)
	{
		RION[i] = 0;
		TotalIonization[i] = 0;
		RVAC[i] = 0;
		IVAC[i] = 0;

		for (j = 0; j < 51; j++)
		{
			MTOT[i][j] = 0;	 // �������ֲ�����ά����
			MION[i][j] = 0;	 // ���ӵ��������ֲ�
			MPART[i][j] = 0; // �������ӷֲ�
			MVAC[i][j] = 0;	 // ��Ѩ�ֲ�
		}
	}

	for (i = 0; i < 4; i++)
	{
		A[i] = 0;
		EPSDG[i] = 0;

		for (j = 0; j < 1001; j++)
		{
			SE[i][j] = 0;
		}
	}

	for (i = 0; i < 1001; i++)
	{
		SEO[i] = 0;
	}

	for (i = 0; i < 9; i++)
	{
		PCOEF[i] = 0;
	}

	for (i = 0; i < 102; i++)
	{
		for (j = 0; j < 31; j++)
		{
			M[i][j] = 0; // ͸�����ӵ������ͽǶ�
		}
	}
}