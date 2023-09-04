#include "montecarlo.h"
#include "math.h"
#include <stdlib.h> // define exit() function
#include <time.h>

void PKU_MC::DoOneIon()
{
	int IonFlag;
	// init the orginal paramters for current ions
	InitialIon();
	IonFlag = 1;

	while (IonFlag == 1) // the current ion are under processing
	{
		FreeLength();
		ChoseCollisionAtom();
		EleEnergyLoss();
		EmissionAngle();
		IonNewCondition();
		ProcessRecorder();
		OutOrNot(IonFlag);
	} // end of one ions

	CollisionRecording(); // add the results of this ion to the recording arraies
}

double PKU_MC::ChoseCollisionAtom()
{

	int MaxElement;
	float tempRandValue;

	tempRandValue = PKURandom();
	while (tempRandValue == 0)
	{
		tempRandValue = PKURandom();
	}
	P = PMAX[CurrentLayer] * sqrt(tempRandValue);
	MaxElement = LayerElementNo[CurrentLayer];

	// decided which kinds of atoms to be knocked!

	tempRandValue = PKURandom();

	for (CurrentElement = 1; CurrentElement <= MaxElement; CurrentElement++)
	{
		tempRandValue = tempRandValue - TargetAtomConc[CurrentLayer][CurrentElement];
		if (tempRandValue < 0)
			break;
	}
	if (tempRandValue >= 0)
		CurrentElement = LayerElementNo[CurrentLayer];

	ReducedCollEnergy = FItemp[CurrentLayer][CurrentElement] * E;
	B = P / ScreenLen[CurrentLayer][CurrentElement];

	return TargetAtomNumber[CurrentLayer][CurrentElement];
}

double PKU_MC::EleEnergyLoss()
{
	double Ion_Se;

	IE = int(E / IonEnergyKeV + .5);
	if (IE != 0)
		SEE = SE[CurrentLayer][IE];
	if (E < IonEnergyKeV)
		SEE = SE[CurrentLayer][1] * sqrt(E / IonEnergyKeV);
	DEE = IonFreeLength * SEE;
	Ion_Se = DEE;
	return Ion_Se;
}

double PKU_MC::EmissionAngle()
{
	double AngleValue, R, Q;
	// using the magic method

	if (ReducedCollEnergy <= 10)
	{
		R = B;
		RR = -2.7 * log(ReducedCollEnergy * B);
		if (RR >= B)
		{
			RR = -2.7 * log(ReducedCollEnergy * RR);
			if (RR >= B)
				R = RR;
		}
		do
		{
			EX1 = 0.18175 * exp(-3.1998 * R);
			EX2 = 0.50986 * exp(-0.94229 * R);
			EX3 = 0.28022 * exp(-0.4029 * R);
			EX4 = 0.028171 * exp(-0.20162 * R);
			V = (EX1 + EX2 + EX3 + EX4) / R;
			V1 = -(V + 3.1998 * EX1 + 0.94229 * EX2 + 0.4092 * EX3 + 0.20162 * EX4) / R;

			FR = B * B / R + V * R / ReducedCollEnergy - R;
			FR1 = -B * B / (R * R) + (V + V1 * R) / ReducedCollEnergy - 1;
			Q = FR / FR1;
			R = R - Q;
		} while (fabs(Q / R) > 0.001);

		ROC = -2.0 * (ReducedCollEnergy - V) / V1;
		SQE = sqrt(ReducedCollEnergy);
		CC = (0.011615 + SQE) / (0.0071222 + SQE);
		AA = 2.0 * ReducedCollEnergy * (1.0 + (0.99229 / SQE)) * pow(B, CC);
		FF = (sqrt(AA * AA + 1.0) - AA) * ((9.3066 + ReducedCollEnergy) / (14.813 + ReducedCollEnergy));
		DELTA = (R - B) * AA * FF / (FF + 1.0);
		CO = (B + DELTA + ROC) / (R + ROC);
		C2 = CO * CO;
		S2 = 1.0 - C2;
		if (S2 < -1)
			printf("\n%d  %lf\n", ElapsedIonNo, B);
		CT = 2.0 * C2 - 1.0;
		ST = sqrt(1.0 - CT * CT);
	}
	// using the RusefScattering methods
	else
	{
		S2 = 1.0 / (1.0 + (1.0 + B * (1.0 + B)) * (2.0 * ReducedCollEnergy * B) * (2.0 * ReducedCollEnergy * B));
		if (S2 < -1)
			printf("\n%d  %lf\n", ElapsedIonNo, B);
		C2 = 1.0 - S2;
		CT = 2.0 * C2 - 1.0;
		ST = sqrt(1.0 - CT * CT);
	}
	AngleValue = acos(CT);
	return AngleValue;
}

double PKU_MC::FreeLength()
{
	float tempvalue;

	// when the ion energy is high
	IonCollisionNo++;
	ReducedCollEnergy = E * F[CurrentLayer];
	EEG = sqrt(ReducedCollEnergy * EPSDG[CurrentLayer]);
	PMAX[CurrentLayer] = A[CurrentLayer] / (EEG + sqrt(EEG) + 0.125 * pow(EEG, 0.1));
	IonFreeLength = 1.0 / (PI * PMAX[CurrentLayer] * PMAX[CurrentLayer] * ARHO[CurrentLayer]);

	// when the ion enrgy is low

	tempvalue = PKURandom();

	if (IonCollisionNo == 1)
		IonFreeLength = tempvalue * AMIN1(IonFreeLength, channelwidth);

	return IonFreeLength;
}

double PKU_MC::IonNewCondition()
{
	// int I,J;
	double MAX(0), X1; // MAX is not OK!!!!!!!!!!!!!!!!!!!!!!!!!
	float tempvalue;

	DEN = recoilfactor[CurrentLayer][CurrentElement] * S2 * E;
	E = E - DEN - DEE;
	if (DEE > MAX)
		MAX = DEE;
	IonWayLength = IonWayLength + IonFreeLength - DistanceNuclear;

	X = X + (IonFreeLength - DistanceNuclear) * COSIN;
	Y = Y + (IonFreeLength - DistanceNuclear) * COSY;

	I = AMIN1(fabs(X / channelwidth) + 1.0, 100.0);
	J = AMIN1(fabs(Y / channelwidth) + 1.0, 50.0);

	// if the porous materials are used, then the new position should consider the bubbles
	// the method is also simple: if yes, then let the ion fly to the new edge of solid materisl
	// be careful, this is only a 2D problem!

	/*

		BubbleEffect();  //X,Y is gloabal!

	*/

	if ((ALPHA != 0) && (Y <= 0))
		J = 1;

	tempvalue = PKURandom();
	PHI = 2.0 * PI * tempvalue;

	PSI = atan(ST / (CT + M1_to_M2[CurrentLayer][CurrentElement]));
	if (PSI < 0)
		PSI = PSI + PI;
	X1 = -COSIN * COSY / (SINE * SINY + pow(10, -8));
	if (fabs(X1) > 1.0)
		X1 = X1 / fabs(X1);
	DELTA = PHI - acos(X1);

	COSIN = COSIN * cos(PSI) + SINE * sin(PSI) * cos(PHI);
	COSY = COSY * cos(PSI) + SINY * sin(PSI) * cos(DELTA);

	SINY = sqrt(1.0 - COSY * COSY);
	SINE = sqrt(1.0 - COSIN * COSIN);
	return 1;
}

void PKU_MC::ProcessRecorder()
{

	double EPSD, EN;

	// total energy
	MTOT[I - 1][J] = MTOT[I - 1][J] + DEN + DEE;

	// the nuclear energy loss treamting process
	if (DEN < DisEnergy)
		PhoneEDis[I] = PhoneEDis[I] + DEN;
	else
	{
		EPSD = DamageEFract[CurrentLayer] * DEN;

		// EN ����ȱ�ݲ���������Modified Kinchin-Pease  model
		EN = DEN / (1.0 + ElectronEFract[CurrentLayer] * (EPSD + 0.4 * pow(EPSD, 0.75) + 3.4 * pow(EPSD, (1.0 / 6.0))));

		if (EN < DisEnergy)
			PhoneEDis[I] = PhoneEDis[I] + DEN;
		else
		{

			MVAC[I - 1][J] = MVAC[I - 1][J] + 1;
			IVAC[I] = IVAC[I] + 1;
			if (EN > 0) //!!!!
				RPHON[I] = RPHON[I] + EN - DisEnergy;

			// Multi-defect production!
			if (EN >= 2.5 * DisEnergy)
			{
				MVAC[I - 1][J] = MVAC[I - 1][J] - 1.0 + 0.4 * EN / DisEnergy;
				RPHON[I] = RPHON[I] + DisEnergy - 0.4 * EN;
				RVAC[I] = RVAC[I] - 1.0 + 0.4 * EN / DisEnergy;
			}

			MION[I - 1][J] = MION[I - 1][J] + DEN - EN;
			RION[I] = RION[I] + DEN - EN;
		}
	}

	// eletric energy loss treating process
	IonizatinEDis[I] = IonizatinEDis[I] + DEE;

	MION[I - 1][J] = MION[I - 1][J] + DEE;
}

// the ion is out of the sample?
// transmitted or RBS
int PKU_MC::OutOrNot(int &Ionflag)
{
	if (X < 0)
	{
		Ionflag = 0;
		return 1;
	}
	else
	{
		if (X <= LayertoSurface[1])
			CurrentLayer = 1;
		else if (X <= LayertoSurface[2])
			CurrentLayer = 2;
		else
			CurrentLayer = 3;
		if (X >= LayertoSurface[CurrentLayer])
		{
			Ionflag = 0;
			return 1;
		}
	}
	if (E <= StoppingEnergy)
	{
		Ionflag = 0;
		return 1;
	}
}

// this will add the results of ions at different collision to the recording arraies
void PKU_MC::CollisionRecording()
{
	int j, IP;
	AtomType LocalAtom;

	if (X < 0)
	{
		RBSIonNo = RBSIonNo + 1; // �ص���
		RBSIonE = RBSIonE + E;

		//***************    recording the transmitted ions  **************
		LocalAtom.flag = RBSION;
		LocalAtom.x = X;
		LocalAtom.y = Y;
		LocalAtom.Z = IonZ;
		LocalAtom.mass = IonM;
		LocalAtom.alpha = COSIN; // here only recording the cos of alpha
		LocalAtom.energy = E;

		RecordingAtom(LocalAtom);
	}
	else if (X > LayertoSurface[MaxLayer])
	{
		TransIonNo = TransIonNo + 1; // ���
		TransIonE = TransIonE + E;

		//***************    recording the transmitted ions  **************
		LocalAtom.flag = TRANSMITTEDION;
		LocalAtom.x = X;
		LocalAtom.y = Y;
		LocalAtom.Z = IonZ;
		LocalAtom.mass = IonM;
		LocalAtom.alpha = COSIN; // here only recording the cos of alpha
		LocalAtom.energy = E;
		RecordingAtom(LocalAtom);
	}
	else
	{
		//***************    recording the normal ions  **************
		LocalAtom.flag = NORMALION;
		LocalAtom.x = X;
		LocalAtom.y = Y;
		LocalAtom.Z = IonZ;
		LocalAtom.mass = IonM;
		LocalAtom.alpha = COSIN; // here only recording the cos of alpha
		LocalAtom.energy = E;
		RecordingAtom(LocalAtom);

		// insity anayliss
		ionrange[I] = ionrange[I] + 1;
		MPART[I - 1][J] = MPART[I - 1][J] + 1;

		IP = int(IonWayLength / channelwidth + 1.0);
		if (IP > 100)
			IP = 100;
		iondistance[IP] = iondistance[IP] + 1;
		XSUM = XSUM + X;
		X2SUM = X2SUM + pow(X, 2);
		X3SUM = X3SUM + pow(X, 3);
		X4SUM = X4SUM + pow(X, 4);
		Y2SUM = Y2SUM + Y * Y;
		XY2SUM = XY2SUM + X * Y * Y;
		X2Y2SUM = X2Y2SUM + X * X * Y * Y;
		Y4SUM = Y4SUM + pow(Y, 4);
		PLSUM = PLSUM + IonWayLength;
		PL2SUM = PL2SUM + IonWayLength * IonWayLength;
		ICSUM = ICSUM + IonCollisionNo;
	}
}