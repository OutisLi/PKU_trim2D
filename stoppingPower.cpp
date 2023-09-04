#include "stoppingpower.h"
#include "math.h"

void EnergyLoss::Se(NORMALTYPE *RESULT)
{
	if (Z1 - 2 < 0)
		H_Se(RESULT);
	else if (Z1 - 2 == 0)
		He_Se(RESULT);
	else
		Heavyion_Se(RESULT);
}

void EnergyLoss::SetParameters(int az1, NORMALTYPE am1, int az2, NORMALTYPE aee)
{
	NORMALTYPE temp1, temp2, temp3, temp4, temp5, uuu;

	Z1 = az1;
	M1 = am1;
	Z2 = az2;
	EE = aee;

	SCOEF(Z1, &temp1, &MM1, &temp2, &temp3, &temp4, &temp5, &LFCTR, &uuu);

	SCOEF(Z2, &temp1, &temp2, &temp3, &temp4, &ATRHO, &VFERMI, &temp5, &uuu);

	SCOEFH(Z2, PCOEF); // һ����ȡ����
}

void EnergyLoss::Sn()
{
}

void EnergyLoss::Range()
{
}

// reslut is in
void EnergyLoss::RPSTOP(int Z2, NORMALTYPE E, NORMALTYPE PCOEF[8], NORMALTYPE *SP)
{

	NORMALTYPE PEO, PE, Result;

	PEO = 25.0;
	PE = AMAX1(PEO, E);

	SL = PCOEF[1] * pow(PE, PCOEF[2]) + PCOEF[3] * pow(PE, PCOEF[4]);
	SH = PCOEF[5] / pow(PE, PCOEF[6]) * log(PCOEF[7] / PE + PCOEF[8] * PE);
	Result = SL * SH / (SL + SH);

	if (E <= PEO)
	{
		if (Z2 <= 6)
			VELPWR = 0.25;
		else
			VELPWR = 0.45;
		Result = Result * pow(E / PEO, VELPWR);
	}
	*SP = Result; // se of P at energy of E
}

void EnergyLoss::H_Se(NORMALTYPE *result)
{

	NORMALTYPE E0, E;
	int I;
	E0 = EE / 1000.0 / M1; // EE the max energy,MM1 is the target atom mass

	for (I = 1; I <= 1000; I++)
	{
		E = E0 * I;
		RPSTOP(Z2, E, PCOEF, &result[I]);
		result[I] = result[I] * 10; // 10 is to change the uini,unit of eV-A2
	}

	/*
	   A=0.5*result[1000]*ATRHO*pow(10,-5);  //the unit changed to eV/A

   for(I=1;I<=1000;I++)
	   {
		   result[I]=result[I]*A;  //10 is to change the uini,unit of eV-A2
	   }
   */

	//	B=result[1000]/10;
	//	printf("Z1=%4d,LAMBDA=%5.2f;Z2=%4d,VFERMI=%5.2f\n",Z1,LFCTR,Z2,VFERMI);
	//	printf("FOR ION ENERGY=%6.0f,KEV,STOPPING=%7,2f,EV/A=%7.2f,KEV-CM2\n",EE,A,B);
}

void EnergyLoss::He_Se(NORMALTYPE *SE)
{
	NORMALTYPE HEO, A, B, HEN, E0, E;
	int I;

	E0 = EE / M1 / 1000.0; // EE the max energy,MM1 is the target atom mass

	for (I = 1; I <= 1000; I++)
	// 8888888888888888888
	{
		E = I * E0;

		HEO = 1.0;
		HE = AMAX1(HEO, E);
		B = log(HE);

		A = 0.2865 + 0.1266 * B - 0.001429 * B * B + 0.02402 * B * B * B - 0.01135 * B * B * B * B + 0.001475 * B * B * B * B * B;
		HEN = 1.0 - exp(-AMIN1(30, A));

		A = 1.0 + (0.007 + 0.00005 * Z2) * exp(-(7.6 - AMAX1(0.0, log(HE))) * (7.6 - AMAX1(0.0, log(HE))));
		HEN = HEN * A * A;

		// get that for proton
		RPSTOP(Z2, HE, PCOEF, &SP);

		// cal the ratio
		//		SE[I]=SP*HEN*4;
		SE[I] = SP * HEN * Z1 * Z1;
		if (E <= HE0)
			SE[I] = SE[I] * sqrt(E / HEO);
		SE[I] = SE[I] * 10;
	}
	// 8888888888888888888888888888888888888888888
	//	A=SE[1000]*ATRHO*pow(10,-24);
	//	B=SE[1000]/10;
	//	printf("Z1=%4d,LAMBDA=%5.2f;Z2=%4d,VFERMI=%5.2f\n",Z1,LFCTR,Z2,VFERMI);
	//	printf("FOR ION ENERGY=%6.0f,KEV,STOPPING=%7,2f,EV/A=%7.2f,KEV-CM2\n",EE,A,B);
}

void EnergyLoss::Heavyion_Se(NORMALTYPE *SE)
{
	NORMALTYPE YRMIN = 0.13, VRMIN = 1.0;
	NORMALTYPE A, tempB, YR, L1, VR, EEE, LO, POWER;
	NORMALTYPE E, E0, Q;
	int I;

	E0 = EE / M1 / 1000.0; // EE the max energy,MM1 is the target atom mass

	for (I = 1; I <= 1000; I++)
	// 8888888888888888888
	{
		E = I * E0;

		Velocity = sqrt(E / 25.0) / VFERMI;

		if (Velocity < 1.0) //!!!
			VR = (3 * VFERMI / 4.0) * (1 + (2 * Velocity * Velocity / 3.0) - Velocity * Velocity * Velocity * Velocity / 15.0);
		else
			VR = Velocity * VFERMI * (1 + 1 / (5 * Velocity * Velocity));

		YR = AMAX1(YRMIN, VR / pow(Z1, 0.6667));
		YR = AMAX1(YR, VRMIN / pow(Z1, 0.6667));

		A = -0.803 * pow(YR, 0.3) + 1.3167 * pow(YR, 0.6) + 0.38157 * YR + 0.008983 * YR * YR;
		Q = AMIN1(1.0, AMAX1(0.0, 1.0 - exp(-AMIN1(A, 50.0))));

		tempB = (AMIN1(0.43, AMAX1(0.32, 0.12 + 0.025 * Z1))) / pow(Z1, 0.3333);
		LO = (0.8 - Q * AMIN1(1.2, 0.6 + Z1 / 30.0)) / pow(Z1, 0.3333);

		if (Q < 0.2)
			L1 = 0.0;
		else if (Q < (AMAX1(0.0, 9 - 0.025 * Z1)))
		{
			// Q1=0.2;  ///!!!!
			L1 = tempB * (Q - 0.2) / fabs(AMAX1(0.0, 0.9 - 0.025 * Z1) - 0.2000001);
		}
		else if (Q < (AMAX1(0.0, 1.0 - 0.025 * AMIN1(16.0, 1.0 * Z1))))
			L1 = tempB;
		else
			L1 = tempB * (1.0 - Q) / (0.025 * AMIN1(16.0, 1.0 * Z1));

		L = AMAX1(L1, LO * LFCTR);

		ZETA = Q + (1.0 / (2.0 * VFERMI * VFERMI)) * (1.0 - Q) * log(1 + pow((4 * L * VFERMI / 1.919), 2));
		A = -pow((7.6 - AMAX1(0.0, log(E))), 2);

		ZETA = ZETA * (1.0 + (1.0 / (Z1 * Z1)) * (0.18 + 0.0015 * Z2) * exp(A));

		if (YR > AMAX1(YRMIN, VRMIN / pow(Z1, 0.6667)))
		{
			RPSTOP(Z2, E, PCOEF, &SP);
			SE[I] = SP * (ZETA * Z1) * (ZETA * Z1);
		}
		else
		{
			VRMIN = AMAX1(VRMIN, YRMIN * pow(Z1, 0.6667));
			VMIN = 0.5 * (VRMIN + sqrt(AMAX1(0.0, VRMIN * VRMIN - 0.8 * VFERMI * VFERMI)));

			EEE = 25 * VMIN * VMIN;

			RPSTOP(Z2, EEE, PCOEF, &SP);
			POWER = 0.5;
			if ((Z2 == 6) || (((Z2 == 14) || (Z2 == 32)) && (Z1 <= 19)))

				POWER = 0.375;
			SE[I] = SP * (ZETA * Z1) * (ZETA * Z1) * pow(E / EEE, POWER);
		}

		SE[I] = SE[I] * 10;

	} // 888888888888888888888888888888888888888888

	//		A=SE[1000]*ATRHO*pow(10,-24);
	//		B=SE[1000]/10;
	//		printf("Z1=%4d,LAMBDA=%5.2f;Z2=%4d,VFERMI=%5.2f\n",Z1,LFCTR,Z2,VFERMI);
	//		printf("FOR ION ENERGY=%6.0f,KEV,STOPPING=%7,2f,EV/A=%7.2f,KEV-CM2\n",EE,A,B);
}

void EnergyLoss::RStop(int Z1a, NORMALTYPE M1a, int Z2a, NORMALTYPE EEa, NORMALTYPE *SEa, NORMALTYPE LFCTRa, NORMALTYPE &VFERMIa)
{

	NORMALTYPE temp1, temp2, temp3, temp4, temp5, uuu;

	Z1 = Z1a;
	M1 = M1a;
	Z2 = Z2a;
	EE = EEa;
	SCOEF(Z1, &temp1, &MM1, &temp2, &temp3, &temp4, &temp5, &LFCTR, &uuu);

	SCOEF(Z2, &temp1, &temp2, &temp3, &temp4, &ATRHO, &VFERMI, &temp5, &uuu);

	SCOEFH(Z2, PCOEF); // һ����ȡ����
	if (M1a == 0)
		M1a = MM1;
	// E0=0.001*EE/M1;

	Se(SEa);
	return;
}

void EnergyLoss::InitialSE()
{
	int LayerNo, ElementNo, II, I;
	NORMALTYPE IZT;
	for (LayerNo = 1; LayerNo <= MaxLayer; LayerNo++)
	{
		// ARHO[LayerNo]=RHO[LayerNo]*0.6022/M2[LayerNo];
		// MU[LayerNo]=IonM/M2[LayerNo];

		II = LayerElementNo[LayerNo];
		for (ElementNo = 1; ElementNo <= II; ElementNo++)
		{
			IZT = TargetAtomNumber[LayerNo][ElementNo];
			RStop(IonZ, IonM, IZT, IonEnergyKeV, SEO, XLAMB, VFERMI);

			VF[LayerNo][ElementNo] = VFERMI;

			for (I = 1; I <= 1000; I++)
				SE[LayerNo][I] = SE[LayerNo][I] + SEO[I] * TargetAtomConc[LayerNo][ElementNo] * ARHO[LayerNo]; //!!!!
		}																									   // element

	} // layer loop
}