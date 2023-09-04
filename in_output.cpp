#include "in_output.h"
#include "io.h"
#include "stdio.h"
#include <iostream>

// appead the atom to the AtomList.txt file for analysis
void in_output::RecordingAtom(AtomType temp)
{

	FILE *handle;

	if ((handle = fopen("AtomList.txt", "ab")) != NULL)
	{
		fwrite(&temp, sizeof(AtomType), 1, handle);
		fclose(handle);
	}
	else
	{
		printf("****************    the atomlist.txt is missing    **********************");
	}
};

void in_output::ResutlsOut()
{

	FILE *handle;
	int i;
	// range

	if ((handle = fopen("ionrange.txt", "w+")) != NULL)
	{
		for (i = 1; i <= 100; i++)
			fprintf(handle, "%d %lf\n", i, ionrange[i]);
		fclose(handle);
		printf("\n\n*****   the ion range has been output!    *******\n");
	}

	// vacancy
	if ((handle = fopen("vacancy.txt", "w+")) != NULL)
	{
		for (i = 1; i <= 100; i++)
			fprintf(handle, "%d %lf\n", i, IVAC[i] + RVAC[i]);
		fclose(handle);
		printf("*****   the vacancy distribution has been output! *******\n");
	}
	// ionization
	if ((handle = fopen("ionization.txt", "w+")) != NULL)
	{
		for (i = 1; i <= 100; i++)
			fprintf(handle, "%d %lf\n", i, IonizatinEDis[i] + RION[i]);
		fclose(handle);
		printf("*****   the ionization has been output! *******\n");
	}

	// phonons
	if ((handle = fopen("phonons.txt", "w+")) != NULL)
	{
		for (i = 1; i <= 100; i++)
			fprintf(handle, "%d %lf\n", i, PhoneEDis[i] + RPHON[i]);
		fclose(handle);
		printf("*****   the phonons distribution has been output! *******\n");
	}

	printf("\n\n*************   END of THE SIMULATION       ***************\n");
}

// this will displaying all the inforamtion needed for the simulation for checking
// li 20230824
void in_output::DisplayInformation()
{
	int i, j; // temp for looping
	char flag;

	std::cout << "\n*****************  ion information *******************" << std::endl
			  << std::endl;
	std::cout << "Ion energy (KeV)  : " << IonEnergyKeV << std::endl;
	std::cout << "Ion Mass          : " << IonM << std::endl;
	std::cout << "Ion Z             : " << IonZ << std::endl;
	std::cout << "Ion incident angle: " << ALPHA << std::endl;
	std::cout << "Incident Ion number:" << ionnumber << std::endl;

	std::cout << std::endl
			  << "****************  sample information  *******************" << std::endl;

	for (i = 1; i <= num_layers; ++i)
	{
		std::cout << std::endl
				  << "layer " << i << "  :" << std::endl;

		for (j = 1; j <= LayerElementNo[i]; ++j)
		{
			std::cout << "element Z: " << TargetAtomNumber[i][j]
					  << " mass is: " << TargetAtomMass[i][j]
					  << " concentration is: " << TargetAtomConc[i][j]
					  << " density: " << layerdensity[i] << std::endl;
		}
	}

	std::cout << std::endl
			  << "************  end of simulation information  ************" << std::endl;

	// If you do not want to ask for confirmation, comment below
	std::cout << std::endl
			  << "Are these parameters OK?  y(n)" << std::endl;
	std::cin >> flag;
	if (flag != 'y')
	{
		std::cout << std::endl
				  << "^^^^^^^^^^^^The parameters are not correct, reset them!^^^^^^^^^^^^^^^^" << std::endl;
		exit(0);
	}
}

void in_output::init()
{
	// Declare local variables
	FILE *filehandle;
	int LayerNo, ElementNo, L, II, I, num_elements;
	double TMIN, X, tempmass;
	NORMALTYPE YY[MAX_ELEMENT];

	// current, the formating should be obeyed, later will change to a more friendly interface
	// Open 'parameters.txt' to read simulation parameters
	if ((filehandle = fopen("parameters.txt", "r")) != NULL)
	{
		// Read various ion properties and simulation parameters from the file
		fscanf(filehandle, "ionenergy(kev):%lf  \n", &IonEnergyKeV);
		fscanf(filehandle, "ionmass:%lf  \n", &IonM);
		fscanf(filehandle, "ionz:%lf  \n", &IonZ);
		fscanf(filehandle, "ionangle:%lf  \n", &ALPHA); // in degree!! alpa is in pai values!
		fscanf(filehandle, "ionnumber:%d  \n", &ionnumber);
		fscanf(filehandle, "displaceenergy:%lf  \n", &DisEnergy);
		fscanf(filehandle, "channelwidth:%lf  \n", &channelwidth);
		fscanf(filehandle, "num_layers:%d\n", &num_layers); // Read the number of layers
		L = num_layers;

		// Loop to read each layer's properties and elements
		for (int i = 1; i <= num_layers; i++)
		{
			fscanf(filehandle, "layerthickness:%lf density:%lf num_elements:%d\n", &LayerThickness[i], &layerdensity[i], &num_elements);
			LayerElementNo[i] = num_elements;

			for (int j = 1; j <= num_elements; j++)
			{
				fscanf(filehandle, "elementmass:%lf elementZ:%lf elementconcentration:%lf\n", &TargetAtomMass[i][j], &TargetAtomNumber[i][j], &TargetAtomConc[i][j]);
			}
		}

		fclose(filehandle);
	}
	else
	{
		std::cout << "The parameters file has an error!" << std::endl;
		exit(3);
	}

	MaxLayer = L; // meanwhile, the layer number and the element number in each layer
				  // is also obtained!!!

	// setup the program constants

	IZ = IonZ;
	SCOEF(IZ, &X, &tempmass, &X, &X, &X, &X, &X, YY); // B is the most abundant isotop of Z1
													  //  X and YY are meaningless, only B is needed
	if (IonM == 0)
		IonM = tempmass;

	IonEnergyEV = IonEnergyKeV * 1000; // now is in unit of eV
	StoppingEnergy = 5.0;			   // AMAX1(5.,IonEnergyKeV*.1);     not neceasray to chose,fix at 5

	TMIN = 5;
	DistanceNuclear = 0; // ion travel distance during the collsion
	DA = 3;

	// Calculate layer-to-surface distance for each layer
	LayertoSurface[1] = LayerThickness[1]; // max lenght from surface to the layer I
	for (LayerNo = 2; LayerNo <= 3; LayerNo++)
		LayertoSurface[LayerNo] = LayerThickness[LayerNo] + LayertoSurface[LayerNo - 1];

	if (channelwidth == 0)
		channelwidth = 0.01 * LayertoSurface[3]; // for recording, the total sample are dividied into 100 segment

	//	L0=channelwidth;

	// Initialize layer average mass and z (atomic number)
	for (LayerNo = 1; LayerNo <= L; LayerNo++) // layer
	{
		II = LayerElementNo[LayerNo]; // II is the element number in layer LL
		for (ElementNo = 1; ElementNo <= II; ElementNo++)

			H[LayerNo] = H[LayerNo] + TargetAtomConc[LayerNo][ElementNo]; // the sum of fraction could be not equal to 100%!!
	}

	for (LayerNo = 1; LayerNo <= L; LayerNo++)
	{
		II = LayerElementNo[LayerNo];
		for (ElementNo = 1; ElementNo <= II; ElementNo++)
		{
			TargetAtomConc[LayerNo][ElementNo] = TargetAtomConc[LayerNo][ElementNo] / H[LayerNo];

			// M2: average Mass; Z2: average z
			M2[LayerNo] = M2[LayerNo] + TargetAtomConc[LayerNo][ElementNo] * TargetAtomMass[LayerNo][ElementNo];
			Z_2[LayerNo] = Z_2[LayerNo] + TargetAtomConc[LayerNo][ElementNo] * TargetAtomNumber[LayerNo][ElementNo];
		}
	}
	printf("\n***** GETTING ZBL STOPPING POWERS FOR TARGET *****");

	// Calculate material densities and other layer-related parameters
	for (LayerNo = 1; LayerNo <= MaxLayer; LayerNo++)
	{
		ARHO[LayerNo] = layerdensity[LayerNo] * 0.6022 / M2[LayerNo];
		MU[LayerNo] = IonM / M2[LayerNo];
	} // layer loop

	// caulculate the mean flight path os that the scatting angle is less than 5
	for (LayerNo = 1; LayerNo <= L; LayerNo++)
	{

		// for free length calculation
		A[LayerNo] = 0.5292 * 0.8853 / (pow(IonZ, 0.23) + pow(Z_2[LayerNo], 0.23)); // screen constant
		F[LayerNo] = A[LayerNo] * M2[LayerNo] / (IonZ * Z_2[LayerNo] * 14.4 * (IonM + M2[LayerNo]));
		//	ESP0=IonEnergyEV*F[LayerNo];
		EPSDG[LayerNo] = TMIN * F[LayerNo] * pow((1.0 + MU[LayerNo]), 2) / (4.0 * MU[LayerNo]);

		// the following decided howmuc of the recoid erngy goes into se and vacancyies
		DamageEFract[LayerNo] = 0.01 * pow(Z_2[LayerNo], (-7.0 / 3.0));						   // se
		ElectronEFract[LayerNo] = 0.1334 * pow(Z_2[LayerNo], (2.0 / 3.0)) / sqrt(M2[LayerNo]); // vacancy
	}

	for (LayerNo = 1; LayerNo <= L; LayerNo++)
	{
		II = LayerElementNo[LayerNo];
		for (I = 1; I <= II; I++)
		{

			//
			M1_to_M2[LayerNo][I] = IonM / TargetAtomMass[LayerNo][I];
			recoilfactor[LayerNo][I] = 4.0 * M1_to_M2[LayerNo][I] / pow((1.0 + M1_to_M2[LayerNo][I]), 2);
			ScreenLen[LayerNo][I] = 0.5292 * 0.8853 / (pow(IonZ, 0.23) + pow(TargetAtomNumber[LayerNo][I], 0.23));
			FItemp[LayerNo][I] = ScreenLen[LayerNo][I] * TargetAtomMass[LayerNo][I] / (IonZ * TargetAtomNumber[LayerNo][I] * 14.4 * (IonM + TargetAtomMass[LayerNo][I]));
		}
	}
	printf("******");
	printf("\nSETUP FINISHED. NOW STARTING MONTE-CARLO LOOPS.");
}

//*************************    xue  OK   ************************************

void in_output::SCOEF(int z, NORMALTYPE *massN, NORMALTYPE *weight, NORMALTYPE *mass, NORMALTYPE *density, NORMALTYPE *atomicdensity, NORMALTYPE *FermiV, NORMALTYPE *Screenactor, NORMALTYPE *HeatSub)
{
	FILE *filehandle;
	int flag = 0, z1;

	if ((filehandle = fopen("scoef.88", "r")) == NULL)
		printf("missing of the data file");
	else
	{
		while (flag++ < 92)
		{
			fscanf(filehandle, "   %d %lf  %lf %lf %lf %lf %lf %lf %lf\n", &z1, massN, weight, mass, density, atomicdensity, FermiV, Screenactor, HeatSub);
			if (z1 == z)
				flag = 100;
		}
	}
}

//*************************    xue  OK   ************************************
void in_output::SCOEFH(int z, NORMALTYPE a[])
{
	FILE *filehandle;
	int flag = 0, z1;

	if ((filehandle = fopen("h.88", "r")) == NULL)
		printf("missing of the data file");
	else
	{
		while (flag++ < 92)
		{
			fscanf(filehandle, "   %d %lf  %lf %lf %lf %lf %lf %lf %lf\n", &z1, &a[1], &a[2], &a[3], &a[4], &a[5], &a[6], &a[7], &a[8]);
			//		printf(" %d %lf  %lf %lf %lf %lf %lf %lf\n\n\n",z1,massN,weight,mass,density,FermiV,Screenactor,HeatSub);
			if (z1 == z)
				flag = 100;
		}
	}
}

void in_output::ResutlsAnalysis()
{
	double AverageRange, Straggling;
	int n;
	n = (ElapsedIonNo - RBSIonNo - TransIonNo);
	AverageRange = (XSUM / n);
	Straggling = sqrt(X2SUM / n - AverageRange * AverageRange);

	printf("\nTransIonNo= %d    RBSIonNo= %d \n", TransIonNo, RBSIonNo);
	printf("AverageRange= %lf", AverageRange);
	printf("  Straggling= %lf\n", Straggling);

	AverageRange = (PLSUM / n);
	printf("Ion average traveled length: %lf\n", AverageRange);

	// 打开文件用于输出结果
	FILE *filehandle = fopen("result.txt", "w");
	if (filehandle == NULL)
	{
		printf("Error: Could not open file for writing.\n");
		return;
	}
	// 文件输出
	fprintf(filehandle, "\nTransIonNo= %d    RBSIonNo= %d \n", TransIonNo, RBSIonNo);
	fprintf(filehandle, "AverageRange= %lf", AverageRange);
	fprintf(filehandle, "  Straggling= %lf\n", Straggling);
	fprintf(filehandle, "\n\nthe energies deposited into diferent ways by ions and by recoiled atoms:(%)");
	fprintf(filehandle, "\nI_ION= %lf   R_ION= %lf\nI_VAC= %lf   R_VAC= %lf\nI_PHON= %lf   R_PHON= %lf\n", I_ION, R_ION, I_VAC, R_VAC, I_PHON, R_PHON);
	// 关闭文件
	fclose(filehandle);

	int i;
	double E_total;
	I_PHON = 0;
	R_PHON = 0;
	I_VAC = 0;
	R_VAC = 0;
	I_ION = 0;
	R_ION = 0;
	for (i = 1; i < MAX_ARRAY; i++)
	{
		// since vac[] is the number of the vacancy, not energy, so they should times DIsenergy

		I_PHON = I_PHON + PhoneEDis[i];
		R_PHON = R_PHON + RPHON[i];

		I_VAC = I_VAC + IVAC[i] * DisEnergy;
		R_VAC = R_VAC + RVAC[i] * DisEnergy;

		I_ION = I_ION + IonizatinEDis[i];
		R_ION = R_ION + RION[i];
	}

	E_total = I_PHON + R_PHON + I_VAC + R_VAC + I_ION + R_ION;

	I_PHON = 100 * I_PHON / E_total;
	R_PHON = 100 * R_PHON / E_total;
	I_VAC = 100 * I_VAC / E_total;
	R_VAC = 100 * R_VAC / E_total;
	I_ION = 100 * I_ION / E_total;
	R_ION = 100 * R_ION / E_total;

	printf("\n\nthe energies deposited into diferent ways by ions and by recoiled atoms:(%)");
	printf("\nI_ION= %lf   R_ION= %lf\nI_VAC= %lf   R_VAC= %lf\nI_PHON= %lf   R_PHON= %lf\n", I_ION, R_ION, I_VAC, R_VAC, I_PHON, R_PHON);
}
