#include "montecarlo.h"
#include "io.h"
#include "stdlib.h"
#include "malloc.h"
#include <time.h>
extern void Get_Information(struct Information *In_Data);

int main()
{
	long ioncounter = 1;
	PKU_MC MCRUNNING;

#ifdef MONIORING
	MCRUNNING.TestRand();
#endif

	// init all the things for the simulation except the ion incident condition
	MCRUNNING.SimulationInit();
	// displaying all the inforamtion for the simulation for correction checking
	MCRUNNING.DisplayInformation();

	//******************   begin the simualtion  *********************************

	//******************        ion loop         *********************************
	while (ioncounter <= MCRUNNING.ionnumber) // for incident ions
	{
#ifdef MONIORING
		MCRUNNING.CurrentCondition(ioncounter);
#endif

		MCRUNNING.DoOneIon();
		MCRUNNING.ElapsedIonNo = ++ioncounter;
	}
	//******************   END the simualtion  *********************************

	MCRUNNING.ResutlsAnalysis();
	MCRUNNING.ResutlsOut();
	printf("\n************     ENDING MONTE-CARLO SIMULATION     **************\n");
	system("pause");
}