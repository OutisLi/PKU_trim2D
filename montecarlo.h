#include "stoppingpower.h"

// structure for the recoild atom recording

// typedef struct
// {
// 	double x, y, layerno, angle, mass, z;
// } CascadeAtom;

class PKU_MC : public EnergyLoss
{
public:
	void SimulationInit(); // all the init are finished here
	void InitialPara();
	void InitialIon();
	void DoOneIon();			// this will finished one ion simulation
	void CurrentCondition(int); // ������ʾ����ʾ��ǰ״��

private:
	double ChoseCollisionAtom(); // ѡ��Ŀ��ԭ�ӵ����࣬������ԭ������
	double FreeLength();
	double EleEnergyLoss();
	double EmissionAngle();
	double IonNewCondition();
	void ProcessRecorder();
	int OutOrNot(int &);	   //  �ж������Ƿ����
	void CollisionRecording(); //  �ۼ�ͳ������
	void BubbleEffect();
	void BubbleList();

	void CascadedAtom(); // this will finish the high energy recoiled atoms

private:
	int IE, I, J;
	int CurrentElement; // which atom is been knocked

	NORMALTYPE MVA, AEX, VARI, SIGMA, V, V2, GAMMA, BETA, Y, X, AVEEPL, SIGPL;
	NORMALTYPE SEE, DEE, EEG, DEN, E;
	NORMALTYPE B, P, RR, EX1, EX2, EX3, EX4, V1, FR, FR1, CC, ROC, AA, FF, C2, S2, CT, ST, SQE, CO;

	NORMALTYPE TransIonE,
		RBSIonE,
		IonFreeLength,
		IonWayLength;
};