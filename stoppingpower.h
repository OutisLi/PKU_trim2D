#include "in_output.h"
class EnergyLoss : public in_output
{

public:
	void RStop(int Z1, NORMALTYPE M1, int Z2, NORMALTYPE EE, NORMALTYPE *SE, NORMALTYPE LFCTR, NORMALTYPE &VFERMI);
	void RPSTOP(int, NORMALTYPE, NORMALTYPE *, NORMALTYPE *);

	void SetParameters(int az1, NORMALTYPE am1, int az2, NORMALTYPE aee);
	void InitialSE();
	void Se(NORMALTYPE *);
	void Sn();
	void Range();
	void H_Se(NORMALTYPE *);
	void He_Se(NORMALTYPE *);
	void Heavyion_Se(NORMALTYPE *);

private:
	// �������ڱ��ؼ����õģ���
	int Z1, Z2;
	NORMALTYPE E0KEV, DA, Velocity;
	NORMALTYPE MM1, M1, EE, L; // ION
	NORMALTYPE ATRHO;
	NORMALTYPE PEO, PE, SL, SH, VELPWR, HE, HEN, SP, HE0, VMIN;
};