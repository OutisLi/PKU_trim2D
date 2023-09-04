#include "definition.h"

class in_output : public Data
{

public:
	void init();

	void SCOEF(int, NORMALTYPE *, NORMALTYPE *, NORMALTYPE *, NORMALTYPE *, NORMALTYPE *, NORMALTYPE *, NORMALTYPE *, NORMALTYPE *); // һ����ȡ����
	void SCOEFH(int, NORMALTYPE *);																									 // һ����ȡ����

	void DisplayInformation();
	void ResutlsAnalysis();
	void ResutlsOut();
	void RecordingAtom(AtomType);

private:
	NORMALTYPE I_PHON, R_PHON, I_VAC, R_VAC, I_ION, R_ION;
};