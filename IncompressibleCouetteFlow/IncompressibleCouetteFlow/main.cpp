#include "IncompressibleCouetteFlow.h"

int main() {
	IncompressibleCouetteFlow ICF;
	Data data;
	ICF.pressureCorrection(data);
	return 0;
}
