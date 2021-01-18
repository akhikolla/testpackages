
#include "ObliqueVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	ObliqueVariable::ObliqueVariable(integer n, integer num)
	{
		SphereVariable SV(n);

		Element **SVs = new Element *[num];
		for (integer i = 0; i < num; i++)
		{
			SVs[i] = &SV;
		}
		integer *powsintev = new integer[2];
		powsintev[0] = 0;
		powsintev[1] = num;

		ProductElementInitialization(SVs, num, powsintev, 1);

		delete[] powsintev;
		delete[] SVs;
	};

	ObliqueVariable::~ObliqueVariable(void)
	{
	};

	ObliqueVariable *ObliqueVariable::ConstructEmpty(void) const
	{
		return new ObliqueVariable(elements[0]->Getlength(), numofelements);
	};

	//void ObliqueVariable::RandInManifold(void)
	//{
	//	SphereVariable *SV = nullptr;
	//	for (integer i = 0; i < numofelements; i++)
	//	{
	//		SV = dynamic_cast<SphereVariable *> (elements[i]);
	//		SV->RandInManifold();
	//	}
	//};

	void ObliqueVariable::Print(const char *name, bool isonlymain) const
	{
		if (isonlymain)
		{
			if (Space == nullptr)
			{
				if (size == nullptr)
				{
					OUTSTREAM << name << " is an empty data with size 0";
				}
				else
				{
					OUTSTREAM << name << " is an empty data with size " << size[0];
				}
				for (integer i = 1; i < ls; i++)
					OUTSTREAM << " x " << size[i];
				OUTSTREAM << std::endl;
				return;
			}
			OUTSTREAM << name << ", shared times:" << *sharedtimes << ", shared times address:" << sharedtimes << std::endl;
			integer n = elements[0]->Getlength();
			integer num = numofelements;
			for (integer i = 0; i < n; i++)
			{
				for (integer j = 0; j < num; j++)
				{
					OUTSTREAM << elements[j]->GetSpace()[i] << "\t";
				}
				OUTSTREAM << std::endl;
			}
			return;
		}
		ProductElement::Print(name, isonlymain);
	};
} /*end of ROPTLIB namespace*/
