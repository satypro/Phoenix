#include "ClassifiersFactory.h"
#include "DiffusedNormalVoting.h"
#include "CovarianceMatrix.h"
#include "Tensor3dvoting.h"

std::map<ClassifierTypes, Classifier*> ClassifiersFactory::_classifierInstanceMap;

ClassifiersFactory::ClassifiersFactory()
{
}

ClassifiersFactory::~ClassifiersFactory()
{
}

Classifier* ClassifiersFactory::GetInstance(ClassifierTypes classifierType)
{
	std::map<ClassifierTypes, Classifier*>::iterator it;
	it = ClassifiersFactory::_classifierInstanceMap.find(classifierType);

	if (it != ClassifiersFactory::_classifierInstanceMap.end())
		return it->second;

	Classifier* instance = ClassifiersFactory::GetClassifier(classifierType);
	ClassifiersFactory::_classifierInstanceMap.insert(
		std::pair<ClassifierTypes, Classifier*>(classifierType, instance)
	);

	return instance;
}

Classifier* ClassifiersFactory::GetClassifier(ClassifierTypes classifierType)
{
	Classifier* object;
	switch (classifierType)
	{
	case C_3DVTGET:
		object = new DiffusedNormalVoting();
		break;
	case C_3DVT:
		object = new DiffusedNormalVoting();
		break;
	case C_3DCM:
		object = new CovarianceMatrix();
		break;
	case C_3DMCM:
		object = new CovarianceMatrix();
		break;
	case C_2DGET:
		object = new CovarianceMatrix();
		break;
	case C_Hessian:
		object = new CovarianceMatrix();
		break;
	case C_2DCM:
		object = new CovarianceMatrix();
		break;
	default:
		object = new CovarianceMatrix();
	}

	return object;
}