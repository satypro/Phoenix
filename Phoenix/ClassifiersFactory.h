#pragma once
#include <map>
#include "Classifier.h"
#include "ClassifierTypes.h"

class ClassifiersFactory
{
public:
	ClassifiersFactory();
	~ClassifiersFactory();
	static Classifier* GetClassifier(ClassifierTypes classifierType);
private:
	static std::map<ClassifierTypes, Classifier*> _classifierInstanceMap;
	static Classifier* GetInstance(ClassifierTypes classifierType);
};

