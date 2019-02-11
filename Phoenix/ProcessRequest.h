#pragma once
#include<string>

class ProcessRequest
{
public:
	float rmin;
	float rmax;
	float rmaxpt;
	float epsi;
	float scale;
	std::string cloudFileName;
	std::string optimalFileName;
	std::string labelFileName;
	int dataStructure; // Octreee or Kdtree
	int searchType; // Radius or KNearest or RadiusOptimal or KNearestOptimal
	int classifierType;
};