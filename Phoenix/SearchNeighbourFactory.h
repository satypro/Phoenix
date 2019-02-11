#pragma once
#include<map>
#include "SearchNeighbour.h"
#include "SearchStructure.h"

class SearchNeighbourFactory
{
public:
	SearchNeighbourFactory();
	~SearchNeighbourFactory();
	static SearchNeighbour* GetNeighbourSearchDataStructure(SearchStructure option);
private:
	static std::map<SearchStructure, SearchNeighbour*> _neighbourSearchMap;
	static SearchNeighbour* GetInstance(SearchStructure structure);
};