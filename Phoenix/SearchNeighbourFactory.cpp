#include "SearchNeighbourFactory.h"
#include "SearchNeighbourOctTree.h"

std::map<SearchStructure, SearchNeighbour*> SearchNeighbourFactory::_neighbourSearchMap;

SearchNeighbourFactory::SearchNeighbourFactory()
{
}

SearchNeighbourFactory::~SearchNeighbourFactory()
{
}

SearchNeighbour * SearchNeighbourFactory::GetNeighbourSearchDataStructure(SearchStructure structure)
{
	std::map<SearchStructure, SearchNeighbour*>::iterator it;
	it = SearchNeighbourFactory::_neighbourSearchMap.find(structure);

	if (it != SearchNeighbourFactory::_neighbourSearchMap.end())
		return it->second;

	SearchNeighbour* instance = SearchNeighbourFactory::GetInstance(structure);

	SearchNeighbourFactory::_neighbourSearchMap.insert(
		std::pair<SearchStructure, SearchNeighbour*>(structure, instance)
	);

	return instance;
}

SearchNeighbour* SearchNeighbourFactory::GetInstance(SearchStructure structure)
{
	SearchNeighbour* object;
	switch (structure)
	{
	case KdTree:
		object = nullptr;
		break;
	case OctTree:
		object = new SearchNeighbourOctTree();
	default:
		object = new SearchNeighbourOctTree();
	}

	return object;
}