#pragma once
#include "SearchNeighbour.h"
#include <pcl/octree/octree_search.h>

class SearchNeighbourOctTree :
	public SearchNeighbour
{
public:
	SearchNeighbourOctTree();
	~SearchNeighbourOctTree();

	// Inherited via SearchNeighbour
	virtual void build(pcl::PointCloud<pcl::PointXYZ>::Ptr & cloud, Configuration* _config) override;
	virtual bool radiusSearch(pcl::PointXYZ searchPoint, float radius, std::vector<int> pointIdxRadiusSearch, std::vector<float> pointRadiusSquaredDistance) override;
	virtual bool nearestKSearch(pcl::PointXYZ searchPoint, int K, std::vector<int> pointIdxNKNSearch, std::vector<float> pointNKNSquaredDistance) override;
	virtual bool voxelSearch(pcl::PointXYZ searchPoint, std::vector<int> pointIdxVec) override;

private:
	pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>* _octree;
	pcl::PointCloud<pcl::PointXYZ>::Ptr _cloud;
	Configuration* _config;
};