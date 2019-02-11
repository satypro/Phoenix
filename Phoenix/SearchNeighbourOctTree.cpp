#include "SearchNeighbourOctTree.h"

SearchNeighbourOctTree::SearchNeighbourOctTree()
{
}

SearchNeighbourOctTree::~SearchNeighbourOctTree()
{
}

void SearchNeighbourOctTree::build(pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud, Configuration* _config)
{
	_cloud = cloud;
	this->_config = _config;

	char* pEnd;
	float resolution = ::strtof(_config->GetValue("Resolution").c_str(), &pEnd);

	_octree = new pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>(resolution);
	_octree->setInputCloud(cloud);
	_octree->addPointsFromInputCloud();
}

bool SearchNeighbourOctTree::radiusSearch(
	pcl::PointXYZ searchPoint, 
	float radius, 
	std::vector<int> pointIdxRadiusSearch, 
	std::vector<float> pointRadiusSquaredDistance)
{
	if (_octree == NULL)
		return false;

	int result = _octree->radiusSearch(searchPoint, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance);
	return result > 0 ? true : false;
}

bool SearchNeighbourOctTree::nearestKSearch(
	pcl::PointXYZ searchPoint, 
	int K, 
	std::vector<int> pointIdxNKNSearch, 
	std::vector<float> pointNKNSquaredDistance)
{
	if (_octree == NULL)
		return false;

	int result = _octree->nearestKSearch(searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance);
	return result > 0 ? true : false;
}

bool SearchNeighbourOctTree::voxelSearch(pcl::PointXYZ searchPoint, std::vector<int> pointIdxVec)
{
	if (_octree == NULL)
		return false;

	bool result = _octree->voxelSearch(searchPoint, pointIdxVec);
	return result;
}