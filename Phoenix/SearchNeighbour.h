#pragma once
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <vector>
#include "Configuration.h"

class SearchNeighbour
{
public:
	virtual void build(pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud, Configuration* _config) = 0;
	virtual bool radiusSearch(pcl::PointXYZ searchPoint, float radius, std::vector<int> pointIdxRadiusSearch, std::vector<float> pointRadiusSquaredDistance) = 0;
	virtual bool nearestKSearch(pcl::PointXYZ searchPoint, int K, std::vector<int> pointIdxNKNSearch, std::vector<float> pointNKNSquaredDistance ) = 0;
	virtual bool voxelSearch(pcl::PointXYZ searchPoint, std::vector<int> pointIdxVec) = 0;
};

