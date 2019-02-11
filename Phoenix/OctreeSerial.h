#ifndef _MY_OCTREE_H
#define _MY_OCTREE_H

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/octree/octree_impl.h>
#include <pcl/octree/octree_iterator.h>
#include <pcl/octree/octree_pointcloud_pointvector.h>
#include <pcl/octree/octree_container.h>
#include <iostream>
#include <vector>
#include <ctime>

template<typename PointT = pcl::PointXYZ, typename LeafContainerT = pcl::octree::OctreeContainerPointIndices, typename BranchContainerT = pcl::octree::OctreeContainerEmpty >
class OctreeXYZ : public pcl::octree::OctreePointCloudSearch<PointT, LeafContainerT, BranchContainerT>
{
public:
	OctreeXYZ(const double resolution) : pcl::octree::OctreePointCloudSearch < PointT, LeafContainerT, BranchContainerT>(resolution)
	{
	}

	~OctreeXYZ() {}
private:

};

#endif