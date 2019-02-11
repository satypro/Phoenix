#pragma once

#include <stdio.h>
#include <vector>
#include "OctreeSerial.h"

#include "data_type.h"

class Utils
{
public:
	typedef pcl::PointXYZ PointT;
	typedef pcl::octree::OctreeContainerPointIndices LeafContainerT;
	typedef pcl::octree::OctreeContainerEmpty BranchContainerT;

	Utils();

	void setInputCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud);
	void setOctree(OctreeXYZ<PointT, LeafContainerT, BranchContainerT > * octree);
	void setParams(float rmin, float rmax, float rmaxpt, float epsi, float scale);

	pcl::PointCloud<pcl::PointXYZ>::Ptr _inCloud;
	OctreeXYZ<PointT, LeafContainerT, BranchContainerT > * _octree;
	float _rmin, _rmax, _rmaxpt, _epsi, _scale;
	float max_cl, min_cl;
	int numTriangles;

	void triangulate(float radius, std::vector<probfeatnode>&probfeatVal);
	void pruneTriangles(float radius, std::vector<probfeatnode>&probfeatVal);

	void generateTensorLines1(vector<probfeatnode>&probfeatVal);
	void generateTensorLines2(vector<probfeatnode>&probfeatVal);
	void generateTensorLines3(vector<probfeatnode>&probfeatVal);
	void simplifyTensorLines(vector<probfeatnode>&probfeatVal);

	void computeDONS(vector<probfeatnode>&probfeatVal, float radius);

	void computeDoEs(vector<probfeatnode>&probfeatVal);

	void computeDoNs_MSO(vector<probfeatnode>&probfeatVal, float rmin, float rmax);

	void ground_density(vector<probfeatnode>&probfeatVal);

	void label_classification(vector<probfeatnode>&probfeatVal);

	void writeSalToFile(vector<probfeatnode>&probfeatVal);

	void writeDONToFile(vector<probfeatnode>&probfeatVal);

	void writeNeighborsToFile(float radius);

	void writeTriangulationToFile(vector<probfeatnode>&probfeatVal);

	void completeTriangulation();

	void normalize(float(&vec)[3]);

	void correctMisclassifiedPoints(ftType eigenvalues[3], ftType evecs[9],
		float Lp, float Ll, tensorType &out_tensor);

	void eigen_decomposition_to_file(vector<tensorType> &aggregated_tensors);
};