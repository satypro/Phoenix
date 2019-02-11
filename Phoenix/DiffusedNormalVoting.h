#ifndef DIFFUSEDNORMALVOTING_H
#define DIFFUSEDNORMALVOTING_H

#include <stdio.h>
#include <vector>
#include "OctreeSerial.h"
#include "data_type.h"
#include "SearchNeighbour.h"
#include "SearchStructure.h"
#include "SearchNeighbourFactory.h"
#include "Classifier.h"

using namespace std;
using namespace Eigen;

class DiffusedNormalVoting : public Classifier
{
public:
	typedef Eigen::Matrix<double, 3, 1> ColumnVector;
	typedef Eigen::Matrix<double, 1, 3> RowVector;
	typedef pcl::PointXYZ PointT;
	typedef pcl::octree::OctreeContainerPointIndices LeafContainerT;
	typedef pcl::octree::OctreeContainerEmpty BranchContainerT;

	DiffusedNormalVoting() {}

	void setInputCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud);

	void setOctree(OctreeXYZ<PointT, LeafContainerT, BranchContainerT > * octree);

	void setParams(float rmin, float rmax, float rmaxpt, float epsi, float scale);

	void process(ProcessRequest processRequest, vector <probfeatnode> *probval, vector <tensorType> *aggregated_tensor);
	void processPoint(ProcessRequest request, probfeatnode* _probval, tensorType* _accum, int index);

	void prob_measure(vector <probfeatnode> *probval, vector <tensorType> *aggregated_tensor);
	//New Prob measure
	void prob_measure(probfeatnode* probval, tensorType* tensor, int index);

	void filterPointTypes(float radius, std::vector<probfeatnode>&probfeatVal, vector <tensorType> *aggregated_tensor);

	void calculatePartialDerivative(float radius) {
		//Method not Implemented
	}

	void calculateSecondDerivative(float radius) {
		//Method not Implemented
	}

	void calculateThirdDerivative(float radius) {
		//Method not Implemented
	}

private:

	void getdiffusionvelocity(Eigen::Vector3f evals, metaVelData *diffVel);

	void getnormalizeTensor(size_t idx, float weight);

	bool setAccumzero();

	tensorType compute3DBallVote(ColumnVector V, float *weight);

	bool addBallVote(size_t idx, tensorType ballTensorVote, float lambdaN);

	bool makeVector(pcl::PointXYZ p_i, pcl::PointXYZ p_j, ColumnVector *V);

	bool voteAccum(float radius, vector<probfeatnode>&probfeatVal);

	// New Method
	bool voteAccum(float radius, probfeatnode probVal, int index);
	// New Method
	bool probMeasure(probfeatnode &probfeatVal, float radius, int index);

	bool probMeasure(vector<probfeatnode>&probfeatVal, float radius);

	void writeGlyphVars(float radius);

	void computeSaliencyVals(glyphVars &glyphtemp);

	void glyphAnalysis(glyphVars &glyphtemp);

	bool eigendecomposition(glyphVars &glyphtemp, size_t idx);

	float _rmin, _rmax, _rmaxpt, _epsi, _scale;

	float _sigma, _delta;

	pcl::PointCloud<pcl::PointXYZ>::Ptr _inCloud;

	OctreeXYZ<PointT, LeafContainerT, BranchContainerT > * _octree;
	vector < tensorType	> _accum;
	vector<glyphVars> _glyphData;
	vector < tensorType	> averaged_tensors;
	vector < tensorType	> averaged_accum;

	tensorType _n_accum;
	glyphVars _n_glyphData;
	tensorType _n_avg_acum;
	tensorType _n_avg_tensor;

	float max_cl, min_cl;
};

#endif