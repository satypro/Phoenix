#include "stdafx.h"
#include "SearchNeighbour.h"
#include "SearchNeighbourFactory.h"
#include "SearchStructure.h"
#include "Processing.h"
#include "ProcessRequest.h"
#include "ClassifierTypes.h"
#include "SearchStructure.h"
#include "SearchType.h"
#include <pcl/io/pcd_io.h>

int TestCase3DVT_Process_Cloud()
{
	/*Request Test Cases*/
	ProcessRequest request;
	request.cloudFileName = "data/label/d4-alg2l.txt";
	request.optimalFileName = "data/label/optimal.txt";
	request.labelFileName = "data/3d_4.txt";

	request.classifierType = static_cast<int>(C_3DVT);
	request.dataStructure = static_cast<int>(OctTree);
	request.epsi = 0.69;
	request.rmin = 0.03;
	request.rmax = 0.09;
	request.rmaxpt = 0.03;
	request.scale = 3;
	request.searchType = static_cast<int>(Radius);

	Processing* p = new Processing();
	p->buildStructure(request);
	p->processCloud(request);
	return 0;
}

int TestCase3DVT_Process_Point_Wise()
{
	/*Request Test Cases*/
	ProcessRequest request;
	request.cloudFileName = "data/label/d4-alg2l.txt";
	request.optimalFileName = "data/label/optimal.txt";
	request.labelFileName = "data/3d_4.txt";
	request.classifierType = static_cast<int>(C_3DVT);
	request.dataStructure = static_cast<int>(OctTree);
	request.epsi = 0.69;
	request.rmin = 0.03;
	request.rmax = 0.09;
	request.rmaxpt = 0.03;
	request.scale = 3;
	request.searchType = static_cast<int>(Radius);

	Processing* p = new Processing();
	// Build the Structre --- like whether to use the KdTree or Octree
	p->buildStructure(request);
	
	pcl::PointCloud<pcl::PointXYZ>::Ptr _cloud = p->getCloud();
	
	for (int index = 0; index < _cloud->size(); index++)
	{
		// We may change the request here to put the classifier type 
		request.classifierType = static_cast<int>(C_3DVT);
		// Also like  serach type from Radius to KNeasers
		// request.searchType = static_cast<int>(KNearest);
		// This way we can perform the analysis against each point.
		p->processPoint(request, index);
	}

	return 0;
}

int main()
{
	TestCase3DVT_Process_Cloud();
	TestCase3DVT_Process_Point_Wise();
}