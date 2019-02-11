#include "Utils.h"
#include <algorithm>
#include <vector>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/search/organized.h>
#include <pcl/search/kdtree.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/filters/conditional_removal.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/features/don.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Segment_3.h>
#include <pcl/common/common.h>
#include <pcl/common/eigen.h>
#include <pcl/common/centroid.h>
#include "eig3.h"

#define GROUND_THRESHOLD 0.0132672
#define CONTOUR_VALUE 0.4
#define DOUGLASPEUCKER_EPSILON 0.005
#define NUM_SEED_POINTS 3500
#define SEED_POINTS_INTERVAL 80
#define STEP_SIZE 0.003
#define MAX_LENGTH 100
#define DON_LOW 0.3
#define DON_HIGH 0.9

struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};

//typedef CGAL::Triangulation_euclidean_traits_xy_3<K>  Gt;
//typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;

typedef K::Point_3   Point;
typedef K::Triangle_3  Triangle_3;
typedef K::Vector_3  Vector_3;
typedef K::Segment_3  Segment_3;

void Barycentric(Point p, Point a, Point b, Point c, float &u, float &v, float &w);
int LineTriangleIntersection(Triangle_3 tri, Point p, Point q, Point& intersectionP);
std::vector<Point> DouglasPeucker(std::vector<Point> pointList, float epsilon, int start, int end);

Utils::Utils()
{
}

void Utils::setInputCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud)
{
	_inCloud = cloud;
}

void Utils::setOctree(OctreeXYZ<PointT, LeafContainerT, BranchContainerT > * octree)
{
	_octree = octree;
}

void Utils::setParams(float rmin, float rmax, float rmaxpt, float epsi, float scale)
{
	_rmin = rmin;
	_rmax = rmax;
	_rmaxpt = rmaxpt;
	_epsi = epsi;
	_scale = scale;
}

void Utils::normalize(float(&vec)[3])
{
	float mag = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
	vec[0] /= mag;
	vec[1] /= mag;
	vec[2] /= mag;
}

/*****************************************************************************************/
void Utils::triangulate(float radius, std::vector<probfeatnode>&probfeatVal)
{
	std::vector<int> pointIdxRadiusSearch;
	std::vector<float> pointRadiusSquaredDistance;
	pcl::PointXYZ searchPoint, tempPoint;

	max_cl = -1;
	min_cl = 100;

	numTriangles = 0;
	int countDONlines = 0;

	float don_low = DON_LOW;
	float don_high = DON_HIGH;

	// Populate truthVal
	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		//        probfeatVal[i].truthVal = ( ((probfeatVal[i].csclcp[1] >= 0.3 && probfeatVal[i].csclcp[1] <= 0.5) || (probfeatVal[i].csclcp[1] >= 0.8 && probfeatVal[i].csclcp[1] <= 1.0)) || (probfeatVal[i].don>=don_low && probfeatVal[i].don<=don_high) || (probfeatVal[i].anisotropy <= .3)) ? 1 : 0;
		//        probfeatVal[i].truthVal = (probfeatVal[i].don>=don_low && probfeatVal[i].don<=don_high) ? 1 : 0;
		//        probfeatVal[i].truthVal = ( (probfeatVal[i].anisotropy <= .3)) ? 1 : 0;
	}

	//    for(size_t i =0; i < _inCloud->points.size(); i++)
	//    {
	//        searchPoint = _inCloud->points[i];
	//        std::vector<Point> triangle_points;
	//        std::vector<int> triangle_points_ndx;
	//
	//        if(probfeatVal[i].don>=don_low && probfeatVal[i].don<=don_high)
	//            countDONlines++;
	////        bool triPoints = (probfeatVal[i].don>=don_low && probfeatVal[i].don<=don_high) || (probfeatVal[i].csclcp[1] > probfeatVal[i].csclcp[0] && probfeatVal[i].csclcp[1] > probfeatVal[i].csclcp[2]);
	////        bool triPoints = (probfeatVal[i].csclcp[1] > probfeatVal[i].csclcp[0] && probfeatVal[i].csclcp[1] > probfeatVal[i].csclcp[2]);
	//        bool triPoints = (probfeatVal[i].csclcp[1]>=.15 && probfeatVal[i].csclcp[1]<=0.5) || (probfeatVal[i].don>=don_low && probfeatVal[i].don<=don_high) || (probfeatVal[i].csclcp[1] > probfeatVal[i].csclcp[0] && probfeatVal[i].csclcp[1] > probfeatVal[i].csclcp[2]);
	//        if( triPoints )
	//        {
	//            if(max_cl<probfeatVal[i].csclcp[1])
	//                max_cl = probfeatVal[i].csclcp[1];
	//            if(min_cl>probfeatVal[i].csclcp[1])
	//                min_cl = probfeatVal[i].csclcp[1];
	//
	//            if(_octree->radiusSearch(searchPoint, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0)
	//            {
	//                for (size_t j = 0; j < pointIdxRadiusSearch.size (); j++)
	//                {
	////                    if((probfeatVal[pointIdxRadiusSearch[j]].don>=don_low && probfeatVal[pointIdxRadiusSearch[j]].don<=don_high) || (probfeatVal[pointIdxRadiusSearch[j]].csclcp[1] > probfeatVal[pointIdxRadiusSearch[j]].csclcp[0] && probfeatVal[pointIdxRadiusSearch[j]].csclcp[1] > probfeatVal[pointIdxRadiusSearch[j]].csclcp[2]))
	////                    if((probfeatVal[pointIdxRadiusSearch[j]].csclcp[1] > probfeatVal[pointIdxRadiusSearch[j]].csclcp[0] && probfeatVal[pointIdxRadiusSearch[j]].csclcp[1] > probfeatVal[pointIdxRadiusSearch[j]].csclcp[2]))
	//                    if((probfeatVal[pointIdxRadiusSearch[j]].csclcp[1]>=.15 && probfeatVal[pointIdxRadiusSearch[j]].csclcp[1]<=.65 ) || (probfeatVal[pointIdxRadiusSearch[j]].don>=don_low && probfeatVal[pointIdxRadiusSearch[j]].don<=don_high) || (probfeatVal[pointIdxRadiusSearch[j]].csclcp[1] > probfeatVal[pointIdxRadiusSearch[j]].csclcp[0] && probfeatVal[pointIdxRadiusSearch[j]].csclcp[1] > probfeatVal[pointIdxRadiusSearch[j]].csclcp[2]))
	//                    {
	//                        tempPoint = _inCloud->points[pointIdxRadiusSearch[j]];
	//                        triangle_points.push_back(Point(tempPoint.x,tempPoint.y,tempPoint.z));
	//                        triangle_points_ndx.push_back(pointIdxRadiusSearch[j]);
	//                    }
	//                }
	//                triangle_points.push_back(Point(searchPoint.x,searchPoint.y,searchPoint.z));
	//                triangle_points_ndx.push_back(i);
	//                Delaunay dt;
	//                dt.insert(triangle_points.begin(), triangle_points.end());
	//                std::vector<Triangle> triangles;
	//                for(Delaunay::Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); ++it)
	//                {
	//                    Point p1 = it->vertex(0)->point();
	//                    Point p2 = it->vertex(1)->point();
	//                    Point p3 = it->vertex(2)->point();
	//                    int i1=-1;
	//                    int i2=-1;
	//                    int i3=-1;
	//                    for(int k=0;k<triangle_points.size();k++){
	//                        if(p1==triangle_points[k]){
	//                            i1 = triangle_points_ndx[k];
	//                        }
	//                        if(p2==triangle_points[k]){
	//                            i2 = triangle_points_ndx[k];
	//                        }
	//                        if(p3==triangle_points[k]){
	//                            i3 = triangle_points_ndx[k];
	//                        }
	//                    }
	//                    Triangle temp;
	//                    temp.p[0] = p1.x();
	//                    temp.p[1] = p1.y();
	//                    temp.p[2] = p1.z();
	//                    temp.q[0] = p2.x();
	//                    temp.q[1] = p2.y();
	//                    temp.q[2] = p2.z();
	//                    temp.r[0] = p3.x();
	//                    temp.r[1] = p3.y();
	//                    temp.r[2] = p3.z();
	//                    temp.pid  = i1;
	//                    temp.qid  = i2;
	//                    temp.rid  = i3;
	//                    temp.average_cl = (probfeatVal[i1].csclcp[1] + probfeatVal[i2].csclcp[1] + probfeatVal[i3].csclcp[1])/3.0;
	//                    temp.average_anisotropy = (probfeatVal[i1].anisotropy + probfeatVal[i2].anisotropy + probfeatVal[i3].anisotropy)/3.0;
	//
	//                    //score=modulus(dot product (unit vector-edge,unit vector-eigenvector at vertex)
	//
	//                    Vector_3 pq(p1,p2);
	//                    Vector_3 qr(p2,p3);
	//                    Vector_3 rp(p3,p1);
	//                    pq = pq/sqrt(pq.squared_length());
	//                    qr = qr/sqrt(qr.squared_length());
	//                    rp = rp/sqrt(rp.squared_length());
	//
	//                    float score_pq = abs(pq.x()*probfeatVal[i1].evecs[0] + pq.y()*probfeatVal[i1].evecs[1] + pq.z()*probfeatVal[i1].evecs[2]);
	//                    float score_pr = abs(rp.x()*probfeatVal[i1].evecs[0] + rp.y()*probfeatVal[i1].evecs[1] + rp.z()*probfeatVal[i1].evecs[2]);
	//
	//                    float score_qp = abs(pq.x()*probfeatVal[i2].evecs[0] + pq.y()*probfeatVal[i2].evecs[1] + pq.z()*probfeatVal[i2].evecs[2]);
	//                    float score_qr = abs(qr.x()*probfeatVal[i2].evecs[0] + qr.y()*probfeatVal[i2].evecs[1] + qr.z()*probfeatVal[i2].evecs[2]);
	//
	//                    float score_rp = abs(rp.x()*probfeatVal[i1].evecs[0] + rp.y()*probfeatVal[i2].evecs[1] + rp.z()*probfeatVal[i3].evecs[2]);
	//                    float score_rq = abs(qr.x()*probfeatVal[i1].evecs[0] + qr.y()*probfeatVal[i2].evecs[1] + qr.z()*probfeatVal[i3].evecs[2]);
	//
	//                    temp.tlScore = max(score_pq,max(score_pr,max(score_qp,max(score_qr,max(score_rp,score_rq)))));
	//
	//                    if(temp.tlScore == score_pq || temp.tlScore == score_qp)
	//                        temp.tlEdge = 0;
	//                    else if(temp.tlScore == score_qr || temp.tlScore == score_rq)
	//                        temp.tlEdge = 1;
	//                    else
	//                        temp.tlEdge = 2;
	//
	//                    float average_height_tri = (temp.p[2] + temp.q[2] + temp.r[2])/3;
	//
	//                    // Marching Triangles
	//                    float scalar1 = probfeatVal[i1].csclcp[1];
	//                    float scalar2 = probfeatVal[i2].csclcp[1];
	//                    float scalar3 = probfeatVal[i3].csclcp[1];
	//
	//                    if(scalar1 > CONTOUR_VALUE && scalar2 > CONTOUR_VALUE && scalar3 > CONTOUR_VALUE ) {
	//                        //+++
	//                        temp.hasCL = false;
	//                    } else if(scalar1 > CONTOUR_VALUE && scalar2 > CONTOUR_VALUE && scalar3 <= CONTOUR_VALUE) {
	//                        //++-
	//                        temp.hasCL = true;
	//                        float per1 = (scalar1 - CONTOUR_VALUE) / (scalar1 - scalar3);
	//                        float per2 = (scalar2 - CONTOUR_VALUE) / (scalar2 - scalar3);
	//
	//                        temp.CL_p1[0] = p1.x()*(1-per1) + p3.x()*per1;
	//                        temp.CL_p1[1] = p1.y()*(1-per1) + p3.y()*per1;
	//                        temp.CL_p1[2] = p1.z()*(1-per1) + p3.z()*per1;
	//
	//                        temp.CL_p2[0] = p2.x()*(1-per2) + p3.x()*per2;
	//                        temp.CL_p2[1] = p2.y()*(1-per2) + p3.y()*per2;
	//                        temp.CL_p2[2] = p2.z()*(1-per2) + p3.z()*per2;
	//
	//                    } else if(scalar1 > CONTOUR_VALUE && scalar2 <= CONTOUR_VALUE && scalar3 > CONTOUR_VALUE) {
	//                        //+-+
	//                        temp.hasCL = true;
	//                        float per1 = (scalar1 - CONTOUR_VALUE) / (scalar1 - scalar2);
	//                        float per2 = (scalar3 - CONTOUR_VALUE) / (scalar3 - scalar2);
	//
	//                        temp.CL_p1[0] = p1.x()*(1-per1) + p2.x()*per1;
	//                        temp.CL_p1[1] = p1.y()*(1-per1) + p2.y()*per1;
	//                        temp.CL_p1[2] = p1.z()*(1-per1) + p2.z()*per1;
	//
	//                        temp.CL_p2[0] = p3.x()*(1-per2) + p2.x()*per2;
	//                        temp.CL_p2[1] = p3.y()*(1-per2) + p2.y()*per2;
	//                        temp.CL_p2[2] = p3.z()*(1-per2) + p2.z()*per2;
	//
	//                    } else if(scalar1 > CONTOUR_VALUE && scalar2 <= CONTOUR_VALUE && scalar3 <= CONTOUR_VALUE) {
	//                        //+--
	//                        temp.hasCL = true;
	//                        float per1 = (scalar1 - CONTOUR_VALUE) / (scalar1 - scalar2);
	//                        float per2 = (scalar1 - CONTOUR_VALUE) / (scalar1 - scalar3);
	//
	//                        temp.CL_p1[0] = p1.x()*(1-per1) + p2.x()*per1;
	//                        temp.CL_p1[1] = p1.y()*(1-per1) + p2.y()*per1;
	//                        temp.CL_p1[2] = p1.z()*(1-per1) + p2.z()*per1;
	//
	//                        temp.CL_p2[0] = p1.x()*(1-per2) + p3.x()*per2;
	//                        temp.CL_p2[1] = p1.y()*(1-per2) + p3.y()*per2;
	//                        temp.CL_p2[2] = p1.z()*(1-per2) + p3.z()*per2;
	//
	//                    } else if(scalar1 <= CONTOUR_VALUE && scalar2 > CONTOUR_VALUE && scalar3 > CONTOUR_VALUE) {
	//                        //-++
	//                        temp.hasCL = true;
	//                        float per1 = (scalar2 - CONTOUR_VALUE) / (scalar2 - scalar1);
	//                        float per2 = (scalar3 - CONTOUR_VALUE) / (scalar3 - scalar1);
	//
	//                        temp.CL_p1[0] = p2.x()*(1-per1) + p1.x()*per1;
	//                        temp.CL_p1[1] = p2.y()*(1-per1) + p1.y()*per1;
	//                        temp.CL_p1[2] = p2.z()*(1-per1) + p1.z()*per1;
	//
	//                        temp.CL_p2[0] = p3.x()*(1-per2) + p1.x()*per2;
	//                        temp.CL_p2[1] = p3.y()*(1-per2) + p1.y()*per2;
	//                        temp.CL_p2[2] = p3.z()*(1-per2) + p1.z()*per2;
	//
	//                    } else if(scalar1 <= CONTOUR_VALUE && scalar2 > CONTOUR_VALUE && scalar3 <= CONTOUR_VALUE) {
	//                        //-+-
	//                        temp.hasCL = true;
	//                        float per1 = (scalar2 - CONTOUR_VALUE) / (scalar2 - scalar1);
	//                        float per2 = (scalar2 - CONTOUR_VALUE) / (scalar2 - scalar3);
	//
	//                        temp.CL_p1[0] = p2.x()*(1-per1) + p1.x()*per1;
	//                        temp.CL_p1[1] = p2.y()*(1-per1) + p1.y()*per1;
	//                        temp.CL_p1[2] = p2.z()*(1-per1) + p1.z()*per1;
	//
	//                        temp.CL_p2[0] = p2.x()*(1-per2) + p3.x()*per2;
	//                        temp.CL_p2[1] = p2.y()*(1-per2) + p3.y()*per2;
	//                        temp.CL_p2[2] = p2.z()*(1-per2) + p3.z()*per2;
	//                    } else if(scalar1 <= CONTOUR_VALUE && scalar2 <= CONTOUR_VALUE && scalar3 > CONTOUR_VALUE) {
	//                        //--+
	//                        temp.hasCL = true;
	//                        float per1 = (scalar3 - CONTOUR_VALUE) / (scalar3 - scalar1);
	//                        float per2 = (scalar3 - CONTOUR_VALUE) / (scalar3 - scalar2);
	//
	//                        temp.CL_p1[0] = p3.x()*(1-per1) + p1.x()*per1;
	//                        temp.CL_p1[1] = p3.y()*(1-per1) + p1.y()*per1;
	//                        temp.CL_p1[2] = p3.z()*(1-per1) + p1.z()*per1;
	//
	//                        temp.CL_p2[0] = p3.x()*(1-per2) + p2.x()*per2;
	//                        temp.CL_p2[1] = p3.y()*(1-per2) + p2.y()*per2;
	//                        temp.CL_p2[2] = p3.z()*(1-per2) + p2.z()*per2;
	//                    } else if(scalar1 <= CONTOUR_VALUE && scalar2 <= CONTOUR_VALUE && scalar3 <= CONTOUR_VALUE) {
	//                        //---
	//                        temp.hasCL = false;
	//                    }
	//
	//                    if(average_height_tri > 0.02)
	//                        triangles.push_back(temp);
	//                }
	//                probfeatVal[i].triangles = triangles;
	//                numTriangles += triangles.size();
	//
	//                pointIdxRadiusSearch.clear();
	//                pointRadiusSquaredDistance.clear();
	////                cout << probfeatVal[i].triangle_edges.size() << " " << endl;
	//            }
	//        }
	//    }
	//    std::cout << "count dons = " << countDONlines << std::endl;
	//    std::cout << "# Triangles = " << numTriangles << std::endl;
}


// Delete triangles with less line saliency than a certain threshold
void Utils::pruneTriangles(float radius, std::vector<probfeatnode>&probfeatVal)
{
	float threshold = (max_cl + min_cl)*.25;
	numTriangles = 0;

	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		for (vector<Triangle>::iterator it = (probfeatVal[i].triangles).begin();
			it != (probfeatVal[i].triangles).end(); )
		{
			if (it->average_cl < threshold || (it->average_anisotropy < 0.17 || it->average_anisotropy > 0.35))
				it = (probfeatVal[i].triangles).erase(it);
			else
				++it;
		}
		numTriangles += probfeatVal[i].triangles.size();
	}
}

// Fill ground point density in cylindrical neighbourhood
void Utils::ground_density(vector<probfeatnode>&probfeatVal)
{
	//    pcl::PointXYZ searchPoint;
	//    pcl::PointXYZ targetPoint;

	//    for(size_t i =0; i < _inCloud->points.size(); i++)
	//    {
	//        searchPoint = _inCloud->points[i];
	//        int numGroundPoints = 0;
	//        int totalPoints = 0;

	//        for(size_t j =0; j < _inCloud->points.size(); j++)
	//        {
	//            targetPoint = _inCloud->points[j];
	//            float squared_distance = (searchPoint.x-targetPoint.x)*(searchPoint.x-targetPoint.x) + (searchPoint.y-targetPoint.y)*(searchPoint.y-targetPoint.y);
	//            if(squared_distance < _rmin*_rmin)
	//            {
	//                totalPoints++;
	//                if(targetPoint.z < GROUND_THRESHOLD)
	//                    numGroundPoints++;
	//            }
	//        }
	//        if(totalPoints!=0)
	//            probfeatVal[i].groundDensity = numGroundPoints/totalPoints;
	//        else
	//            probfeatVal[i].groundDensity = 0;
	//    }
}

// Tensor lines algorithm 1
void Utils::generateTensorLines1(vector<probfeatnode>&probfeatVal)
{
	// Get some random seed point
	int seed_point = -1;
	for (size_t j = 0; j < _inCloud->points.size(); j++)
	{
		if (probfeatVal[j].triangles.size() > 0) {
			seed_point = j;
			break;
		}
	}

	float step_size = 0.008;
	float max_length = 10.0;
	float current_length = 0.0;

	std::vector<Line> tensor_line;
	float current_point[3] = { _inCloud->points[seed_point].x,_inCloud->points[seed_point].y,_inCloud->points[seed_point].z }; //=seed point
	float next_point[3] = { 0,0,0 };
	float major_evec[3] = { probfeatVal[seed_point].evecs[0],probfeatVal[seed_point].evecs[1],probfeatVal[seed_point].evecs[2] };
	std::vector<Triangle> neighbouring_triangles;
	neighbouring_triangles.insert(neighbouring_triangles.end(), probfeatVal[seed_point].triangles.begin(), probfeatVal[seed_point].triangles.end());

	while (current_length < max_length)
	{
		next_point[0] = step_size * major_evec[0] + current_point[0];
		next_point[1] = step_size * major_evec[1] + current_point[1];
		next_point[2] = step_size * major_evec[2] + current_point[2];

		//        std::cout << "me= " << major_evec[0] << "," << major_evec[1] << "," << major_evec[2] << std::endl;
		//        std::cout << "cp= " << current_point[0] << "," << current_point[1] << "," << current_point[2] << std::endl;
		//        std::cout << "np= " << next_point[0] << "," << next_point[1] << "," << next_point[2] << std::endl;

		bool np_valid = false;
		Point cp(current_point[0], current_point[1], current_point[2]);
		Point np(next_point[0], next_point[1], next_point[2]);
		Point intersectionP;
		for (int i = 0; i < neighbouring_triangles.size(); i++)
		{
			Point p(neighbouring_triangles[i].p[0], neighbouring_triangles[i].p[1], neighbouring_triangles[i].p[2]);
			Point q(neighbouring_triangles[i].q[0], neighbouring_triangles[i].q[1], neighbouring_triangles[i].q[2]);
			Point r(neighbouring_triangles[i].r[0], neighbouring_triangles[i].r[1], neighbouring_triangles[i].r[2]);
			Triangle_3 tri(p, q, r);
			if (LineTriangleIntersection(tri, cp, np, intersectionP) == 1) {
				//                std::cout << "ip" << intersectionP << std::endl;
				//                std::cout << "np" << np << std::endl;
				//                std::cout << "cp" << cp << std::endl;
				next_point[0] = intersectionP.x();
				next_point[1] = intersectionP.y();
				next_point[2] = intersectionP.z();
				np = intersectionP;
			}
			if (tri.has_on(np)) {
				np_valid = true;
				std::cout << "np valid" << endl;
				Line step;
				step.p[0] = current_point[0];
				step.p[1] = current_point[1];
				step.p[2] = current_point[2];
				step.q[0] = next_point[0];
				step.q[1] = next_point[1];
				step.q[2] = next_point[2];
				tensor_line.push_back(step);

				float u, v, w;
				Barycentric(np, p, q, r, u, v, w);
				current_length += step_size;
				int pid = neighbouring_triangles[i].pid;
				int qid = neighbouring_triangles[i].qid;
				int rid = neighbouring_triangles[i].rid;
				major_evec[0] = probfeatVal[pid].evecs[0] * u + probfeatVal[qid].evecs[0] * v + probfeatVal[rid].evecs[0] * w;
				major_evec[1] = probfeatVal[pid].evecs[1] * u + probfeatVal[qid].evecs[1] * v + probfeatVal[rid].evecs[1] * w;
				major_evec[2] = probfeatVal[pid].evecs[2] * u + probfeatVal[qid].evecs[2] * v + probfeatVal[rid].evecs[2] * w;
				current_point[0] = next_point[0];
				current_point[1] = next_point[1];
				current_point[2] = next_point[2];

				std::cout << "uvw= " << u << "," << v << "," << w << std::endl;
				std::cout << "me= " << major_evec[0] << "," << major_evec[1] << "," << major_evec[2] << std::endl;
				std::cout << "cp= " << current_point[0] << "," << current_point[1] << "," << current_point[2] << std::endl;
				std::cout << "np= " << next_point[0] << "," << next_point[1] << "," << next_point[2] << std::endl;

				neighbouring_triangles.clear();
				neighbouring_triangles.insert(neighbouring_triangles.end(), probfeatVal[pid].triangles.begin(), probfeatVal[pid].triangles.end());
				neighbouring_triangles.insert(neighbouring_triangles.end(), probfeatVal[qid].triangles.begin(), probfeatVal[qid].triangles.end());
				neighbouring_triangles.insert(neighbouring_triangles.end(), probfeatVal[rid].triangles.begin(), probfeatVal[rid].triangles.end());

				break;
			}
		}

		if (!np_valid)
			break;
	}
	probfeatVal[seed_point].tensor_line = tensor_line;
}

// Tensor lines algorithm 2
void Utils::generateTensorLines2(vector<probfeatnode>&probfeatVal)
{
	// Get some random seed point
	std::vector<int> seed_points;
	for (size_t j = 0; j < _inCloud->points.size(); j++)
	{
		if (probfeatVal[j].triangles.size() > 0) {
			seed_points.push_back(j);
			j += SEED_POINTS_INTERVAL;
			if (seed_points.size() >= NUM_SEED_POINTS)
				break;
		}
	}
	std::cout << "# seed points = " << seed_points.size() << std::endl;
	for (int sp = 0; sp < seed_points.size(); sp++)
	{
		int seed_point = seed_points[sp];
		probfeatVal[seed_point].tl_visited = true;

		float step_size = STEP_SIZE;
		float max_length = MAX_LENGTH;
		float current_length = 0.0;

		std::vector<Line> tensor_line;
		float current_point[3] = { _inCloud->points[seed_point].x,_inCloud->points[seed_point].y,_inCloud->points[seed_point].z }; //=seed point
		float next_point[3] = { 0,0,0 };
		float major_evec[3] = { probfeatVal[seed_point].evecs[0],probfeatVal[seed_point].evecs[1],probfeatVal[seed_point].evecs[2] };
		std::vector<Triangle> neighbouring_triangles;
		neighbouring_triangles.insert(neighbouring_triangles.end(), probfeatVal[seed_point].triangles.begin(), probfeatVal[seed_point].triangles.end());

		while (current_length < max_length)
		{
			next_point[0] = step_size * major_evec[0] + current_point[0];
			next_point[1] = step_size * major_evec[1] + current_point[1];
			next_point[2] = step_size * major_evec[2] + current_point[2];

			float minDis = 1000;
			float np_id = -1;
			Point cp(current_point[0], current_point[1], current_point[2]);
			Point np(next_point[0], next_point[1], next_point[2]);
			for (int i = 0; i < neighbouring_triangles.size(); i++)
			{
				Point p(neighbouring_triangles[i].p[0], neighbouring_triangles[i].p[1], neighbouring_triangles[i].p[2]);
				Point q(neighbouring_triangles[i].q[0], neighbouring_triangles[i].q[1], neighbouring_triangles[i].q[2]);
				Point r(neighbouring_triangles[i].r[0], neighbouring_triangles[i].r[1], neighbouring_triangles[i].r[2]);
				Segment_3 cpp(np, p);
				Segment_3 cpq(np, q);
				Segment_3 cpr(np, r);
				if (!probfeatVal[neighbouring_triangles[i].pid].tl_visited && p != cp && cpp.squared_length() < minDis) {
					minDis = cpp.squared_length();
					np_id = neighbouring_triangles[i].pid;
					probfeatVal[np_id].tl_visited = true;
				}
				if (!probfeatVal[neighbouring_triangles[i].qid].tl_visited && q != cp && cpq.squared_length() < minDis) {
					minDis = cpq.squared_length();
					np_id = neighbouring_triangles[i].qid;
					probfeatVal[np_id].tl_visited = true;
				}
				if (!probfeatVal[neighbouring_triangles[i].rid].tl_visited && r != cp && cpr.squared_length() < minDis) {
					minDis = cpr.squared_length();
					np_id = neighbouring_triangles[i].rid;
					probfeatVal[np_id].tl_visited = true;
				}
			}
			if (np_id == -1)
				break;
			else {
				Line step;
				step.p[0] = current_point[0];
				step.p[1] = current_point[1];
				step.p[2] = current_point[2];
				step.q[0] = _inCloud->points[np_id].x;
				step.q[1] = _inCloud->points[np_id].y;
				step.q[2] = _inCloud->points[np_id].z;
				step.ndxq = np_id;
				tensor_line.push_back(step);
				current_length += step_size;
				current_point[0] = _inCloud->points[np_id].x;
				current_point[1] = _inCloud->points[np_id].y;
				current_point[2] = _inCloud->points[np_id].z;
				neighbouring_triangles.clear();
				neighbouring_triangles.insert(neighbouring_triangles.end(), probfeatVal[np_id].triangles.begin(), probfeatVal[np_id].triangles.end());
			}
		}
		probfeatVal[seed_point].tensor_line = tensor_line;
	}
}

// Tensor lines algorithm 3 (using this one)
void Utils::generateTensorLines3(vector<probfeatnode>&probfeatVal)
{
	// Get some random seed point
	std::vector<int> seed_points;
	for (size_t j = 0; j < _inCloud->points.size(); j++)
	{
		if (probfeatVal[j].triangles.size() > 0) {
			seed_points.push_back(j);
			j += SEED_POINTS_INTERVAL;
			if (seed_points.size() >= NUM_SEED_POINTS)
				break;
		}
	}
	std::cout << "# seed points = " << seed_points.size() << std::endl;
	for (int sp = 0; sp < seed_points.size(); sp++)
	{
		int seed_point = seed_points[sp];
		probfeatVal[seed_point].tl_visited = true;

		float step_size = STEP_SIZE;
		float max_length = MAX_LENGTH;
		float current_length = 0.0;

		std::vector<Line> tensor_line;
		float current_point[3] = { _inCloud->points[seed_point].x,_inCloud->points[seed_point].y,_inCloud->points[seed_point].z }; //=seed point
		float next_point[3] = { 0,0,0 };
		float major_evec[3] = { probfeatVal[seed_point].evecs[0],probfeatVal[seed_point].evecs[1],probfeatVal[seed_point].evecs[2] };
		std::vector<Triangle> neighbouring_triangles;
		neighbouring_triangles.insert(neighbouring_triangles.end(), probfeatVal[seed_point].triangles.begin(), probfeatVal[seed_point].triangles.end());

		while (current_length < max_length)
		{
			next_point[0] = step_size * major_evec[0] + current_point[0];
			next_point[1] = step_size * major_evec[1] + current_point[1];
			next_point[2] = step_size * major_evec[2] + current_point[2];

			float minDis = 1000;
			float np_id = -1;
			Point cp(current_point[0], current_point[1], current_point[2]);
			Point np(next_point[0], next_point[1], next_point[2]);
			for (int i = 0; i < neighbouring_triangles.size(); i++)
			{
				Point p(neighbouring_triangles[i].p[0], neighbouring_triangles[i].p[1], neighbouring_triangles[i].p[2]);
				Point q(neighbouring_triangles[i].q[0], neighbouring_triangles[i].q[1], neighbouring_triangles[i].q[2]);
				Point r(neighbouring_triangles[i].r[0], neighbouring_triangles[i].r[1], neighbouring_triangles[i].r[2]);
				Segment_3 cpp(np, p);
				Segment_3 cpq(np, q);
				Segment_3 cpr(np, r);
				if (!probfeatVal[neighbouring_triangles[i].pid].tl_visited && p != cp && cpp.squared_length() < minDis) {
					minDis = cpp.squared_length();
					np_id = neighbouring_triangles[i].pid;
					probfeatVal[np_id].tl_visited = true;
				}
				if (!probfeatVal[neighbouring_triangles[i].qid].tl_visited && q != cp && cpq.squared_length() < minDis) {
					minDis = cpq.squared_length();
					np_id = neighbouring_triangles[i].qid;
					probfeatVal[np_id].tl_visited = true;
				}
				if (!probfeatVal[neighbouring_triangles[i].rid].tl_visited && r != cp && cpr.squared_length() < minDis) {
					minDis = cpr.squared_length();
					np_id = neighbouring_triangles[i].rid;
					probfeatVal[np_id].tl_visited = true;
				}
			}
			if (np_id == -1)
				break;
			else {
				Line step;
				step.p[0] = current_point[0];
				step.p[1] = current_point[1];
				step.p[2] = current_point[2];
				step.q[0] = _inCloud->points[np_id].x;
				step.q[1] = _inCloud->points[np_id].y;
				step.q[2] = _inCloud->points[np_id].z;
				step.ndxq = np_id;
				tensor_line.push_back(step);
				current_length += step_size;
				current_point[0] = _inCloud->points[np_id].x;
				current_point[1] = _inCloud->points[np_id].y;
				current_point[2] = _inCloud->points[np_id].z;
				neighbouring_triangles.clear();
				neighbouring_triangles.insert(neighbouring_triangles.end(), probfeatVal[np_id].triangles.begin(), probfeatVal[np_id].triangles.end());
			}
		}
		probfeatVal[seed_point].tensor_line = tensor_line;
	}
}

float Dot(Vector_3 p, Vector_3 q)
{
	return p.x()*q.x() + p.y()*q.y() + p.z()*q.z();
}

// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
void Barycentric(Point p, Point a, Point b, Point c, float &u, float &v, float &w)
{
	Vector_3 v0(a, b), v1(a, c), v2(a, p);
	float d00 = Dot(v0, v0);
	float d01 = Dot(v0, v1);
	float d11 = Dot(v1, v1);
	float d20 = Dot(v2, v0);
	float d21 = Dot(v2, v1);
	float denom = d00 * d11 - d01 * d01;
	v = (d11 * d20 - d01 * d21) / denom;
	w = (d00 * d21 - d01 * d20) / denom;
	u = 1.0f - v - w;
}

int LineTriangleIntersection(Triangle_3 tri, Point cp, Point np, Point& intersectionP)
{
	Segment_3 seg(cp, np);
	CGAL::cpp11::result_of<K::Intersect_3(Segment_3, Triangle_3)>::type
		result = CGAL::intersection(seg, tri);
	if (result) {
		if (const Point* r = boost::get<Point>(&*result)) {
			if (*r == cp)
				return 0;
			intersectionP = *r;
			return 1;
		} /*else {
			intersectionP = q;
			return 1;
		}*/
	}
	return 0;
}

/*****************************************************************************************
 * Difference of normals (Beena's Code)
*****************************************************************************************/

void Utils::computeDONS(vector<probfeatnode>&probfeatVal, float radius)
{
	// Recreating arguments from Beena's code
	std::vector<int> neighborIndices;
	std::vector<int> neighborSize;

	std::vector<int> pointIdxRadiusSearch;
	std::vector<float> pointRadiusSquaredDistance;
	pcl::PointXYZ searchPoint;
	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		searchPoint = _inCloud->points[i];

		if (_octree->radiusSearch(searchPoint, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0)
		{
			neighborIndices.insert(neighborIndices.end(), pointIdxRadiusSearch.begin(), pointIdxRadiusSearch.end());

			if (i == 0)
				neighborSize.push_back(pointIdxRadiusSearch.size());
			else
				neighborSize.push_back(pointIdxRadiusSearch.size() + neighborSize[i - 1]);
			pointIdxRadiusSearch.clear();
			pointRadiusSquaredDistance.clear();
		}

	}

	pcl::search::Search<pcl::PointXYZ>::Ptr tree = boost::shared_ptr<pcl::search::Search<pcl::PointXYZ> >(new pcl::search::KdTree<pcl::PointXYZ>);

	pcl::PointCloud <pcl::Normal>::Ptr normals(new pcl::PointCloud <pcl::Normal>);
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normal_estimator;


	normal_estimator.setSearchMethod(tree);
	normal_estimator.setInputCloud(_inCloud);
	normal_estimator.setRadiusSearch(radius);
	normal_estimator.compute(*normals);

	size_t totCount, loc;
	double len;

	float max = -1.0;
	float min = 10000.0;
	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		if (i == 0)
		{
			totCount = neighborSize[0];
			loc = 0;
		}
		else
		{
			totCount = neighborSize[i] - neighborSize[i - 1];
			loc = neighborSize[i - 1];
		}

		double ncurv = 0.0;

		for (size_t j = 0; j < totCount; j++)
		{
			idxType indx = neighborIndices[loc + j];

			float n1 = sqrt((*normals)[i].normal_x * (*normals)[i].normal_x +
				(*normals)[i].normal_y * (*normals)[i].normal_y +
				(*normals)[i].normal_z * (*normals)[i].normal_z);

			float n2 = sqrt((*normals)[indx].normal_x * (*normals)[indx].normal_x +
				(*normals)[indx].normal_y * (*normals)[indx].normal_y +
				(*normals)[indx].normal_z * (*normals)[indx].normal_z);


			double nmag = (*normals)[i].normal_x * (*normals)[indx].normal_x +
				(*normals)[i].normal_y * (*normals)[indx].normal_y +
				(*normals)[i].normal_z * (*normals)[indx].normal_z;

			if (n1 != 0.0 && n2 != 0.0)
			{
				nmag = abs(nmag / (n1 * n2));
				if (nmag >= 1.0)
					nmag = 1.0;
				nmag = acos(nmag);
				ncurv += abs(nmag);
			}
		}

		if (totCount > 0)
		{
			len = totCount;
			probfeatVal[i].don = ncurv / len;
			if (probfeatVal[i].don > max)
				max = probfeatVal[i].don;
			if (probfeatVal[i].don < min)
				min = probfeatVal[i].don;
		}
	}
	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		probfeatVal[i].don = (probfeatVal[i].don - min) / (max - min);
	}
}

/*****************************************************************************************
 * Difference of normals using multi scale operator
*****************************************************************************************/


void Utils::computeDoNs_MSO(vector<probfeatnode>&probfeatVal, float rmin, float rmax)
{
	pcl::search::Search<pcl::PointXYZ>::Ptr tree = boost::shared_ptr<pcl::search::Search<pcl::PointXYZ> >(new pcl::search::KdTree<pcl::PointXYZ>);

	pcl::PointCloud <pcl::Normal>::Ptr normals1(new pcl::PointCloud <pcl::Normal>);
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normal_estimator1;

	pcl::PointCloud <pcl::Normal>::Ptr normals2(new pcl::PointCloud <pcl::Normal>);
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normal_estimator2;


	normal_estimator1.setSearchMethod(tree);
	normal_estimator1.setInputCloud(_inCloud);
	normal_estimator1.setRadiusSearch(rmin);
	normal_estimator1.compute(*normals1);

	normal_estimator2.setSearchMethod(tree);
	normal_estimator2.setInputCloud(_inCloud);
	normal_estimator2.setRadiusSearch(rmax);
	normal_estimator2.compute(*normals2);

	float max = -1.0;
	float min = 10000.0;
	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{

		float n1 = sqrt(((*normals1)[i].normal_x - (*normals2)[i].normal_x) * ((*normals1)[i].normal_x - (*normals2)[i].normal_x) +
			((*normals1)[i].normal_y - (*normals2)[i].normal_y) * ((*normals1)[i].normal_y - (*normals2)[i].normal_x) +
			((*normals1)[i].normal_z - (*normals2)[i].normal_z) * ((*normals1)[i].normal_z - (*normals2)[i].normal_x));

		probfeatVal[i].don = (n1) / 2;

		if (probfeatVal[i].don > max)
			max = probfeatVal[i].don;
		if (probfeatVal[i].don < min)
			min = probfeatVal[i].don;
	}
	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		probfeatVal[i].don = (probfeatVal[i].don - min) / (max - min);
	}
}

/*****************************************************************************************
 * Difference of normals using minor eigenvector
*****************************************************************************************/

void Utils::computeDoEs(vector<probfeatnode>&probfeatVal)
{
	float max = -1.0;
	float min = 10000.0;
	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		float minor_evec_rmin[] = { probfeatVal[i].evecs[3],probfeatVal[i].evecs[4],probfeatVal[i].evecs[5] };
		float minor_evec_rmax[] = { probfeatVal[i].evecs[6],probfeatVal[i].evecs[7],probfeatVal[i].evecs[8] };
		float orient[] = { -_inCloud->points[i].x, -_inCloud->points[i].y, 1 - _inCloud->points[i].z };
		int sign1 = (minor_evec_rmax[0] * orient[0] +
			minor_evec_rmax[1] * orient[1] +
			minor_evec_rmax[2] * orient[2] >= 0) ? 1 : -1;
		int sign2 = (minor_evec_rmin[0] * orient[0] +
			minor_evec_rmin[1] * orient[1] +
			minor_evec_rmin[2] * orient[2] >= 0) ? 1 : -1;

		normalize(minor_evec_rmax);
		normalize(minor_evec_rmin);

		minor_evec_rmin[0] *= sign2;
		minor_evec_rmin[1] *= sign2;
		minor_evec_rmin[2] *= sign2;
		minor_evec_rmax[0] *= sign1;
		minor_evec_rmax[1] *= sign1;
		minor_evec_rmax[2] *= sign1;

		float don[3];
		don[0] = minor_evec_rmin[0] - minor_evec_rmax[0];
		don[1] = minor_evec_rmin[1] - minor_evec_rmax[1];
		don[2] = minor_evec_rmin[2] - minor_evec_rmax[2];

		probfeatVal[i].don = sqrt(don[0] * don[0] + don[1] * don[1] + don[2] * don[2]);

		if (probfeatVal[i].don > max)
			max = probfeatVal[i].don;
		if (probfeatVal[i].don < min)
			min = probfeatVal[i].don;
	}
	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		probfeatVal[i].don = (probfeatVal[i].don - min) / (max - min);
	}
}

/*****************************************************************************************
 * Tensor line smoothening
*****************************************************************************************/

void Utils::simplifyTensorLines(vector<probfeatnode>&probfeatVal)
{
	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		if (probfeatVal[i].tensor_line.size() > 0)
		{
			std::vector<Point> pointList;
			for (int j = 0; j < probfeatVal[i].tensor_line.size(); j++)
				pointList.push_back(Point(probfeatVal[i].tensor_line[j].p[0], probfeatVal[i].tensor_line[j].p[1], probfeatVal[i].tensor_line[j].p[2]));
			std::vector<Point> resultSet = DouglasPeucker(pointList, DOUGLASPEUCKER_EPSILON, 0, pointList.size() - 1);
			probfeatVal[i].tensor_line.clear();
			for (int j = 0; j < resultSet.size() - 1; j++)
			{
				Line step;
				step.p[0] = resultSet[j].x();
				step.p[1] = resultSet[j].y();
				step.p[2] = resultSet[j].z();
				step.q[0] = resultSet[j + 1].x();
				step.q[1] = resultSet[j + 1].y();
				step.q[2] = resultSet[j + 1].z();
				probfeatVal[i].tensor_line.push_back(step);
			}
		}
	}
}

std::vector<Point> DouglasPeucker(std::vector<Point> pointList, float epsilon, int start, int end)
{
	// Find the point with the maximum distance
	float dmax = 0;
	int index = -1;
	for (int i = start + 1; i < end; i++) {
		float d = CGAL::squared_distance(pointList[i], Segment_3(pointList[start], pointList[end]));
		if (d > dmax && d != 0) {
			index = i;
			dmax = d;
		}
	}
	cout << "dmax = " << dmax << endl;
	// If max distance is greater than epsilon, recursively simplify
	std::vector<Point> resultList;
	if (dmax > epsilon*epsilon) {
		// Recursive call
		std::vector<Point> recResults1 = DouglasPeucker(pointList, epsilon, start, index);
		std::vector<Point> recResults2 = DouglasPeucker(pointList, epsilon, index, end);

		// Build the result list
		resultList.insert(resultList.end(), recResults1.begin(), recResults1.end() - 1);
		resultList.insert(resultList.end(), recResults2.begin(), recResults2.end());
	}
	else {
		resultList.push_back(pointList[start]);
		resultList.push_back(pointList[end]);
	}
	// Return the result
	return resultList;
}

void Utils::label_classification(vector<probfeatnode>&probfeatVal)
{
	int classification[3][9];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 9; j++) {
			classification[i][j] = 0;
		}
	}

	for (size_t j = 0; j < _inCloud->points.size(); j++)
	{
		if (probfeatVal[j].csclcp[0] >= probfeatVal[j].csclcp[1] && probfeatVal[j].csclcp[0] >= probfeatVal[j].csclcp[2]) {
			// PointType
			classification[0][probfeatVal[j].label]++;
		}
		else if (probfeatVal[j].csclcp[1] >= probfeatVal[j].csclcp[0] && probfeatVal[j].csclcp[1] >= probfeatVal[j].csclcp[2]) {
			// LineType
			classification[1][probfeatVal[j].label]++;
		}
		else {
			// SurfaceType
			classification[2][probfeatVal[j].label]++;
		}
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 9; j++) {
			std::cout << classification[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

void Utils::writeSalToFile(vector<probfeatnode>&probfeatVal)
{
	ofstream fout("clcscp.txt");
	fout << probfeatVal.size() << endl;
	for (int i = 0; i < probfeatVal.size(); i++) {
		fout << i << " " << probfeatVal[i].csclcp[1] << " " << probfeatVal[i].csclcp[2] << " " << probfeatVal[i].csclcp[0] << endl;
	}
	fout.close();
}

void Utils::writeDONToFile(vector<probfeatnode>&probfeatVal)
{
	ofstream fout("don.txt");
	fout << probfeatVal.size() << endl;
	for (int i = 0; i < probfeatVal.size(); i++) {
		fout << i << " " << probfeatVal[i].don << endl;
	}
	fout.close();
}

void Utils::writeNeighborsToFile(float radius)
{
	const char* filenames[6] = { "1_neighbours.txt","2_neighbours.txt","3_neighbours.txt","4_neighbours.txt","5_neighbours.txt" };
	for (int scale = 0; scale < 5; scale++)
	{
		ofstream fout(filenames[scale]);
		fout << _inCloud->points.size() << endl;

		std::vector<int> pointIdxRadiusSearch;
		std::vector<float> pointRadiusSquaredDistance;
		pcl::PointXYZ searchPoint;
		for (size_t i = 0; i < _inCloud->points.size(); i++)
		{
			searchPoint = _inCloud->points[i];
			fout << i << " ";

			if (_octree->radiusSearch(searchPoint, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0)
			{
				fout << pointIdxRadiusSearch.size() << " ";
				for (size_t m = 0; m < pointIdxRadiusSearch.size(); m++)
				{
					fout << pointIdxRadiusSearch[m] << " ";
				}
				pointIdxRadiusSearch.clear();
				pointRadiusSquaredDistance.clear();
			}
			fout << endl;
		}
		fout.close();
		radius += 0.001;
	}
}

void  Utils::writeTriangulationToFile(vector<probfeatnode>&probfeatVal) {
	cout << "Writing triangulation to file" << endl;
	ofstream fout("triangulation.txt");
	int count = 0;
	for (int i = 0; i < probfeatVal.size(); i++)
		count += probfeatVal[i].triangles.size();

	fout << count << endl;
	int ndx = 0;
	for (int i = 0; i < probfeatVal.size(); i++) {
		for (int j = 0; j < probfeatVal[i].triangles.size(); j++) {
			fout << ndx << " " << probfeatVal[i].triangles[j].pid << " " << probfeatVal[i].triangles[j].qid << " " << probfeatVal[i].triangles[j].rid << endl;
			ndx++;
		}
	}
	fout.close();
}

/*****************************************************************************************
 * Triangulation of entire point set
*****************************************************************************************/

void Utils::completeTriangulation()
{
	ofstream fout("complete_triangulation.txt");
	std::vector<Point> triangle_points;
	std::vector<int> triangle_points_ndx;
	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		triangle_points.push_back(Point(_inCloud->points[i].x, _inCloud->points[i].y, _inCloud->points[i].z));
		triangle_points_ndx.push_back(i);
	}
	/*Delaunay dt;
	dt.insert(triangle_points.begin(), triangle_points.end());
	int tndx = 0;
	for(Delaunay::Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); ++it)
	{
		Point p1 = it->vertex(0)->point();
		Point p2 = it->vertex(1)->point();
		Point p3 = it->vertex(2)->point();
		int i1=-1;
		int i2=-1;
		int i3=-1;
		for(int k=0;k<triangle_points.size();k++){
			if(p1==triangle_points[k]){
				i1 = triangle_points_ndx[k];
			}
			if(p2==triangle_points[k]){
				i2 = triangle_points_ndx[k];
			}
			if(p3==triangle_points[k]){
				i3 = triangle_points_ndx[k];
			}
		}

		fout << tndx  << " " << i1 << " " << i2  << " " << i3 << endl;
		tndx++;
	}*/
	fout.close();
}

void Utils::correctMisclassifiedPoints(ftType eigenvalues[3], ftType evecs[9], float Lp, float Ll, tensorType &out_tensor)
{
	float factor = (1 + 1.5 * (Lp / Ll)) * eigenvalues[2];

	out_tensor.evec0[0] = factor * (evecs[0] * evecs[0] + evecs[3] * evecs[3]);
	out_tensor.evec0[1] = factor * (evecs[0] * evecs[1] + evecs[3] * evecs[4]);
	out_tensor.evec0[2] = factor * (evecs[0] * evecs[2] + evecs[3] * evecs[5]);

	out_tensor.evec1[0] = factor * (evecs[0] * evecs[1] + evecs[3] * evecs[4]);
	out_tensor.evec1[1] = factor * (evecs[1] * evecs[1] + evecs[4] * evecs[4]);
	out_tensor.evec1[2] = factor * (evecs[1] * evecs[2] + evecs[4] * evecs[5]);

	out_tensor.evec2[0] = factor * (evecs[0] * evecs[2] + evecs[3] * evecs[5]);
	out_tensor.evec2[1] = factor * (evecs[1] * evecs[2] + evecs[4] * evecs[5]);
	out_tensor.evec2[2] = factor * (evecs[2] * evecs[2] + evecs[5] * evecs[5]);

	out_tensor.evec0[0] += evecs[6] * evecs[6];
	out_tensor.evec0[1] += evecs[6] * evecs[7];
	out_tensor.evec0[2] += evecs[6] * evecs[8];

	out_tensor.evec1[0] += evecs[6] * evecs[7];
	out_tensor.evec1[1] += evecs[7] * evecs[7];
	out_tensor.evec1[2] += evecs[7] * evecs[8];

	out_tensor.evec2[0] += evecs[6] * evecs[8];
	out_tensor.evec2[1] += evecs[7] * evecs[8];
	out_tensor.evec2[2] += evecs[8] * evecs[8];

	out_tensor.evec0[0] *= eigenvalues[2];
	out_tensor.evec0[1] *= eigenvalues[2];
	out_tensor.evec0[2] *= eigenvalues[2];

	out_tensor.evec1[0] *= eigenvalues[2];
	out_tensor.evec1[1] *= eigenvalues[2];
	out_tensor.evec1[2] *= eigenvalues[2];

	out_tensor.evec2[0] *= eigenvalues[2];
	out_tensor.evec2[1] *= eigenvalues[2];
	out_tensor.evec2[2] *= eigenvalues[2];
}

void Utils::eigen_decomposition_to_file(vector <tensorType> &aggregated_tensors)
{
	cout << "writing to eigen_decomposition.txt" << endl;
	ofstream fout("eigen_decomposition.txt");
	fout << aggregated_tensors.size() << endl;
	for (size_t idx = 0; idx < aggregated_tensors.size(); idx++)
	{
		float A[3][3], V[3][3], d[3];
		for (int i = 0; i < 3; i++)
		{
			A[i][0] = (aggregated_tensors)[idx].evec0[i];
			A[i][1] = (aggregated_tensors)[idx].evec1[i];
			A[i][2] = (aggregated_tensors)[idx].evec2[i];
		}

		eigen_decomposition(A, V, d); //d[2] > d[1] > d[0]

		fout << d[2] << " " << d[1] << " " << d[0] << " "
			<< V[0][2] << " " << V[1][2] << " " << V[2][2] << " "
			<< V[0][1] << " " << V[1][1] << " " << V[2][1] << " "
			<< V[0][0] << " " << V[1][0] << " " << V[2][0] << endl;

	}
	fout.close();
}