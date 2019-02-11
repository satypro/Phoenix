#pragma once
#ifndef DATATYPE_COMMON_H
#define DATATYPE_COMMON_H

#include <stdio.h>
#include <vector>

using namespace std;

typedef float ftType;
typedef int idxType;

const double MIN_X = -1.0;
const double MIN_Y = -1.0;
const double MIN_Z = -1.0;
const double MAX_X = 1.0;
const double MAX_Y = 1.0;
const double MAX_Z = 1.0;
const double DEL_X = MAX_X - MIN_X;
const double DEL_Y = MAX_Y - MIN_Y;
const double DEL_Z = MAX_Z - MIN_Z;


struct SLasHeaderBlock
{
	unsigned char ucVerMaj;   //1 byte
	unsigned char ucVerMin;  //1 byte
	unsigned char ucFormatID;  //1 byte
	unsigned short usnRecLength;  //2byte
	unsigned int usnNumOfPoint;    //4 byte	
	unsigned int usnNumOfPointReturn[5];    //20 byte	
	unsigned int ulnOffset;   //4 byte
	double dXScale;
	double dYScale;
	double dZScale;
	double dXOffset;
	double dYOffset;
	double dZOffset;

	double xmax;
	double xmin;

	double ymax;
	double ymin;

	double zmax;
	double zmin;
};

struct Line
{
	float p[3];
	float q[3];
	int ndxp;
	int ndxq;
};

struct Rect
{
	float p[3];
	float q[3];
	float r[3];
	float s[3];
	int ndxp;
	int ndxq;
	int ndxr;
	int ndxs;
};

struct Triangle
{
	float p[3];
	float q[3];
	float r[3];
	float average_cl;
	bool hasCL;
	float CL_p1[3]; //Isoline p
	float CL_p2[3]; //Isoline q

	float average_anisotropy;

	int tlEdge;
	float tlScore;

	int pid;
	int qid;
	int rid;
};

struct probfeatnode
{
	ftType prob[3]; //prob 0 = spher, 1 = curve, 2 = disc
	ftType featStrength[3]; // 0 =spher, 1 = curve, 2 = disc
	ftType csclcp[3];
	ftType eigenvalues[3];

	/*
	 * evec[0,1,2] = major evec at r_max
	 * evec[3,4,5] = minor evec at r_min
	 * evec[6,7,8] = minor evec at r_max
	 *
	 * ** USED this array for storing r_min minor evec and not separately because of struct's memory limitation
	 * */
	ftType evecs[9];
	ftType planarity;
	ftType anisotropy;
	ftType sphericity;
	ftType sum_eigen;
	ftType linearity;
	ftType omnivariance;
	ftType eigenentropy;
	double don;

	int label;
	int optimalScale;
	int numScale;

	std::vector<Triangle>  triangles;
	std::vector<Line>  tensor_line;
	bool tl_visited = false;
};

struct featProps
{
	int type; // 0 - line, 1 - surface, 2   - surface and line - transition points
	float n_ps; // point type density  - 4
	float n_ls; // line type density - 0 
	float n_surfs; //surface type density - 1 
	float n_pls; // line and  point density 0 and 4  
	float n_psusfs;  // surface and  surface density 0 and 4
};

struct glyphVars
{
	float evals[3];
	float evecs[9];
	float uv[2];
	float abc[3];
	float csclcp[3];
};

struct tensorType
{
	float evec0[3];
	float evec1[3];
	float evec2[3];
};

struct metaVelData
{
	float vel[3];
	int index[3];
};

struct tempminval
{
	float mval;
	idxType idx;
};

struct node
{
	idxType idx;
	bool status;
};


struct nodeInfo
{
	int idx; // corresponding to  original point cloud
	bool type; // true = mot critical point, false = critical point
	int seedIdx;
};

struct gnode
{
	nodeInfo nd;
	std::vector<nodeInfo> edge;
};

struct myfloat3
{
	float x;
	float y;
	float z;
};
#endif