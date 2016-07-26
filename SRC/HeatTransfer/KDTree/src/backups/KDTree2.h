/*Copyright (c) 2005, Steven Michael@MIT 
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met:

    * Redistributions of source code must retain the above copyright 
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright 
      notice, this list of conditions and the following disclaimer in 
      the documentation and/or other materials provided with the distribution
      
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.
*/

//Date:March 2012
//Modified:Yaqiang Jiang @ The University of Edinburgh

#ifndef _KDTREE_H_
#define _KDTREE_H_
#include <vector>

#pragma pack(1)
struct _Node {
	double key;
	int pntidx;
	int leftIdx;
	int rightIdx;
	};
#pragma pack()


typedef double *Range;
typedef int IRange[2];


class KDTree {

public:
  KDTree(double *setpoints, int N, int setndim);
  KDTree();
  virtual ~ KDTree();

  const std::vector<int>& get_points_in_range(Range *range);

  int create(double *setpoints, int setnpoints, int setndim,
	  bool setCopy = false, struct _Node *setNodeMem = (struct _Node *)0);

  // The following functions allow all the information in the class
  // to be serialized and unserialized.  This is convenient, for example,
  // for writing the tree to a disk or to a MATLAB variable
  //static int get_serialize_length(int npoints,int ndim);
  //static KDTree *unserialize(void *mem);

  //int set_verbosity(int v){verbosity=v;return 0;}

private:

	//int create(double *setpoints,int setnpoints,int setndim,
	//	void *mem);

	// Do we copy the points or not
	// if we do, this is where they go
	bool copyPoints;

	int range_search(int nodeIdx, Range *range, int dim);
	int heapsort(int dim, int *idx, int len);

	int build_kdtree(int **sortidx, int dim, int *pidx, int len);

	int *workArr;
	struct _Node *nodeMem;
	int nodeMemCnt;
	bool nodeMemAlloc;
	struct _Node *node_alloc();

	int *intArrMem;
	int  intArrMemCnt;
	int *int_alloc(int len);

	static int (*logmsg)(const char *,...);
	int verbosity;

	int nPntsInRange;
	int *pntsInRange;
	int ndim;
	int npoints;
	double *points;
	std::vector<int>* PTSInRange;
};

#endif
