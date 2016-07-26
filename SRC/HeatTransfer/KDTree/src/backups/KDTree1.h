#ifndef _KDTREE_H_
#define _KDTREE_H_

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

  // Search for the nearest neighbor to pnt and 
  // return its index
  //int closest_point(double *pnt, int &idx, bool approx=false);

  int get_points_in_range(Range *range);

  int create(double *setpoints, int setnpoints, int setndim,
	  bool setCopy = false, struct _Node *setNodeMem = (struct _Node *)0);

  // Return distance squared between points
  // between points
  //inline double distsq(double *pnt1, double *pnt2);

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

	//int check_border_distance(int nodeIdx, int dim,
	//	double *pnt, double &dist, int &idx);

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
};

#endif
