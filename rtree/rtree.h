/* rtree.h
   this file defines the class RTree*/

#ifndef __RTREE
#define __RTREE
//------------------------------------------------------------
#include "../func/gendef.h"
#include "../heap/heap.h"
//Added by Sarowar
#include "../rtree_time.h"
//-----------------------------
//Added by Tanzima
//#include "../global.h"
//...............

//Added by Tanzima
#include <vector>
using namespace std;


//added for TP KNN------------------------------------------------
#define SEQUENCE_SENSITIVE false
  //set it true if you want a new partition point every time the order of the NNs change
#define PAST_PAIR 1000

//Added for kGNN by Tanzima


const double MAXDOUBLE = numeric_limits<double>::infinity();
const int kMAX=100; // changed by Eunus : it was set to 500
const int MAXDATALIMIT=20000;
//const int MAXGROUP = 128;
//const int MINSUBGROUP = MAXGROUP/2;
//------------------------------------------------------------
class LinList;
class SortedLinList;
class Cache;
class RTNode;
class Entry;
//Added by Tanzima
struct Rectangle1;
struct DistfromPOI;
//------------------------------------------------------------
class RTree : public Cacheable
{
public:
//--===on disk===--
	int dimension;                       
	int num_of_data;	                 
    int num_of_dnodes;	                 
    int num_of_inodes;	                 
	int root;                            
	bool root_is_data;                   
//--===others===--
	RTNode *root_ptr;
    bool *re_level;  
    LinList *re_data_cands; 
	LinList *deletelist;

//--===Sarowar : KEYWORDS : FRUITS===--
	int bitmask_array_length;
	string fruits[128][20];
	vector<string> vector_string;
	//debug -- print
	int visit_order;
	int leaf_count;
	int negative_bitmask_count;
	int query_count;
	int internal_node;
	//---------------

//--===added for TP KNN===--
	int last_pair[PAST_PAIR][2]; //records the pairs of points that produce the minimum inflence time
	int lastcnt; //next last pair to be replaced
	Heap *tpheap;

//--===functions===--
	RTree(char *fname, int _b_length, Cache* c, int _dimension);
    RTree(char *fname, Cache* c);
    RTree(char *inpname, char *fname, int _blength, Cache* c, int _dimension);
	RTree(char *inpname, char *fname,char *keywords, int _blength, Cache* c, int _dimension);
    ~RTree();
	void del_root();
	bool delete_entry(Entry *d);
	bool FindLeaf(Entry *e);
    int get_num() { return num_of_data; }
	void insert(Entry *d);
	void load_root();  
	void rangeQuery(float *mbr, SortedLinList *res);
	//	Sarowar : debug
	void rangeQuery_with_time(float *mbr,RTIME *mbt, SortedLinList *res);
	//--------------------------------------------------------------------------
	void read_header(char *buffer);      
	void write_header(char *buffer);
	int update_rslt(Entry *_e, float _dist, Entry *_rslt, 
					 float *_key, int _k);

	// This function was added to perform TP-kNN queries by Bobby
	void TPNN_TP(float *_qline, int _k, Entry *_nn, Entry *_rslt, float _max_trvl);
	
	//--added for valdity region queries---
	void rect_win_query(float *mbr, LinList *in_objs, LinList *out_objs_so_far);
	void rect_win_query(float *mbr, float *exclmbr, LinList *c_inf_objs);
	void BFNN(float *_qpt, int _k, Entry *_rslt);
	//Nasim
	void BFNN_with_keywords_time_uncertainty(float *_qpt, RTIME *_qtime,int *keywords, char *BFNN, int _k, int threshold, Entry *_rslt, float *probability, int *total_result_found);
	//	Sarowar : debug
	void BFNN_with_time(float *_qpt,RTIME *_qtime, int _k, Entry *_rslt);
	void BFNN_with_time_keywords(float *_qpt,RTIME *_qtime,int *keywords, char *BFNN, int _k, Entry *_rslt);
	int load_keywords_from_file(char *keywords_file);
	void print_rtree_nodes(char *file);
	void dfs_visit_print(FILE *fp, int son);
	void print_node_during_insert();
	//--------------------------------------------------------------------------
	void BFNNCont(float *qmbr, float *qmbr2, Heap *heap, HeapEntry *e, int k);
	//Added by Tanzima
	float KMin(Point2D m, int c, float cl, int k, Pointlocation _rslt[], int *num_of_data,float d_safe_corner,DistfromPOI _cornertoPOI[]);
	void UpdateCount(Rectangle1 R, float cl, int k, Pointlocation _rslt[], int *num_of_data, float r, int *count, DistfromPOI _cornertoPOI[]);
	int UpdateStatus(Rectangle1 R, float cl, int k, Pointlocation _rslt[], int *num_of_data, int *count, DistfromPOI _cornertoPOI[]);
	void Rect_kNNQ(Rectangle1 R, float cl, int k, Pointlocation *_rslt, int *num_of_data);
	void Point_BFN_NNQ(Point2D o, double *_rslt);

	//Added by Tanzima for kGNN
	float RTree::KMax(int k, Pointlocation _rslt[], int *num_of_data);
	void private_kGNN_sum(Rectangle1 R[], int g_size, int k, Pointlocation _rslt[], int *num_of_data);
	void private_kGNN_max(Rectangle1 R[], int g_size, int k, Pointlocation rslt[], int *num_of_data);

	//Added by Eunus for Consensus
	void RTree::kGNN(Rectangle1 R[], int g_size, int k, int f, Pointlocation rslt[], int *num_of_data);
	void RTree::Naive_Consensus_kGNN(Rectangle1 R[], int g_size, int sg_size, int k, int f, SG L[][16], int *num_of_data);
	void RTree::consensus_kGNN(Rectangle1 R[], int g_size, int sg_size, int k, int f, SG L[][16], int *num_of_data, int tb_lb);
	void RTree::MQ_consensus_kGNN(Rectangle1 R[], int g_size, int sg_size, int k, int f, SG L[][16], int *num_of_data);
	void RTree::FANN_consensus_kGNN(Rectangle1 R[], int g_size, int sg_size, int k, int f, SG L[][16], int lind, int *num_of_data, int tb_lb); //FANN each category at a time

};

#endif // __RTREE
