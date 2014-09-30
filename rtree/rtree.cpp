/*rtree.cpp
  this file implements the RTree class*/
/*
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
*/
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include "rtree.h"
#include "entry.h"
#include "rtnode.h"
#include "distance.h"
#include "../blockfile/cache.h"
#include "../blockfile/blk_file.h"
#include "../linlist/linlist.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <istream>
#include <iostream>


//Added by Tanzima
extern int io_access;
extern int disktime;
extern int updatetime;
extern int counttime;
extern int kmintime;

extern int ind_comb;
//......................
//------------------------------------------------------------
RTree::RTree(char *fname, int _b_length, Cache *c, int _dimension)
  //use this constructor to build a new tree
{
    file = new BlockFile(fname, _b_length);
    cache = c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    dimension = _dimension;
    root = 0;
    root_ptr = NULL;
    root_is_data = TRUE;
    num_of_data = num_of_inodes = num_of_dnodes = 0;

    root_ptr = new RTNode(this);
	  //note that when a tree is constructed, the root is automatically created
	  //though at this time there is no entry in the root yet.
    num_of_dnodes++;
    root_ptr -> level = 0;
    root = root_ptr -> block;

	//added for TP KNN----------------------------------------
	tpheap=NULL;
}

//------------------------------------------------------------
RTree::RTree(char *fname, Cache *c)
  //use this constructor to restore a tree from a file
{
    file = new BlockFile(fname, 0);
    cache =c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    char *header = new char [file->get_blocklength()];
    file -> read_header(header);
    read_header(header);
	delete [] header;

    root_ptr = NULL;

	//added for TP KNN----------------------------------------
	tpheap=NULL;
}
//------------------------------------------------------------
RTree::RTree(char *inpname, char *fname, int _b_length, Cache *c, int _dimension)
  // construct new R-tree from a specified input textfile with rectangles
{		/*	Sarowar debug : RTree constructor */
	printf("\n\nSarowar debug : In the Constructor with the *inpname, *fname.\n");
	Rtree_time rtime;
	time_t t1,t2;
	bitmask_array_length = (int) (TOTAL_FRUIT_NUMBER / 32);
	//-----------------------------------------
    Entry *d;
    FILE *fp;
	file = new BlockFile(fname, _b_length);
    cache =c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    char *header = new char [file->get_blocklength()];
    read_header(header);// reads header from this char * header.
	delete [] header;

    dimension = _dimension;
    root = 0;
    root_ptr = NULL;
    root_is_data = TRUE;
    num_of_data = num_of_inodes = num_of_dnodes = 0;

    root_ptr = new RTNode(this);
    num_of_dnodes++;
    root_ptr -> level = 0;
    root = root_ptr->block;

	//added for TP KNN----------------------------------------
	tpheap=NULL;

	int record_count = 0;

    if((fp = fopen(inpname,"r")) == NULL)
    {
      delete this;
      error("Cannot open R-Tree text file", TRUE);
    }
    else
    {
      // Sarowar : debug
	  int ij=1;
	  int x,y,a;
	  //---------------
      while (!feof(fp))
      {
		record_count ++;

		d = new Entry(dimension, NULL);

    	fscanf(fp, "%d", &(d -> son));
		
		fscanf(fp, " %f %f %f %f", &(d->bounces[0]),&(d->bounces[1]), &(d->bounces[2]), &(d->bounces[3]));
		
    	/*
			Sarowar debug : RTree constructor
		*/
		/*
			Sarowar time : start : 9

			printf("time add : Linlist : linkable : clone :9\n");
		*/
			//fscanf(fp, " %lf %lf", &(d->times[0]),&(d->times[1])	);
			rtime.time_input_from_file(fp, &(d->times[0]), &(d->times[1]) );
			
			d->fruit_bitmask[0]=d->fruit_bitmask[1]=d->fruit_bitmask[2]=d->fruit_bitmask[3]=0;
			y=4;
			while(y--){
				fscanf(fp,"%d",&a);
				x = a/32 ;
				a = a%32 ;
				d->fruit_bitmask[x] |= ( 1<<a );
			}	
			//fscanf(fp, " %d %d %d %d ", &(d -> fruit_bitmask[0]), &(d->fruit_bitmask[1]), &(d->fruit_bitmask[2]), &(d->fruit_bitmask[3]) );
		/*
			Sarowar time : end : 9
		*/

		printf("\n\n\n%d'th data : Sarowar debug : rtree constructor : before insert() : \n d->son= %d \n d->bounces[0]= %f \n d->bounces[1]= %f \n d->bounces[2]= %f \n d->bounces[3]= %f \n d->times[0] = %lld \n d->times[1] = %lld \nd->fruit_bitmask[0]=%d d->fruit_bitmask[1]= %d \n d->fruit_bitmask[2] =%d d->fruit_bitmask[3] =%d \n\n",ij++, (d -> son), (d->bounces[0]), (d->bounces[1]), (d->bounces[2]), (d->bounces[3]),(d->times[0]),(d->times[1]), (d->fruit_bitmask[0]), (d->fruit_bitmask[1]), (d->fruit_bitmask[2]), (d->fruit_bitmask[3]) );
		printf("//getchar() : Press Enter 2 Continue.\n");
		//getchar();
		insert(d);

		  //d will be deleted in insert()

		if (record_count % 100 == 0)
		{
			for (int i = 0; i < 79; i ++)  //clear a line
				printf("\b");
			// prints for : 50, 5000
			printf("inserting object %d", record_count);
		}
      }
    }
	//TANZIMA
	// prints for : 5000
	printf("inserting object %d", record_count);
	fclose(fp);

	printf("\n");
	delete root_ptr;
	root_ptr = NULL;
}
//------------------------------------------------------------
//---------------------------------------------------------------------------------------
/*--------------------------- Sarowar : Constructor : Keywords ------------------------*/
RTree::RTree(char *inpname, char *fname,char *keywords, int _blength, Cache* c, int _dimension)
  // construct new R-tree from a specified input textfile with rectangles
{		/*	Sarowar debug : RTree constructor */
	printf("\n\nSarowar debug : In the Constructor with the *keywords, *inpname, *fname.\n");
	Rtree_time rtime;
	time_t t1,t2;
	bitmask_array_length = (int) (TOTAL_FRUIT_NUMBER / 32);
	//-----------------------------------------
    Entry *d;
    FILE *fp;
	file = new BlockFile(fname, _blength);
    cache =c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    char *header = new char [file->get_blocklength()];
    read_header(header);// reads header from this char * header.
	delete [] header;

    dimension = _dimension;
    root = 0;
    root_ptr = NULL;
    root_is_data = TRUE;
    num_of_data = num_of_inodes = num_of_dnodes = 0;

    root_ptr = new RTNode(this);
    num_of_dnodes++;
    root_ptr -> level = 0;
    root = root_ptr->block;

	//added for TP KNN----------------------------------------
	tpheap=NULL;

	int record_count = 0;

    if((fp = fopen(inpname,"r")) == NULL)
    {
      delete this;
      error("Cannot open R-Tree text file", TRUE);
    }
    else
    {
      // Sarowar : debug
	  int ij=1;
	  int x,y,a;
	  //---------------
      while (!feof(fp))
      {
		record_count ++;

		d = new Entry(dimension, NULL);

    	fscanf(fp, "%d", &(d -> son));
		
		fscanf(fp, " %f %f %f %f", &(d->bounces[0]),&(d->bounces[1]), &(d->bounces[2]), &(d->bounces[3]));
		/*	Sarowar debug : RTree constructor	*/
		/*	Sarowar time : start : 9
			printf("time add : Linlist : linkable : clone :9\n");	*/
			//fscanf(fp, " %lf %lf", &(d->times[0]),&(d->times[1])	);
			rtime.time_input_from_file(fp, &(d->times[0]), &(d->times[1]) );
			
			d->fruit_bitmask[0]=d->fruit_bitmask[1]=d->fruit_bitmask[2]=d->fruit_bitmask[3]=0;
			y=4;
			while(y--){
				fscanf(fp,"%d",&a);
				x = a/32 ;
				a = a%32 ;
				d->fruit_bitmask[x] |= ( 1<<a );
			}
			//fscanf(fp, " %d %d %d %d ", &(d -> fruit_bitmask[0]), &(d->fruit_bitmask[1]), &(d->fruit_bitmask[2]), &(d->fruit_bitmask[3]) );
		/*	int a,b,c,d;
			fscanf(fp," %d %d %d %d ", &a,&b,&c,&d );
			------------------------------
			int j=4;
			while(j--){
				fscanf(fp," %d ",&a);
				i = a/32 ;
				a = a%32 ;
				d->fruit_bitmask[i] |= ( 1<<a );
			}
			------------------------------
			Sarowar time : end : 9
		*/
			

		printf("\n\n\n%d'th data : Sarowar debug : rtree constructor : before insert() : \n d->son= %d \n d->bounces[0]= %f \n d->bounces[1]= %f \n d->bounces[2]= %f \n d->bounces[3]= %f \n d->times[0] = %lld \n d->times[1] = %lld \nd->fruit_bitmask[0]=%d d->fruit_bitmask[1]= %d \n d->fruit_bitmask[2] =%d d->fruit_bitmask[3] =%d \n\n",ij++, (d -> son), (d->bounces[0]), (d->bounces[1]), (d->bounces[2]), (d->bounces[3]),(d->times[0]),(d->times[1]), (d->fruit_bitmask[0]), (d->fruit_bitmask[1]), (d->fruit_bitmask[2]), (d->fruit_bitmask[3]) );
		printf("//getchar() : Press Enter 2 Continue.\n");
		//getchar();
		insert(d);

		  //d will be deleted in insert()

		if (record_count % 100 == 0)
		{
			for (int i = 0; i < 79; i ++)  //clear a line
				printf("\b");
			// prints for : 50, 5000
			printf("inserting object %d", record_count);
		}
      }
    }
	//TANZIMA
	// prints for : 5000
	printf("inserting object %d", record_count);
	fclose(fp);

	printf("\n");
	delete root_ptr;
	root_ptr = NULL;
}
//---------------------------------------------------------------------------------------
RTree::~RTree()
{
	char *header = new char[file -> get_blocklength()];
    write_header(header);
    file->set_header(header);
    delete [] header;

    if (root_ptr != NULL)
    {
        delete root_ptr;
        root_ptr = NULL;
    }

	if (cache)
      cache -> flush();

    delete file;

    delete re_data_cands;
	delete deletelist;

    //printf("This R-Tree contains %d internal, %d data nodes and %d data\n",
	//   num_of_inodes, num_of_dnodes, num_of_data);
}
//------------------------------------------------------------
void RTree::del_root()
{
	printf("\n root is deleted.\n");
	delete root_ptr;
	root_ptr = NULL;
}
//------------------------------------------------------------
bool RTree::delete_entry(Entry *d)
{
	load_root();

	R_DELETE del_ret;
	del_ret=root_ptr->delete_entry(d);

	if (del_ret == NOTFOUND) return false;
	if (del_ret == ERASED)
		error("RTree::delete_entry--The root has been deleted\n",true);

	if (root_ptr -> level > 0 && root_ptr -> num_entries == 1)
		//there is only one entry in the root but the root
		//is not leaf.  in this case, the child of the root is exhalted to root
	{
		root = root_ptr -> entries[0].son;
		delete root_ptr;
		root_ptr = NULL;
		load_root();
		num_of_inodes--;
	}

	//Now will reinsert the entries
	while (deletelist -> get_num() > 0)
	{
		Linkable *e;
		e = deletelist -> get_first();
		Entry *new_e = new Entry(dimension, NULL);
		new_e -> set_from_Linkable(e);
		deletelist -> erase();
		insert(new_e);
	}

	delete root_ptr;
	root_ptr = NULL;

	return true;
}
//------------------------------------------------------------
void RTree::insert(Entry* d)
{
    int i, j;
    RTNode *sn;
    RTNode *nroot_ptr;
    int nroot;
    Entry *de;
    R_OVERFLOW split_root;
    Entry *dc;
    float *nmbr;
	// Sarowar : added new field for time. used when node is to be splited.
	//
	RTIME *new_minimum_bounding_time;
	int	  *new_minimum_bitmask;
    //----------------------------------------------------------------------
    // load root into memory
    load_root();

    // no overflow occured until now
    re_level = new bool[root_ptr -> level + 1];
    for (i = 0; i <= root_ptr -> level; i++)
        re_level[i] = FALSE;

    // insert d into re_data_cands as the first entry to insert
    // make a copy of d because it should be erased later
    Linkable *new_link;
	new_link = d -> gen_Linkable();
	re_data_cands -> insert(new_link);

	delete d;  //we follow the convention that the entry will be deleted when insertion finishes

    j = -1;
    while (re_data_cands -> get_num() > 0)
    {
        // first try to insert data, then directory entries
	    Linkable *d_cand;
		d_cand = re_data_cands -> get_first();
        if (d_cand != NULL)
        {
            // since "erase" deletes the data itself from the
            // list, we should make a copy of the data before
            // erasing it
			dc = new Entry(dimension, NULL);
            dc -> set_from_Linkable(d_cand);
            re_data_cands -> erase();

            // start recursive insert with root

			/*
				Sarowar debug :

			*/
			printf("\n\nSarowar debug : root: insert(): before :split_root = root_ptr -> insert(dc, &sn);\n");
			printf("//getchar() : Press Enter 2 Continue.\n");
			//getchar();

			split_root = root_ptr -> insert(dc, &sn);
        }
        else
	        error("RTree::insert: inconsistent list re_data_cands", TRUE);

    	if (split_root == SPLIT)
    	// insert has lead to split --> new root-page with two sons (i.e. root and sn)
    	{
			/*
				Sarowar debug
			*/
			printf("\n\nSarowar debug : if (split_root == SPLIT)=true \n");
			printf("//getchar() : Press Enter 2 Continue.\n");
			//getchar();

    	    nroot_ptr = new RTNode(this);
    	    nroot_ptr -> level = root_ptr -> level + 1;
    	    num_of_inodes++;
    	    nroot = nroot_ptr -> block;

			/*
				Sarowar debug
			*/
			printf("\n\nSarowar debug : \n nroot_ptr -> level=%d \n num_of_inodes=%d \n nroot's block=%d \n",nroot_ptr -> level, num_of_inodes, nroot);
			printf("//getchar() : Press Enter 2 Continue.\n");
			//getchar();


			/*
				Sarowar debug : 1 partition : the root_ptr
			*/
    	    de = new Entry(dimension, this);
    	    nmbr = root_ptr -> get_mbr();
			memcpy(de->bounces, nmbr, 2*dimension*sizeof(float));
			delete [] nmbr;

			// Sarowar : time : returning minimum_bounding_time of root node.
			new_minimum_bounding_time = root_ptr -> get_minimumboundingtime();
			memcpy( de->times, new_minimum_bounding_time, RTIME_SIZE * sizeof(RTIME) );
			delete [] new_minimum_bounding_time;

			new_minimum_bitmask       = root_ptr -> get_minimumFruitBitmask();
			memcpy( de->fruit_bitmask, new_minimum_bitmask, bitmask_array_length * sizeof(int));
			delete [] new_minimum_bitmask;
			// Sarowar


    	    de->son = root_ptr->block;
    	    de->son_ptr = root_ptr;
    	    nroot_ptr -> enter(de);

			/*
				Sarowar debug : 2nd partition
			*/

    	    de = new Entry(dimension, this);
    	    nmbr = sn -> get_mbr();
    	    memcpy(de -> bounces, nmbr, 2*dimension*sizeof(float));
    	    delete [] nmbr;
    	    de -> son = sn -> block;
    	    de -> son_ptr = sn;
			// Sarowar : time : returning minimum_bounding_time of sn node the 2nd partition.
			new_minimum_bounding_time = sn -> get_minimumboundingtime();//root_ptr -> get_minimuboundingtime();
			memcpy( de->times, new_minimum_bounding_time, RTIME_SIZE * sizeof(RTIME) );
			delete [] new_minimum_bounding_time;

			new_minimum_bitmask       = sn -> get_minimumFruitBitmask();
			memcpy( de->fruit_bitmask, new_minimum_bitmask, bitmask_array_length * sizeof(int));
			delete [] new_minimum_bitmask;
			// Sarowar

    	    nroot_ptr->enter(de);

			/*
				Sarowar debug : root assingment.
			*/
    	    root = nroot;
            root_ptr = nroot_ptr;

            root_is_data = FALSE;
        }
        j++;
		/*
			Sarowar debug:
		*/
		printf("\n\nwhile (re_data_cands -> get_num()=%d > 0) : loop count var j=%d \n",re_data_cands -> get_num(),j);
    }

    num_of_data++;

    delete [] re_level;

	delete root_ptr;
	root_ptr = NULL;
}
//------------------------------------------------------------
void RTree::load_root()
{
	if (root_ptr == NULL){
        root_ptr = new RTNode(this, root);
		printf("\nload_root() : (root_ptr == NULL)\n-------------\nroot=%d",root);
	}
	else
		printf("\nload_root() : (root_ptr != NULL)\n-------------\nroot=%d",root);
}
//------------------------------------------------------------
void RTree::rangeQuery(float *mbr, SortedLinList *res)
{
    load_root();
	//Added by Tanzima
	io_access++;
	//..............
    root_ptr -> rangeQuery(mbr,res);

	delete root_ptr;
	root_ptr = NULL;
}
//------------------------------------------------------------
//	Sarowar : debug
void RTree::rangeQuery_with_time(float *mbr,RTIME *mbt, SortedLinList *res)
{
	load_root();
	//Added by Tanzima
	io_access++;
	//..............
	root_ptr -> rangeQuery_with_time(mbr, mbt, res); //rangeQuery(mbr,res);

	delete root_ptr;
	root_ptr = NULL;
}
//--------------------------------------------------------------------------
//-------------------------------------------------------------
void RTree::read_header(char *buffer)
{
    int i;

    memcpy(&dimension, buffer, sizeof(dimension));
    i = sizeof(dimension);

    memcpy(&num_of_data, &buffer[i], sizeof(num_of_data));
    i += sizeof(num_of_data);

    memcpy(&num_of_dnodes, &buffer[i], sizeof(num_of_dnodes));
    i += sizeof(num_of_dnodes);

    memcpy(&num_of_inodes, &buffer[i], sizeof(num_of_inodes));
    i += sizeof(num_of_inodes);

    memcpy(&root_is_data, &buffer[i], sizeof(root_is_data));
    i += sizeof(root_is_data);

    memcpy(&root, &buffer[i], sizeof(root));
    i += sizeof(root);
}
//------------------------------------------------------------
int RTree::update_rslt(Entry *_e, float _dist, Entry *_rslt, float *_key, int _k)
{
	for (int i = 0; i < _k; i ++)
	{
		if (_dist < _key[i])
		{
			for (int j = _k - 1; j > i; j --)
			{
				_rslt[j] = _rslt[j - 1];
				_key[j] = _key[j - 1];
			}
			_rslt[i] = *_e;
			_key[i] = _dist;
			return i;
		}
	}
	error("Error in update_rslt\n", true);
	return -1;
}
//------------------------------------------------------------
void RTree::write_header(char *buffer)
{
    int i;

    memcpy(buffer, &dimension, sizeof(dimension));
    i = sizeof(dimension);

    memcpy(&buffer[i], &num_of_data, sizeof(num_of_data));
    i += sizeof(num_of_data);

    memcpy(&buffer[i], &num_of_dnodes, sizeof(num_of_dnodes));
    i += sizeof(num_of_dnodes);

    memcpy(&buffer[i], &num_of_inodes, sizeof(num_of_inodes));
    i += sizeof(num_of_inodes);

    memcpy(&buffer[i], &root_is_data, sizeof(root_is_data));
    i += sizeof(root_is_data);

    memcpy(&buffer[i], &root, sizeof(root));
    i += sizeof(root);
}
//------------------------------------------------------------
//---added for valdity region queries-------------------------
// perform a window query mbr to get the query result into in_objs, put the outer objects retrieved into out_objs_so_far
void RTree::rect_win_query(float *mbr, LinList *in_objs, LinList *out_objs_so_far)
{
    load_root();

    root_ptr->rect_win_query(mbr, in_objs, out_objs_so_far);

	delete root_ptr;
	root_ptr = NULL;
}

// perform a window query mbr (excluding the window excr) to get the query result into c_inf_objs
void RTree::rect_win_query(float *mbr, float *exclmbr, LinList *c_inf_objs)
{
    load_root();

    root_ptr->rect_win_query(mbr, exclmbr, c_inf_objs);

	delete root_ptr;
	root_ptr = NULL;
}

void RTree::BFNN(float *_qpt, int _k, Entry *_rslt)
{
	//first init an array for storing the keys of retrieve objects
	float *key = new float [_k];
	for (int i = 0; i < _k; i ++)
		key[i] = (float) MAXREAL; //initially they are infinity
	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------

	int son = root; //this entry is to be visited next
	while (son != -1)
	{
		RTNode *rtn = new RTNode(this, son);
		for (int i = 0; i < rtn -> num_entries; i ++)
		{
			float edist = MINDIST(_qpt, rtn->entries[i].bounces, dimension);
			if (rtn->level == 0)
			{
				if (edist < key[_k - 1])
					update_rslt(&(rtn->entries[i]), edist, _rslt, key, _k);
			}
			else
			{
				if (edist<key[_k - 1])
					//meaning that edist is valid and we insert it to heap
				{
					HeapEntry *he = new HeapEntry();
					he -> key = edist;
					he -> level = rtn -> level;
					he -> son1 = rtn->entries[i].son;
					heap -> insert(he);
					delete he;
				}
			}
		}

		delete rtn;

		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{
				if (he->key>key[_k - 1]) //the algorithm terminates
					son = -1;
				else if (he->level == 0)
						//protection. if you see this message, debug
						error("testing... leaf entries found in heap\n", true);
					 else
						son=he->son1;
			}
		}
		delete he;
		//--------------------------------------------------------
	}

	delete [] key;
	delete heap;
}



/*
Nasim: BFNN_with_time_uncertainty:similar to BFNN_with_time except that, after obtaining
k neighbours, a probabilistic score is attached with the results
*/
void RTree::BFNN_with_keywords_time_uncertainty(float *_qpt, RTIME *_qtime,int *keywords, char *BFNN, int _k, int threshold, Entry *_rslt, float *probability, int *total_result_found)
{
	//first init an array for storing the keys of retrieve objects
	float *key = new float [_k];
	for (int i = 0; i < _k; i ++)
		key[i] = (float) MAXREAL; //initially they are infinity
	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	
	int heap_entry_count=0;
	int heap_remove_count=0;
	int update_rslt_call_count=0;
	
	//use threshold value given as hour, to adjust upper and lower time
	//in paper, threshold is denoted as alpha
	_qtime[0]-=(threshold*3600);
	_qtime[1]+=(threshold*3600);

	
	//FILE *fp = fopen(BFNN,"a");
	//fprintf(fp,"				query_no=%d\n",query_count++);
	//fprintf(fp,"----------------------------------------------------------------\n");
	//fprintf(fp,"keywords[0]=%d keywords[1]=%d keywords[2]=%d keywords[3]=%d\n",keywords[0],keywords[1],keywords[2],keywords[3]);
	io_access = 0;
	int son = root; //this entry is to be visited next
	while (son != -1)
	{
		//fprintf(fp,"----------------------------------------------------------------\n");
		RTNode *rtn = new RTNode(this, son);
		io_access++;
		//fprintf(fp,"son=%d io_access=%d rtn->num_entries=%d\n\n",son,io_access,rtn->num_entries);		
		for (int i = 0; i < rtn -> num_entries; i ++)
		{
			//fprintf(fp,"\t\t\tson=%d's i=%d'th entry :\n",son,i);
			//fprintf(fp,"==============================\n");
			float edist = MINDIST(_qpt, rtn->entries[i].bounces, dimension);
			if (rtn->level == 0)
			{
				//fprintf(fp,"leaf_node:\n");
				if (edist < key[_k - 1])
				{
					//fprintf(fp,"(edist<key[_k-1])==true\n");
					// Sarowar :	
					// time overlap is checked here : only difference from normal BFNN
					if( rtn->entries[i].section_time( _qtime ) != S_NONE   && (rtn->entries[i].section_keywords(keywords) != S_NONE) )
					{
						//fprintf(fp,"\n(time && keywords)==true\nupdate_rslt()::called.\n");
						update_rslt(&(rtn->entries[i]), edist, _rslt, key, _k);
						update_rslt_call_count++;
					}
				}
			}
			else
			{
				if (edist<key[_k - 1])
				{	
					//fprintf(fp,"internal_node:\n\n");
					//fprintf(fp,"(edist<key[_k-1])==true\n");			
					// Sarowar : 
					// time overlap is checked here : only difference from normal BFNN
					if( rtn->entries[i].section_time( _qtime ) != S_NONE && (rtn->entries[i].section_keywords(keywords) != S_NONE) )
						//meaning that edist is valid and we insert it to heap
					{
						//fprintf(fp,"\n(time && keywords)==true\n\n");
						HeapEntry *he = new HeapEntry();
						he -> key = edist;
						he -> level = rtn -> level;
						he -> son1 = rtn->entries[i].son;
						heap -> insert(he);
						heap_entry_count++;
						//fprintf(fp,"Entered in heap:\nhe->key=%f he->level=%d heap->son=%d\n",he->key,he->level,he->son1);
						//fprintf(fp,"(//--------------------------------------------------------------------//\n\n");
						delete he;
					}
				}
			}
		}

		delete rtn;

		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;
			if (!heap->remove(he))  //heap is empty
			{
				//fprintf(fp,"***************************************************\n");
				//fprintf(fp,"\t\t\t HEAP IS EMPTY.\n");
				//fprintf(fp,"***************************************************\n");
				son = -1;
			}
			else
			{
				heap_remove_count++;
				if (he->key>key[_k - 1]) //the algorithm terminates
				{
					//fprintf(fp,"***************************************************\n");
					//fprintf(fp,"\t\t\t HEAP->KEY IS LARGER THAN RESULT.\n");
					//fprintf(fp,"***************************************************\n");
					son = -1;
				}
				else if (he->level == 0)
				{
						//protection. if you see this message, debug
						error("testing... leaf entries found in heap\n", true);
				}
				else
				{
					son=he->son1;
					//fprintf(fp,"***************************************************\n");
					//fprintf(fp,"\t\tHEAP->REMOVE()=%d.\n\t\tHEAP VALUE IS USED FOR NEXT SON.\n",son);
					//fprintf(fp,"***************************************************\n");
				}
			}
			
		}
		delete he;
		//--------------------------------------------------------
	}
	//fprintf(fp,"############################################################\n");
	//fprintf(fp,"heap_entry_count=%d heap_remove_count=%d update_rslt_call_count=%d\n",heap_entry_count,heap_remove_count,update_rslt_call_count);
	//fclose(fp);

	
	
	//k-nearest neighbour has been found, now we calculate probability 
	//score for each of them
	
	//nasim:debug----------------------------------------------------------------------
	//FILE *nasim = fopen("nasim.txt","w");
	//fprintf(nasim,"q0=%lld,q1=%lld\n",_qtime[0],_qtime[1]);

	*total_result_found=update_rslt_call_count;
	for(int i=0;i<update_rslt_call_count;i++)
	{
		RTIME q;
		float total_prob=0.0,ind_prob=0.0;
		//fprintf(nasim,"inside probability calculation for object %d->\n",i);
		//fprintf(nasim,"%s",asctime(localtime(&_rslt[i].times[0])));
		//fprintf(nasim,"%s\n",asctime(localtime(&_rslt[i].times[1])));
		for(q=_qtime[0];q<=_qtime[1];q+=3600)
		{
			int d_to_init=difftime(q,_rslt[i].times[0])/3600;
			int d_to_end=difftime(q,_rslt[i].times[1])/3600;

            int max_radius= max(d_to_end,d_to_init);
            int min_radius=min(d_to_end,d_to_init);

            if(min_radius>0)min_radius++;
            if(max_radius<0)max_radius--;
			for(int j=min_radius;j<=max_radius;j++)
            {
                if(j==0)continue;
                int curRadius = abs(j);

                ind_prob=1;
				for(int k=0;k<update_rslt_call_count;k++)
                {
                    if(i==k)continue;
                    int dist_to_init_time=difftime(q,_rslt[k].times[0])/3600;
                    int dist_to_end_time=difftime(q,_rslt[k].times[1])/3600;

					if((dist_to_init_time<0)==(dist_to_end_time<0))
                    {
                        double mx=max(abs(dist_to_end_time),abs(dist_to_init_time));
                        double mn=min(abs(dist_to_end_time),abs(dist_to_init_time));

                        if(mx<=curRadius)
                        {
                            ind_prob=0;
                            break;
                        }
                        else if(mn>=curRadius)
                        {

                        }
                        else
                            ind_prob*=((mx-curRadius)/((_rslt[k].times[1]-_rslt[k].times[0])/3600));
                    }
                    else
                    {
                        float m1=abs(dist_to_end_time);
                        float m2=abs(dist_to_init_time);
                        ind_prob*=(((m1>curRadius?(m1-curRadius):0)+(m2>curRadius?(m2-curRadius):0))/((_rslt[k].times[1]-_rslt[k].times[0])/3600));
                    }

                }
                total_prob+=ind_prob;
				
            }
		}

		total_prob/=(((_rslt[i].times[1]-_rslt[i].times[0])/3600)*((_qtime[1]-_qtime[0])/3600+1));
		//we will assign probability score here
		probability[i]=total_prob;
	}

	delete [] key;
	delete heap;
}

//Nasim: BFNN_with_time_uncertainty: ends here -----------------------------------

// Sarowar : BFNN_with_time(). here node's time needs to be overlapped with the _qtime.
// another approach can be time proximity rather than equality. 
// another problem is time is non periodic. i.e. rather than sat or sunday, it is y-m-d specific.
void RTree::BFNN_with_time_keywords(float *_qpt, RTIME *_qtime,int *keywords, char *BFNN, int _k, Entry *_rslt)
{
	//first init an array for storing the keys of retrieve objects
	float *key = new float [_k];
	for (int i = 0; i < _k; i ++)
		key[i] = (float) MAXREAL; //initially they are infinity
	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	//char *BFNN = "C:/CSGLBSQ/Datasets/experiment/BFNN_50.txt";
	//char *BFNN = "C:/CSGLBSQ/Datasets/experiment/BFNN_20000_k_20.txt";
	int heap_entry_count=0;
	int heap_remove_count=0;
	int update_rslt_call_count=0;
	//------------------------------------------------------------
	FILE *fp = fopen(BFNN,"a");
	fprintf(fp,"				query_no=%d\n",query_count++);
	fprintf(fp,"----------------------------------------------------------------\n");
	fprintf(fp,"keywords[0]=%d keywords[1]=%d keywords[2]=%d keywords[3]=%d\n",keywords[0],keywords[1],keywords[2],keywords[3]);
	io_access = 0;
	int son = root; //this entry is to be visited next
	while (son != -1)
	{
		fprintf(fp,"----------------------------------------------------------------\n");
		RTNode *rtn = new RTNode(this, son);
		io_access++;
		fprintf(fp,"son=%d io_access=%d rtn->num_entries=%d\n\n",son,io_access,rtn->num_entries);		
		for (int i = 0; i < rtn -> num_entries; i ++)
		{
			fprintf(fp,"\t\t\tson=%d's i=%d'th entry :\n",son,i);
			fprintf(fp,"==============================\n");
			float edist = MINDIST(_qpt, rtn->entries[i].bounces, dimension);
			if (rtn->level == 0)
			{
				fprintf(fp,"leaf_node:\n");
				if (edist < key[_k - 1])
				{
					fprintf(fp,"(edist<key[_k-1])==true\n");
					// Sarowar :	
					// time overlap is checked here : only difference from normal BFNN
					if( rtn->entries[i].section_time( _qtime ) != S_NONE   && (rtn->entries[i].section_keywords(keywords) != S_NONE) )
					{
						fprintf(fp,"\n(time && keywords)==true\nupdate_rslt()::called.\n");
						update_rslt(&(rtn->entries[i]), edist, _rslt, key, _k);
						update_rslt_call_count++;
					}
				}
			}
			else
			{
				if (edist<key[_k - 1])
				{	
					fprintf(fp,"internal_node:\n\n");
					fprintf(fp,"(edist<key[_k-1])==true\n");			
					// Sarowar : 
					// time overlap is checked here : only difference from normal BFNN
					if( rtn->entries[i].section_time( _qtime ) != S_NONE && (rtn->entries[i].section_keywords(keywords) != S_NONE) )
						//meaning that edist is valid and we insert it to heap
					{
						fprintf(fp,"\n(time && keywords)==true\n\n");
						HeapEntry *he = new HeapEntry();
						he -> key = edist;
						he -> level = rtn -> level;
						he -> son1 = rtn->entries[i].son;
						heap -> insert(he);
						heap_entry_count++;
						fprintf(fp,"Entered in heap:\nhe->key=%f he->level=%d heap->son=%d\n",he->key,he->level,he->son1);
						fprintf(fp,"(//--------------------------------------------------------------------//\n\n");
						delete he;
					}
				}
			}
		}

		delete rtn;

		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;
			if (!heap->remove(he))  //heap is empty
			{
				fprintf(fp,"***************************************************\n");
				fprintf(fp,"\t\t\t HEAP IS EMPTY.\n");
				fprintf(fp,"***************************************************\n");
				son = -1;
			}
			else
			{
				heap_remove_count++;
				if (he->key>key[_k - 1]) //the algorithm terminates
				{
					fprintf(fp,"***************************************************\n");
					fprintf(fp,"\t\t\t HEAP->KEY IS LARGER THAN RESULT.\n");
					fprintf(fp,"***************************************************\n");
					son = -1;
				}
				else if (he->level == 0)
				{
						//protection. if you see this message, debug
						error("testing... leaf entries found in heap\n", true);
				}
				else
				{
					son=he->son1;
					fprintf(fp,"***************************************************\n");
					fprintf(fp,"\t\tHEAP->REMOVE()=%d.\n\t\tHEAP VALUE IS USED FOR NEXT SON.\n",son);
					fprintf(fp,"***************************************************\n");
				}
			}
			
		}
		delete he;
		//--------------------------------------------------------
	}
	fprintf(fp,"############################################################\n");
	fprintf(fp,"heap_entry_count=%d heap_remove_count=%d update_rslt_call_count=%d\n",heap_entry_count,heap_remove_count,update_rslt_call_count);
	fclose(fp);
	delete [] key;
	delete heap;
}


void RTree::print_rtree_nodes(char *file)
{
	//char *file="C:/CSGLBSQ/Datasets/experiment/printed_in_dfs_preorder_d_50.txt";
	FILE *fp;
	if((fp = fopen(file,"w")) == NULL) 
		error("Cannot open R-Tree text file", TRUE);
    else
    {
		visit_order=0;
		leaf_count=0;
		negative_bitmask_count=0;
		fprintf(fp,"			ROOT");
		dfs_visit_print(fp,root);
	}
	fprintf(fp,"\n\n------------------------------------------------\n--------------------------------------------\n");
	fprintf(fp,"leaf_count=%d\nnegative_bitmask_count=%d\n",leaf_count,negative_bitmask_count);
	fprintf(fp,"internal_node=%d\n",internal_node);
	internal_node=0;
	visit_order=0;
	leaf_count = 0;
	negative_bitmask_count=0;
	fclose(fp);
}

void RTree::dfs_visit_print(FILE *fp, int son)
{
	RTNode *rtn = new RTNode(this, son);
	
	fprintf(fp,"\n---------------------------------------------------------------------\n");
	fprintf(fp,"\n---------------------------------------------------------------------\n");
	fprintf(fp,"visit_order=%d son=%d rtn -> num_entries=%d level=%d\n\n",visit_order++,son,rtn -> num_entries,rtn->level);
	for (int i = 0; i < rtn -> num_entries; i ++)
	{
		fprintf(fp,"i=%d\nx1=%f x2=%f y1=%f y2=%f\n",i,rtn->entries[i].bounces[0],rtn->entries[i].bounces[1],rtn->entries[i].bounces[2],rtn->entries[i].bounces[3]);
		//fprintf(fp,"\nentries[ %d ].son=%d\n--------------\n",i,rtn->entries[i].son);
		fprintf(fp,"\nentries[ %d ].son=%d\n\n",i,rtn->entries[i].son);
		fprintf(fp,"TIME:   t1=%lld t2=%lld\n\n",rtn->entries[i].times[0],rtn->entries[i].times[1]);
		fprintf(fp,"FRUITS: b0=%d   b1=%d   b2=%d   b3=%d\n",rtn->entries[i].fruit_bitmask[0],rtn->entries[i].fruit_bitmask[1],rtn->entries[i].fruit_bitmask[2],rtn->entries[i].fruit_bitmask[3]);
		fprintf(fp,"//------------------------------------------------------//\n\n");
	}
	if(rtn->level == 0 )
	{
		leaf_count++;// += rtn->num_entries;
		for(int i=0; i<rtn->num_entries ; i++)
		{
			for(int j=0; j<4; j++)
			{
				if( rtn->entries[i].fruit_bitmask[j] < 0 )
				{
					negative_bitmask_count++;
				}
			}
		}
		delete rtn;
		return;
	}
	internal_node++;
	for (int i = 0; i < rtn -> num_entries; i ++)
	{
		//visit
		dfs_visit_print(fp,rtn->entries[i].son);
	}	
	delete rtn;		
}

void RTree::BFNN_with_time(float *_qpt, RTIME *_qtime, int _k, Entry *_rslt)
{
	//first init an array for storing the keys of retrieve objects
	float *key = new float [_k];
	for (int i = 0; i < _k; i ++)
		key[i] = (float) MAXREAL; //initially they are infinity
	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------

	int son = root; //this entry is to be visited next
	while (son != -1)
	{
		RTNode *rtn = new RTNode(this, son);
		for (int i = 0; i < rtn -> num_entries; i ++)
		{
			float edist = MINDIST(_qpt, rtn->entries[i].bounces, dimension);
			if (rtn->level == 0)
			{
				if (edist < key[_k - 1])
					// time overlap is checked here : only difference from normal BFNN
					if( rtn->entries[i].section_time( _qtime ) != S_NONE )
					{
						update_rslt(&(rtn->entries[i]), edist, _rslt, key, _k);
					}
			}
			else
			{
				if (edist<key[_k - 1])
				{	
					// time overlap is checked here : only difference from normal BFNN
					if( rtn->entries[i].section_time( _qtime ) != S_NONE )
						//meaning that edist is valid and we insert it to heap
					{
						HeapEntry *he = new HeapEntry();
						he -> key = edist;
						he -> level = rtn -> level;
						he -> son1 = rtn->entries[i].son;
						heap -> insert(he);
						delete he;
					}
				}
			}
		}

		delete rtn;

		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{
				if (he->key>key[_k - 1]) //the algorithm terminates
					son = -1;
				else if (he->level == 0)
						//protection. if you see this message, debug
						error("testing... leaf entries found in heap\n", true);
					 else
						son=he->son1;
			}
		}
		delete he;
		//--------------------------------------------------------
	}

	delete [] key;
	delete heap;
}

//Sarowar : BFNN_with_time() -->  End ---------------------------------------------------

//RTree::BFNN --- Modified by Tanzima for Rectangle

float RTree::KMin(Point2D m, int c, float cl, int k, Pointlocation _rslt[], int *num_of_data, float d_safe_corner,DistfromPOI _cornertoPOI[])
{
	float *kmindist = new float[k];
	float tempdist;
	for(int i=0; i<k; i++)
	{
			kmindist[i]=(float) MAXREAL;
	}

	for(int i=0; i<*num_of_data;i++)
	{
		if(cl* _cornertoPOI[i].d[c] <= d_safe_corner)
		{
			tempdist=Dist(_rslt[i],m);

			for(int j=0; j<k; j++)
			{

				if(kmindist[j]>tempdist)
				{
					for(int l=k-1;l>j;l--)
					{
						kmindist[l]= kmindist[l-1];
					}
					kmindist[j]= tempdist;
					break;
				}
			}
		}

	}
	return kmindist[k-1];
}


void RTree::UpdateCount(Rectangle1 R, float cl, int k, Pointlocation _rslt[], int *num_of_data, float d_safe_corner, int *count, DistfromPOI _cornertoPOI[])
{
	Point2D c[4];
	c[0][0] = R.x1; c[0][1]=R.y1;
	c[1][0] = R.x1; c[1][1]=R.y2;
	c[2][0] = R.x2; c[2][1]=R.y2;
	c[3][0] = R.x2; c[3][1]=R.y1;


	_cornertoPOI[*num_of_data-1].d[0]= Dist(_rslt[*num_of_data-1],c[0]);
	_cornertoPOI[*num_of_data-1].d[1]= Dist(_rslt[*num_of_data-1],c[1]);
	_cornertoPOI[*num_of_data-1].d[2]= Dist(_rslt[*num_of_data-1],c[2]);
	_cornertoPOI[*num_of_data-1].d[3]= Dist(_rslt[*num_of_data-1],c[3]);

	for(int j=0; j<4; j++)
	{
		if(count[j]<k)
		{
			count[j]=0;
			for(int i=0; i<*num_of_data;i++)
			{
				if(cl*_cornertoPOI[i].d[j] <= d_safe_corner)	count[j]++;

			}
		}
	}
	/*
	for(int i=0; i<_rslt.size();i++)
	{
		if(i<_cornertoPOI.size())
		{
			if(cl*_cornertoPOI[i].d[0] <= d_safe_corner)	count[0]++;
			if(cl*_cornertoPOI[i].d[1] <= d_safe_corner)	count[1]++;
			if(cl*_cornertoPOI[i].d[2] <= d_safe_corner)	count[2]++;
			if(cl*_cornertoPOI[i].d[3] <= d_safe_corner)	count[3]++;

		}
		else
		{
			Pointlocation p = _rslt[i];
			dist_c_P.d[0]= Dist(p,c[0]);
			if(cl*dist_c_P.d[0] <= d_safe_corner)	count[0]++;
			dist_c_P.d[1]= Dist(p,c[1]);
			if(cl*dist_c_P.d[1] <= d_safe_corner)	count[1]++;
			dist_c_P.d[2]= Dist(p,c[2]);
			if(cl*dist_c_P.d[2] <= d_safe_corner)	count[2]++;
			dist_c_P.d[3]= Dist(p,c[3]);
			if(cl*dist_c_P.d[3] <= d_safe_corner)	count[3]++;
			_cornertoPOI.push_back(dist_c_P);

		}
	}
	*/
}

int RTree::UpdateStatus(Rectangle1 R, float cl, int k, Pointlocation _rslt[], int *num_of_data, int *count, DistfromPOI _cornertoPOI[])
{
	Point2D o;
	o[0] = (R.x1+R.x2)/2;
	o[1] = (R.y1+R.y2)/2;
	float r = Dist(_rslt[*num_of_data-1],o);

	Point2D c[4];
	c[0][0] = R.x1; c[0][1]=R.y1;
	c[1][0] = R.x1; c[1][1]=R.y2;
	c[2][0] = R.x2; c[2][1]=R.y2;
	c[3][0] = R.x2; c[3][1]=R.y1;

	float d1=Dist(c[0],c[1]);
	float d2=Dist(c[1],c[2]);
	float d_safe_corner = r-0.5*sqrt(d1*d1+d2*d2);




	UpdateCount(R,cl,k,_rslt,num_of_data,d_safe_corner,count,_cornertoPOI);


	for(int i=0; i<4; i++)
	{
		if (count[i]<k)	return 0;
	}



	float d_i, d_j,d_max_m,d_max = 0;
	int i,j;
	Point2D m;
	for (i=0; i<4; i++)
	{
		j=(i+1)%4;

		m[0]= (c[i][0]+c[j][0])/2;
		m[1]= (c[i][1]+c[j][1])/2;
		d_i= KMin(m,i,cl,k,_rslt,num_of_data,d_safe_corner,_cornertoPOI);
		d_j= KMin(m,j,cl,k,_rslt,num_of_data,d_safe_corner,_cornertoPOI);
		if(d_i>d_j) d_max_m = d_i;
		else d_max_m = d_j;
		if(d_max_m > d_max)
			d_max = d_max_m;

	}
	float c_dmax=d1;
	if(d2>d1)	c_dmax=d2;
	float d_safe = r- 0.5 * c_dmax;

	if(cl*d_max>d_safe)
		return ceil(r+cl*(d_max-d_safe));
	else
		return -1;
}
void RTree::Rect_kNNQ(Rectangle1 R, float cl, int k, Pointlocation _rslt[], int *num_of_data)
{
	int io=0;
	//init status
	int count[4];
	for(int i=0; i<4; i++)	count[i]=0;
	int status=0;
	Point2D o;
	o[0] = (R.x1+R.x2)/2;
	o[1] = (R.y1+R.y2)/2;
	DistfromPOI cornertoPOI[10000];

	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------

	int son = root; //this entry is to be visited next
	while (son != -1 && status != -1)
	{


		RTNode *rtn = new RTNode(this, son);
		//Experiment
		io_access++;
		for (int i = 0; i < rtn -> num_entries; i++)
		{
			float o1[2];
			o1[0]=(float)o[0];
			o1[1]=(float)o[1];
			float edist = MINDIST(o1, rtn->entries[i].bounces, dimension);

			HeapEntry *he = new HeapEntry();
			he -> key = edist;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;

		}

		delete rtn;


		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{
				if (he->level == 0) //p is an object
				{
					if(*num_of_data==10000)
						//printf("\nGreater than 10000\n");
						error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);

					_rslt[*num_of_data].x=he->x1;
					_rslt[*num_of_data].y=he->y1;
					*num_of_data=*num_of_data+1;
					//if(*num_of_data%1000==0) printf("\n%d",*num_of_data);



					if(status==0)
					{
						status = UpdateStatus(R,cl,k,_rslt,num_of_data,count,cornertoPOI);

						if(status != -1)	again=true;
					}
					else if (status < Dist(_rslt[*num_of_data-1],o))
					{
						status = -1;
					}
					else
					{
						again = true;
					}
				}
				else
				{
					if(status>0)
					{
						if(status<he->key)
							status = -1;
						else
							son=he->son1;
					}
					else
					{
						son=he->son1;
					}
				}
			}
		}

		delete he;
	}
	delete heap;
}

void RTree::Point_BFN_NNQ(Point2D o, double *_rslt)
{

	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------

	int son = root; //this entry is to be visited next
	while (son != -1)
	{
		RTNode *rtn = new RTNode(this, son);
		io_access++;
		for (int i = 0; i < rtn -> num_entries; i ++)
		{
			float o1[2];
			o1[0]=(float)o[0];
			o1[1]=(float)o[1];
			float edist = MINDIST(o1, rtn->entries[i].bounces, dimension);

			HeapEntry *he = new HeapEntry();
			he -> key = edist;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;

		}

		delete rtn;


		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{
				if (he->level == 0) //p is an object
				{
					_rslt[0] = he->x1;
					_rslt[1] = he->y1;
					son=-1;

				}
				else
				{
					son=he->son1;

				}
			}
		}

		delete he;
	}
	delete heap;
}
//END




// The following code was copied from the implementation of TP-kNN queries from Tony
void RTree::TPNN_TP(float *_qline, int _k, Entry *_nn, Entry *_rslt, float _max_trvl)
{
	float key = _max_trvl;
	  // the minimum distance that the query point must travel

//we comment these lines to avoid initing the heap everytime--
//this function is called
//Heap *heap = new Heap();
//heap->init(dimension);
//------------------------------------------------------------
	if (tpheap==NULL)
		error("tpheap is not initialized\n", true);
	tpheap->used=0;

	int son = root;
	while (son != -1)
	{
		RTNode *rtn = new RTNode(this, son);
		for (int i = 0; i < rtn -> num_entries; i ++)
		{
			//first create an array m2 for e to cal dist----------
			float *m2 = new float [4 * dimension];
			memset(m2, 0, 4 * dimension * sizeof(float));
			memcpy(m2, rtn -> entries[i].bounces, 2 * dimension * sizeof(float));
			//----------------------------------------------------
//testing--------------------------
//if (rtn->entries[i].son==573673 && cnt==84-1)
//	printf("testing...\n");
//---------------------------------

			float edist = (float) MAXREAL;
			if (rtn -> level == 0)
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^leaf node case^^^^^^^^^^^^^^^^^^^^^
			{
				int eNNsub=-1;
				//if e (i.e., m2) is one of the current NN points-
				//eNNsub stores its subsript in _nn; otherwise,
				//eNNsub=-1
				for (int j=0; j<_k; j++)
					if (rtn->entries[i].son==_nn[j].son)
					{ eNNsub=j; j=_k;}
				//------------------------------------------------

				if (eNNsub==-1 || SEQUENCE_SENSITIVE)
					//namely, if sequence insensitive and e is indeed a NN
					//then we do not handle this entry
				{

					float *m1 = new float [4 * dimension];
					//find the NN that leads to the minimum-------
					//influence time
					int nn_num=-1; //the subsript of the NN to be found
					for (int j = 0; j < _k; j ++)
						//for each of the NN found in the 1st step
					{
						bool yesdo=true; //whether to compute
						  //the inflence time of nn[j] and e
						if (j==eNNsub) yesdo=false;

						//check if this pair has produced --------
						//influence time before
						for (int l=0; l<PAST_PAIR; l++)
							if (min(_nn[j].son, rtn->entries[i].son)==min(last_pair[l][0], last_pair[l][1]) &&
								max(_nn[j].son, rtn->entries[i].son)==max(last_pair[l][0], last_pair[l][1]))
							{	yesdo=false; l=PAST_PAIR; }
						//----------------------------------------

						if (yesdo)
						{
//these codes use NNinf===========================================
/*
							//first create an array m1 (for nn[j])
							//to compute dist
							memset(m1, 0, 4 * dimension * sizeof(float));
							memcpy(m1, _nn[j].bounces, 2 * dimension * sizeof(float));
							//get the influence time of m2--------
							//(for entry e) with respect to(nn[j])
							float this_inf = NNinf(m1, m2, _qline, dimension);
							//------------------------------------
*/
//================================================================

//these codes use NNinf2==========================================
							//first create an array m1 (for nn[j])
							//to compute dist
							memset(m1, 0, 4 * dimension * sizeof(float));
							memcpy(m1, _nn[j].bounces, 2 * dimension * sizeof(float));
							m1[1]=m1[2]; m2[1]=m2[2];
							//create an arry m3 for _qline----------------
							float *m3=new float[2*dimension];
							m3[0]=_qline[0]; m3[1]=_qline[2];
							m3[2]=_qline[4]; m3[3]=_qline[6];
							//--------------------------------------------
							//get the influence time of m2--------
							//(for entry e) with respect to(nn[j])
							float this_inf = NNinf2(m1, m2, m3);
							//------------------------------------
							delete []m3;
//================================================================

							if (this_inf>0 && this_inf<edist)
								//this_inf=0 means that there is another point that has the same distance
								//to the current query position as the current NN. In this implementation,
								//we choose to ignore handling such special cases, which, however, may cause
								//problems for datasets with big cardinality
//							if (this_inf>=0 && this_inf<edist)
							{
								edist=this_inf; nn_num=j;
							}
						}  //END if (yesdo)
					}//END checking all neighbors
					//-------------------------------------------------
					//if (edist<key && edist!=0)
					if (edist<key)
					{
						update_rslt(&(rtn->entries[i]), edist, _rslt, &key, 1);
						_rslt->nn_num=nn_num;
					}
					delete []m1;
				}
			}
//^^^^^^^^^^^^^^^^^^^^^^^^^non-leaf node case^^^^^^^^^^^^^^^^^^^^^
			else
				//Next handle non-leaf node case
			{
				float *m1 = new float [4 * dimension];
				for (int j = 0; j < _k; j ++)
				{
					//first create an array m1 to cal dist--------
					memset(m1, 0, 4 * dimension * sizeof(float));
					memcpy(m1, _nn[j].bounces, 2 * dimension * sizeof(float));
					//--------------------------------------------
					float this_mininf = NNmininf(m1, m2, _qline, dimension);
					if (this_mininf < edist)
						edist = this_mininf;
				}
				delete [] m1;

				if (edist < key)
				{
					HeapEntry *he = new HeapEntry();
					he -> key = edist;
					he -> level = rtn -> level;
					he -> son1 = rtn->entries[i].son;
					tpheap -> insert(he);
					delete he;
				}
			}
			delete [] m2;

		}
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		delete rtn;

		//get next entry from the heap
		bool again = true;
		while (again)
		{
			again = false;
			HeapEntry *he = new HeapEntry();
			if (!tpheap -> remove(he))  //heap is empty
				son = -1;
			else
				if (he -> key > key)
					//the algorithm can terminate
					son = -1;
				else
					son = he -> son1;
			delete he;
		}
	}
//delete heap;
}



//RTree::private_kGNN --- Created by Tanzima for kGNN



void RTree::private_kGNN_sum(Rectangle1 R[], int g_size, int k, Pointlocation rslt[], int *num_of_data)
{
	//variable initialization
	int end=0;

	//Variables
	int i, j;

	double maxdist[kMAX];
	for(i=0; i<k; i++)
		maxdist[i] = MAXDOUBLE;

	//Find the MBR from provided rectangles
	Rectangle1 mbr;
	mbr.x1=R[0].x1;
	mbr.x2=R[0].x2;
	mbr.y1=R[0].y1;
	mbr.y2=R[0].y2;

	for(i=1; i<g_size; i++)
	{
		if(R[i].x1 < mbr.x1)	mbr.x1=R[i].x1;
		if(R[i].x2 > mbr.x2)	mbr.x2=R[i].x2;
		if(R[i].y1 < mbr.y1)	mbr.y1=R[i].y1;
		if(R[i].y2 > mbr.y2)	mbr.y2=R[i].y2;
	}
	//end

	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------

	int son = root; //this entry is to be visited next

	while (son != -1 && end==0)
	{
		RTNode *rtn = new RTNode(this, son);

		for (int i = 0; i < rtn -> num_entries; i++)
		{
			double edist1=0, edist2=0;
			edist1 = MINRECTDIST(mbr, rtn->entries[i].bounces);
			if((g_size*edist1)>maxdist[k-1])	continue;
			edist1 = MINRECTDIST(R[0], rtn->entries[i].bounces);
			for (j=1; j < g_size; j++)
			{
				edist1 += MINRECTDIST(R[j], rtn->entries[i].bounces);
			}

			if(edist1>maxdist[k-1])	continue;

			edist2 = MAXRECTDIST(R[0], rtn->entries[i].bounces);
			for (j=1; j < g_size; j++)
			{
				edist2 += MAXRECTDIST(R[j], rtn->entries[i].bounces);
			}

			//update maxdistk
			for(int l=0; l<k; l++)
			{
				if (edist2 < maxdist[l])
				{
					for(int j=k-1; j>l; j--)
					{
						maxdist[j]=maxdist[j-1];
					}
					maxdist[l]=edist2;
					break;
				}
			}

			HeapEntry *he = new HeapEntry();
			he -> key = edist1;
			he -> key1 = edist2;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;

		}

		delete rtn;

		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else if (he->key>maxdist[k-1])
			{
				end = 1;

			}
			else
			{
				if (he->level == 0) //p is an object
				{
						//enter into result set
						if(*num_of_data==MAXDATALIMIT)
							//printf("\nGreater than 10000\n");
							error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);

						rslt[*num_of_data].x=he->x1;
						rslt[*num_of_data].y=he->y1;
						rslt[*num_of_data].dmin=he->key;
						rslt[*num_of_data].dmax=he->key1;


						*num_of_data=*num_of_data+1;

						//get next data  from heap
						again =true;
				}
				else //not leaf node
				{
					son=he->son1;
				}
			}
		}

		delete he;
	}
	delete heap;
}


//END


void RTree::private_kGNN_max(Rectangle1 R[], int g_size, int k, Pointlocation rslt[], int *num_of_data)
{


	//variable initialization
	int end=0;

	//Variables
	int i, j;

	double maxdist[kMAX];
	for(i=0; i<k; i++)
		maxdist[i] = MAXDOUBLE;

	//Find the MBR from provided rectangles
	Rectangle1 mbr;
	mbr.x1=R[0].x1;
	mbr.x2=R[0].x2;
	mbr.y1=R[0].y1;
	mbr.y2=R[0].y2;

	for(i=1; i<g_size; i++)
	{
		if(R[i].x1 < mbr.x1)	mbr.x1=R[i].x1;
		if(R[i].x2 > mbr.x2)	mbr.x2=R[i].x2;
		if(R[i].y1 < mbr.y1)	mbr.y1=R[i].y1;
		if(R[i].y2 > mbr.y2)	mbr.y2=R[i].y2;
	}
	//end


	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------

	int son = root; //this entry is to be visited next


	while (son != -1 && end==0)
	{


		RTNode *rtn = new RTNode(this, son);

		for (int i = 0; i < rtn -> num_entries; i++)
		{
			double edist1=0, edist2=0, tdist=0;
			edist1 = MINRECTDIST(mbr, rtn->entries[i].bounces);
			if((edist1)>maxdist[k-1])	continue;


			edist1 = MINRECTDIST(R[0], rtn->entries[i].bounces);
			for (j=1; j < g_size; j++)
			{
				tdist = MINRECTDIST(R[j], rtn->entries[i].bounces);
				if(tdist>edist1)
					edist1=tdist;
			}

			if(edist1>maxdist[k-1])	continue;

			edist2 = MAXRECTDIST(R[0], rtn->entries[i].bounces);
			for (j=1; j < g_size; j++)
			{
				tdist= MAXRECTDIST(R[j], rtn->entries[i].bounces);
				if(tdist>edist2)
					edist2=tdist;
			}



			//update maxdistk
			for(int l=0; l<k; l++)
			{
				if (edist2 < maxdist[l])
				{
					for(int j=k-1; j>l; j--)
					{
						maxdist[j]=maxdist[j-1];
					}
					maxdist[l]=edist2;
					break;
				}
			}

			HeapEntry *he = new HeapEntry();
			he -> key = edist1;
			he -> key1 = edist2;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;

		}

		delete rtn;


		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else if (he->key>maxdist[k-1])
			{
				end = 1;

			}
			else
			{

				if (he->level == 0) //p is an object
				{
						//enter into result set
						if(*num_of_data==MAXDATALIMIT)
							//printf("\nGreater than 10000\n");
							error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);

						rslt[*num_of_data].x=he->x1;
						rslt[*num_of_data].y=he->y1;
						rslt[*num_of_data].dmin=he->key;
						rslt[*num_of_data].dmax=he->key1;

						*num_of_data=*num_of_data+1;
						//if(*num_of_data%1000==0) printf("\n%d",*num_of_data);

						//get next data  from heap
						again =true;
				}
				else //not leaf node
				{
					son=he->son1;
				}
			}
		}

		delete he;
	}
	delete heap;

}



//added by Eunus:


void RTree::kGNN(Rectangle1 R[], int g_size, int k, int f, Pointlocation rslt[], int *num_of_data)
{


	//variable initialization
	int end=0;

	//Variables
	int i, j;

	double maxdist[kMAX];
	for(i=0; i<k; i++)
		maxdist[i] = MAXDOUBLE;

	//Find the MBR from provided query points
	Rectangle1 mbr;
	mbr.x1=R[0].x1;
	//mbr.x2=R[0].x2;
	mbr.x2=R[0].x1;
	mbr.y1=R[0].y1;
	mbr.y2=R[0].y1;
	//mbr.y2=R[0].y2;

	for(i=1; i<g_size; i++)
	{
		if(R[i].x1 < mbr.x1)	mbr.x1=R[i].x1;
		//if(R[i].x2 > mbr.x2)	mbr.x2=R[i].x2;
		if(R[i].x1 > mbr.x2)	mbr.x2=R[i].x1;
		if(R[i].y1 < mbr.y1)	mbr.y1=R[i].y1;
		//if(R[i].y2 > mbr.y2)	mbr.y2=R[i].y2;
		if(R[i].y1 > mbr.y2)	mbr.y2=R[i].y1;
	}
	//end


	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------

	int son = root; //this entry is to be visited next


	while (son != -1 && end==0)
	{


		RTNode *rtn = new RTNode(this, son);

		for (int i = 0; i < rtn -> num_entries; i++)
		{
			double edist1=0, edist2=0, tdist=0;
			edist1 = MINRECTDIST(mbr, rtn->entries[i].bounces);
			if((edist1)>maxdist[k-1])	continue;

			Pointlocation pl;
			pl.x = R[0].x1;
			pl.y = R[0].y1;


			Rectangle1 nd;



			nd.x1 = rtn->entries[i].bounces[0];
			nd.x2 = rtn->entries[i].bounces[1];
			nd.y1 = rtn->entries[i].bounces[2];
			nd.y2 = rtn->entries[i].bounces[3];


			edist1 = MINRECTDIST1(pl, nd);
			for (j=1; j < g_size; j++)
			{
				pl.x = R[j].x1;
				pl.y = R[j].y1;
				tdist = MINRECTDIST1(pl, nd);
				//if(tdist>edist1)

				edist1=fDIST(edist1, tdist, f);
			}

			if(edist1>maxdist[k-1])	continue;



			HeapEntry *he = new HeapEntry();
			he -> key = edist1;
			//he -> key1 = edist2;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;

		}

		delete rtn;


		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else if (he->key>maxdist[k-1])
			{
				end = 1;

			}
			else
			{

				if (he->level == 0) //p is an object
				{
						//enter into result set
						if(*num_of_data==MAXDATALIMIT)
							//printf("\nGreater than 10000\n");
							error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);


					double aggdist = 0.0;


					Pointlocation xp;
					xp.x  = he->x1;
					xp.y = he->y1;
					aggdist = MINRECTDIST1(xp,mbr);

					if(aggdist>maxdist[k-1])	continue;

					Point2D dp;
					dp[0] = he->x1;
					dp[1] = he->y1;
					for(i=0; i<g_size;i++){

							//HeapEntry *hq = new HeapEntry();

							Point2D qp; //query point
							qp[0] = R[i].x1;
							qp[1] = R[i].y1;

							double cdist = Dist(dp,qp);//edist1; distance between q and p

							if(i==0){
								aggdist = cdist;
							}else{
								aggdist = fDIST(aggdist,cdist,f);
							}
					}

					if((aggdist)>maxdist[k-1]){
						end = 1;
						continue;
					}

				 //update maxdistk
					for(int l=0; l<k; l++)
					{
						if (aggdist < maxdist[l])
						{
							for(int j=k-1; j>l; j--)
							{
								maxdist[j]=maxdist[j-1];
							}
							maxdist[l]=aggdist;
							break;
						}
					}


						if(aggdist != he->key)
						{
							//printf("something is wrong!!");
						}

						rslt[*num_of_data].x=he->x1;
						rslt[*num_of_data].y=he->y1;
						rslt[*num_of_data].dmin=aggdist;//he->key; //aggregate distance
						rslt[*num_of_data].dmax=he->key1;

						rslt[*num_of_data].oid=he->son1; //oid

						*num_of_data=*num_of_data+1;
						//if(*num_of_data%1000==0) printf("\n%d",*num_of_data);

						//get next data  from heap
						again =true;
				}
				else //not leaf node
				{
					son=he->son1;
				}
			}
		}

		delete he;
	}
	delete heap;

}





void RTree::Naive_Consensus_kGNN(Rectangle1 R[], int g_size, int sg_size, int k, int f, SG L[][16], int *num_of_data)
{

	//Variables
	int i, j;

	for(i=0; i<g_size-sg_size+1; i++){
		//BestAggDist[i] = MAXDOUBLE;
		//status[i] = false; //temporary; true-> final
		for(j=0; j<k; j++){
			//need to initialize
			//SG sg;
			//L[i][j] = sg; //RESULT List

			L[i][j].sgDist = MAXDOUBLE;
			L[i][j].m = sg_size+i; //subgrop size
			L[i][j].o = -1;
			L[i][j].gDist = MAXDOUBLE;
			L[i][j].status = false;
		}

	}

	//group index
	int n[32];
	for(i=0;i<g_size;i++){
		n[i] = i;//+'0';
	}
	//n[i] = '\0';

	//for 8 : 5000
	//for 16: 20000
	// for 20: 170000
	// for 24: 2500000
	//for 32: set 130000000

	int **comb;

	comb = new int*[2500000];
	for(int j = 0; j < 2500000; j++)
		comb[j] = new int[32];


	//evaluate for each subgroup and populate L[i]
	for(i=sg_size; i <= g_size; i++){

		//int comb[2000][32];


		int r[32];


		//r[i+1] = '\0';

		ind_comb = 0;

		if(i<g_size){

			combination(n,0,r,i,0,g_size-i,comb);
		}
		else{
			//strcpy(comb[0],n);
			for(j=0;j<g_size;j++){
				comb[0][j]=n[j];//+'0';
			}
		}

		int cnr; //number of combination
		calcCNR( g_size, i, &cnr );

		Rectangle1 GR[32];

		printf("g_size: %d sg_size: %d cnr: %d \n",g_size, i, cnr);

		for(int c = 0; c<cnr; c++){
			//comb[c][i]='\0';
			//printf("%s\n",comb[c]);

			for(int l =0; l < i; l++){
				int ind = comb[c][l];//-'0';
				//printf("%d ",comb[c][l]);

				GR[l].x1 = R[ind].x1;
				GR[l].y1 = R[ind].y1;
				GR[l].x2 = R[ind].x2;
				GR[l].y2 = R[ind].y2;
			}

			//printf("\n");

			//call kGNN
			Pointlocation rslt[100];
			int no_of_data = 0;

			kGNN(GR, i, k, f, rslt, &no_of_data);

			*num_of_data=no_of_data+1;

			for(int no = 0; no<no_of_data; no++)
			{
				//if the min aggdist is greater than the bestaggdis
				if(rslt[no].dmin >  L[i-sg_size][k-1].sgDist)
					break;

				for(int l=0; l<k; l++)
				{
					if (rslt[no].dmin < L[i-sg_size][l].sgDist)
					{
						for(int j=k-1; j>l; j--)
						{
							//maxdist[j]=maxdist[j-1];
							L[i-sg_size][j].sgDist = L[i-sg_size][j-1].sgDist;
							L[i-sg_size][j].m = L[i-sg_size][j-1].m;
							L[i-sg_size][j].o = L[i-sg_size][j-1].o;
							//L[i-1][j].aggDist = L[i-1][j-1].aggDist;
							//L[i-sg_size][j].gDist = L[i-sg_size][j-1].gDist;
							//L[i-sg_size][j].status = L[i-sg_size][j-1].status;

							strcpy(L[i-sg_size][j].sgList,L[i-sg_size][j-1].sgList);

						}
						//maxdist[l]=edist2;
						L[i-sg_size][l].sgDist = rslt[no].dmin;
						L[i-sg_size][l].m = sg_size;
						L[i-sg_size][l].o = rslt[no].oid;



						//L[i-sg_size][l].gDist = sgi.gDist; //total aggregate distance

						//L[i-sg_size][l].status = sgi.status;
						//strcpy(L[i-sg_size][l].sgList,comb[c]);

						//calculate total aggregate distance
						double aggdist = 0.0;
						Point2D dp;
						dp[0] = rslt[no].x;
						dp[1] = rslt[no].y;
						for(int w =0; w < g_size; w++){
							Point2D qp; //query point
							qp[0] = R[w].x1;
							qp[1] = R[w].y1;

							double cdist = Dist(dp,qp);//edist1; distance between q and p

							if(w==0){
								aggdist = cdist;
							}else{
								aggdist = fDIST(aggdist,cdist,f);
							}
						}
						L[i-sg_size][l].gDist = aggdist; //total aggregate distance




						break;
					} //if
				} //for
			}

		}//for all combination

		printf("group size : %d finished",i);
	} // for all grp size


	for(j = 0; j < 2500000; ++j)
			delete [] comb[j];
		delete [] comb;

}




//I AM using Rectangle1 as query points: fann for

void RTree::FANN_consensus_kGNN(Rectangle1 R[], int g_size, int sg_size, int k, int f, SG L[][16], int lind, int *num_of_data, int tb_lb)
{


	//variable initialization
	int end=0;

	//Variables
	int i, j;

	double maxdist[kMAX];
	for(i=0; i<k; i++)
		maxdist[i] = MAXDOUBLE;



	double BestAggDist[128];
	bool status[128];

	//SG L[g_size-sg_size+1][kMAX]; //result list

	//for(i=0; i<g_size-sg_size+1; i++){
	i = lind;
	BestAggDist[i] = MAXDOUBLE;
		status[i] = false; //temporary; true-> final
		for(j=0; j<k; j++){
			//need to initialize
			//SG sg;
			//L[i][j] = sg; //RESULT List

			L[i][j].sgDist = MAXDOUBLE;
			L[i][j].m = sg_size; //+i; //subgrop size
			L[i][j].o = -1;
			L[i][j].gDist = MAXDOUBLE;
			L[i][j].status = false;


		}

	//}

	Point2D x;
	GEOCENTROID(R,g_size,x);
	//calculate maximum distance to query points
	float maxQDist = MAXGROUPDIST(R,g_size,x);


	//Find the MBR from provided rectangles
	Rectangle1 mbr;
	mbr.x1=R[0].x1;
	mbr.x2=R[0].x2;
	mbr.y1=R[0].y1;
	mbr.y2=R[0].y2;

	for(i=1; i<g_size; i++)
	{
		if(R[i].x1 < mbr.x1)	mbr.x1=R[i].x1;
		if(R[i].x2 > mbr.x2)	mbr.x2=R[i].x2;
		if(R[i].y1 < mbr.y1)	mbr.y1=R[i].y1;
		if(R[i].y2 > mbr.y2)	mbr.y2=R[i].y2;
	}
	//end


	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------

	int son = root; //this entry is to be visited next

	//son == -1 means empty && end = 0 means exit

	while (son != -1 && end==0)
	{


		RTNode *rtn = new RTNode(this, son);

		for (i = 0; i < rtn -> num_entries; i++)
		{
			double edist1=0;//edist2=0, tdist=0;
			//edist1 = MINRECTDIST(mbr, rtn->entries[i].bounces);

			Pointlocation xp;
			xp.x = x[0];
			xp.y = x[1];

			Rectangle1 nd;



			nd.x1 = rtn->entries[i].bounces[0];
			nd.x2 = rtn->entries[i].bounces[1];
			nd.y1 = rtn->entries[i].bounces[2];
			nd.y2 = rtn->entries[i].bounces[3];

			if(rtn->level == 0){ //child node is a object
				nd.x2 = nd.x1;
				nd.y2 = nd.y1;
			}

			edist1 = MINRECTDIST1(xp, nd); // mindist from x to node



			HeapEntry *h = new HeapEntry();
			h -> key = edist1;
			//he -> key1 = edist2;
			h -> level = rtn -> level;
			h -> son1 = rtn->entries[i].son;


			h-> x1 = nd.x1; //rtn->entries[i].bounces[0];
			h-> x2 = nd.x2; //rtn->entries[i].bounces[1]; //UNCOMMENTED BY EUNUS ; WAS commented IN sHETU
			h-> y1 = nd.y1; //rtn->entries[i].bounces[2];
			h-> y2 = nd.y2;//rtn->entries[i].bounces[3]; //UNCOMMENTED BY eUNUS
			heap -> insert(h);
			delete h;

		}

		delete rtn;


		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;

			bool isempty = !heap->remove(he);

			Rectangle1 p;
			p.x1 = he-> x1;
			p.y1 = he-> y1;
			p.x2 = he-> x2;
			p.y2 = he-> y2;




			bool isterminate = false;

			if(tb_lb ==1){ //Tight bound  method
				isterminate = !TERMINATE_TB(R, g_size,sg_size,x,p,maxQDist,f,BestAggDist,status);
			}else{ //lower bound method
				isterminate = !FANN_TERMINATE_LB(R, g_size,sg_size,x,p,maxQDist,f,BestAggDist,lind,status);
			}

			if (isempty)  //heap is empty
				son = -1;
			else if (isterminate)//(he->key>maxdist[k-1]) //check terminating condition
			{
				//bool TERMINATE_LB(Rectangle1 R[], int g_size, int sg_size, Point2D x, Rectangle1 p, float maxQDist, int f, double BestAggDist[], bool status[])

				end = 1;

			}
			else
			{

				if (he->level == 0) //p is an object ; code for data object
				{

					//enter into result set
						//if(*num_of_data==MAXDATALIMIT)
							//printf("\nGreater than 10000\n");
							//error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);
						/*
						rslt[*num_of_data].x=he->x1;
						rslt[*num_of_data].y=he->y1;
						rslt[*num_of_data].dmin=he->key;
						rslt[*num_of_data].dmax=he->key1;
						*/

						*num_of_data=*num_of_data+1;
						//if(*num_of_data%1000==0) printf("\n%d",*num_of_data);


						//Eunus Code

						Heap *heapQ = new Heap();
						heapQ->init(dimension);
						double totalDist = 0.0;
						//float p[2];
						//p[0] = he->x1;
						//p[1] = he->y1;

						Point2D dp; //data point
						dp[0] = he->x1;
						dp[1] = he->y1;

						for(i=0; i<g_size;i++){

							HeapEntry *hq = new HeapEntry();

							Point2D qp; //query point
							qp[0] = R[i].x1;
							qp[1] = R[i].y1;

							hq -> key = Dist(dp,qp);//edist1; distance between q and p
							//he -> key1 = edist2;
							//he -> level = rtn -> level;
							hq -> son1 = i;//rtn->entries[i].son;
							hq-> x1 = qp[0];//rtn->entries[i].bounces[0];
							//he-> x2 = rtn->entries[i].bounces[1];
							hq-> y1 = qp[1];//rtn->entries[i].bounces[2];
							//he-> y2 = rtn->entries[i].bounces[3];
							heapQ -> insert(hq);


							if(i==0){
								totalDist = hq -> key;
							}else{
								totalDist = fDIST(totalDist,hq -> key,f);
							}

							delete hq;
						} //end of creating heapQ for query based on the distnace from data point

						//popping elements from HeapQ

						double aggDist = 0.0;


						char sgList[129];

						for(i=0;i<129;i++)
							sgList[i] = '0';
						sgList[128] = '\0';

						i = 0;
						HeapEntry *hq = new HeapEntry();
						while (heapQ->remove(hq))
						{
							if(i==0){
								aggDist = hq -> key;
							}else{
								aggDist = fDIST(aggDist,hq -> key,f);
							}

							sgList[hq -> son1] = '1';

							i++;
							if(i>sg_size) break;
							//if(i >= sg_size && BestAggDist[i-sg_size] > aggDist){
							if(i == sg_size && BestAggDist[lind] > aggDist){

								SG sgi;
								sgi.m = i;
								sgi.o = he -> son1;
								sgi.sgDist = aggDist;
								sgi.gDist = totalDist;
								sgi.status = false;


								//sgi.sgList[64] = '\0';
								strcpy(sgi.sgList,sgList);



								//order the result list RL and update BestAggDist


								for(int l=0; l<k; l++)
								{


									if (sgi.sgDist < L[lind][l].sgDist)
									{
										for(int j=k-1; j>l; j--)
										{
											//maxdist[j]=maxdist[j-1];
											L[lind][j].sgDist = L[lind][j-1].sgDist;
											L[lind][j].m = L[lind][j-1].m;
											L[lind][j].o = L[lind][j-1].o;
											//L[i-1][j].aggDist = L[i-1][j-1].aggDist;
											L[lind][j].gDist = L[lind][j-1].gDist;
											L[lind][j].status = L[lind][j-1].status;

											strcpy(L[lind][j].sgList,L[lind][j-1].sgList);

										}
										//maxdist[l]=edist2;
										L[lind][l].sgDist = sgi.sgDist;
										L[lind][l].m = sgi.m;
										L[lind][l].o = sgi.o;
										L[lind][l].gDist = sgi.gDist;
										L[lind][l].status = sgi.status;
										strcpy(L[lind][l].sgList,sgi.sgList);

										//update BestAggDist
										//BestAggDist[i-sg_size] = sgi.sgDist;

										break;
									} //if
								} //for

								//update BestAggDist with the kth best
								BestAggDist[lind] = L[lind][k-1].sgDist;

							} //if



						} //while

						delete hq;
						//get next data  from heap
						again =true;
						delete heapQ;
				}
				else //not leaf node
				{
					son=he->son1;
				}
			}
		}

		delete he;
	}
	delete heap;

}




//I AM using Rectangle1 as query points

void RTree::consensus_kGNN(Rectangle1 R[], int g_size, int sg_size, int k, int f, SG L[][16], int *num_of_data, int tb_lb)
{


	//variable initialization
	int end=0;

	//Variables
	int i, j;

	double maxdist[kMAX];
	for(i=0; i<k; i++)
		maxdist[i] = MAXDOUBLE;



	double BestAggDist[128];
	bool status[128];

	//SG L[g_size-sg_size+1][kMAX]; //result list

	for(i=0; i<g_size-sg_size+1; i++){
		BestAggDist[i] = MAXDOUBLE;
		status[i] = false; //temporary; true-> final
		for(j=0; j<k; j++){
			//need to initialize
			//SG sg;
			//L[i][j] = sg; //RESULT List

			L[i][j].sgDist = MAXDOUBLE;
			L[i][j].m = sg_size+i; //subgrop size
			L[i][j].o = -1;
			L[i][j].gDist = MAXDOUBLE;
			L[i][j].status = false;


		}

	}

	Point2D x;
	GEOCENTROID(R,g_size,x);
	//calculate maximum distance to query points
	float maxQDist = MAXGROUPDIST(R,g_size,x);


	//Find the MBR from provided rectangles
	Rectangle1 mbr;
	mbr.x1=R[0].x1;
	mbr.x2=R[0].x2;
	mbr.y1=R[0].y1;
	mbr.y2=R[0].y2;

	for(i=1; i<g_size; i++)
	{
		if(R[i].x1 < mbr.x1)	mbr.x1=R[i].x1;
		if(R[i].x2 > mbr.x2)	mbr.x2=R[i].x2;
		if(R[i].y1 < mbr.y1)	mbr.y1=R[i].y1;
		if(R[i].y2 > mbr.y2)	mbr.y2=R[i].y2;
	}
	//end


	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------

	int son = root; //this entry is to be visited next

	//son == -1 means empty && end = 0 means exit

	while (son != -1 && end==0)
	{


		RTNode *rtn = new RTNode(this, son);

		for (i = 0; i < rtn -> num_entries; i++)
		{
			double edist1=0;//edist2=0, tdist=0;
			//edist1 = MINRECTDIST(mbr, rtn->entries[i].bounces);

			Pointlocation xp;
			xp.x = x[0];
			xp.y = x[1];

			Rectangle1 nd;



			nd.x1 = rtn->entries[i].bounces[0];
			nd.x2 = rtn->entries[i].bounces[1];
			nd.y1 = rtn->entries[i].bounces[2];
			nd.y2 = rtn->entries[i].bounces[3];

			if(rtn->level == 0){ //child node is a object
				nd.x2 = nd.x1;
				nd.y2 = nd.y1;
			}

			edist1 = MINRECTDIST1(xp, nd); // mindist from x to node

			/*
			if((edist1)>maxdist[k-1])	continue;
			edist1 = MINRECTDIST(R[0], rtn->entries[i].bounces);
			for (j=1; j < g_size; j++)
			{
				tdist = MINRECTDIST(R[j], rtn->entries[i].bounces);
				if(tdist>edist1)
					edist1=tdist;
			}

			if(edist1>maxdist[k-1])	continue;

			edist2 = MAXRECTDIST(R[0], rtn->entries[i].bounces);
			for (j=1; j < g_size; j++)
			{
				tdist= MAXRECTDIST(R[j], rtn->entries[i].bounces);
				if(tdist>edist2)
					edist2=tdist;
			}

			*/

			//update maxdistk
			/*
			for(int l=0; l<k; l++)
			{
				if (edist2 < maxdist[l])
				{
					for(int j=k-1; j>l; j--)
					{
						maxdist[j]=maxdist[j-1];
					}
					maxdist[l]=edist2;
					break;
				}
			}
			*/

			HeapEntry *h = new HeapEntry();
			h -> key = edist1;
			//he -> key1 = edist2;
			h -> level = rtn -> level;
			h -> son1 = rtn->entries[i].son;


			h-> x1 = nd.x1; //rtn->entries[i].bounces[0];
			h-> x2 = nd.x2; //rtn->entries[i].bounces[1]; //UNCOMMENTED BY EUNUS ; WAS commented IN sHETU
			h-> y1 = nd.y1; //rtn->entries[i].bounces[2];
			h-> y2 = nd.y2;//rtn->entries[i].bounces[3]; //UNCOMMENTED BY eUNUS
			heap -> insert(h);
			delete h;

		}

		delete rtn;


		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;

			bool isempty = !heap->remove(he);

			Rectangle1 p;
			p.x1 = he-> x1;
			p.y1 = he-> y1;
			p.x2 = he-> x2;
			p.y2 = he-> y2;




			bool isterminate = false;

			if(tb_lb ==1){ //Tight bound  method
				isterminate = !TERMINATE_TB(R, g_size,sg_size,x,p,maxQDist,f,BestAggDist,status);
			}else{ //lower bound method
				isterminate = !TERMINATE_LB(R, g_size,sg_size,x,p,maxQDist,f,BestAggDist,status);
			}

			if (isempty)  //heap is empty
				son = -1;
			else if (isterminate)//(he->key>maxdist[k-1]) //check terminating condition
			{
				//bool TERMINATE_LB(Rectangle1 R[], int g_size, int sg_size, Point2D x, Rectangle1 p, float maxQDist, int f, double BestAggDist[], bool status[])

				end = 1;

			}
			else
			{

				if (he->level == 0) //p is an object ; code for data object
				{

					//enter into result set
						//if(*num_of_data==MAXDATALIMIT)
							//printf("\nGreater than 10000\n");
							//error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);
						/*
						rslt[*num_of_data].x=he->x1;
						rslt[*num_of_data].y=he->y1;
						rslt[*num_of_data].dmin=he->key;
						rslt[*num_of_data].dmax=he->key1;
						*/

						*num_of_data=*num_of_data+1;
						//if(*num_of_data%1000==0) printf("\n%d",*num_of_data);


						//Eunus Code

						Heap *heapQ = new Heap();
						heapQ->init(dimension);
						double totalDist = 0.0;
						//float p[2];
						//p[0] = he->x1;
						//p[1] = he->y1;

						Point2D dp; //data point
						dp[0] = he->x1;
						dp[1] = he->y1;

						for(i=0; i<g_size;i++){

							HeapEntry *hq = new HeapEntry();

							Point2D qp; //query point
							qp[0] = R[i].x1;
							qp[1] = R[i].y1;

							hq -> key = Dist(dp,qp);//edist1; distance between q and p
							//he -> key1 = edist2;
							//he -> level = rtn -> level;
							hq -> son1 = i;//rtn->entries[i].son;
							hq-> x1 = qp[0];//rtn->entries[i].bounces[0];
							//he-> x2 = rtn->entries[i].bounces[1];
							hq-> y1 = qp[1];//rtn->entries[i].bounces[2];
							//he-> y2 = rtn->entries[i].bounces[3];
							heapQ -> insert(hq);


							if(i==0){
								totalDist = hq -> key;
							}else{
								totalDist = fDIST(totalDist,hq -> key,f);
							}

							delete hq;
						} //end of creating heapQ for query based on the distnace from data point

						//popping elements from HeapQ

						double aggDist = 0.0;


						char sgList[129];

						for(i=0;i<129;i++)
							sgList[i] = '0';
						sgList[128] = '\0';

						i = 0;
						HeapEntry *hq = new HeapEntry();
						while (heapQ->remove(hq))
						{
							if(i==0){
								aggDist = hq -> key;
							}else{
								aggDist = fDIST(aggDist,hq -> key,f);
							}

							sgList[hq -> son1] = '1';

							i++;
							if(i >= sg_size && BestAggDist[i-sg_size] > aggDist){

								SG sgi;
								sgi.m = i;
								sgi.o = he -> son1;
								sgi.sgDist = aggDist;
								sgi.gDist = totalDist;
								sgi.status = false;


								//sgi.sgList[64] = '\0';
								strcpy(sgi.sgList,sgList);

								//check with a exisiting subgroup and remove if it exists:

								/*
								int isduplicate = 0;
								for(int l=0; l<k; l++)
								{
									if(strcmp(L[i-sg_size][l].sgList,sgi.sgList)==0){
										if(L[i-sg_size][l].sgDist > sgi.sgDist){ //remove existing and procede
											for(int j = l; j<k; j++){ //shift everythign one cell up
												L[i-sg_size][j].sgDist = L[i-sg_size][j+1].sgDist;
												L[i-sg_size][j].m = L[i-sg_size][j+1].m;
												L[i-sg_size][j].o = L[i-sg_size][j+1].o;

												L[i-sg_size][j].gDist = L[i-sg_size][j+1].gDist;
												L[i-sg_size][j].status = L[i-sg_size][j+1].status;

												strcpy(L[i-sg_size][j].sgList,L[i-sg_size][j+1].sgList);

											}


										}else{ // not a candidate - continue while
											isduplicate = 1;
										}
										break;
									}
								}

								if(isduplicate ==1 ) continue;
								*/

								//order the result list RL and update BestAggDist


								for(int l=0; l<k; l++)
								{


									if (sgi.sgDist < L[i-sg_size][l].sgDist)
									{
										for(int j=k-1; j>l; j--)
										{
											//maxdist[j]=maxdist[j-1];
											L[i-sg_size][j].sgDist = L[i-sg_size][j-1].sgDist;
											L[i-sg_size][j].m = L[i-sg_size][j-1].m;
											L[i-sg_size][j].o = L[i-sg_size][j-1].o;
											//L[i-1][j].aggDist = L[i-1][j-1].aggDist;
											L[i-sg_size][j].gDist = L[i-sg_size][j-1].gDist;
											L[i-sg_size][j].status = L[i-sg_size][j-1].status;

											strcpy(L[i-sg_size][j].sgList,L[i-sg_size][j-1].sgList);

										}
										//maxdist[l]=edist2;
										L[i-sg_size][l].sgDist = sgi.sgDist;
										L[i-sg_size][l].m = sgi.m;
										L[i-sg_size][l].o = sgi.o;
										L[i-sg_size][l].gDist = sgi.gDist;
										L[i-sg_size][l].status = sgi.status;
										strcpy(L[i-sg_size][l].sgList,sgi.sgList);

										//update BestAggDist
										//BestAggDist[i-sg_size] = sgi.sgDist;

										break;
									} //if
								} //for

								//update BestAggDist with the kth best
								BestAggDist[i-sg_size] = L[i-sg_size][k-1].sgDist;

							} //if



						} //while

						delete hq;
						//get next data  from heap
						again =true;
						delete heapQ;
				}
				else //not leaf node
				{
					son=he->son1;
				}
			}
		}

		delete he;
	}
	delete heap;

}



void RTree::MQ_consensus_kGNN(Rectangle1 R[], int g_size, int sg_size, int k, int f, SG L[][16], int *num_of_data)
{


	//variable initialization
	int end=0;

	//Variables
	int i, j;

	double maxdist[kMAX];
	for(i=0; i<k; i++)
		maxdist[i] = MAXDOUBLE;



	double BestAggDist[128];
	bool status[128];

	//SG L[g_size-sg_size+1][kMAX]; //result list

	for(i=0; i<g_size-sg_size+1; i++){
		BestAggDist[i] = MAXDOUBLE;
		status[i] = false; //temporary; true-> final
		for(j=0; j<k; j++){
			//need to initialize
			//SG sg;
			//L[i][j] = sg; //RESULT List

			L[i][j].sgDist = MAXDOUBLE;
			L[i][j].m = sg_size+i; //subgrop size
			L[i][j].o = -1;
			L[i][j].gDist = MAXDOUBLE;
			L[i][j].status = false;


		}

	}

	Point2D x;
	GEOCENTROID(R,g_size,x);
	//calculate maximum distance to query points
	float maxQDist = MAXGROUPDIST(R,g_size,x);


	//Find the MBR from provided rectangles
	Rectangle1 mbr;
	mbr.x1=R[0].x1;
	mbr.x2=R[0].x2;
	mbr.y1=R[0].y1;
	mbr.y2=R[0].y2;

	for(i=1; i<g_size; i++)
	{
		if(R[i].x1 < mbr.x1)	mbr.x1=R[i].x1;
		if(R[i].x2 > mbr.x2)	mbr.x2=R[i].x2;
		if(R[i].y1 < mbr.y1)	mbr.y1=R[i].y1;
		if(R[i].y2 > mbr.y2)	mbr.y2=R[i].y2;
	}
	//end


	int rtind[128];
	float smq[128]; // mindist from q to explored area for each q
	bool isemp[128];


	Heap *heap[128];
	//init a heap that stores the non-leaf entries to be accessed-
	for(i=0; i<g_size; i++)
	{
		heap[i] = new Heap();
		heap[i]->init(dimension);
		smq[i] = 0.0;
		rtind[i] = root;
		isemp[i] =false;
	}

	int hind = 0; //heap index


	//------------------------------------------------------------

	int son = rtind[hind];//root; //this entry is to be visited next

	//son == -1 means empty && end = 0 means exit

	//while (son != -1 && end==0)
	//{

	//read children of root nodes


		RTNode *rtn = new RTNode(this, son);

		for (i = 0; i < rtn -> num_entries; i++)
		{
			double edist1=0;//edist2=0, tdist=0;
			//edist1 = MINRECTDIST(mbr, rtn->entries[i].bounces);

			Pointlocation xp;
			xp.x = x[0];
			xp.y = x[1];

			Rectangle1 nd;



			nd.x1 = rtn->entries[i].bounces[0];
			nd.x2 = rtn->entries[i].bounces[1];
			nd.y1 = rtn->entries[i].bounces[2];
			nd.y2 = rtn->entries[i].bounces[3];

			if(rtn->level == 0){ //child node is a object
				nd.x2 = nd.x1;
				nd.y2 = nd.y1;
			}

			//insert nodes to every heap
			for(j=0;j<g_size; j++)
			{

				//every query point
				xp.x = R[j].x1;
				xp.y = R[j].y1;

				edist1 = MINRECTDIST1(xp, nd); // mindist from x to node



				HeapEntry *h = new HeapEntry();
				h -> key = edist1;
				//he -> key1 = edist2;
				h -> level = rtn -> level;
				h -> son1 = rtn->entries[i].son;


				h-> x1 = nd.x1; //rtn->entries[i].bounces[0];
				h-> x2 = nd.x2; //rtn->entries[i].bounces[1]; //UNCOMMENTED BY EUNUS ; WAS commented IN sHETU
				h-> y1 = nd.y1; //rtn->entries[i].bounces[2];
				h-> y2 = nd.y2;//rtn->entries[i].bounces[3]; //UNCOMMENTED BY eUNUS

				heap[j] -> insert(h);
				delete h;
			}

		}

		delete rtn;

		//end of exploring the childrens of root nodes




		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		end = 0;

		while (!end)
		{
			//again = false;

			bool isempty = !heap[hind]->remove(he);

			Rectangle1 p;
			p.x1 = he-> x1;
			p.y1 = he-> y1;
			p.x2 = he-> x2;
			p.y2 = he-> y2;


			if (isempty){  //heap is empty
				son = -1; // should be all empty????
				isemp[hind] = true;
				end = 1; //assume all empty
				for(int t=0; t<g_size; t++)
				{
					if(isemp[t] ==false){
						end = 0;
					}
				}
				//one heap is empty: go to next
				hind = (hind+1) % g_size;

			}
			else if (!MQ_TERMINATE_LB(R, g_size,sg_size,smq,maxQDist,f,BestAggDist,status))//(he->key>maxdist[k-1]) //check terminating condition
			{
				//bool TERMINATE_LB(Rectangle1 R[], int g_size, int sg_size, Point2D x, Rectangle1 p, float maxQDist, int f, double BestAggDist[], bool status[])

				end = 1;

			}
			else
			{

				if (he->level == 0) //p is an object ; code for data object
				{


					//enter into result set
						//if(*num_of_data==MAXDATALIMIT)
							//printf("\nGreater than 10000\n");
						//	error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);

						*num_of_data=*num_of_data+1;
						//if(*num_of_data%1000==0) printf("\n%d",*num_of_data);


						//Eunus Code

						Heap *heapQ = new Heap();
						heapQ->init(dimension);
						double totalDist = 0.0;
						//float p[2];
						//p[0] = he->x1;
						//p[1] = he->y1;

						Point2D dp; //data point
						dp[0] = he->x1;
						dp[1] = he->y1;


						//update explored known area : to be used in the TERMINATE_LB

							Point2D cqp; //current query point
							cqp[0] = R[hind].x1;
							cqp[1] = R[hind].y1;
							smq[hind] = Dist(dp,cqp);
						//end updating expored area
							//if(smq[0] > 840.21){
								//printf("bestagg: %f ",BestAggDist[0]);
							//}

						for(i=0; i<g_size;i++){

							HeapEntry *hq = new HeapEntry();

							Point2D qp; //query point
							qp[0] = R[i].x1;
							qp[1] = R[i].y1;

							hq -> key = Dist(dp,qp);//edist1; distance between q and p
							//he -> key1 = edist2;
							//he -> level = rtn -> level;
							hq -> son1 = i;//rtn->entries[i].son;
							hq-> x1 = qp[0];//rtn->entries[i].bounces[0];
							//he-> x2 = rtn->entries[i].bounces[1];
							hq-> y1 = qp[1];//rtn->entries[i].bounces[2];
							//he-> y2 = rtn->entries[i].bounces[3];
							heapQ -> insert(hq);


							if(i==0){
								totalDist = hq -> key;
							}else{
								totalDist = fDIST(totalDist,hq -> key,f);
							}

							delete hq;
						} //end of creating heapQ for query based on the distnace from data point

						//popping elements from HeapQ

						double aggDist = 0.0;


						char sgList[129];


						for(i=0;i<129;i++)
							sgList[i] = '0';
						sgList[128] = '\0';

						i = 0;

						HeapEntry *hq = new HeapEntry();
						while (heapQ->remove(hq))
						{
							if(i==0){
								aggDist = hq -> key;
							}else{
								aggDist = fDIST(aggDist,hq -> key,f);
							}

							sgList[hq -> son1] = '1';

							i++;
							if(i >= sg_size && BestAggDist[i-sg_size] > aggDist){

								SG sgi;
								sgi.m = i;
								sgi.o = he -> son1;
								sgi.sgDist = aggDist;
								sgi.gDist = totalDist;
								sgi.status = false;


								//sgi.sgList[64] = '\0';
								strcpy(sgi.sgList,sgList);

								//check with a exisiting subgroup and remove if it exists:


								//order the result list RL and update BestAggDist


								for(int l=0; l<k; l++)
								{
									//duplicate check
									if(sgi.sgDist ==  L[i-sg_size][l].sgDist && sgi.o ==  L[i-sg_size][l].o)
										break;

									if (sgi.sgDist < L[i-sg_size][l].sgDist)
									{
										for(int j=k-1; j>l; j--)
										{
											//maxdist[j]=maxdist[j-1];
											L[i-sg_size][j].sgDist = L[i-sg_size][j-1].sgDist;
											L[i-sg_size][j].m = L[i-sg_size][j-1].m;
											L[i-sg_size][j].o = L[i-sg_size][j-1].o;
											//L[i-1][j].aggDist = L[i-1][j-1].aggDist;
											L[i-sg_size][j].gDist = L[i-sg_size][j-1].gDist;
											L[i-sg_size][j].status = L[i-sg_size][j-1].status;

											strcpy(L[i-sg_size][j].sgList,L[i-sg_size][j-1].sgList);

										}
										//maxdist[l]=edist2;
										L[i-sg_size][l].sgDist = sgi.sgDist;
										L[i-sg_size][l].m = sgi.m;
										L[i-sg_size][l].o = sgi.o;
										L[i-sg_size][l].gDist = sgi.gDist;
										L[i-sg_size][l].status = sgi.status;
										strcpy(L[i-sg_size][l].sgList,sgi.sgList);

										//update BestAggDist
										//BestAggDist[i-sg_size] = sgi.sgDist;

										break;
									} //if
								} //for

								//update BestAggDist with the kth best
								BestAggDist[i-sg_size] = L[i-sg_size][k-1].sgDist;

							} //if



						} //while

						delete heapQ;
						//get next data  from heap
						//again =true;
					//next object from next heap
					hind = (hind+1) % g_size;


				}
				else //not leaf node
				{
					son=he->son1;

					RTNode *rtn = new RTNode(this, son);

					for (i = 0; i < rtn -> num_entries; i++)
					{
						double edist1=0;//edist2=0, tdist=0;
						//edist1 = MINRECTDIST(mbr, rtn->entries[i].bounces);

						Pointlocation xp;

						Rectangle1 nd;

						nd.x1 = rtn->entries[i].bounces[0];
						nd.x2 = rtn->entries[i].bounces[1];
						nd.y1 = rtn->entries[i].bounces[2];
						nd.y2 = rtn->entries[i].bounces[3];

						if(rtn->level == 0){ //child node is a object
							nd.x2 = nd.x1;
							nd.y2 = nd.y1;
						}

							//every query point
							xp.x = R[hind].x1;
							xp.y = R[hind].y1;

							edist1 = MINRECTDIST1(xp, nd); // mindist from x to node


							HeapEntry *h = new HeapEntry();
							h -> key = edist1;
							//he -> key1 = edist2;
							h -> level = rtn -> level;
							h -> son1 = rtn->entries[i].son;


							h-> x1 = nd.x1; //rtn->entries[i].bounces[0];
							h-> x2 = nd.x2; //rtn->entries[i].bounces[1]; //UNCOMMENTED BY EUNUS ; WAS commented IN sHETU
							h-> y1 = nd.y1; //rtn->entries[i].bounces[2];
							h-> y2 = nd.y2;//rtn->entries[i].bounces[3]; //UNCOMMENTED BY eUNUS

							heap[hind] -> insert(h);
							delete h;

					}

					delete rtn;


				}
			}
		}

		delete he;
	//}
		for(i=0;i<sg_size;i++)
			delete heap[i];

}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
int RTree::load_keywords_from_file(char * keywords_file)
{	
    string STRING;
	ifstream infile;
	//infile.open(keywords_file);
	infile.open ("fruits.txt");
	int i,j;
	i=j=0;
        while(!infile.eof()) // To get you all the lines.
        {
	        //getline(infile,STRING); // Saves the line in STRING.
			//infile>>STRING;
	        //cout<<STRING<<endl; // Prints our STRING.
	        vector_string.push_back(STRING);
	        i++;
        }
        while( j<i)
        {
           // cout<<vector_string[j++]<<endl;
        }
	infile.close();
	system ("pause");

	return 0;
}
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
/* Sarowar
	in rtree there are total 32 functions and 15 fields

	1# 3 constructors and 1 destructor
	2# CRUD -> 4 functions
				 load_root()
		CREATE : insert()
		UPDATE :
		DELETE : del_root()
				 del_entry()
	 3# left 24 = 32 - 8
		QUERY  : total - 13

		-----  QUERY - 3 -----
		rangeQuery()
		rect_win_query() - 2
		-----  NNQ - 5   -----
		BFNN()
		KMin()
		Rect_kNNQ()
		Point_BFN_NNQ()
		TPNN_TP()
		-----  KGNN - 6  -----
		private_kGNN_max()
		kGNN()
		Naive_Consensus_kGNN()
		FANN_consensus_kGNN()
		consensus_kGNN()
		MQ_consensus_kGNN()

	 4# left 11 = 24 - 13
		----- there r 3 functions with update - i don't know what they do ? -----
		read_header()
		write_header()

		update_rslt()

		UpdateCount()
		UpdateStatus()

----------------------------------------------
				FIELDS
----------------------------------------------
	1	cache
		file

		deletelist

		last_pair
		lastcnt

		------- loaded from header -----
		dimension
		num_of_dnodes
		num_of_data
		num_of_inodes
		root
		root_is_data
		------- loaded from header -----

		re_data_cands
		re_level

		root_ptr
		tpheap

*/
