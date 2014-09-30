/*rtnode.cpp
  this file implements class RTNode*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "rtnode.h"
#include "rtree.h"
#include "entry.h"
#include "../blockfile/blk_file.h"
#include "../blockfile/cache.h"
#include "../linlist/linlist.h"

extern int io_access;
//------------------------------------------------------------
RTNode::RTNode(RTree *rt)
  //use this ructor to create a new node on disk.
{
    char *b;
    int header_size;
    Entry * d;
    int i;

    my_tree = rt;
    dimension = rt->dimension;
    num_entries = 0;
	dirty = TRUE;
	//Sarowar
	bitmask_array_length = (int) (TOTAL_FRUIT_NUMBER / 32);
	//----------------------------------------------------------------------------------
    d = new Entry();
	d -> init_entry(dimension, NULL);
    header_size = sizeof(char) + sizeof(int);  // level + num_entries
    capacity = (rt -> file -> get_blocklength() - header_size) / d -> get_size();

	/*
		Sarowar debug : RTnode insert()
	*/
	printf("Sarowar debug : RTnode(Rtree *rt) : capacity=%d = (rt -> file -> get_blocklength()=%d - header_size=%d) / d -> get_size()=%d \n",capacity, rt -> file -> get_blocklength(), header_size, d -> get_size());
	printf("//getchar() : Press Enter 2 Continue.\n");
	//getchar();
	//---------------------------------

    delete d;

    entries = new Entry[capacity];
    for (i = 0; i < capacity; i++)
		entries[i].init_entry(dimension, rt);

	//assign a new block on the disk
    b = new char[rt -> file -> get_blocklength()];
    block = rt -> file -> append_block(b);
    delete [] b;
}
//------------------------------------------------------------
RTNode::RTNode(RTree *rt, int _block)//Rui Zhang 002
  //use this ructor to restore a node from the disk.
{
    char *b;
    int header_size;
    Entry * d;
    int i;
	// sarowar
	bitmask_array_length = (int) (TOTAL_FRUIT_NUMBER / 32);
	//----------------------------------------------------
    my_tree = rt;
    dimension = rt->dimension;
    num_entries = 0;
	dirty = FALSE;

    d = new Entry();
	d -> init_entry(dimension, NULL);
    header_size = sizeof(char) + sizeof(int);
    capacity = (rt -> file -> get_blocklength() - header_size) / d -> get_size();
    delete d;

	/*
		Sarowar debug : RTnode insert()
	*/
	printf("Sarowar debug : RTNode(RTree *rt, int _block) : capacity=%d = (rt -> file -> get_blocklength()=%d - header_size=%d) / d -> get_size()=20\n_block=%d  \n",capacity, rt -> file -> get_blocklength(), header_size,_block);
	printf("//getchar() : Press Enter 2 Continue.\n");
		//getchar();

    entries = new Entry[capacity];
    for (i = 0; i < capacity; i++)
		entries[i].init_entry(dimension, rt);

    block = _block;
    b = new char[rt -> file -> get_blocklength()];
    if (rt -> cache == NULL) // no cache
        rt -> file -> read_block(b, block);
    else
        rt -> cache -> read_block(b, block, rt);

    read_from_buffer(b);
    delete [] b;
}
//------------------------------------------------------------
RTNode::~RTNode()
{
    char *b;

    if (dirty)
    {
        b = new char[my_tree->file->get_blocklength()];
        write_to_buffer(b);

        if (my_tree->cache == NULL) // no cache
            my_tree->file->write_block(b, block);
        else
            my_tree->cache->write_block(b, block, my_tree);

        delete [] b;
    }

    delete [] entries;
}
//------------------------------------------------------------
/* Sarowar
	used before calling insert()
	it is mainly used for choosing best subtree depending on mbr
	even if any entry doesn't contain, it produces best result
	by checking minimum enlargement.also considers overlapping.
*/
int RTNode::choose_subtree(float *mbr)
{
    int i, j, follow, minindex, *inside, inside_count, *over;
    float *bmbr, old_o, o, omin, a, amin, f, fmin;

    inside_count = 0;
    inside = new int[num_entries];
    over = new int[num_entries];
	// checks if mbr falls within any entries of this rtnode.
	// inside_count keeps track of no of them
    for (i = 0; i < num_entries; i++)
    {
    	switch (entries[i].section(mbr))
    	{
        	case INSIDE:
        	    inside[inside_count++] = i;
        	    break;
        }
    }

    if (inside_count == 1)
        // Case 1: There is exactly one dir_mbr that contains mbr
    	follow = inside[0];
    else if (inside_count > 1)
    // Case 2: There are many dir_mbrs that contain mbr
    // choose the one with the minimum area
    {
    	fmin = MAXREAL;
    	for (i = 0; i < inside_count; i++)
    	{
    	    f = area(dimension, entries[inside[i]].bounces);
    	    if (f < fmin)
    	    {
    	    	minindex = i;
          		fmin = f;
       	    }
       	}
    	follow = inside[minindex];
    }
    else
    // Case 3: There are no dir_mbrs that contain mbr
    // choose the one for which insertion causes the minimun overlap if son_is_data
    // else choose the one for which insertion causes the minimun area enlargement
    {
       	if (level == 1) // son_is_data
    	{
            omin = MAXREAL;
    	    fmin = MAXREAL;
    	    amin = MAXREAL;
    	    for (i = 0; i < num_entries; i++)
    	    {
        		enlarge(dimension, &bmbr, mbr, entries[i].bounces);

        		// calculate area and area enlargement
        		a = area(dimension, entries[i].bounces);
        		f = area(dimension, bmbr) - a;

        		// calculate overlap before enlarging entry_i
        		old_o = o = 0.0;
        		for (j = 0; j < num_entries; j++)
        		{
        		    if (j != i)
        		    {
    			        old_o += overlap(dimension,
    					 entries[i].bounces,
    					 entries[j].bounces);
    			        o += overlap(dimension,
    				     bmbr,
    				     entries[j].bounces);
    		        }
    	        }
    	        o -= old_o;

    	        // is this entry better than the former optimum ?
    	        if ((o < omin) ||
    		    (o == omin && f < fmin) ||
    		    (o == omin && f == fmin && a < amin))
    	        {
    	       	    minindex = i;
        		    omin = o;
        		    fmin = f;
        		    amin = a;
        	    }
    	        delete [] bmbr;
    	    }
        }
        else // son is not a data node
        {
    	    fmin = MAXREAL;
    	    amin = MAXREAL;
    	    for (i = 0; i < num_entries; i++)
    	    {
    	        enlarge(dimension, &bmbr, mbr, entries[i].bounces);

    	        // calculate area and area enlargement
    	        a = area(dimension, entries[i].bounces);
    	        f = area(dimension, bmbr) - a;

    	        // is this entry better than the former optimum ?
    	        if ((f < fmin) || (f == fmin && a < amin))
    	        {
    	       	    minindex = i;
    		        fmin = f;
    	            amin = a;
    	        }
	            delete [] bmbr;
	        }
        }

    	follow = minindex;

    	dirty = TRUE;
    }

    delete [] inside;
    delete [] over;

    return follow;
}
//------------------------------------------------------------
R_DELETE RTNode::delete_entry(Entry *e)
{
	RTNode *succ;
	float *tmp;
	if (level > 0)
	{
		if (this == my_tree->root_ptr)
			//i.e. this is the root
		{
			for (int i = 0; i < num_entries; i++)
			{
				tmp = overlapRect(dimension, entries[i].bounces, e -> bounces);
				if (tmp != NULL)
				{
					delete [] tmp;
					succ = entries[i].get_son();
					R_DELETE del_ret;
					del_ret = succ -> delete_entry(e);
					if (del_ret != NOTFOUND)
					{
						switch (del_ret)
						{
						case NORMAL:

							float *mbr;
							mbr = succ -> get_mbr();
							memcpy(entries[i].bounces, mbr, sizeof(float) * 2 * dimension);
							dirty = true;
							delete [] mbr;

							delete entries[i].son_ptr;
							entries[i].son_ptr = NULL;

							return NORMAL;
							break;

						case ERASED:
							delete entries[i].son_ptr;
							entries[i].son_ptr = NULL;

							int j;
							for (j = i; j < num_entries - 1; j++)
								entries[j] = entries[j+1];
							for (j = num_entries - 1; j < capacity; j++)
								entries[j].son_ptr = NULL;

							num_entries--;

							dirty = true;
							return NORMAL;
							break;
						}
					}
				}
			}
			//Not found;
			return NOTFOUND;
		}
		else//is not root and not leaf
		{
			for (int i = 0; i < num_entries; i++)
			{
				tmp = overlapRect(dimension, entries[i].bounces, e -> bounces);
				if (tmp != NULL)
				{
					delete [] tmp;
					succ = entries[i].get_son();
					R_DELETE del_ret;
					del_ret = succ->delete_entry(e);
					if (del_ret != NOTFOUND)
					{
						switch (del_ret)
						{
						case NORMAL:

							float *mbr;
							mbr = succ -> get_mbr();
							memcpy(entries[i].bounces, mbr, sizeof(float) * 2 * dimension);
							dirty = true;
							delete [] mbr;

							entries[i].del_son();

							return NORMAL;
							break;

						case ERASED:

							entries[i].del_son();

							int j;
							for (j = i; j < num_entries - 1; j++)
								entries[j] = entries[j+1];
							for (j = num_entries - 1; j < capacity; j++)
								entries[j].son_ptr = NULL;

							num_entries--;

							dirty = true;
							delete succ;

							if (num_entries < (int)ceil(0.4 * capacity))
							{
								for (int j = 0; j < num_entries; j++)
								{
									Linkable *e;
									e = entries[j].gen_Linkable();
									my_tree -> deletelist -> insert(e);
								}

								my_tree -> num_of_inodes --;
								return ERASED;
							}
							else
								return NORMAL;
							break;
						}
					}
				}
			}
		}
	}
	else//it's a leaf
	{
		for (int i = 0; i < num_entries; i++)
		{
			if (entries[i] == (*e))
			{
				my_tree -> num_of_data --;

				for (int j = i; j < num_entries-1; j++)
					entries[j] = entries[j+1];

				num_entries--;
				dirty = true;

				if (this != my_tree -> root_ptr && num_entries < (int)ceil(0.4 * capacity))
				{
					for (int k = 0; k < num_entries; k++)
					{
						Linkable *en;
					    en = entries[k].gen_Linkable();
						en -> level = 0;
						my_tree -> deletelist -> insert(en);
					}

					my_tree -> num_of_dnodes --;
					return ERASED;
				}
				else
					return NORMAL;
			}
		}
		return NOTFOUND;
	}
}
//------------------------------------------------------------
void RTNode::enter(Entry *de)
  //note that de will be deleted after being entered.
{
    if (num_entries > (capacity-1))
        error("RTNode::enter: called, but node is full", TRUE);

    entries[num_entries] = *de;

    num_entries++;

	dirty = true;

    de->son_ptr = NULL;
	//sarowar
	/*
	FILE *fp=fopen("C:/CSGLBSQ/Datasets/experiment/enter_node_values.txt","a");
	fprintf(fp,"Parameter :Entry *de: %f %f %f %f \n",de->bounces[0],de->bounces[1],de->bounces[2],de->bounces[3]);
	fprintf(fp,"field: %lf %lf %lf %lf \n", entries[num_entries].bounces[0],entries[num_entries].bounces[1],entries[num_entries].bounces[2],entries[num_entries].bounces[3]);
	fprintf(fp,"-----------------------------------\n");
	fclose(fp);
	*///sarowar
    delete de;
}
//------------------------------------------------------------
/* Sarowar
	finds the leaf node that contains this entry starting from
	current node. checks all the entries of the current node in
	a loop and if an overlap is found checks recursively untill
	the leaf i.e. level == 0 node.
*/
bool RTNode::FindLeaf(Entry *e)
{
	RTNode *succ;
	if (level > 0)
	{
		for (int i = 0; i < num_entries; i++)
		{
			float *f;
			f = overlapRect(my_tree -> dimension,
				  entries[i].bounces, e -> bounces);
			if (f != NULL)
			{
				delete [] f;
				succ = entries[i].get_son();
				bool find;
				find = succ->FindLeaf(e);
				entries[i].del_son();
				if (find)
					return true;
			}
		}
		return false;
	}
	else
	{
		for (int i = 0; i < num_entries; i++)
		{
			if (entries[i] == (*e))
				return true;
		}
		return false;
	}
	return false;
}
//------------------------------------------------------------
/* Sarowar
	returns the minimum mbr of all the entries of a rtnode.
*/
float* RTNode::get_mbr()
{
    int i, j;
    float *mbr;

    mbr = new float[2*dimension];
    for (i = 0; i < 2*dimension; i ++ )
        mbr[i] = entries[0].bounces[i];

    for (j = 1; j < num_entries; j++)
    {
    	for (i = 0; i < 2*dimension; i += 2)
    	{
    	    mbr[i]   = min(mbr[i],   entries[j].bounces[i]);
    	    mbr[i+1] = max(mbr[i+1], entries[j].bounces[i+1]);
        }
    }

    return mbr;
}
/**
	Sarowar implementation of get_minimuboundingtime();  -- start
	implemented just as a similar copy of get_mbr() -- function above
*/
//--------------------------------------------------------------------------------------//

RTIME* RTNode::get_minimumboundingtime(){
	int i, j;
    //double *mbt;
	RTIME *mbt;

	mbt = new RTIME[RTIME_SIZE];
	for (i = 0; i < RTIME_SIZE; i ++ )
        mbt[i] = entries[0].times[i];

    for (j = 1; j < num_entries; j++)
    {
    	mbt[0]   = min(mbt[0],   entries[j].times[0]);
		mbt[1]	 = max(mbt[1], entries[j].times[1]);
		/*
		for (i = 0; i < 2*dimension; i += 2)
    	{
    	    mbt[i]   = min(mbr[i],   entries[j].bounces[i]);
    	    mbt[i+1] = max(mbr[i+1], entries[j].bounces[i+1]);
        }
		*/
    }

    return mbt;

}

int * RTNode::get_minimumFruitBitmask()
{
	int *fruit_bitmask;
	int size;
	size = (TOTAL_FRUIT_NUMBER / 32);
	int i,j;
	fruit_bitmask = new int[size];
	for (i = 0; i < size; i++ )
        fruit_bitmask[i] = 0;

    for (j = 0; j < num_entries; j++)
    {
		for (i = 0; i < size ; i ++)
    	{
			fruit_bitmask[i] |= entries[j].fruit_bitmask[i];
        }		
    }
    return fruit_bitmask;
}
//--------------------------------------------------------------------------------------
//------------------------------------------------------------
/* Sarowar
	returns total no of data under a rtnode. it root calls
	it will return total no.  of leaf nodes in the tree.
*/
int RTNode::get_num_of_data()
{
    int i, sum;
    RTNode* succ;

    if (level == 0)
        return num_entries;

    sum = 0;
    for (i = 0; i < num_entries ; i++)
    {
        succ = entries[i].get_son();
        sum += succ->get_num_of_data();
		entries[i].del_son();
    }

    return sum;
}
//------------------------------------------------------------
R_OVERFLOW RTNode::insert(Entry *d, RTNode **sn)
{
    int follow;
    RTNode *succ, *new_succ;
    RTNode *brother;
    Entry *de;
    R_OVERFLOW ret;
    float *mbr,*nmbr;
	//Sarowar : time
	RTIME *minimum_bounding_time;
	int *minimum_fruit_bitmask;
	//--------------

    int i, last_cand;
    float *center;
    SortMbr *sm;
    Entry *new_entries;

    if (level > 0) // direcrtory node
    {
	   /*	Sarowar debug : RTnode insert()		*/
		printf("\n\nSarowar debug : RTnode insert() : if( level=%d > 0 )=TRUE \n",level);
		printf("//getchar() : Press Enter 2 Continue.\n");
		//getchar();
	  if (level > d -> level)
	  {
		/*	Sarowar debug : RTnode insert()	*/
		  printf("\n\nSarowar debug : RTnode insert() : (level > d->level) : if( level=%d > d->level=%d )=TRUE \n", level, d->level);
		  printf("//getchar() : Press Enter 2 Continue.\n");
		  //getchar();

        follow = choose_subtree(d -> bounces);

        succ   = entries[follow].get_son();

		/*	Sarowar debug : RTnode insert()	*/
		printf("\n\nSarowar debug : choose_subtree() returns => follow=%d'th son_ptr (*rtnode) is used for recursive insert() call. \n",follow);
		printf("//getchar() : Press Enter 2 Continue.\n");
		//getchar();

        ret = succ -> insert(d, &new_succ);

        mbr = succ -> get_mbr(); // region bounding all the entries of succ is the mbr of this entry input.
        memcpy(entries[follow].bounces, mbr, sizeof(float) * 2 * dimension);
        delete [] mbr;

		// Sarowar : rtnode : insert() : 1
		minimum_bounding_time = succ -> get_minimumboundingtime();
		memcpy( entries[follow].times ,minimum_bounding_time , sizeof(RTIME) * RTIME_SIZE );
		delete [] minimum_bounding_time;

		minimum_fruit_bitmask		  = succ -> get_minimumFruitBitmask();
		//int size = entries[follow].bitmask_array_length;
		memcpy(entries[follow].fruit_bitmask,minimum_fruit_bitmask, bitmask_array_length * sizeof(int));
		delete [] minimum_fruit_bitmask;
		//-------------------------

		entries[follow].del_son();
		// 3 possible action : enum R_OVERFLOW {SPLIT, REINSERT, NONE};
        if (ret == SPLIT)
        // node has split into itself and *new_succ
        {
			/*
			Sarowar debug : RTnode insert()
			*/
			printf("\n\nSarowar debug : Rtnode insert() : SPLIT : if( ret == SPLIT )=TRUE \n");
			printf("//getchar() : Press Enter 2 Continue.\n");
			//getchar();

            if (num_entries == capacity)
         	    error("RTNode::insert: maximum capacity violation", TRUE);

            de = new Entry(dimension, my_tree);
    	    nmbr = new_succ -> get_mbr();
            memcpy(de -> bounces, nmbr, 2 * dimension * sizeof(float));
    	    delete [] nmbr;

			// Sarowar : rtnode : insert() : 1
			minimum_bounding_time = new_succ -> get_minimumboundingtime();
			memcpy( de->times ,minimum_bounding_time , sizeof(RTIME) * RTIME_SIZE );
			delete [] minimum_bounding_time;

			minimum_fruit_bitmask = new_succ -> get_minimumFruitBitmask();
			memcpy( de->fruit_bitmask , minimum_fruit_bitmask, bitmask_array_length * sizeof(int));
			delete [] minimum_fruit_bitmask;
			//-------------------------

            de -> son = new_succ -> block;
			delete new_succ;
            de -> son_ptr = NULL;
            enter(de);

			/*	Sarowar debug : RTnode insert()	*/
			//printf("Sarowar debug : entered Entry *de in Entry[] : de->son= %d \n de->bounces[0]= %f \n de->bounces[1]= %f \n de->bounces[2]= %f \n de->bounces[3]= %f \n .\n", (de -> son), (de->bounces[0]), (de->bounces[1]), (de->bounces[2]), (de->bounces[3]) );
			//printf("//getchar() : Press Enter 2 Continue.\n");
			////getchar();

            if (num_entries == (capacity - 1))
            {
        	    brother = new RTNode(my_tree);
        	    my_tree -> num_of_inodes++;
        	    brother -> level = level;
        	    split(brother);
                *sn = brother;
                ret = SPLIT;
        	}
            else
          	    ret = NONE;
        }
        dirty = TRUE;

        return ret;
	  }
	  else //level==d->level
	  {
		  /*	Sarowar debug : RTnode insert()	*/
		  printf("\n\nSarowar debug : RTnode insert() : level==d->level : if( level=%d == d->level=%d )=TRUE \n", level, d->level);
		  printf("Sarowar debug :entered Entry *d in Entry[] : d->son= %d \n d->bounces[0]= %f \n d->bounces[1]= %f \n d->bounces[2]= %f \n d->bounces[3]= %f \n d->times[0]=%lld \n d->times[1]=%lld \n\n before insert() from rtree constructor.\n", (d -> son), (d->bounces[0]), (d->bounces[1]), (d->bounces[2]), (d->bounces[3]), (d->times[0]), (d->times[1]) );
		  printf("//getchar() : Press Enter 2 Continue.\n");
		  //getchar();


		  enter(d);    //note that d will be deleted on return

		  if (num_entries == (capacity - 1))
            // maximun no of entries --> Split
            // this happens already if the node is nearly filled
            // for the algorithms are more easy then
		  {
			/*
			Sarowar debug : RTnode insert()
			*/
			  printf("\n\nSarowar debug : RTnode insert() : if(num_entries=%d == (capacity=%d - 1))=TRUE : return SPLIT : called SPLIT.\n", num_entries, capacity);
			  printf("//getchar() : Press Enter 2 Continue.\n");
			  //getchar();

        	brother = new RTNode(my_tree);
        	my_tree -> num_of_inodes++;
        	brother -> level = level;
        	split(brother);
            *sn = brother;
            ret = SPLIT;
		  }
          else
          	ret = NONE;

		  dirty=true;
		  return ret;
	  }
    }
    else // data (leaf) node i.e. level == 0
    {


        if (num_entries == capacity)
        	error("RTDataNode::insert: maximum capacity violation", TRUE);
		/*	Sarowar debug : RTnode insert()	*/
		printf("\n\nSarowar debug : Rtnode insert() : data (leaf) node : if (level=%d > 0)=false \n",level);
		printf("Sarowar debug :entered Entry *d in Entry[] : \nd->son=%d \n d->bounces[0]=%f \n d->bounces[1]=%f \n d->bounces[2]=%f \n d->bounces[3]=%f .\n d->times[0]=%lld \n d->times[1]=%lld \n\n", (d -> son), (d->bounces[0]), (d->bounces[1]), (d->bounces[2]), (d->bounces[3]),(d->times[0]),(d->times[1]) );
		printf("//getchar() : Press Enter 2 Continue.\n");
		//getchar();

        enter(d);



        dirty = TRUE;

		// during insertion in root and capacity reached.
        if ( num_entries == (capacity - 1) )
        // maximum # of entries --> Split
        // this happens already if the node is nearly filled
        // for the algorithms are more easy then
        {
			/*	Sarowar debug : RTnode insert()	*/
			  printf("\n\nSarowar debug : RTnode insert() : if(num_entries=%d == (capacity=%d - 1))=TRUE :.\n", num_entries, capacity);
			  printf("//getchar() : Press Enter 2 Continue.\n");
			  //getchar();
			  // this wasn't re_leveled && this level != root->level
            if (my_tree->re_level[0] == FALSE && my_tree -> root_ptr -> level != level)
    	    // there was no reinsert on level 0 during this insertion
			//-------------------------------------------------------
			//Here I changed the condition as if it is already root, no need
			//to reinsert.  Split directly in this case
			//-----------By TAO Yufei
            {
				/*
				Sarowar debug : RTnode insert()
					Here :
					1 mbr is calculated
					2 from mbr by averaging center x & y is got
					3 in SortMbr[] entries[] are kept
					4 a qSort on SortMbr[] is done
					5 according to SortMbr[i] after sorting, last_cand no. of data are entered in new_entry
					6 rest of entries are converted to linkable and inserted in to re_data_cand linklist.
					7 my_tree->re_level[0] i.e. for root is made true.dirty=true, num_entries-=last_cand
					8 returned REINSERT
				*/
				  printf("\n\nSarowar debug : RTnode insert() : if (my_tree->re_level[0] == FALSE && my_tree -> root_ptr -> level=%d != level=%d)=TRUE .\n returned REINSERT\n", my_tree -> root_ptr -> level, level);
				  printf("//getchar() : Press Enter 2 Continue.\n");
					//getchar();


                // calculate center of page
                mbr = get_mbr();
                center = new float[dimension];
                for (i = 0; i < dimension; i++)
                     center[i] = (mbr[2*i] + mbr[2*i+1]) / 2.0;

                new_entries = new Entry[capacity];

				for (i = 0; i < capacity; i ++)
					new_entries[i].init_entry(dimension, my_tree);

        	    sm = new SortMbr[num_entries];
        	    for (i = 0; i < num_entries; i++)
        	    {
            		sm[i].index = i;
            		sm[i].dimension = dimension;
            		sm[i].mbr = entries[i].bounces;
            		sm[i].center = center;
                }

                qsort(sm, num_entries, sizeof(SortMbr), sort_center_mbr);

                last_cand = (int) ((float)num_entries * 0.30);

                // copy the nearest candidates to new array
                for (i = 0; i < num_entries - last_cand; i++)
    	            new_entries[i] = entries[sm[i].index];

                // insert candidates into reinsertion list
                for ( ; i < num_entries; i++)
                {
					Linkable *nd = entries[sm[i].index].gen_Linkable();
                    my_tree -> re_data_cands -> insert(nd);
                }

                // free and copy data array
                delete [] entries;
        	    entries = new_entries;

        	    delete sm;
        	    delete [] mbr;
        	    delete [] center;
        	    my_tree -> re_level[0] = TRUE;

        	    // correct # of entries
        	    num_entries -= last_cand;

        	    // must write page
        	    dirty = TRUE;

                return REINSERT;
        	}
           	else  //there has been reinsertion on this level
           	{

				/*
				Sarowar debug : RTnode insert()
				*/
				printf("\n\nSarowar debug : RTnode insert() : if (my_tree->re_level[0] == FALSE && my_tree -> root_ptr -> level=%c != level=%c)=false : \n called split((RTNode *) *sn); return SPLIT \n", my_tree -> root_ptr -> level, level);
				printf("//getchar() : Press Enter 2 Continue.\n");
				//getchar();

        	    *sn = new RTNode(my_tree);
        	    (*sn) -> level = level;
        	    my_tree -> num_of_dnodes++;
        	    split((RTNode *) *sn);
    	    }
    	    return SPLIT;
        }
        else
            return NONE;
    }
}
//------------------------------------------------------------
void RTNode::print()
{
    int i;

	printf("level %d  Block: %d\n", level, block);

    for (i = 0; i < num_entries ; i++)
    {
        printf("(%4.1lf, %4.1lf, %4.1lf, %4.1lf)\n(%lld, %lld)\n",
	       entries[i].bounces[0],
	       entries[i].bounces[1],
	       entries[i].bounces[2],
	       entries[i].bounces[3],
		   entries[i].times[0],
		   entries[i].times[1]);
    }
}
//------------------------------------------------------------
void RTNode::rangeQuery(float *mbr, SortedLinList *res)
{
    int i, n;
    SECTION s;
    RTNode *succ;

    n = num_entries;
    for (i = 0; i < n; i++)
    {
        s = entries[i].section(mbr);
        if (s == INSIDE || s == OVERLAP)
        {
            if (level == 0)
            {
                Linkable *copy;
				copy = entries[i].gen_Linkable();
        		res -> insert(copy);
            }
            else
            {
                succ = entries[i].get_son();
				//Added by Tanzima
				io_access++;
				//..............
                succ -> rangeQuery(mbr, res);
				entries[i].del_son();
            }
        }
    }
}
//------------------------------------------------------------
//	Sarowar : debug	: 
void RTNode::rangeQuery_with_time(float *mbr,RTIME *mbt, SortedLinList *res)
{
	printf("\n inside : rangeQuery_with_time()\n");
    int i, n;
    SECTION s,t;
    RTNode *succ;

    n = num_entries;
    for (i = 0; i < n; i++)
    {
        s = entries[i].section(mbr);
        if (s == INSIDE || s == OVERLAP)
        {

			printf("\nrangeQuery_with_time :: Rectangle inside \nEnter 2 continue.\n");
			//getchar();

			t = entries[i].section_time(mbt);
			
			if (t == INSIDE || t == OVERLAP)
			{	
				printf("\nrangeQuery_with_time :: time inside \nEnter 2 continue.\n");
				//getchar();
				if (level == 0)
				{
					printf("\nrangeQuery_with_time :: leaf node. Enter to result. \nEnter 2 continue.\n");
					//getchar();
					Linkable *copy;
					copy = entries[i].gen_Linkable();
        			res -> insert(copy);
				}
				else
				{
					printf("\nrangeQuery_with_time :: Directory node. Recursive query. \nEnter 2 continue.\n");
					//getchar();
					succ = entries[i].get_son();
					//Added by Tanzima
					io_access++;
					//..............
					succ ->rangeQuery_with_time(mbr, mbt, res);
					//succ -> rangeQuery(mbr, res);
					entries[i].del_son();
					printf("\nrangeQuery_with_time :: Recursive query. Ends. \nEnter 2 continue.\n");
					//getchar();
				}
			}
        }
    }
}

//------------------------------------------------------------
void RTNode::read_from_buffer(char *buffer)
{
    int i, j, s;

    // Level
    memcpy(&level, buffer, sizeof(char));
    j = sizeof(char);

    // num_entries
    memcpy(&num_entries, &(buffer[j]), sizeof(int));
    j += sizeof(int);
	//----------------
	printf("\nRTnode :load_root() : RTnode(rtree, block) : read_from_buffer() :num_entries=%d.\nenter 2 continue.\n\n",num_entries);
	//getchar();
	//----------------
    s = entries[0].get_size();
    for (i = 0; i < num_entries; i++)
    {
    	entries[i].read_from_buffer(&buffer[j]);
    	j += s;
    }
}
//------------------------------------------------------------
int RTNode::split(float **mbr, int **distribution)
{
    bool lu;
    int i, j, k, l, s, n, m1, dist, split_axis;
    SortMbr *sml, *smu;
    float minmarg, marg, minover, mindead, dead, over, *rxmbr, *rymbr;

    n = num_entries;

    m1 = (int) ceil((float)n * 0.40);

    sml = new SortMbr[n];
    smu = new SortMbr[n];
    rxmbr = new float[2*dimension];
    rymbr = new float[2*dimension];

    // choose split axis
    minmarg = MAXREAL;
    for (i = 0; i < dimension; i++)
    // for each axis
    {
        for (j = 0; j < n; j++)
        {
            sml[j].index = smu[j].index = j;
            sml[j].dimension = smu[j].dimension = i;
            sml[j].mbr = smu[j].mbr = mbr[j];
        }

        // Sort by lower and upper value perpendicular axis_i
      	qsort(sml, n, sizeof(SortMbr), sort_lower_mbr);
        qsort(smu, n, sizeof(SortMbr), sort_upper_mbr);

        marg = 0.0;
        // for all possible distributions of sml
        for (k = 0; k < n - 2 * m1 + 1; k++)
        {
			for (s = 0; s < 2 * dimension; s += 2)
			{
				rxmbr[s] =    MAXREAL;
				rxmbr[s+1] = -MAXREAL;
			}
            for (l = 0; l < m1 + k; l++)
            {
				for (s = 0; s < 2*dimension; s += 2)
				{
					rxmbr[s]   =   min(rxmbr[s],   sml[l].mbr[s]);
					rxmbr[s+1] = max(rxmbr[s+1], sml[l].mbr[s+1]);
				}
			}
			marg += margin(dimension, rxmbr);

			for (s = 0; s < 2 * dimension; s += 2)
			{
				rxmbr[s]   =    MAXREAL;
				rxmbr[s+1] = -MAXREAL;
			}
            for ( ; l < n; l++)
            {
				for (s = 0; s < 2 * dimension; s += 2)
				{
					rxmbr[s] =   min(rxmbr[s],   sml[l].mbr[s]);
					rxmbr[s+1] = max(rxmbr[s+1], sml[l].mbr[s+1]);
				}
            }
			marg += margin(dimension, rxmbr);
        }

        // for all possible distributions of smu
       	for (k = 0; k < n - 2 * m1 + 1; k++)
        {
            // now calculate margin of R1
			// initialize mbr of R1
			for (s = 0; s < 2 * dimension; s += 2)
			{
				rxmbr[s] =    MAXREAL;
				rxmbr[s+1] = -MAXREAL;
			}
            for (l = 0; l < m1+k; l++)
            {
                // calculate mbr of R1
				for (s = 0; s < 2 * dimension; s += 2)
				{
					rxmbr[s] =   min(rxmbr[s],   smu[l].mbr[s]);
					rxmbr[s+1] = max(rxmbr[s+1], smu[l].mbr[s+1]);
				}
            }
			marg += margin(dimension, rxmbr);

            // now calculate margin of R2
			// initialize mbr of R2
			for (s = 0; s < 2 * dimension; s += 2)
			{
				rxmbr[s] =    MAXREAL;
				rxmbr[s+1] = -MAXREAL;
			}
            for ( ; l < n; l++)
            {
                // calculate mbr of R1
				for (s = 0; s < 2 * dimension; s += 2)
				{
					rxmbr[s] =   min(rxmbr[s],   smu[l].mbr[s]);
					rxmbr[s+1] = max(rxmbr[s+1], smu[l].mbr[s+1]);
				}
            }
			marg += margin(dimension, rxmbr);
        }

        if (marg < minmarg)
        {
            split_axis = i;
            minmarg = marg;
        }
    }

    // choose best distribution for split axis
    for (j = 0; j < n; j++)
    {
		sml[j].index = smu[j].index = j;
		sml[j].dimension = smu[j].dimension = split_axis;
		sml[j].mbr = smu[j].mbr = mbr[j];
    }

    // Sort by lower and upper value perpendicular split axis
    qsort(sml, n, sizeof(SortMbr), sort_lower_mbr);
    qsort(smu, n, sizeof(SortMbr), sort_upper_mbr);

    minover = MAXREAL;
    mindead = MAXREAL;
    // for all possible distributions of sml and snu
    for (k = 0; k < n - 2 * m1 + 1; k++)
    {
        dead = 0.0;
		for (s = 0; s < 2 * dimension; s += 2)
		{
			rxmbr[s] =    MAXREAL;
			rxmbr[s+1] = -MAXREAL;
		}
		for (l = 0; l < m1 + k; l++)
		{
			for (s = 0; s < 2*dimension; s += 2)
			{
				rxmbr[s] =   min(rxmbr[s],   sml[l].mbr[s]);
				rxmbr[s+1] = max(rxmbr[s+1], sml[l].mbr[s+1]);
			}
			dead -= area(dimension, sml[l].mbr);
		}
        dead += area(dimension, rxmbr);
		  //**************note**************
		  //this does not compute the dead space for all the cases.  some overlapping
		  //area may be subtrated twice.
		  //********************************

		for (s = 0; s < 2 * dimension; s += 2)
		{
			rymbr[s] =    MAXREAL;
       		rymbr[s+1] = -MAXREAL;
		}
		for ( ; l < n; l++)
		{
			for (s = 0; s < 2*dimension; s += 2)
			{
				rymbr[s] =   min(rymbr[s],   sml[l].mbr[s]);
				rymbr[s+1] = max(rymbr[s+1], sml[l].mbr[s+1]);
			}
			dead -= area(dimension, sml[l].mbr);
		}
        dead += area(dimension, rymbr);

		over = overlap(dimension, rxmbr, rymbr);

        if ((over < minover) || (over == minover) && dead < mindead)
        {
            minover = over;
            mindead = dead;
            dist = m1+k;
            lu = TRUE;
        }

		//Now we do the same thing for smu
        dead = 0.0;
		for (s = 0; s < 2*dimension; s += 2)
		{
			rxmbr[s] =    MAXREAL;
			rxmbr[s+1] = -MAXREAL;
		}
		for (l = 0; l < m1+k; l++)
		{
			for (s = 0; s < 2*dimension; s += 2)
			{
				rxmbr[s] =   min(rxmbr[s],   smu[l].mbr[s]);
				rxmbr[s+1] = max(rxmbr[s+1], smu[l].mbr[s+1]);
			}
			dead -= area(dimension, smu[l].mbr);
		}
        dead += area(dimension, rxmbr);

		for (s = 0; s < 2*dimension; s += 2)
		{
			rymbr[s] =    MAXREAL;
			rymbr[s+1] = -MAXREAL;
		}
		for ( ; l < n; l++)
		{
			for (s = 0; s < 2*dimension; s += 2)
			{
				rymbr[s] =   min(rymbr[s],   smu[l].mbr[s]);
				rymbr[s+1] = max(rymbr[s+1], smu[l].mbr[s+1]);
			}
			dead -= area(dimension, smu[l].mbr);
		}
		//correcting errors
        dead += area(dimension, rymbr);

		over = overlap(dimension, rxmbr, rymbr);

        if ((over < minover) || (over == minover) && dead < mindead)
        {
            minover = over;
            mindead = dead;
            dist = m1+k;
            lu = FALSE;
        }
    }

    // calculate best distribution
	// the array distribution is deleted in split(rtnode *sn);
    *distribution = new int[n];
    for (i = 0; i < n; i++)
    {
        if (lu)
            (*distribution)[i] = sml[i].index;
        else
            (*distribution)[i] = smu[i].index;
    }

    delete [] sml;
    delete [] smu;
    delete [] rxmbr;
    delete [] rymbr;

    return dist;
}
//------------------------------------------------------------
void RTNode::split(RTNode *sn)
{
    int i, *distribution, dist, n;
    float **mbr_array;
    Entry *new_entries1, *new_entries2;

    n = num_entries;

    mbr_array = new floatptr[n];
    for (i = 0; i < n; i++)
       	mbr_array[i] = entries[i].bounces;

    dist = split(mbr_array, &distribution);

    new_entries1 = new Entry[capacity];
    new_entries2 = new Entry[capacity];

	for (i = 0; i < capacity; i ++)
	{
		new_entries1[i].init_entry(dimension, my_tree);
		new_entries2[i].init_entry(dimension, my_tree);
	}

    for (i = 0; i < dist; i++)
       	new_entries1[i] = entries[distribution[i]];

    for (i = dist; i < n; i++)
       	new_entries2[i-dist] = entries[distribution[i]];

    for (i = 0; i < n; i++)
    {
       	entries[i].son_ptr = NULL;
       	sn->entries[i].son_ptr = NULL;
    }
    delete [] entries;
    delete [] sn->entries;

    entries = new_entries1;
    sn->entries = new_entries2;

    num_entries = dist;
    sn->num_entries = n - dist;

    delete [] mbr_array;
	delete [] distribution;
}
//------------------------------------------------------------
void RTNode::write_to_buffer(char *buffer)
{
    int i, j, s;

    // Level
    memcpy(buffer, &level, sizeof(char));
    j = sizeof(char);

    // num_entries
    memcpy(&buffer[j], &num_entries, sizeof(int));
    j += sizeof(int);

    s = entries[0].get_size();
    for (i = 0; i < num_entries; i++)
    {
    	entries[i].write_to_buffer(&buffer[j]);
       	j += s;
    }
}

//------------------------------------------------------------
void RTNode::rect_win_query(float *mbr, LinList *in_objs, LinList *out_objs_so_far)
{
    int i, n;
    SECTION s;
    RTNode *succ;

    n = num_entries;
    for (i = 0; i < n; i++)
    {
		s = section_new(2, mbr, entries[i].bounces);

		if (level == 0)
        {
			if (s == INSIDE || s == OVERLAP)
			{
				Linkable *copy;
				copy = entries[i].gen_Linkable();
        		in_objs -> insert(copy);
	        }

/*  we remarked this because we assume there is no buffer available to put the outer objects found so far
			else
			{
				Linkable *copy;
				copy = entries[i].gen_Linkable();
        		out_objs_so_far -> insert(copy);
			}
*/
		}
        else
        {
			if (s == INSIDE || s == OVERLAP)
			{
				succ = entries[i].get_son();
				succ -> rect_win_query(mbr, in_objs, out_objs_so_far);
				entries[i].del_son();
	        }
        }
    }
}

//------------------------------------------------------------
void RTNode::rect_win_query(float *mbr, float *exclmbr, LinList *c_inf_objs)
{
    int i, n;
    SECTION s, s1;
    RTNode *succ;

    n = num_entries;
    for (i = 0; i < n; i++)
    {
		s = section_new(2, mbr, entries[i].bounces);
		s1 = section_new(2, exclmbr, entries[i].bounces);

		if (level == 0)
        {
			if ((s == INSIDE || s == OVERLAP) && s1 != INSIDE && s1 != OVERLAP)
			{
				Linkable *copy;
				copy = entries[i].gen_Linkable();
        		c_inf_objs -> insert(copy);
	        }
		}

/* we remarked this because we assume there is no buffer available to put the outer objects found so far
		else if (level == 1)
        {
			if ((s == INSIDE || s == OVERLAP) && s1 != INSIDE && s1 != OVERLAP)
			{
				succ = entries[i].get_son();
				succ -> rect_win_query(mbr, exclmbr, c_inf_objs);
				entries[i].del_son();
	        }
        }
*/
		else
        {
			if ((s == INSIDE || s == OVERLAP) && s1 != INSIDE)
			{
				succ = entries[i].get_son();
				succ -> rect_win_query(mbr, exclmbr, c_inf_objs);
				entries[i].del_son();
	        }
        }
    }
}


/* Sarowar
	in rtree there are total 19 functions and 8 fields

	1# 2 constructors and 1 destructor
	2# CRUD -> 3 functions

		CREATE : insert()
				 enter()
		UPDATE :
		DELETE : delete_entry()
	 3# left 13 = 19 - 6
		QUERY  : total - 2

		-----  QUERY - 3 -----

		rangeQuery()
		rect_win_query() - 2


	 4# left 6 = 13 - 7

		choose_subtree()
		FindLeaf()
		get_mbr()
		get_num_of_data()
		print()
		split() - 2

	 5# - 2

		read_from_buffer()
		write_to_buffer()


-----------------------------------------------------
					FIELD
-----------------------------------------------------
//--===on disk===--
	char level;
	int block;
	int num_entries;
	Entry *entries;
//--===others===--
	bool dirty;
	int capacity;
    int dimension;
	RTree *my_tree;

*/
