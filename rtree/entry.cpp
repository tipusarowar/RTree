/* entry.cpp
   implementation of class Entry */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "entry.h"
#include "rtnode.h"
#include "../linlist/linlist.h"
//------------------------------------------------------------
Entry::Entry()
  //this constructor does nothing.  remember you should call init_entry
  //to initialize if you use this ructor
{
	son_ptr = NULL;
	bounces = NULL;
	// Sarowar -- # 
	times   = NULL;
	fruit_bitmask = NULL;
	//-------------
	nn_num=-1;
}
//------------------------------------------------------------
Entry::Entry(int _dimension, RTree *rt)
{
    dimension = _dimension;
    my_tree = rt;
    bounces = new float[2*dimension];
	// --===Sarowar===-- #
	times = new RTIME[RTIME_SIZE];
	bitmask_array_length = (int)TOTAL_FRUIT_NUMBER / 32 ;
	fruit_bitmask = new int[ bitmask_array_length ];
	//-------------
    son_ptr = NULL;
    son = 0;
	level = 0;
	nn_num=-1;
}
//------------------------------------------------------------
Entry::~Entry()
{
    if (bounces)
		delete [] bounces;
    if (son_ptr != NULL)
	    delete son_ptr;
	// Sarowar -- #
	if( times )
		delete [] times;
	if( fruit_bitmask )
		delete [] fruit_bitmask;
	// ------------
}
//------------------------------------------------------------
void Entry::del_son()
{
	if (son_ptr != NULL)
	{
		delete son_ptr;
		son_ptr = NULL;
	}
}
//------------------------------------------------------------
Linkable* Entry::gen_Linkable()
{
	Linkable *new_link = new Linkable(dimension);
	new_link -> son = son;
	//memcpy(new_link -> bounces, bounces, 2 * dimension * sizeof(float));
	for (int i = 0; i < 2 * dimension; i ++)
		new_link -> bounces[i] = bounces[i];
	new_link -> level = level;

	/*	Sarowar time : start : 2	*/
	new_link -> times[0] = times[0];
	new_link -> times[1] = times[1];

	int size = (TOTAL_FRUIT_NUMBER /32 );
	for(int i = 0 ; i< bitmask_array_length ; i++)
	{
		new_link ->fruit_bitmask[i]=fruit_bitmask[i];
	}
	/*	Sarowar time : end	*/
	return new_link;
}
//------------------------------------------------------------
int Entry::get_size()
{
    //return 2 * dimension * sizeof(float) + sizeof(int);
	  //for bounces and son
	//  before fruit_bimask was added.
	//return 2 * dimension * sizeof(float) + sizeof(int) + ( 2 * sizeof(RTIME)) ;	 
	/*	Sarowar time : start : 4
		printf("time add : 4\n");	*/
				// bounces + son + time + fruit_bitmask
	int size = TOTAL_FRUIT_NUMBER / 32 ;
	return 2 * dimension * sizeof(float) + sizeof(int) + ( 2 * sizeof(RTIME)) + (bitmask_array_length * sizeof(int)) ;	 

	/*	Sarowar time : end : 4	*/
}
//------------------------------------------------------------
RTNode* Entry::get_son()
{
	// Sarowar debug 
	printf("\n\n**************************\n");
	printf("Sarowar : RTNode* Entry::get_son() : son=%d\n ",son);
    
	if (son_ptr == NULL)
	    son_ptr = new RTNode(my_tree, son);
	// disk access is made
	if( son_ptr == NULL )
		printf("son_ptr == NULL");
	else
		printf("son_ptr != NULL");
	printf("\n**************************\n");
    return son_ptr;
}
//------------------------------------------------------------
void Entry::init_entry(int _dimension, RTree *_rt)
{
	dimension = _dimension;
    my_tree = _rt;
    bounces = new float[2 * dimension];
	// Sarowar -- #
	times   = new RTIME[RTIME_SIZE];
	int size  = (TOTAL_FRUIT_NUMBER / 32);
	bitmask_array_length = (TOTAL_FRUIT_NUMBER / 32);
	fruit_bitmask = new int[bitmask_array_length];
	//-------------
    son_ptr = NULL;
    son = 0;
	level = 0;
	//sarowar : debug
	printf("void Entry::init_entry(int _dimension, RTree *_rt) ---- ENDS \n");
}
//------------------------------------------------------------
void Entry::read_from_buffer(char *buffer)
{
    int i;

    i = 2 * dimension * sizeof(float);
    memcpy(bounces, buffer, i);

    memcpy(&son, &buffer[i], sizeof(int));
    i += sizeof(int);

	/*	Sarowar time : start : 5
		printf("time add : 5\n");	*/

	int size = RTIME_SIZE * sizeof(RTIME);
	memcpy(times, &buffer[i], size);
	i+=size;
	// fruit_bitmask
	//byte 
	size = bitmask_array_length * sizeof(int);
	memcpy(fruit_bitmask, &buffer[i],size );
	i+=size;
	/*	Sarowar time : end : 5	*/
}
//------------------------------------------------------------
SECTION Entry::section(float *mbr)
{
    bool inside;
    bool overlap;

    overlap = TRUE;
    inside = TRUE;

    for (int i = 0; i < dimension; i++)
    {
		// Sarowar : mbr.x1 > bounces.x2 || mbr.x2 < bounces.x1
		// Sarowar : mbr.y1 > bounces.y2 || mbr.y2 < bounces.y1
		if (mbr[2 * i] > bounces[2 * i + 1] ||  mbr[2 * i + 1] < bounces[2 * i])
			overlap = FALSE;
		if (mbr[2 * i] < bounces[2 * i] ||
			mbr[2 * i + 1] > bounces[2 * i + 1])
			inside = FALSE;
    }
    if (inside)
		return INSIDE;
    else if (overlap)
		return OVERLAP;
    else
		return S_NONE;
}
// Sarowar :: ----  ####### ---------
//------------------------------------------------------------
SECTION Entry::section_time(RTIME *mbt)
{
	/*
	bool inside,overlap,within_mbt;
    within_mbt = inside = overlap = false;
	*/	
	bool none=false;
	RTIME mbt_t1,mbt_t2,entry_t1,entry_t2;
	mbt_t1 = mbt[0];
	mbt_t2 = mbt[1];
	entry_t1 = times[0];
	entry_t2 = times[1];	
	/*
	if( entry_t1 <= mbt_t1  && mbt_t2 <= entry_t2 ){inside = true;}
	//else if( ( mbt_t1 >= entry_t1 && mbt_t2 > entry_t2 ) || ( mbt_t1 <  entry_t1 && mbt_t2 <= entry_t2 ) )//here a mistake is done: case mbt_t1 & mbt_t2 both less than e_t1 & e_t2 will be shown as inside.{overlap = true;}
	else if( (mbt_t1 <= entry_t1 && mbt_t2 >= entry_t1) || ((mbt_t1 <= entry_t2 && mbt_t2 >= entry_t2)) )
	{// not implemented. what will it return upon return.//within_mbt = true;overlap = true;}	
	if(inside)return INSIDE;else if(overlap)return OVERLAP;else return S_NONE;
	*/
	if( (mbt_t1 < entry_t1 && mbt_t2 < entry_t1) || (mbt_t1 > entry_t2 && mbt_t2 > entry_t2) )
	{
		none = true;
	}
	if(none==true)
		return S_NONE;
	else
		return OVERLAP;
}


SECTION Entry::section_keywords(int *keywords)
{
	if(		 (fruit_bitmask[0] & keywords[0]) != 0)
		return OVERLAP;
	else if( (fruit_bitmask[1] & keywords[1]) != 0)
		return OVERLAP;
	else if( (fruit_bitmask[2] & keywords[2]) != 0)
		return OVERLAP;
	else if( (fruit_bitmask[3] & keywords[3]) != 0)
		return OVERLAP;
	else 
		return S_NONE;
			
}
//------------------------------------------------------------
bool Entry::section_circle(float *center, float radius)
{
	float r2;

	r2 = radius * radius;

	if ((r2 - MINDIST(center,bounces,dimension)) < FLOATZERO)
		return TRUE;
	else
		return FALSE;
}
//------------------------------------------------------------
void Entry::set_from_Linkable(Linkable *link)
{
	son = link -> son;
	dimension = link -> dimension;
	memcpy(bounces, link -> bounces, 2 * dimension * sizeof(float));
	level = link -> level;
	
	/*	Sarowar time : start : 3	*/
	 
	memcpy(times, link -> times , RTIME_SIZE * sizeof(RTIME) );
	
	memcpy(fruit_bitmask, link->fruit_bitmask, bitmask_array_length * sizeof(int));
	/*	Sarowar time : end : 3	*/

	my_tree = NULL;
	son_ptr = NULL;
}
//------------------------------------------------------------
void Entry::write_to_buffer(char *buffer)
{
    int i;

    i = 2 * dimension * sizeof(float);
    memcpy(buffer, bounces, i);

    memcpy(&buffer[i], &son, sizeof(int));
    i += sizeof(int);

	/*	Sarowar time : start : 6
		printf("time add : 6\n");	*/
	int size = RTIME_SIZE * sizeof(RTIME);
	memcpy(&buffer[i], times, size);
	i += size;	

	size = bitmask_array_length * sizeof(int);
	memcpy(&buffer[i], fruit_bitmask, size);
	i+=size;
	/*	Sarowar time : end : 6	*/
}
//------------------------------------------------------------
bool Entry::operator == (Entry &_d)
  //this function compares two entries based on (1)son (2)dimension (3)extents
{
	if (son != _d.son) return false;
	if (dimension != _d.dimension) return false;
	for (int i = 0; i < 2 * dimension; i++)
		if (fabs(bounces[i] - _d.bounces[i]) > FLOATZERO) return false;
	/*	Sarowar time : start : 7

		printf("time add : 7\n");
	*/
	for(int i = 0; i < 2 ; i++){
		//if( _abs64( times[i] - _d.times[i] ) > FLOATZERO ) return false;
		if( ( times[i] - _d.times[i] ) > FLOATZERO  ) return false;
		else if( ( _d.times[i] - times[i] ) > FLOATZERO  ) return false;
	}	

	for(int i = 0 ; i< bitmask_array_length ; i++)
	{
		if( fruit_bitmask[i] != _d.fruit_bitmask[i] )
			return false;
	}
	/*	Sarowar time : end : 7	*/
	return true;
}
//------------------------------------------------------------
Entry& Entry::operator = (Entry &_d)
  //this function assigns all fieds of _d with the same values of this entry
{
    dimension = _d.dimension;
    son = _d.son;
    son_ptr = _d.son_ptr;
    memcpy(bounces, _d.bounces, sizeof(float) * 2 * dimension);
    my_tree = _d.my_tree;
	level = _d.level;

	nn_num=_d.nn_num;

	/*	Sarowar time : start : 7

		printf("time add : 7\n");
	*/
	memcpy(times, _d.times , sizeof(RTIME) * RTIME_SIZE );

	memcpy(fruit_bitmask, _d.fruit_bitmask, sizeof(int) * bitmask_array_length);
	
	/*	Sarowar time : end : 7	*/

    return *this;
}
//------------------------------------------------------------