#pragma once
#include <iostream>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <cstdio>

using namespace std;

class Rtree_time
{
public:
	//-------------------
	struct tm *time_tm;
	//-------------------
	Rtree_time(void);
	~Rtree_time(void);
	int  time_input_from_file(FILE *&fp, time_t *rt1,time_t *rt2);
	struct tm *int64_to_tm( time_t time_int64);
	void print_time_onscreen(time_t time_int64);
	void print_time_onscreen(struct tm  struct_tm);
};
