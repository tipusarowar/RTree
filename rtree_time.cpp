#include "rtree_time.h"

Rtree_time::Rtree_time(void)
{
}

Rtree_time::~Rtree_time(void)
{
}


int Rtree_time::time_input_from_file(FILE *&fp, time_t *rt1,time_t *rt2){
	
	cout<<"\n\n------------ Function starts -------------"<<endl;
    int a;
    struct tm date1,date2;
    time_t t1,t2;
    char s,p,z,x;s=p='-';z=x=':';
	
	// FORMAT : YY MM DD HH MIN SEC
    fscanf( fp," %d %d %d %d %d %d ",&a,&date1.tm_mon,&date1.tm_mday,&date1.tm_hour,&date1.tm_min,&date1.tm_sec);
	cout<<"1st y-m-d = "<<a<<"-"<<date1.tm_mon<<"-"<<date1.tm_mday<<" || h-m-s = "<<date1.tm_hour<<" : "<<date1.tm_min<<" : "<<date1.tm_sec<<endl;
	
	date1.tm_year = a - 1900;
	t1 = mktime( &date1 );
	cout<<"t1 ="<<t1<<endl;
	*rt1=t1;
	
	// debug : conversion & reverse conversion works :
    /*	struct tm *newDate1 = localtime( &t1 );cout<<"Reverse y m d ="<<newDate1->tm_year<<" "<<newDate1->tm_mon << " " << newDate1->tm_mday <<" : h m s =" << newDate1->tm_hour<<" "<<newDate1->tm_min<<" "<<newDate1->tm_sec<<endl;time_t new_t1 = mktime( newDate1 );cout<<"new_t1 = "<<new_t1<<endl;
	*/
	// FORMAT : YY MM DD HH MIN SEC
    fscanf( fp," %d %d %d %d %d %d ",&a,&date2.tm_mon,&date2.tm_mday,&date2.tm_hour,&date2.tm_min,&date2.tm_sec);
	cout<<"2nd y-m-d = "<<a<<"-"<<date2.tm_mon<<"-"<<date2.tm_mday<<" || h-m-s = "<<date2.tm_hour<<":"<<date2.tm_min<<":"<<date2.tm_sec<<endl;

    date2.tm_year = a - 1900;
	t2 = mktime( &date2);
	*rt2=t2;
    cout<<"t2 = "<<t2<<endl;

    cout<<"-------------Function ends----------\n\n "<<endl;
	return 1;
}



struct tm *Rtree_time::int64_to_tm(time_t time_int64)
{

	return localtime( &time_int64 );
}

void Rtree_time::print_time_onscreen(time_t time_int64)
{
	struct tm *struct_tm = localtime( &time_int64);
	cout<<"y-m-d="<<struct_tm->tm_year<<"-"<<struct_tm->tm_mon<<"-"<<struct_tm->tm_mday<<" || h:m:s = "<<struct_tm->tm_hour<<" : "<<struct_tm->tm_min<<" : "<<struct_tm->tm_sec<<endl;
	cout<<"time __int64 format:"<<time_int64<<endl;

}

void Rtree_time::print_time_onscreen(tm struct_tm)
{
	time_t time_int64 = mktime( &struct_tm);
	cout<<"y-m-d="<<struct_tm.tm_year<<"-"<<struct_tm.tm_mon<<"-"<<struct_tm.tm_mday<<" || h:m:s = "<<struct_tm.tm_hour<<" : "<<struct_tm.tm_min<<" : "<<struct_tm.tm_sec<<endl;
	cout<<"time __int64 format:"<<time_int64<<endl;
}

