#ifndef __LINLIST
#define __LINLIST

#include <stdio.h>
#include <ctime>
#include "../func/gendef.h"
////////////////////////////////////////////////////////////////////////
// LinList  (SLink)
////////////////////////////////////////////////////////////////////////
// Sarowar  --- ##  ----
typedef time_t RTIME;
static int RTIME_SIZE = 2;
//---------------------

class RTree;

struct Linkable
{
public:
	int son;
	int dimension;
	int level;
	float *bounces;
	/*	Sarowar time : start : 8

		printf("time add : 8\n");	*/
	//double *times;
	RTIME   *times;
	int		*fruit_bitmask;
	//int fruit_bitmask[4];
	/*	Sarowar time : end : 8	*/
	float distanz;

	Linkable(int dim)
	{ dimension = dim;
	  bounces = new float[2 * dim];
	  // Sarowarn  --- ## ---
	  times	  = new RTIME[RTIME_SIZE];
	  int size = TOTAL_FRUIT_NUMBER/32;
	  fruit_bitmask = new int[size];
	  //---------------------
    }

	~Linkable()
	{
		delete [] bounces;
		// Sarowar --- ## ---
		delete [] times;
		delete [] fruit_bitmask;
		// ------------------
	}

	void clone(Linkable *old);
};

struct SLink
{
    Linkable *d;          // Zeiger auf Element-Daten
    SLink *next;          // Zeiger auf naechstes Element
    SLink *prev;          // Zeiger auf vorhergehendes Element

    SLink();
    ~SLink();
};

////////////////////////////////////////////////////////////////////////
// LinList
////////////////////////////////////////////////////////////////////////


class LinList
{
protected:
    SLink *first;         // Rootzeiger des Datenbestands
    SLink *last;          // Zeiger auf letztes Element
    int anz;                    // Anzahl der belegten Elemente in der Liste
    SLink *akt;           // zeigt auf aktuelles Element
    int akt_index;              // Index des zuletzt mit get geholten Elements
public:
    LinList();
    virtual ~LinList();
    int get_num()               // gibt Anzahl der im Index belegten Elements
        { return anz; }         // zurueck

    void check();               // ueberprueft Konsistenz der Liste
    void print();
	void printnew();

    void insert(Linkable *f);       // haengt ein Element vorne an die Liste an
    bool erase();               // loescht aktuelles Element aus der Liste

    Linkable * get(int i);          // liefert i-tes Element
    Linkable * get_first();         // liefert erstes Element im Index
    Linkable * get_last();          // liefert erstes Element im Index
	Linkable * get_lastnew();		//added by Shirley, get the last linkable but the current poistion is unchanged
	Linkable * get_next();          // liefert naechstes Element im Index
    Linkable * get_prev();          // liefert vorhergehendes Element im Index
};


////////////////////////////////////////////////////////////////////////
// SortedLinList
////////////////////////////////////////////////////////////////////////


class SortedLinList : public LinList
{
    bool increasing;
public:
    SortedLinList();

    void set_sorting(bool _increasing); // wenn increasing gleich TRUE, wird
                                // aufsteigend einsortiert
                                // DIESE FUNKTION MUSS VOR DEM ERSTEN EINFUEGEN
                                // GERUFEN WERDEN !!!!!!!!!!!
    void insert(Linkable *f);       // fuegt ein Element durch direktes Einfuegen ein

    void sort(bool _increasing);// sortiert die Liste
};

#endif  // __LINLIST
