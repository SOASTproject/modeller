#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <math.h>

// in getparcel locs ONLY SETDOWN ALLOWED
#define REALMAXLOCS 100000
#define MAXLOCS 10000
#define MAXPACKAGES 20000
#define MAXBUSFLEETS 20

#define OUTSIDE 123
#define INSIDE  234

#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)

typedef struct point {
  double x,y;
} POINT;

POINT *points;
int npoints = 0;



double HDMULTIPLIER = 1.2;

float tmat[MAXLOCS][MAXLOCS];
float dmat[MAXLOCS][MAXLOCS];

#define __BUS 0
#define __PARCEL 1
#define __LOCKERPARCEL 11
#define __LASTMILE 15
#define __BUSPLUS 20
#define __PARCELPLUS 30
#define __BUSDEPOT 334

float parcel_speedfactor = 1, bus_speedfactor = 0.9;
float _speedfactor;  // speed in comparison to parcel van

float busdepot_arrival = 540;  // when bundle carried by parcel van must get to bus depot
float bundle_available = 560; //  latest that a parcel bundle gets to the relevant bus depot. So, bus can pick it up then or later
float bundle_drop = 900; //  latest that a parcel bundle can get to the village drop point, for lastmile pckup
float bundle_available_lastmile = 930; //  latest that a parcel bundle can get to the village drop point, for lastmile pckup

int _parcel_vans;

int bus_maxparcels = 20, bus_maxpassengers = 24, parcelvan_maxparcels = 120;

int _parcelvan_maxparcels = 120;

int doing = 0;

/****************************************************************************************************\

 location types

  0  bus can set down
  1  bus can pick up
  2  bus can set down or pick up
  3  bus passenger pickup raw
  4  bus passenger setdown raw
 23  parcel locker for pickup
 24  raw setdown position for parcel
 28  raw setdown for parcel that starts at parcel depot
  
  parcel types

  -1  no parcel activity
  10 parcel setdown
  20 parcel pickup
  30 parcel pickup and setdown

\****************************************************************************************************/

typedef struct location {
  float lat, lon;
  int type, partype;
  char str[100];

  int bef;
  float v, w;  // useful to have these here ??

  int nparcels;

  float arrtime, nolater, noearlier;


  //loc temporarily holds lastmile fleet details
  int nvans ; // or whatever other vehicle
  int maxparcels;
  float maxmins, earlystart, loadtime, speedfactor, permile, perhour;



} LOCATION;


int nlocs = 0, slocs = 0; // slocs are the snapped locs - the matrix will only use them
LOCATION locs[REALMAXLOCS];


typedef struct package {

  int pu_raw, pu_snap, sd_raw, sd_snap; // actual PU, snapped PU, 
  float pu_open, pu_close, sd_open, sd_close;
  int inplay;
  int type; // 0 bus   1 parcel

  int sd_assigned;
  int nparcels;
  int pu_bindex, sd_bindex;
  
  float v, w, arrtime, nolater, noearlier;
  int id;
  int isbundle;
  int bundle_index;

  int bus_pickup_lindex;
  int bus_setdown_lindex;

} PACKAGE;

PACKAGE packages[MAXPACKAGES];
int npackages = 0;

typedef struct busdemand {

  int bindex ;// index into bloc
  int before;  // which busdemand it must be before
  float preftime;// preferred pickup time
  float latesttime;
  int type; //  0 means pickup / 1 is setdown
  int nparcels;
} BUSDEMAND;


BUSDEMAND *busdemands;
int nbusdemand;

typedef struct depot {

  int loc;
  int type; 

} DEPOT;
DEPOT *depots;
int ndepots = 0;

int addloc(int t, int pt, float lat,float lon, char *str);
int addbusdem(float lat, float lon, char *tstr, float lat2, float lon2);
double lldistance(double lat1, double lon1, double lat2, double lon2, char unit);
double deg2rad(double deg);
double rad2deg(double rad);
double pi = 3.141592653589793;


float bus_earlystart = 360, bus_latestart = 900, bus_maxmins = 480;
float bus_permile = 5, bus_perhour = 15;


float _permile, _perhour;

float passenger_early = 0, passenger_late = 20;

float lastmile_permile = 4, lastmile_perhour = 12, lastmile_speed = 5;

float parcel_lat, parcel_lon;
int parcel_vans;
float parcel_w, parcel_v;
float parcel_permile, parcel_perhour, parcel_visit_time = 1, parcel_maxmins = 480, parcel_earlystart = 360, parcel_latestart = 540, bus_stop_time = 1;
float parcel_loadtime = 30, bundle_visit_time = 2;


float _parcel_permile, _parcel_perhour, _parcel_visit_time = 1, _parcel_maxmins = 480, _parcel_earlystart = 360, _parcel_latestart = 540;
float _parcel_loadtime = 30;

typedef struct busfleet {
  float lat, lon;
  int n;

  int lindex, bindex; // locations of the fleet in locs and buslocs respectively

} BUSFLEET;
int nbusfleets = 0;

BUSFLEET busfleets[MAXBUSFLEETS];
char lastmile[100];
int addbusfleet(int n, float lat, float lon);
int addobpdem(float lat, float lon, float v, float w);


LOCATION *uselocs;
int nuselocs;

float parcelvan_speed = 35; // mph


typedef struct vanplan {

  int ndels;
  int *sequence;
  double objectives[5];
  
} VANPLAN;

typedef struct busplan {

  int fleet; // determines start and end point - also times etc
  int ndels;
  int *sequence;
  double objectives[5];
  
} BUSPLAN;




typedef struct candidate {

  int nvanplans, nbusplans;
  VANPLAN *vanplans;
  BUSPLAN *busplans;  
  int maxon;

  double objectives[5];
  
} CANDIDATE;
CANDIDATE *population;

int popsize = 20;


typedef struct best3 {

  float q[3];  // lower the better
  int i1[3], i2[3], i3[3], i4[3];
  int n;
} BEST3;

float addity(CANDIDATE *cp, int ploc, int vi, int posi, int type);
double travtime(int a, int b);
int applyq3(float q,int i1,int i2,int i3,int i4,int to, BEST3 *b3);
int addbest3(float q, int i1, int i2, int i3, int i4, BEST3 *b3);


int parcel_iterations = 3000000;
int bus_iterations = 3000000;
int closest_busloc(float lat, float lon);

// parcel lastmile site - bus drops here, lastmile picks up here
typedef struct sdlocker {
  int lindex, bindex, fleet, nassigned, id;

  // the lastmile fleet parameters at this site
  
  int nvans ; // or whatever other vehicle
  int maxparcels;
  float maxmins, earlystart, loadtime, speedfactor, permile, perhour;

  
} SDLOCKER;

int nsdlockers = 0, npulockers = 0;

SDLOCKER *sdlockers, *pulockers;


int allfleet = 0;
float addity_bus(CANDIDATE *cp, int bpu, int bsd, int vi, int posi, int posi2);

float lastmile_radius = 1; // miles


typedef struct bundle {

  int sdlocker;
  int contents[100];
  int ncontents;

  int bus_pu_lindex, bus_sd_lindex, van_sd_lindex; 
} BUNDLE;

BUNDLE bundles[1000];
int nbundles = 0;

int uselocs_alloced = 0;
int ndeliveries = 0;
double travdist(int a, int b);
char mapfilebase[100];
FILE *mf;

char areapoly[50];
char pclls[50];

int gotareapoly = 0, gotpclls = 0;

int nparcellocs = 0;

int lms_density = 0;

double *aplats, *aplons;
int naplls = 0;
int apllmalloc = 1000;
int napllmalloced = 0;


double *pclats, *pclons;
int pcmalloc = 1000;
int npcmalloced = 0;
int npcs = 0;

double filter_pclls(void);
void add_apll(double lat, double lon);
void add_pcll(double lat, double lon);

int num_lms = 0;

typedef struct lms_pop {
  int *centres;
  double f;
} LMS_POP;
LMS_POP *lms_pop;
int lms_popsize = 20;

char *pcllarray;


int nclusters = 0;
int inferred_lms = 0;

  typedef struct gridcell {
    double lon, lat;  // top left;
    int npoints;
    int x, y;
  } GRIDCELL;

void main(int argc, char **argv)
{
  int arg = 1;

  getbuslocs(argv[arg++]);  // potential places the bus could stop to pickup or setdown

  getparcellocs(argv[arg++]);  // places where a bus can drop parcels for later take-away
                               // (and/OR could contain 'automated' signal
  
  getbusdem(argv[arg++]);    // normal bus passenger demand - pickups and setdowns
  getobpdem(argv[arg++]);    // parcels only
  getconfig(argv[arg++]);    // config

  strcpy(mapfilebase,argv[arg++]);

  printf("... picked up base name for mapfiles:  %s\n", mapfilebase);
  printf("... arg is %d nd argc is %d\n", arg, argc);
  
  if(argc > arg)
    {
      printf("... expecting area polygon ... \n");
      strcpy(areapoly,argv[arg++]); gotareapoly = 1;
      printf("... got it:  %s\n", areapoly);
    }
  if(argc > arg)
    {
      printf("... expecting full set of UK postcode latlons ... \n");
      strcpy(pclls,argv[arg++]); gotpclls = 1;
      printf("... hopefully that's them in %s \n", pclls);
    }  

  printlocs(10);
  printparams();

  if(lms_density > 0)
    {
      identify_lms();

    }


  
  organise();

}


int identify_lms(void)
{

  pclats = (double *)malloc(pcmalloc*sizeof(double));
  pclons = (double *)malloc(pcmalloc*sizeof(double));  
  npcmalloced = pcmalloc;
  
  // get the data
  printf("starting to identify lastmile sites\n");
  
  get_areapoly();

  print_areapoly();

  filter_pclls();
  printf("there are %d postcodes in the region\n", npcs);
  
  find_lms_centres();
  //  exit(0);
  
}


int find_lms_centres(void)
{

  // first, grid up the area so that it contains num_lms squares
  num_lms = npcs / lms_density;

  printf("there will be %d lastmile sites\n", num_lms);
  

  // get bounding box

  int i, j;
  double lonlow, lonhigh, latlow, lathigh;
  lonlow = lonhigh = pclons[0];
  latlow = lathigh = pclats[0];
  
  for(i=1;i<npcs;i++)
    {
      if(pclons[i] < lonlow) lonlow = pclons[i];
      if(pclons[i] > lonhigh) lonhigh = pclons[i];
      if(pclats[i] < latlow) latlow = pclats[i];
      if(pclats[i] > lathigh) lathigh = pclats[i];                  
    }

  // add s small margin
  lonlow -= 0.0001;
  lonhigh += 0.0001;
  latlow -= 0.0001;
  lathigh += 0.0001;  
  
  // now grid this box until we have num_lms cluster centres
  double ndivs = 3;



  GRIDCELL **grid = (GRIDCELL **)malloc(1000 * sizeof(GRIDCELL *));
  for(i=0;i<1000;i++) grid[i] = (GRIDCELL *)malloc(1000 * sizeof(GRIDCELL));
  double lonsize, latsize, lat, lon;
  GRIDCELL *gc;
  int x, y;
  
  while(1)
    {
      // set up the grid with ndivs points
      // what is the lon cell size?
      lonsize = (lonhigh - lonlow)/ndivs;
      latsize = (lathigh - latlow)/ndivs;      
      
      for(i=0;i<ndivs;i++)
	for(j=0;j<ndivs;j++)
	  {
	    gc = &grid[i][j];
	    gc->x = i; gc->y = j;
	    gc->lon = lonlow + i * lonsize;
	    gc->lat = latlow + j * latsize;	    
	    gc->npoints = 0;
	  }

      // now assign the points

      for(i=0;i<npcs; i++)
	{
	  lat = pclats[i];  lat -= latlow;
	  lon = pclons[i];  lon -= lonlow;
	  x = (int) floor (lon / lonsize);
	  y = (int) floor (lat / latsize);	  
	  
	  gc = &grid[x][y];
	  gc->npoints++;	  
	}
      // how many cells are occupied?

      int occ = 0;
      for(i=0;i<ndivs;i++)
	for(j=0;j<ndivs;j++)
	  {
	    gc = &grid[i][j];
	    if(gc->npoints > 0) occ++;
	  }     
      printf("with %f x %f cells (%f), occupied are: %d  (num_lms: %d)\n", ndivs, ndivs, ndivs*ndivs, occ, num_lms);

      ndivs+=1.0;
      if(ndivs > 1000) break;
      if(occ >= num_lms) break;
      
    }

  // okey dokey,  the occupied parts of the current grid now represent our lastmile regions
  // we will take the centremost pcll of each one  to represent it.

  system("rm __lmsites.txt");

  for(i=0;i<ndivs;i++)
    for(j=0;j<ndivs;j++)
      {
	gc = &grid[i][j];
	if(gc->npoints > 0)
	  {
	    establish_lmspoint(gc,latlow,lonlow,latsize,lonsize);
	  }
      }

  //  exit(0);
}

int establish_lmspoint(GRIDCELL *gc, double latlow, double lonlow,double latsize,double lonsize)
{

  double clats[1000], clons[10000], lat, lon, clat, clon;

  int i, x, y, n=0;

  for(i=0;i<npcs;i++)
    {
      if(gc->npoints < 1) continue;
      
      lat = pclats[i];  lat -= latlow;
      lon = pclons[i];  lon -= lonlow;
      x = (int) floor (lon / lonsize);
      y = (int) floor (lat / latsize);

      if(x!=gc->x) continue;
      if(y!=gc->y) continue;      

      // OK we have a point in this cell  - saveit

      clats[n] = lat+latlow;
      clons[n++] = lon+lonlow;      

    }

  // ok, now what is the centre point
  clat = clon = 0;
  for(i=0;i<n;i++)
    {
      clat += clats[i];
      clon += clons[i];
    }
  clat /= (double)(n);
  clon /= (double)(n);  

  // which point is closest to this??

  int cl = 0;
  double thisdist, cldist = lldistance(clats[0], clons[0],clat,clon,'M');
  for(i=1;i<n;i++)
    {
      thisdist = lldistance(clats[i], clons[i],clat,clon,'M');
      if(thisdist < cldist) {cldist = thisdist; cl = i;}
    }

  // OK our  lms site is i

  char msg[20];
  sprintf(msg,"LMS_%05d", ++inferred_lms);
  addloc(2, 10, clats[cl],clons[cl],msg);

  addlastmileout(clats[cl],clons[cl]);
  nparcellocs++;

  //now add default settings for the set
  LOCATION *lp = &locs[nlocs-1];

  lp->nvans = 3;
  lp->maxparcels = 100;
  lp->maxmins = 600;
  lp->earlystart = 720;
  lp->loadtime = 30;      
  lp->speedfactor = 0.3;
  lp->permile = 5;
  lp->perhour = 15;


  printf("LMS  established at %f %f\n", lp->lat, lp->lon);
}

int addlastmileout(double lat, double lon)
{
  FILE *f = fopen("__lmsites.txt","a");

  fprintf(f,"%f %f\n", lat,lon);
  fclose(f);

}

  


/****************************************************************************************************\
 
  we first go for the cth entry in u. If u says it is unused (0) then fine.  
  else we keep incrementing by 1 until we find first unused one.

\****************************************************************************************************/
int get_indicated_centre(int c, char *u, int n)
{
  int a = c;
  
  while ( u[a]==1 )
    {
      // increment a
      a++;
      if(a==n) {a = 0;} // wraparound
    }
  return a;
}


int print_areapoly(void)
{
  POINT *p;
  int i;

  for(i=0;i<naplls;i++)
    {
      p = &points[i];
      printf(" --- area poly %d: %f %f\n", i, p->x, p->y);
    }
}


int get_areapoly(void)
{

  FILE *f = fopen(areapoly,"r");

  char instr[100];
  int r;
  double lat, lon;
  
  while(1)
    {
      r = fscanf(f,"%s",instr);

      if(r<1) break;

      lat = atof(instr);
      fscanf(f,"%s",instr);
      lon = atof(instr);
      add_apll(lat,lon);
    }

  fclose(f);
}


double filter_pclls(void)
{

  FILE *f = fopen(pclls,"r");

  char instr[100];
  int r;
  double lat, lon;

  while(1)  
    {
      r = fscanf(f,"%s",instr);

      if(r<1) break;

      lat = atof(instr);
      fscanf(f,"%s",instr);
      lon = atof(instr);
      add_pcll(lat,lon);

    }
  
  fclose(f);
  
}
  


void add_apll(double lat, double lon)
{

  if(naplls==0)
    {
      aplats = (double *)malloc(apllmalloc*sizeof(double));
      aplons = (double *)malloc(apllmalloc*sizeof(double));
      points = (POINT *)malloc(apllmalloc*sizeof(POINT));

      napllmalloced = apllmalloc;
    }
  else if ((naplls + 1) >= napllmalloced)
    {
      napllmalloced += apllmalloc;
      aplats = (double *)realloc(aplats,napllmalloced*sizeof(double));
      aplons = (double *)realloc(aplons,napllmalloced*sizeof(double));
      points = (POINT *)realloc(points,napllmalloced*sizeof(POINT));      
    }
  
  aplats[naplls] = lat;
  aplons[naplls++] = lon;  

  POINT *p = &points[naplls-1];
  p->x = lat;
  p->y = lon;
  npoints = naplls;
}


void add_pcll(double lat, double lon)
{


  if(inpoly(lat,lon)==0) return;
  
  if(npcs==0)
    {
      pclats = (double *)malloc(pcmalloc*sizeof(double));
      pclons = (double *)malloc(pcmalloc*sizeof(double));
      npcmalloced = pcmalloc;
    }
  else if ((npcs + 1) >= npcmalloced)
    {
      npcmalloced += pcmalloc;
      pclats = (double *)realloc(pclats,npcmalloced*sizeof(double));
      pclons = (double *)realloc(pclons,npcmalloced*sizeof(double));
    }
  
  pclats[npcs] = lat;
  pclons[npcs++] = lon;  

  printf("includified:   %f %f\n", pclats[npcs-1],  pclons[npcs-1]);
  
}



// WORKING RIGHT HERE
int inpoly(double lat, double lon)
{

  POINT q, *qp;

  qp = &q;
  qp->x = lat;
  qp->y = lon;
  
  int r = InsidePolygon(points, npoints, qp);

  if(r==INSIDE) return 1;
  return 0;

}


int organise(void)
{

  prelims();

  
  doing = __PARCEL;
  parcel_solve();
  // what we need for van-only delivery plan

  doing = __BUS;  
  bus_solve();  
  // solve the bus delivery problem

  doing = __PARCELPLUS;  
  combined_parcel_solve();

  doing = __BUSPLUS;    
  combined_bus_solve();
 
  doing = __LASTMILE;   
  lastmile_solve();    
}

int prelims(void)
{
  int i;
  for(i=0;i<nbusfleets;i++) allfleet += busfleets[i].n;
}



int print_packages(void)
{

  int i, np=0;
  LOCATION *lp;
  PACKAGE *pp;


  printf("------------------packages\n");
  for(i=0;i<npackages;i++)
    {
      pp = &packages[i];
      if(pp->type==__BUS) printf("%d  BUSP ", i);
      else       if(pp->type==__PARCEL) printf("%d  PARC ", i);
      else       if(pp->type==__LOCKERPARCEL) printf("%d  LOPA ", i);
      else       printf("%d  ???? ", i);

      printf(" %d (%d)  %f %f\n", pp->sd_raw, locs[pp->sd_raw].type, locs[pp->sd_raw].lat,  locs[pp->sd_raw].lon);
      
    }

  printf("-----------------\n");

}



int parcel_solve(void)
{
  int i, np=0, ul;
  LOCATION *lp;
  PACKAGE *pp;

  char fname[100];
  sprintf(fname,"%s_parcels_alone.txt",mapfilebase);
  
  mf = fopen(fname,"w");
  
  // show all packages

  print_packages();

  
  for(i=0;i<npackages;i++)
    {
      pp = &packages[i];
      if(pp->type==__BUS) continue; // only bus

      np++;
      // this is definitely a parcel
    }

  if(uselocs_alloced==0)
    {
      uselocs = (LOCATION *)malloc((np + 1) * sizeof(LOCATION)); // includes depot
      uselocs_alloced = 1;
    }
  else
    {
      uselocs = (LOCATION *)realloc(uselocs,(np + 1) * sizeof(LOCATION)); // includes depot

    }

  nuselocs = np+1; // 1st loc is depot, remainder are raw parel delivery sites

  ndeliveries = nuselocs - 1;
  
  lp = &uselocs[0];
  lp->lat = parcel_lat;
  lp->lon = parcel_lon;  

  ul = 1;
  for(i=0;i<npackages;i++)
    {
      pp = &packages[i];
      if(pp->type==__BUS) continue; // only bus
      if(pp->inplay==0) continue;
      lp = &uselocs[ul++];

      lp->lat = locs[pp->sd_raw].lat;
      lp->lon = locs[pp->sd_raw].lon;

      lp->nparcels = 1;
      lp->arrtime = -1;  // no worries on the time
      lp->nolater = lp->noearlier = -1;  // no constraints in play for this run without bundles
      pp->nolater = pp->noearlier = -1;  // no constraints in play for this run without bundles      
      // lp->isbundle = 0;
      
      lp->v = pp->v;
      lp->w =pp->w;
    }


  // setup

  _parcel_vans = parcel_vans;
  _parcel_maxmins = parcel_maxmins;
  _parcel_earlystart = parcel_earlystart;
  _parcel_loadtime = parcel_loadtime;  
  _parcelvan_maxparcels = parcelvan_maxparcels;
  _speedfactor = 1;

  _permile = parcel_permile;
  _perhour = parcel_perhour;
  
  get_parcel_matrix();

  print_matrix();
  
  parcel_optimize();

  fclose(mf);
}

int get_bus_snaps(PACKAGE *pp)
{

  LOCATION *pu = &locs[pp->pu_raw], *sd = &locs[pp->sd_raw];

  float lat = pu->lat, lon = pu->lon;

  int cloc = closest_busloc(lat, lon);

  pp->pu_snap = cloc;

  lat = sd->lat; lon = sd->lon;
  cloc = closest_busloc(lat, lon);

  pp->sd_snap = cloc;

}


int closest_busloc(float lat, float lon)
{

  int i, cl = -1;
  float dist, bestdist;
  LOCATION *lp;
  
  for(i=0;i<nlocs;i++)
    {
      if(locs[i].type == 3) continue;
      if(locs[i].type == 4) continue;
      if(locs[i].type == 24) continue;            
      if(locs[i].type == 28) continue;

      lp = &locs[i];
      dist = lldistance(lat,lon,lp->lat,lp->lon,'M');
      if(cl<0) {cl = i; bestdist = dist;}
      else if(dist < bestdist)  {cl = i; bestdist = dist;}
    }
  
  return cl;

}

int buslocations[10000], nbuslocations=0;

int addbuslocation(int locindex)
{
  int i;

  for(i=0;i<nbuslocations;i++)
    if(buslocations[i]==locindex) return i;

  buslocations[nbuslocations++] = locindex;
  return(nbuslocations -1);
}

typedef struct busdem {

  int puloc, sdloc; // indices into locs
  float putime;
  
} BUSDEM;


int bus_solve(void)
{
  int i, np=0, ul, pu, sd, nbd=0;
  LOCATION *lp;
  PACKAGE *pp;

  char fname[100];

  sprintf(fname,"%s_bus_alone.txt",mapfilebase);
  
  mf = fopen(fname,"w");
    
  for(i=0;i<npackages;i++)
    {
      pp = &packages[i];
      if (!(pp->type == __BUS)) continue;

      printf("snapping before  %d %d   %d %d  -- ", pp->pu_raw, pp->sd_raw, pp->pu_snap, pp->sd_snap);   
      get_bus_snaps(pp);
      printf("after   %d %d\n ", pp->pu_snap, pp->sd_snap);         

      pu = addbuslocation(pp->pu_snap);
      sd = addbuslocation(pp->sd_snap);      

      pp->pu_bindex = pu;
      pp->sd_bindex = sd;      
      nbd++;
    }

  // sort out busdemand structure
  nbusdemand = nbd*2;
  
  busdemands = (BUSDEMAND *)malloc(nbusdemand * sizeof(BUSDEMAND));

  nbd = 0;
  BUSDEMAND *bdp;

  for(i=0;i<npackages;i++)
    {
      pp = &packages[i];
      if (!(pp->type == __BUS)) continue;

      bdp = &busdemands[nbd++];
      bdp->bindex = pp->pu_bindex;
      bdp->before = nbd;
      bdp->preftime = pp->pu_open; 
      bdp->latesttime = 2880; // default
      bdp->type = 0;

      bdp = &busdemands[nbd++];
      bdp->bindex = pp->sd_bindex;
      bdp->before = -1;
      bdp->preftime = -1; 
      bdp->latesttime = 2880; // default
      bdp->type = 1;
    }
  
  
  // add the bus fleet locations
  for(i=0;i< nbusfleets;i++)
    {
      BUSFLEET *bf = &busfleets[i];
      pu = addbuslocation(bf->lindex);
      bf->bindex = pu;
    }

  // and add the lockers for parcel setdown

  np = 0;
  for(i=0;i<nlocs;i++)
    {
      lp = &locs[i];
      if (!((lp->type == 2) && (lp->partype==10))) continue;
      np++;      
    }
  
  sdlockers = (SDLOCKER *)malloc(np *sizeof(SDLOCKER));
  
  for(i=0;i<nlocs;i++)
    {
      lp = &locs[i];
      if (!((lp->type == 2) && (lp->partype==10))) continue;

      pu = addbuslocation(i);
      addsetdownlocker(i,pu,lp);
      // this is a place where buses can setdown parcels      
    }

  /////////////////////////////
  // and add the lockers for parcel pickup

  
  uselocs = (LOCATION *)realloc(uselocs, nbuslocations * sizeof(LOCATION)); // includes depot

  nuselocs = nbuslocations; // 1st loc is depot, remainder are raw parel delivery sites

  
  for(i=0;i<nbuslocations;i++)
    {

      lp = &uselocs[i];

      lp->lat = locs[buslocations[i]].lat;
      lp->lon = locs[buslocations[i]].lon;

    }

  // setup
  
  _speedfactor = bus_speedfactor;  
  _permile = bus_permile;
  _permile = bus_perhour;  
  
  get_parcel_matrix();

  print_matrix();

  bus_optimize();
  fclose(mf);

  // stash 

}



int combined_bus_solve(void)
{
  int i, np=0, ul, pu, sd, nbd=0;
  LOCATION *lp;
  PACKAGE *pp;
  BUNDLE *bp;

  int newnbusdemand = nbusdemand;

  char fname[100];

  sprintf(fname,"%s_bus_withparcels.txt",mapfilebase);
  
  mf = fopen(fname,"w");

  

  // sort out busdemand structure

  // NOOOO we need to update this 

  //  nbusdemand = nbd*2;
  
  //  busdemands = (BUSDEMAND *)realloc(busdemands, nbusdemand * sizeof(BUSDEMAND));


  BUSDEMAND *bdp;

  int newbd = 0;
  for(i=0;i<npackages;i++)
    {
      pp = &packages[i];
      if (pp->isbundle==0) continue;
      newbd +=2;
    }

  nbusdemand += newbd;
  busdemands = (BUSDEMAND *)realloc(busdemands, nbusdemand * sizeof(BUSDEMAND));
  
  for(i=0;i<npackages;i++)
    {
      pp = &packages[i];
      if (pp->isbundle==0) continue;

      bdp = &busdemands[newnbusdemand++];
      bdp->bindex = addbuslocation(pp->bus_pickup_lindex) ; // this right???  
      bdp->before = newnbusdemand;

      bdp->preftime  = bundle_available; // no earlier
      bdp->latesttime = 2880; // default
      bdp->type = 100;
      bdp->nparcels = pp->nparcels;

      bdp = &busdemands[newnbusdemand++];
      bdp->bindex = addbuslocation(pp->bus_setdown_lindex) ; // this right???  
      bdp->before = -1;

      bdp->type = 200;
      bdp->preftime  = bundle_available; // no earlier
      bdp->latesttime = bundle_drop; // no later
      bdp->nparcels = pp->nparcels;
    }

  
  
  uselocs = (LOCATION *)realloc(uselocs, nbuslocations * sizeof(LOCATION)); // includes depot

  nuselocs = nbuslocations; // 1st loc is depot, remainder are raw parel delivery sites

  
  for(i=0;i<nbuslocations;i++)
    {

      lp = &uselocs[i];

      lp->lat = locs[buslocations[i]].lat;
      lp->lon = locs[buslocations[i]].lon;

    }

  _speedfactor = bus_speedfactor;  
  _permile = bus_permile;
  _permile = bus_perhour;  

  
  get_parcel_matrix();

  print_matrix();

  bus_optimize();

  fclose(mf);
}


int assign_package(PACKAGE *pp, float lastmile_radius)
{
  // get closest sdloc
  int sdi, csdi;
  float dist, bestdist, plat = locs[pp->sd_raw].lat, plon = locs[pp->sd_raw].lon;
  SDLOCKER *sdp; 

  for(sdi = 0; sdi < nsdlockers; sdi++)
    {
      sdp = &sdlockers[sdi];
      
      dist = lldistance(plat, plon, uselocs[sdp->bindex].lat, uselocs[sdp->bindex].lon,'M');
      if(sdi==0) {csdi = sdi; bestdist = dist;}
      else if( dist < bestdist) {csdi = sdi; bestdist = dist;}
    }
  if(bestdist > lastmile_radius) return 0;

  // if we are here, we can assign this package to csdi

  pp->sd_assigned = csdi;

  sdp = &sdlockers[csdi];
  sdp->nassigned++;
}

/****************************************************************************************************\
   put the bus-assigned parcels into bundles that are at most bus_maxparcels
\****************************************************************************************************/
int bundle_packages(void)
{
  int i,j,  np=0, ul, pu, sd, nbd=0;
  LOCATION *lp;
  PACKAGE *pp;
  SDLOCKER *sdp;
  
  for(i=0;i<npackages;i++)
    {
      pp = &packages[i];
      if (pp->type == __BUS) continue;

      if(pp->sd_assigned >=0)
	bundle_package(pp);
    }

  // now show bundles

  for(i=0;i<nbundles;i++)
    {
      BUNDLE *bp = &bundles[i];
      printf("b %d  sd %d  n %d : ", i, bp->sdlocker, bp->ncontents);
      for(j=0;j<bp->ncontents;j++) printf("%d ", bp->contents[j]); printf("\n");

    }

}

int bundle_package(PACKAGE *pp)
{

  int i;
  BUNDLE *bp;

  for(i=0;i<nbundles;i++)
    {
      if((bundles[i].sdlocker == pp->sd_assigned) && (bundles[i].ncontents < bus_maxparcels))
	{bp = &bundles[i]; bp->contents[bp->ncontents++] = pp->id; return 1;}
    }

  // if we're still here, we need to create a new bundle

  bp = &bundles[nbundles++];
  bp->sdlocker = pp->sd_assigned;
  bp->ncontents = 1;
  bp->contents[0] = pp->id;

  // what are these?
  int bus_pu_lindex, bus_sd_lindex, van_sd_lindex;   

  // bus picks it up from -- and van sets it down at -- the fleet
  bp->bus_pu_lindex = bp->van_sd_lindex =  busfleets[sdlockers[bp->sdlocker].fleet].lindex;

  // bus sd_lindex is the bus drop point
  bp->bus_sd_lindex = sdlockers[bp->sdlocker].lindex;
  
}


/****************************************************************************************************\
  reconstitute the packages array to contain the bus-assigned packages only as bundles

\****************************************************************************************************/
int reshape_packages(void)
{
  
  PACKAGE *pp;
  int i;

  for(i=0;i<npackages;i++)
    {
      pp = &packages[i];
      if(pp->sd_assigned >=0) pp->inplay = 0;
    }

  for(i=0;i<nbundles;i++)
    {
      pp = &packages[npackages++];
      pp->id = npackages - 1;

      pp->inplay = 1;

      pp->type = __PARCEL; // parcel demand
      pp->isbundle = 1;
      pp->bundle_index = i;
      
      pp->nparcels = bundles[i].ncontents;
      
      pp->arrtime = pp->nolater = busdepot_arrival;
      pp->noearlier = -1;
      
      pp->sd_raw = pp->sd_snap = bundles[i].van_sd_lindex; // the one we just added -- already a  snap //  .......................

      pp->bus_pickup_lindex = pp->sd_raw;
      pp->bus_setdown_lindex = bundles[i].bus_sd_lindex;
    }


}

/****************************************************************************************************\
 
 solve the parcel van delivery problem with the bus-assigned parcels now bundled up for particular
 drop sites, and due to   arrive at their drop sites by lastmile_arrival time 

\****************************************************************************************************/
int combined_parcel_solve(void)
{
  int i, np=0, ul, pu, sd, nbd=0;
  LOCATION *lp;
  PACKAGE *pp;


  char fname[100];
  sprintf(fname,"%s_parcels_withbus.txt",mapfilebase);  
  mf = fopen(fname,"w");

  
  //  parcel drops are the sdlockers
  //  we want to assign parcels to these


  for(i=0;i<npackages;i++)
    {
      pp = &packages[i];
      if (pp->type == __BUS) continue;

      pp->sd_assigned = -1;  // by default, not in a locker

      assign_package(pp,lastmile_radius);
    }

  for(i=0;i<nsdlockers;i++)
    printf("assigned to locker %d:  %d\n", i+1, sdlockers[i].nassigned);

  bundle_packages();

  reshape_packages();

  
  for(i=0;i<npackages;i++)
    {
      pp = &packages[i];
      if(pp->type==__BUS) continue; // only bus
      if(pp->inplay == 0) continue; // excluded the individuals now assigned to bundles
      np++;
    }

  nuselocs = np  + 1 + nbusfleets;
  uselocs = (LOCATION *)realloc(uselocs, nuselocs * sizeof(LOCATION)); // includes depot
  ndeliveries = np;
  
  // 1st loc is depot, last locs are bus fleet depots, remainder are raw parcel delivery sites

  lp = &uselocs[0];
  lp->lat = parcel_lat;
  lp->lon = parcel_lon;  

  ul = 1;
  for(i=0;i<npackages;i++)
    {
      pp = &packages[i];
      if(pp->type==__BUS) continue; // only bus
      if(pp->inplay==0) continue;

      lp = &uselocs[ul++];
      // WORKING RIGHT HERE

      lp->lat = locs[pp->sd_raw].lat;
      lp->lon = locs[pp->sd_raw].lon;

      lp->nparcels = pp->nparcels;

      if(pp->nparcels > 1)
	{ // change lat and lon to those of the drop fleet
	  ; // no change needed we already set the latlon correctly
	 }
      lp->arrtime = pp->arrtime;
      lp->noearlier = pp->noearlier;
      lp->nolater = pp->nolater;
      
      lp->v = pp->v;
      lp->w =pp->w;
    }

  BUSFLEET *b;
  
  for(i=0;i<nbusfleets;i++)
    {
      lp = &uselocs[ul++];
      b = &busfleets[i];
      b->lindex = ul - 1;
      lp->lat = locs[b->lindex].lat;
      lp->lon = locs[b->lindex].lon;
    }
  
 // setup

  _parcel_vans = parcel_vans;
  _parcel_maxmins = parcel_maxmins;
  _parcel_earlystart = parcel_earlystart;
  _parcel_loadtime = parcel_loadtime;
  _parcelvan_maxparcels = parcelvan_maxparcels;
  _speedfactor = 1;

  _permile = parcel_permile;
  _perhour = parcel_perhour;
  
  
  get_parcel_matrix();

  print_matrix();
  
  parcel_optimize();


  fclose(mf);
}


/****************************************************************************************************\
 
 the bundle drop sites now have all of their bundles delivered.
 solve the mini delivery problem at each one

\****************************************************************************************************/
int lastmile_solve(void)
{
  int i, np=0, ul, pu, sd, nbd=0;
  LOCATION *lp;
  PACKAGE *pp;

  char fname[100];
  sprintf(fname,"%s_lastmile.txt",mapfilebase);  
  mf = fopen(fname,"w");

  
  //  parcel drops are the sdlockers
  //  we want to assign parcels to these

  for(i=0;i<nsdlockers;i++)
    lastmile_solve_site(&sdlockers[i]);

}

/****************************************************************************************************\
 
 solve the mini delivery problem at this setdown locker site

// to do this, we build up a list of locations that we need to deliver to

\****************************************************************************************************/

int lastmile_solve_site(SDLOCKER *sdp)
{
  LOCATION *lp;
  BUNDLE *bup;
  PACKAGE *pp;
  int i, p, ppid, ul;


  // first, how many locs do we need for the bundles?
  ul = 0;
  
  for(i=0;i<nbundles;i++)
    {
      bup = &bundles[i];
      if(bup->sdlocker != sdp->id) continue;

      ul += bup->ncontents;
    }

  if(ul<1) return 0;
  
  nuselocs = ul+1; // including the site itself
  ul = 0;

  uselocs = (LOCATION *)realloc(uselocs, nuselocs * sizeof(LOCATION)); // includes depot
  ndeliveries = nuselocs - 1;  
  
  
  // first location is the 'depot' in this case, which is the setdown locher site.

  lp = &uselocs[0];
  lp->lat = locs[sdp->lindex].lat;
  lp->lon = locs[sdp->lindex].lon;
  ul = 1;
  // now, pickup the bundles at this site and put together their locations  


  for(i=0;i<nbundles;i++)
    {
      bup = &bundles[i];
      if(bup->sdlocker != sdp->id) continue;

      for(p=0;p< bup->ncontents; p++)
	{
	  ppid = bup->contents[p];
	  pp = &packages[ppid];

	  lp = &uselocs[ul++];
	  // so, where do we need to set this fellow down?
	  lp->lat = locs[pp->sd_raw].lat;
	  lp->lon = locs[pp->sd_raw].lon;     

	  //timing

	  pp->arrtime = bundle_available_lastmile;
	  lp->noearlier = pp->noearlier = bundle_available_lastmile;
	  lp->nolater = pp->nolater = -1;
	}
    }


 // setup

  _parcel_vans = sdp->nvans;
  _parcel_maxmins = sdp->maxmins;
  _parcel_earlystart = sdp->earlystart;
  _parcel_loadtime = sdp->loadtime;
  _parcelvan_maxparcels = sdp->maxparcels;
  _speedfactor = sdp->speedfactor;


  
  _permile = sdp->permile;
  _perhour = sdp->perhour;


  
  printf("\n-----------------------\n solving lastmile for site %d with %d locs speedfactor %f earlystart %f loadtime %f\n-------------------------\n",
	 sdp->id,  nuselocs, _speedfactor, _parcel_earlystart, _parcel_loadtime);
  
  // MAKE COMBINED VERSYIONS BELOW

 
  get_parcel_matrix();

  print_matrix();
  
  parcel_optimize();

}










int printbusdemand(void)
{

  int i;
  BUSDEMAND *bdp;

  printf("busdemand:  %d  ---------------------\n",nbusdemand);
  
  for(i=0;i<nbusdemand;i++)
    {
      bdp = &busdemands[i];
      printf("%d  %d %d %d  %f\n", i, bdp->bindex, bdp->before, bdp->type, bdp->preftime);
    }
}

int print_matrix(void)
{

  int i, j;

  printf("skipping matrix (%d locations) ...\n", nuselocs);
  return 1;
  
  for(i=0;i<nuselocs;i++)
    for(j=0;j<nuselocs;j++)
      {
	printf(" %d %d   %f   %f  \n", i, j, tmat[i][j],   dmat[i][j]);

      }


}


/****************************************************************************************************\
   solve with a simple but sensible genetic algorithm
\****************************************************************************************************/
int parcel_optimize(void)
{

  allocate_parcel();
  initialize_parcel();

  run_parcel();

}


/****************************************************************************************************\
   solve with a simple but sensible genetic algorithm
\****************************************************************************************************/
int bus_optimize(void)
{

  //  printf("   ===   COST    allfleet %d  nbusdemand %d  nbuslocations %d\n", allfleet, nbusdemand, nbuslocations);
  allocate_bus();

  initialize_bus();
  
  run_bus();

}

int allocate_bus(void)
{
  // we reuse the population already alloced for the parcel problem

  int i, j, bi, thisb;
  BUSPLAN *bp;
  CANDIDATE *cp;
  
  
  for(i=0;i< (popsize+3); i++)
    {
      cp = &population[i];

      cp->nvanplans = 0;  // mainly this stops the copying
      
      cp->nbusplans = allfleet;
      
      cp->busplans = (BUSPLAN *)malloc(allfleet * sizeof(BUSPLAN));

      for(j=0;j< allfleet; j++)
	{
	  bp = &cp->busplans[j];
	  bp->ndels = 0;
	  bp->sequence = (int *)malloc(nbusdemand * sizeof(int));
	  zapobjectives(bp->objectives);
	}

      thisb=0;
      for(j=0;j<nbusfleets;j++)
	for(bi=0;bi<busfleets[j].n; bi++)
	  {
	    bp = &cp->busplans[thisb++];
	    bp->fleet = j;
	  }
	  zapobjectives(cp->objectives);
    }

  printbusdemand();
}



int allocate_parcel(void)
{

  population = (CANDIDATE *)malloc((popsize+3) * sizeof(CANDIDATE));
  CANDIDATE *cp;
  VANPLAN *vp;
  
  int i, j;


  for(i=0;i< (popsize+3); i++)
    {
      cp = &population[i];
      cp->nvanplans = _parcel_vans;
      
      cp->vanplans = (VANPLAN *)malloc(_parcel_vans * sizeof(VANPLAN));

      for(j=0;j< _parcel_vans; j++)
	{
	  vp = &cp->vanplans[j];
	  vp->ndels = 0;
	  vp->sequence = (int *)malloc(nuselocs * sizeof(int));
	  zapobjectives(vp->objectives);
	}
	  zapobjectives(cp->objectives);
    }

}


/****************************************************************************************************\
   there are nuselocs-1 parcels and parcel_van vans
\****************************************************************************************************/
int initialize_parcel(void)
{
  int i, j, done = 0;

  int *upp = (int *)malloc(nuselocs * sizeof(int)), nupp = nuselocs - 1;;
  CANDIDATE *cp;
  
  // build unplanned parcel list

  for(i=0;i<popsize;i++)
    {
      cp = &population[i];

      for(j=0;j<ndeliveries;j++) upp[j] = j+1; // loc of parcel
      init_build_parcels(cp,nupp,upp);
      eval_parcel(cp,0);
      print_parcel(cp);
    }

  free(upp);

}


/****************************************************************************************************\
   there are nuselocs-1 parcels and parcel_van vans
\****************************************************************************************************/
int initialize_bus(void)
{
  int i, j, done = 0;

  int *upp = (int *)malloc(nbusdemand * sizeof(int)), nupp = nbusdemand;;
  CANDIDATE *cp;
  
  // build unplanned parcel list

  
  for(i=0;i<popsize;i++)
    {
      cp = &population[i];

      for(j=0;j<nbusdemand;j++){ upp[j] = j; }// loc of parcel
      init_build_bus(cp,nupp,upp);
      //      eval_parcel(cp);
      print_bus(cp);
    }

  free(upp);

}



/****************************************************************************************************\
   there are nuselocs-1 parcels and parcel_van vans
\****************************************************************************************************/
int run_parcel(void)
{
  int i, j, done = 0, iter, s;

  int *upp = (int *)malloc(nuselocs * sizeof(int)), nupp = nuselocs - 1;;
  CANDIDATE *cp;
  
  // build unplanned parcel list

  for(iter = 0; iter < parcel_iterations; iter++)
    {
      s = select_candidate(2);  // standard binary tournament

      mutate_parcel(s); // make a copy in pop[popsize] and mutate it

      eval_parcel(&population[popsize],0);
      
      replace_candidate(2,__PARCEL);  // if the mutation is good enough, put it in pop

      if(iter%100000==0)
	print_progress(iter,2);

    }

  print_best_candidate(2,__PARCEL);

}

/****************************************************************************************************\
   there are ...
\****************************************************************************************************/
int run_bus(void)
{
  int i, j, done = 0, iter, s;

  int *upp = (int *)malloc(nuselocs * sizeof(int)), nupp = nuselocs - 1;;
  CANDIDATE *cp;
  
  // build unplanned parcel list

  for(iter = 0; iter < bus_iterations; iter++)
    {
      s = select_candidate(3);  // standard binary tournament

      mutate_bus(s); // make a copy in pop[popsize] and mutate it

      eval_bus(&population[popsize],0);
      
      replace_candidate(3,__BUS);  // if the mutation is good enough, put it in pop

      if(iter%100000==0)
	print_progress(iter,3);

    }

  print_best_candidate(3,__BUS);

}

int replace_candidate(int selobj, int what)
{

  int i, w=0;

  for(i=1;i<popsize;i++)
    if(population[i].objectives[selobj] > population[w].objectives[selobj]) {w = i;};

  if(population[popsize].objectives[selobj] < population[w].objectives[selobj])
    { // overwrite

      if(what==__PARCEL)
	copy_parcel(&population[popsize],&population[w]);
      else  if(what==__BUS)
	copy_bus(&population[popsize],&population[w]);	
    }
}


int print_progress(int iter, int selobj)
{

  int i, w=0;

  for(i=1;i<popsize;i++)
    if(population[i].objectives[selobj] < population[w].objectives[selobj]) {w = i;};
  CANDIDATE *cp = &population[w];

  double *o = cp->objectives;
  printf("  %d  best so far  %f %f %f %f\n", iter, o[0], o[1],o[2], o[3]);


}


int print_best_candidate(int selobj, int what)
{

  int i, w=0;

  for(i=1;i<popsize;i++)
    if(population[i].objectives[selobj] < population[w].objectives[selobj]) {w = i;};

  CANDIDATE *cp = &population[w];
  if(what==__PARCEL)  printmore_parcel(cp);
  else   if(what==__BUS)  printmore_bus(cp);

}


int select_candidate(int selobj)
{


  int a, b;

  a = rand()%popsize;
  b = rand()%popsize;

  if(population[a].objectives[selobj] < population[b].objectives[selobj]) return a;
  else return b;

  // down to dist and time  
}


int mutate_parcel(int candi)
{

  CANDIDATE *cp = &population[candi], *m = &population[popsize];

  copy_parcel(cp,m);

  operate_parcel(m);
  
}


int mutate_bus(int candi)
{

  CANDIDATE *cp = &population[candi], *m = &population[popsize];

  copy_bus(cp,m);

  operate_bus(m);
  
}



int operate_parcel(CANDIDATE *m)
{

  int choice = rand()%5;
  
  if(choice==0)  invert_parcel(m,1); // invert
  else if (choice == 1) invert_parcel(m,2); // double invert
  else if (choice == 2) invert_parcel(m,3); // double invert  
  else if (choice == 3) move_parcel(m);
  else if (choice == 4) move_parcel(m);  


}



int operate_bus(CANDIDATE *m)
{

  int choice = rand()%4;
  // choice = 3;

  copy_bus(m,&population[popsize+1]);
	   
  if(choice==0)  invert_bus(m,1); // invert
  else if (choice == 1) invert_bus(m,2); // double invert
  else if (choice == 2) invert_bus(m,3); // double invert  
  else if (choice == 3) move_bus(m);
  else if (choice == 4) move_bus(m);  

  check_integrity(m);
}

int check_integrity(CANDIDATE *m)
{

  int i,j, x, ndt = 0, stop = 0, prob;
  int vans[1000], pos[1000], count[1000];
  int diffnds[1000], ndiffnds=0;
  BUSPLAN *bp;
  CANDIDATE *base = &population[popsize+1];
  
  for(i=0;i<nbusdemand;i++)
    {vans[i] = pos[i] = -1; count[i] = 0;}
  
  for(i=0;i<allfleet;i++)
    {
      bp = &m->busplans[i];

      ndt += bp->ndels;
      if(bp->ndels != base->busplans[i].ndels) diffnds[ndiffnds++] = i;
      
      for(j=0;j<bp->ndels;j++)
	{
	  x = bp->sequence[j];
	  vans[x] = i;
	  pos[x] = j;
	  count[x]++;
	}

    }
  if(ndt != nbusdemand)
    {
      stop = 1;
    }

  for(i=0;i<nbusdemand;i++)
    {
      if(count[x] != 1) {prob = x; stop = 1;}
	
    }


  if(stop)
    {
      print_bus(base);
      printf("demand %d has count %d\n", x, count[x]);
      print_bus(m);

      printf("buses with changed ndel numbers: ");
      for(i=0;i<ndiffnds;i++) printf("%d ", diffnds[i]); printf("\n");
      
      exit(1);
    }
  
}




int invert_parcel(CANDIDATE *m, int n)
{

  int i, cands[100], ncands = 0;

  // which have dels?
  for(i=0;i<_parcel_vans; i++) if(m->vanplans[i].ndels>1) cands[ncands++] = i;

  if(ncands==0) return 1;
    
  int v = cands[rand()%ncands];

  VANPLAN *vp = &m->vanplans[v];

  invert_van(vp,--n);

}


int invert_bus(CANDIDATE *m, int n)
{

  int i, cands[100], ncands = 0;

  // which have dels?
  for(i=0;i<allfleet; i++) if(m->busplans[i].ndels>2) cands[ncands++] = i;

  if(ncands==0) return 1;
    
  int v = cands[rand()%ncands];

  BUSPLAN *vp = &m->busplans[v];

  invert_busplan(vp,--n);

}



int move_parcel(CANDIDATE *m)
{

  int i, cands[100], ncands = 0;

  // which have dels?
  for(i=0;i<_parcel_vans; i++) if(m->vanplans[i].ndels>0) cands[ncands++] = i;

  if(ncands<2) return 1;
    
  int v1 = cands[rand()%ncands];
  int v2 = cands[rand()%ncands];  

  while(v2==v1)  v2 = cands[rand()%ncands];  

  VANPLAN *vp1 = &m->vanplans[v1];
  VANPLAN *vp2 = &m->vanplans[v2];  

  int p = extract_parcel(vp1);
  insert_parcel_rand(vp2,p);

}



int move_bus(CANDIDATE *m)
{

  int i, cands[100], ncands = 0;

  // which have dels?
  for(i=0;i< allfleet; i++) if(m->busplans[i].ndels>0) cands[ncands++] = i;

  if(ncands<2) return 1;
    
  int v1 = cands[rand()%ncands];
  int v2 = cands[rand()%ncands];  

  while(v2==v1)  v2 = cands[rand()%ncands];  

  BUSPLAN *vp1 = &m->busplans[v1];
  BUSPLAN *vp2 = &m->busplans[v2];  

  //  printf(" vp1 before: ");
  //  for(i=0;i<vp1->ndels;i++) printf("%d ", vp1->sequence[i]); printf("\n");
  int p = extract_busdem(vp1);
  //  printf(" vp1 after extract: ");
  //  for(i=0;i<vp1->ndels;i++) printf("%d ", vp1->sequence[i]); printf("\n");
  fixbusplan(vp1);

  //  printf(" vp2 before insert: ");
  //  for(i=0;i<vp2->ndels;i++) printf("%d ", vp2->sequence[i]); printf("\n");
  insert_busdem_rand(vp2,p);
  //  printf(" vp2 after insert: ");
  //  for(i=0;i<vp2->ndels;i++) printf("%d ", vp2->sequence[i]); printf("\n");
  


  fixbusplan(vp2);  
  
}




int extract_parcel(VANPLAN *vp)
{

  int n = vp->ndels;
  if(n<1) return -1;

  int pos = rand()%n;

  int p = vp->sequence[pos];

  int i;
  for(i=pos;i<(n-1);i++)
    vp->sequence[i] = vp->sequence[i+1];
  vp->ndels--;

  return p;

}


int extract_busdem(BUSPLAN *vp)
{

  int n = vp->ndels;
  if(n<1) return -1;

  int pos, p, s, x;

  do
    {
      pos = rand()%n;
      p = vp->sequence[pos];
    }
  while (busdemands[p].before < 0);

  s = busdemands[p].before;

  int newp[1000], nnew=0, i;
  
  for(i=0;i<n;i++)
    {
      x = vp->sequence[i];
      if( (x != p) && (x != s)) newp[nnew++] = x;
    }

  for(i=0;i<nnew;i++) vp->sequence[i] = newp[i];
  vp->ndels -=2;

  return p;

}


int invert_van(VANPLAN *vp, int n)
{

  int i, chunk, a, b, tmp, nd = vp->ndels, news[500], ibit[4];

  if(nd<2) return 1;

  if(nd==2)
    {
      tmp = vp->sequence[0];
      vp->sequence[0]= vp->sequence[1]; 
      vp->sequence[1]= tmp;
      return 1;
    }
  
  a = rand()%nd-1;

  chunk = 1 + rand()%2;

  a = rand()%(nd-chunk);
  
  b = (a + chunk);
  
  if(b>= nd) return 1;


  for(i=0;i<a;i++)  news[i] = vp->sequence[i];
  tmp = chunk;
  for(i=a;i<=b;i++)  ibit[tmp--] = vp->sequence[i];
  tmp = 0;
  for(i=a;i<=b;i++)  news[i] = ibit[tmp++];
  for(i=b+1;i<nd;i++) news[i] = vp->sequence[i];

  for(i=0;i<nd;i++) vp->sequence[i] = news[i];

  if(n>1)
    invert_van(vp,1);
}


int invert_busplan(BUSPLAN *vp, int n)
{

  int i, chunk, a, b, tmp, tmp2, nd = vp->ndels, news[500], ibit[4];

  if(nd<4) return 1;

  if(nd==4)
    {
      tmp = vp->sequence[0];
      tmp2 = vp->sequence[1];      
      vp->sequence[0]= vp->sequence[2];
      vp->sequence[1]= vp->sequence[3];       
      vp->sequence[2]= tmp;
      vp->sequence[3]= tmp2;      
      fixbusplan(vp);
      return 1;
    }
  
  a = rand()%nd-2;

  chunk = 1 + rand()%2;

  a = rand()%(nd-chunk);
  
  b = (a + chunk);
  
  if(b>= nd) return 1;


  for(i=0;i<a;i++)  news[i] = vp->sequence[i];
  tmp = chunk;
  for(i=a;i<=b;i++)  ibit[tmp--] = vp->sequence[i];
  tmp = 0;
  for(i=a;i<=b;i++)  news[i] = ibit[tmp++];
  for(i=b+1;i<nd;i++) news[i] = vp->sequence[i];

  for(i=0;i<nd;i++) vp->sequence[i] = news[i];

  if(n>1)
    invert_busplan(vp,1);
  else fixbusplan(vp);
    
}

int fixbusplan(BUSPLAN *bp)
{
  int i, nd = bp->ndels;
  int *seq = bp->sequence;
  int pos[5000], sds[5000];
  int bd, thisdem, tmp;

  //  printf("before: ");
  //  for(i=0;i<bp->ndels;i++) printf("%d ", seq[i]); printf("\n");
  
  for(i=0;i<nbusdemand;i++) pos[i]= -1;
  
  for(i=0;i<nd;i++)
    {
      bd = seq[i];
      pos[bd] = i;
    }

  for(i=0;i<nbusdemand;i+=2)
    {
      if(pos[i]<0) continue;  // not contained here
      
      if(pos[i] > pos[i+1])  // they need swapping

	{
	  tmp = seq[pos[i]];
	  seq[pos[i]] = seq[pos[i+1]];
	  seq[pos[i+1]] = tmp;
	}
    }

  //  printf("after: ");
  //  for(i=0;i<bp->ndels;i++) printf("%d ", seq[i]); printf("\n");
}





int copy_vanplan(VANPLAN *f, VANPLAN *t)
{
  t->ndels = f->ndels;
  int i;
  for(i=0;i<f->ndels;i++)  t->sequence[i] = f->sequence[i];
  for(i=0;i<5;i++) t->objectives[i] = f->objectives[i];
}

int copy_busplan(BUSPLAN *f, BUSPLAN *t)
{
  t->ndels = f->ndels;
  t->fleet = f->fleet;
  
  int i;
  for(i=0;i<f->ndels;i++)  t->sequence[i] = f->sequence[i];
  for(i=0;i<5;i++) t->objectives[i] = f->objectives[i];
}


int copy_parcel(CANDIDATE *f, CANDIDATE *t)
{

  t->nvanplans = f->nvanplans;
  t->maxon = f->maxon;
  int i;

  for(i=0;i<f->nvanplans; i++)
    {
      VANPLAN *vf = &f->vanplans[i], *vt = &t->vanplans[i];
      copy_vanplan(vf,vt);
    }

  for(i=0;i<5;i++) t->objectives[i] = f->objectives[i];
}


int copy_bus(CANDIDATE *f, CANDIDATE *t)
{

  t->nbusplans = f->nbusplans;
  t->maxon = f->maxon;
  
  int i;

  for(i=0;i<f->nbusplans; i++)
    {
      BUSPLAN *vf = &f->busplans[i], *vt = &t->busplans[i];
      copy_busplan(vf,vt);
    }

  for(i=0;i<5;i++) t->objectives[i] = f->objectives[i];
}



int print_parcel(CANDIDATE *cp)
{
  VANPLAN *vp;

  int i, j;
  
  printf("\n\n--------------------------------------------\n");
  for(i=0;i<_parcel_vans;i++)
    {
      vp = &cp->vanplans[i];
      printf("van %02d: ", i+1);
      for(j=0;j<vp->ndels;j++)
	printf("%d ", vp->sequence[j]);
      printf(" |  %f   %f   %f\n", vp->objectives[0],vp->objectives[1],vp->objectives[2]);
    }
  printf("%f    %f     %f--------------------------------------------\n",cp->objectives[0],cp->objectives[1],cp->objectives[2]);
}


int printmore_parcel(CANDIDATE *cp)
{
  VANPLAN *vp;
  int i, j;
  int nvans = 0, npars = 0;
  
  eval_parcel(cp,1);

  printf("\n\n---------------------------------------------------\n");
  for(i=0;i<_parcel_vans;i++)
    {
      vp = &cp->vanplans[i];
      if(vp->ndels>0) {nvans++; npars += vp->ndels;}
      mapoutvan(vp,mf);
      printf("van %02d based at %f %f - loadtime %f  ", i+1, uselocs[0].lat,uselocs[0].lon, _parcel_loadtime);
      for(j=0;j<vp->ndels;j++)
	printf("%d ", vp->sequence[j]);
      printf(" |  %f   %f   %f\n", vp->objectives[0],vp->objectives[1],vp->objectives[2]);
    }
  printf("%f    %f     %f  CST:  %f--------------------------------------------\n",cp->objectives[0],cp->objectives[1],cp->objectives[2], cp->objectives[4]);
  printf(" %d COST:  %f   vans %d  parcels %d --------------------------------------\n",doing,cp->objectives[4],nvans,npars);  
}


int print_bus(CANDIDATE *cp)
{
  BUSPLAN *vp;

  eval_bus(cp,1);
  int i, j, nvans = 0, npars = 0;
  
  printf("\n\n--------------------------------------------\n");
  for(i=0;i<allfleet;i++)
    {
      vp = &cp->busplans[i];
      printf("bus %02d fleet %d: ", i+1, vp->fleet+1);

      if(vp->ndels>0) {nvans++; npars += vp->ndels;}
      
      for(j=0;j<vp->ndels;j++)
	printf("%d ", vp->sequence[j]);
      printf(" |  %f   %f   %f  %f\n", vp->objectives[0],vp->objectives[1],vp->objectives[2],vp->objectives[3]);
    }
  printf("%f    %f     %f   %f--------------------------------------------\n",cp->objectives[0],cp->objectives[1],cp->objectives[2],cp->objectives[3]);
  // printf(" %d COST:  %f--------------------------------------------\n",doing,cp->objectives[4]);  
}

int printmore_bus(CANDIDATE *cp)
{
  BUSPLAN *vp;

  int i, j, nvans = 0, npars = 0;  
  printf("\n\n--------------------------------------------\n");
  for(i=0;i<allfleet;i++)
    {
      vp = &cp->busplans[i];
      printf("bus %02d fleet %d: ", i+1, vp->fleet+1);


      if(vp->ndels>0) {nvans++; npars += vp->ndels;}
      
      for(j=0;j<vp->ndels;j++)
	printf("%d ", vp->sequence[j]);
      printf(" |  %f   %f   %f  %f\n", vp->objectives[0],vp->objectives[1],vp->objectives[2],vp->objectives[3]);
    }
  printf("%f    %f     %f   %f--------------------------------------------\n",cp->objectives[0],cp->objectives[1],cp->objectives[2],cp->objectives[3]);
  printf(" %d COST:  %f    vans  %d  parcs %d maxon %d-----------------------------------\n",doing,cp->objectives[4], nvans, npars,cp->maxon);  
}



int mapoutvan(VANPLAN *vp, FILE *f)
{


  
  if(vp->ndels < 1) return 1;

  int i, lastloc = 0, thisloc;

  char col[20];
  if(doing == __PARCEL) strcpy(col,"blue");
  else strcpy(col,"green");
  
  for(i=0;i<vp->ndels;i++)
    {
      thisloc = vp->sequence[i];
      fprintf(f,"line %f %f  %f %f 3 %s\n", uselocs[lastloc].lat,uselocs[lastloc].lon,uselocs[thisloc].lat,uselocs[thisloc].lon, col);
      lastloc = thisloc;
    }
  thisloc = 0;
  fprintf(f,"line %f %f  %f %f 3 %s\n", uselocs[lastloc].lat,uselocs[lastloc].lon,uselocs[thisloc].lat,uselocs[thisloc].lon,col);

}

int eval_van(VANPLAN *vp, int o) //
{

  zapobjectives(vp->objectives);

  float tt = 0, km = 0; 
  int i, lastloc = 0, packages_onboard = 0, sdp; // 0 is the parcel depot location
  int thisone;
  
  LOCATION *thisl;
  
  float tstart = _parcel_earlystart, thist, tnow, late=0, vwait = 0;

  tnow = tstart + _parcel_loadtime;
  
  for(i=0;i<vp->ndels; i++)
    {
     thisone  = vp->sequence[i];
     thisl = &uselocs[thisone];
     
     packages_onboard += thisl->nparcels;

     if(packages_onboard > _parcelvan_maxparcels)
       { // reload 
	 thist = travtime(thisone,0); // go back to base
	 km += travdist(thisone,0);
	 tt += thist;
         tnow += thist;
	 tnow += _parcel_loadtime;  // reload
	 packages_onboard = thisl->nparcels;
	 lastloc = 0;
       }
     // travel there     

     thist = travtime(lastloc,thisone);
     km += travdist(lastloc,thisone);
     
     tt += thist;
     tnow += thist;

     
     if(thisl->noearlier >=0 )
       {
	 if(i==0) // first trip
	   {
	     if(tnow < thisl->noearlier)
	       { tnow = thisl->noearlier; tstart = tnow - (thist + _parcel_loadtime);}
	   }
	 else
	   {
	     if(tnow < thisl->noearlier) {vwait +=  thisl->noearlier - tnow; tnow = thisl->noearlier;}
	   }
       }
     if(thisl->nolater >= 0)
       {
	 if(tnow > thisl->nolater) {late +=   (tnow - thisl->nolater);}
       }
     
     // visit
     if(thisl->nparcels > 1) tnow += bundle_visit_time;
     else tnow += parcel_visit_time;

     if(o)
       {
	 if(i==0)
	   {printf("TSTART is %f \n",tstart);
	     printf("AFTER LOADING is %f \n",tstart + _parcel_loadtime);
	   }
	 printf("-->TNOW is %f  *noe %f  nol %f  acclate %f \n",tnow, thisl->noearlier, thisl->nolater,late  );

       }
    }

    
  // back to base
  
    thist = travtime(vp->sequence[vp->ndels-1],0);
    km += travdist(vp->sequence[vp->ndels-1],0);
    
    tt += thist;
    tnow += thist;
    if(o)   printf("FINAL TNOW is %f \n",tnow);

  // OK what are the pens?

    float overwas, over = (tnow - tstart) - _parcel_maxmins;
    overwas = over;
  over += late;

  
  if(over < 0) over = 0;

  if(o)   printf("FINAL TNOW is %f   over %f late %f (shift %f vs max %f) \n",tnow,overwas,late,tnow-tstart,_parcel_maxmins);
  
  vp->objectives[0] = tnow - tstart;
  vp->objectives[1] = over; // must reduce this
  vp->objectives[2] = tt + over*1000; // must reduce this  


  // cost
  float cost =   _permile * km;
  cost += (_perhour * (tnow - tstart)/60);
  vp->objectives[4] = cost;
}




int eval_parcel(CANDIDATE *cp, int o)
{

  int i;

  
  for(i=0;i<_parcel_vans;i++)
    eval_van(&cp->vanplans[i],o);

  zapobjectives(cp->objectives);  

  for(i=0;i<_parcel_vans;i++)
    {
      cp->objectives[0] += cp->vanplans[i].objectives[0];
      cp->objectives[1] += cp->vanplans[i].objectives[1];
      cp->objectives[4] += cp->vanplans[i].objectives[4];        
    }

  cp->objectives[2] = cp->objectives[0] + 1000 * cp->objectives[1];


}

int eval_bus(CANDIDATE *cp, int o)
{

  int i, j;

  cp->maxon = 0;
  for(i=0;i<allfleet;i++)
    eval_busplan(cp,&cp->busplans[i], o);

  zapobjectives(cp->objectives);  


  for(i=0;i<allfleet;i++)
    for(j=0;j<5;j++)
      {
	cp->objectives[j] += cp->busplans[i].objectives[j];
      }
}

int init_build_parcels(CANDIDATE *cp, int np, int *parcels)
{

  // add these parcels to this candidate in best way possible

  int nparc = np;

  int v[MAXLOCS]; // van
  int p[MAXLOCS];   //pos
  //  float q[MAXLOCS]; // quality

  int pos, pi, vi, posi, ploc;
  float t, q;
  BEST3 b3, *pb3 = &b3;

  
  while(1)
    {
      // find best place for a parcel in existing group

      if(nparc==0) break;      

      pb3->n = 0;      
      for(pi= 0; pi < nparc; pi++)
	{
	  ploc = parcels[pi];
	  
	  for(vi=0;vi<_parcel_vans; vi++)
	    for(posi=0;posi<=cp->vanplans[vi].ndels; posi++)
	      {		
		// SOMETHING ABOUT PARCEL CONSTRAINT !!!!
		q = addity(cp,ploc,vi,posi,__PARCEL);
		//		printf("added %f  for vi %d posi %d\n", q, vi+1, posi);
		addbest3(q,ploc,vi,posi,0,pb3);
	      }	  
	}
      // take one of the best 3 and add it
      //      printstate(pb3);

      int x = rand()%pb3->n;

      add_parcel_from_b3(cp,pb3,x);
      eval_parcel(cp,0);
      for(pi=0;pi<nparc;pi++) if(parcels[pi] == pb3->i1[x]) break;
      parcels[pi] = parcels[nparc - 1];
      nparc--;
      
    }

}


int init_build_bus(CANDIDATE *cp, int np, int *dempairs)
{

  // add these parcels to this candidate in best way possible

  int nbd = np;

  int v[MAXLOCS]; // van
  int p[MAXLOCS];   //pos
  //  float q[MAXLOCS]; // quality

  int pos, bdi, vi, posi, posi2, dembase, bpu, bsd;
  float t, q, senspos, lastpos, caft;
  BEST3 b3, *pb3 = &b3;

  
  while(1)
    {
      // find best place for a parcel in existing group

      if(nbd==0) break;      

      pb3->n = 0;      
      for(bdi= 0; bdi < nbd; bdi+=2)
	{
	  bpu = dempairs[bdi];
	  bsd = dempairs[bdi+1];
	  
	  for(vi=0;vi<allfleet; vi++)
	    {
	      senspos = first_sensible_bus_addition(&cp->busplans[vi],bpu,30);
	      lastpos = last_sensible_bus_addition(&cp->busplans[vi],bpu,60);

	      if(lastpos < senspos) {senspos = 0; lastpos = cp->busplans[vi].ndels;}
	      
	      for(posi=senspos; posi<=lastpos; posi++)
		{
		  // what is the best place for bsd? 
		  caft = closest_after(&cp->busplans[vi],bsd,posi+1);

		      q = addity_bus(cp,bpu,bsd,vi,posi,caft);		  
		      addbest3(q,bpu,vi,posi,caft,pb3);
		      /*
		      for(posi2=caft;posi2<=(cp->busplans[vi].ndels + 1); posi2++)
			{
			  q = addity_bus(cp,bpu,bsd,vi,posi,posi2);
		    //		    printf("added %f  for bpu %d  vi %d posi %d %d  (ndels %d)\n", bpu, q, vi, posi, posi2, cp->busplans[vi].ndels);
			  addbest3(q,bpu,vi,posi,posi2,pb3);
			  }*/
		}
	    }
	}
      // take one of the best 3 and add it
      //      printstate(pb3);
      // REDO BELOW .....
      int x = rand()%pb3->n;
      //      printf("  ------------------------------>   going with   %f   bpu %d b %d   %d %d  \n", pb3->q[x], pb3->i1[x], pb3->i2[x], pb3->i3[x], pb3->i4[x]);
      add_bus_from_b3(cp,pb3,x);
      //  print_bus(cp);
      eval_bus(cp,0);
      //      print_bus(cp);
      // finally remove these two from dempairs

      for(bdi=0;bdi<nbd;bdi++)
	if(dempairs[bdi] == pb3->i1[x]) break;

      dempairs[bdi] = dempairs[nbd-2];
      dempairs[bdi+1] = dempairs[nbd-1];
      nbd -= 2;
    }

}


int printstate(BEST3 *b3)
{

  printf("STATE:  nst: %d |", b3->n);
  int i;

  for(i=0;i<b3->n;i++) printf(" %d %d %d %d %.2f | ", b3->i1[i],b3->i2[i],b3->i3[i],b3->i4[i],b3->q[i]);
  printf("\n");

}



int add_parcel_from_b3(CANDIDATE *cp, BEST3 *b3, int x)
{

  VANPLAN *vp = &cp->vanplans[b3->i2[x]];
  insert_parcel(vp,b3->i1[x],b3->i3[x]); // vp  loc  posi

}


int add_bus_from_b3(CANDIDATE *cp, BEST3 *b3, int x)
{

  BUSPLAN *vp = &cp->busplans[b3->i2[x]];
  insert_busdem(vp,b3->i1[x],b3->i3[x],b3->i4[x]); // vp  loc  posi

}

int insert_busdem(BUSPLAN *vp, int bpu, int pos1, int pos2)
{
  int newseq[1000], i, nnseq = 0; //
  int bsd =  busdemands[bpu].before;

  insert_bus(vp, bpu, bsd, pos1, pos2);

}

int insert_busdem_rand(BUSPLAN *vp, int bpu)
{
  int newseq[1000], i, nnseq = 0; //
  int bsd =  busdemands[bpu].before;

  int pos1, pos2, thepos;

  int  senspos = first_sensible_bus_addition(vp,bpu,60);
  int  lastpos = last_sensible_bus_addition(vp,bpu,60);

  if(lastpos < senspos)
    {pos2 = lastpos; senspos = lastpos; lastpos = pos2; pos1 = senspos;}

  if(lastpos == senspos) pos1 = lastpos;
  else
    {
      pos1 = senspos + rand()%((lastpos + 1) - senspos);
    }

  pos2 = closest_after(vp,busdemands[bpu].before,pos1+1);

  int adjust = rand()%7;
  adjust -= 3;
  pos2 += adjust;

  if(pos2 < (pos1+1)) pos2 = pos1+1;
  if(pos2 > (vp->ndels + 1)) pos2 = vp->ndels + 1;
  insert_bus(vp, bpu, bsd, pos1, pos2);
}




int insert_parcel(VANPLAN *vp, int parc, int pos)
{
  int newseq[1000], i, nnseq = 0;

  for(i=0;i<pos;i++) newseq[nnseq++] = vp->sequence[i];
  newseq[nnseq++] = parc;
  for(i=pos;i<vp->ndels;i++) newseq[nnseq++] = vp->sequence[i];  
    
  vp->ndels++;

  for(i=0;i<vp->ndels;i++) vp->sequence[i] = newseq[i];

}

int insert_parcel_rand(VANPLAN *vp, int parc)
{
  int newseq[1000], i, nnseq = 0;

  int pos = rand()%(vp->ndels+1);
  
  for(i=0;i<pos;i++) newseq[nnseq++] = vp->sequence[i];
  newseq[nnseq++] = parc;
  for(i=pos;i<vp->ndels;i++) newseq[nnseq++] = vp->sequence[i];  
    
  vp->ndels++;

  for(i=0;i<vp->ndels;i++) vp->sequence[i] = newseq[i];

}

  
float addity(CANDIDATE *cp, int ploc, int vi, int posi, int type)
{

  VANPLAN *vp = &cp->vanplans[vi];

  float deptop = travtime(0,ploc), ptodep = travtime(ploc,0), add;

  add = parcel_visit_time;
  
  if(vp->ndels == 0) {add +=  deptop + ptodep;}

  else  if(posi==0) 
    {add +=  deptop + travtime(ploc, vp->sequence[0]);}

  else if(posi == vp->ndels)
    { add +=  travtime( vp->sequence[vp->ndels-1], ploc) + ptodep; }    
  else 
    add +=  travtime( vp->sequence[posi-1],ploc) + travtime(ploc,vp->sequence[posi]);


  float tsf = vp->objectives[0];
  float over = (tsf + add) - _parcel_maxmins;

  if(over < 0) over = 0;


  //  printf("  deptop %f  ptodep %f  add %f  tsf  %f  over %f \n", deptop, ptodep, add, tsf, over);
  return add + 1000 * over;
}


int closest_after(BUSPLAN *bp, int bsd, int pos)
{
  int loc  = busdemands[bsd].bindex;
  int nd = bp->ndels;
  int i, spos = nd, thisbp, bestpos = pos;
  float dist, bestdist = -1;
  
  for(i=pos;i<nd; i++)
    {
      thisbp = bp->sequence[i];
      //      if(busdemands[thisbp].before >= 0) continue;
      dist = travtime(loc,busdemands[thisbp].bindex);
      if((i == pos)  || (dist < bestdist))
	{bestdist = dist; bestpos = i;}
    }
  return bestpos;
}  

int first_sensible_bus_addition(BUSPLAN *bp, int bpu, int mins)
{
  float ptime = busdemands[bpu].preftime;
  int nd = bp->ndels;
  int i, spos = nd, thisbp;

  for(i=0;i<nd; i++)
    {
      thisbp = bp->sequence[i];
      if(busdemands[thisbp].before < 1) continue;
      
      if((ptime  - busdemands[thisbp].preftime) <= mins) {spos = i; break;}
    }
  return spos;
}

int last_sensible_bus_addition(BUSPLAN *bp, int bpu, int mins)
{
  float ptime = busdemands[bpu].preftime;
  int nd = bp->ndels;
  int i, spos = 0, thisbp;

  for(i=nd-1;i>=0; i--)
    {
      thisbp = bp->sequence[i];
      if(busdemands[thisbp].before < 1) continue;
      
      if((busdemands[thisbp].preftime - ptime) <= mins) {spos = i; break;}
    }
  return spos;
}

float addity_bus(CANDIDATE *cp, int bpu, int bsd, int vi, int posi, int posi2)
{

  CANDIDATE *cc = &population[popsize];
  BUSPLAN *vp = &cp->busplans[vi];
  BUSPLAN *mp = &cc->busplans[vi];  
  
  
  int floc = busfleets[vp->fleet].bindex;
  int puloc = busdemands[bpu].bindex;
  int sdloc = busdemands[bsd].bindex;  
  
  float deptop = travtime(floc,puloc), ptodep = travtime(sdloc,floc), atob = travtime(puloc,sdloc), add;
  float prepen, postpen;
  
  add = 2 * bus_stop_time;
  
  if(vp->ndels == 0) {add +=  deptop + atob + ptodep;} // add is just the time between

  else
    {
      eval_busplan(cp,vp,0);
      prepen = vp->objectives[3];

      copy_busplan(vp, mp);
      
      insert_bus(mp, bpu, bsd, posi, posi2);	      
      eval_busplan(cp,mp,0);
      
      postpen = mp->objectives[3];      
      
      add =  postpen - prepen;
    }

  return add;
}


int insert_bus(BUSPLAN *vp, int bpu, int bsd, int pos1, int pos2)
{

  int i, newp[1000], nnew = 0;

  for(i=0;i<pos1;i++) newp[nnew++] = vp->sequence[i];
  newp[nnew++] = bpu;
  for(i=pos1;i<vp->ndels;i++) newp[nnew++] = vp->sequence[i];
  for(i=0;i<nnew;i++) vp->sequence[i] = newp[i];
  vp->ndels++;

  nnew = 0;
  for(i=0;i<pos2;i++) newp[nnew++] = vp->sequence[i];
  newp[nnew++] = bsd;
  for(i=pos2;i<vp->ndels;i++) newp[nnew++] = vp->sequence[i];
  for(i=0;i<nnew;i++) vp->sequence[i] = newp[i];
  vp->ndels++;  
  

}



int eval_busplan(CANDIDATE *cp, BUSPLAN *vp, int o)
{

  int bd, fl = vp->fleet;
  int floc = busfleets[fl].bindex;
  int puloc;
  int sdloc;
  
  float add; 
  float prepen, postpen,  tt, ptime, passwait = 0, buswait = 0, passlate;


  float tstart, tride = 0, tnow = bus_earlystart, latest, latedrop = 0, km = 0;

  zapobjectives(vp->objectives);

  tstart = tnow;

  int i, ptype, loc, lastloc = floc, parcelson = 0, overload = 0, folkon = 0;

  
  for(i=0; i < vp->ndels; i++)
    { bd = vp->sequence[i];
      loc = busdemands[bd].bindex;
      tt = travtime(lastloc, loc);
      km += travdist(lastloc,loc);
      // now figure time
      tnow += tt; // but will we have to wait  ?

      ptime = busdemands[bd].preftime;
      ptype = busdemands[bd].type;
      latest = busdemands[bd].latesttime;

      if(i==0)
	{ // my need to adjust tstart / tnow based on bus start time constraints 
	  tstart = ptime - tt; tnow = ptime;
	  // OK, but this may be too early
	  if(tstart  < bus_earlystart)  {tstart = bus_earlystart; tnow = tstart + tt;}
	  else if(tstart > bus_latestart) { tstart = bus_latestart; tnow = tstart + tt;}
	}

      // OK what's the damage?

      if(ptime > 0) // shows it's a consideration 
	{
	  if(tnow < ptime) {buswait += (ptime - tnow); tnow = ptime;}
	  else if (tnow > ptime)
	    {
	      passlate = (tnow - ptime) - passenger_late; if(passlate <0) passlate = 0;

	      if(ptype < 50) {passwait += passlate;}
	    }
	}

      // now look at latesttime
      if(tnow > latest)
	{
	  latedrop += (tnow - latest);
	}

      if(ptype < 100) tnow += bus_stop_time;
      else tnow += bundle_visit_time;
      
      // add km here
      lastloc = loc;

      if(ptype == 100)
	{parcelson += busdemands[bd].nparcels;
	  if(parcelson > bus_maxparcels) overload += (parcelson - bus_maxparcels);}
      else if(ptype==200) {parcelson -= busdemands[bd].nparcels;}
      else if (ptype == 0) {folkon++; if(folkon>cp->maxon) {cp->maxon = folkon;} if(folkon > bus_maxpassengers) overload += (folkon - bus_maxpassengers); }
      else if (ptype == 1) {folkon--; }

     if(o)
       {
	 if(i==0)
	   {printf("TSTART is %f \n",tstart);
	   }
	 printf("-->TNOW is %f type %d *noe %f  nol %f -- %f %f \n",tnow, ptype, ptime, latest, passlate, tnow - latest);

       }
    }
  // add the back to base
  tt = travtime(lastloc, floc);
  km += travdist(lastloc,floc);
  
  tnow += tt;

  
  float shift = tnow - tstart, ot;


  vp->objectives[0] = shift;
  vp->objectives[1] = passwait;
  vp->objectives[2] = latedrop;
  if(shift > bus_maxmins)
    vp->objectives[2] += (ot = (shift - bus_maxmins));
  vp->objectives[3] = vp->objectives[0] + 50 * vp->objectives[1] + 100 * vp->objectives[2] + 1000 * overload; 

  float cost = km * _permile;
  cost += (shift * _perhour)/60;
  vp->objectives[4] = cost;

  if(o)   printf("FINAL TNOW is %f   over %f (%f %f)  passwait %f latedrop %f overload %d \n",tnow,ot,shift, bus_maxmins,passwait,latedrop, overload);

  
}




int addbest3(float q, int i1, int i2, int i3, int i4, BEST3 *b3)
{

  if(b3->n == 0)
    {
      applyq3(q,i1,i2,i3,i4,0,b3);
      b3->n = 1;
      return;
    }

  if(b3->n == 1)
    {
      if(q < b3->q[0])
	{
	  shiftb3(0,1,b3);	  
	  applyq3(q,i1,i2,i3,i4,0,b3);
	}
      else 
	{
	  applyq3(q,i1,i2,i3,i4,1,b3);
	}
      b3->n = 2;
      return;
    }

  if(b3->n == 2)
    {
      if(q < b3->q[0])
	{
	  shiftb3(1,2,b3);
	  shiftb3(0,1,b3);	  
	  applyq3(q,i1,i2,i3,i4,0,b3);
	}
      else if (q < b3->q[1])
	{	
	  shiftb3(1,2,b3);  
	  applyq3(q,i1,i2,i3,i4,1,b3);
	}
      else 
	{	
	  applyq3(q,i1,i2,i3,i4,2,b3);
	}
      
      b3->n = 3;
      return;
    }

  // if here, we already have 3

      if(q < b3->q[0])
	{
	  shiftb3(1,2,b3);
	  shiftb3(0,1,b3);	  
	  applyq3(q,i1,i2,i3,i4,0,b3);
	}
      else if (q < b3->q[1])
	{	
	  shiftb3(1,2,b3);  
	  applyq3(q,i1,i2,i3,i4,1,b3);
	}
      else if (q < b3->q[2])
	{	
	  applyq3(q,i1,i2,i3,i4,2,b3);
	}
}

int applyq3(float q,int i1,int i2,int i3,int i4, int to, BEST3 *b3)
{

  b3->q[to] = q;  b3->i1[to] = i1;       b3->i2[to] = i2;       b3->i3[to] = i3;     b3->i4[to] = i4; 
}

int shiftb3(int from, int to, BEST3 *b3)
{

  b3->q[to] = b3->q[from];
  b3->i1[to] = b3->i1[from];
  b3->i2[to] = b3->i2[from];  
  b3->i3[to] = b3->i3[from];
  b3->i4[to] = b3->i4[from];

}


int get_parcel_matrix(void)
{
  int i, j;
  LOCATION *li, *lj;
  
  for(i=0;i<nuselocs;i++)
    for(j=0;j<nuselocs;j++)
      {
	if(i==j)
	  {dmat[i][j] = dmat[j][i] = tmat[j][i]= tmat[i][j] = 0; continue;}

	li = &uselocs[i];
	lj = &uselocs[j];	
	
	dmat[i][j] = dmat[j][i] =   HDMULTIPLIER * lldistance(li->lat, li->lon, lj->lat, lj->lon, 'M');

	tmat[i][j] = tmat[j][i] =  60 *  (dmat[i][j] / parcelvan_speed); // parcel speed is in mph - so we mult 60 to get minutes
      }
}

int printparams(void)
{

  int i;

  for(i=0;i<nbusfleets;i++)
    {
      printf("bus fleet %d:  %f %f   %d buses   %f %f (%f)  cpm %f cph %f\n", i+1, busfleets[i].lat,busfleets[i].lon,busfleets[i].n,
	     bus_earlystart, bus_latestart, bus_maxmins, bus_permile, bus_perhour);
    }

}

/****************************************************************************************************\
    print first and last n, or all if n < 0
\****************************************************************************************************/
int printlocs(int n)
{
  int i;

  for(i=0;i<nlocs;i++)
    {
      if((n>0) && ((i < n) || (i > (nlocs - n))))
	printf(" %d  %d %f %f  ->  %s\n", locs[i].type, locs[i].partype, locs[i].lat,locs[i].lon, locs[i].str);
    }

}


int getbuslocs(char *fname)
{

  FILE *f = fopen(fname,"r");

  int r, t;
  float lat, lon;
  char field[1000];
  
  while(1)
    {
      r = getfield(f,field);      if(r==0) break;       
      t = atoi(field);
      r = getfield(f,field);      if(r==0) break; // finished
      lat = atof(field);
      r = getfield(f,field);      if(r==0) break; // finished      
      lon = atof(field);
      addloc(t, 5, lat,lon,"");
    }
  fclose(f);
}




int getparcellocs(char *fname)
{

  FILE *f = fopen(fname,"r");

  int r, t;
  float lat, lon;
  char field[1000];
  
  while(1)
    {
      r = getfield(f,field);      if(r==0) break;       

      if(!strcmp(field,"automate"))
	{
	  r = getfield(f,field);
	  printf(" after automate : got %s\n", field);
	  
	  lms_density  = atoi(field); // postcodes per lm
	  lastmile_radius = 1000;
	  
	  printf("  .. OK, I have been instructed to automatically determine lastmile sites at a density of %d postcodes per site\n",
		 lms_density);
	  return 1;
	}
      
      t = atoi(field);
      r = getfield(f,field);      if(r==0) break; // finished
      lat = atof(field);
      r = getfield(f,field);      if(r==0) break; // finished      
      lon = atof(field);
      r = getstrfield(f,field);      if(r==0) break; // finished      

      addloc(2, 10, lat,lon, field);
      nparcellocs++;

      
      LOCATION *lp = &locs[nlocs-1];

      r = getfield(f,field);      if(r==0) break;       
      t = atoi(field);
      lp->nvans = t;
      r = getfield(f,field);      if(r==0) break;       
      t = atoi(field);
      lp->maxparcels = t;
      r = getfield(f,field);      if(r==0) break;       
      lon = atof(field);
      lp->maxmins = lon;
      r = getfield(f,field);      if(r==0) break;       
      lon = atof(field);
      lp->earlystart = lon;      
      r = getfield(f,field);      if(r==0) break;       
      lon = atof(field);
      lp->loadtime = lon;      
      r = getfield(f,field);      if(r==0) break;       
      lon = atof(field);
      lp->speedfactor = lon;
      r = getfield(f,field);      if(r==0) break;       
      lon = atof(field);
      lp->permile = lon;
      r = getfield(f,field);      if(r==0) break;       
      lon = atof(field);
      lp->perhour = lon;                  
    }
  fclose(f);
}



int getbusdem(char *fname)
{

  FILE *f = fopen(fname,"r");

  int r, t;
  float lat, lon, lat2, lon2;
  char field[1000];
  char tstr[100];
  
  while(1)
    {
      r = getfield(f,field);      if(r==0) break; // finished
      lat = atof(field);
      r = getfield(f,field);      if(r==0) break; // finished      
      lon = atof(field);

      r = getfield(f,field);      if(r==0) break; // finished      
      strcpy(tstr,field);

      r = getfield(f,field);      if(r==0) break; // finished
      lat2 = atof(field);
      r = getfield(f,field);      if(r==0) break; // finished      
      lon2 = atof(field);      

      addbusdem(lat, lon, tstr, lat2, lon2);
    }
  fclose(f);
}




int getobpdem(char *fname)
{

  FILE *f = fopen(fname,"r");

  int r, t;
  float lat, lon, lat2, lon2, v, w;
  char field[1000];
  char tstr[100];
  
  while(1)
    {
      // time lat lon  lat lon v w
       r = getfield(f,field);      if(r==0) break; // finished
      lat = atof(field);
      r = getfield(f,field);      if(r==0) break; // finished      
      lon = atof(field);

      r = getfield(f,field);      if(r==0) break; // finished      
      v = atof(field);

      r = getfield(f,field);      if(r==0) break; // finished      
      w = atof(field);                  

      addobpdem(lat, lon, v, w);
    }
  fclose(f);
}

/*
bus_fleet,  22, 51.12345, 0.321
bus_fleet,  14, 50.42345, 1.321
bus_times, 06:00, 09:00, 08:00
bus_costs,  3.20, 10.50
passenger_window, 40, 5

parcel_depot,  50.4567, -1.45623
parcel_vans,  20, 500, 6
parcel_costs,  12,  4

last_mile_method , trolley 
last_mile_costs, 3, 8
last_mile_speed, 5  -- HERE  

*/

int getconfig(char *fname)
{

  FILE *f = fopen(fname,"r");

  int r, t, nbuses;
  float lat, lon, lat2, lon2, v, w;
  char field[1000];
  char tstr[100];
  
  while(1)
    {
      // time lat lon  lat lon v w
       r = getfield(f,field);      if(r==0) break; // finished
       if(!strcmp(field,"bus_fleet"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   nbuses = atoi(field);

	   r = getfield(f,field);      if(r==0) break; // finished      
	   lat = atof(field);

	   r = getfield(f,field);      if(r==0) break; // finished      
	   lon = atof(field);

	   addbusfleet(nbuses,lat, lon);
	   continue;
	 }

       if(!strcmp(field,"bus_times"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   bus_earlystart = getmins(field);
	   r = getfield(f,field);      if(r==0) break; // finished      
	   bus_latestart = getmins(field);	   
	   r = getfield(f,field);      if(r==0) break; // finished      
	   bus_maxmins = getmins(field);	   	  
	   continue;
	 }

       if(!strcmp(field,"bus_costs"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   bus_permile = atof(field);
	   r = getfield(f,field);      if(r==0) break; // finished      
	   bus_perhour = atof(field);	   
	   continue;
	 }

       if(!strcmp(field,"bus_speedfactor"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   bus_speedfactor = atof(field);
	   continue;
	 }
       if(!strcmp(field,"lastmile_radius"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      


	   lastmile_radius = atof(field);
	   // otherwise we ignore this because there is no fixed radius
	   continue;
		      
	 }

       if(!strcmp(field,"parcel_speedfactor"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   parcel_speedfactor = atof(field);
	   continue;
	 }
       if(!strcmp(field,"parcel_visit_time"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   parcel_visit_time = atof(field);
	   continue;
	 }
       if(!strcmp(field,"bundle_visit_time"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   bundle_visit_time = atof(field);
	   continue;
	 }                     
       if(!strcmp(field,"passenger_window"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   passenger_early = atoi(field);
	   r = getfield(f,field);      if(r==0) break; // finished      
	   passenger_late = atoi(field);	   

	   continue;
	 }              

       if(!strcmp(field,"parcel_depot"))
	 {

	   r = getfield(f,field);      if(r==0) break; // finished      
	   parcel_lat = atof(field);

	   r = getfield(f,field);      if(r==0) break; // finished      
	   parcel_lon = atof(field);

	   continue;
	 }

       if(!strcmp(field,"parcel_vans"))
	 {

	   r = getfield(f,field);      if(r==0) break; // finished      
	   parcel_vans = atoi(field);

	   r = getfield(f,field);      if(r==0) break; // finished      
	   parcel_w = atof(field);

	   r = getfield(f,field);      if(r==0) break; // finished      
	   parcel_v = atof(field);	   

	   continue;
	 }
       if(!strcmp(field,"parcel_costs"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   parcel_permile = atof(field);
	   r = getfield(f,field);      if(r==0) break; // finished      
	   parcel_perhour = atof(field);	   
	   continue;
	 }

       if(!strcmp(field,"last_mile_method"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   strcpy(lastmile,field);
	   continue;
	 }       

       if(!strcmp(field,"last_mile_costs"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   lastmile_permile = atof(field);
	   r = getfield(f,field);      if(r==0) break; // finished      
	   lastmile_perhour = atof(field);	   
	   continue;
	 }

       if(!strcmp(field,"last_mile_speed"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   lastmile_speed = atof(field);
	   continue;
	 }
       if(!strcmp(field,"parcelvan_maxparcels"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   parcelvan_maxparcels = atoi(field);
	   continue;
	 }
       if(!strcmp(field,"bus_maxpassengers"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   bus_maxpassengers = atoi(field);
	   continue;
	 }
       if(!strcmp(field,"bus_maxparcels"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   bus_maxparcels = atoi(field);
	   continue;
	 }
       if(!strcmp(field,"bundle_drop"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   bundle_drop = atof(field);
	   continue;
	 }
       if(!strcmp(field,"bundle_available"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   bundle_available = atof(field);
	   continue;
	 }
       if(!strcmp(field,"bundle_available_lastmile"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   bundle_available_lastmile = atof(field);
	   continue;
	 }              
       if(!strcmp(field,"busdepot_arrival"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   busdepot_arrival = atof(field);
	   continue;
	 }
       if(!strcmp(field,"parcel_earlystart"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   parcel_earlystart = atof(field);
	   continue;
	 }                                   
       if(!strcmp(field,"parcel_loadtime"))
	 {
	   r = getfield(f,field);      if(r==0) break; // finished      
	   parcel_loadtime = atof(field);
	   continue;
	 }                                   
	        
    }
  fclose(f);
}

int addbusfleet(int n, float lat, float lon)
{
  BUSFLEET *b = &busfleets[nbusfleets++];
  b->n = n;
  b->lat = lat;
  b->lon = lon;


  addloc(__BUSDEPOT, __BUSDEPOT, lat, lon, "");
  b->lindex = nlocs-1;
}

int addloc(int t, int pt, float lat,float lon, char *str)
{
  LOCATION *lp = &locs[nlocs++];

  lp->type = t;
  lp->partype = pt;  
  lp->lat = lat;
  lp->lon = lon;  

  strcpy(lp->str,str);
}


int addbusdem(float lat, float lon, char *tstr, float lat2, float lon2)
{
  PACKAGE *p = &packages[npackages++];


  p->inplay = 1;
  p->id = npackages-1;
  p->type = __BUS; // bus demand
  p->arrtime = -1;
  p->nparcels = 1;
  p->isbundle = 0;
  
  addloc(3,22,lat,lon,"");
  p->pu_raw = nlocs - 1; // the one we just added
  
  addloc(4,23,lat2,lon2,"");  
  p->sd_raw = nlocs - 1; // the one we just added

  p->pu_open = p->pu_close = getmins(tstr);

}


int addsetdownlocker(int lin, int bin, LOCATION *lp)
{
  SDLOCKER *sp;

  sp = &sdlockers[nsdlockers++];
  sp->lindex = lin;
  sp->bindex = bin;

  sp->fleet = closest_busfleet(bin);
  sp->nassigned = 0;
  sp->id = nsdlockers - 1;

  sp->nvans = lp->nvans;
  sp->maxparcels = lp->maxparcels;
  sp->maxmins = lp->maxmins;    
  sp->earlystart = lp->earlystart;
  sp->loadtime = lp->loadtime;
  sp->speedfactor = lp->speedfactor;
  sp->perhour = lp->perhour;
  sp->permile = lp->permile;  

}

int closest_busfleet(int bindex)
{
  int i, c = -1;
  float dist, bestdist;
  
  for(i=0;i<nbusfleets;i++)
    {
      dist = travtime(bindex,busfleets[i].bindex);
      if(i==0) {c = i; bestdist = dist;}
      else if(dist < bestdist) {c = i; bestdist = dist;}
    }
  return c;
}


  
int addpickuplocker(int lin, int bin)
{
  SDLOCKER *sp;

  sp = &pulockers[npulockers++];
  sp->lindex = lin;
  sp->bindex = bin;  
}


// add a parcel demand 

int addobpdem(float lat, float lon, float v, float w)
{
  PACKAGE *p = &packages[npackages++];


  p->inplay = 1;
  p->id = npackages-1;
  p->arrtime = -1;
  p->nparcels = 1;
  p->isbundle = 0;
 
  p->type = __PARCEL; // parcel demand

  addloc(28,100,lat,lon,""); // parcel demand
  p->sd_raw = nlocs - 1; // the one we just added -- already a  snap
  
  p->v = v; p->w = w;

}


int getmins(char *hhmm)
{

  char h[3],m[3];
  int min;
  
  h[0] = hhmm[0];  h[1] = hhmm[1]; h[2] = '\0';

  min = 60 *atoi(h);

  m[0] = hhmm[3];  m[1] = hhmm[4]; m[2] = '\0';  

  min += atoi(m);

  return min;
}


   

int getfield(FILE *f, char *field)
{

  char fld[1000], i=0;
  int state = 0;  // moving over whitespace looking for non-comma content
  int r;
  char c;
  
  while(1)  // keep going over whitespace
    {
      r = fscanf(f,"%c", &c);      if(r<0) return 0;
      if(c==' ') continue;
      if(c=='\t') continue;
      if(c=='\n') continue;
      if(c=='\r') continue;            

      if(c==',') continue; // empty field carry on
      fld[i++] = c;  state = 1; break;  // picking up content
    }

  while(1)
    {
      r = fscanf(f,"%c", &c);      if(r<0) break; // was RETURN 0
      if( (c==' ')  ||
          (c=='\t') || 
          (c=='\n') || 
          (c=='\r') ||
          (c==',')) break; //  end of field
      fld[i++] = c;   // picking up content
    }
  fld[i++] = '\0';

  strcpy(field,fld);
}



int getstrfield(FILE *f, char *field)
{

  char fld[1000], i=0;
  int state = 0;  // moving over whitespace looking for non-comma content
  int r;
  char c;
  
  while(1)  // keep going over whitespace
    {
      r = fscanf(f,"%c", &c);      if(r<0) return 0;
      if(c==' ') continue;
      if(c=='\t') continue;
      if(c=='\n') continue;
      if(c=='\r') continue;            

      if(c==',') continue; // empty field carry on
      if(c=='"') {  state = 1; break;} // now to get content 
    }

  while(1)
    {
      r = fscanf(f,"%c", &c);      if(r<0) return 0;
      if (c=='"') {  fld[i++] = '\0';  break;}
      fld[i++] = c;   // picking up content
    }
  strcpy(field,fld);
}


double travtime(int a, int b)
{

  double r = tmat[a][b];

  return  r/_speedfactor;

}

double travdist(int a, int b)
{

  double r = dmat[a][b];

  return r;
}


double lldistance(double lat1, double lon1, double lat2, double lon2, char unit)
{
  double theta, dist;
  if ((lat1 == lat2) && (lon1 == lon2)) {
    return 0;
  }
  else {
    theta = lon1 - lon2;
    dist = sin(deg2rad(lat1)) * sin(deg2rad(lat2)) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * cos(deg2rad(theta));
    dist = acos(dist);
    dist = rad2deg(dist);
    dist = dist * 60 * 1.1515;
    switch(unit) {
    case 'M':
      break;
    case 'K':
      dist = dist * 1.609344;
      break;
    case 'N':
      dist = dist * 0.8684;
      break;
    }
    return (dist);
  }
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  This function converts decimal degrees to radians             :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

 double deg2rad(double deg) {
  return (deg * pi / 180);
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  This function converts radians to decimal degrees             :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
double rad2deg(double rad) {
  return (rad * 180 / pi);
}

int zapobjectives(double *o)
{
  int i;
  for (i = 0; i < 5 ; i++) o[i] = 0;

}



int InsidePolygon(POINT *polygon, int N, POINT *p)
{
  int counter = 0;
  int i;
  double xinters;
  POINT *p1,*p2;

  p1 = &polygon[0];
  for (i=1;i<=N;i++) {
    p2 = &polygon[i % N];
    if (p->y > MIN(p1->y,p2->y)) {
      if (p->y <= MAX(p1->y,p2->y)) {
        if (p->x <= MAX(p1->x,p2->x)) {
          if (p1->y != p2->y) {
            xinters = (p->y-p1->y)*(p2->x-p1->x)/(p2->y-p1->y)+p1->x;
            if (p1->x == p2->x || p->x <= xinters)
              counter++;
          }
        }
      }
    }
    p1 = p2;
  }

  if (counter % 2 == 0)
    return(OUTSIDE);
  else
    return(INSIDE);
}
