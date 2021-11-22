#ifndef _MSHMET_H
#define _MSHMET_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>

#include "chrono.h"
#include "memory.h"
#include "libmesh5.h"
#include "eigenv.h"

#define MS_VER   "2.3c"
#define MS_REL   "Aug. 25, 2009"
#define MS_CPY   "Copyright (c) LJLL, 2007-09"
#define MS_STR   "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"

#define PRECI    1.0
#define EPS      1.e-6
#define EPS1     1.e-15
#define EPST    -1.e-2
#define EPSR     1.e+2
#define LONMAX   4096

#define CTE2D    2.0 / 9.0
#define CTE3D    9.0 / 32.0
  
#define MS_MAX(a,b)   ( ((a) < (b)) ? (b) : (a) )
#define MS_MIN(a,b)   ( ((a) < (b)) ? (a) : (b) )

char idir[5];

typedef struct {
  double         aire,rins;
  double         c[3];
  int            s,nv,mark;
  unsigned char  b,h;
} Point;
typedef Point * pPoint;

typedef struct {
  double  aire;
  int     v[3];
  int     mark;
} Tria;
typedef Tria * pTria;

typedef struct {
  int     v[4];
  int     mark;
} Tetra;
typedef Tetra * pTetra;

typedef struct {
  double   delta;
  double   min[3],max[3];
  float    eps,hmin,hmax,width;
  int      nnu,nsol,nlis;
  char     imprim,ddebug,iso,bin,metric,ls;
} Info;

typedef struct {
  int      np,nt,ne,ver,dim;
  int     *adja,mark;
  char    *name,*mname;

  pPoint   point;
  pTria    tria;
  pTetra   tetra;
  Info     info;
} Mesh;
typedef Mesh * pMesh;

typedef struct {
  int         np,ver,dim,type,size,typtab[GmfMaxTyp];
  double     *sol,*met,umin,umax;
  char       *name,*outn;
} Sol;
typedef Sol * pSol;

typedef struct {
  double         grd[3];
  double         hes[6];
} Deriv;
typedef Deriv * pDeriv;


/* prototypes */
int eigen_3d(int ,double *,double *,double v[3][3]);
int eigen_2d(double *,double *,double vp[2][2]);

int loadMesh(pMesh ,char *);
int loadSol(pSol ,Info *,char *);
int loadMetric(pSol ,Info *,char *);
int saveMet(pSol ,Info *,char *);
int saveSol(pSol ,Info *,char *);
int outder(pMesh ,pDeriv ,int );
int scaleMesh(pMesh ,pSol );
int unscaleMesh(pMesh ,pSol );

int mshme1(pMesh ,pSol );

int    zaldy(pMesh );
int    zaldy1(pSol );
pDeriv zaldy2(int );
int    zaldy3(Info *,pSol );
void   freeMesh(pMesh mesh);

/* function pointers */
int   boulep_3d(pMesh ,int ,int ,int *);
int   boulep_2d(pMesh ,int ,int ,int *);
int   hashel_3d(pMesh );
int   hashel_2d(pMesh );
int   gradLS_3d(pMesh ,pSol ,int ,int ,double *);
int   gradLS_2d(pMesh ,pSol ,int ,int ,double * );
int   gradLS_s(pMesh ,pSol ,int ,int ,double * );
int   hessLS_3d(pMesh ,pSol ,int ,int ,double *,double *);
int   hessLS_2d(pMesh ,pSol ,int ,int ,double *,double *);
int   hessLS_s(pMesh ,pSol ,int ,int ,double *,double *);
int   avgval_3d(pMesh ,pDeriv ,int ,double *);
int   avgval_2d(pMesh ,pDeriv ,int ,double *);
int   avgval_s(pMesh ,pDeriv ,int ,double *);
int   clsval_3d(pMesh ,pDeriv ,int ,double *);
int   clsval_2d(pMesh ,pDeriv ,int ,double *);
int   clsval_s(pMesh ,pDeriv ,int ,double *);
double getSol_3d(pSol ,int ,int );
double getSol_2d(pSol ,int ,int );
int   nrmhes_3d(pMesh ,pSol ,pDeriv ,int );
int   nrmhes_2d(pMesh ,pSol ,pDeriv ,int );
int   laplac_3d(pMesh ,pDeriv );
int   laplac_2d(pMesh ,pDeriv );
int   laplac_s(pMesh ,pDeriv );
int   redsim_3d(double *,double *,double *);
int   redsim_2d(double *,double *,double *);
int   defmet_3d(pMesh ,pSol ,pDeriv ,int );
int   defmet_2d(pMesh ,pSol ,pDeriv ,int );
int   metrLS_3d(pMesh mesh,pSol ,pDeriv der);
int   metrLS_2d(pMesh mesh,pSol ,pDeriv der);
int   lissag_3d(pMesh mesh,pSol sol,int is);
int   lissag_2d(pMesh mesh,pSol sol,int is);

int   (*boulep)(pMesh ,int ,int ,int *);
int   (*hashel)(pMesh );
int   (*gradLS)(pMesh ,pSol ,int ,int ,double *);
int   (*hessLS)(pMesh ,pSol ,int ,int ,double *,double *);
int   (*avgval)(pMesh ,pDeriv ,int ,double *);
int   (*clsval)(pMesh ,pDeriv ,int ,double *);
int   (*nrmhes)(pMesh ,pSol ,pDeriv ,int );
int   (*laplac)(pMesh mesh,pDeriv der);
int   (*redsim)(double *,double *,double *);
int   (*defmet)(pMesh ,pSol ,pDeriv ,int );
double (*getSol)(pSol ,int ,int );
int   (*metrLS)(pMesh mesh,pSol ,pDeriv der);
int   (*lissag)(pMesh ,pSol , int );

#endif
