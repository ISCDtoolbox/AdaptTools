/* mshmet.c
 * compute metric based on hessian
 *
 * Pascal FREY, LJLL
 * Copyright (c) LJLL, 2004.
*/

#include "mshmet.h"
#include "compil.date"

char     idir[5]     = {0,1,2,0,1};
mytime   ctim[TIMEMAX];


static void excfun(int sigid) {
  fprintf(stdout,"\n Unexpected error:");  fflush(stdout);
  switch(sigid) {
    case SIGABRT:
      fprintf(stdout,"  Abnormal stop\n");  exit(1);
    case SIGFPE:
      fprintf(stdout,"  Floating-point exception\n"); exit(1);
    case SIGILL:
      fprintf(stdout,"  Illegal instruction\n"); exit(1);
    case SIGSEGV:
      fprintf(stdout,"  Segmentation fault\n");  exit(1);
    case SIGTERM:
    case SIGINT:
      fprintf(stdout,"  Program killed\n");  exit(1);
  }
  exit(1);
}


static void usage(char *prog) {
  fprintf(stdout,"\n usage: %s filein[.mesh] [solin[.sol]] [fileout.sol] -eps x -hmin y -hmax z -v -iso -norm\n",prog);
  
  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"-d      Turn on debug mode\n");
  fprintf(stdout,"-h      Print this message\n");
  fprintf(stdout,"-ls     Build levelset metric\n");
  fprintf(stdout,"-v [n]  Tune level of verbosity\n");
  fprintf(stdout,"-m file Use metric file\n");

  fprintf(stdout,"\n** Specific options : \n");
  fprintf(stdout,"  -eps :  tolerance\n");
  fprintf(stdout,"  -hmin:  min size\n");
  fprintf(stdout,"  -hmax:  max size\n");
  fprintf(stdout,"  -iso :  isotropic\n");
  fprintf(stdout,"  -w   :  relative width for LS (0<w<1)\n");
  fprintf(stdout,"  -n[i]:  normalization (level i), default=0\n");
  fprintf(stdout,"  -s n :  consider solution n (only)\n");

  fprintf(stdout,"\n** DEFAULT.hmet file allows to store parameter values\n");
  fprintf(stdout,"eps   x\n");
  fprintf(stdout,"hmin  y\n");
  fprintf(stdout,"hmax  z\n");
  fprintf(stdout,"iso\n");
  exit(1);
}


static int parsar(int argc,char *argv[],pMesh mesh,pSol sol) {
  Info    *info;
  int      i;
  char    *ptr;

  /* default */
  info = &mesh->info;
  info->hmin   = 0.01;
  info->hmax   = 1.0;
  info->eps    = 0.01;
  info->nnu    = 1;
  info->iso    = 0;
  info->bin    = 1;
  info->metric = 0;
  info->ls     = 0;
  info->nsol   = -1;
  info->width  = 0.05;
  info->ddebug = 0;
  info->imprim = -99;
	info->nlis   = 0;

  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case 'h':  /* on-line help */
        if ( !strcmp(argv[i],"-hmin") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            info->hmin = atof(argv[i]);
        }
        else if ( !strcmp(argv[i],"-hmax") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            info->hmax = atof(argv[i]);
        }
        else  
          usage(argv[0]);
        break;
      case '?':
        usage(argv[0]);
        break;

      case 'd':  /* debug */
        info->ddebug = 1;
        break;

      case 'e':
        if ( ++i < argc && isdigit(argv[i][0]) )
          info->eps = atof(argv[i]);
        else
          i--;
        break;

      case 'i':
        if ( !strcmp(argv[i],"-iso") )
          info->iso = 1;
        else if ( !strcmp(argv[i],"-in") ) {
		      mesh->name = argv[i];
	      }
        break;

      case 'l':
        if ( !strcmp(argv[i],"-ls") )
          info->ls = 1;
        else if ( !strcmp(argv[i],"-lis") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            info->nlis = atoi(argv[i]);
        }
        break;

      case 'm':
        if ( ++i < argc ) {
          mesh->mname  = argv[i];
          info->metric = 1;
        }
        else {
          usage(argv[0]);
        }
        break;

      case 'n':
        info->nnu = 1;
        if ( ++i < argc && isdigit(argv[i][0]) )
          info->nnu = atoi(argv[i]);
        else
          i--;
        break;

      case 's':
        if ( ++i < argc && isdigit(argv[i][0]) )
          info->nsol = atoi(argv[i]) - 1;
        else
          i--;
        break;

      case 'v':
        if ( ++i < argc ) {
          if ( argv[i][0] == '-' || isdigit(argv[i][0]) )
            info->imprim = atoi(argv[i]);
          else 
            i--;
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          usage(argv[0]);
        }
        break;

      case 'w':
        if ( ++i < argc ) {
          if ( argv[i][0] == '-' || isdigit(argv[i][0]) )
            info->width = MS_MIN(1,MS_MAX(0,atof(argv[i])));
          else 
            i--;
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          usage(argv[0]);
        }
        break;
      default:
        fprintf(stderr,"  Unrecognized option %s\n",argv[i]);
        usage(argv[0]);
      }
    }

    else {
      if ( mesh->name == NULL ) {
        mesh->name = argv[i];
        if ( info->imprim == -99 )  info->imprim = 5;
      }
      else if ( sol->outn == NULL ) {
		sol->outn = argv[i];
	  }
      else {
        fprintf(stdout,"  Argument %s ignored\n",argv[i]);
        usage(argv[0]);
      }
    }
    i++;
  }

  /* check params */
  if ( mesh->info.imprim == -99 ) {
    fprintf(stdout,"\n  -- PRINT (0 10(advised) -10) ?\n");
    fflush(stdin);
    fscanf(stdin,"%d",&i);
    info->imprim = i;
  }

  if ( mesh->name == NULL ) {
    mesh->name = (char *)calloc(128,sizeof(char));
    assert(mesh->name);
    fprintf(stdout,"  -- MESH BASENAME ?\n");
    fflush(stdin);
    fscanf(stdin,"%s",mesh->name);
  }
  sol->name = (char *)calloc(128,sizeof(char));
  assert(sol->name);
  strcpy(sol->name,mesh->name);
  ptr = strstr(sol->name,".mesh");
  if ( ptr ) *ptr = '\0';
  if ( sol->outn == NULL) {
    sol->outn = (char *)calloc(128,sizeof(char));
    assert(sol->outn);
    strcpy(sol->outn,sol->name);
    strcat(sol->outn,".new");
  }

  return(1);
}


int parsop(pMesh mesh,pSol sol) {
  char    *ptr,data[256],key[256];
  float    dummy;
  int      i,ret;
  FILE    *in;

  strcpy(data,sol->name);
  ptr = strstr(data,".sol");
  if ( ptr )  *ptr = '\0';
  strcat(data,".mhes");

  in = fopen(data,"r");
  if ( !in ) {
    strcpy(data,"DEFAULT.hmet");
    in = fopen(data,"r");
    if ( !in )  {
      if ( mesh->info.imprim < 0 )
        fprintf(stdout,"  %%%% DEFAULT VALUES (%g %g %g)\n",
                mesh->info.eps,mesh->info.hmin,mesh->info.hmax);
      return(1);
    }
  }
  fprintf(stdout,"  %%%% %s FOUND\n",data);

  while ( !feof(in) ) {
    ret = fscanf(in,"%s",key);
    if ( !ret || feof(in) )  break;
    for (i=0; i<strlen(key); i++) key[i] = tolower(key[i]);

    if ( !strcmp(key,"hmin") ) {
      fscanf(in,"%f",&dummy);
      mesh->info.hmin = dummy;
    }
    else if ( !strcmp(key,"hmax") ) {
      fscanf(in,"%f",&dummy);
      mesh->info.hmax = dummy;
    }
    else if ( !strcmp(key,"eps") ) {
      fscanf(in,"%f",&dummy);
      mesh->info.eps = dummy;
    }
    else if ( !strcmp(key,"iso") ) {
      mesh->info.iso = 1;
    }
    else if ( !strcmp(key,"norm") ) {
      fscanf(in,"%d",&mesh->info.nnu);
    }
    else if ( key[0] == '#' ) {
      fgets(key,255,in);
    }
    else fprintf(stderr,"  unrecognized keyword : %s\n",key);
  }

  fclose(in);
  return(1);
}


static void stats(pMesh mesh,pSol sol) {
  fprintf(stdout,"     NUMBER OF GIVEN VERTICES   %8d\n",mesh->np);
  if ( mesh->nt )
    fprintf(stdout,"     NUMBER OF GIVEN TRIANGLES  %8d\n",mesh->nt);
  if ( mesh->ne )
    fprintf(stdout,"     NUMBER OF GIVEN TETRAHEDRA %8d\n",mesh->ne);
  fprintf(stdout,"     NUMBER OF GIVEN DATA       %8d\n",sol->np);
}


static void endcod() {
  double   ttot,ttim[TIMEMAX];
  int      k,call[TIMEMAX];

  chrono(OFF,&ctim[0]);
  for (k=0; k<TIMEMAX; k++) {
    call[k] = ctim[k].call;
    ttim[k] = ctim[k].call ? gttime(ctim[k]) : 0.0;
  }
  ttot    = ttim[1]+ttim[2]+ttim[3]+ttim[4];
  ttim[0] = MS_MAX(ttim[0],ttot);

  if ( ttot > 0.01 ) {
    fprintf(stdout,"\n  -- CPU REQUIREMENTS\n");
    fprintf(stdout,"  in/out    %8.2f %%    %3d. calls,   %7.2f sec/call\n",
        100.*ttim[1]/ttim[0],call[1],ttim[1]/(float)call[1]);
    fprintf(stdout,"  analysis  %8.2f %%    %3d. calls,   %7.2f sec/call\n",
        100.*ttim[2]/ttim[0],call[2],ttim[2]/(float)call[2]);
    fprintf(stdout,"  metric    %8.2f %%    %3d. calls,   %7.2f sec/call\n",
        100.*ttim[3]/ttim[0],call[3],ttim[3]/(float)call[3]);
    fprintf(stdout,"  total     %8.2f %%    %3d. calls,   %7.2f sec/call\n",
        100.*ttot/ttim[0],call[0],ttot/(float)call[0]);
  }
  fprintf(stdout,"\n   ELAPSED TIME  %.2f SEC.  (%.2f)\n",ttim[0],ttot);
}


/* set function pointers */
void setfunc(pMesh mesh) {
  if ( !mesh->ne ) {
    if ( mesh->dim == 2 ) {
      boulep = boulep_2d;
      hashel = hashel_2d;
      gradLS = gradLS_2d;
      hessLS = hessLS_2d;
      getSol = getSol_2d;
      avgval = avgval_2d;
      clsval = clsval_2d;
      nrmhes = nrmhes_2d;
      laplac = laplac_2d;
      defmet = defmet_2d;
      redsim = redsim_2d;
      metrLS = metrLS_2d;
      lissag = lissag_2d;
    }
    else {
      boulep = boulep_2d;
      hashel = hashel_2d;
      gradLS = gradLS_s;
      hessLS = hessLS_s;
      getSol = getSol_3d;
      avgval = avgval_s;
      clsval = clsval_s;
      nrmhes = nrmhes_3d;
      laplac = laplac_s;
      defmet = defmet_3d;
      redsim = redsim_3d;
      metrLS = metrLS_3d;
      lissag = lissag_2d;
    }
  }
  else {
    boulep = boulep_3d;
    hashel = hashel_3d;
    gradLS = gradLS_3d;
    hessLS = hessLS_3d;
    getSol = getSol_3d;
    avgval = avgval_3d;
    clsval = clsval_3d;
    nrmhes = nrmhes_3d;
    laplac = laplac_3d;
    defmet = defmet_3d;
    redsim = redsim_3d;
    metrLS = metrLS_3d;
    lissag = lissag_3d;
  }
}


int main(int argc,char **argv) {
  pMesh       mesh;
  pSol        sol;

  fprintf(stdout,"  -- MSHMET, Release %s (%s) \n",MS_VER,MS_REL);
  fprintf(stdout,"     %s\n",MS_CPY);
  fprintf(stdout,"    %s\n",COMPIL);

  /* trap exceptions */
  signal(SIGABRT,excfun);
  signal(SIGFPE,excfun);
  signal(SIGILL,excfun);
  signal(SIGSEGV,excfun);
  signal(SIGTERM,excfun);
  signal(SIGINT,excfun);
  atexit(endcod);

  tminit(ctim,TIMEMAX);
  chrono(ON,&ctim[0]);

  /* mem alloc */
  mesh = (pMesh)M_calloc(1,sizeof(Mesh),"main");
  assert(mesh);
  sol  = (pSol)M_calloc(1,sizeof(Sol),"main");
  assert(sol);

  if ( !parsar(argc,argv,mesh,sol) )  exit(1);

  /* load data */
  if ( mesh->info.imprim )   fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&ctim[1]);
  if ( !loadMesh(mesh,mesh->name) )           exit(1);
  if ( !loadSol(sol,&mesh->info,sol->name) )  exit(1);
  if ( sol->np != mesh->np ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    exit(1);
  }
  if ( mesh->info.metric && mesh->mname )
    if ( !loadMetric(sol,&mesh->info,mesh->mname) )  exit(1);
  if ( !parsop(mesh,sol) )  exit(1);
  setfunc(mesh);
  chrono(OFF,&ctim[1]);
  if ( mesh->info.imprim )  stats(mesh,sol);
  fprintf(stdout,"  -- DATA READING COMPLETED.     %.2f sec.\n",gttime(ctim[1]));

  fprintf(stdout,"\n  %s\n   MODULE MSHMET-LJLL : %s (%s)\n  %s\n",
          MS_STR,MS_VER,MS_REL,MS_STR);

  /* analysis */
  chrono(ON,&ctim[2]);
  if ( mesh->info.imprim )   fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  if ( abs(mesh->info.imprim) > 4 ) {
    fprintf(stdout,"  ** SETTING ADJACENCIES\n");
    fflush(stdout);
  }
  if ( !scaleMesh(mesh,sol) )  return(1);
  if ( !hashel(mesh) )         exit(1);
  chrono(OFF,&ctim[2]);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %.2f sec.\n",gttime(ctim[2]));

  /* metric */
  chrono(ON,&ctim[3]);
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- PHASE 2 : METRIC\n");
  if ( !mshme1(mesh,sol) )  exit(1);
  chrono(OFF,&ctim[3]);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %.2f sec.\n",gttime(ctim[3]));

  fprintf(stdout,"\n  %s\n   END OF MODULE MSHMET \n  %s\n",MS_STR,MS_STR);

  /* save file */
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh->name);
  chrono(ON,&ctim[1]);
  if ( !saveMet(sol,&mesh->info,sol->outn) )  exit(1);
  chrono(OFF,&ctim[1]);
  if ( mesh->info.imprim )  fprintf(stdout,"  -- WRITING COMPLETED\n");

  /* free mem */
  M_free(mesh->point);
  if ( mesh->nt )  M_free(mesh->tria);
  if ( mesh->ne )  M_free(mesh->tetra);
  M_free(mesh->adja);
  M_free(sol->sol);
  M_free(sol->met);
  M_free(sol);

  if ( mesh->info.imprim < -4 || mesh->info.ddebug )  M_memDump();
  return(0);
}
