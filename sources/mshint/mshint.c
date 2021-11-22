#include "mshint.h"
#include "compil.date"

Param    par;


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
  fprintf(stdout,"usage: %s [-v[n]] [-h] file1[.mesh] file2[.mesh]\n",prog);
  
  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"-d      Turn on debug mode\n");
  fprintf(stdout,"-h      Print this message\n");
	fprintf(stdout,"-r      Scaling factor on ROI\n"); 
  fprintf(stdout,"-v[n]  Tune level of verbosity\n");
  fprintf(stdout,"-e      Return element number in mesh1 containing vertex in mesh2\n");

  exit(1);
}


static int parsar(int argc,char *argv[],pMesh mesh1,pSol sol1,pMesh mesh2,pSol sol2) {
  int      i;
  char    *ptr;

  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case 'h':
      case '?':
	    usage(argv[0]);
	    break;

      case 'd':
        par.ddebug = 1;
        break;

      case 'e':
        par.option = 2;
        break;

      case 'r':
        if ( ++i < argc ) {
	        if ( argv[i][0] == '-' || isdigit(argv[i][0]) )
	          par.ray = I_MAX(0.,atof(argv[i]));
	        else 
	          i--;
	      }
	      else {
	        fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
	        usage(argv[0]);
	      }
				break;
      case 'v':
        if ( ++i < argc ) {
	        if ( argv[i][0] == '-' || isdigit(argv[i][0]) )
	          par.imprim = atoi(argv[i]);
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
      if ( mesh1->name == NULL ) {
        mesh1->name = argv[i];
        if ( par.imprim == -99 )  par.imprim = 5;
      } 
      else if ( mesh2->name == NULL )
        mesh2->name = argv[i];
      else {
        fprintf(stdout,"  Argument %s ignored\n",argv[i]);
	    usage(argv[0]);
      }
    }
    i++;
  }

  /* check params */
  if ( mesh1->name == NULL  || par.imprim == -99 ) {
    fprintf(stdout,"\n  -- PRINT (0 10(advised) -10) ?\n");
    fflush(stdin);
    fscanf(stdin,"%d",&i);
    par.imprim = i;
  }

  if ( mesh1->name == NULL ) {
    mesh1->name = (char *)calloc(128,sizeof(char));
    assert(mesh1->name);
    fprintf(stdout,"  -- MESH1 BASENAME ?\n");
    fflush(stdin); 
    fscanf(stdin,"%s",mesh1->name);
  }
  sol1->name = (char *)calloc(128,sizeof(char));
  assert(sol1->name);
  strcpy(sol1->name,mesh1->name);
  ptr = strstr(sol1->name,".mesh");
  if ( ptr ) *ptr = '\0';

  if ( mesh2->name == NULL ) {
    mesh2->name = (char *)calloc(128,sizeof(char));
    assert(mesh2->name);
    fprintf(stdout,"  -- MESH2 BASENAME ?\n");
    fflush(stdin); 
    fscanf(stdin,"%s",mesh2->name);
  }
  sol2->name = (char *)calloc(128,sizeof(char));
  assert(sol2->name);
  strcpy(sol2->name,mesh2->name);

  return(1);
}
 

static void stats(pMesh mesh1,pSol sol1,pMesh mesh2,pSol sol2) {
  fprintf(stdout,"     NUMBER OF GIVEN VERTICES   %8d  %8d\n",mesh1->np,mesh2->np);
  if ( mesh1->nt )
    fprintf(stdout,"     NUMBER OF GIVEN TRIANGLES  %8d  %8d\n",mesh1->nt,mesh2->nt);
  if ( mesh1->ne )
    fprintf(stdout,"     NUMBER OF GIVEN TETRAHEDRA %8d  %8d\n",mesh1->ne,mesh2->ne);
  fprintf(stdout,"     NUMBER OF GIVEN DATA       %8d  %8d\n",sol1->np+sol1->ne,0);
}


static void endcod() {
  double   ttot,ttim[TIMEMAX];
  int      k,call[TIMEMAX];

  chrono(OFF,&par.ctim[0]);
  for (k=0; k<TIMEMAX; k++) {
    call[k] = par.ctim[k].call;
    ttim[k] = par.ctim[k].call ? gttime(par.ctim[k]) : 0.0;
  }
  ttot    = ttim[1]+ttim[2]+ttim[3]+ttim[4];
  ttim[0] = I_MAX(ttim[0],ttot);

  if ( abs(par.imprim) > 5 ) {
    fprintf(stdout,"\n  -- CPU REQUIREMENTS\n");
    fprintf(stdout,"  in/out    %8.2f %%    %3d. calls,   %7.2f sec/call\n",
	    100.*ttim[1]/ttim[0],call[1],ttim[1]/(float)call[1]);
    fprintf(stdout,"  analysis  %8.2f %%    %3d. calls,   %7.2f sec/call\n",
	    100.*ttim[2]/ttim[0],call[2],ttim[2]/(float)call[2]);
    fprintf(stdout,"  interpol  %8.2f %%    %3d. calls,   %7.2f sec/call\n",
	    100.*ttim[3]/ttim[0],call[3],ttim[3]/(float)call[3]);
    fprintf(stdout,"  total     %8.2f %%    %3d. calls,   %7.2f sec/call\n",
	    100.*ttot/ttim[0],call[0],ttot/(float)call[0]);
  }
  fprintf(stdout,"\n   ELAPSED TIME  %.2f SEC.  (%.2f)\n",ttim[0],ttot);
}


/* set function pointers */
int setfunc(int dim) {
  if ( dim == 2 ) {
    newBucket = newBucket_2d;
    buckin    = buckin_2d;
    intpp1    = intpp1_2d;
		intpp0    = intpp0_2d;
		closept   = closept_2d;
		locelt    = locelt_2d;
		hashelt   = hashelt_2d;
		boulep    = boulep_2d;
  }
  else {
    newBucket = newBucket_3d;
    buckin    = buckin_3d;
    intpp1    = intpp1_3d;
		intpp0    = intpp0_3d;
		closept   = closept_3d;
		locelt    = locelt_3d;
		hashelt   = hashelt_3d;
		boulep    = boulep_3d;
  }
  return 1;
}


int main(int argc,char **argv) {
  Mesh       mesh1,mesh2;
  Sol        sol1,sol2;
  int        i;

  fprintf(stdout,"  -- MSHINT, Release %s (%s) \n",I_VER,I_REL);
  fprintf(stdout,"     %s\n",I_CPY);
  fprintf(stdout,"    %s\n",COMPIL);

  /* trap exceptions */
  signal(SIGABRT,excfun);
  signal(SIGFPE,excfun);
  signal(SIGILL,excfun);
  signal(SIGSEGV,excfun);
  signal(SIGTERM,excfun);
  signal(SIGINT,excfun);
  atexit(endcod);

  tminit(par.ctim,TIMEMAX);
  chrono(ON,&par.ctim[0]);

  /* default values */
  memset(&mesh1,0,sizeof(Mesh));
  memset(&sol1,0,sizeof(Sol));
  memset(&mesh2,0,sizeof(Mesh));
  memset(&sol2,0,sizeof(Sol));
  par.imprim = -99;
  par.ddebug = 0;
  par.option = 1;
	par.ray    = 1.2;

  if ( !parsar(argc,argv,&mesh1,&sol1,&mesh2,&sol2) )  return(1);

  /* load data */
  if ( par.imprim )   fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&par.ctim[1]);
  if ( !loadMesh(&mesh1,mesh1.name) )  return(1);
  if ( !loadMesh(&mesh2,mesh2.name) )  return(1);
  if ( !setfunc(mesh1.dim) )  return(1);

  if ( par.option == 1 ) {
    if ( !loadSol(&sol1,sol1.name) )  return(1);
    sol2.dim  = mesh2.dim;
    if ( sol1.np )
      sol2.np = mesh2.np;
    if ( sol1.ne ) {
			sol2.ne = mesh2.dim == 2 ? mesh2.nt : mesh2.ne;
		}
    sol2.ver  = sol1.ver;
    sol2.iter = sol1.iter;
    sol2.time = sol1.time;
    memcpy(sol2.type,sol1.type,2*sizeof(int));
    memcpy(sol2.size,sol1.size,2*sizeof(int));
    memcpy(&sol2.typtab[0],&sol1.typtab[0],sol1.type[0]*sizeof(int));
    memcpy(&sol2.typtab[1],&sol1.typtab[1],sol1.type[1]*sizeof(int));
  }
  else {
    sol2.dim  = mesh2.dim;
    sol2.np   = mesh2.np;
    sol2.ver  = GmfFloat;
    sol2.type[0]   = 1;
    sol2.size[0]   = 1;
    sol2.typtab[0][0] = GmfSca;
  }
  if ( sol2.np ) {
    sol2.valp1 = (double*)M_calloc(sol2.np+1,sol2.size[0]*sizeof(double),"mshint");
    assert(sol2.valp1);
  }
  if ( sol2.ne ) {
    sol2.valp0 = (double*)M_calloc(sol2.ne+1,sol2.size[1]*sizeof(double),"mshint");
    assert(sol2.valp0);
	}
  chrono(OFF,&par.ctim[1]);
  if ( par.imprim )  stats(&mesh1,&sol1,&mesh2,&sol2);
  fprintf(stdout,"  -- DATA READING COMPLETED.     %.2f sec.\n",
          gttime(par.ctim[1]));

  fprintf(stdout,"\n  %s\n   MODULE MSHINT-LJLL : %s (%s)\n  %s\n",
          I_STR,I_VER,I_REL,I_STR);

  /* analysis */
  chrono(ON,&par.ctim[2]);
  if ( par.imprim )   fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  for (i=0; i<mesh1.dim; i++) {
    mesh1.info.min[i] = I_MIN(mesh1.info.min[i],mesh2.info.min[i]);
    mesh1.info.max[i] = I_MAX(mesh1.info.max[i],mesh2.info.max[i]);
  }
  mesh1.info.delta = I_MAX(mesh1.info.delta,mesh2.info.delta);
  memcpy(&mesh2.info,&mesh1.info,sizeof(Info));
  if ( !scaleMesh(&mesh1,&sol1) )  return(1);
  if ( !scaleMesh(&mesh2,0) )      return(1);
  if ( !hashelt(&mesh1) )          return(1);
  chrono(OFF,&par.ctim[2]);
  if ( par.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %.2f sec.\n",gttime(par.ctim[2]));

  /* interpolation */
  if ( par.imprim )   fprintf(stdout,"\n  -- PHASE 2 : INTERPOLATION\n");
  chrono(ON,&par.ctim[3]);
  if ( !mshin1(&mesh1,&sol1,&mesh2,&sol2) )  return(1);
  chrono(OFF,&par.ctim[3]);
  if ( par.imprim )
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %.2f sec.\n",gttime(par.ctim[3]));

  /* save file */
  if ( par.imprim )  fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh2.name);
  chrono(ON,&par.ctim[1]);
  if ( !unscaleMesh(&mesh2,&sol2) )  return(1);
  if ( !saveSol(&sol2,sol2.name) )   return(1);
  chrono(OFF,&par.ctim[1]);
  if ( par.imprim )  fprintf(stdout,"  -- WRITING COMPLETED\n");

  M_free(mesh1.point);
	M_free(mesh1.adja);
	M_free(mesh2.point);
	M_free(mesh2.adja);
  if ( mesh1.dim == 3 ) {
	  M_free(mesh1.tetra);
    M_free(mesh2.tetra);
	}
  else {
    M_free(mesh1.tria);
    M_free(mesh2.tria);
	}
  if ( sol2.np )
    M_free(sol2.valp1);
	if ( sol2.ne )
		M_free(sol2.valp0);
  if ( par.imprim < -4 || par.ddebug )  M_memDump();

  return(0);
}

