#include "mshmet.h"


int mshme1(pMesh mesh,pSol sol) {
  pPoint   ppt;
  pDeriv   der;
  int      j,k,ret,nc,nex,ney;

  /* mem alloc */
  der = (Deriv *)M_calloc(sol->np+1,sizeof(Deriv),"deriv");
  assert(der);

  if ( mesh->info.metric == 0 ) {
    if ( mesh->info.iso )
      sol->met = (double*)M_calloc(sol->np+1,sizeof(double),"zaldy3.met");
    else {
      if ( mesh->dim == 2 )
        sol->met = (double*)M_calloc(sol->np+1,3*sizeof(double),"zaldy3.met");
      else
        sol->met = (double*)M_calloc(sol->np+1,6*sizeof(double),"zaldy3.met");
    }
    assert(sol->met);
  }
  if ( mesh->info.ls )  sol->type = 1;

  /* solution smoothing */
  if ( mesh->info.nlis ) {
		for (j=0; j<mesh->info.nlis; j++)
			lissag(mesh,sol,0);
		/*saveSol(sol,&mesh->info,"titi.sol");*/
  }

  for (j=0; j<sol->type; j++) {
    if ( mesh->info.nsol > -1 && mesh->info.nsol != j )  continue;
    if ( mesh->info.imprim < 0 )
      fprintf(stdout,"     Solution %d: %d vertices\n",j+1,mesh->np);

    /* compute gradient */
    nex = ney = 0;
    for (k=1; k<=mesh->np; k++) {
      if ( !gradLS(mesh,sol,k,j,der[k].grd) )  return(0);
    }

    /* compute hessian */
    for (k=1; k<=mesh->np; k++) {
      ret = hessLS(mesh,sol,k,j,der[k].grd,der[k].hes);
      if ( !ret )  return(0);
      else if ( ret < 0 )  nex++;
    }
    /* correction */
    if ( nex ) {
      ney = nex;
      do {
	      nc  = 0;
 	      nex = 0;
        for (k=1; k<=mesh->np; k++) {
          ppt = &mesh->point[k];
          if ( ppt->h )  continue;
	        nex++;
	        nc += avgval(mesh,der,k,der[k].hes);
        }
      }
      while (nc > 0 && nex > 0);

      /* value of closest vertex */
      if ( nex ) {
	      do {
	        nc = 0;
	        for (k=1; k<=mesh->np; k++) {
            ppt = &mesh->point[k];
	          if ( ppt->h )  continue;
	          nc += clsval(mesh,der,k,der[k].hes);
          }
          nex -= nc;
	      }
	      while (nc > 0 && nex > 0);
      }
    }

    if ( mesh->info.imprim < 0 && nex+ney )
      fprintf(stdout,"  %%%% %d corrected  %d unknowns\n",ney,nex);

    /* norm + avg hessian */
    if ( !nrmhes(mesh,sol,der,j) )  return(0);
    if ( !laplac(mesh,der) )        return(0);

    /* define metric */
    if ( mesh->info.ddebug )  outder(mesh,der,j);
	  if ( mesh->info.ls )  break;
    else if ( !defmet(mesh,sol,der,j) )  return(0);
  }

  /* levelset metric */
  if ( mesh->info.ls )
    if ( !metrLS(mesh,sol,der) )  return(0);

  M_free(der);
  return(1);
}
