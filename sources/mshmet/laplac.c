#include "mshmet.h"

#define lambda 0.04
#define mu     0.03


int laplac_3d(pMesh mesh,pDeriv der) {
  printf(" Warning: laplace not implemented in 3D.\n");
  return(1);
}

int laplac_s(pMesh mesh,pDeriv der) {
  printf(" Warning: laplace not implemented for 3D surface meshes.\n");
  return(1);
}


int laplac_2d(pMesh mesh,pDeriv der) {
  pTria       pt;
  pPoint      p0,p1,p2;
  double     *hes,dd,d3;
  int         i,k,l,i0,i1,i2,it,maxtou;

  /* mem alloc */
  hes = (double*)M_malloc((3*mesh->nt+1)*sizeof(double),"laplac");
  assert(hes);

  /* average hessian */
  maxtou = 3;
  d3 = 1.0 / 3.0;
  for (it=0; it<maxtou; it++) {
    l = 1;

    /* set 0 to boundary */
    for (k=1; k<=mesh->nt; k++) {
      pt  = &mesh->tria[k];
      i0 = pt->v[0];
      i1 = pt->v[1];
      i2 = pt->v[2];
      p0 = &mesh->point[i0];
      p1 = &mesh->point[i1];
      p2 = &mesh->point[i2];

      hes[l++] = d3 * (der[i0].hes[0] + der[i1].hes[0] + der[i2].hes[0]);
      hes[l++] = d3 * (der[i0].hes[1] + der[i1].hes[1] + der[i2].hes[1]);
      hes[l++] = d3 * (der[i0].hes[2] + der[i1].hes[2] + der[i2].hes[2]);
      if ( p0->b )  memset(der[i0].hes,0,3*sizeof(double));
      if ( p1->b )  memset(der[i1].hes,0,3*sizeof(double));
      if ( p2->b )  memset(der[i2].hes,0,3*sizeof(double));
    }

    /* avg value to boundary */
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      for (i=0; i<3; i++) {
        i0 = pt->v[i];
        p0 = &mesh->point[i0];
        if ( !p0->b )  continue;
        dd  = pt->aire / p0->aire;
        dd *= d3 * (p0->nv-1);
        der[i0].hes[0] += dd * hes[3*(k-1)+1];
        der[i0].hes[1] += dd * hes[3*(k-1)+2];
        der[i0].hes[2] += dd * hes[3*(k-1)+3];
      }
    }
  }
  
  M_free(hes);
  return(1);
}


int lissag_3d(pMesh mesh,pSol sol,int is) {
  pPoint    p0,p1;
	pTetra    pt;
  double   *nu,u,v;
  int       list[LONMAX+2],ilist,i,k,iadr,nsdep,ip1,nb;

	nu = calloc(mesh->np+1,sizeof(double));
	assert(nu);

  /* 1st stage: new value */
  for (k=1; k<=mesh->np; k++) {
		p0    = &mesh->point[k];
    nsdep = p0->s;
    assert(nsdep);

	  pt = &mesh->tetra[nsdep];
    i  = 0;
    if ( pt->v[1] == k )      i = 1;
    else if ( pt->v[2] == k ) i = 2;
    else if ( pt->v[3] == k ) i = 3;

		ilist = boulep(mesh,nsdep,i,list);
    u  = getSol(sol,k,is);
		v  = 0.0;
		nb = 0;
	
	  /*mesh->info.ddebug = (k == 13149);
	  if ( mesh->info.ddebug )  printf("ATTENTION point %d,  INIVAL %f\n",k,u);*/
    for (i=1; i<=ilist; i++) {
      ip1 = list[i];
      p1  = &mesh->point[ip1];
      v  += getSol(sol,ip1,is);
			nb++;
    }
		v     = v / nb;
		nu[k] = u + lambda*(v-u);
    /*if ( mesh->info.ddebug)  printf("   MOY  %f\n",v);*/
  }

  /* 2nd stage: update value */
  for (k=1; k<=mesh->np; k++) {
		p0    = &mesh->point[k];
    nsdep = p0->s;
    assert(nsdep);

	  pt = &mesh->tetra[nsdep];
    i  = 0;
    if ( pt->v[1] == k )      i = 1;
    else if ( pt->v[2] == k ) i = 2;
    else if ( pt->v[3] == k ) i = 3;

		ilist = boulep(mesh,nsdep,i,list);
    u  = nu[k];
		v  = 0.0;
		nb = 0;
    for (i=1; i<=ilist; i++) {
      ip1 = list[i];
      p1  = &mesh->point[ip1];
      v  += nu[ip1];
			nb++;
    }
		v = v / nb;

    iadr = (k-1)*sol->size + is+1;
		sol->sol[iadr] = u - mu*(u-v);
	  /*mesh->info.ddebug = (k == 13149);
	  if ( mesh->info.ddebug )  printf("ATTENTION point %d,  FINVAL %f\n",k,sol->sol[iadr]);*/
  }

	free(nu);
	return(1);
}


/* regularisation of sol is */
int lissag_2d(pMesh mesh,pSol sol,int is) {
  pPoint    p0,p1;
  pTria     pt;
  double   *nu,u,v;
  int      *adja,adj,i,k,iadr,nsdep,ip1,nb;
  unsigned char voy,i1;

	nu = calloc(mesh->np+1,sizeof(double));
	assert(nu);

  /* 1st stage: new value */
  for (k=1; k<=mesh->np; k++) {
		p0    = &mesh->point[k];
    nsdep = p0->s;
    assert(nsdep);

    pt = &mesh->tria[nsdep];
    u  = getSol(sol,k,is);
		v  = 0.0;
    i  = 0;
    if ( pt->v[1] == k )      i = 1;
    else if ( pt->v[2] == k ) i = 2;
    iadr = 3*(nsdep-1)+1;
    adja = &mesh->adja[iadr];
    adj  = nsdep;
    voy  = i;
		nb   = 0;
	
	  /*mesh->info.ddebug = (k == 13149);
	  if ( mesh->info.ddebug )  printf("ATTENTION point %d,  INIVAL %f\n",k,u);*/
	
    do {
      pt  = &mesh->tria[adj];
      i1  = idir[voy+1];
      ip1 = pt->v[i1];
      p1  = &mesh->point[ip1];
      v  += getSol(sol,ip1,is);
			nb++;

      iadr = 3*(adj-1)+1;
      adja = &mesh->adja[iadr];
      adj  = adja[i1] / 3;
      voy  = adja[i1] % 3;
      voy  = idir[voy+1];
    }
    while ( adj && adj != nsdep );

    /* check open ball */
    if ( !adj ) {
			i = 0;
      if ( pt->v[1] == k )      i = 1;
      else if ( pt->v[2] == k ) i = 2;
      voy = idir[i+2];
      ip1 = pt->v[voy];
      p1  = &mesh->point[ip1];
      v   += getSol(sol,ip1,is);
			nb++;
    }
		v = v / nb;
    /*if ( mesh->info.ddebug)  printf("   MOY  %f\n",v);*/

		nu[k] = u + lambda*(v-u);
  }

  /* 2nd stage: update value */
  for (k=1; k<=mesh->np; k++) {
		p0    = &mesh->point[k];
    nsdep = p0->s;
    assert(nsdep);

    pt = &mesh->tria[nsdep];
    u  = nu[k];
		v  = 0.0;
    i  = 0;
    if ( pt->v[1] == k )      i = 1;
    else if ( pt->v[2] == k ) i = 2;
    iadr = 3*(nsdep-1)+1;
    adja = &mesh->adja[iadr];
    adj  = nsdep;
    voy  = i;
		nb   = 0;
    do {
      pt  = &mesh->tria[adj];
      i1  = idir[voy+1];
      ip1 = pt->v[i1];
      p1  = &mesh->point[ip1];
      v  += nu[ip1];
			nb++;
  
      iadr = 3*(adj-1)+1;
      adja = &mesh->adja[iadr];
      adj  = adja[i1] / 3;
      voy  = adja[i1] % 3;
      voy  = idir[voy+1];
    }
    while ( adj && adj != nsdep );

    /* check open ball */
    if ( !adj ) {
			i = 0;
      if ( pt->v[1] == k )      i = 1;
      else if ( pt->v[2] == k ) i = 2;
      voy = idir[i+2];
      ip1 = pt->v[voy];
      p1  = &mesh->point[ip1];
      v   += nu[ip1];
			nb++;
    }
		v = v / nb;

    iadr = (k-1)*sol->size + is+1;
		sol->sol[iadr] = u - mu*(u-v);
	  /*mesh->info.ddebug = (k == 13149);
	  if ( mesh->info.ddebug )  printf("ATTENTION point %d,  FINVAL %f\n",k,sol->sol[iadr]);*/
  }

	free(nu);
  return(1);
}

