#include "mshmet.h"


/* return unique sol is at vertex ip,
   convert vector-tensor into scalar */
double getSol_3d(pSol sol,int ip,int is) {
  double  ux,uy,uz,u,m[6],lambda[3],vp[3][3];
  int     i,iadr;

  u = 0.0;
  iadr = (ip-1)*sol->size + is+1;
  switch(sol->typtab[is]) {
  case GmfSca:
    u = sol->sol[iadr];
    break;
  case GmfVec:
    ux = sol->sol[iadr];
    uy = sol->sol[iadr+1];
    uz = sol->sol[iadr+2];
    u  = sqrt(ux*ux + uy*uy + uz*uz);
    break;
  case GmfSymMat:
    for (i=0; i<6; i++)  m[i] = sol->sol[iadr+i];
    eigenv(1,m,lambda,vp);
    u = MS_MAX(lambda[0],MS_MAX(lambda[1],lambda[2]));
  }

  return(u);
}


double getSol_2d(pSol sol,int ip,int is) {
  double  ux,uy,u,m[3],lambda[2],vp[2][2];
  int     i,iadr;

  u = 0.0;
  iadr = (ip-1)*sol->size + is+1;
  switch(sol->typtab[is]) {
  case GmfSca:
    u = sol->sol[iadr];
    break;
  case GmfVec:
    ux = sol->sol[iadr];
    uy = sol->sol[iadr+1];
    u  = sqrt(ux*ux + uy*uy);
    break;
  case GmfSymMat:
    for (i=0; i<3; i++)  m[i] = sol->sol[iadr+i];
    eigen2(m,lambda,vp);
    u = MS_MAX(lambda[0],lambda[1]);
  }

  return(u);
}


/* compute gradient via least square approx 
   in:  ip = vertex num   is = sol num 
   out: grd[3] gradient */
int gradLS_3d(pMesh mesh,pSol sol,int ip,int is,double *grd) {
  pPoint    p0,p1;
	pTetra    pt;
  double    a[6],b[3],ax,ay,az,du,u,u1,aa,bb,cc,dd,ee,ff;
  int       i,nsdep,ip1,list[LONMAX+2],ilist;

  p0     = &mesh->point[ip];
  p0->nv = 0;
  nsdep  = p0->s;
	assert(nsdep);

	pt = &mesh->tetra[nsdep];
	i  = 0;
  if ( pt->v[1] == ip )      i = 1;
  else if ( pt->v[2] == ip ) i = 2;
  else if ( pt->v[3] == ip ) i = 3;

  ilist = boulep(mesh,nsdep,i,list);
	if ( ilist < 1 )  return(0);

  memset(a,0,6*sizeof(double));
  memset(b,0,3*sizeof(double));
  u = getSol(sol,ip,is);

  /* Gradient: ui = u + <grad u,PPi> */
	for (i=1; i<=ilist; i++) {
		ip1 = list[i];
		p1  = &mesh->point[ip1];
    u1  = getSol(sol,ip1,is);
    du  = u1 - u;

    /* M = At*A symmetric definite positive */
    ax = p1->c[0] - p0->c[0];
    ay = p1->c[1] - p0->c[1];
    az = p1->c[2] - p0->c[2];
    a[0] += ax * ax;
    a[1] += ax * ay;
    a[2] += ax * az;
    a[3] += ay * ay;
    a[4] += ay * az;
    a[5] += az * az;

    /* b = At*du */
    b[0] += ax * du;
    b[1] += ay * du;
    b[2] += az * du;
    p0->nv++;
	}

  /* solution of A(3,3)*grad(1,3) = b(1,3) */ 
	aa = a[3]*a[5] - a[4]*a[4];
  bb = a[4]*a[2] - a[1]*a[5];
  cc = a[1]*a[4] - a[2]*a[3];
  du = aa*a[0] + bb*a[1] + cc*a[2];
  if ( fabs(du) == 0.0 ) {
    fprintf(stdout," Invalid matrix (%E). Exit\n",du);
    exit(1);
  }
  du = 1.0 / du;
  dd = a[0]*a[5] - a[2]*a[2];
  ee = a[1]*a[2] - a[0]*a[4];
  ff = a[0]*a[3] - a[1]*a[1];

  grd[0] = (aa*b[0] + bb*b[1] + cc*b[2]) * du;
  grd[1] = (bb*b[0] + dd*b[1] + ee*b[2]) * du;
  grd[2] = (cc*b[0] + ee*b[1] + ff*b[2]) * du;

  /* normalize */
  if ( mesh->info.ls ) {
    du = grd[0]*grd[0] + grd[1]*grd[1] + grd[2]*grd[2];
    if ( du > EPS1 ) {
      du      = 1.0 / sqrt(du);
      grd[0] *= du;
      grd[1] *= du;
      grd[2] *= du;
    }
  }
  return(1);
}


int gradLS_2d(pMesh mesh,pSol sol,int ip,int is,double *grd) {
  pPoint    p0,p1;
  pTria     pt;
  double    a[3],b[2],ax,ay,du,u,u1;
  int      *adja,adj,i,iadr,nsdep,ip1;
  unsigned char voy,i1;

  p0      = &mesh->point[ip];
  p0->nv  = 0;
  nsdep   = p0->s;
  if ( !nsdep ) {
    fprintf(stdout," No simplex stored. Exit\n");
    exit(EXIT_FAILURE);
  }
  pt   = &mesh->tria[nsdep];
  iadr = 3*(nsdep-1)+1;
  adja = &mesh->adja[iadr];

  memset(a,0,3*sizeof(double));
  memset(b,0,2*sizeof(double));
  u = getSol(sol,ip,is);

  if ( pt->v[0] == ip )      i = 0;
  else if ( pt->v[1] == ip ) i = 1;
  else                       i = 2;

  /* Gradient: ui = u + <grad u,PPi> */
  adj = nsdep;
  voy = i;
  do {
    pt   = &mesh->tria[adj];
    i1   = idir[voy+1];
    ip1  = pt->v[i1];
    p1   = &mesh->point[ip1];
    u1   = getSol(sol,ip1,is);

    /* M = At*A symmetric definite positive */
    ax = p1->c[0] - p0->c[0];
    ay = p1->c[1] - p0->c[1];
    a[0] += ax * ax;
    a[1] += ax * ay;
    a[2] += ay * ay;

    /* b = At*du */
    du    = u1 - u;
    b[0] += ax * du;
    b[1] += ay * du;
    p0->nv++;

    iadr = 3*(adj-1)+1;
    adja = &mesh->adja[iadr];
    adj  = adja[i1] / 3;
    voy  = adja[i1] % 3;
    voy  = idir[voy+1];
  }
  while ( adj && adj != nsdep );

  /* check open ball */
  if ( !adj ) {
    if ( pt->v[0] == ip )       i = 0;
    else if ( pt->v[1] == ip )  i = 1;
    else                        i = 2;
    voy   = idir[i+2];
    ip1   = pt->v[voy];
    p1    = &mesh->point[ip1];
    u1    = getSol(sol,ip1,is);

    ax = p1->c[0] - p0->c[0];
    ay = p1->c[1] - p0->c[1];
    a[0] += ax * ax;
    a[1] += ax * ay;
    a[2] += ay * ay;

    /* b = At*du */
    du    = u1 - u;
    b[0] += ax * du;
    b[1] += ay * du;
    p0->nv++;
    p0->b = 1;
  }

  /* solution of A(2,2)*grad(1,2) = b(1,2) */ 
  du = a[0]*a[2] - a[1]*a[1];
  if ( fabs(du) == 0.0 ) {
    fprintf(stdout," Invalid matrix (%E). Exit\n",du);
    exit(1);
  }
  du  = 1.0 / du;
  grd[0] = (a[2]*b[0] - a[1]*b[1]) * du;
  grd[1] = (a[0]*b[1] - a[1]*b[0]) * du;

  /* normalize */
  if ( mesh->info.ls ) {
    du = grd[0]*grd[0] + grd[1]*grd[1];
    if ( du > EPS1 ) {
      du      = 1.0 / sqrt(du);
      grd[0] *= du;
      grd[1] *= du;
    }
    else {
      grd[0] = 1;
      grd[1] = 0;
    }
  }

  return(1);
}

int gradLS_s(pMesh mesh,pSol sol,int ip,int is,double *grd) {
  pPoint    p0,p1;
  pTria     pt;
  double    a[6],b[3],ax,ay,az,du,u,u1,aa,bb,cc,dd,ee,ff;
  int      *adja,adj,i,iadr,nsdep,ip1;
  unsigned char voy,i1;

  p0      = &mesh->point[ip];
  p0->nv  = 0;
  nsdep   = p0->s;
  if ( !nsdep ) {
    // Algiane: Try to ignore those point do not work for now
    fprintf(stdout," Warning: No simplex stored at point %d.\n",ip);
    grd[0] = 1;
    grd[1] = 0;
    grd[2] = 0;
    return -1;
  }
  pt   = &mesh->tria[nsdep];
  iadr = 3*(nsdep-1)+1;
  adja = &mesh->adja[iadr];

  memset(a,0,6*sizeof(double));
  memset(b,0,3*sizeof(double));
  u = getSol(sol,ip,is);

  if ( pt->v[0] == ip )      i = 0;
  else if ( pt->v[1] == ip ) i = 1;
  else                       i = 2;

  /* Gradient: ui = u + <grad u,PPi> */
  adj = nsdep;
  voy = i;
  do {
    pt   = &mesh->tria[adj];
    i1   = idir[voy+1];
    ip1  = pt->v[i1];
    p1   = &mesh->point[ip1];

    if ( p1->s ) {
      u1   = getSol(sol,ip1,is);
      du  = u1 - u;

      /* M = At*A symmetric definite positive */
      ax = p1->c[0] - p0->c[0];
      ay = p1->c[1] - p0->c[1];
      az = p1->c[2] - p0->c[2];
      a[0] += ax * ax;
      a[1] += ax * ay;
      a[2] += ax * az;
      a[3] += ay * ay;
      a[4] += ay * az;
      a[5] += az * az;

      /* b = At*du */
      b[0] += ax * du;
      b[1] += ay * du;
      b[2] += az * du;
      p0->nv++;
    }

    iadr = 3*(adj-1)+1;
    adja = &mesh->adja[iadr];
    adj  = adja[i1] / 3;
    voy  = adja[i1] % 3;
    voy  = idir[voy+1];
  }
  while ( adj && adj != nsdep );

  /* check open ball */
  if ( !adj ) {
    if ( pt->v[0] == ip )       i = 0;
    else if ( pt->v[1] == ip )  i = 1;
    else                        i = 2;
    voy   = idir[i+2];
    ip1   = pt->v[voy];
    p1    = &mesh->point[ip1];

    if ( p1->s ) {
      u1    = getSol(sol,ip1,is);
      du  = u1 - u;

      ax = p1->c[0] - p0->c[0];
      ay = p1->c[1] - p0->c[1];
      az = p1->c[2] - p0->c[2];
      a[0] += ax * ax;
      a[1] += ax * ay;
      a[2] += ax * az;
      a[3] += ay * ay;
      a[4] += ay * az;
      a[5] += az * az;

      /* b = At*du */
      b[0] += ax * du;
      b[1] += ay * du;
      b[2] += az * du;
      p0->nv++;
      p0->b = 1;
    }
  }

  /* solution of A(3,3)*grad(1,3) = b(1,3) */
  aa = a[3]*a[5] - a[4]*a[4];
  bb = a[4]*a[2] - a[1]*a[5];
  cc = a[1]*a[4] - a[2]*a[3];
  du = aa*a[0] + bb*a[1] + cc*a[2];
  if ( fabs(du) < EPS1 ) {
    /* Points of the ball are in the same plane */
    /* solution of A(2,2)*grad(1,2) = b(1,2) */ 
    du = a[0]*a[3] - a[1]*a[1];
    if ( fabs(du) < EPS1 ) {
      fprintf(stdout," Surface 2d Invalid matrix (%E). Exit\n",du);
      return(0);
    }
    du  = 1.0 / du;
    grd[0] = (a[3]*b[0] - a[1]*b[1]) * du;
    grd[1] = (a[0]*b[1] - a[1]*b[0]) * du;
    grd[2] = 0;

    return 1;
  }
  du = 1.0 / du;
  dd = a[0]*a[5] - a[2]*a[2];
  ee = a[1]*a[2] - a[0]*a[4];
  ff = a[0]*a[3] - a[1]*a[1];

  grd[0] = (aa*b[0] + bb*b[1] + cc*b[2]) * du;
  grd[1] = (bb*b[0] + dd*b[1] + ee*b[2]) * du;
  grd[2] = (cc*b[0] + ee*b[1] + ff*b[2]) * du;

  /* normalize */
  if ( mesh->info.ls ) {
    du = grd[0]*grd[0] + grd[1]*grd[1] + grd[2]*grd[2];
    if ( du > EPS1 ) {
      du      = 1.0 / sqrt(du);
      grd[0] *= du;
      grd[1] *= du;
      grd[2] *= du;
    }
    else {
      grd[0] = 1;
      grd[1] = 0;
      grd[2] = 0;
    }
  }

  return(1);
}
