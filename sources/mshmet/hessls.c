#include "mshmet.h"


static int gauss(int n,double m[n][n],double *x,double *b,char dbg) {
  double    nn,dd;
  int       i,j,k,ip;

  nn = m[0][0];
  for (i=0; i<n; i++) 
    for (j=0; j<n; j++)
      nn = MS_MAX(nn,m[i][j]);

  if ( fabs(nn) < EPS1 ) {
    fprintf(stdout,"  %%%% Null matrix\n");
    return(0);
  }
  nn = 1.0 / nn;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++)
      m[i][j] *= nn;
    b[i] *= nn; 
  }

  /* partial pivoting, column j */
  for (j=0; j<n-1; j++) {
    /* find line pivot */
    ip = j;
    for (i=j+1; i<n; i++)
      if ( fabs(m[i][j]) > fabs(m[ip][j]) )  ip = i;

    /* swap i <-> ip */
    if ( i != ip ) {
      for (k=j; k<n; k++) {
        dd       = m[j][k];
        m[j][k]  = m[ip][k];
        m[ip][k] = dd;
      }
      dd    = b[j];
      b[j]  = b[ip];
      b[ip] = dd;
    }
    if ( fabs(m[j][j]) < EPS1 ) {
      if ( dbg )  fprintf(stdout,"  %%%% Null pivot: %e.\n",m[j][i]);
      return(0);
    }

    /* elimination */
    for (i=j+1; i<n; i++) {
      dd = m[i][j] / m[j][j];
      m[i][j] = 0.0;
      for (k=j+1; k<n; k++)
        m[i][k] -= dd * m[j][k]; 
      b[i] -= dd * b[j];
    }
  }

  if ( fabs(m[n-1][n-1]) < EPS1 ) {
    if ( dbg )  fprintf(stdout,"  %%%% Null pivot.\n");
    return(0);
  }

  x[n-1] = b[n-1] / m[n-1][n-1];
  for (i=n-2; i>=0; i--) {
    dd = 0.0;
    for (j=i+1; j<n; j++)
      dd += m[i][j] * x[j];
    x[i] = (b[i] - dd) / m[i][i];
  }

  if ( dbg ) {
    for (i=0; i<n; i++) {
      dd = 0.;
      for (j=0; j<n; j++)
        dd += m[i][j] * x[j];
      if ( fabs(dd-b[i]) > EPS ) {
        fprintf(stdout,"  Ax[%d] = %f   b[%d] = %f\n",i,dd,i,b[i]);
        exit(1);
      }
    }
  }

  return(1);
}


int hessLS_3d(pMesh mesh,pSol sol,int ip,int is,double *grd,double *hes) {
  pPoint    p0,p1;
  pTetra    pt;
  double    a[6],b[6],m[6][6],ax,ay,az;
  double    du,u,u1;
  int       i,j,l,ier,nsdep,ip1,list[LONMAX+2],ilist;

  /*mesh->info.ddebug = (ip == 117);
  if ( mesh->info.ddebug )  printf("ATTENTION point %d\n",ip);*/

  p0 = &mesh->point[ip];
  if ( p0->nv < 6 )  return(-1);
  nsdep = p0->s;
  assert(nsdep);

  pt = &mesh->tetra[nsdep];
  i  = 0;
  if ( pt->v[1] == ip )      i = 1;
  else if ( pt->v[2] == ip ) i = 2;
  else if ( pt->v[3] == ip ) i = 3;

  ilist = boulep(mesh,nsdep,i,list);
  if ( ilist < 1 )  return(0);
	
  memset(b,0,6*sizeof(double));
  memset(m,0,36*sizeof(double));
  u = getSol(sol,ip,is);

  /* Hessian: Ui = U + <gradU,PPi> + 0.5*(tPPi.Hess.PPi) */
  for (i=1; i<=ilist; i++) {
    ip1 = list[i];
    p1  = &mesh->point[ip1];
    u1  = getSol(sol,ip1,is);
	
    ax = p1->c[0] - p0->c[0];
    ay = p1->c[1] - p0->c[1];
    az = p1->c[2] - p0->c[2];
 
    a[0] = 0.5 * ax*ax;
    a[1] =       ax*ay;
    a[2] =       ax*az;
    a[3] = 0.5 * ay*ay;
    a[4] =       ay*az;
    a[5] = 0.5 * az*az;
    du   = (u1-u) - (ax*grd[0] + ay*grd[1] + az*grd[2]);

    /*if ( mesh->info.ddebug ) {
      printf("  %d  u1 %f  a: %f %f %f\n",ip1,u1,ax,ay,az);
      printf("  A: %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f\n",a[0],a[1],a[2],a[3],a[4],a[5]);
    }*/

    /* M = At*A symmetric positive definite */
    for (j=0; j<6; j++) {
      for (l=j; l<6; l++) {
        m[j][l] += a[j]*a[l];
        m[l][j] += a[l]*a[j];
      }
    }

    /* c = At*b */
    for (j=0; j<6; j++)
      b[j] += a[j]*du;
  }

  /* solve m(6,6)*x(6) = b(6) */
  ier = gauss(6,m,hes,b,mesh->info.ddebug);
  if ( !ier && mesh->info.ddebug ) {
    fprintf(stdout," Ill cond'ed matrix (%d, %d).\n",ip,ier);
    return(-1);
  }
  /*if ( mesh->info.ddebug )  
    printf("hessien: %f %f %f  kappa %f\n",hes[0],hes[1],hes[2],fabs(hes[0])+fabs(hes[2]));*/

  p0->h = 1;
  return(1);
}


/* least square approximation of hessian */
int hessLS_2d(pMesh mesh,pSol sol,int ip,int is,double *grd,double *hes) {
  pPoint    p0,p1;
  pTria     pt;
  double    a[3],b[2],ma[6],mb[3],ax,ay,wgt;
  double    det,aa,bb,cc,dd,ee,ff;
  double    u,u1;
  int      *adja,adj,i,iadr,nsdep,ip1;
  unsigned char voy,i1;

  /*mesh->info.ddebug = (ip == 13149);
  if ( mesh->info.ddebug )  printf("ATTENTION point %d\n",ip);*/

  p0    = &mesh->point[ip];
  nsdep = p0->s;
  if ( p0->nv < 3 )  return(-1);
  assert(nsdep);
  pt   = &mesh->tria[nsdep];
  iadr = 3*(nsdep-1)+1;
  adja = &mesh->adja[iadr];

  memset(ma,0,6*sizeof(double));
  memset(mb,0,3*sizeof(double));
  memset(a,0,3*sizeof(double));
  memset(b,0,2*sizeof(double));
  u = getSol(sol,ip,is);

  /*if ( mesh->info.ddebug ) 
    printf("point %d: tria %d  sol %f  grad %f %f\n",ip,nsdep,u,grd[0],grd[1]);*/

  i = 0;
  if ( pt->v[1] == ip )      i = 1;
  else if ( pt->v[2] == ip ) i = 2;

  /* Hessian: Ui = U + <gradU,PPi> + 0.5*(tPPi.Hess.PPi) */
  adj = nsdep;
  voy = i;
  wgt = 1.0;
  do {
    pt    = &mesh->tria[adj];
    i1    = idir[voy+1];
    ip1   = pt->v[i1];
    p1    = &mesh->point[ip1];
    u1    = getSol(sol,ip1,is);
    /*if ( mesh->info.ddebug) 
      printf("  tr: %d : %d %d %d   sol %f / %d\n",adj,pt->v[0],pt->v[1],pt->v[2],u1,ip1);*/

    ax    = p1->c[0] - p0->c[0];
    ay    = p1->c[1] - p0->c[1];
    a[0]  = 0.5 * ax*ax;
    a[1]  =       ax*ay;
    a[2]  = 0.5 * ay*ay;
    dd    = (u1-u) - (ax*grd[0] + ay*grd[1]);

    /* M = At*A symmetric definite positive */
    ma[0] += a[0]*a[0] * wgt;
    ma[1] += a[0]*a[1] * wgt;
    ma[2] += a[1]*a[1] * wgt;
    ma[3] += a[0]*a[2] * wgt;
    ma[4] += a[1]*a[2] * wgt;
    ma[5] += a[2]*a[2] * wgt;

    /* c = At*b */
    mb[0] += a[0]*dd * wgt;
    mb[1] += a[1]*dd * wgt;
    mb[2] += a[2]*dd * wgt;

    iadr = 3*(adj-1)+1;
    adja = &mesh->adja[iadr];
    adj  = adja[i1] / 3;
    voy  = adja[i1] % 3;
    voy  = idir[voy+1];
  }
  while ( adj && adj != nsdep );

  /* check open ball */
  if ( !adj ) {
    if ( pt->v[0] == ip )      i = 0;
    else if ( pt->v[1] == ip ) i = 1;
    else                       i = 2;
    voy   = idir[i+2];
    ip1   = pt->v[voy];
    p1    = &mesh->point[ip1];
    u1    = getSol(sol,ip1,is);

    ax    = p1->c[0] - p0->c[0];
    ay    = p1->c[1] - p0->c[1];
    a[0]  = 0.5 * ax*ax;
    a[1]  =       ax*ay;
    a[2]  = 0.5 * ay*ay;
    dd    = (u1-u) - (ax*grd[0] + ay*grd[1]);

    /* M = At*A symmetric definite positive */
    ma[0] += a[0]*a[0] * wgt;
    ma[1] += a[0]*a[1] * wgt;
    ma[2] += a[1]*a[1] * wgt;
    ma[3] += a[0]*a[2] * wgt;
    ma[4] += a[1]*a[2] * wgt;
    ma[5] += a[2]*a[2] * wgt;

    /* c = At*b */
    mb[0] += a[0]*dd * wgt;
    mb[1] += a[1]*dd * wgt;
    mb[2] += a[2]*dd * wgt;
  }

  if ( fabs(ma[0]) == 0.0 ) {
    if ( mesh->info.ddebug )  fprintf(stdout," Ill cond'ed matrix (%d, %E).\n",ip,ma[0]);
    exit(1);
  }

  /* direct solving */
  aa = ma[2]*ma[5] - ma[4]*ma[4];
  bb = ma[4]*ma[3] - ma[1]*ma[5];
  cc = ma[1]*ma[4] - ma[2]*ma[3];
  dd = ma[0]*aa + ma[1]*bb + ma[3]*cc;

  /* singular matrix */
  if ( fabs(dd) == 0.0 ) {
    /*if ( mesh->info.ddebug )  fprintf(stderr," Singular matrix (%d, %E).\n",ip,dd);*/
    return(-1);
  }
  det = 1.0 / dd;
  dd = ma[0]*ma[5] - ma[3]*ma[3];
  ee = ma[1]*ma[3] - ma[0]*ma[4];
  ff = ma[0]*ma[2] - ma[1]*ma[1];

  hes[0] = (mb[0]*aa + mb[1]*bb + mb[2]*cc) * det;
  hes[1] = (mb[0]*bb + mb[1]*dd + mb[2]*ee) * det;
  hes[2] = (mb[0]*cc + mb[1]*ee + mb[2]*ff) * det;

  /*if ( mesh->info.ddebug )  
    printf("hessien: %f %f %f  kappa %f\n",hes[0],hes[1],hes[2],fabs(hes[0])+fabs(hes[2]));*/
  
  if ( !p0->b && p0->nv == 4 )  return(-1);
  p0->h = 1;

  return(1);
}

/* least square approximation of hessian */
int hessLS_s(pMesh mesh,pSol sol,int ip,int is,double *grd,double *hes) {
  pPoint    p0,p1;
  pTria     pt;
  double    a[6],b[6],m[6][6],ax,ay,az,wgt;
  double    det;
  double    du,u,u1;
  int      *adja,adj,i,j,l,iadr,nsdep,ip1;
  unsigned char voy,i1;

  p0    = &mesh->point[ip];
  nsdep = p0->s;
  if ( p0->nv < 3 ) {
    /* Undetermined boundaries */
    return(-1);
  }

  if ( !nsdep ) {
    p0->h = 1;
    // Algiane: Try to ignore those point do not work for now
    return -1;
  }

  pt   = &mesh->tria[nsdep];
  iadr = 3*(nsdep-1)+1;
  adja = &mesh->adja[iadr];

  memset(b,0,6*sizeof(double));
  memset(m,0,36*sizeof(double));
  u = getSol(sol,ip,is);

  i = 0;
  if ( pt->v[1] == ip )      i = 1;
  else if ( pt->v[2] == ip ) i = 2;

  /* Hessian: Ui = U + <gradU,PPi> + 0.5*(tPPi.Hess.PPi) */
  adj = nsdep;
  voy = i;
  do {
    pt    = &mesh->tria[adj];
    i1    = idir[voy+1];
    ip1   = pt->v[i1];
    p1    = &mesh->point[ip1];
    if ( p1->s ) {
      u1    = getSol(sol,ip1,is);
      /*if ( mesh->info.ddebug)
        printf("  tr: %d : %d %d %d   sol %f / %d\n",adj,pt->v[0],pt->v[1],pt->v[2],u1,ip1);*/

      ax    = p1->c[0] - p0->c[0];
      ay    = p1->c[1] - p0->c[1];
      az    = p1->c[2] - p0->c[2];

      a[0] = 0.5 * ax*ax;
      a[1] =       ax*ay;
      a[2] =       ax*az;
      a[3] = 0.5 * ay*ay;
      a[4] =       ay*az;
      a[5] = 0.5 * az*az;
      du   = (u1-u) - (ax*grd[0] + ay*grd[1] + az*grd[2]);

      /* M = At*A symmetric positive definite */
      for (j=0; j<6; j++) {
        m[j][j] += a[j]*a[j];
        for (l=j+1; l<6; l++) {
          m[j][l] += a[j]*a[l];
          m[l][j] += a[l]*a[j];
        }
      }

      /* c = At*b */
      for (j=0; j<6; j++)
        b[j] += a[j]*du;
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
    if ( pt->v[0] == ip )      i = 0;
    else if ( pt->v[1] == ip ) i = 1;
    else                       i = 2;
    voy   = idir[i+2];
    ip1   = pt->v[voy];
    p1    = &mesh->point[ip1];
    if ( p1->s ) {
      u1    = getSol(sol,ip1,is);

      ax    = p1->c[0] - p0->c[0];
      ay    = p1->c[1] - p0->c[1];
      az    = p1->c[2] - p0->c[2];

      a[0] = 0.5 * ax*ax;
      a[1] =       ax*ay;
      a[2] =       ax*az;
      a[3] = 0.5 * ay*ay;
      a[4] =       ay*az;
      a[5] = 0.5 * az*az;
      du   = (u1-u) - (ax*grd[0] + ay*grd[1] + az*grd[2]);

      /* M = At*A symmetric positive definite */
      for (j=0; j<6; j++) {
        m[j][j] += a[j]*a[j];
        for (l=j+1; l<6; l++) {
          m[j][l] += a[j]*a[l];
          m[l][j] += a[l]*a[j];
        }
      }

      /* c = At*b */
      for (j=0; j<6; j++)
        b[j] += a[j]*du;
    }
  }

  /* solve m(6,6)*x(6) = b(6) */
  int ier = gauss(6,m,hes,b,mesh->info.ddebug);

  if ( !ier ) {
    /* Points of the ball lie probably along the same plane */
    assert( fabs(grd[2]) < EPS1 );

    double ma[6];
    ma[0] = m[0][0];
    ma[1] = m[0][1];
    ma[2] = m[1][1];
    ma[3] = m[0][3];
    ma[4] = m[1][3];
    ma[5] = m[3][3];


    if ( fabs(ma[0]) < EPS1 ) {
      fprintf(stdout,"Surface 2d Ill cond'ed matrix (%d, %E).\n",ip,ma[0]);
      return(-1);
    }
        /* direct solving */
    double aa = ma[2]*ma[5] - ma[4]*ma[4];
    double bb = ma[4]*ma[3] - ma[1]*ma[5];
    double cc = ma[1]*ma[4] - ma[2]*ma[3];
    double dd = ma[0]*aa + ma[1]*bb + ma[3]*cc;

    /* singular matrix */
    if ( fabs(dd) < EPS1 ) {
      /* Constant solution: put hmax */
      hes[0] = hes[3] = hes[5] = 0;
      hes[1] = hes[4] = hes[2] = 0;
      p0->h = 1;
    }
    else {
      det = 1.0 / dd;
      double dd = ma[0]*ma[5] - ma[3]*ma[3];
      double ee = ma[1]*ma[3] - ma[0]*ma[4];
      double ff = ma[0]*ma[2] - ma[1]*ma[1];
      hes[0] = (b[0]*aa + b[1]*bb + b[3]*cc) * det;
      hes[1] = (b[0]*bb + b[1]*dd + b[3]*ee) * det;
      hes[3] = (b[0]*cc + b[1]*ee + b[3]*ff) * det;
      hes[2] = hes[4] = 0;
      /* Put 0 along z-dir to impose hmax at the end */
      hes[5] = 0;
    }
  }

  if ( !ier && mesh->info.ddebug ) {
    fprintf(stdout," Ill cond'ed matrix (%d, %d).\n",ip,ier);
    return(-1);
  }

  if ( !p0->b && p0->nv == 4 )  return(-1);
  p0->h = 1;

  return(1);
}


/* average values of neighbors */
int avgval_3d(pMesh mesh,pDeriv der,int ip,double *hes) {
  pPoint    p0,p1;
  pTetra    pt;
  double    dd;
  int       list[LONMAX+2],ilist,i,j,nsdep,ip1,nb;

  p0    = &mesh->point[ip];
  nsdep = p0->s;
  assert(nsdep);

  pt = &mesh->tetra[nsdep];
  i  = 0;
  if ( pt->v[1] == ip )      i = 1;
  else if ( pt->v[2] == ip ) i = 2;
  else if ( pt->v[3] == ip ) i = 3;

  ilist = boulep(mesh,nsdep,i,list);
  if ( ilist < 1 )  return(0);

  memset(hes,0,6*sizeof(double));
  nb = 0;
  for (i=1; i<=ilist; i++) {
    ip1 = list[i];
    p1  = &mesh->point[ip1];
    if ( p1->h > 0 ) {
      for (j=0; j<6; j++)
        hes[j] += der[ip1].hes[j];
      nb++;
    }
  }

  if ( nb > 0 ) {
    dd = 1.0 / (double)nb;
    for (j=0; j<6; j++)
      hes[j] = hes[j] * dd;
    p0->h = 1;
    return(1);
  }

  return(0);
}


int avgval_2d(pMesh mesh,pDeriv der,int ip,double *hes) {
  pPoint    p0,p1;
  pTria     pt;
  double    hxx,hxy,hyy,dd;
  int      *adja,adj,i,iadr,nsdep,ip1,nb;
  unsigned char voy,i1;

  p0    = &mesh->point[ip];
  nsdep = p0->s;
  assert(nsdep);
  pt   = &mesh->tria[nsdep];
  iadr = 3*(nsdep-1)+1;
  adja = &mesh->adja[iadr];

  if ( pt->v[0] == ip )      i = 0;
  else if ( pt->v[1] == ip ) i = 1;
  else                       i = 2;

  adj = nsdep;
  voy = i;
  hxx = 0.0;
  hxy = 0.0;
  hyy = 0.0;
  nb  = 0;
  do {
    pt    = &mesh->tria[adj];
    i1    = idir[voy+1];
    ip1   = pt->v[i1];
    p1    = &mesh->point[ip1];
    if ( p1->h > 0 ) {
      hxx += der[ip1].hes[0];
      hxy += der[ip1].hes[1];
      hyy += der[ip1].hes[2];
      nb++;
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
    if ( p1->h > 0 ) {
      hxx += der[ip1].hes[0];
      hxy += der[ip1].hes[1];
      hyy += der[ip1].hes[2];
      nb++;
    }
  }

  if ( nb > 0 ) {
    dd = 1.0 / (double)nb;
    hes[0] = hxx * dd;
    hes[1] = hxy * dd;
    hes[2] = hyy * dd;
    p0->h = 1;
    return(1);
  }

  return(0);
}

int avgval_s(pMesh mesh,pDeriv der,int ip,double *hes) {
  pPoint    p0,p1;
  pTria     pt;
  double    dd;
  int      *adja,adj,i,j,iadr,nsdep,ip1,nb;
  unsigned char voy,i1;

  p0    = &mesh->point[ip];
  nsdep = p0->s;
  if ( !nsdep ) {
    p0->h = 1;
    // Algiane: Try to ignore those point do not work for now
    return -1;
  }

  pt   = &mesh->tria[nsdep];
  iadr = 3*(nsdep-1)+1;
  adja = &mesh->adja[iadr];

  if ( pt->v[0] == ip )      i = 0;
  else if ( pt->v[1] == ip ) i = 1;
  else                       i = 2;

  adj = nsdep;
  voy = i;
  memset(hes,0,6*sizeof(double));
  nb  = 0;
  do {
    pt    = &mesh->tria[adj];
    i1    = idir[voy+1];
    ip1   = pt->v[i1];
    p1    = &mesh->point[ip1];
    if ( p1->s && p1->h > 0 ) {
      for (j=0; j<6; j++)
        hes[j] += der[ip1].hes[j];
      nb++;
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
    if ( p1->s && p1->h > 0 ) {
      for (j=0; j<6; j++)
        hes[j] += der[ip1].hes[j];
      nb++;
    }
  }

  if ( nb > 0 ) {
    dd = 1.0 / (double)nb;
    for (j=0; j<6; j++)
      hes[j] = hes[j] * dd;
    p0->h = 1;
    return(1);
  }

  return(0);
}


/* assign value of closest vertex */
int clsval_3d(pMesh mesh,pDeriv der,int ip,double *hes) {
  pPoint    p0,p1;
  pTetra    pt;
  int       i,j,nsdep,ip1;

  p0    = &mesh->point[ip];
  nsdep = p0->s;
  assert(nsdep);
  pt   = &mesh->tetra[nsdep];

  pt = &mesh->tetra[nsdep];
  i  = 1;
  if ( pt->v[1] == ip )      i = 2;
  else if ( pt->v[2] == ip ) i = 3;
  else if ( pt->v[3] == ip ) i = 0;

  ip1 = pt->v[i];
  p1  = &mesh->point[ip1];
  if ( !p1->h )  return(0);
  for (j=0; j<6; j++)
    hes[j] = der[ip1].hes[j];

  p0->h = 1;
  return(1);
}


int clsval_2d(pMesh mesh,pDeriv der,int ip,double *hes) {
  pPoint    p0,p1;
  pTria     pt;
  int       nsdep,ip1;
  unsigned char i;

  p0    = &mesh->point[ip];
  nsdep = p0->s;
  if ( !nsdep ) {
    fprintf(stdout," No simplex stored. Exit\n");
    exit(EXIT_FAILURE);
  }
  pt = &mesh->tria[nsdep];
  i  = 1;
  if ( pt->v[1] == ip )      i = 2;
  else if ( pt->v[2] == ip ) i = 0;

  ip1 = pt->v[i];
  p1  = &mesh->point[ip1];
  if ( !p1->h )  return(0);
  hes[0] = der[ip1].hes[0];
  hes[1] = der[ip1].hes[1];
  hes[2] = der[ip1].hes[2];

  p0->h = 1;
  return(1);
}

int clsval_s(pMesh mesh,pDeriv der,int ip,double *hes) {
  pPoint    p0,p1;
  pTria     pt;
  int       nsdep,ip1;
  unsigned char i,j;

  p0    = &mesh->point[ip];
  nsdep = p0->s;
  if ( !nsdep ) {
    fprintf(stdout," Warning: No simplex stored at point %d.\n",ip);
    for (j=0; j<6; j++)
      hes[j] = 0;
    // Algiane: Try to ignore those point do not work for now
    return -1;
  }
  pt = &mesh->tria[nsdep];
  i  = 1;
  if ( pt->v[1] == ip )      i = 2;
  else if ( pt->v[2] == ip ) i = 0;

  ip1 = pt->v[i];
  p1  = &mesh->point[ip1];
  if ( !p1->h )  return(0);
  for (j=0; j<6; j++)
    hes[j] = der[ip1].hes[j];

  p0->h = 1;
  return(1);
}


int nrmhes_3d(pMesh mesh,pSol sol,pDeriv der,int is) {
  pPoint   p0;
  Info     info;
  double   err,err1,errs,u,norm;
  double   maxu,maxg;
  int      i,k;

  info = mesh->info;

  switch(info.nnu) {  
  /* no norm */
  case 0:
    err = CTE3D / info.eps;
    for (k=1; k<=mesh->np; k++) {
      for (i=0; i<6; i++)  der[k].hes[i] *= err;
    }
    break;

  /* relative value: M(u)= |H(u)| / (err*||u||_inf) */
  case 1:
    maxu = 0.0;
    for (k=1; k<=mesh->np; k++) {
      u = fabs(getSol(sol,k,is));
      if ( u > maxu )  maxu = u;
    }
    if ( fabs(maxu) == 0.0 )  return(1);
    maxu =  CTE3D / (info.eps * maxu);
    for (k=1; k<=mesh->np; k++) {
      for (i=0; i<6; i++)  der[k].hes[i] *= maxu;
    }
    break;

  /* local norm: M(u)= |H(u)| / err*|u| */
  case 2:
    maxu = 0.0;
    for (k=1; k<=mesh->np; k++) {
      u = fabs(getSol(sol,k,is));
      if ( u > maxu )  maxu = u;
    }
    if ( maxu == 0.0 )
      errs = 0.01;
    else
      errs = maxu * 0.01;
    for (k=1; k<=mesh->np; k++) {
      u     = fabs(getSol(sol,k,is));
      err1  = MS_MAX(errs,u);
      maxu  = CTE3D / err1;
      for (i=0; i<6; i++)  der[k].hes[i] *= maxu;
    }
    break;

  /* local norm: M(u)= |H(u)| / err*(e1*|u|+e2*h*|du|) */
  case 3:
    maxu = 0.0;
    maxg = 0.0;
    for (k=1; k<=mesh->np; k++) {
      u    = fabs(getSol(sol,k,is));
      maxu = MS_MAX(maxu,u);
      norm = der[k].grd[0]*der[k].grd[0] + der[k].grd[1]*der[k].grd[1] + der[k].grd[2]*der[k].grd[2];
      norm = sqrt(norm);
      if ( norm > maxg )  maxg = norm;
    }
    if ( maxu == 0.0 )
      errs = 0.01;
    else
      errs = maxu * 0.01;
    err1 = maxg * 0.01;

    for (k=1; k<=mesh->np; k++) {
      p0    = &mesh->point[k];
      u     = getSol(sol,k,is);
      norm  = der[k].grd[0]*der[k].grd[0] + der[k].grd[1]*der[k].grd[1] + der[k].grd[2]*der[k].grd[2];
      norm  = sqrt(norm);
      norm  = MS_MAX(err1,norm);
      norm *= p0->rins / p0->nv;
      err1  = MS_MAX(errs,fabs(u));
      
      /* variable weight */
      err1  = 0.01*err1 + 0.01*norm;
      maxu  = CTE3D / (info.eps*err1);
      for (i=0; i<6; i++)  der[k].hes[i] *= maxu;
    } 
    break;
  }

  return(1);
}


int nrmhes_2d(pMesh mesh,pSol sol,pDeriv der,int is) {
  pPoint   p0;
  Info     info;
  double   err,err1,errs,u,norm;
  double   maxu,maxg;
  int      k;

  info = mesh->info;

  if ( info.nnu > 0 ) {
    for (k=1; k<=mesh->np; k++) {
      u = fabs(getSol(sol,k,is));
      sol->umax = MS_MAX(u,sol->umax);
      sol->umin = MS_MIN(u,sol->umin);
    }
	}
	
  switch(info.nnu) {  
  /* no norm */
  case 0:
    err = CTE2D / info.eps;
    for (k=1; k<=mesh->np; k++) {
      der[k].hes[0] *= err;
      der[k].hes[1] *= err;
      der[k].hes[2] *= err;
    }
    break;

  /* relative value: M(u)= |H(u)| / err*||u||_inf */
  case 1:
    if ( fabs(sol->umax) < 1e-30 )  return(1);
    maxu =  CTE2D / (info.eps * sol->umax);
    for (k=1; k<=mesh->np; k++) {
      der[k].hes[0] *= maxu; 
      der[k].hes[1] *= maxu;
      der[k].hes[2] *= maxu;
    }
    break;

  /* local norm: M(u)= |H(u)| / err*|u| */
  case 2:
    if ( sol->umax < 1.0e-30 )
      errs = 0.01;
    else
      errs = sol->umax * 0.01;
    for (k=1; k<=mesh->np; k++) {
      u     = fabs(getSol(sol,k,is));
      err1  = MS_MAX(errs,u);
      maxu  = CTE2D / err1;
      der[k].hes[0] *= maxu; 
      der[k].hes[1] *= maxu;
      der[k].hes[2] *= maxu;
    }
    break;

  /* local norm: M(u)= |H(u)| / err*(e1*|u|+e2*h*|du|) */
  case 3:
    maxg = 0.0;
    for (k=1; k<=mesh->np; k++) {
      u    = fabs(getSol(sol,k,is));
      norm = der[k].grd[0]*der[k].grd[0] + der[k].grd[1]*der[k].grd[1];
      norm = sqrt(norm);
      if ( norm > maxg )  maxg = norm;
    }
    if ( sol->umax < 1.0e-30 )
      errs = 0.01;
    else
      errs = sol->umax * 0.01;
    err1 = maxg * 0.01;

    for (k=1; k<=mesh->np; k++) {
      p0    = &mesh->point[k];
      u     = getSol(sol,k,is);
      norm  = der[k].grd[0]*der[k].grd[0] + der[k].grd[1]*der[k].grd[1];
      norm  = sqrt(norm);
      norm  = MS_MAX(err1,norm);
      norm *= p0->rins / p0->nv;
      err1  = MS_MAX(errs,fabs(u));
      
      /* variable weight */
      err1  = 0.01*err1 + 0.01*norm;
      maxu  = CTE2D / (info.eps*err1);
      der[k].hes[0] *= maxu;
      der[k].hes[1] *= maxu;
      der[k].hes[2] *= maxu;
    } 
    break;
  }

  return(1);
}


/* build metric for levelsets */
static double fsize(pMesh mesh,double u) {
  double   dd;

  if ( u < mesh->info.width )
    dd = mesh->info.hmin;
  else
	  dd = mesh->info.hmin + u * (mesh->info.hmax-mesh->info.hmin);
	
//  dd = MS_MAX(mesh->info.hmin,dd);
//	dd = MS_MIN(mesh->info.hmax,dd);
  return(dd);
}


/* compute ansotropic metric for levelset */
int metrLS_3d(pMesh mesh,pSol sol,pDeriv der) {
  double   *tmet,mat[6],mr[6],e1[3],e2[3],e3[3],dd,u,urel;
  double    tail,kappa,hh,hhmin,hhmax,lambda1,lambda2,lambda3,dx2,dy2,dz2,dxy,dxz,dyz;
  int       i,j,k,l,ias;

  hhmin = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  hhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);
	
  for (k=1; k<=mesh->np; k++) {
    u       = fabs(getSol(sol,k,0));
    urel    = u / sol->umax;
    hh      = fsize(mesh,urel);
    lambda3 = 1.0 / (hh*hh);
		kappa   = 0.0;

    /* curvature */
    dx2 = der[k].grd[0]*der[k].grd[0];
    dxy = der[k].grd[0]*der[k].grd[1];
		dxz = der[k].grd[0]*der[k].grd[2];
    dy2 = der[k].grd[1]*der[k].grd[1];
		dyz = der[k].grd[1]*der[k].grd[2];
		dz2 = der[k].grd[2]*der[k].grd[2];
		dd  = dx2 + dy2 + dz2;
    if ( dd > EPS1 ) {
      kappa = dx2*der[k].hes[3] + dz2*der[k].hes[3] \
            + dy2*der[k].hes[0] + dz2*der[k].hes[0] \
            + dx2*der[k].hes[5] + dy2*der[k].hes[5] \
            - 2.0 * (dxy*der[k].hes[1] + dxz*der[k].hes[2] + dyz*der[k].hes[4]);
      kappa = fabs(kappa) / 2.0 * sqrt(dd*dd*dd);
    }
    if ( urel < mesh->info.width ) {
      lambda1 = kappa / mesh->info.eps; //kappa*kappa*eps; //(eps*hh);
			lambda1 = MS_MIN(lambda1,hhmin);
			lambda1 = MS_MAX(lambda1,hhmax);
			lambda2 = lambda1;
    }
    else 
			lambda1 = lambda2 = hhmax;
    //lambda3 = MS_MIN(lambda3,hhmin);
    //lambda3 = MS_MAX(lambda3,hhmax);

    if ( mesh->info.iso ) {
      if ( lambda3 > lambda1 ) 
        tail = 1.0 / sqrt(lambda3);
      else
        tail = 1.0 / sqrt(lambda1);
      if ( !mesh->info.metric ) {
        sol->met[k] = tail;
      }
      else
        sol->met[k] = MS_MIN(sol->met[k],tail);
    }
    else {
      ias  = (k-1)*6 + 1;
      tmet = &sol->met[ias];

      /* compute local basis */
			memcpy(e3,der[k].grd,3*sizeof(double));
			dd = sqrt(dx2+dy2+dz2);
			e3[0] /= dd;
			e3[1] /= dd;
			e3[2] /= dd;
      if ( fabs(der[k].grd[0]) > EPS1 ) {
        e2[0] = -(der[k].grd[1]+der[k].grd[2]) / der[k].grd[0];
        e2[1] = e2[2] = 1.0;
      }
      else if ( fabs(der[k].grd[1]) > EPS1 ) {
        e2[1] = -(der[k].grd[0]+der[k].grd[2]) / der[k].grd[1];
        e2[0] = e2[2] = 1.0;
      }
      else if ( fabs(der[k].grd[2]) > EPS1 ) {
        e2[2] = -(der[k].grd[0]+der[k].grd[1]) / der[k].grd[2];
        e2[0] = e2[1] = 1.0;
      }
      else {
				e2[0] = 0.0; e2[1] = 0.; e2[2] = 1.0;
	    }
	    dd = sqrt(e2[0]*e2[0]+e2[1]*e2[1]+e2[2]*e2[2]);
			e2[0] /= dd;
			e2[1] /= dd;
			e2[2] /= dd;

      e1[0] = e2[1]*e3[2] - e2[2]*e3[1];
      e1[1] = e2[2]*e3[0] - e2[0]*e3[2];
      e1[2] = e2[0]*e3[1] - e2[1]*e3[0];

      if ( !mesh->info.metric ) {
        /* set coeffs directly */
        for (l=0,i=0; i<3; i++)
          for (j=i; j<3; j++)
           tmet[l++] = lambda1*e1[i]*e1[j] + lambda2*e2[i]*e2[j] + lambda3*e3[i]*e3[j];
      }
      else {
        for (l=0,i=0; i<3; i++)
          for (j=i; j<3; j++)
            mat[l++] = lambda1*e1[i]*e1[j] + lambda2*e2[i]*e2[j] + lambda3*e3[i]*e3[j];
        if ( !redsim(tmet,mat,mr) )  exit(0);
				for (i=0; i<6 ; i++)  tmet[i] = mr[i];
      }
    }
  }

  return(1);
}


int metrLS_2d(pMesh mesh,pSol sol,pDeriv der) {
  double   *tmet,mat[3],mr[3],u,urel;
  double    tail,kappa,hh,hhmin,hhmax,lambda1,lambda2,dd,dx2,dy2,dxy;
  int       k,ias;

  hhmin = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  hhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);

  for (k=1; k<=mesh->np; k++) {
		//mesh->info.ddebug = k == 207030;
    u       = fabs(getSol(sol,k,0));
		urel    = u / sol->umax;
    hh      = fsize(mesh,urel);
    lambda1 = 1.0 / (hh*hh);
    //kappa   = fabs(der[k].hes[0])+fabs(der[k].hes[2]);
    
    dx2 = der[k].grd[0]*der[k].grd[0];
    dxy = der[k].grd[0]*der[k].grd[1];
    dy2 = der[k].grd[1]*der[k].grd[1];
		dd  = dx2 + dy2;
    kappa = dx2*der[k].hes[2] - 2.0*dxy*der[k].hes[1] + dy2*der[k].hes[0];
    kappa = fabs(kappa) / 2.0 * sqrt(dd*dd*dd);
		
    if ( u < mesh->info.width*mesh->info.hmin ) {
      lambda2 = kappa / mesh->info.eps;
			lambda2 = MS_MAX(lambda2,hhmin);
			lambda2 = MS_MIN(lambda2,hhmax);
    }
    else {
      lambda2 = hhmax;
    }
    lambda2 = MS_MAX(lambda2,hhmin);   /* <-- chgt le 16/03/09 */
    lambda2 = MS_MIN(lambda2,hhmax);

    if ( mesh->info.iso ) {
      if ( lambda1 > lambda2 ) 
        tail = 1.0 / sqrt(lambda1);
      else
        tail = 1.0 / sqrt(lambda2);
      if ( !mesh->info.metric )
        sol->met[k] = tail;
      else
        sol->met[k] = MS_MIN(sol->met[k],tail);
    }
    else {
      ias  = (k-1)*3 + 1;
      tmet = &sol->met[ias];
      if ( !mesh->info.metric ) {
        tmet[0] = lambda1*dx2 + lambda2*dy2;
        tmet[1] = (lambda1-lambda2)*dxy;
        tmet[2] = lambda2*dx2 + lambda1*dy2;
      }
      else {
        mat[0] = lambda1*dx2 + lambda2*dy2;
        mat[1] = (lambda1-lambda2)*dxy;
        mat[2] = lambda2*dx2 + lambda1*dy2;
        if ( !redsim(tmet,mat,mr) )  return(0);
        tmet[0] = mr[0];
        tmet[1] = mr[1];
        tmet[2] = mr[2];
      }
    }
  }

  return(1);
}
