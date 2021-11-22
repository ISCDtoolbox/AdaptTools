#include "mshmet.h"


/* read mesh data */
int loadMesh(pMesh mesh,char *filename) {
  pPoint       ppt,p0,p1,p2;
  pTetra       pt;
  pTria        pt1;
  double       ux,uy,h1,h2,h3,pe,rins;
  float        fp1,fp2,fp3;
  int          i,k,inm,ref;
  char        *ptr,data[128];

  strcpy(data,filename);
  ptr = strstr(data,".mesh");
  mesh->info.bin = 0;
  if ( !ptr ) {
    strcat(data,".meshb");
    mesh->info.bin = 1;
    if (!(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if (!(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(0);
      }
      mesh->info.bin = 0;
    }
  }
  else if (!(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
    fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"  -- READING DATA FILE %s\n",data);

  mesh->np = GmfStatKwd(inm,GmfVertices);
  mesh->nt = GmfStatKwd(inm,GmfTriangles);
  mesh->ne = GmfStatKwd(inm,GmfTetrahedra);
  if ( !mesh->np || mesh->ne+mesh->nt == 0 ) {
    fprintf(stdout,"  ** MISSING DATA\n");
    return(0);
  }

  /* mem alloc */
  mesh->point = (pPoint)M_calloc(mesh->np+1,sizeof(Point),"point");
  assert(mesh->point);
  if ( mesh->ne ) {
    mesh->tetra = (pTetra)M_calloc(mesh->ne+1,sizeof(Tetra),"tetra");
    assert(mesh->tetra);
    mesh->adja = (int*)M_calloc(4*mesh->ne+5,sizeof(int),"adja");
    assert(mesh->adja);
    mesh->nt = 0;
  }
  else if ( mesh->nt ) {
    mesh->tria  = (pTria)M_calloc(mesh->nt+1,sizeof(Tria),"tria");
    assert(mesh->tria);
    mesh->adja = (int*)M_calloc(3*mesh->nt+5,sizeof(int),"adja");
    assert(mesh->adja);
  }

  /* read mesh vertices */
  GmfGotoKwd(inm,GmfVertices);
  if ( mesh->dim == 2 ) {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( mesh->ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ref);
    }
  }
  else {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( mesh->ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
        ppt->c[2] = fp3;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ref);
    }
  }

  /* read mesh triangles */
  GmfGotoKwd(inm,GmfTriangles);
  if ( mesh->dim == 2 ) {
    for (k=1; k<=mesh->nt; k++) {
      pt1 = &mesh->tria[k];
      GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&ref);
      for (i=0; i<3; i++) {    
              ppt = &mesh->point[pt1->v[i]];
        if ( !ppt->s )  ppt->s = k;
      }
      p0 = &mesh->point[pt1->v[0]];
      p1 = &mesh->point[pt1->v[1]];
      p2 = &mesh->point[pt1->v[2]];

      ux = p1->c[0] - p0->c[0];
      uy = p1->c[1] - p0->c[1];
      h1 = sqrt(ux*ux + uy*uy);

      ux = p2->c[0] - p0->c[0];
      uy = p2->c[1] - p0->c[1];
      h2 = sqrt(ux*ux + uy*uy);

      ux = p2->c[0] - p1->c[0];
      uy = p2->c[1] - p1->c[1];
      h3 = sqrt(ux*ux + uy*uy);

      pe        = 0.5 * (h1 + h2 + h3);
      pt1->aire = pe * (pe-h1) * (pe-h2) * (pe-h3);
      pt1->aire = sqrt(pt1->aire);
      rins      = 2.0 * pt1->aire / pe;

      p0->aire += pt1->aire;
      p0->rins += rins;

      p1->aire += pt1->aire;
      p1->rins += rins;

      p2->aire += pt1->aire;
      p2->rins += rins;
    }
  }
  else {
    for (k=1; k<=mesh->nt; k++) {
      pt1 = &mesh->tria[k];
      GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&ref);
    }
  }

  /* read mesh tetrahedra */
  GmfGotoKwd(inm,GmfTetrahedra);
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    GmfGetLin(inm,GmfTetrahedra,&pt->v[0],&pt->v[1],&pt->v[2],&pt->v[3],&ref);
    for (i=0; i<4; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( mesh->dim == 3 && !ppt->s )  ppt->s = k;
    }
  }

  GmfCloseMesh(inm);
  return(1);
}


/* load solution (metric) */
int loadSol(pSol sol,Info *info,char *filename) {
  double       dbuf[ GmfMaxTyp ];
  float        fbuf[ GmfMaxTyp ];
  int          k,i,ia,inm;
  char        *ptr,data[128];

  strcpy(data,filename);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  strcat(data,".solb");
  if (!(inm = GmfOpenMesh(data,GmfRead,&sol->ver,&sol->dim)) ) {
    ptr  = strstr(data,".solb");
    *ptr = '\0';
    strcat(data,".sol");
    if (!(inm = GmfOpenMesh(data,GmfRead,&sol->ver,&sol->dim)) ) {
      fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
      return(0);
    }
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  if ( abs(info->imprim) > 3 )
    fprintf(stdout,"  -- READING DATA FILE %s\n",data);

  sol->np = GmfStatKwd(inm,GmfSolAtVertices,&sol->type,&sol->size,sol->typtab);
  if ( !sol->np ) {
    fprintf(stdout,"  ** MISSING DATA.\n");
    return(0);
  }

  /* mem alloc */
  sol->sol = (double*)M_calloc(sol->np+1,sol->size*sizeof(double),"inout");
  assert(sol->sol);

  /* read mesh solutions */
  GmfGotoKwd(inm,GmfSolAtVertices);
  for (k=1; k<=sol->np; k++) {
    if ( sol->ver == GmfFloat ) {
      GmfGetLin(inm,GmfSolAtVertices,fbuf);
      ia = (k-1)*sol->size + 1;
      for (i=0; i<sol->size; i++)
        sol->sol[ia+i] = fbuf[i];   // < -0.0001 ? -1.0 : fbuf[i];
    }
    else {
      GmfGetLin(inm,GmfSolAtVertices,dbuf);
      ia = (k-1)*sol->size + 1;
      for (i=0; i<sol->size; i++)
        sol->sol[ia+i] = dbuf[i];
    }
  }

  GmfCloseMesh(inm);
  return(1);  
}


/* load solution (metric) */
int loadMetric(pSol sol,Info *info,char *filename) {
  double   dbuf[ GmfMaxTyp ];
  float    fbuf[ GmfMaxTyp ];
  int      i,k,inm,ia,np,ver,dim,type,size,typtab[ GmfMaxTyp ];

  if ( !(inm = GmfOpenMesh(filename,GmfRead,&ver,&dim)) ) {
    fprintf(stderr,"  ** %s  NOT FOUND.\n",filename);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",filename);

  if ( abs(info->imprim) > 3 )
    fprintf(stdout,"  -- READING METRIC FILE %s\n",filename);

  np = GmfStatKwd(inm,GmfSolAtVertices,&type,&size,typtab);
  if ( !np || np != sol->np ) {
    fprintf(stdout,"  ** INCOMPLETE DATA.\n");
    return(0);
  }

  /* mem alloc */
  sol->met = (double*)M_calloc(np+1,size*sizeof(double),"inout");
  assert(sol->met);

  /* read mesh solutions */
  GmfGotoKwd(inm,GmfSolAtVertices);
  for (k=1; k<=np; k++) {
    if ( ver == GmfFloat ) {
      GmfGetLin(inm,GmfSolAtVertices,fbuf);
      ia = (k-1)*size + 1;
      for (i=0; i<size; i++)
        sol->met[ia+i] = fbuf[i];
    }
    else {
      GmfGetLin(inm,GmfSolAtVertices,dbuf);
      ia = (k-1)*size + 1;
      for (i=0; i<size; i++)
        sol->met[ia+i] = dbuf[i];
    }
  }

  GmfCloseMesh(inm);
  return(1);  
}


int saveMet(pSol sol,Info *info,char *filename) {
  double       dbuf[ GmfMaxTyp ],tmpd;
  float        fbuf[ GmfMaxTyp ],tmpf;
  int          k,i,ia,outm;
  char        *ptr,data[128];

  strcpy(data,filename);
  ptr = strstr(data,".mesh");
  if ( ptr )  {
    *ptr = '\0';
    strcat(data,".solb");
  }
  else {
    ptr = strstr(data,".sol");
    if ( !ptr )  strcat(data,".sol");
  }

  if (!(outm = GmfOpenMesh(data,GmfWrite,sol->ver,sol->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }

  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /* write sol */
  sol->type = 1;
  if ( info->iso ) {
    sol->typtab[0] = GmfSca;
    sol->size = 1;
  }
  else {
    sol->typtab[0] = GmfSymMat;
    sol->size = sol->dim*(sol->dim+1) / 2;
  }
  GmfSetKwd(outm,GmfSolAtVertices,sol->np,sol->type,sol->typtab);
  for (k=1; k<=sol->np; k++) {
    ia = (k-1)*sol->size + 1;
    if ( sol->ver == GmfFloat ) {
      for (i=0; i<sol->size; i++)
        fbuf[i] = sol->met[ia+i];
      if ( sol->dim == 3 ) {
                                tmpf    = fbuf[3];
                                fbuf[3] = fbuf[2];
                                fbuf[2] = tmpf;
            }
      GmfSetLin(outm,GmfSolAtVertices,fbuf);
    }
    else {
      for (i=0; i<sol->size; i++)
        dbuf[i] = sol->met[ia+i];
      if ( sol->dim == 3 ) {
                                tmpd    = dbuf[3];
                                dbuf[3] = dbuf[2];
                                dbuf[2] = tmpd;
            }
      GmfSetLin(outm,GmfSolAtVertices,dbuf);
    }
  }

  GmfCloseMesh(outm);
  return(1);
}


int saveSol(pSol sol,Info *info,char *filename) {
  double       dbuf[ GmfMaxTyp ];
  float        fbuf[ GmfMaxTyp ];
  int          k,i,ia,outm;
  char        *ptr,data[128];

  strcpy(data,filename);
  ptr = strstr(data,".mesh");
  if ( ptr )  {
    *ptr = '\0';
    strcat(data,".solb");
  }
  else {
    ptr = strstr(data,".sol");
    if ( !ptr )  strcat(data,".sol");
  }

  if (!(outm = GmfOpenMesh(data,GmfWrite,sol->ver,sol->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }

  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /* write sol */
  sol->type = 1;
  sol->typtab[0] = GmfSca;
  sol->size = 1;
  GmfSetKwd(outm,GmfSolAtVertices,sol->np,sol->type,sol->typtab);
  for (k=1; k<=sol->np; k++) {
    ia = (k-1)*sol->size + 1;
    if ( sol->ver == GmfFloat ) {
      for (i=0; i<sol->size; i++)
        fbuf[i] = sol->sol[ia+i];
      GmfSetLin(outm,GmfSolAtVertices,fbuf);
    }
    else {
      for (i=0; i<sol->size; i++)
        dbuf[i] = sol->sol[ia+i];
      GmfSetLin(outm,GmfSolAtVertices,dbuf);
    }
  }

  GmfCloseMesh(outm);
  return(1);
}

int outder(pMesh mesh,pDeriv der,int j) { 
        double       dd;
  float        fbuf[ GmfMaxTyp ];
  int          i,k,outm,ver,siz,type,typtab[ GmfMaxTyp ];

  ver = GmfFloat;
  if ( !(outm = GmfOpenMesh("gradient.sol",GmfWrite,ver,mesh->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN gradient\n");
    return(0);
  }
  type = 1;
  typtab[0] = GmfVec;
  GmfSetKwd(outm,GmfSolAtVertices,mesh->np,type,typtab);
  if ( mesh->dim == 2) {
    for (k=1; k<=mesh->np; k++) {
                  dd = 1.0 / sqrt(der[k].grd[0]*der[k].grd[0] + der[k].grd[1]*der[k].grd[1]);
      fbuf[0] = dd * der[k].grd[0];
      fbuf[1] = dd * der[k].grd[1];
      GmfSetLin(outm,GmfSolAtVertices,fbuf);
    }
  }
        else {
    for (k=1; k<=mesh->np; k++) {
                  dd = 1.0 / sqrt(der[k].grd[0]*der[k].grd[0] + der[k].grd[1]*der[k].grd[1] \
                                + der[k].grd[2]*der[k].grd[2]);
      fbuf[0] = dd * der[k].grd[0];
      fbuf[1] = dd * der[k].grd[1];
      fbuf[2] = dd * der[k].grd[2];
      GmfSetLin(outm,GmfSolAtVertices,fbuf);
    }
        }
  GmfCloseMesh(outm);

  if ( !(outm = GmfOpenMesh("hessien.sol",GmfWrite,ver,mesh->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN hessien\n");
    return(0);
  }
  type = 1;
  typtab[0] = GmfSymMat;
  GmfSetKwd(outm,GmfSolAtVertices,mesh->np,type,typtab);
        siz = mesh->dim == 2 ? 3 : 6;
  for (k=1; k<=mesh->np; k++) {
                for (i=0; i<siz; i++)
      fbuf[i] = der[k].hes[i];
    GmfSetLin(outm,GmfSolAtVertices,fbuf);
  }
  GmfCloseMesh(outm);

  if ( mesh->info.ls ) {
    if ( !(outm = GmfOpenMesh("curvature.sol",GmfWrite,ver,mesh->dim)) ) {
      fprintf(stderr,"  ** UNABLE TO OPEN hessien\n");
      return(0);
    }
    type = 1;
    typtab[0] = GmfSca;
    GmfSetKwd(outm,GmfSolAtVertices,mesh->np,type,typtab);
    for (k=1; k<=mesh->np; k++) {      
      fbuf[0] = fabs(der[k].hes[0]) + fabs(der[k].hes[2]);
                  GmfSetLin(outm,GmfSolAtVertices,fbuf);
    }   
    GmfCloseMesh(outm);
  }
  return 1;
}


