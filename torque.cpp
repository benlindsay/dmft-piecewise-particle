#include "globals.h"
void rot_noise(double** quatern_noise){
  int i,j,k,l;
  double tmp_gn[3];
  for(i=0; i<nP ;i++){
    for(j=0 ;j<3 ;j++){


      tmp_gn[j] = sqrt(2.0*Diff_rot*delt)* gasdev2();
    }

    for(j=0;j<4;j++){
      quatern_noise[i][j] = 0;
      for(k=0 ;k<3 ;k++){
        quatern_noise[i][j] += euler_B[i][j][k+1]*tmp_gn[k];
      }//k
    }//j=4 quaternions

  }//i = nP

}

void torque(void){

  int i,j,m,k,l, gind, t1,t2, center_ind,ind;
  double *tmp_trp,mt, mdr2, mdrp2, mdr, dr[Dim], drp[Dim],delr ;
  int sites_per_gnp = 1 + ng_per_partic * ( 1 + Ng );
  // set all the torque to zero//

#pragma omp parallel for private(j)
  for(i=0;i<nP;i++){
    for(j = 0; j<4 ;j++)
      trq[i][j] = 0;
    for(j = 0; j<3; j++)
      real_trq[i][j] = 0;
  }


  if(sites_per_gnp >0){
#pragma omp parallel for private(dr,drp,m,j,k,l,ind,gind,center_ind,mdrp2,mdr2)\
    reduction(+:Ubond)
    for(i=0;i<nP;i++){
      center_ind = nD * ( Nda + Ndb ) + nA * Nha + nB * Nhb + sites_per_gnp * i ;
      //real torque on each nP
      for( m=0 ; m<ng_per_partic ; m++){
        ind = center_ind + m * ( Ng + 1 ) + 1 ;
        gind = ng_per_partic*i + m;
        mdr2 = pbc_mdr2( x[ind+1] , x[ind] , dr ) ;//dist between r' and r
        mdrp2 = pbc_mdr2(x[ind],x[center_ind],drp) ; //dist between rcm and r
        if(tp[ind] != -1) { cout<<"wrong trq calc !!"<<endl; exit(1);}

        Ubond += 1.5*mdr2;
        for(j=0;j<3;j++){

          for(k=0 ; k<3 ; k++){
            for(l=0; l<3 ; l++){
              real_trq[i][j] += 3.0*epslon[j][k][l]*drp[k]*dr[l];
              //	real_trq[i][j] += epslon[j][k][l]*drp[k]*f[ind][l];
            }//l



          }//k

        }//j, torque dim




      }//m=ng_per_oartic

      //get the quaternions' langevin force term, trq//
      for( m = 0; m<4; m++){
        for(j= 0 ; j<3 ; j++){
          for(k= 0 ; k<3 ;k++){
            trq[i][m] += euler_B[i][m][j+1]*euler_A[i][j][k]*real_trq[i][k];

          }//k
        }//j
      }//m = 4 quaternions

    }//nP
  }//site_np>0?


}
