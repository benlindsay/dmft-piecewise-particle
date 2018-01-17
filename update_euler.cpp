#include "globals.h"
void calc_A();
void rot_noise(double**);
double crct_q(double* , double* );

void update_euler(){

  int i, j,k, cent_ind,gind, ind;
  int sites_per_gnp = 1 + ng_per_partic * ( 1 + Ng );
  double **tmp_q,lamb;

  tmp_q = (double**) calloc(nP, sizeof(double*));
  for(i=0 ; i<nP ; i++)
    tmp_q[i] = (double*) calloc(4, sizeof(double));

  rot_noise(q_noise);

#pragma omp parallel for private(j,lamb)
  for(i=0 ; i<nP ; i++){
    for(j=0;j<4;j++){

      //get the predicted quaterns, q~//	
      tmp_q[i][j] =  euler_q[i][j] + Diff_rot*delt*trq[i][j] + q_noise[i][j];	


    }//j
    lamb = crct_q(euler_q[i],tmp_q[i]); // Lagrange multiplier , lambda


    for(j=0;j<4;j++){

      euler_q[i][j] = tmp_q[i][j] + lamb*euler_q[i][j]; //unit length correction 
    }//j

  }//i=nP

  calc_A();




#pragma omp parallel for private(k,j,cent_ind,ind,gind)
  for(i = 0 ; i<nP ;i++){
    cent_ind = nD * ( Nda + Ndb ) + nA * Nha + nB * Nhb + sites_per_gnp * i ;
    for(j = 0; j<ng_per_partic; j++){
      ind = cent_ind + j * ( Ng + 1 ) + 1 ; 
      gind = ng_per_partic*i + j;
      matrix_trs_vect(euler_A[i],grf_bf_x[gind],x[ind]);

      if(tp[ind] != -1) { cout<<"wrong em update !!"<<endl; exit(1);} 

      for(k=0 ;k<Dim; k++){
        x[ind][k] += x[cent_ind][k] ;
        if ( x[ind][k] > L[k] )
          x[ind][k] -= L[k] ;
        else if ( x[ind][k] < 0.0 )
          x[ind][k] += L[k] ;


      }
    }//ng_per_partic


  }//nP

  for(i=0 ;i<nP;i++)
    free(tmp_q[i]);
  free(tmp_q);

}


double crct_q(double* q, double* q_tld){

  double delta,q_qtld=0, qtld2=0;
  int i;
  for(i=0 ; i<4 ;i++){

    q_qtld += q[i]*q_tld[i];
    qtld2 += q_tld[i]*q_tld[i];

  }

  delta = q_qtld*q_qtld - qtld2 +1;

  if(delta <0 ){
    cout<<"can not find lamb_a !"<<endl;
    exit(1);
  }

  return (-q_qtld+sqrt(delta));
}
