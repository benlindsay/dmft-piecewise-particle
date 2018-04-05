#include "globals.h"

void update_positions() {
  /*update the position using SV-GJF*/
  int i,j ;
  double tmp_d,tmp_pm,tmp_x,tmp_gn,delt_2 = delt*delt;
  //#pragma openmp parallel for \
  //  private(i, j)
  for ( i=0 ; i<nstot ; i++ ) {

    //check if tp[i] ==2
    if(tp[i]==2)
      tmp_pm = p_m;
    else
      tmp_pm = 1;

    if(tp[i] ==-1)
      continue;


    for ( j=0 ; j<Dim ; j++ ) {

      tmp_x = x[i][j];
      tmp_gn = gasdev2();

      tmp_d = x[i][j] - x_bac[i][j];

      if(tmp_d >Lh[j])
        x_bac[i][j] += L[j];
      else if(tmp_d <-Lh[j])
        x_bac[i][j] -= L[j];
      else{

      }



      x[i][j] = 2.0*verlet_b[ tp[i] ]*x[i][j] - verlet_a[ tp[i] ]*x_bac[i][j] + verlet_b[ tp[i] ]*delt_2*f[i][j]/tmp_pm + verlet_b[ tp[i] ]*delt*sqrt( 2.0/Diff[ tp[i] ] * delt )*(tmp_gn+gn_bac[i][j])/2.0/tmp_pm;

      if ( x[i][j] > L[j] ) {
        x[i][j] -= L[j] ;
        if (x[i][j] > L[j]) {
          printf("Particle %d is out of bounds (x[%d][%d] = %lf).\n",
                 i, i, j, x[i][j]);
          die("Try reducing timestep.\n");
        }
      }

      else if ( x[i][j] < 0.0 ) {
        x[i][j] += L[j] ;
        if (x[i][j] < 0) {
          printf("Particle %d is out of bounds (x[%d][%d] = %lf).\n",
                 i, i, j, x[i][j]);
          die("Try reducing timestep.\n");
        }
      }

      x_bac[i][j] = tmp_x;
      gn_bac[i][j] = tmp_gn;
    }

  }//for ( i=0 ; ...


}

void update_positions_init() {
  /*initialize the position using using EM */
  int i,j ;

  double tmp_pm ,tmp_gn,delt_2 = delt*delt;


  //#pragma openmp parallel for \
  //  private(i, j)
  for ( i=0 ; i<nstot ; i++ ) {

    if(tp[i] ==2)
      tmp_pm = p_m;
    else
      tmp_pm =1 ;

    if(tp[i] == -1 )
      continue;

    for ( j=0 ; j<Dim ; j++ ) {
      x_bac[i][j] = x[i][j] ;

      tmp_gn=gasdev2();

      gn_bac[i][j] = tmp_gn;

      x[i][j] = x[i][j] + verlet_b[ tp[i] ]*delt_2*f[i][j]/tmp_pm/2.0 + verlet_b[ tp[i] ]*delt*sqrt( 2.0/Diff[ tp[i] ] * delt )*tmp_gn/2.0/tmp_pm;

      // x[i][j] = x[i][j] + Diff[ tp[i] ] * delt * f[i][j]
      //       + sqrt( 2.0 * Diff[ tp[i] ] * delt ) * gasdev2() ;

      if ( x[i][j] > L[j] )
        x[i][j] -= L[j] ;

      else if ( x[i][j] < 0.0 )
        x[i][j] += L[j] ;
    }

  }//for ( i=0 ; ...


}
