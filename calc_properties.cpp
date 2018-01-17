#include "globals.h"

void calc_properties( ) {
  int i, j, k, m, t1, t2 ;
  
/*

  // Evaluate rho(r) * u(r-r') * rho(r')
  // t1: loops over distinct species for rho(r')
  
  Unb = 0.0 ;
  for ( t1=0 ; t1<ntypes ; t1++ ) {

    convolve_fields( rho[t1] , u , tmp ) ;

    for ( t2=t1 ; t2<ntypes ; t2++ ) {

      for ( i=0 ; i<M ; i++ ) 
        tmp[i] *= rho[t2][i] ;

      Unb += integ_trapPBC( tmp ) ;

    }
  }

  
  Pvir = 0.0 ;

  */
}
