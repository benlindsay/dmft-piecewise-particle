#include "globals.h"

void calc_Unb() {

  int i ;
  // Take Gaussian potential to k-space
  fftw_fwd( uG , ktmp2 ) ;
  
  
  // Polymer chi A-B contribution //
  fftw_fwd( rho[0] , ktmp ) ;

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    ktmp[i] *= ktmp2[i] ;

  fftw_back( ktmp , tmp ) ;

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    tmp[i] *= rho[1][i] ;

  U_chi_gg = integrate( tmp ) * chiAB / rho0 ;


  // Polymer kappa contribution //
#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    tmp2[i] = rhot[i]  ;

  fftw_fwd( tmp2 , ktmp ) ;

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    ktmp[i] *= ktmp2[i] ;

  fftw_back( ktmp , tmp ) ;

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    tmp[i] *= tmp2[i] ;

  U_kappa_gg = integrate( tmp ) * kappa / 2.0 / rho0 ;

  // Polymer  particle kappa contribution //

fftw_fwd( uPG , ktmp2 ) ;


#pragma omp parallel for
    for ( i=0 ; i<M ; i++ )
        tmp2[i] = rhot[i]  ;
    fftw_fwd( tmp2 , ktmp ) ;

#pragma omp parallel for
    for ( i=0 ; i<M ; i++ )
        ktmp[i] *= ktmp2[i] ;
    
    fftw_back( ktmp , tmp ) ;

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    tmp[i] *= rho[2][i] ;




  U_kappa_pg = integrate( tmp ) * kappa / rho0 ;

// Polymer  particle chi  contribution //

#pragma omp parallel for
    for ( i=0 ; i<M ; i++ )
        tmp2[i] = rho[1][i]  ;
    fftw_fwd( tmp2 , ktmp ) ;

#pragma omp parallel for
    for ( i=0 ; i<M ; i++ )
        ktmp[i] *= ktmp2[i] ;
    
    fftw_back( ktmp , tmp ) ;

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    tmp[i] *= rho[2][i] ;




  U_chi_pg = integrate( tmp ) * chiAB/ rho0 ;

//Particle particle kappa  contribution //


  fftw_fwd( uP , ktmp2 ) ;

 #pragma omp parallel for
    for ( i=0 ; i<M ; i++ )
        tmp2[i] = rho[2][i]  ;
    fftw_fwd( tmp2 , ktmp ) ;

#pragma omp parallel for
    for ( i=0 ; i<M ; i++ )
        ktmp[i] *= ktmp2[i] ;
    
    fftw_back( ktmp , tmp ) ;

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    tmp[i] *= rho[2][i] ;

 U_kappa_pp = integrate( tmp ) * kappa_p/ rho0/2.0 ;
 
  
  Utt = U_kappa_pp + U_kappa_pg+ U_kappa_gg + U_chi_gg + U_chi_pg +Ubond;

}

