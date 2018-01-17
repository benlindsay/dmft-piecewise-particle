#include "globals.h"
void charge_grid( void ) ;
void bonds( void ) ;


void forces() {

  int i,j, m, gind, t1, t2 ;

  charge_grid() ;

  ///////////////////////////////////////////////
  // Reset the particle forces and grid grad w //
  ///////////////////////////////////////////////
#pragma omp parallel for private(j)
  for ( i=0 ; i<nstot ; i++ )
    for ( j=0 ; j<Dim ; j++ )
      f[i][j] = 0.0 ;

#pragma omp parallel for private(j)
  for ( i=0 ; i<M ; i++ )
    for ( j=0 ; j<Dim ; j++ ) 
      gradwA[j][i] = gradwB[j][i] = gradwP[j][i] = 0.0 ;



  //////////////////////////////////////////////////
  // Accumulate the monomer-monomer contributions //
  //////////////////////////////////////////////////

  if ( chiAB != 0.0 ) {
    // A acting on B //
    fftw_fwd( rho[0] , ktmp ) ;

    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwB[j][i] += tmp[i] * chiAB / rho0 ;
    }


    // B acting on A //
    fftw_fwd( rho[1] , ktmp ) ;

    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwA[j][i] += tmp[i] * chiAB / rho0 ;
    }
  }

  if ( eps != 0.0 ) {
    // P pulling on A //
    fftw_fwd( rho[2] , ktmp ) ;

    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uAG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwA[j][i] -= tmp[i] * eps/rho0 ;
    }
  } 

  if (eps != 0.0) {
    // A pulling on P //
    fftw_fwd( rho[0] , ktmp ) ;

    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uAG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwP[j][i] -= tmp[i] * eps/rho0 ;
    }
  }

  // Compressibility contribution //
#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    tmp[i] = rhot[i] -rho0 ;

  fftw_fwd( tmp , ktmp ) ;

  for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
    for ( i=0 ; i<M ; i++ )
      ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;

    fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
    for ( i=0 ; i<M ; i++ ) {
      gradwA[j][i] += tmp[i] * kappa / rho0 ;
      gradwB[j][i] += tmp[i] * kappa / rho0 ;
    }

  }


  ///////////////////////////////////////////
  // Accumulate the particle contributions //
  ///////////////////////////////////////////
  if ( nP > 0 ) {

    // Particle-particle //
    fftw_fwd( rho[2] , ktmp ) ;


    if(nP >1){
      for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
        for ( i=0 ; i<M ; i++ )
          ktmp2[i] = grad_uP_hat[j][i] * ktmp[i] ;

        fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
        for ( i=0 ; i<M ; i++ )
          gradwP[j][i] += tmp[i] * kappa_p / rho0 ;
      }
    }//nP>1
    // Particles acting on monomers //
    for ( j=0 ; j<Dim ; j++ ) {
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ ) {
        gradwA[j][i] += tmp[i] * kappa / rho0 ;
        gradwB[j][i] += tmp[i] * ( kappa + ( A_partics ? chiAB : 0.0 ) ) / rho0 ;
      }
    }

    // A Monomers acting on particles //
#pragma omp parallel for
    for ( i=0 ; i<M ; i++ )
      tmp[i] = rho[0][i] -rho0 ;


    fftw_fwd( tmp , ktmp ) ;
    for ( j=0 ; j<Dim ; j++ ) {

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwP[j][i] += tmp[i] * kappa / rho0 ;

    }

    // B Monomers acting on particles //
    fftw_fwd( rho[1] , ktmp ) ;
    for ( j=0 ; j<Dim ; j++) {

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        gradwP[j][i] += tmp[i] * ( kappa + ( A_partics ? chiAB : 0.0 ) ) / rho0 ;

    }


  }// if ( nP > 0 )




  //////////////////////////////////////////////////////
  // Accumulate the nonbonded forces on each particle //
  //////////////////////////////////////////////////////
#pragma omp parallel for private(m,j)
  for ( i=0 ; i<nstot ; i++ ) {
    for ( m=0 ; m < grid_per_partic ; m++ ) {

      gind = grid_inds[ i ][ m ] ;

      for ( j=0 ; j<Dim ; j++ ) {
        if ( tp[i] == 0 )
          f[i][j] -= gradwA[ j ][ gind ] * grid_W[i][m] ;
        else if ( tp[i] == 1 )
          f[i][j] -= gradwB[ j ][ gind ] * grid_W[i][m] ;
        else if ( tp[i] == 2 )
          f[i][j] -= gradwP[ j ][ gind ] * grid_W[i][m] ;
      }
    }

    for ( j=0 ; j<Dim ; j++ )
      f[i][j] *= gvol ;
  }


  ////////////////////////////
  // Call the bonded forces //
  ////////////////////////////
  bonds() ;

}
