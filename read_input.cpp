#include "globals.h"

void read_input( void ) {

  FILE *inp ;
  int i,j;
  double d1 ;

  char tt[80] ;

  inp = fopen( "dyft.input" , "r" ) ;


  fscanf( inp , "%d %d" , &Nda , &Ndb ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf %d" , &phiHA , &Nha ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf %d" , &phiHB , &Nhb ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &CG_ratio ) ;
  fgets( tt , 80 , inp ) ;


  // Blank line //
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &phiP ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf %d %d" , &sigma, &Ng, &uni_sig ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &Rp ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &Xi ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &a ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%d" , &A_partics ) ;
  fgets( tt , 80 , inp ) ;


  // Blank line //
  fgets( tt , 80 , inp ) ;


  fscanf( inp , "%lf" , &C ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &chiAB ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &eps ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &well_shift ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf %lf" , &kappa ,&kappa_p) ;
  fgets( tt , 80 , inp ) ;
  cout<<"kappa "<<kappa<<" "<<kappa_p<<endl;
  Diff = ( double* ) calloc( 3 , sizeof( double ) ) ;
  fscanf( inp , "%lf %lf %lf %lf" , &Diff[0] , &Diff[1], &Diff[2],&Diff_rot ) ;
  fgets( tt , 80 , inp ) ;


  // Blank line //
  fgets( tt , 80 , inp ) ;


  for ( i=0 ; i<Dim ; i++ )
    fscanf( inp , "%lf" , &L[i] ) ;
  fgets( tt , 80 , inp ) ;

  for ( i=0 ; i<Dim ; i++ )
    fscanf( inp , "%d" , &Nx[i] ) ;
  fgets( tt , 80 , inp ) ;


  fscanf( inp , "%lf", &delt ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%d" , &pmeorder ) ;
  fgets( tt , 80 , inp ) ;


  // Blank line //
  fgets( tt , 80 , inp ) ;


  fscanf( inp , "%d" , &nsteps ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d" , &print_freq ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d" , &pre_equil_steps ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d" , &avg_freq ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d %d" , &sample_wait , &sample_freq ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d" , &stress_freq ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d" , &RAN_NUM ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d" , &charge_freq ) ;
  fclose( inp ) ;


}
