#include "globals.h"

void write_stress( ) {

  int i, j, k ;

  FILE *otp,*otp_pp,*otp_ng ;

  if ( step <= print_freq )
    otp = fopen( "stress_euler.dat" , "w" ) ;
  else
    otp = fopen( "stress_euler.dat" , "a" ) ;

  if ( step <= print_freq )
    otp_pp = fopen( "stress_PP.dat" , "w" ) ;
  else
    otp_pp = fopen( "stress_PP.dat" , "a" ) ;

  if ( step <= print_freq )
    otp_ng = fopen( "stress_ga.dat" , "w" ) ;
  else
    otp_ng = fopen( "stress_ga.dat" , "a" ) ;



  for ( i=0 ; i<buff_ind ; i++ ) {
    // Write diagonals first //
    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%lf " , sts_buf[i][j][j] ) ;

    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp_pp , "%lf " , sts_buf_pp[i][j][j] ) ;

    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp_ng , "%lf " , sts_buf_ng[i][j][j] ) ;



    for ( j=0 ; j<Dim  ; j++ )
      for ( k=j+1 ; k<(Dim) ; k++ )
        fprintf( otp , "%lf " , sts_buf[i][j][k] ) ;

    for ( j=0 ; j<Dim ; j++ )
      for ( k=j+1 ; k<(Dim) ; k++ )
        fprintf( otp_pp , "%lf " , sts_buf_pp[i][j][k] ) ;

    for ( j=0 ; j<Dim  ; j++ )
      for ( k=j+1 ; k<(Dim) ; k++ )
        fprintf( otp_ng , "%lf " , sts_buf_ng[i][j][k] ) ;


    // for ( j=0 ; j<3  ; j++ )
    /* for(k=Dim ; k<(Dim+nP) ; k++)
       for ( j=0 ; j<4  ; j++ )
       fprintf( otp , "%lf " , sts_buf[i][j][k] ) ;
       */
    fprintf( otp , "\n" ) ;	
    fprintf( otp_pp , "\n" ) ;
    fprintf( otp_ng , "\n" ) ;
  }

  fclose( otp ) ;
  fclose( otp_pp ) ;
  fclose( otp_ng ) ;

  buff_ind = 0 ;

}




void write_kspace_data( const char *nm , complex<double> *kdt ) {
  int i, j , nn[Dim] ;
  FILE *otp ;
  double kv[Dim], k2 ;

  otp = fopen( nm , "w" ) ;

  for ( i=1 ; i<M ; i++ ) {
    unstack( i , nn ) ;

    k2 = get_k( i , kv ) ;

    for ( j=0 ; j<Dim ; j++ ) 
      fprintf( otp , "%lf " , kv[j] ) ;

    fprintf( otp , "%1.5e %1.5e %1.5e %1.5e\n" , abs(kdt[i]), sqrt(k2), 
        real(kdt[i]) , imag(kdt[i]) ) ;

    if ( Dim == 2 && nn[0] == Nx[0]-1 )
      fprintf( otp , "\n" ) ;
  }

  fclose( otp ) ;


}

void write_avg_grid_data( const char *nm , double *dat ) {

  int i, j, nn[Dim] ;
  FILE *otp ;
  double r[Dim] ;
  otp = fopen( nm , "w" ) ;

  for ( i=0 ; i<M ; i++ ) {
    unstack( i , nn ) ;

    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%lf " , double(nn[j]) * dx[j] ) ;

    fprintf( otp , "%1.16e \n" , dat[i]/n_samples ) ;

    if ( Dim == 2 && nn[0] == Nx[0]-1 )
      fprintf( otp , "\n" ) ;
  }

  fclose( otp ) ;

}

void write_avg_grid_data_iter( const char *nmi , double *dat ) {

  int i, j, nn[Dim] ;
  FILE *otp ;
  double r[Dim] ;
  char nm[50];
  sprintf(nm, "%s_%d.dat", nmi, avg_iter);
  otp = fopen( nm , "w" ) ;

  for ( i=0 ; i<M ; i++ ) {
    unstack( i , nn ) ;

    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%lf " , double(nn[j]) * dx[j] ) ;

    fprintf( otp , "%1.16e \n" , dat[i]/n_samples_iter ) ;

    if ( Dim == 2 && nn[0] == Nx[0]-1 )
      fprintf( otp , "\n" ) ;
  }

  fclose( otp ) ;

}

void write_grid_data( const char *nm , double *dat ) {

  int i, j, nn[Dim] ;
  FILE *otp ;
  double r[Dim] ;
  otp = fopen( nm , "w" ) ;

  for ( i=0 ; i<M ; i++ ) {
    unstack( i , nn ) ;

    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%lf " , double(nn[j]) * dx[j] ) ;

    fprintf( otp , "%1.16e \n" , dat[i] ) ;

    if ( Dim == 2 && nn[0] == Nx[0]-1 )
      fprintf( otp , "\n" ) ;
  }

  fclose( otp ) ;

}


