#include "globals.h"
void calc_A(void) ;
void fibonacci_u(double *, int );
void unif_sig (double *, int );
void charge_grid( void ) ;
void write_grid_data( const char* , double* ) ;
void random_u( double* ) ;
void random_config( void ) ;
void write_gro( void ) ;
void read_gro( FILE* ) ;
void random_config( void ) ;
void read_quaternions(FILE *);


void initialize_configuration( ) {

  int i, j ;

  FILE *inp ,*inp2;
  inp = fopen( "input.gro" , "r" ) ;
  inp2 = fopen( "input_q.gro" , "r" ) ;
  //cout<<"here reach "<<endl;
  if ( inp == NULL || inp2 ==NULL){
    random_config() ;
    //calc_A();
    rst_para = 0 ;
    eps_flag = 1 ;
  }
  else {
    read_gro( inp ) ;
    read_quaternions( inp2 );
    fclose( inp ) ;
    fclose( inp2 ) ;
    printf("input.gro read!\n" ) ;
    printf("input_q.gro read!\n" ) ;
    rst_para =1 ;
    eps_flag = 0 ;
  }


}

void read_gro( FILE *ip ) {

  int i, m, j, k, di , ind ;
  char tt[80] ;
  fgets( tt , 80 , ip ) ;
  fscanf( ip , "%d\n" , &di ) ;

  if ( di != nstot )
    die("Number of sites in input.gro does not match!\n");

  ind = 0 ;

  int a = 0;
  for ( k=0 ; k<nD ; k++ ) {
    for ( m=0 ; m<Nda + Ndb ; m++ ) {
      a++;
      fscanf( ip , "%5d" , &di ) ;
      fscanf( ip , "%s" , tt ) ;
      fscanf( ip , "%s" , tt ) ;
      if (a < 10000 )
        fscanf( ip , "%d" , &di ) ;
      if (a % 100000 == 0)
        a = 0;

      for ( j=0 ; j<Dim ; j++ ) {
        fscanf( ip , "%lf" , &x[ind][j] ) ;
      }
      for ( j=0 ; j<Dim ; j++ ) {
        fscanf( ip , "%lf" , &x_bac[ind][j] ) ;
      }
      for ( j=0 ; j<Dim ; j++ ) {
        fscanf( ip , "%lf" , &gn_bac[ind][j] ) ;
      }




      fgets( tt, 80 , ip ) ;

      tp[ ind ] = ( m < Nda ? 0 : 1 ) ;

      ind++ ;
    }
  }

  for ( k=0 ; k<nA ; k++ ) {
    for ( m=0 ; m<Nha ; m++ ) {
      a++;
      fscanf( ip , "%5d" , &di ) ;
      fscanf( ip , "%s" , tt ) ;
      fscanf( ip , "%s" , tt ) ;
      if (a < 10000 )
        fscanf( ip , "%d" , &di ) ;
      if (a % 100000 == 0)
        a = 0;

      for ( j=0 ; j<Dim ; j++ ) {
        fscanf( ip , "%lf" , &x[ind][j] ) ;
      }
      for ( j=0 ; j<Dim ; j++ ) {
        fscanf( ip , "%lf" , &x_bac[ind][j] ) ;
      }
      for ( j=0 ; j<Dim ; j++ ) {
        fscanf( ip , "%lf" , &gn_bac[ind][j] ) ;
      }


      fgets( tt, 80 , ip ) ;

      tp[ ind ] = 0 ;

      ind++ ;
    }
  }

  for ( k=0 ; k<nB ; k++ ) {
    for ( m=0 ; m<Nhb ; m++ ) {
      a++;
      fscanf( ip , "%5d" , &di ) ;
      fscanf( ip , "%s" , tt ) ;
      fscanf( ip , "%s" , tt ) ;
      if (a < 10000 )
        fscanf( ip , "%d" , &di ) ;
      if (a % 100000 == 0)
        a = 0;


      for ( j=0 ; j<Dim ; j++ ) {
        fscanf( ip , "%lf" , &x[ind][j] ) ;
      }
      for ( j=0 ; j<Dim ; j++ ) {
        fscanf( ip , "%lf" , &x_bac[ind][j] ) ;
      }
      for ( j=0 ; j<Dim ; j++ ) {
        fscanf( ip , "%lf" , &gn_bac[ind][j] ) ;
      }


      fgets( tt, 80 , ip ) ;

      tp[ ind ] = 1 ;

      ind++ ;
    }
  }

  for ( k=0 ; k<nP ; k++ ) {

    a++;
    fscanf( ip , "%5d" , &di ) ;
    fscanf( ip , "%s" , tt ) ;
    fscanf( ip , "%s" , tt ) ;
    if (a < 10000 )
      fscanf( ip , "%d" , &di ) ;
    if (a % 100000 == 0)
      a = 0;


    for ( j=0 ; j<Dim ; j++ ) {
      fscanf( ip , "%lf" , &x[ind][j] ) ;
    }
    for ( j=0 ; j<Dim ; j++ ) {
      fscanf( ip , "%lf" , &x_bac[ind][j] ) ;
    }
    for ( j=0 ; j<Dim ; j++ ) {
      fscanf( ip , "%lf" , &gn_bac[ind][j] ) ;
    }


    fgets( tt, 80 , ip ) ;

    tp[ ind ] = 2 ;



    ind++ ;


    for(i=0; i<(ng_per_partic); i++){
      for(m = 0 ; m<(Ng+1); m++){
        a++;
        fscanf( ip , "%5d" , &di ) ;
        fscanf( ip , "%s" , tt ) ;
        fscanf( ip , "%s" , tt ) ;
        if (a < 10000 )
          fscanf( ip , "%d" , &di ) ;
        if (a % 100000 == 0)
          a = 0;


        for ( j=0 ; j<Dim ; j++ ) {
          fscanf( ip , "%lf" , &x[ind][j] ) ;
        }

        for ( j=0 ; j<Dim ; j++ ) {
          fscanf( ip , "%lf" , &x_bac[ind][j] ) ;
        }
        for ( j=0 ; j<Dim ; j++ ) {
          fscanf( ip , "%lf" , &gn_bac[ind][j] ) ;
        }

        fgets( tt, 80 , ip ) ;

        tp[ ind ] = (m == 0 ? -1 :0) ;

        ind++ ;

      }//Ng+1
    }//ng_p_np

  }//nP

  for (i=0; i<nstot; i++) {
    for (j=0; j<Dim; j++) {
      x[i][j] *= 10.0;
      x_bac[i][j] *= 10.0;
    }
  }

  // Assign the labels //
  for ( i=0 ; i<nstot ; i++ ) {
    if ( tp[i] == 0 )
      xc[i] = "H" ;
    else if ( tp[i] == 1 )
      xc[i] = "He" ;
    else if ( tp[i] == 2 )
      xc[i] = "O" ;
    else if ( tp[i] == 3 )
      xc[i] = "S" ;
    else if ( tp[i] == 4 )
      xc[i] = "N" ;
    else if ( tp[i] == 5 )
      xc[i] = "Br" ;
    else if ( tp[i] == 6 )
      xc[i] = "C" ;
    else if ( tp[i] == 7 )
      xc[i] = "Na" ;
    else if ( tp[i] == 8 )
      xc[i] = "P" ;
    else if ( tp[i] == 9 )
      xc[i] = "Ca" ;
    else if ( tp[i] == -1 )
      xc[i] = "X" ;
  }

}



void read_quaternions(FILE *ip){

  int center_ind, gind , i, m, j, k, di , ind ;
  int sites_per_gnp = 1 + ng_per_partic * ( 1 + Ng );
  double dr[Dim],mdr;
  char tt[80] ;

  fgets( tt , 80 , ip ) ;
  fscanf( ip , "%d" , &di ) ;
  fgets( tt, 80 , ip ) ;
  if ( di != nP)
    die("Number of nP in input_q.gro does not match!\n");


  for(i=0; i<nP ;i++){

    fscanf( ip , "%d" , &ind ) ;

    for ( j=0 ; j<4 ; j++ ) {
      fscanf( ip , "%lf" , &euler_q[ind][j] ) ;
    }


    fgets( tt, 80 , ip ) ;


    //cout<<"i "<<i<<endl;

  }

  //iniitialize A and B matrixes
  calc_A();

  //initialize x_grft at bf frames
  for(i=0; i<nP;i++){

    center_ind = nD * ( Nda + Ndb ) + nA * Nha + nB * Nhb + sites_per_gnp * i ;

    for(m=0; m<ng_per_partic ; m++){

      ind = center_ind + m * ( Ng + 1 ) + 1 ;

      gind = ng_per_partic*i + m;

      mdr = pbc_mdr2(x[ind],x[center_ind],dr);

      for(j=0; j<Dim; j++){

        grf_bf_x[gind][j] = 0.0;

        for(k=0 ;k<Dim ; k++){
          grf_bf_x[gind][j] += euler_A[i][j][k]*dr[k];

        }//k

      }//j

    }//ng_per_partic

  }//nP

}

void write_rst_quaternions(){

  FILE *otp ;
  int i, j, center_ind,gind,ind, k, m, resind , sites_per_gnp = 1 + ng_per_partic * ( 1 + Ng );

  otp = fopen( "rst_quaternions.gro" , "w" ) ;

  fprintf( otp , "%lf = steps * delt, Writing frame %d\n" , double(step) * delt , step ) ;
  fprintf( otp , "%d np\n" , nP );
  for(i= 0; i<nP;i++){
    center_ind = nD * ( Nda + Ndb ) + nA * Nha + nB * Nhb + sites_per_gnp * i ;

    fprintf( otp ,"%d %lf %lf %lf %lf\n",i,euler_q[i][0],euler_q[i][1],euler_q[i][2],euler_q[i][3]);

  }
  fprintf( otp , "%d %d\n" , Ng,nP  );
  fclose( otp ) ;
}

void write_quaternions(){

  FILE *otp ;
  int i, j, center_ind,gind,ind, k, m, resind , sites_per_gnp = 1 + ng_per_partic * ( 1 + Ng );

  if ( step == 0 )
    otp = fopen( "quaternions.gro" , "w" ) ;
  else
    otp = fopen( "quaternions.gro" , "a" ) ;

  fprintf( otp , "%lf = steps * delt, Writing frame %d\n" , double(step) * delt , step ) ;
  fprintf( otp , "%d np\n" , nP );
  for(i= 0; i<nP;i++){
    center_ind = nD * ( Nda + Ndb ) + nA * Nha + nB * Nhb + sites_per_gnp * i ;

    fprintf( otp ,"%d %lf %lf %lf %lf\n",i,euler_q[i][0],euler_q[i][1],euler_q[i][2],euler_q[i][3]);

  }
  fprintf( otp , "%d %d\n" , Ng,nP  );
  fclose( otp ) ;
}

void write_np() {

  FILE *otp ;
  int gind,i, j, ind, k, m, resind , center_ind;
  int real_ind,sites_per_gnp = 1 + ng_per_partic * ( 1 + Ng );

  if ( step == 0 )
    otp = fopen( "traj_np.gro" , "w" ) ;
  else
    otp = fopen( "traj_np.gro" , "a" ) ;

  fprintf( otp , "%lf = steps * delt, Writing frame %d\n" , double(step) * delt , step ) ;
  fprintf( otp , "%d\n" , nP*sites_per_gnp);

  ind = 0 ;
  resind = 0 ;

  for ( k=0 ; k<nP ; k++ ) {
    fprintf( otp , "%5d" , (resind + 1) % 100000 ) ;
    fprintf( otp , "%-5s" , "GP" ) ;
    fprintf( otp , "%5s" , "O" ) ;

    fprintf( otp , "%5d" , (ind + 1) % 100000 ) ;

    center_ind = nD * ( Nda + Ndb ) + nA * Nha + nB * Nhb + sites_per_gnp * k ;
    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%8.3lf" , x[center_ind][j] / 10.0 ) ;

    for ( j=Dim ; j<3 ; j++ )
      fprintf( otp , "%8.3lf" , 0.0 );

    fprintf( otp , "\n" ) ;

    real_ind  = center_ind +1;
    ind++ ;

    for ( i=0 ; i<ng_per_partic * ( Ng+1 ) ; i++ ){
      fprintf( otp , "%5d" , (resind + 1) % 100000 ) ;
      fprintf( otp , "%-5s" , "GP" ) ;
      fprintf( otp , "%5s" , xc[real_ind] ) ;

      fprintf( otp , "%5d" , (ind + 1) % 100000 ) ;

      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.3lf" , x[real_ind][j] / 10.0 ) ;
      for ( j=Dim ; j<3 ; j++ )
        fprintf( otp , "%8.3lf" , 0.0 );

      fprintf( otp , "\n" ) ;
      ind++ ;
      real_ind ++;
    }//ng_per_partic



    resind++ ;
  }//nP

  for ( j=0 ; j<Dim ; j++ )
    fprintf( otp , "%lf " , L[j] / 10.0 ) ;
  for ( j=Dim ; j<3 ; j++ )
    fprintf( otp , "1.0" ) ;
  fprintf( otp , "\n" ) ;


  fclose( otp ) ;

}

void write_rst_gro() {

  FILE *otp ;
  int i, j, ind, k, m, resind ;

  otp = fopen( "rst_coord.gro" , "w" ) ;

  fprintf( otp , "%lf = steps * delt, Writing frame %d\n" , double(step) * delt , step ) ;
  fprintf( otp , "%d\n" , nstot );

  ind = 0 ;
  resind = 0 ;
  for ( k=0 ; k<nD ; k++ ) {
    for ( m=0 ; m<Nda + Ndb ; m++ ) {
      fprintf( otp , "%5d" , (resind + 1) % 100000 ) ;
      fprintf( otp , "%-5s" , "BCP" ) ;
      fprintf( otp , "%5s" , xc[ind] ) ;

      fprintf( otp , "%5d" , (ind + 1) % 100000 ) ;

      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.4lf" , x[ind][j] / 10.0 ) ;

      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.4lf" ,  x_bac[ind][j] / 10.0 );

      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.4lf" ,  gn_bac[ind][j]  );


      fprintf( otp , "\n" ) ;

      ind++;
    }
    resind++ ;
  }

  for ( k=0 ; k<nA ; k++ ) {
    for ( m=0 ; m<Nha ; m++ ) {
      fprintf( otp , "%5d" , (resind + 1) % 100000 ) ;
      fprintf( otp , "%-5s" , "HA" ) ;
      fprintf( otp , "%5s" , xc[ind] ) ;

      fprintf( otp , "%5d" , (ind + 1) % 100000 ) ;

      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.3lf" , x[ind][j] / 10.0 ) ;

      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.4lf" ,  x_bac[ind][j] / 10.0 );

      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.4lf" ,  gn_bac[ind][j]  );


      fprintf( otp , "\n" ) ;

      ind++;
    }
    resind++ ;
  }

  for ( k=0 ; k<nB ; k++ ) {
    for ( m=0 ; m<Nhb ; m++ ) {
      fprintf( otp , "%5d" , (resind + 1) % 100000 ) ;
      fprintf( otp , "%-5s" , "HB" ) ;
      fprintf( otp , "%5s" , xc[ind] ) ;

      fprintf( otp , "%5d" , (ind + 1) % 100000 ) ;

      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.3lf" , x[ind][j] / 10.0 ) ;
      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.4lf" ,  x_bac[ind][j] / 10.0 );

      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.4lf" ,  gn_bac[ind][j]  );



      fprintf( otp , "\n" ) ;

      ind++;
    }
    resind++ ;
  }

  for ( k=0 ; k<nP ; k++ ) {
    fprintf( otp , "%5d" , (resind + 1) % 100000 ) ;
    fprintf( otp , "%-5s" , "GP" ) ;
    fprintf( otp , "%5s" , xc[ind] ) ;

    fprintf( otp , "%5d" , (ind + 1) % 100000 ) ;

    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%8.3lf" , x[ind][j] / 10.0 ) ;
    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%8.4lf" ,  x_bac[ind][j] / 10.0 );

    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%8.4lf" ,  gn_bac[ind][j]  );



    fprintf( otp , "\n" ) ;

    ind++ ;

    for ( i=0 ; i<ng_per_partic * ( Ng+1 ) ; i++ ) {
      fprintf( otp , "%5d" , (resind + 1) % 100000 ) ;
      fprintf( otp , "%-5s" , "GP" ) ;
      fprintf( otp , "%5s" , xc[ind] ) ;

      fprintf( otp , "%5d" , (ind + 1) % 100000 ) ;

      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.3lf" , x[ind][j] / 10.0 ) ;
      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.4lf" ,  x_bac[ind][j] / 10.0 );

      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.4lf" ,  gn_bac[ind][j]  );



      fprintf( otp , "\n" ) ;
      ind++ ;
    }


    resind++ ;
  }

  for ( j=0 ; j<Dim ; j++ )
    fprintf( otp , "%lf " , L[j] / 10.0 ) ;
  for ( j=Dim ; j<3 ; j++ )
    fprintf( otp , "1.0" ) ;
  fprintf( otp , "\n" ) ;


  fclose( otp ) ;

}


void write_gro() {

  FILE *otp ;
  int i, j, ind, k, m, resind ;
  if ( step == 0 )
    otp = fopen( "output.gro" , "w" ) ;
  else
    otp = fopen( "output.gro" , "a" ) ;

  fprintf( otp , "%lf = steps * delt, Writing frame %d\n" , double(step) * delt , step ) ;
  fprintf( otp , "%d\n" , nstot );

  ind = 0 ;
  resind = 0 ;
  for ( k=0 ; k<nD ; k++ ) {
    for ( m=0 ; m<Nda + Ndb ; m++ ) {
      fprintf( otp , "%5d" , (resind + 1) % 100000 ) ;
      fprintf( otp , "%-5s" , "BCP" ) ;
      fprintf( otp , "%5s" , xc[ind] ) ;

      fprintf( otp , "%5d" , (ind + 1) % 100000 ) ;

      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.3lf" , x[ind][j] / 10.0 ) ;

      for ( j=Dim ; j<3 ; j++ )
        fprintf( otp , "%8.3lf" , 0.0 );

      fprintf( otp , "\n" ) ;

      ind++;
    }
    resind++ ;
  }

  for ( k=0 ; k<nA ; k++ ) {
    for ( m=0 ; m<Nha ; m++ ) {
      fprintf( otp , "%5d" , (resind + 1) % 100000 ) ;
      fprintf( otp , "%-5s" , "HA" ) ;
      fprintf( otp , "%5s" , xc[ind] ) ;

      fprintf( otp , "%5d" , (ind + 1) % 100000 ) ;

      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.3lf" , x[ind][j] / 10.0 ) ;

      for ( j=Dim ; j<3 ; j++ )
        fprintf( otp , "%8.3lf" , 0.0 );

      fprintf( otp , "\n" ) ;

      ind++;
    }
    resind++ ;
  }

  for ( k=0 ; k<nB ; k++ ) {
    for ( m=0 ; m<Nhb ; m++ ) {
      fprintf( otp , "%5d" , (resind + 1) % 100000 ) ;
      fprintf( otp , "%-5s" , "HB" ) ;
      fprintf( otp , "%5s" , xc[ind] ) ;

      fprintf( otp , "%5d" , (ind + 1) % 100000 ) ;

      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.3lf" , x[ind][j] / 10.0 ) ;

      for ( j=Dim ; j<3 ; j++ )
        fprintf( otp , "%8.3lf" , 0.0 );

      fprintf( otp , "\n" ) ;

      ind++;
    }
    resind++ ;
  }

  for ( k=0 ; k<nP ; k++ ) {
    fprintf( otp , "%5d" , (resind + 1) % 100000 ) ;
    fprintf( otp , "%-5s" , "GP" ) ;
    fprintf( otp , "%5s" , xc[ind] ) ;

    fprintf( otp , "%5d" , (ind + 1) % 100000 ) ;

    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%8.3lf" , x[ind][j] / 10.0 ) ;

    for ( j=Dim ; j<3 ; j++ )
      fprintf( otp , "%8.3lf" , 0.0 );

    fprintf( otp , "\n" ) ;

    ind++ ;

    for ( i=0 ; i<ng_per_partic * ( Ng+1 ) ; i++ ) {
      fprintf( otp , "%5d" , (resind + 1) % 100000 ) ;
      fprintf( otp , "%-5s" , "GP" ) ;
      fprintf( otp , "%5s" , xc[ind] ) ;

      fprintf( otp , "%5d" , (ind + 1) % 100000 ) ;

      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , "%8.3lf" , x[ind][j] / 10.0 ) ;

      for ( j=Dim ; j<3 ; j++ )
        fprintf( otp , "%8.3lf" , 0.0 );

      fprintf( otp , "\n" ) ;
      ind++ ;
    }


    resind++ ;
  }

  for ( j=0 ; j<Dim ; j++ )
    fprintf( otp , "%lf " , L[j] / 10.0 ) ;
  for ( j=Dim ; j<3 ; j++ )
    fprintf( otp , "1.0" ) ;
  fprintf( otp , "\n" ) ;


  fclose( otp ) ;

}

void random_config( void ) {

  int i, m, j, k, ind = 0 ;
  int tmp_dr[Dim];
  double sqrt_nP;

  sqrt_nP = sqrt(nP);

  for ( k=0 ; k<nD ; k++ ) {
    x[ind][0] = ran2() * L[0];
    x[ind][1] = ran2() * L[1];
    x[ind][2] = ran2() * L[2] ;
    
    if ( x[ind][0] > L[0] )
      x[ind][0] -= L[0] ;
    else if ( x[ind][0] < 0.0 )
      x[ind][0] += L[0] ;

    if ( x[ind][1] > L[1] )
      x[ind][1] -= L[1] ;
    else if ( x[ind][1] < 0.0 )
      x[ind][1] += L[1] ;
    
    if ( x[ind][2] > L[2] )
      x[ind][2] -= L[2] ;
    else if ( x[ind][2] < 0.0 )
      x[ind][2] += L[2] ;

    tp[ ind ] = ( Nda > 0 ? 0 : 1 ) ;

    ind += 1 ;

    for ( m=1 ; m<Nda + Ndb ; m++ ) {

      for ( j=0 ; j<Dim ; j++ ) {
        x[ind][j] = x[ ind-1 ][ j ] + gasdev2() ;

        if ( x[ind][j] > L[j] )
          x[ind][j] -= L[j] ;
        else if ( x[ind][j] < 0.0 )
          x[ind][j] += L[j] ;
      }

      tp[ ind ] = ( m < Nda ? 0 : 1 ) ;

      ind++ ;

    }
  }

  /*
  // Random A homopolymers //
  for ( k=0 ; k<nA ; k++ ) {
  for ( j=0 ; j<Dim ; j++ )
  x[ind][j] = ran2() * L[j] ;

  tp[ ind ] = 0 ;

  ind += 1 ;

  for ( m=1 ; m<Nha ; m++ ) {

  for ( j=0 ; j<Dim ; j++ ) {
  x[ind][j] = x[ ind-1 ][ j ] + gasdev2() ;

  if ( x[ind][j] > L[j] )
  x[ind][j] -= L[j] ;
  else if ( x[ind][j] < 0.0 )
  x[ind][j] += L[j] ;
  }

  tp[ ind ] = 0 ;

  ind++ ;

  }
  }

  // Random B homopolymers //
  for ( k=0 ; k<nB ; k++ ) {
  for ( j=0 ; j<Dim ; j++ )
  x[ind][j] = ran2() * L[j] ;

  tp[ ind ] = 1 ;

  ind += 1 ;

  for ( m=1 ; m<Nhb ; m++ ) {

  for ( j=0 ; j<Dim ; j++ ) {
  x[ind][j] = x[ ind-1 ][ j ] + gasdev2() ;

  if ( x[ind][j] > L[j] )
  x[ind][j] -= L[j] ;
  else if ( x[ind][j] < 0.0 )
  x[ind][j] += L[j] ;
  }

  tp[ ind ] = 1 ;

  ind++ ;

  }
  }
  */


  for ( k=0 ; k<nP ; k++ ) {
    int center_ind , gind,prev_graft ;
    double u[Dim] ;


    double center_box_x = agg_frac * (1 - phi_M);
    x[ind][0] = L[0] * ran2();
    x[ind][1] = L[1] * ran2();
    if (double(k) / double(nP) < agg_frac) {
      // place aggregated particles in center of box on z axis
      x[ind][2] = L[2] * (0.5 + center_box_x * (ran2() - 0.5));
    } else {
      // place non-aggregated particles away from center of box on z axis
      x[ind][2] = L[2] * (1.0 - center_box_x) * (ran2() - 0.5);
    }

    if ( x[ind][0] > L[0] )
      x[ind][0] -= L[0] ;
    else if ( x[ind][j] < 0.0 )
      x[ind][0] += L[0] ;

    if ( x[ind][1] > L[1] )
      x[ind][1] -= L[1] ;
    else if ( x[ind][j] < 0.0 )
      x[ind][1] += L[1] ;
    
    if ( x[ind][2] > L[2] )
      x[ind][2] -= L[2] ;
    else if ( x[ind][j] < 0.0 )
      x[ind][2] += L[2] ;

    if(Dim ==3 ) {  for(j=0 ; j<3 ; j++){

      euler_ang[k][j] = 0;
      euler_q[k][j+1] = 0;
    }

    euler_q[k][0] = 1;
    // euler_ang[k][0] = PI/2.0;
    }
    tp[ ind ] = 2 ;

    center_ind = ind ;

    ind += 1 ;


    for ( m=0 ; m<(ng_per_partic) ; m++ ) {

      // Place the grafting site //
      if(Dim == 3 and uni_sig ==1){ //fibonacci grafting sites
        fibonacci_u ( u , m );
        gind = k*(ng_per_partic) + m;
        if(k==0 and m==0)cout<<"fibonacci grafting is used"<<endl;
      }
      else if( Dim ==3 and  uni_sig ==2){//uniform grafting
        unif_sig ( u , m );
        gind = k*(ng_per_partic) + m;
        if(k==0  and m==0)cout<<"uniform grafting is used"<<endl;
      }
      else{ //random grafting sites
        random_u( u ) ;
        if(Dim == 3) gind = k*(ng_per_partic) + m;
        if(k==0 and m==0)cout<<"random grafting is used"<<endl;
      }


      for ( j=0 ; j<Dim ; j++ ) {
        x[ind][j] = x[center_ind][j] + (Rp+Xi) * u[j] ;
        if(Dim == 3)  grf_bf_x[gind][j] = (Rp+Xi) * u[j];

        if ( x[ind][j] > L[j] )
          x[ind][j] -= L[j] ;
        else if ( x[ind][j] < 0.0 )
          x[ind][j] += L[j] ;
      }
      //  cout<<" ng "<<m<<" "<<x[ind][0]<<" "<<x[ind][1]<<" "<<x[ind][2]<<endl;

      tp[ ind ] = -1 ;

      if ( Dim == 3 and m > 0 ) {
        double mdr2, dr[Dim], mdr ;
        mdr2 = pbc_mdr2( x[ prev_graft ] , x[ ind ] , dr ) ;
        mdr = sqrt( mdr2 ) ;

        graft_req[k][m-1] = mdr ;
      }

      prev_graft = ind ;

      ind++ ;

      // Place the chain grafted to this site
      if(m <ng_per_partic){
        for ( i=0 ; i<Ng ; i++ ) {

          for ( j=0 ; j<Dim ; j++ ) {
            x[ind][j] = x[ind-1][j] + 0.5 * u[j] ;

            if ( x[ind][j] > L[j] )
              x[ind][j] -= L[j] ;
            else if ( x[ind][j] < 0.0 )
              x[ind][j] += L[j] ;
          }

          tp[ ind ] = 0 ;
          ind++ ;
        }
      }//m<ng_per

    }//m
  }//nP


  // Assign the labels //
  for ( i=0 ; i<nstot ; i++ ) {
    //cout << "i: " << i << "tp[i]: " << tp[i] << endl;
    if ( tp[i] == 0 )
      xc[i] = "H" ;
    else if ( tp[i] == 1 )
      xc[i] = "He" ;
    else if ( tp[i] == 2 )
      xc[i] = "O" ;
    else if ( tp[i] == 3 )
      xc[i] = "S" ;
    else if ( tp[i] == 4 )
      xc[i] = "N" ;
    else if ( tp[i] == 5 )
      xc[i] = "Br" ;
    else if ( tp[i] == 6 )
      xc[i] = "C" ;
    else if ( tp[i] == 7 )
      xc[i] = "Na" ;
    else if ( tp[i] == 8 )
      xc[i] = "P" ;
    else if ( tp[i] == 9 )
      xc[i] = "Ca" ;
    else if ( tp[i] == -1 )
      xc[i] = "X" ;
  }

  calc_A();

  printf("config generated!\n") ; fflush( stdout ) ;
}

