#define MAIN
#include "globals.h"
void update_positions( void ) ;
void update_positions_init( void ) ;
void update_euler(void);
void initialize( void ) ;
void write_np(void) ;
void write_gro( void ) ;
void write_rst_gro(void );
void write_quaternions(void) ;
void write_rst_quaternions(void) ;
void write_grid( void ) ;
void forces( void ) ;
void torque(void);
double integrate( double* ) ;
void write_stress( void ) ;
void bond_stress( void ) ;
void calc_Unb( void ) ;
void calc_stress();


int main( int argc , char** argv ) {

  int i,j,k,l;

  n_samples = 0;

  if ( argc == 3 && !strcmp( "-nt" , argv[1] ) ) {
    nthreads = atoi( argv[2] ) ;
    omp_set_num_threads( nthreads ) ;
    printf("\nNumber of threads set to %d!\n" , nthreads ) ;
  }
  else {
    nthreads = 1 ;
    omp_set_num_threads( nthreads ) ;
    printf("\nNumber of threads set to %d!\n" , nthreads ) ;
  }

  n_samples_iter = 0;
  initialize() ;

  write_gro( ) ;
  write_quaternions( ) ;
  // Save chi for pre-equilibration steps //
  double tmp_ang,chi_bkp = chiAB ;
  double eps_bkp = eps ;
  double delt_bkp = delt ;
  double kappa_bkp = kappa ;
  double kappa_p_bkp = kappa_p ;

  FILE *otp ;
  //otp = fopen( "data.dat" , "w" ) ;

  for (i=0; i<M; i++) {
    avg_rhoga_iter[i] = 0.0;
    avg_rhoda_iter[i] = 0.0;
    avg_rhodb_iter[i] = 0.0;
    avg_rhop_iter[i] = 0.0;
  }

  printf("Entering main loop!\n") ; fflush( stdout ) ;
  for ( step = 0 ; step < nsteps ; step++ )  {


    if ( step < sample_wait && eps_flag == 1 ) {
      delt = 0.0002 ;
      eps = 0.0 ;
      kappa = kappa_p = 10;
    }
    else {
      eps = eps_bkp ;
      delt = delt_bkp ;
      kappa = kappa_bkp ;
      kappa_p = kappa_p_bkp ;
    }

    forces() ;
    if(sigma>0){
      torque();
    }
    //  exit(1);

    if(step >0 || rst_para == 1)
      update_positions() ;
    else
      update_positions_init() ;

    if(sigma>0){
      update_euler();
    }

    /*
       if ( stress_freq > 0 && step % stress_freq == 0 ) {
       calc_stress() ;

       for ( j=0 ; j<Dim ; j++ ){
       for ( k=0 ; k<Dim ; k++ ) {
       sts_buf[buff_ind][j][k]= Rg3*Ptens[j][k];//( j<Dim ? Stress_bonds[j][k]:0.0) ;
       sts_buf_pp[buff_ind][j][k] = Rg3*Stress_PP[j][k];
       sts_buf_ng[buff_ind][j][k] = Rg3*Stress_Ng[j][k];

    //for ( k=0 ; k<nP;  k++ ){
    //euler_adot[k][j] = euler_q[k][j];
    //sts_buf[buff_ind][j][k+Dim] = euler_q[k][j];
    //}
    }
    }
    buff_ind++ ;
    }*/

    if(step % sample_freq == 0){
      write_np();
    }


    /*
       if ( step > sample_wait && step % sample_freq == 0 ) {
       fftw_fwd( rho[0] , ktmp ) ;
       for ( i=0 ; i<M ; i++ ) {
       avg_sk[0][i] += ktmp[i] * conj(ktmp[i]) ;
       }

       if ( nP > 0 ) {

       fftw_fwd( rho[2] , ktmp ) ;
       for ( i=0 ; i<M ; i++ ) {
       avg_sk[2][i] += ktmp[i] * conj( ktmp[i] ) ;
       }
       }

       num_averages += 1.0 ;
       }*/


    if ( step % print_freq == 0 || step == nsteps-1 ) {
      printf("step %d of %d  Ubond: %lf\n" , step , nsteps , Ubond ) ;
      fflush( stdout ) ;
      write_gro() ;
      write_rst_gro();
      write_quaternions();
      write_rst_quaternions();

      //if ( stress_freq > 0 )
      //  write_stress() ;

      if ( step > sample_wait && step % sample_freq == 0 ) {
        if (step % avg_freq == 0) {
          n_samples_iter = 0.0;
          avg_iter += 1;
          for (i=0; i<M; i++) {
            avg_rhoga_iter[i] = 0.0;
            avg_rhoda_iter[i] = 0.0;
            avg_rhodb_iter[i] = 0.0;
            avg_rhop_iter[i] = 0.0;
          }
        }
        n_samples += 1.0 ;
        n_samples_iter += 1.0 ;
        for ( i=0 ; i<M ; i++ ) {
          avg_rhoga[i] += rhoga[i] ;
          avg_rhoda[i] += rhoda[i] ;
          avg_rhodb[i] += rhodb[i] ;
          avg_rhop[i] += rhop[i] ;
          avg_rhoga_iter[i] += rhoga[i] ;
          avg_rhoda_iter[i] += rhoda[i] ;
          avg_rhodb_iter[i] += rhodb[i] ;
          avg_rhop_iter[i] += rhop[i] ;
        }
        if (ng_per_partic > 0.0) {
          write_avg_grid_data( "avg_rhoga.dat" , avg_rhoga ) ;
          write_avg_grid_data_iter( "avg_rhoga" , avg_rhoga_iter ) ;
        }
        if (Nda > 0.0) {
          write_avg_grid_data( "avg_rhoda.dat" , avg_rhoda ) ;
          write_avg_grid_data_iter( "avg_rhoda" , avg_rhoda_iter ) ;
        }
        if (Ndb > 0.0) {
          write_avg_grid_data( "avg_rhodb.dat" , avg_rhodb ) ;
          write_avg_grid_data_iter( "avg_rhodb" , avg_rhodb_iter ) ;
        }
        if (nP > 0.0) {
          write_avg_grid_data( "avg_rhop.dat" , avg_rhop ) ;
          write_avg_grid_data_iter( "avg_rhop" , avg_rhop_iter ) ;
        }
      }


      write_grid_data( "rhoda.dat" , rhoda ) ;
      write_grid_data( "rhodb.dat" , rhodb ) ;

      if ( nA > 0.0 )
        write_grid_data( "rhoha.dat" , rhoha ) ;

      if ( nB > 0.0 )
        write_grid_data( "rhohb.dat" , rhohb ) ;

      if ( nP > 0.0 )
        write_grid_data( "rhop.dat" , rhop ) ;

      /*
         if ( step > sample_wait ) {
         for ( i=0 ; i<M ; i++ )
         ktmp2[i] = avg_sk[0][i] / num_averages ;
         write_kspace_data( "avg_sk_A.dat" , ktmp2 ) ;

         if ( nP > 0 ) {
         for ( i=0 ; i<M ; i++ )
         ktmp2[i] = avg_sk[2][i] / num_averages ;
         write_kspace_data( "avg_sk_np.dat" , ktmp2 ) ;
         }
         }*/

      //calc_Unb() ;

      //fprintf( otp , "%d %lf %lf %lf %lf %lf %lf %lf\n" , step , Ubond , U_chi_gg, U_kappa_gg,U_chi_pg,U_kappa_pg,U_kappa_pp , Utt) ;
      //fflush( otp ) ;

    }// if step % print_Freq == 0

  }

  //fclose( otp ) ;

  return 0 ;

}
