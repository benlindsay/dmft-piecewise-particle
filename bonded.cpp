#include "globals.h"
void gnp_bonds( void ) ;

void gnp_bonds( ) {

  int i, j, k, m, ind , center_ind ;
  double mf, mdr2, mdr, dr[Dim], delr ;

  int sites_per_gnp = 1 + ng_per_partic * ( 1 + Ng ), prev_graft_site ;

  kgraft = 20.0 ;
#pragma omp parallel for \
  private( center_ind, m, ind, mdr2, dr, j, k, prev_graft_site ) \
  reduction(+:Ubond)
  for ( i=0 ; i<nP ; i++ ) {

    center_ind = nD * ( Nda + Ndb ) + nA * Nha + nB * Nhb + sites_per_gnp * i ;

    for ( m=0 ; m<ng_per_partic ; m++ ) {
      ind = center_ind + m * ( Ng + 1 ) + 1 ;


      prev_graft_site = ind ;
      ind++ ;

      for ( k=0 ; k<Ng ; k++ ) {
        mdr2 = pbc_mdr2( x[ind] , x[ind-1] , dr ) ;
        Ubond += mdr2 * 1.5 ;

        for ( j=0 ; j<Dim ; j++ ) {
          f[ind][j] -= 3.0 * dr[j] ;
          f[ind-1][j] += 3.0 * dr[j] ;

          if( k==0 ){
            f[center_ind][j] += 3.0 * dr[j] ;
          }
        }

        ind++ ;
      }

    }// for ( m=0 ; m<ng_per_partic ;
  }//for ( i=0 ; i<nP
}

void bonds( ) {

  int i, j, m, ind ;
  double mdr2, dr[Dim] ;

  Ubond = 0.0 ;

  // Diblock bonds //
#pragma omp parallel for \
  private(ind, i, j, m, dr, mdr2) \
  reduction(+:Ubond)
  for ( i=0 ; i<nD ; i++ ) {
    for ( m=0 ; m<Nda + Ndb - 1 ; m++ ) {

      ind = i * (Nda + Ndb) + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;

      Ubond += mdr2 * 1.5 ;

      for ( j=0 ; j<Dim ; j++ ) {
        f[ind][j] -= 3.0 * dr[j] ;
        f[ind+1][j] += 3.0 * dr[j] ;
      }

    }

  } // for ( i=0 ; i<nT[k]


  // Homopolymer A bonds //
#pragma omp parallel for \
  private(ind, j, m, dr, mdr2) \
  reduction(+:Ubond)
  for ( i=0 ; i<nA ; i++ ) {
    for ( m=0 ; m<Nha - 1 ; m++ ) {

      ind = nD * (Nda + Ndb) + i * Nha + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;

      Ubond += mdr2 * 1.5 ;

      for ( j=0 ; j<Dim ; j++ ) {
        f[ind][j] -= 3.0 * dr[j] ;
        f[ind+1][j] += 3.0 * dr[j] ;
      }

    }

  } // for ( i=0 ; i<nT[k]

  // Homopolymer B bonds //
#pragma omp parallel for \
  private(ind, i, j, m, dr, mdr2) \
  reduction(+:Ubond)
  for ( i=0 ; i<nB ; i++ ) {
    for ( m=0 ; m<Nhb - 1 ; m++ ) {

      ind = nD * (Nda + Ndb) + nA * Nha + i * Nhb + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;

      Ubond += mdr2 * 1.5 ;

      for ( j=0 ; j<Dim ; j++ ) {
        f[ind][j] -= 3.0 * dr[j] ;
        f[ind+1][j] += 3.0 * dr[j] ;
      }

    }

  } // for ( i=0 ; i<nT[k]


  gnp_bonds() ;


}


void bond_stress() {
  int center_ind,i,k,j1, j2, m, ind ;
  double mdr,mdr2, dr[Dim] ;
  int sites_per_gnp = 1 + ng_per_partic * ( 1 + Ng ), prev_graft_site ;

  ind = 0 ;



  for (j1=0 ; j1<Dim ; j1++ )
    for ( j2=0 ; j2<Dim ; j2++ ){
      Stress_bonds[j1][j2] = 0.0 ;
      for(m=0; m<nthreads; m++)
        Stress_bond_t[j1][j2][m] = 0.0 ;
    }


  // Diblock bonds //
#pragma omp parallel for private(ind,m,mdr2,dr,j1,j2)
  for ( i=0 ; i<nD ; i++ ) {
    for ( m=0 ; m<Nda + Ndb - 1 ; m++ ) {

      ind = i * (Nda + Ndb) + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;

      int tid = omp_get_thread_num() ;
      for ( j1=0 ; j1<Dim ; j1++ )
        for ( j2=0 ; j2<Dim ; j2++ )
          Stress_bond_t[j1][j2][tid] +=  dr[j1] * dr[j2] ;

    }

  } // for ( i=0 ; i<nT[k]


  // Homopolymer A bonds //
#pragma omp parallel for private(ind,m,dr,mdr2,j1,j2)
  for ( i=0 ; i<nA ; i++ ) {
    for ( m=0 ; m<Nha - 1 ; m++ ) {

      ind = nD * (Nda + Ndb) + i * Nha + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;
      int tid = omp_get_thread_num() ;

      for ( j1=0 ; j1<Dim ; j1++ )
        for ( j2=0 ; j2<Dim ; j2++ )
          Stress_bond_t[j1][j2][tid] +=  dr[j1] * dr[j2] ;

    }

  } // for ( i=0 ; i<nT[k]

  // Homopolymer B bonds //
#pragma omp parallel for private(ind,m,dr,mdr2,j1,j2)
  for ( i=0 ; i<nB ; i++ ) {
    for ( m=0 ; m<Nhb - 1 ; m++ ) {

      ind = nD * (Nda + Ndb) + nA * Nha + i * Nhb + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;

      int tid = omp_get_thread_num() ;

      for ( j1=0 ; j1<Dim ; j1++ )
        for ( j2=0 ; j2<Dim ; j2++ )
          Stress_bond_t[j1][j2][tid] +=  dr[j1] * dr[j2] ;

    }

  } // for ( i=0 ; i<nT[k]

  // HomoA gaft chain bons //
#pragma omp parallel for \
  private( dr,center_ind, m, ind, mdr2, j1,j2, k, prev_graft_site )
  for ( i=0 ; i<nP ; i++ ) {

    center_ind = nD * ( Nda + Ndb ) + nA * Nha + nB * Nhb + sites_per_gnp * i ;

    for ( m=0 ; m<ng_per_partic ; m++ ) {
      ind = center_ind + m * ( Ng + 1 ) + 1 ;

      int tid = omp_get_thread_num() ;

      prev_graft_site = ind ;
      ind += 1 ;
      for ( k=0 ; k<Ng ; k++ ) {
        mdr2 = pbc_mdr2( x[ind-1] , x[ind] , dr ) ;

        for ( j1=0 ; j1<Dim ; j1++ )
          for ( j2=0 ; j2<Dim ; j2++ )
            Stress_bond_t[j1][j2][tid] +=  dr[j1] * dr[j2] ;

        ind++;
      }

    }//m

  }//for ( i=0 ; i<nP

  for ( j1=0 ; j1<Dim ; j1++ )
    for ( j2=0 ; j2<Dim ; j2++ )
      for ( m=0 ; m<nthreads ; m++ )
        Stress_bonds[j1][j2] += 3.0 * Stress_bond_t[j1][j2][m] /V;

}
