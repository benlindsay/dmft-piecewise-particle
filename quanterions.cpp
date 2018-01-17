#include "globals.h"



void  calc_dAdang();




void calc_A(){

	int i ,j,k,l;
	double tmp_qt[4],tmp_4norm;
#pragma omp parallel for private(j,k,tmp_4norm)
  for(i =0 ; i<nP ;i++){
/*  	   euler_q[i][0] = cos(euler_ang[i][1]/2.0)*cos( 0.5*(euler_ang[i][0]+euler_ang[i][2]) );
	     euler_q[i][1] = sin(euler_ang[i][1]/2.0)*cos( 0.5*(euler_ang[i][0]-euler_ang[i][2]) );
	  euler_q[i][2] = sin(euler_ang[i][1]/2.0)*sin( 0.5*(euler_ang[i][0]-euler_ang[i][2]) );
	  euler_q[i][3] = cos(euler_ang[i][1]/2.0)*sin( 0.5*(euler_ang[i][0]+euler_ang[i][2]) );
*/	
	for(j=0 ; j<3; j++){
		euler_A[i][j][j] =  euler_q[i][0]*euler_q[i][0] ;
		for(k=1 ; k<4;k++)
		  euler_A[i][j][j] += (k ==(j+1) ? 1 :-1) * euler_q[i][k]*euler_q[i][k] ;
	}//j=j
	
	euler_A[i][0][1] = 2*(euler_q[i][1]*euler_q[i][2] +euler_q[i][0]*euler_q[i][3]);
	euler_A[i][0][2] = 2*(euler_q[i][1]*euler_q[i][3] -euler_q[i][0]*euler_q[i][2]);
	euler_A[i][1][0] = 2*(euler_q[i][1]*euler_q[i][2] -euler_q[i][0]*euler_q[i][3]);
  	euler_A[i][2][0] = 2*(euler_q[i][1]*euler_q[i][3] +euler_q[i][0]*euler_q[i][2]);
  	euler_A[i][2][1] = 2*(euler_q[i][2]*euler_q[i][3] - euler_q[i][0]*euler_q[i][1]);
  	euler_A[i][1][2] = 2*(euler_q[i][2]*euler_q[i][3] +euler_q[i][0]*euler_q[i][1]);

      tmp_4norm =0 ;

      for(j=0;j<4;j++)
      	 tmp_4norm += euler_q[i][j]*euler_q[i][j];
      
      tmp_4norm = tmp_4norm*tmp_4norm;


      for(j=0 ;j<4 ;j++){
        for(k=0;k<4;k++){
	 if(j ==k)
	 	euler_B[i][j][k] = euler_q[i][0];
	 else if( (j+k)==3)
	  	euler_B[i][j][k] = (j<k ? -1:1.0)*euler_q[i][3];
	 else if ((j+k) ==2)
	 	euler_B[i][j][k] = (j<k ? -1:1.0)*euler_q[i][2];
	 else if ((j+k) ==4)
	 	euler_B[i][j][k] = (j<k ?  1: -1.0)*euler_q[i][2];
	 else if ((j+k) ==1)
	 	euler_B[i][j][k] = (j<k ?  -1 : 1.0)*euler_q[i][1];
	 else{
		euler_B[i][j][k] = (j<k ?  -1 : 1.0)*euler_q[i][1];
	 }
	 euler_B[i][j][k] *= 0.5/tmp_4norm;
	//   cout<<euler_B[0][j][k]<<" ";
	}//k
//	cout<<endl;
      }//j
   }//nP


}



void calc_dAdang(){

  int i,j,k,l ;
   double tmp_qdot[4];


//calcu dAdphi //

 #pragma omp parallel for private(j,k)
  for(i =0 ; i<nP ;i++){
  	  tmp_qdot[0] = -0.5*cos(euler_ang[i][1]/2.0)*sin( 0.5*(euler_ang[i][0]+euler_ang[i][2]) );
	  tmp_qdot[1] = -0.5*sin(euler_ang[i][1]/2.0)*sin( 0.5*(euler_ang[i][0]-euler_ang[i][2]) );
	  tmp_qdot[2] = 0.5*sin(euler_ang[i][1]/2.0)*cos( 0.5*(euler_ang[i][0]-euler_ang[i][2]) );
	  tmp_qdot[3] = 0.5*cos(euler_ang[i][1]/2.0)*cos( 0.5*(euler_ang[i][0]+euler_ang[i][2]) );

	for(j=0 ; j<3; j++){
		dAdphi[i][j][j] = 2.0* euler_q[i][0]*tmp_qdot[0] ;
		for(k=1 ; k<4;k++)
		  dAdphi[i][j][j] += (k ==(j+1) ? 1 :-1) * 2.0*euler_q[i][k]*tmp_qdot[k];
	}//j=j
	
	dAdphi[i][0][1] = 2*(euler_q[i][1]*tmp_qdot[2]+tmp_qdot[1]*euler_q[i][2] + euler_q[i][0]*tmp_qdot[3] + tmp_qdot[0]*euler_q[i][3]);

	dAdphi[i][0][2] = 2*(euler_q[i][1]*tmp_qdot[3]+tmp_qdot[1]*euler_q[i][3] - tmp_qdot[0]*euler_q[i][2]- euler_q[i][0]*tmp_qdot[2]);
	
	dAdphi[i][1][0] = 2*(euler_q[i][1]*tmp_qdot[2]+euler_q[i][2]*tmp_qdot[1] -euler_q[i][0]*tmp_qdot[3]  - tmp_qdot[0]*euler_q[i][3]);
  	
	dAdphi[i][2][0] = 2*(euler_q[i][1]*tmp_qdot[3]+tmp_qdot[1]*euler_q[i][3] +euler_q[i][0]*tmp_qdot[2]+ tmp_qdot[0]*euler_q[i][2]);

	dAdphi[i][2][1] = 2*(euler_q[i][2]*tmp_qdot[3]+tmp_qdot[2]*euler_q[i][3] - euler_q[i][0]*tmp_qdot[1] -tmp_qdot[0]*euler_q[i][1]);
  	
	dAdphi[i][1][2] = 2*(euler_q[i][2]*tmp_qdot[3]+tmp_qdot[2]*euler_q[i][3] +euler_q[i][0]*tmp_qdot[1] +tmp_qdot[0]*euler_q[i][1]);
  }//nP

 
//calcu dAdtheta //

 #pragma omp parallel for private(j,k)
  for(i =0 ; i<nP ;i++){
  	  tmp_qdot[0] = -0.5*sin(euler_ang[i][1]/2.0)*cos( 0.5*(euler_ang[i][0]+euler_ang[i][2]) );
	  tmp_qdot[1] = 0.5*cos(euler_ang[i][1]/2.0)*cos( 0.5*(euler_ang[i][0]-euler_ang[i][2]) );
	  tmp_qdot[2] = 0.5*cos(euler_ang[i][1]/2.0)*sin( 0.5*(euler_ang[i][0]-euler_ang[i][2]) );
	  tmp_qdot[3] = -0.5*sin(euler_ang[i][1]/2.0)*sin( 0.5*(euler_ang[i][0]+euler_ang[i][2]) );

	for(j=0 ; j<3; j++){
		dAdtheta[i][j][j] = 2.0* euler_q[i][0]*tmp_qdot[0] ;
		for(k=1 ; k<4;k++)
		  dAdtheta[i][j][j] += (k ==(j+1) ? 1 :-1) * 2.0*euler_q[i][k]*tmp_qdot[k];
	}//j=j
	
	dAdtheta[i][0][1] = 2*(euler_q[i][1]*tmp_qdot[2]+tmp_qdot[1]*euler_q[i][2] + euler_q[i][0]*tmp_qdot[3] + tmp_qdot[0]*euler_q[i][3]);

	dAdtheta[i][0][2] = 2*(euler_q[i][1]*tmp_qdot[3]+tmp_qdot[1]*euler_q[i][3] - tmp_qdot[0]*euler_q[i][2]- euler_q[i][0]*tmp_qdot[2]);
	
	dAdtheta[i][1][0] = 2*(euler_q[i][1]*tmp_qdot[2]+euler_q[i][2]*tmp_qdot[1] -euler_q[i][0]*tmp_qdot[3]  - tmp_qdot[0]*euler_q[i][3]);
  	
	dAdtheta[i][2][0] = 2*(euler_q[i][1]*tmp_qdot[3]+tmp_qdot[1]*euler_q[i][3] +euler_q[i][0]*tmp_qdot[2]+ tmp_qdot[0]*euler_q[i][2]);

	dAdtheta[i][2][1] = 2*(euler_q[i][2]*tmp_qdot[3]+tmp_qdot[2]*euler_q[i][3] - euler_q[i][0]*tmp_qdot[1] -tmp_qdot[0]*euler_q[i][1]);
  	
	dAdtheta[i][1][2] = 2*(euler_q[i][2]*tmp_qdot[3]+tmp_qdot[2]*euler_q[i][3] +euler_q[i][0]*tmp_qdot[1] +tmp_qdot[0]*euler_q[i][1]);
  
 
  }//nP

//calcu dAdpsi //

 #pragma omp parallel for private(l,j,k)
  for(i =0 ; i<nP ;i++){
  	  tmp_qdot[0] = -0.5*cos(euler_ang[i][1]/2.0)*sin( 0.5*(euler_ang[i][0]+euler_ang[i][2]) );
	  tmp_qdot[1] = 0.5*sin(euler_ang[i][1]/2.0)*sin( 0.5*(euler_ang[i][0]-euler_ang[i][2]) );
	  tmp_qdot[2] = -0.5*sin(euler_ang[i][1]/2.0)*cos( 0.5*(euler_ang[i][0]-euler_ang[i][2]) );
	  tmp_qdot[3] = 0.5*cos(euler_ang[i][1]/2.0)*cos( 0.5*(euler_ang[i][0]+euler_ang[i][2]) );

	for(j=0 ; j<3; j++){
		dAdpsi[i][j][j] = 2.0* euler_q[i][0]*tmp_qdot[0] ;
		for(k=1 ; k<4;k++)
		  dAdpsi[i][j][j] += (k ==(j+1) ? 1 :-1) * 2.0*euler_q[i][k]*tmp_qdot[k];
	}//j=j
	
	dAdpsi[i][0][1] = 2*(euler_q[i][1]*tmp_qdot[2]+tmp_qdot[1]*euler_q[i][2] + euler_q[i][0]*tmp_qdot[3] + tmp_qdot[0]*euler_q[i][3]);

	dAdpsi[i][0][2] = 2*(euler_q[i][1]*tmp_qdot[3]+tmp_qdot[1]*euler_q[i][3] - tmp_qdot[0]*euler_q[i][2]- euler_q[i][0]*tmp_qdot[2]);
	
	dAdpsi[i][1][0] = 2*(euler_q[i][1]*tmp_qdot[2]+euler_q[i][2]*tmp_qdot[1] -euler_q[i][0]*tmp_qdot[3]  - tmp_qdot[0]*euler_q[i][3]);
  	
	dAdpsi[i][2][0] = 2*(euler_q[i][1]*tmp_qdot[3]+tmp_qdot[1]*euler_q[i][3] +euler_q[i][0]*tmp_qdot[2]+ tmp_qdot[0]*euler_q[i][2]);

	dAdpsi[i][2][1] = 2*(euler_q[i][2]*tmp_qdot[3]+tmp_qdot[2]*euler_q[i][3] - euler_q[i][0]*tmp_qdot[1] -tmp_qdot[0]*euler_q[i][1]);
  	
	dAdpsi[i][1][2] = 2*(euler_q[i][2]*tmp_qdot[3]+tmp_qdot[2]*euler_q[i][3] +euler_q[i][0]*tmp_qdot[1] +tmp_qdot[0]*euler_q[i][1]);
 
       /*       if(i ==0 ){ 
       	for(j=0; j<3;j++){
 		for(k=0;k<3;k++){
			cout<<dAdpsi[i][j][k]<<" ";
		}
		cout<<endl;
	}
      }*/
 

  }//nP






}
