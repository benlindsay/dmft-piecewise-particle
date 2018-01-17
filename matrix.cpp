#include "globals.h"

void matrix_vect(double** mat,double* vect,double* rlst_vect){

int i , j ;
   
   for(i=0;i<Dim;i++)
   	rlst_vect[i] = 0.0;


   for(i=0;i<Dim;i++)
      for(j=0;j<Dim;j++){
	 rlst_vect[i] += mat[i][j]*vect[j];
      
      }


}

void matrix_trs_vect(double** mat,double* vect,double* rlst_vect){

int i , j ;
   
   for(i=0;i<Dim;i++)
   	rlst_vect[i] = 0.0;


   for(i=0;i<Dim;i++)
      for(j=0;j<Dim;j++){
	 rlst_vect[i] += mat[j][i]*vect[j];
      
      }


}

void fibonacci_u(double *u, int m){


    double offset = 2./double(ng_per_partic) ; 

    double incrmt = PI*(3.0 - sqrt(5.0));

    u[1] = (m)*offset -1.0 +offset/2.0 ;

    double r = sqrt(1.0-pow(u[1],2.0));

    double phi = double((m+1)%ng_per_partic)*incrmt ; 

    u[0] = cos(phi)*r;
   
    u[2] = sin(phi)*r;

//    cout<<u[0]<<" "<<u[1]<<" "<<u[2]<<endl;


}

void unif_sig (double *u, int m){

   if(ng_per_partic ==4){

	if(m==0){
		u[0] = 0;
		u[1] = 0;
		u[2] = -1;
	}
	else if (m==1){
		u[0] = 2.0*sqrt(2)/3.0;
		u[1] = 0;
		u[2] = 1.0/3.0;
	}
	else if(m==2){
		u[0] = -sqrt(2.0)/3.0;
		u[1] = sqrt(2.0/3.0);
		u[2] = 1.0/3.0;
	}
	else {

		u[0] = -sqrt(2.0)/3.0;
		u[1] = -sqrt(2.0/3.0);
		u[2] = 1.0/3.0;
	}
   }
   
   if(ng_per_partic ==6){

	if(m==0){
		u[0] = 1;
		u[1] = 0;
		u[2] = 0;
	}
	else if (m==1){
		u[0] = 0;
		u[1] = 1;
		u[2] = 0;
		
	}
	else if(m==2){
		u[0] = -1;
		u[1] = 0;
		u[2] = 0;	
	}
	else if(m==3){

		u[0] = 0;//-sqrt(2.0)/3.0;
		u[1] = -1;//-sqrt(2.0/3.0);
		u[2] = 0;//1.0/3.0;
	}
	else if(m==4){
		u[0] = 0;
		u[1] = 0;
		u[2] = 1;

	}
	else {
	       	u[0] = 0;
		u[1] = 0;
		u[2] = -1;
	}
   }
 
   if(ng_per_partic ==8){

	if(m==0){
		u[0] = sqrt(1.0/3.0);
		u[1] = sqrt(1.0/3.0);
		u[2] = -sqrt(1.0/3.0);
	}
	else if (m==1){
		u[0] = -sqrt(1.0/3.0);
		u[1] = sqrt(1.0/3.0);
		u[2] = -sqrt(1.0/3.0);
		
	}
	else if(m==2){
		u[0] = -sqrt(1.0/3.0);
		u[1] = -sqrt(1.0/3.0);
		u[2] = -sqrt(1.0/3.0);	
	}
	else if(m==3){

		u[0] = sqrt(1.0/3.0);//-sqrt(2.0)/3.0;
		u[1] = -sqrt(1.0/3.0);//-sqrt(2.0/3.0);
		u[2] = -sqrt(1.0/3.0);//1.0/3.0;
	}
	else if(m==4){
		u[0] = sqrt(1.0/3.0);
		u[1] = sqrt(1.0/3.0);
		u[2] = sqrt(1.0/3.0);

	}
	else if(m==5){
	       	u[0] = -sqrt(1.0/3.0);
		u[1] = sqrt(1.0/3.0);
		u[2] = sqrt(1.0/3.0);
	}
	else if(m==6){
	       	u[0] = -sqrt(1.0/3.0);
		u[1] = -sqrt(1.0/3.0);
		u[2] = sqrt(1.0/3.0);
	}
	else if(m==7){
	       	u[0] = sqrt(1.0/3.0);
		u[1] = -sqrt(1.0/3.0);
		u[2] = sqrt(1.0/3.0);
	}
	else{

	}


   }
 
    //cout<<u[0]<<" "<<u[1]<<" "<<u[2]<<endl;

}


