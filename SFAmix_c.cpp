#include <stdlib.h>
#include <unistd.h>
#include <ctime>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>

#include <iostream>
#include <Eigen/Dense>
#include "myHeader.cpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

#include <dirent.h>
#include <errno.h>

using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;

//const double PI  =3.141592653589793238462;
/*
  Chuan Gao C++
*/

/* 
   usage
   ./SFAmix_conditional --nf 50 --y ../data/sim_data_SFA/Y_1.txt --out result --sep tab --a 0.1 --b 0.1 --c 0.5 --d 0.5 --e 1 --f 1
   ./SFAmix_conditional --nf 50 --y ../data/sim_data_SFA/Y_1.txt --out result --sep tab --a 0.5 --b 0.5 --c 0.5 --d 0.5 --e 0.5 --f 0.5
   ./SFAmix_conditional --nf 50 --y ../data/sim_data_SFA/Y_1.txt --out result --sep tab --a 1 --b 0.5 --c 1 --d 0.5 --e 1 --f 0.5
   ./SFAmix_conditional --nf 50 --y ../data/sim_data_SFA/Y_1.txt --out result --sep tab --a 1 --b 1 --c 1 --d 1 --e 1 --f 1
   ./SFAmix --nf 20 --y ../data/sim_data_SFA/Y_1.txt --out result --sep tab --a 1 --b 1 --c 1 --d 1 --e 1 --f 1
   ./SFAmix --nf 20 --y ../data/sim_data_SFA/Y_1.txt --out result --sep tab --a 0.5 --b 0.5 --c 0.5 --d 0.5 --e 0.5 --f 0.5
   ./SFAmix --y /gpfs/fs0/data/engelhardtlab/cg148/data/CAP/S480RjQ2N_t_noCrossHyb.dat --nf 50 --sep space --out result
*/

int main(int argc,char *argv[]){
    
    // declare variables 
	int nf=100,nf2=0,s_n=0,d_y=0,i_seed=0,n_itr=5001,write_itr=50;
    double a=0.5,b=0.5,c=0.5,d=0.5,g=0.5,h=0.5,alpha=1,beta=1,nf_keep=10;
    
    string file_y,dir_out,sep;
    stringstream ss;

    sep="tab";
    int interval=200;
    
    // read in argument
    string s_df="--nf",s_y="--y",s_out="--out",s_sep="--sep",s_a="--a",s_b="--b",s_c="--c",s_d="--d",s_e="--e",s_f="--f",s_seed="--seed",s_itr="--itr",s_interval="--interval",s_write_itr="--write_itr",s_df_keep="--nf_keep";
    for(int i=0;i<argc;i++){
        if(s_df.compare(argv[i])==0){nf=atoi(argv[i+1]);}
        if(s_y.compare(argv[i])==0){file_y=argv[i+1];}
        if(s_out.compare(argv[i])==0){dir_out=argv[i+1];}
        if(s_sep.compare(argv[i])==0){sep=argv[i+1];}
        if(s_a.compare(argv[i])==0){a=atof(argv[i+1]);}
        if(s_b.compare(argv[i])==0){b=atof(argv[i+1]);}
        if(s_c.compare(argv[i])==0){c=atof(argv[i+1]);}
        if(s_d.compare(argv[i])==0){d=atof(argv[i+1]);}
		if(s_e.compare(argv[i])==0){g=atof(argv[i+1]);}
        if(s_f.compare(argv[i])==0){h=atof(argv[i+1]);}
        if(s_interval.compare(argv[i])==0){interval=atoi(argv[i+1]);}
		if(s_write_itr.compare(argv[i])==0){write_itr=atoi(argv[i+1]);}
        //if(s_seed.compare(argv[i])==0){i_seed=atoi(argv[i+1]);}
        if(s_itr.compare(argv[i])==0){n_itr=atoi(argv[i+1]);}
        if(s_df_keep.compare(argv[i])==0){nf_keep=atoi(argv[i+1]);}
    }
    // convert directory name to char_array
    int n_char = dir_out.length();
    // declaring character array
    char char_array[n_char+1];
    // copying the contents of the
    // string to char array
    strcpy(char_array, dir_out.c_str());
    DIR* dir = opendir(char_array);
    if (dir)
    {
        /* Directory exists. */
        closedir(dir);
    }
    else
    {
        /* Directory does not exist. */
        printf("Can't open results directory, stopping. \n");exit(0);
    }
    //else
    //{
    //    /* opendir() failed for some other reason. */
    //}
    
    // calculate the sample size and the gene numbers 
    string line;
    string field;
    ifstream f;
    
    f.open(file_y.c_str());
    if (! f.is_open())
		{
			printf("Gene expression file open failed\n");exit(0);
		}
    getline(f,line);
    s_n++;
    istringstream iss(line);
    if(sep.compare("space")==0){
        while(getline(iss,field,' ')){d_y++;}
    }else if(sep.compare("tab")==0){
        while(getline(iss,field,'\t')){d_y++;}
    }else{
        cout << "Please specify a valid separator." << endl << endl;
    }
    while(getline(f,line)){s_n++;}

    // write command into file for later references, also write the dimension of the gene expression matrix 
    ss.str("");
    ss.clear();
    ss << dir_out << "/command.txt";
    ofstream f_com (ss.str().c_str());
    if (f_com.is_open()){
        for(int i=0;i<argc;i++){
            f_com << argv[i] << " ";
        }
        f_com << endl;
    }
    f_com << "Y matrix has dimension of " << s_n << " by " << d_y << endl << endl;
    f_com << "Starting analysis using factor number of " << nf << endl;
    //f_com << "Command is written in command.txt" << endl;
    f_com << "The total number of runs is set to "  << n_itr << endl;
    f_com << "Results will be written for every " << write_itr << " iterations" << endl;
    f_com << "Convergence is reached if the total number of sparse elements remain unchanged for " << interval << " iterations" << endl;
    f_com << "The algorithm impose strong shrinkage, to avoid over shrinkage, it is set to keep approximately " << nf_keep << " factors at least" << endl;
    f_com.close();
    
    cout << "Starting analysis using factor number of " << nf << endl;
    cout << "Details of the run parameters can be found in command.txt" << endl;
    cout << "The total number of runs is set to "  << n_itr << endl;
    cout << "Results will be written for every " << write_itr << " iterations" << endl;
    cout << "Convergence is reached if the total number of sparse elements remain unchanged for " << interval << " iterations" << endl;
    cout << "The algorithm impose strong shrinkage, to avoid over shrinkage, it is set to keep approximately " << nf_keep << " factors at least" << endl;
    
    // read in the Y matrix 
    MatrixXd Y_TMP=MatrixXd::Constant(s_n,d_y,0);
    
    f.clear();
    f.seekg (0, ios_base::beg);
    int i=0,j=0;
    while(getline(f,line)){
        istringstream iss(line);
        j=0;
        if(sep.compare("space")==0){
            while(getline(iss,field,' ')){
                Y_TMP(i,j)=atof(field.c_str());
                j++;
            }
            i++;
        }else{
            while(getline(iss,field,'\t')){
                Y_TMP(i,j)=atof(field.c_str());
                j++;
            }
            i++;
        }
    }
    f.close();
    
    // chunk Y to a smaller sub matrix to test if necessary
	//s_n=500;
    MatrixXd Y=MatrixXd::Constant(s_n,d_y,0);
    Y=Y_TMP.block(0,0,s_n,d_y);
    
    
    // Declare variables independent of factor number to prepare for the EM algorithm
    
    VectorXd psi_v = VectorXd::Constant(s_n,1);
    VectorXd PSI=VectorXd::Constant(s_n,1);
    VectorXd PSI_INV=VectorXd::Constant(s_n,1);
    MatrixXd LX=MatrixXd::Constant(s_n,d_y,0);
    MatrixXd YLX=MatrixXd::Constant(s_n,d_y,0);
    MatrixXd YLXi=VectorXd::Constant(s_n,0);
    MatrixXd YLX2=MatrixXd::Constant(s_n,d_y,0);
    VectorXd TOP=VectorXd::Constant(s_n,0);
    
    // Declare variables that are dependent of the factor number
    nf2 = nf;
    int nt=nf;
    
    MatrixXd EX=MatrixXd::Constant(nt,d_y,0);
    MatrixXd TEX=MatrixXd::Constant(d_y,nt,0);
    MatrixXd VX=MatrixXd::Constant(nt,nt,0);
    MatrixXd EXX=MatrixXd::Constant(nt,nt,0);
    MatrixXd LP=MatrixXd::Constant(nt,s_n,0);
    MatrixXd ID=MatrixXd::Constant(nt,nt,0);
    VectorXd id_v = VectorXd::Constant(nt,1);
    
    MatrixXd LAM=MatrixXd::Constant(s_n,nt,0);
	//MatrixXd LAM_T=MatrixXd::Constant(nt,s_n,0);
    MatrixXd LAM_BAK=MatrixXd::Constant(s_n,nt,0);
    MatrixXd THETA=MatrixXd::Constant(s_n,nf,0);
    MatrixXd DELTA=MatrixXd::Constant(s_n,nf,0);
    VectorXd PHI = VectorXd::Constant(nf,0);
    VectorXd TAU = VectorXd::Constant(nf,0);
    double nu = 1;
    double ETA = 1;
    double GAMMA = 1;
    
    VectorXd count_lam = VectorXd::Constant(nf,0);
    VectorXd index = VectorXd::Constant(nf,0);
    
    
    MatrixXd LAM_TOP=MatrixXd::Constant(s_n,nf,0);
    MatrixXd LAM_BOT=MatrixXd::Constant(s_n,nf,0);
    
    double nmix=2;
    double zi = double(1)/nmix;
    MatrixXd Z = MatrixXd::Constant(nmix,nf,zi);
    MatrixXd logZ = MatrixXd::Constant(nmix,nf,log(zi));
    MatrixXd LOGV = MatrixXd::Constant(nmix,1,log(0.5));
    MatrixXd Zm = MatrixXd::Constant(nmix,nf,log(zi));

	MatrixXd LPL = MatrixXd::Constant(nf,nf,0);

	//MatrixXd VXM=MatrixXd::Constant(d_y,nf,0);

	long seed;
	gsl_rng *r;  // random number generator
	r=gsl_rng_alloc(gsl_rng_mt19937);
	
	seed = time (NULL) * getpid();
	//seed=20877551893284;
	gsl_rng_set (r, seed);                  // set seed

	ss.str("");
    ss.clear();
    ss << dir_out << "/seed";
    ofstream f_seed (ss.str().c_str());
    f_seed << seed << endl;
    f_seed.close();
       
	//gsl_rng_set (r, i_seed);
	
    // Initialize parameters 
    ID.diagonal()=id_v;
   
    //PSI.diagonal() = psi_v;
    //inv_psi(PSI,PSI_INV,s_n);

    // fill in the lambda matrix  
    for (int i=0; i<s_n; i++) {
        for(int j=0;j<nt;j++){
            LAM(i,j)=gsl_ran_gaussian(r,1);
        }
    }
	
    for(int i=0;i<s_n;i++){
        for(int j=0;j<nf;j++){
            THETA(i,j)=1;
            DELTA(i,j)=1;
        }
    }
    
    // continue initializing 
    for(int i=0;i<nf;i++){
        PHI(i)=1;
        TAU(i)=1;
    }
    VectorXd lam_count_v = VectorXd::Constant(n_itr,0);

	MatrixXd LAM_T=LAM.transpose();
    for(int i=0;i<s_n;i++){
        LP.col(i)=LAM_T.col(i)*PSI_INV(i);
    }
	VX=(LP*LAM+ID).inverse();
	EX=VX*LP*Y;
	EXX=EX*EX.transpose()+VX*d_y;
	
	TEX=EX.transpose();

	// int gene_start=0;
	// int gene_stop=0;
	
    for(int itr=0;itr<(n_itr-1);itr++){

        
        // Expectation step

     	MatrixXd LAM_T=LAM.transpose();
		for(int i=0;i<s_n;i++){
			LP.col(i)=LAM_T.col(i)*PSI_INV(i);
		}
        //VX=(LPL+ID).lu().solve(ID);
		VX=(LP*LAM+ID).lu().solve(ID);
        EX=VX*LP*Y;
        EXX=EX*EX.transpose()+VX*d_y;		
        TEX=EX.transpose();
        
        // Mazimization step 
        double alpha_sum=TAU.sum();
        ETA=double((d*nf+g-1))/(GAMMA+alpha_sum);
        GAMMA=double(g+h)/(ETA+nu);
        for(int i=0;i<nf;i++){
            TAU(i)=double(c+d)/(PHI(i)+ETA);
			
		}
		for(int i=0;i<nf;i++){
			double sum_c=s_n*b*Z(0,i)+c-1-0.5*s_n*Z(1,i);
			double at = 2*(TAU(i)+Z(0,i)*(DELTA.col(i).sum()));
			double lam_sum=LAM.col(i).dot(LAM.col(i));
			double bt = Z(1,i)*lam_sum;
			PHI(i)=double(sum_c+sqrt(sum_c*sum_c+at*bt))/at;
            //if(nf <= nf_keep){
            //    PHI(i) = PHI(i) + 1e-280;
            //}
			// if(PHI(i)<1e-300){
			// 	PHI(i)=1e-300;
			// }
			
		}

		double min_phi=1e100;
		for(int i=1;i<nf;i++){
			if(PHI(i)!=0&PHI(i)<min_phi){
				min_phi=PHI(i);
			}
		}
		for(int i=1;i<nf;i++){
			if(PHI(i)==0){
				PHI(i)=min_phi;
			}
		}
		
		// count the number of loadings that are set to all zero
		count_lam.setZero();
		for(int i=0;i<nf;i++){
			for(int j=0;j<s_n;j++){
				if(LAM(j,i)!=0){
					count_lam(i) +=  1;
				}
			}
		}

	
		// Count the number of loadings that are active, either all zero or PHI_k zero will kill 
		int count_nonzero = 0;
		for(int i=0;i<nf;i++){
			if(count_lam(i)>1&&PHI(i)!=0){
				index(count_nonzero)=i;
				count_nonzero++;
			}
		}
		// remove inactive factors, loadings etc. and assign to new matrix
		//cout << "number of factors before " << nf << endl;
		if(count_nonzero != nf){

			nf=count_nonzero;
			nt=nf;
			MatrixXd EX2=MatrixXd::Constant(nt,d_y,0);
			MatrixXd TEX2=MatrixXd::Constant(d_y,nt,0);
			MatrixXd VX2=MatrixXd::Constant(nt,nt,0);
			MatrixXd EXX2=MatrixXd::Constant(nt,nt,0);
			MatrixXd LP2=MatrixXd::Constant(nt,s_n,0);
			MatrixXd ID2=MatrixXd::Constant(nt,nt,0);
			VectorXd id_v2 = VectorXd::Constant(nt,1);
			MatrixXd LAM2=MatrixXd::Constant(s_n,nt,0);
			MatrixXd LAM_BAK2=MatrixXd::Constant(s_n,nt,0);
			MatrixXd THETA2=MatrixXd::Constant(s_n,nf,0);
			MatrixXd DELTA2=MatrixXd::Constant(s_n,nf,0);
			VectorXd PHI2 = VectorXd::Constant(nf,0);
			VectorXd TAU2 = VectorXd::Constant(nf,0);
			MatrixXd Z2 = MatrixXd::Constant(nmix,nf,zi);
			MatrixXd logZ2 = MatrixXd::Constant(nmix,nf,log(zi));
			VectorXd count_lam2 = VectorXd::Constant(nf,0);
			VectorXd index2 = VectorXd::Constant(nf,0);
			MatrixXd LAM_TOP2=MatrixXd::Constant(s_n,nt,0);
			MatrixXd LAM_BOT2=MatrixXd::Constant(s_n,nt,0);
			MatrixXd LPL2 = MatrixXd::Constant(nf,nf,0);
			ID2.diagonal()=id_v2;
                     
			for(int i=0;i<nf;i++){
				EX2.row(i)=EX.row(index(i));
				//TEX2=EX2.transpose();
				for(int j=0;j<nf;j++){
					VX2(i,j)=VX(index(i),index(j));
					EXX2(i,j)=EXX(index(i),index(j));
					LPL2(i,j)=LPL(index(i),index(j));
				}
				LP2.row(i)=LP.row(index(i));
				LAM2.col(i)=LAM.col(index(i));
				LAM_BAK2.col(i)=LAM_BAK.col(index(i));
				THETA2.col(i)=THETA.col(index(i));
				DELTA2.col(i)=DELTA.col(index(i));
				PHI2(i)=PHI(index(i));
				TAU2(i)=TAU(index(i));
				Z2.col(i)=Z.col(index(i));
				logZ2.col(i)=logZ.col(index(i));
				LAM_TOP2.col(i)=LAM_TOP.col(index(i));
				LAM_BOT2.col(i)=LAM_BOT.col(index(i));
				count_lam2(i)=count_lam(index(i));
				index2(i)=index(i);
			}
            
			// Assign the new parameters back 
			EX=EX2;
			//TEX=TEX2;
			VX=VX2;
			EXX=EXX2;
			LP=LP2;
			ID=ID2;
			LAM=LAM2;
			LAM_BAK=LAM_BAK2;
			THETA=THETA2;
			DELTA=DELTA2;
			PHI=PHI2;
			TAU=TAU2;
			Z=Z2;
			logZ=logZ2;
			LAM_TOP=LAM_TOP2;
			LAM_BOT=LAM_BOT2;
			count_lam=count_lam2;
			index=index2;
			LPL=LPL2;
            
		}
		//cout << "number of factors after" << nf << endl;
        
		// continue updating parameters 
		for(int i=0;i<s_n;i++){
			for(int j=0;j<nf;j++){
				DELTA(i,j)=double((a+b))/(THETA(i,j)+PHI(j));
			}
		}
		
		for(int i=0;i<s_n;i++){
			for(int j=0;j<nf;j++){
				double a23=(2*a-3);
				THETA(i,j)=double(a23+sqrt(a23*a23+8*LAM(i,j)*LAM(i,j)*DELTA(i,j)))/4/DELTA(i,j);
			}
		}

        
		// Need to update separately because 0/0=NA
		/*
		  VectorXd indexALL = VectorXd::Constant(nf,0);
		  MatrixXd partV=MatrixXd::Constant(nf,nf,0);
		  MatrixXd partL = MatrixXd::Constant(s_n,nf,0);

		  cal_lam(LAM,indexALL, partL, Y, EX, PSI_INV, EXX, Z, LPL, THETA,PHI, partV, s_n,nf);
		*/

				
		YLX=Y-LAM*EX;
		for(int i=0;i<nf;i++){
			YLX=YLX+LAM.col(i)*(EX.row(i));
			//YLX=YLX+LAM.col(i)*(TEX.col(i).transpose());
			for(int j=0;j<s_n;j++){
				TOP(j) = YLX.row(j).dot(EX.row(i));
			}
                    
			for(int j=0;j<s_n;j++){
				if(Z(0,i)==0){
					LAM(j,i) = TOP(j)/(EXX(i,i)+PSI(j)*Z(1,i)/PHI(i));
				}
				else if(Z(1,i)==0){
					LAM(j,i) = TOP(j)/(EXX(i,i)+PSI(j)*Z(0,i)/THETA(j,i));
				}
				else{
					LAM(j,i) = TOP(j)/(EXX(i,i)+PSI(j)*(Z(1,i)/PHI(i)+Z(0,i)/THETA(j,i)));
				}
			}
			YLX=YLX-LAM.col(i)*EX.row(i);                   
		}
			

		/*
		  YLX=Y*TEX-LAM*EXX;
		  for(int i=0;i<nf;i++){
		  YLX.col(i)=YLX.col(i)+LAM.col(i)*(EXX(i,i));
		  for(int j=0;j<s_n;j++){
		  if(Z(0,i)==0){
		  LAM(j,i) = YLX(j,i)/(EXX(i,i)+PSI(j,j)*Z(1,i)/PHI(i));
		  }
		  else if(Z(1,i)==0){
		  LAM(j,i) = YLX(j,i)/(EXX(i,i)+PSI(j,j)*Z(0,i)/THETA(j,i));
		  }
		  else{
		  LAM(j,i) = YLX(j,i)/(EXX(i,i)+PSI(j,j)*(Z(1,i)/PHI(i)+Z(0,i)/THETA(j,i)));
		  }
		  }
		  YLX.col(i)=YLX.col(i)-LAM.col(i)*(EXX(i,i));
		  }
		*/		
		// If theta=0, then LAM=0

		/*
		  VectorXd maxEX = VectorXd::Constant(nf,0);
		
		  for(int k=0;k<nf;k++){
		  maxEX(k)=EX(k,0);
		  for(int i=1;i<d_y;i++){
		  if(abs(maxEX(k))<abs(EX(k,i))){
		  maxEX(k)=EX(k,i);
		  }
				
		  }
	
		  for(int i=1;i<s_n;i++){
			
		  if(abs(maxEX(k)*LAM(i,k))<1e-10){
		  LAM(i,k)=0;
		  }
		  }
	
		  }
		
		*/
		for(int i=0;i<s_n;i++){
			for(int j=0;j<nf;j++){
				if(THETA(i,j)==0){
					LAM(i,j)=0;
				}
				if(abs(LAM(i,j))<1e-10){
					LAM(i,j)=0;
				}
			}
		}
		
		// logZ 
		for(int i=0;i<nf;i++){
			logZ(0,i)=LOGV(0,0);
			logZ(1,i)=LOGV(1,0);
			for(int j=0;j<s_n;j++){
				logZ(0,i)=logZ(0,i)+log_norm(LAM(j,i),0,THETA(j,i))+log_gamma(THETA(j,i),a,DELTA(j,i))+log_gamma(DELTA(j,i),b,PHI(i));
				logZ(1,i)=logZ(1,i)+log_norm(LAM(j,i),0,PHI(i));
			}

		}	
        
		// Z 
		for(int i=0;i<nf;i++){
			Z(0,i)=double(1)/(1+exp(logZ(1,i)-logZ(0,i)));
			Z(1,i)=1-Z(0,i);
		}
        
		double ps1=alpha;
		double ps2=beta;
        
		// Probability 
		for(int i=0;i<nf;i++){
			ps1=ps1+Z(0,i);
			ps2=ps2+Z(1,i);
		}
		double dgama = gsl_sf_psi(ps1+ps2);
		LOGV(0,0)=gsl_sf_psi(ps1)-dgama;
		LOGV(1,0)=gsl_sf_psi(ps2)-dgama;
               
		// PSI 
		LX=LAM*EX;
		for(int i=0;i<s_n;i++){
			PSI(i)=Y.row(i).dot(Y.row(i))-2*(LX.row(i)).dot(Y.row(i))+(LAM.row(i)*EXX).dot(LAM.row(i));
		}
		PSI=PSI/d_y;
		inv_psi_vec(PSI,PSI_INV,s_n);
        
		// LAM active 
		int lam_count=0;
		for(int i=0;i<s_n;i++){
			for(int j=0;j<nf;j++){
				if(LAM(i,j)!=0){
					lam_count++;
				}
			}
		}

		MatrixXd range_LAM = MatrixXd::Constant(2,nf,0);
		MatrixXd range_EX = MatrixXd::Constant(2,nf,0);

		//MatrixXd range_THETA = MatrixXd::Constant(2,nf,0);

		//MatrixXd range_LAM_cov = MatrixXd::Constant(2,nfc,0);
		//MatrixXd range_EX_cov = MatrixXd::Constant(2,nfc,0);	
		 
		//range_colwise(range_LAM,LAM,s_n,nf);
		//range_rowwise(range_EX,EX,nf,d_y);
		//range_colwise(range_THETA,THETA,s_n,nf);
	
		if(itr%10==0){
			cout << "after reduction" << endl;
			cout << "itr " << itr << endl;
			cout << "number of factors " << nf << endl;
			
			cout << "count_lam" << endl << count_lam.transpose() << endl;
			//cout << "PHI" << endl << PHI.transpose() << endl;
			//cout << a << b << c << d << g << h << endl;

			//cout << "range LAM\n" << std::setprecision(3) << range_LAM << endl;
			//cout << "range EX\n" << std::setprecision(3) << range_EX << endl;
			//cout << "range THETA\n" << std::setprecision(3) << range_THETA << endl;
			//cout << "range Z\n" << std::setprecision(3) << Z << endl;
			//cout << "range logZ\n" << std::setprecision(3) << logZ << endl;
		}
		
		lam_count_v(itr+1)=lam_count;
		//int interval=50;
		//cout << "lam_count " << lam_count << endl;

		// claim convergence if the number of values is stable for 200 iterations.
		if(itr>=interval){
			//if(lam_count_v(itr+1)!=(s_n*nf)&&(lam_count_v(itr-interval)-lam_count_v(itr+1))<(0.05*lam_count_v(itr-interval))){
			if(lam_count_v(itr)!=(s_n*nf)&&(lam_count_v(itr-interval)-lam_count_v(itr))==0){
				ss.str("");
				ss.clear();
				ss << dir_out << "/itr";
				ofstream f_itr (ss.str().c_str());
				if (f_itr.is_open()){
					f_itr << itr << endl;
				}
				f_itr.close();
    
				ss.str("");
				ss.clear();
				ss << dir_out << "/LAM";
				ofstream f_lam (ss.str().c_str());
				if (f_lam.is_open()){
					f_lam << LAM << endl;
				}
				f_lam.close();
                
				ss.str("");
				ss.clear();
				ss << dir_out << "/Z";
				ofstream f_Z (ss.str().c_str());
				if (f_Z.is_open()){
					f_Z << Z << endl;
				}
				f_Z.close();
                
				// ss.str("");
				// ss.clear();
				// ss << dir_out << "/LOGV";
				// ofstream f_V (ss.str().c_str());
				// if (f_V.is_open()){
				// 	f_V << LOGV << endl;
				// }
				// f_V.close();
                            
				ss.str("");
				ss.clear();
				ss << dir_out << "/EX";
				ofstream f_EX (ss.str().c_str());
				if (f_EX.is_open()){
					f_EX << EX << endl;
				}
				f_EX.close();
                
				ss.str("");
				ss.clear();
				ss << dir_out << "/EXX";
				ofstream f_EXX (ss.str().c_str());
				if (f_EXX.is_open()){
					f_EXX << EXX << endl;
				}
				f_EXX.close();

				ss.str("");
				ss.clear();
				ss << dir_out << "/PSI";
				ofstream f_PSI (ss.str().c_str());
				if (f_PSI.is_open()){
					f_PSI << PSI << endl;
				}
				f_PSI.close();

				ss.str("");
				ss.clear();
				ss << dir_out << "/lam_count";
				ofstream f_lam_count (ss.str().c_str());
				if (f_lam_count.is_open()){
					f_lam_count << count_lam.transpose() << endl;
				}
				f_lam_count.close();

                
				
				exit(0);
 
			}
	
                
		}

		if((itr%write_itr==0)){
              
			ss.str("");
			ss.clear();
			ss << dir_out << "/itr";
			ofstream f_itr (ss.str().c_str());
			if (f_itr.is_open()){
				f_itr << itr << endl;
			}
			f_itr.close();
    
			ss.str("");
			ss.clear();
			ss << dir_out << "/LAM_" << itr;
			ofstream f_lam (ss.str().c_str());
			if (f_lam.is_open()){
				f_lam << LAM << endl;
			}
			f_lam.close();
                
			ss.str("");
			ss.clear();
			ss << dir_out << "/Z_" << itr;
			ofstream f_Z (ss.str().c_str());
			if (f_Z.is_open()){
				f_Z << Z << endl;
			}
			f_Z.close();
                
			// ss.str("");
			// ss.clear();
			// ss << dir_out << "/LOGV_" <<  itr;
			// ofstream f_V (ss.str().c_str());
			// if (f_V.is_open()){
			// 	f_V << LOGV << endl;
			// }
			// f_V.close();
                
                
			ss.str("");
			ss.clear();
			ss << dir_out << "/EX_" << itr;
			ofstream f_EX (ss.str().c_str());
			if (f_EX.is_open()){
				f_EX << EX << endl;
			}
			f_EX.close();
                
			ss.str("");
			ss.clear();
			ss << dir_out << "/EXX_" << itr;
			ofstream f_EXX (ss.str().c_str());
			if (f_EXX.is_open()){
				f_EXX << EXX << endl;
			}
			f_EXX.close();

			ss.str("");
			ss.clear();
			ss << dir_out << "/PSI_" << itr;
			ofstream f_PSI (ss.str().c_str());
			if (f_PSI.is_open()){
				f_PSI << PSI << endl;
			}
			f_PSI.close();

            ss.str("");
            ss.clear();
            ss << dir_out << "/Phi_" << itr;
            ofstream f_Phi (ss.str().c_str());
            if (f_Phi.is_open()){
                f_Phi << PHI << endl;
            }
            f_PSI.close();
            
			ss.str("");
			ss.clear();
			ss << dir_out << "/lam_count_" << itr;
			ofstream f_lam_count (ss.str().c_str());
			if (f_lam_count.is_open()){
				f_lam_count << count_lam.transpose() << endl;
			}
			f_lam_count.close();
    
		}
	
		// if(itr>interval && lam_count_v(itr+1)-lam_count_v(itr-interval)==0){
		//     ss.str("");
		//     ss.clear();
		//     ss << dir_out << "/final";
		//     ofstream f_FINAL (ss.str().c_str());
		//     if (f_FINAL.is_open()){
		//         f_FINAL << "done" << endl;
		//     }
		//     f_FINAL.close();
		//     exit(0);
		// }
        
	}

}


