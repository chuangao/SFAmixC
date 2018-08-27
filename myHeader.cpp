#include <stdlib.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <Eigen/Dense>
#include <limits>
#include <stdio.h>
#include <ostream>
#include <string>
#include <iostream>
#include <float.h>
//#include <boost/math/special_functions/gamma.hpp>

using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;

const double PI  =3.141592653589793238462;

void range_colwise(MatrixXd& range,MatrixXd& x,int n,int m){
	int i,j;
	for(i=0;i<m;i++){
		//cout << "i" << i << endl;
		int s;
		double min=1E300,max=0;
		// for(s=0;s<n;s++){
		// 	if(x(s,i)!=0){
		// 		min=abs(x(s,i));
		// 		max=abs(x(s,i));
		// 		continue;		
		// 	}
		// 	//cout << "s" << s << endl;
		// }
		for(j=s;j<n;j++){
			if(x(j,i)!=0){
				if(min>abs(x(j,i))){
					min=abs(x(j,i));
				}
				if(max<abs(x(j,i))){
					max=abs(x(j,i));
				}
			}
		}
		range(0,i)=min;
		range(1,i)=max;
	}
	//cout << "range_coleise_finished" << endl;
}

void range_rowwise(MatrixXd& range,MatrixXd& x,int n,int m){
	int i,j;
	for(j=0;j<n;j++){
		int s;
		double min=1E300,max=0;
		// for(s=0;s<m;s++){
		// 	if(x(j,s)!=0){
		// 		min=abs(x(j,s));
		// 		max=abs(x(j,s));
		// 		continue;
		// 	}
			
		// }
		for(i=s;i<m;i++){
			if(x(j,i)!=0){
				if(min>abs(x(j,i))){
					min=abs(x(j,i));
				}
				if(max<abs(x(j,i))){
					max=abs(x(j,i));
				}
			}
		}
		range(0,j)=min;
		range(1,j)=max;
	}
}

void cpy_row_matrix(MatrixXd& M1,MatrixXd& M2,VectorXd& index, int n,int row){
    for(int i=0;i<n;i++){
        M1(0,i)=M2(row,index(i));
    }
}
void cpy_row_matrix_bak(MatrixXd& M1,MatrixXd& M2,VectorXd& index, int n,int row){
    for(int i=0;i<n;i++){
        M1(row,index(i))=M2(0,i);
    }
}

void inv_psi(Eigen::MatrixXd& psi,Eigen::MatrixXd& psi_inv,int n){
    for(int i=0;i<n;i++){
        psi_inv(i,i)=double(1)/psi(i,i);
    }
}
void inv_psi_vec(Eigen::VectorXd& psi,Eigen::VectorXd& psi_inv,int n){
    for(int i=0;i<n;i++){
        psi_inv(i)=double(1)/psi(i);
    }
}
double log_norm(double x,double e, double v){
  //const double PI  =3.141592653589793238462;
    double a = std::numeric_limits<double>::infinity();
    //cout << "a " << a << endl;
	
    if(v==0){
        return a;
    }else{
        return -0.5*(x-e)*(x-e)/v-0.5*log(v)-0.5*log(2*PI);
		//return(log(gsl_ran_gaussian_pdf(x-e,sqrt(v))));
    } 
	
	
}
double log_gamma(double x, double a, double beta){
    //return a*log(beta)-log(gsl_sf_gamma(a))+(a-1)*log(x)-beta*x;
	if(a==1&&x==0){
		return(beta);
	}else{
		return a*log(beta)-lgamma(a)+(a-1)*log(x)-beta*x;
	}
	//return(log(gsl_ran_gamma_pdf(x,a,double(1)/beta)));
    //boost::lambda
}
void cumsum(Eigen::VectorXd& S, Eigen::VectorXd& D, int n){
  for(int i=0;i<n;i++){
    if(i==0){
      D(i)=S(i);
    }
    D(i)=D(i-1)+S(i);
  }
}

double fx(double x,double c){
    return -1*log(x)+0.5/x+c;
}

double dfx(double x){
    return -1*(double(1)/x+0.5/x/x);
}

double NR(double c){
    double x=1e-10;
    for(int i=0;i<500;i++){
        x=x-fx(x,c)/dfx(x);
    }
    return x;
}


void cal_lam(MatrixXd& LAM,VectorXd& indexALL,MatrixXd& partL,MatrixXd& Y,MatrixXd& EX,MatrixXd& PSI_INV,MatrixXd& EXX,MatrixXd& Z,MatrixXd& LPL,MatrixXd& THETA,VectorXd& PHI,MatrixXd& partV, int s_n,int nf){
	indexALL.setZero();
	partL=PSI_INV*Y*EX.transpose();
	//cout << "partL " << endl << partL.block(0,0,5,5) << endl;
	//LPL=LAM.transpose()*PSI_INV*LAM;
	LPL.setZero();
	
	//LAM.setZero();
        
	//MatrixXd partV=MatrixXd::Constant(nf,nf,0);
        
	for(int j=0;j<s_n;j++){
		int count_indexALL=0;
		for(int i=0;i<nf;i++){              
			if(THETA(j,i)!=0 && PHI(i)!=0){
				indexALL(count_indexALL)=i;
				count_indexALL++;
			}
		}
            
		if(count_indexALL==0){
			LAM.row(j).setZero();
			continue;
		}
            
		partV.setZero();
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partV(indexALL(i1),indexALL(i2)) = PSI_INV(j,j)*EXX(indexALL(i1),indexALL(i2));
			}
			partV(indexALL(i1),indexALL(i1)) += Z(0,indexALL(i1))/THETA(j,indexALL(i1))+Z(1,indexALL(i1))/PHI(indexALL(i1));
		}  
       
		MatrixXd partVI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partVI(i1,i2) = partV(indexALL(i1),indexALL(i2));
			}
		}
							
		MatrixXd partLI=MatrixXd::Constant(1,count_indexALL,0);
		MatrixXd LAMI=MatrixXd::Constant(1,count_indexALL,0);
		MatrixXd LLI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		cpy_row_matrix(partLI,partL,indexALL,count_indexALL,j);
            
		MatrixXd IDNF = MatrixXd::Identity(count_indexALL, count_indexALL);

		//LDLT<MatrixXd> ldltOfA(partVI);
		//MatrixXd vl=ldltOfA.solve(IDNF);
			
		MatrixXd vl = partVI.lu().solve(IDNF);
		LAMI=partLI*vl;

		for(int i=0;i<count_indexALL;i++){
			if(THETA(j,indexALL(i))==0){
				LAMI(i)=0;
			}
		}
            
		cpy_row_matrix_bak(LAM,LAMI,indexALL,count_indexALL,j);
			
		LLI=LAMI.transpose()*LAMI;


		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j,j)*(LLI(i1,i2)+vl(i1,i2));
				//ELL(indexALL(i1),indexALL(i2)) += LLI(i1,i2)+vl(i1,i2);
				//LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j,j)*(LLI(i1,i2));
				//ELL(indexALL(i1),indexALL(i2)) += LLI(i1,i2);
				//vLXL(j) += vl(i1,i2)*EXX(indexALL(i1),indexALL(i2));
					
			}
		}
			
			
	}

}


double like(Eigen::MatrixXd& PSI,Eigen::MatrixXd& EXX,Eigen::MatrixXd& LAM,Eigen::MatrixXd& THETA,Eigen::MatrixXd& DELTA,Eigen::VectorXd& PHI,Eigen::VectorXd& TAU,Eigen::MatrixXd& Z,Eigen::MatrixXd& V,double ETA, double GAMMA,double alpha, double beta, int n,int p, int nf,double a, double b, double c, double d, double e, double f, double nu){
  double det_psi=0;
  
  for(int i=0;i<n;i++){
    det_psi = det_psi+log(PSI(i,i));
  }
  
  double like=(-1)*0.5*n*p*log(2*PI)-0.5*p*det_psi;
  
  double sum_x=0;
  for(int i=0;i<nf;i++){
    sum_x = sum_x + (-1)*0.5*EXX(i,i);
  }
  
  like = like - 0.5*nf*p*log(2*PI) + sum_x;
  
  for(int i=0;i<nf;i++){
    like=like + (Z(0,i)+alpha)*V(0,i) + (Z(1,i)+beta)*V(1,i);
  }
  
  for(int i=0;i<n;i++){
    for(int j=0;j<nf;j++){
      if(THETA(i,j)!=0){
        like=like+Z(0,j)*log_norm(LAM(i,j),0,THETA(i,j));
        //if(DELTA(i,j)!=0){
          like=like+Z(0,j)*log_gamma(THETA(i,j),a,DELTA(i,j));
          //}
      }
      like=like+Z(0,j)*log_gamma(DELTA(i,j),b,PHI(j));
      like=like+Z(1,j)*log_norm(LAM(i,j),0,PHI(j));
    }
  }
 
 for(int i=0;i<nf;i++){
   like=like+log_gamma(PHI(i),c,TAU(i));
   like=like+log_gamma(TAU(i),d,ETA);
 }
 
 like=like+log_gamma(ETA,e,GAMMA);
 like=like+log_gamma(GAMMA,f,nu);
  
  return like;
 
}
/*
#include <iostream>     // cout
#include <math.h>       // acos
#include <float.h>      // DBL_MAX
#include <limits>       // numeric_limits
*/
template<typename T>
bool is_infinite( const T &value )
{
    // Since we're a template, it's wise to use std::numeric_limits<T>
    //
    // Note: std::numeric_limits<T>::min() behaves like DBL_MIN, and is the smallest absolute value possible.
    //
 
    T max_value = std::numeric_limits<T>::max();
    T min_value = - max_value;
 
    return ! ( min_value <= value && value <= max_value );
}
 
template<typename T>
bool is_nan( const T &value )
{
    // True if NAN
    return value != value;
}
 
template<typename T>
bool is_valid( const T &value )
{
    return ! is_infinite(value) && ! is_nan(value);
}

/*
int main()
{
    using std::cout;
 
    double a, b, c, d, e;
 
    a = 1.0;
    b = 0.0;
    c = a / c;          // divide by zero
    d = acos(-1.001);   // domain for acos is [-1, 1], anything else is #IND or inf
    e = b / b;          // zero / zero
 
    cout << "Value of a: " << a << " " << is_valid(a) << " " << (is_nan(a) ? " nan " : "") << (is_infinite(a) ? " infinite " : "") << "n";
    cout << "Value of b: " << b << " " << is_valid(b) << " " << (is_nan(b) ? " nan " : "") << (is_infinite(b) ? " infinite " : "") << "n";
    cout << "Value of c: " << c << " " << is_valid(c) << " " << (is_nan(c) ? " nan " : "") << (is_infinite(c) ? " infinite " : "") << "n";
    cout << "Value of d: " << d << " " << is_valid(d) << " " << (is_nan(d) ? " nan " : "") << (is_infinite(d) ? " infinite " : "") << "n";
    cout << "Value of e: " << e << " " << is_valid(e) << " " << (is_nan(e) ? " nan " : "") << (is_infinite(e) ? " infinite " : "") << "n";
 
    return 0;
}
*/
