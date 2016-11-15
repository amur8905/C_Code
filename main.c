

#include <math.h>

//define some element wise operators 
#define for_i for(i = 0; i < h; i++)
#define for_j for(j = 0; j < w; j++)
#define _M double**
#define OPM(name, _op_) \
	void eop_##name(_M a, _M b, _M c, int w, int h){int i,j;\
		for_i for_j c[i][j] = a[i][j] _op_ b[i][j];}
OPM(add, +);OPM(sub, -);OPM(mul, *);OPM(div, /);
 
#define OPS(name, res) \
	void eop_s_##name(_M a, double s, _M b, int w, int h) {double x;int i,j;\
		for_i for_j {x = a[i][j]; b[i][j] = res;}}
OPS(mul, x*s);OPS(div, x/s);OPS(add, x+s);OPS(sub, x-s);OPS(pow, pow(x, s));


#include<stdio.h>

int main()
{    

}

// calculates the value of the given RBF at X with center Xc and radius r
// rbf = rbf type (see case statement in function CalcRBF)
// r = radius, available RBFs are given in the case statement below
// X = point to evaluate rbf at
// Xc = center of rbf

double CalcRBF(int rbf,double r,double* X, double* Xc) {
	
	double Xnorm,Xr,phi;
	double pi = 3.141592653589793;
	Xnorm = sqrt(pow(X[1]-Xc[1],2)+pow(X[2]-Xc[2],2) + pow(X[3]-Xc[3],2));
	Xr = Xnorm/r;

	switch (rbf) {
		case 1: // C0
			phi = (1-Xr)*(1-Xr);
			if (Xnorm > r) {
				phi = 0;
			}
			break;
		case 2: // C2
			phi = pow((1 - Xr),4)*(4*Xr + 1);
			if (Xnorm > r) {
				phi = 0;
			}
			break;
		case 3: // C4
			 phi = pow((1-Xr),6)*(35*Xr*Xr+18*Xr+3)/3;
			if (Xnorm > r) {
				phi = 0;
			}
			break;
		case 4: // C6
			phi = pow((1-Xr),8)*(32*Xr*Xr*Xr+25*Xr*Xr+8*Xr+1);
			if (Xnorm > r) {
				phi = 0;
			}
			break;
		case 5: //Euclid
			phi = pi*((1/12*Xr*Xr*Xr)-0.25*Xr+4/3*0.25)/(pi*(-1*0.25*0+4/3*0.125));			
			if (Xnorm > r) {
				phi = 0;
			}
			break;
		case 6: // Multiquadric
			phi = sqrt(1+Xr*Xr);
			break;
		case 7: // InverseMulti
			phi = 1/sqrt(1+Xr*Xr);
			break;
		case 8: // TPS
			phi = Xr*Xr*log(Xr);
			break;
		case 9: // Gaussian
			phi = exp(-1*Xr*Xr);
			break;
		default :			
		phi = pow((1 - Xr),4)*(4*Xr + 1);
		if (Xnorm > r) {
			phi = 0;
		}
		break;
	}	
	return phi;
}

// constructs the coupling matrix H
// Xa = aerodynamic mesh points
// Xs = structural points
// degP = degree of polynomial term (1 is usually fine)
// rbf = rbf type (see case statement in function CalcRBF)
// r = radius
double** CalcRBFMatrix(double** Xa,double** Xs,int degP,char* rbf,double r) {
			
	// count the number of possible polynomials terms (tetrahedral numbers)
	int pDim = (degP+1)*(degP+2)*(degP+3)/6; 				
	
	int Ns = sizeof(Xs) / sizeof(double);	
	int Na = sizeof(Xa) / sizeof(double);		
	
	double Ps[pDim][Ns];
	double Pa[Na][pDim];
	double Pst[Ns][pDim];
	double Aas[Na][pDim+Ns];
	double F[pDim+Ns][Ns];
	double M[Ns][Ns];
	double Mi[Ns][Ns];
	double Mp[pDim][pDim];
	double right[Na][Ns];
	double H[Na][Ns];	
	
	// build the components of the Css matrix
	int ind = 1;
	int i,j,k,itr;	
	for (k = 0; k <= degP; k++) {
		for (j = 0; j <= degP; j++) {
			for (i = 0; i <= degP; i++) {
				if (i+j+k <= degP) {
					for (itr = 0; itr < Ns; itr++) {
						Ps[ind][itr] = pow(Xs[itr][1],(double)i)*pow(Xs[itr][2],(double)j)*pow(Xs[itr][3],k);
					}
				}
			}
		}
	}
	
	double Xs_i[3],Xs_j[3];
	for (i = 0; i < Ns; i++) {
		for (j = 0; j < Ns; j++) {
			Xs_i[1] = Xs[i][1];
			Xs_i[2] = Xs[i][2];
			Xs_i[3] = Xs[i][3];
			Xs_j[1] = Xs[j][1];
			Xs_j[2] = Xs[j][2];
			Xs_j[3] = Xs[j][3];
			M[i][j] = CalcRBF(rbf,r,Xs_i,Xs_j);
		}
	}
	
	Mi = inverse(M);
	Pst = transpose(Ps);
	Mp = inverse(multiply(Ps,multiply(Mi,Pst)));
	
	// build the Aas matrix		
	for (k = 0; k <= degP; k++) {
		for (j = 0; j <= degP; j++) {
			for (i = 0; i <= degP; i++) {
				if (i+j+k <= degP) {
					for (itr = 0; itr < Ns; itr++) {
						Pa[ind][itr] = pow(Xa[itr][1],(double)i)*pow(Xa[itr][2],(double)j)*pow(Xa[itr][3],k);					
					}
				}
			}
		}
	}
	
	double Xa_i[3],Xa_j[3];
	for (i = 0; i < Na; i++) {
		for (j = 0; j < Na; j++) {
			Xa_i[1] = Xa[i][1];
			Xa_i[2] = Xa[i][2];
			Xa_i[3] = Xa[i][3];
			Xa_j[1] = Xa[j][1];
			Xa_j[2] = Xa[j][2];
			Xa_j[3] = Xa[j][3];
			right[i][j] = CalcRBF(rbf,r,Xa_i,Xa_j);
		}
	}
	
	for (i = 0; i < Na; i++) {
		for (j = 0; j < pDim; j++) {
			Aas[i][j] = Pa[i][j];
		}
	}
	
	// calculate coupling matrix H
	if (degP == 0) {
		Aas = right;
		H = multiply(Aas,Mi);
	} else {
		// matcat_h stack matrices vertically
		F = matcat_v(multiply(Mp,multiply(Ps,Mi),Mi-multiply(Mi,multiply(Pst,multiply(Mp,multiply(Ps,Mi)))));		
		H = multiply(Aas,F);
	}	
	return H;
}





