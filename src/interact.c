#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double Factorial(int n) {
    int i;
    double res;

    if ((n > 0) || (n==0)){
       res = 1.0;
       for(i = 1; i <= n; i++) res *= i;
    }  else {
       printf("invalid n in Factorial...\n" );
       exit(1);
    }
    return(res);
}


/* 	cbintgs(pt,n1+n2+3,fct, b);  */ 
void cbintgs(double x, int k, double *fct, double *b)
{
	int i; 
	int i0; 
	int last;
	int m;  
	double absx; 
	double expx; 
	double expmx;
	double rx; 
	double di;  
	double y; 
	double ytrm; 
	
	absx = fabs(x); 
	i0 = 0; 
	
	if (absx > 3.0) { 
		expx = exp(x); 
		expmx = 1.0/expx; 
		rx = 1.0/x; 
		b[0] = (expx - expmx)*rx; 
		for (i = 1; i <= k; i++) { 
		       di =  (double)(i); 
		       b[i] =  (di*b[i-1]+ pow(-1.0,i)*expx - expmx)*rx ; 
		 }       
		 return; 
	} else 	if (absx > 2.0) { 
		if((k < 10) || (k == 10)) { 
			expx = exp(x);
			expmx = 1.0/expx;
			rx = 1.0/x;
			b[0] = (expx - expmx)*rx;
			for (i = 1; i <= k; i++) {  
			     di =  (double)(i);
			     b[i] = ( di*b[i-1]+ pow(-1.0, i)*expx - expmx)*rx;  
			}      
		} else {  	
			last = 15; 
			for (i = i0; i <= k; i++) { 
				y = 0.0; 
				for (m = i0; m <= last; m++) { 
					ytrm = pow(-x, m-1)*(1.0 - pow(-1.0, m+i+1))/(fct[m]*((double)(m+i+1))   ); 
      					y = y + ytrm*(-x); 
				} 	
				b[i] = y; 
			}	
		 }	
		 return; 
	} else if (absx > 1.0) { 
		if((k < 7) || (k == 7)) { 
			expx = exp(x); 
			expmx = 1.0/expx; 
			rx = 1.0/x; 
			b[0] = (expx - expmx)*rx; 
			for (i = 1; i <= k; i++)  { 
			    di = (double)(i); 
			    b[i] = (di*b[i-1]+pow(-1.0, i)*expx - expmx)*rx; 
			}    
		} else {  	
			last = 12; 
			for (i  = i0; i <= k; i++) { 
				y = 0.0; 
				for (m = i0; m <= last; m++) { 
					ytrm = pow(-x, m-1)*(1.0 - pow(-1.0, m+i+1))/(fct[m]*((double)(m+i+1))); 
      					y = y + ytrm*(-x); 
				} 	
				b[i] = y; 
			}	
		} 	
		return; 		
	} else if (absx > 0.5) { 
		if((k < 5) || (k == 5)) { 
			expx = exp(x); 
			expmx = 1.0/expx; 
			rx = 1.0/x; 
			b[0] = (expx - expmx)*rx; 
			for (i = 1; i <= k; i++) { 
	                       di = (double)(i);  		
			       b[i] = ( di*b[i-1]+pow(-1.0, i)*expx - expmx)*rx ; 
			 }      
		} else {  	
			last = 7; 
			for (i = i0; i <= k; i++) { 
				y = 0.0; 
				for (m = i0; m <= last; m++)  { 
					ytrm = pow(-x, m-1)*(1.0 - pow(-1.0, m+i+1))/(fct[m]*((double)(m+i+1))); 
      					y = y + ytrm*(-x); 
				} 	
				b[i] = y; 
			}	
		} 		
		return; 
	} else if (absx > 1.0e-8) { 
			last = 6; 
			for (i=i0; i <= k; i++) { 
				y = 0.0; 
				for (m=i0; m <= last; m++)  { 
					ytrm = pow(-x, m-1)*(1.0 - pow(-1.0, m+i+1))/(fct[m]*((double)(m+i+1))); 
      					y = y + ytrm*(-x); 
				}	
				b[i] = y; 
			} 	
			return;  
	} else  { 
	     for (i = i0; i <= k; i++)
	            b[i] = (1.0-pow(-1.0, i+1) )/((double)(i+1)); 
	} 
	return; 
}		

/* 	caintgs(p,n1+n2+3, a);  */ 
void caintgs(double x, int k,  double *a)
{
	int i; 
	double dexp; 
	dexp = exp(-x); 
	a[0] = dexp/x; 
	for (i = 1; i <= k; i++) a[i] = ( a[i-1]*((double)(i)) + dexp )/x ; 
	return; 
}	

/*  coeffs(n1, n2,  i1-1, fct);   */ 
double coeffs(int na, int nb, int k, double *fct) 
{
	int i; 
	int ia; 
	int ie; 
	int il;  
	int j; 
	int je; 
	int l; 
	double binm_na_i; 
	double binm_nb_j; 
	double res = 0.0; 
	
	l = na + nb - k; 
	if(l < na) ie = l + 1; 
	else ie = na+1; 
	if(l < nb) je = l; 
	else je = nb; 
	ia = l - je + 1; 
	for (il=ia; il <= ie; il++) {
		i = il - 1; 
		j = l - i; 
		binm_na_i = fct[na]/(fct[na-i]*fct[i]);
		binm_nb_j = fct[nb]/(fct[nb-j]*fct[j]); 
		res = res + binm_na_i*binm_nb_j*pow(-1.0, j); 
	}	
	return(res); 
} 	
	  
	
	  /* css(na2-1,nb2-index,z2ra,z2rb,r);   */
double css(int nn1, int nn2, double alpha, double beta,  double r)
{
	int i;  
	int i1;   
        int k; 
	int n1; 
	int n2;
	int nni1;
	int ulim; 
	double fct[30]; 
	double a[30]; 
	double b[30]; 
	double coff;
	double p;
	double pt; 
	double difzeta; 
	double sumzeta;
	double x; 
	double zeta1; 
	double zeta2; 

	for (i = 0; i < 30; i++) { 
	      fct[i] = Factorial(i); 
	      a[i] = 0.0; 
	      b[i] = 0.0; 
	}   
	   
	n1 = nn1; 
	n2 = nn2; 
	p  = (alpha + beta)*0.5; 
	pt = (alpha - beta)*0.5; 
	zeta1 = alpha/r;
	zeta2 = beta/r;
	sumzeta = zeta1+zeta2;
	difzeta = zeta1-zeta2;

	if (n2 < n1) { 
		k  =  n1;
		n1 =  n2;
		n2 =  k;
		pt = (-1.0)*pt;
		difzeta = (-1.0)*difzeta;
	}	

	if ((p > 86.0) || (pt > 86.0)) return(0.0); 
/****************************
 Find a and b integrals    * 
 ***************************/
	
	caintgs(p,n1+n2+3, a); 
	cbintgs(pt,n1+n2+3,fct, b); 
	ulim = n1+n2+1; 
	x  =  0.0; 
	for (i1= 1; i1 <= ulim; i1++) { 
		nni1 = n1 + n2 - i1 + 2; 
		coff = coeffs(n1, n2,  i1-1, fct); 
		x = x + coff*a[i1-1]*b[nni1-1]; 
	}	
	return(0.5*x); 
}	

/* orb_interaction(iatom_ni, zeta_i, 	jatom_nj, zeta_j, dr); */
double orb_interaction(int na, double za, int nb,  double zb, double r) 
{ 
	int index;
	int i;   
	int na2;
	int nb2; 
	double ctrm; 
	double d2rtrm; 
	double drtrm; 
	double halfr; 
	double rtrm; 
	double ss; 
	double trm1;
	double trm2;
	double trm3; 
	double z2ra; 
	double z2rb; 
	double ztrm;
	
	double gam = 0.0 ; 
	if(fabs(r) < 1.0e-8) return(0.0); 
	
  	z2ra = 2.0*za*r; 
  	z2rb = 2.0*zb*r; 
  	na2 = 2*na; 
  	nb2 = 2*nb; 
  	halfr = 0.5*r;
  	d2rtrm = pow(halfr, na2-2); 
  	drtrm = d2rtrm*halfr; 
  	rtrm = drtrm*halfr; 
	ss =  css(na2-1, 0, z2ra, 0.0, r); 
	gam = rtrm*ss; 

	rtrm = drtrm; 
	drtrm = d2rtrm; 
        ztrm = 0.5/( zb*((double)(nb2))  ); 	
	index = nb2; 
	for (i = 1; i <= nb2; i++) { 
		rtrm = rtrm*halfr; 
		drtrm = drtrm*halfr; 
		ztrm = ztrm*2.0*zb; 
		ctrm = ztrm/Factorial(nb2-index); 
                ss = css(na2-1,nb2-index,z2ra,z2rb,r);  
		trm1 = ctrm*((double)(index)); 
		trm2 = trm1*rtrm;  
		gam = gam - trm2*ss; 		 
		index = index - 1; 
        } 
	
	trm3 = pow(2.0*za, na2+1)/Factorial(na2); 
	gam = gam*trm3; 	
	return(gam) ; 
}	
