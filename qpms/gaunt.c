#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#define ZERO_THRESHOLD 1.e-8
#define BF_PREC 1.e-14
//  "Besides, the determined Real Programmer can write FORTRAN programs in any language."
//  -- Ed Post, Real Programmers Don't Use Pascal, 1982.

// logarithm of factorial (from basicsubs.f90)
double lnf (double z) {
	// expansion parameters
	static const double v_c0[] = {
		0.16427423239836267e5, -0.48589401600331902e5, 0.55557391003815523e5, -0.30964901015912058e5,
	       	0.87287202992571788e4, -0.11714474574532352e4, 0.63103078123601037e2, -0.93060589791758878e0, 
	       	0.13919002438227877e-2,-0.45006835613027859e-8, 0.13069587914063262e-9
	};
	static const double cp =  2.5066282746310005;
	double a = 1.; 
	double b = z + 10.5;
	b = (z + 0.5) * log(b) - b;
	for (int i = 0; i < (sizeof(v_c0) / sizeof(double)); i++) {
		z += 1.;
		a += v_c0[i] / z;
	}

	return b+log(cp*a);
}

// logarithm of Pochhammer function (from basicsubs.f90)
// FIXME replace with standard C99 lgamma functions
double lpoch(double x, double n) {
	if(fabs(n) < 1e-5) // ???
		return 1.;
	double sum = x+n;
	return lnf(sum-1.) - lnf(x-1.);
}

// pochhammer function: substitute for fortran DPOCH
// gamma(a+x) / gamma(a)
double poch(double a, double x) { // FIXME replace with standard C99 lgamma functions
	return exp(lpoch(a,x));
}

double f_a0 (int m, int n, int mu, int nu) {
	double logw = lnf(n+nu-m-mu) - lnf(n-m) - lnf(nu-mu);
	double logp = lpoch(n+1, n) + lpoch(nu+1, nu) - lpoch(n+nu+1, n+nu);
	return exp(logw+logp);
}

// coefficient Ap(m,n,mu,nu,p) (from vec_trans.f90)
int f_Ap(int m, int n, int mu, int nu, int p) {
	return p*(p-1)*(m-mu)-(m+mu)*(n-nu)*(n+nu+1);
}

// coefficient a(m,n,mu,nu,1) normalized by the backward recursion (from vec_trans.f90)

double f_a1norm(int m, int n, int mu, int nu) {
	int n4 = n + nu - m - mu;
	return ((2.*n + 2.*nu - 3.) / 2.) * (1. - ((2.*n + 2.*nu - 1.) / (n4 * (n4-1.))) 
			* ((m - n) * (m - n + 1.) / (2.*n - 1.) + (mu-nu) * (mu-nu+1.)/(2.*nu-1.)));
}

// coefficient a(m,n,mu,nu,2) normalized by the backward recursion (From vec_trans.f90)

double f_a2norm(int m, int n, int mu, int nu) {
	double n4 = n + nu - m - mu;
	double n2nu = 2*n + 2*nu;
	double mn = m - n;
	double munu = mu - nu;

	return ((n2nu-1.)*(n2nu-7.)/4.) * \
                 ( ((n2nu-3.)/(n4*(n4-1.))) * \
                 ( ((n2nu-5.)/(2.*(n4-2.)*(n4-3.))) * \
                 ( mn*(mn+1.)*(mn+2.)*(mn+3.)/((2.*n-1.)*(2.*n-3.)) + \
                   2.*mn*(mn+1.)*munu*(munu+1.)/((2.*n-1.)*(2.*nu-1.)) + \
                   munu*(munu+1.)*(munu+2.)*(munu+3.)/((2.*nu-1.)*(2.*nu-3.)) \
                   ) - mn*(mn+1.)/(2.*n-1.) - munu*(munu+1.)/(2.*nu-1.) ) +0.5);
}

// just for convenience – square of ints
static inline int isq(int x) {return x * x;}

double f_alpha(int n, int nu, int p) {
	return (isq(p) - isq(n-nu))*(isq(p)-isq(n+nu+1)) / (double)(4*isq(p)-1);
}

static inline int min1pow(int pow) { return (pow % 2) ? -1 : 1; }

// starting value of coefficient a(m,n,mu,nu,qmax) for the forward recursion
double f_aqmax(int m, int n, int mu, int nu, int qmax) {
	int pmin = n + nu - 2*qmax;
	// TODO? zde měl int a double varianty téhož – má to nějaký smysl???
	if (pmin == n-nu) {
		double logw = lnf(n+m) + lnf(2*pmin+1) - lnf(nu-mu) - lnf(n-nu) - lnf(pmin+m+mu);
		double logp = lpoch(nu+1, nu) - lpoch(n+1, n+1);
		return min1pow(mu)*exp(logw+logp);
	} else if (pmin == nu-n) {
		double logw = lnf(nu+mu) + lnf(2*pmin+1) - lnf(n-m) - lnf(nu-n) - lnf(pmin+m+mu);
		double logp = lpoch(n+1, n) - lpoch(nu+1, nu+1); // ??? nešel by druhý lpoch nahradit něčím rozumnějším?
		return min1pow(m)*exp(logw+logp);
	} else if (pmin == m+mu) {
        	double logw = lpoch(qmax+1,qmax)+lnf(n+nu-qmax)+lnf(n+m)+lnf(nu+mu) \
                	-lnf(n-qmax)-lnf(nu-qmax)-lnf(n-m)-lnf(nu-mu)-lnf(n+nu+pmin+1);
        	return min1pow(n+m-qmax)*(2*pmin+1)*exp(logw);
	} else if (pmin == -m-mu) {
	        double logw = lpoch(qmax+1,qmax)+lnf(n+nu-qmax)+lnf(pmin-m-mu) \
        	        -lnf(n-qmax)-lnf(nu-qmax)-lnf(n+nu+pmin+1);
		return min1pow(nu+mu-qmax)*(2*pmin+1)*exp(logw);
	} else if (pmin==m+mu+1) {
		int Apmin = f_Ap(m,n,mu,nu,pmin);
		double logw = lpoch(qmax+1,qmax)+lnf(n+nu-qmax)+lnf(n+m)+lnf(nu+mu) \
	                -lnf(n+nu+pmin+1)-lnf(n-qmax)-lnf(nu-qmax)-lnf(n-m)-lnf(nu-mu);
        	return min1pow(n+m-qmax)*Apmin*(2*pmin+1)*exp(logw)/(double)(pmin-1);
	} else if (pmin==-m-mu+1) {
		int Apmin=f_Ap(m,n,mu,nu,pmin);
		double logw=lpoch(qmax+1,qmax)+lnf(n+nu-qmax)+lnf(pmin-m-mu) \
	                -lnf(n+nu+pmin+1)-lnf(n-qmax)-lnf(nu-qmax);
		return min1pow(nu+mu-qmax)*Apmin*(2*pmin+1)*exp(logw) / (double)(pmin-1);
	}
	assert(0);
}

// coeff a(m,n,mu,nu,2) normalizzato per la backward recursion calcolato questa volta per ricorsione
// (from vec_trans.f90)
double f_a2normr(int m, int n, int mu, int nu, double a1norm) {
	int p = n + nu - 4;
	int p1 = p - m - mu;
	int p2 = p + m + mu;
	int Ap4 = f_Ap(m,n,mu,nu,p+4);
	int Ap6 = f_Ap(m,n,mu,nu,p+6);
	
	double alphap1=f_alpha(n,nu,p+1);
	double alphap2=f_alpha(n,nu,p+2);
	
	if (!Ap4) {
		if(!Ap6) {
	
	                double c0=(p+2)*(p1+1)*alphap1;
	                double c1=(p+1)*(p2+2)*alphap2;
	
	                return (c1/c0)*a1norm;
		} else /* Ap6 != 0 */ {
	                int Ap2=f_Ap(m,n,mu,nu,p+2);
	                int Ap3=f_Ap(m,n,mu,nu,p+3);
	                int Ap5=f_Ap(m,n,mu,nu,p+5);
	                double alphap5=f_alpha(n,nu,p+5);
	
	                double c0=(p+2)*(p+3)*(p+5)*(p1+1)*(p1+2)*(p1+4)*Ap6*alphap1;
	                double c1=(p+5)*(p1+4)*Ap6*(Ap2*Ap3+(p+1)*(p+3)*(p1+2)*(p2+2)*alphap2);
	                double c2=(p+2)*(p2+3)*Ap2*(Ap5*Ap6+(p+4)*(p+6)*(p1+5)*(p2+5)*alphap5);
	
	                return (c1/c0)*a1norm+(c2/c0);
		}
	} else /* Ap4 != 0 */ {
	                int Ap2=f_Ap(m,n,mu,nu,p+2);
	                int Ap3=f_Ap(m,n,mu,nu,p+3);
	                double alphap3=f_alpha(n,nu,p+3);
	                double alphap4=f_alpha(n,nu,p+4);
	
	                double c0=(p+2)*(p+3)*(p1+1)*(p1+2)*Ap4*alphap1;
	                double c1=Ap2*Ap3*Ap4+(p+1)*(p+3)*(p1+2)*(p2+2)*Ap4*alphap2 \
	                 + (p+2)*(p+4)*(p1+3)*(p2+3)*Ap2*alphap3;
	                double c2=-(p+2)*(p+3)*(p2+3)*(p2+4)*Ap2*alphap4;
	                
			return (c1/c0)*a1norm+(c2/c0);
	}
}

#define MAX(x,y) (((x) > (y)) ? (x) : (y))

int j3minimum(int j1, int j2, int m1, int m2) {
	return MAX(abs(j1-j2), abs(m1+m2));
}

int nw(int j1, int j2, int m1, int m2) {
	return j1+j2+1-MAX(abs(j1-j2),abs(m1+m2));
}

double wdown(int j1, int j2, int m1, int m2) {
	double logw = .5*(lnf(2.*j1,
TODO TODO TODO ZDE JSEM SKONČIL

//******************************************************************************
//7) subroutine Wigner3jm: calcolo il vettore di simboli 3jm
//******************************************************************************

void  wigner3jm(int j1, int j2, int m1, int m2, int j3min, int j3max, double *v_w3jm){
	// in the original code, the dimension of v_w3jm is (j3min:j3max).
	// In C, this means it has length j3max-j3min+1, and we must
	// always deduct j3min from the indices
	
	// Inizializzo gli indici per la downward recursion
	int j3 = j3max; // we don't use separate j3int as gevero does.
	
	// In questo if separo i casi in cui ho un vettore di lunghezza uno da quelli che
	// necessitano dell'uso della ricorsione
	if (j3min==j3max) // big_if
		v_w3jm[j3max-j3min]=wdown0(j1,j2,m1,m2); // Unico termine da calcolare
	else {
		// Si inizializza la ricorsione
		v_w3jm[j3max-j3min]=wdown0(j1,j2,m1,m2);
		v_w3jm[j3max-1-j3min]=-(dr(j1,j2,j1+j2,m1,m2,-m1-m2)/( (j1+j2+1)*cr(j1,j2,j1+j2,m1,m2,-m1-m2) ))*v_w3jm[j3max-j3min];
	
	
		// Ciclo per il calcolo ricorsivo
		while(j3-2>=j3min){ // down_do
	
				//Primo coeff della ricorsione
				cd1=dr(j1,j2,j3-1,m1,m2,-m1-m2)/(j3*cr(j1,j2,j3-1,m1,m2,-m1-m2))
				cd2=((j3-1)*cr(j1,j2,j3,m1,m2,-m1-m2))/(j3*cr(j1,j2,j3-1,m1,m2,-m1-m2))
				//Ricorsione
				v_w3jm[j3int-2-j3min]=-cd1*v_w3jm[j3int-1-j3min]-cd2*v_w3jm[j3int-j3min]
	
				//Aggiorno gli indici
				--j3;
		} //END DO down_do
	
		// Inizializzo gli indici per la upward recursion
		j3int=j3min
		j3=REAL(j3min,dbl)
	
		// Calcolo del primo termine di wigner dal basso
		v_w3jm[j3int-j3min]=wup0(j1,j2,m1,m2)
	
		// Calcolo del secondo termine di wigner dal basso
		// Pongo anche la condizione sul coefficienti nel caso ci sia signolarita'
		cu3_if: IF (j3min==0) THEN
				cu3=0
		ELSE
				cu3=dr(j1,j2,j3,m1,m2,-m1-m2)/(j3*cr(j1,j2,j3+1,m1,m2,-m1-m2))
		END IF cu3_if
	
		w3jm_temp=-cu3*v_w3jm[j3int-j3min]
	
		// If legato alla monotonia della successione
		up_if: IF (ABS(w3jm_temp)>ABS(v_w3jm[j3min-j3min])) THEN
	
				// Aggiorno gli indici e metto nell'array il secondo valore
				// in questo modo sono pronto per iniziale la upward recursion
				// a tre termini
				j3int=j3int+1
				v_w3jm[j3int-j3min]=w3jm_temp
	
				up_do: DO
						//Aggiorno gli indici
						j3int=j3int+1
						j3=REAL(j3int,dbl)
	
						IF (j3int-1==j3max) EXIT
	
	// 					IF ((INT(-m1)==-1).AND.(INT(j1)==1).AND.(INT(m2)==1).AND.(INT(j2)==2)) THEN
	// 					WRITE(*,*) "j3-1,cr1,cr2",j3-1,cr(j1,j2,j3,m1,m2,-m1-m2),cr(j1,j2,j3,m1,m2,-m1-m2)
	// 					END IF
	
						//Primo e secondo coeff della ricorsione
						cu1=dr(j1,j2,j3-1,m1,m2,-m1-m2)/((j3-1)*cr(j1,j2,j3,m1,m2,-m1-m2))
						cu2=(j3*cr(j1,j2,j3-1,m1,m2,-m1-m2))/((j3-1)*cr(j1,j2,j3,m1,m2,-m1-m2))
	
						//Assegnazione temporanea della ricorsione
						w3jm_temp=-cu1*v_w3jm[j3int-1-j3min]-cu2*v_w3jm[j3int-2-j3min]
	
						IF ((ABS(w3jm_temp)<ABS(v_w3jm[j3int-1-j3min]))  .OR. ((j3int-1)==j3max) ) EXIT // Cond. di uscita
	
						v_w3jm[j3int-j3min]=w3jm_temp	//Assegno perche' e' ok
	
				END DO up_do
	
		END IF up_if
	
	} // big_if
	
	END SUBROUTINE wigner3jm

void gaunt_cz(int m, int n, int mu, int n, int qmaxa, double *v_aq, int *error) {
	*error = 0;
	if (abs(m) > n || abs(mu) > nu) { // error_if
		*error=1;
		assert(0);
		return;
	}

	// calcolo i bounds dei vettori di wigner
	int pmin = j3minimum(n,nu,m,mu);
	int pmax = n+nu;
	int pmin0 = j3minimum(n,nu,0,0);
	// Alloco i vettori di wigner e li calcolo
	ALLOCATE(v_w3jm(pmin:pmax),v_w3jm0(pmin0:pmax),STAT=stat_a)
	CALL wigner3jm(n,nu,m,mu,pmin,pmax,v_w3jm)
	CALL wigner3jm(n,nu,0.0D0,0.0D0,pmin0,pmax,v_w3jm0)
	
	// Entro nel ciclo per il calcolo dei coefficienti di gaunt
	gaunt_do: DO q=0,qmaxa
	
				// Calcolo dell'indice p, sia reale che intero
				p=INT(n+nu,lo)-2*q
				pr=REAL(p,dbl)
	
				//Calcolo del fattoriale
				logw = 0.5D0* (lnf(n+m,r)+lnf(nu+mu,r)+lnf(pr-m-mu,r) - &
					 		   lnf(n-m,r)-lnf(nu-mu,r)-lnf(pr+m+mu,r))
				fac= EXP(logw)
	
				// Calcolo del coefficiente di gaunt
				v_aq(q)=((-1.0D0)**INT(m+mu,lo))*(2.0D0*pr+1.0D0)*fac*v_w3jm(p)*v_w3jm0(p)
	
	END DO gaunt_do
	
	// Disalloco i vettori di wigner a lavoro finito
	DEALLOCATE(v_w3jm,v_w3jm0)
	
	END SUBROUTINE gaunt_cz


// gaunt_xu from vec_trans.f90
// FIXME set some sensible return value
void gaunt_xu(int m, int n, int mu, int nu, int qmax, double *v_aq, int *error) {
	
	*error = 0;
	int v_zero[qmax];
	for (int i = 0; i < qmax; i++) v_zero[i] = 1;
	double v_aq_cz[qmax];
	for (int i = 0; i < qmax; i++) v_aq_cz[i] = 0.;
	int qi = 0;

	if(abs(m)>n || abs(mu)>nu) {
		*error = 1;
		fprintf(stderr, "invalid values for m, n, mu or nu\n");
		return; // FIXME vyřešit chyby
	}

	switch(qmax) { //qmax_case
		case 0:
			v_aq[0]  = f_a0(m,n,mu,nu);
			break;
		case 1:
			v_aq[0] = f_a0(m,n,mu,nu);
			v_aq[1] = v_aq[0] + f_a1norm(m,n,mu,nu);
			// controllo gli zeri
			if (fabs(v_aq[1]/v_aq[0]) < ZERO_THRESHOLD) {
				v_aq[1] = 0.;
				v_zero[1] = 0;
			}
			break;
		case 2:
			v_aq[0] = f_a0(m,n,mu,nu);
			
			v_aq[1] = v_aq[0] + f_a1norm(m,n,mu,nu);
			// controllo gli zeri
			if (fabs(v_aq[1]/v_aq[0]) < ZERO_THRESHOLD) {
				v_aq[1] = 0.;
				v_zero[1] = 0;
			}
			
			v_aq[2] = v_aq[0] * f_a2normr(m,n,mu,nu,v_aq[1]/v_aq[0]);
			// controllo gli zeri
			if (fabs(v_aq[2]/v_aq[0]) < ZERO_THRESHOLD) {
				v_aq[2] = 0.;
				v_zero[2] = 0;
			}
			break;
		default:
			if (m == 0 && mu == 0) { // big_if
				v_aq[0] = f_a0(m,n,mu,nu);

				// backward recursion
				for (int q = 1; q <= qmax; ++q) { // uno_b_do
					int p = n + nu - 2*q;
					double c0 = f_alpha(n,nu,p+1);
					double c1 = f_alpha(n,nu,p+2);
					v_aq[q] = (c1/c0) * v_aq[q-1];

					// Vedo se il q-esimo valore e' zero
					if (v_zero[q-1] == 1) {// v_zero_if_1
						if(fabs(v_aq[q]/v_aq[q-1]) < ZERO_THRESHOLD) { // zg_if_1
							v_aq[q] = 0.;
							v_zero[q] = 0;
						}
					} else if (v_zero[q-1]==0 && v_zero[q-2]) { 
						if(fabs(v_aq[q]/v_aq[q-2]) < ZERO_THRESHOLD) {// zg_if1_1
							v_aq[q] = 0.;
							v_zero[q] = 0;
						}
					} //v_zero_if_1
				} //uno_b_do
			} else if (mu == m && nu == n) {
				v_aq[0] = f_a0(m,n,mu,nu);

				// backward recursion
				for (int q = 1; q <= qmax; ++q) { // due_b_do
					// calculate pre-coefficients
					int p = n + nu - 2*q;
					int p1 = p - m - mu;
					int p2 = p + m + mu;

					// calculate recursion coefficients
					double c0 = (p+2) * (p1+1) * f_alpha(n,nu,p+1);
					double c1 = (p+1) * (p2+2) * f_alpha(n,nu,p+2);

					//recursion
					v_aq[q] = (c1/c0) * v_aq[q-1];

					// Vedo se il q-esimo valore e' zero
					if (v_zero[q-1] == 1) {// v_zero_if_2
						if(fabs(v_aq[q]/v_aq[q-1]) < ZERO_THRESHOLD) { // zg_if_2
							v_aq[q] = 0.;
							v_zero[q] = 0;
						}
					} else if (v_zero[q-1]==0 && v_zero[q-2]) { 
						if(fabs(v_aq[q]/v_aq[q-2]) < ZERO_THRESHOLD) {// zg_if1_2
							v_aq[q] = 0.;
							v_zero[q] = 0;
						}
					} //v_zero_if_2
				} // due_b_do
			} else if (mu == -m) {
				// Primo valore per la backward recursion
				v_aq[0] = f_a0(m,n,mu,nu);
				v_aq[1] = f_a1norm(m,n,mu,nu)*v_aq[0];

				// Controllo gli zeri
				if (fabs(v_aq[1]/v_aq[0]) < ZERO_THRESHOLD) { //zg_if_3_0
					v_aq[1] = 0.;
					v_zero[1] = 0;
				} //zg_if_3_0

				// backward recursion
				for (int q = 2; q <= qmax; ++q) { // tre_b_do
					// calculate pre-coefficient
					int p = n + nu - 2*q;

					// calculate recursion coefficients
					double c0 = f_alpha(n, nu, p+1);
					double c1 = 4*isq(m) + f_alpha(n,nu,p+2) + f_alpha(n,nu,p+3);
					double c2 = - f_alpha(n, nu, p+4);

					// recursion
					v_aq[q] = (c1/c0)*v_aq[q-1] + (c2/c0)*v_aq[q-2];
					
					// Vedo se il q-esimo valore e' zero
					if (v_zero[q-1] == 1) {// v_zero_if_3
						if(fabs(v_aq[q]/v_aq[q-1]) < ZERO_THRESHOLD) { // zg_if_3
							v_aq[q] = 0.;
							v_zero[q] = 0;
						}
					} else if (v_zero[q-1]==0 && v_zero[q-2]) { 
						if(fabs(v_aq[q]/v_aq[q-2]) < ZERO_THRESHOLD) {// zg_if1_3
							v_aq[q] = 0.;
							v_zero[q] = 0;
						}
					} //v_zero_if_3
				} // tre_b_do
				
				// forward recursion	
		                // Primo valore per la forward recursion,errore relativo e suo swap
		                double aq_fwd=f_aqmax(m,n,mu,nu,qmax);
		                double res=fabs(aq_fwd-v_aq[qmax])/fabs(aq_fwd);
		
		                //Se non ho precisione, sostituisco i valori
		                if (res>BF_PREC) { //tre_f_if
		                        v_aq[qmax]=aq_fwd;
		                        int qi=1;
					int zeroswitch = 0; // black magic (gevero's "switch")
		                        //Entro nel ciclo della sostituzione valori
		                        for( int q=qmax-1;q>=0;--q) { // tre_f_do
						switch(qmax-q) {// tre_q_case // FIXME switch -> if/else
							case 1: {// q==qmax-1
		                                        	//Calcolo v_aq[qmax-1]
			                                        int p=n+nu-2*(q+2);
			                                        double c1=4*isq(m)+f_alpha(n,nu,p+2)+f_alpha(n,nu,p+3);
		                                        	double c2=-f_alpha(n,nu,p+4);
			                                        double aq_fwd=-(c1/c2)*v_aq[qmax];
		
		                                        	switch(v_zero[q]) { // z_3_1_case
									case 0:
										v_aq[q] = 0.;
										break;
									case 1:
										res=fabs(aq_fwd-v_aq[q])/fabs(aq_fwd);
										break;
									default:
										assert(0);
								}
								}	
								break;
							default: { //Per tutti gli altri q
		                                       	//Calcolo v_aq[qmax-1]
		                                        int p=n+nu-2*(q+2);
		                                        double c0=f_alpha(n,nu,p+1);
		                                        double c1=4*isq(m)+f_alpha(n,nu,p+2)+f_alpha(n,nu,p+3);
		                                        double c2=-f_alpha(n,nu,p+4);
		                                        aq_fwd=-(c1/c2)*v_aq[q+1]+(c0/c2)*v_aq[q+2];

							switch(v_zero[q]){ // z_3_2_case//FIXME switch -> if/else
								case 0:
									v_aq[q] = 0.;
									break;
		                                        	case 1: //Il valore precedente e' zero
			                                                res=fabs(aq_fwd-v_aq[q])/fabs(aq_fwd);
									break;
								default:
									assert(0);
							}
							}
						} //tre_q_case
		                        	//Adesso se la precisione e' raggiunta esco dal ciclo, 
						//se no sostituisco e rimango
		                        	if (res<BF_PREC || q==0 || fabs(aq_fwd) < fabs(v_aq[q+1])) 
							break; //tre_f_do
		                        	//Sono nel ciclo, allora sostituisco eaggiorno indice e residuo
			                        v_aq[q]=aq_fwd;
		                        	qi=q;
						assert(q); // assert níže přesunut sem
					} // tre_f_do
		                // Check sul ciclo di sostituzione
				//assert(q);
				/*
		                error_if1: IF (q==0) THEN
		                                        WRITE(*,*)
		                                        WRITE(*,*) "Si e' verificato un errore nella subroutine gaunt_xu:"
		                                        WRITE(*,*) "la precisione richiesta per i coefficienti di Gaunt nella backward"
		                                        WRITE(*,*) "e forward recursion non e' stata raggiunta"
		                                        WRITE(*,*)
		                                        error=1
		                                        RETURN
		                END IF error_if1
				*/
				} // tre_f_if
			} else { // caso generale (4)
				// backward
				// Calcolo direttamente i primi tre valori della ricorsione
				v_aq[0]=f_a0(m,n,mu,nu);
				v_aq[1]=v_aq[0]*f_a1norm(m,n,mu,nu);

				// vedo se il secondo valore e' zero
				if (fabs(v_aq[1]/v_aq[0]) < ZERO_THRESHOLD) { // zg1_if
					v_aq[1] = 0.;
					v_zero[1] = 0;
				}

				//...........................................................
				//Calcolo il terzo valore della ricorsione in funzione di Ap4
				//...........................................................
				//Inizializzo i valori comuni per i coefficienti
				int p=n+nu-2*(2);
				int p1=p-m-mu;
				int p2=p+m+mu;
				double alphap1=f_alpha(n,nu,p+1);
				double alphap2=f_alpha(n,nu,p+2);
				int Ap2=f_Ap(m,n,mu,nu,p+2);
				int Ap3=f_Ap(m,n,mu,nu,p+3);
				int Ap4=f_Ap(m,n,mu,nu,p+4);
		
				//Con questo if decido se mi serve la ricorsione a 3 o 4 termini
				if (Ap4==0) { //Ap4_2_if
					//Calcolo i restanti valori preliminari
					int Ap5=f_Ap(m,n,mu,nu,p+5);
					int Ap6=f_Ap(m,n,mu,nu,p+6);
					double alphap5=f_alpha(n,nu,p+5);
					double alphap6=f_alpha(n,nu,p+6);
		
					//Calcolo i coefficienti per la ricorsione ma non c3 perche' qui e solo qui non mi serve
					double c0=(p+2)*(p+3)*(p+5)*(p1+1)*(p1+2)*(p1+4)*Ap6*alphap1;
					double c1=(p+5)*(p1+4)*Ap6*(Ap2*Ap3 + (p+1)*(p+3)*(p1+2)*(p2+2)*alphap2);
					double c2=(p+2)*(p2+3)*Ap2*(Ap5*Ap6 + (p+4)*(p+6)*(p1+5)*(p2+5)*alphap5);
		
					//Calcolo il mio coefficiente
					v_aq[2]=(c1/c0)*v_aq[1]+(c2/c0)*v_aq[0];
		
					//Assegno l'indice segnaposto per Ap4=0
					// q4=2 FIXME UNUSED
				} else {
					//Calcolo i restanti valori preliminari
					double alphap3=f_alpha(n,nu,p+3);
					double alphap4=f_alpha(n,nu,p+4);
		
					//Calcolo coefficienti ricorsione
					double c0=(p+2)*(p+3)*(p1+1)*(p1+2)*Ap4*alphap1;
					double c1=Ap2*Ap3*Ap4+(p+1)*(p+3)*(p1+2)*(p2+2)*Ap4*alphap2+ \
					  (p+2)*(p+4)*(p1+3)*(p2+3)*Ap2*alphap3;
					double c2=-(p+2)*(p+3)*(p2+3)*(p2+4)*Ap2*alphap4;
		
					//Calcolo il mio coefficiente
					v_aq[2]=(c1/c0)*v_aq[1]+(c2/c0)*v_aq[0];
				} // Ap4_2_if
		
				//Vedo se il terzo valore e' zero
				if (v_zero[1]==1) { // v_zero_if1
					 if (fabs(v_aq[2]/v_aq[1])< ZERO_THRESHOLD) { //zg2_if
						v_aq[2]=0;
						v_zero[2]=0;
					 }
				} else if (v_zero[1]==0) {
					 if (fabs(v_aq[2]/v_aq[0])<ZERO_THRESHOLD) { //zg2_if1:
						v_aq[2]=0;
						v_zero[2]=0;
					 }
				} // v_zero_if1
		
		
				//...........................................................
				//Calcolo i restanti valori nel loop
				//...........................................................
				for (int q = 3; q <= qmax; q++ ) { //gen_bwd_do: DO q=3,qmax
		
					//Inizializzo i valori comuni per i coefficienti
					int p=n+nu-2*(q);
					int p1=p-m-mu;
					int p2=p+m+mu;
					double alphap1=f_alpha(n,nu,p+1);
					double alphap2=f_alpha(n,nu,p+2);
					int Ap2=f_Ap(m,n,mu,nu,p+2);
					int Ap3=f_Ap(m,n,mu,nu,p+3);
					int Ap4=f_Ap(m,n,mu,nu,p+4);
		
					//Con questo if decido se mi serve la ricorsione a 3 o 4 termini
					if (Ap4==0) { // Ap4_bwd_if: IF (Ap4==0) THEN
		
						//Calcolo i restanti valori preliminari
						int Ap5=f_Ap(m,n,mu,nu,p+5);
						int Ap6=f_Ap(m,n,mu,nu,p+6);
						double alphap5=f_alpha(n,nu,p+5);
						double alphap6=f_alpha(n,nu,p+6);
		
						//Calcolo i coefficienti per la ricorsione ma non c3 perche' qui e solo qui non mi serve
						double c0=(p+2)*(p+3)*(p+5)*(p1+1)*(p1+2)*(p1+4)*Ap6*alphap1;
						double c1=(p+5)*(p1+4)*Ap6*(Ap2*Ap3 + (p+1)*(p+3)*(p1+2)*(p2+2)*alphap2);
						double c2=(p+2)*(p2+3)*Ap2*(Ap5*Ap6 + (p+4)*(p+6)*(p1+5)*(p2+5)*alphap5);
						double c3=-(p+2)*(p+4)*(p+5)*(p2+3)*(p2+5)*(p2+6)*Ap2*alphap6;
		
						//Calcolo il mio coefficiente
						v_aq[q]=(c1/c0)*v_aq[q-1]+(c2/c0)*v_aq[q-2]+(c3/c0)*v_aq[q-3];
		
						//Assegno l'indice segnaposto per Ap4=0
						//q4=q // FIXME nepoužitá proměnná
					} else {
						//Calcolo i restanti valori preliminari
						double alphap3=f_alpha(n,nu,p+3);
						double alphap4=f_alpha(n,nu,p+4);
		
						//Calcolo coefficienti ricorsione
						double c0=(p+2)*(p+3)*(p1+1)*(p1+2)*Ap4*alphap1;
						double c1=Ap2*Ap3*Ap4+(p+1)*(p+3)*(p1+2)*(p2+2)*Ap4*alphap2+ \
						  (p+2)*(p+4)*(p1+3)*(p2+3)*Ap2*alphap3;
						double c2=-(p+2)*(p+3)*(p2+3)*(p2+4)*Ap2*alphap4;
		
						//Calcolo il mio coefficiente
						v_aq[q]=(c1/c0)*v_aq[q-1]+(c2/c0)*v_aq[q-2];
					} // END IF Ap4_bwd_if
		
					//Vedo se il q-esimo valore e' zero
					if(v_zero[q-1]==1) { //v_zero_ifq: IF (v_zero[q-1]==1) THEN
						if(fabs(v_aq[q]/v_aq[q-1])<ZERO_THRESHOLD) {//zgq_if
							v_aq[q]=0.;
							v_zero[q]=0;
						}
					} else if ((v_zero[q-1]==0) && (v_zero[q-2] !=0)) {
						 if (fabs(v_aq[q]/v_aq[q-2])<ZERO_THRESHOLD) {//zgq_if1:
							v_aq[q]=0.;
							v_zero[q]=0;
						 } // zgq_if1
					} // v_zero_ifq
		
				} //gen_bwd_do
		
				//---------------------------------------------------------------------------------
				//FORWARD
				//---------------------------------------------------------------------------------
		
				//Calcolo pmin,Apmin e la mia variabile logica
				int pmin=n+nu-2*(qmax);
				int Apmin=f_Ap(m,n,mu,nu,pmin);
				int test=(((Apmin)==0) && 
				    (((pmin)==(m+mu+1)) || ((pmin)==(-m-mu+1))));
		
				//........................................................
				//Se la mia variabile logica e' vera, Faccio il mio conto
				//........................................................
				if(test) { //Apmin_if: if (test) THEN
		
					//Il valore per qmax allora e' zero
					v_aq[qmax]=0;
		
					//Calcolo il secondo valore, e se la precisione e' raggiunta esco
					double aq_fwd=f_aqmax_1(m,n,mu,nu,qmax);
					double res=fabs(aq_fwd-v_aq[qmax-1])/fabs(aq_fwd);
					if (res<BF_PREC)
						return;
		
					//Assegno il secondo valore e faccio il ciclo
					v_aq[qmax-1]=aq_fwd;
					qi=1; //FIXME nepoužitá proměnná
		
					for (int q = qmax; q >= 2; --q) { //Apmin_do: DO q=qmax,2,-1
		
						//Calcolo pre-coefficienti
						int p=n+nu-2*(q);
						int p1=p-m-mu;
						int p2=p+m+mu;
						double alphap1=f_alpha(n,nu,p+1);
						double alphap2=f_alpha(n,nu,p+2);
						double alphap3=f_alpha(n,nu,p+3);
						double alphap4=f_alpha(n,nu,p+4);
						int Ap2=f_Ap(m,n,mu,nu,p+2);
						int Ap3=f_Ap(m,n,mu,nu,p+3);
						int Ap4=f_Ap(m,n,mu,nu,p+4);
		
						//Calcolo coefficienti ricorsione
						double c0=(p+2)*(p+3)*(p1+1)*(p1+2)*Ap4*alphap1;
						double c1=Ap2*Ap3*Ap4+(p+1)*(p+3)*(p1+2)*(p2+2)*Ap4*alphap2+ 
						  (p+2)*(p+4)*(p1+3)*(p2+3)*Ap2*alphap3;
						double c2=-(p+2)*(p+3)*(p2+3)*(p2+4)*Ap2*alphap4;
		
						//Ricorsione e residuo
						aq_fwd=-(c1/c2)*v_aq[q-1]+(c0/c2)*v_aq[q];
						res=fabs(aq_fwd-v_aq[q-2])/fabs(aq_fwd);
		
						if (res<BF_PREC) return;
		
						v_aq[q-2]=aq_fwd;
						qi=q-2;
		
					} // END DO Apmin_do
		
					// Check sul ciclo di sostituzione
					assert (qi); // Apmin_error_if1
				       	/*{ 
								WRITE(*,*)
								WRITE(*,*) "Si e' verificato un errore nella subroutine gaunt_xu, caso generale, Apmin=0:"
								WRITE(*,*) "la precisione richiesta per i coefficienti di Gaunt nella backward"
								WRITE(*,*) "e forward recursion non e' stata raggiunta"
								WRITE(*,*)
								error=1
								RETURN
					}*/ // Apmin_error_if1
		
					//Esco dalla subroutine gaunt_xu
					return;
		
				} // Apmin_if
		
				//..........................................................................
				//CASO GENERALE PER LA FORWARD RECURRENCE
				//..........................................................................
		
				//Primo valore per la forward recursion,errore relativo e suo swap
				double aq_fwd=f_aqmax(m,n,mu,nu,qmax);
				double res=fabs(aq_fwd-v_aq[qmax])/fabs(aq_fwd);
				qi=1;
		
				if (res>BF_PREC) { //gen_f_if
				//Se non ho precisione, sostituisco i valori
		
					v_aq[qmax]=aq_fwd;
		
					qi=qmax-1;
		
					//Entro nel ciclo della sostituzione valori
					while(1) { // gen_f_do: DO
		
						switch(qmax-qi) {//gen_q_case:SELECT CASE (qmax-qi)
		
											//$$$$$$$$$$$$$$$$
							case 1: { 			//q=qmax-1
											//$$$$$$$$$$$$$$$$
								//Calcolo Ap4 per qi+2 per vedere quale schema usare
								int p=n+nu-2*(qi+2);
								int Ap4=f_Ap(m,n,mu,nu,p+4);
			
								//Scelgo la ricorrenza a seconda del valore di Ap4
								if (Ap4==0) { // Ap4_q1_if
			
									//Calcolo aq secondo la ricorrenza a 4 termini: uso qi+3 perche' il termine piu' alto e'
									//maggiore di 3 unita' rispetto a qi, pur essendo nullo e non comparendo nella ricorsione
									int p=n+nu-2*(qi+3);
									int p1=p-m-mu;
									int p2=p+m+mu;
									double alphap5=f_alpha(n,nu,p+5);
									double alphap6=f_alpha(n,nu,p+6);
									int Ap2=f_Ap(m,n,mu,nu,p+2);
									int Ap5=f_Ap(m,n,mu,nu,p+5);
									int Ap6=f_Ap(m,n,mu,nu,p+6);
									double c2=(p+2)*(p2+3)*Ap2*(Ap5*Ap6 + (p+4)*(p+6)*(p1+5)*(p2+5)*alphap5);
									double c3=-(p+2)*(p+4)*(p+5)*(p2+3)*(p2+5)*(p2+6)*Ap2*alphap6;
									aq_fwd=-(c2/c3)*v_aq[qi+1];
			
									//A seconda che il mio valore sia 0 o meno confronto i valori
									switch(v_zero[qi]){ //zAp41_case:SELECT CASE (v_zero[qi])
										case 0: 
											v_aq[qi]=0;
											break;
										case 1:
											res=fabs(aq_fwd-v_aq[qi])/fabs(aq_fwd);
											if (res<BF_PREC) {
										       		//EXIT gen_f_do
												//měli bychom breaknout smyčku, ale za 
												//ní už nic „smysluplného“ není
												assert(qi);
												return;
											}
											break;
										default:
											assert(0);
									} // END SELECT zAp41_case
			
									//Qui calcolo il valore successivo dopo aver aggiornato qi:
									//Se v_aq[qi]=0 allora non chiamo cruzan, se no lo chamo e
									//tengo un solo valore
									qi=qi-1;
			
									switch(v_zero[qi]) {//zcz1_case:SELECT CASE (v_zero[qi])
										case 0:
											v_aq[qi]=0;
											qi=qi-1;
											continue; //CYCLE gen_f_do
											break;
										case 1:
											gaunt_cz(m,n,mu,nu,qmax,&(v_aq_cz[qi]),error); // FIXME implementace gaunt_cz
											aq_fwd=(v_aq_cz[qi]);
											res=fabs(aq_fwd-v_aq[qi])/fabs(aq_fwd);
											break;
										default:
											assert(0);
									}
			
										//-----------------
								} else {	//Qui Ap4/=0
										//-----------------
									//Calcolo aq
									int p=n+nu-2*(qi+2);
									int p1=p-m-mu;
									int p2=p+m+mu;
									double alphap2=f_alpha(n,nu,p+2);
									double alphap3=f_alpha(n,nu,p+3);
									double alphap4=f_alpha(n,nu,p+4);
									int Ap2=f_Ap(m,n,mu,nu,p+2);
									int Ap3=f_Ap(m,n,mu,nu,p+3);
									int Ap4=f_Ap(m,n,mu,nu,p+4);
									double c1=Ap2*Ap3*Ap4+(p+1)*(p+3)*(p1+2)*(p2+2)*Ap4*alphap2+ 
									  (p+2)*(p+4)*(p1+3)*(p2+3)*Ap2*alphap3;
									double c2=-(p+2)*(p+3)*(p2+3)*(p2+4)*Ap2*alphap4;
									aq_fwd=-(c1/c2)*v_aq[qi+1]; //E' qui che lo calcolo
			
									//A seconda che il mio valore sia 0 o meno confronto i valori
									switch(v_zero[qi]) { // zAp4d1_case:SELECT CASE (v_zero[qi])
										case 0:
											v_aq[qi]=0;
											qi=qi-1;
											continue; // gen_f_do
											break;
										case 1:
											res=fabs(aq_fwd-v_aq[qi])/fabs(aq_fwd);
											break;
										default:
											assert(0);
									}
			
								} // Ap4_q1_if
								} break;
			
												//$$$$$$$$$$$$$$$$
							case 2: {//CASE(2) gen_q_case	//q=qmax-2
												//$$$$$$$$$$$$$$$$
			
			
			
			
			
			
								//Calcolo Ap4 per qi+2 per vedere quale schema usare
								int p=n+nu-2*(qi+2);
								int Ap4=f_Ap(m,n,mu,nu,p+4);
			
								//Scelgo la ricorrenza a seconda del valore di Ap4
								if (Ap4==0) { // Ap4_q2_if
									//Calcolo aq secondo la ricorrenza a 4 termini: uso qi+3 perche' il termine piu' alto e'
									//maggiore di 3 unita' rispetto a qi, pur essendo nullo e non comparendo nella ricorsione
									int p=n+nu-2*(qi+3);
									int p1=p-m-mu;
									int p2=p+m+mu;
									double alphap2=f_alpha(n,nu,p+2);
									double alphap5=f_alpha(n,nu,p+5);
									double alphap6=f_alpha(n,nu,p+6);
									int Ap2=f_Ap(m,n,mu,nu,p+2);
									int Ap3=f_Ap(m,n,mu,nu,p+3);
									int Ap5=f_Ap(m,n,mu,nu,p+5);
									int Ap6=f_Ap(m,n,mu,nu,p+6);
									double c1=(p+5)*(p1+4)*Ap6*(Ap2*Ap3 + (p+1)*(p+3)*(p1+2)*(p2+2)*alphap2);
									double c2=(p+2)*(p2+3)*Ap2*(Ap5*Ap6 + (p+4)*(p+6)*(p1+5)*(p2+5)*alphap5);
									double c3=-(p+2)*(p+4)*(p+5)*(p2+3)*(p2+5)*(p2+6)*Ap2*alphap6;
									aq_fwd=-(c1/c3)*v_aq[qi+2] -(c2/c3)*v_aq[qi+1];
			
									//A seconda che il mio valore sia 0 o meno confronto i valori
									switch(v_zero[qi]) { // zAp42_case
										case 0:
										v_aq[qi]=0;
										break;
										case 1:
										res=fabs(aq_fwd-v_aq[qi])/fabs(aq_fwd);
										if (res<BF_PREC) { // EXIT gen_f_do
											assert(qi);
											return;
										}
										break;
										default:
										assert(0);
									} // zAp42_case
			
									//Qui calcolo il valore successivo dopo aver aggiornato qi:
									//Se v_aq[qi]=0 allora non chiamo cruzan, se no lo chamo e
									//tengo un solo valore
									qi=qi-1;
			
									switch (v_zero[qi]) {//zcz2_case:SELECT CASE (v_zero[qi])
										case 0:
										v_aq[qi]=0;
										qi=qi-1;
										continue; // gen_f_do
										break;
										case 1:
										//FIXME gaunt_cz !!!!!!!!!!!!!!!!!!
										gaunt_cz(m,n,mu,nu,qmax,&(v_aq_cz[qi]),error);
										aq_fwd=v_aq_cz[qi];
										res=fabs(aq_fwd-v_aq[qi])/fabs(aq_fwd);
										break;
										default:
										assert(0);
									}
								} else { //Qui Ap4!=0
									//Calcolo aq
									int p=n+nu-2*(qi+2);
									int p1=p-m-mu;
									int p2=p+m+mu;
									double alphap2=f_alpha(n,nu,p+2);
									double alphap3=f_alpha(n,nu,p+3);
									double alphap4=f_alpha(n,nu,p+4);
									int Ap2=f_Ap(m,n,mu,nu,p+2);
									int Ap3=f_Ap(m,n,mu,nu,p+3);
									int Ap4=f_Ap(m,n,mu,nu,p+4);
									double c0=(p+2)*(p+3)*(p1+1)*(p1+2)*Ap4*alphap1;
									double c1=Ap2*Ap3*Ap4+(p+1)*(p+3)*(p1+2)*(p2+2)*Ap4*alphap2+ 
									  &(p+2)*(p+4)*(p1+3)*(p2+3)*Ap2*alphap3;
									double c2=-(p+2)*(p+3)*(p2+3)*(p2+4)*Ap2*alphap4;
									aq_fwd=(c0/c2)*v_aq[qi+2]-(c1/c2)*v_aq[qi+1]; //E' qui che lo calcolo
			
									//A seconda che il mio valore sia 0 o meno confronto i valori
									switch(v_zero[qi]) { //zAp4d2_case:SELECT CASE (v_zero[qi])
										case 0:
										v_aq[qi]=0;
										qi=qi-1;
										continue; // gen_f_do
										break;
										case 1:
										res=fabs(aq_fwd-v_aq[qi])/fabs(aq_fwd);
										break;
										default:
										assert(0):
									}
			
								} // Ap4_q2_if
								} break;
			
			
													//$$$$$$$$$$$$$$$$$$$$$$
							default: { //CASE DEFAULT gen_q_case //Per tutti gli altri q
													//$$$$$$$$$$$$$$$$$$$$$$
			
								//Calcolo Ap4 per qi+2 per vedere quale schema usare
								int p=n+nu-2*(qi+2);
								int Ap4=f_Ap(m,n,mu,nu,p+4);
			
								//Scelgo la ricorrenza a seconda del valore di Ap4
								if (Ap4==0) { // Ap4_qq_if
			
									//Calcolo aq secondo la ricorrenza a 4 termini: uso qi+3 perche' il termine piu' alto e'
									//maggiore di 3 unita' rispetto a qi, pur essendo nullo e non comparendo nella ricorsione
									int p=n+nu-2*(qi+3);
									int p1=p-m-mu;
									int p2=p+m+mu;
									double alphap2=f_alpha(n,nu,p+2);
									double alphap5=f_alpha(n,nu,p+5);
									double alphap6=f_alpha(n,nu,p+6);
									int Ap2=f_Ap(m,n,mu,nu,p+2);
									int Ap3=f_Ap(m,n,mu,nu,p+3);
									int Ap5=f_Ap(m,n,mu,nu,p+5);
									int Ap6=f_Ap(m,n,mu,nu,p+6);
									double c0=(p+2)*(p+3)*(p+5)*(p1+1)*(p1+2)*(p1+4)*Ap6*alphap1;
									double c1=(p+5)*(p1+4)*Ap6*(Ap2*Ap3 + (p+1)*(p+3)*(p1+2)*(p2+2)*alphap2);
									double c2=(p+2)*(p2+3)*Ap2*(Ap5*Ap6 + (p+4)*(p+6)*(p1+5)*(p2+5)*alphap5);
									double c3=-(p+2)*(p+4)*(p+5)*(p2+3)*(p2+5)*(p2+6)*Ap2*alphap6;
									aq_fwd=(c0/c3)*v_aq[qi+3]-(c1/c3)*v_aq[qi+2] -(c2/c3)*v_aq[qi+1];
			
									//A seconda che il mio valore sia 0 o meno confronto i valori
									switch(v_zero[qi]) {//zAp4q_case:SELECT CASE (v_zero[qi])
										case 0:
										v_aq[qi]=0;
										break;
										case 1:
										res=fabs(aq_fwd-v_aq[qi])/fabs(aq_fwd);
										if (res<BF_PREC) {// EXIT gen_f_do
											assert(qi);
											return;
										}
										break;
										default:
										assert(0);
									}
			
									//Qui calcolo il valore successivo dopo aver aggiornato qi:
									//Se v_aq[qi]=0 allora non chiamo cruzan, se no lo chiamo e
									//tengo un solo valore.L'if c'e' per non far sballare qi
			
									if (qi>0) { // qi_if
			
										qi=qi-1;
			
										switch(v_zero[qi]) { //zczq_case:SELECT CASE (v_zero[qi])
											case 0:
											v_aq[qi]=0;
											qi=qi-1;
											continue; // CYCLE gen_f_do
											break;
											case 1:
											gaunt_cz(m,n,mu,nu,qmax,&(v_aq_cz[qi]),error); // FIXME je potřeba mít v_aq_cz jako pole?
											aq_fwd=v_aq_cz[qi];
											res=fabs(aq_fwd-v_aq[qi])/fabs(aq_fwd);
											break;
											default:
											assert(0);
										}
			
									} // qi_if
			
										//-----------------
								} else {	//Qui Ap4/=0
										//-----------------
			
									//Calcolo aq
									int p=n+nu-2*(qi+2);
									int p1=p-m-mu;
									int p2=p+m+mu;
									double alphap2=f_alpha(n,nu,p+2);
									double alphap3=f_alpha(n,nu,p+3);
									double alphap4=f_alpha(n,nu,p+4);
									int Ap2=f_Ap(m,n,mu,nu,p+2);
									int Ap3=f_Ap(m,n,mu,nu,p+3);
									int Ap4=f_Ap(m,n,mu,nu,p+4);
									double c0=(p+2)*(p+3)*(p1+1)*(p1+2)*Ap4*alphap1;
									double c1=Ap2*Ap3*Ap4+(p+1)*(p+3)*(p1+2)*(p2+2)*Ap4*alphap2+ 
									  (p+2)*(p+4)*(p1+3)*(p2+3)*Ap2*alphap3;
									double c2=-(p+2)*(p+3)*(p2+3)*(p2+4)*Ap2*alphap4;
									aq_fwd=(c0/c2)*v_aq[qi+2]-(c1/c2)*v_aq[qi+1]; //E' qui che lo calcolo
			
									//A seconda che il mio valore sia 0 o meno confronto i valori
									switch(v_zero[qi]) { //zAp4dq_case:SELECT CASE (v_zero[qi])
										case 0:
										v_aq[qi]=0;
										qi=qi-1;
										continue; // gen_f_do
										break;
										case 1:
										res=fabs(aq_fwd-v_aq[qi])/fabs(aq_fwd);
										break;
										default:
										assert(0);
									}
			
								} // Ap4_qq_if
							} // default
			
						} // END SELECT gen_q_case
		
					//Adesso se la precisione e' raggiunta esco dal ciclo, se no sostituisco e rimango
					if ((res<BF_PREC) || (qi==0) || (fabs(aq_fwd)<fabs(v_aq[qi+1])))  // EXIT
						break; // gen_f_do
		
					//Sono nel ciclo, allora sostituisco eaggiorno indice e residuo
					v_aq[qi]=aq_fwd;
					qi=qi-1;
		
					} // END DO gen_f_do
		
					// Check sul ciclo di sostituzione
					assert(qi); /* gen_error_if1: if (qi==0) {
								WRITE(*,*)
								WRITE(*,*) "Si e' verificato un errore nella subroutine gaunt_xu,caso generale:"
								WRITE(*,*) "la precisione richiesta per i coefficienti di Gaunt nella backward"
								WRITE(*,*) "e forward recursion non e' stata raggiunta"
								WRITE(*,*)
								error=1
								RETURN
					} */ // gen_error_if1
		
		
				} // gen_f_if

		} // big_if
	} // qmax_case
} // gaunt_xu

		
		

