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
int Ap(int m, int n, int mu, int nu, int p) {
	return p*(p-1)*(m-mu)-(m+mu)*(n-nu)*(n+nu+1);
}

// coefficient a(m,n,mu,nu,1) normalized by the backward recursion (from vec_trans.f90)

double f_a1norm(int m, int n, int mu, int nu) {
	int n4 = n + nu - m - mu;
	return ((2.*n + 2.*nu - 3.) / 2.) * (1. - ((2.*n + 2.*nu - 1.) / (n4 * (n4-1.))) 
			* ((m - n) * (m - n + 1.) / (2.*n - 1.) + (mu-nu) * (mu-nu+1.)/(2.nu-1.)));
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
	                -lnf(n+nu+pmin+1)-lnf(n-qmax)-lnf(nu-qmax)-lnf(n-m)-lnf(nu-mu)
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

// gaunt_xu from vec_trans.f90
// btw, what is the difference from gaunt_xu2?
// FIXME set some sensible return value
double gaunt_xu2(int m, int n, int mu, int nu, int qmax, double *v_aq, int *error) {
	
	int v_zero[qmax] = {0};
	*error = 0;

	if(abs(m)>n || abs(mu)=nu) {
		*error = 1;
		fprintf(stderr, "invalid values for m, n, mu or nu\n")
		return NAN;
	}

	switch(qmax) {
		case 0:
			v_aq[0]  = f_a0(m,n,mu,nu);
			break;
		case 1:
			v_aq[0] = f_a0(m,n,mu,nu);
			v_aq[1] = v_aq[0] + f_a1norm(m,n,mu,nu);
			break;
		case 2:
			v_aq[0] = f_a0(m,n,mu,nu);
			v_aq[1] = v_aq[0] + f_a1norm(m,n,mu,nu);
			v_aq[2] = v_aq[0] * f_a2normr(m,n,mu,nu,v_aq[1]/v_aq[0]);
			break;
		default:
			if (m == 0 && mu == 0) {
				v_aq[0] = f_a0(m,n,mu,nu);

				// backward recursion
				for (int q = 1; q <= qmax; ++q) {
					int p = n + nu - 2*q;
					double c0 = f_alpha(n,nu,p+1);
					double c1 = f_alpha(n,nu,p+2);
					v_aq[q] = (c1/c0) * v_aq[q-1];
				}
			} else if (mu == m && nu == n) {
				v_aq[0] = f_a0(m,n,mu,nu);

				// backward recursion
				for (int q = 1; q <= qmax; ++q) {
					// calculate pre-coefficients
					int p = n + nu - 2*q;
					int p1 = p - m - mu;
					int p2 = p + m + mu;

					// calculate recursion coefficients
					double c0 = (p+2) * (p1+1) * f_alpha(n,nu,p+1);
					double c1 = (p+1) * (p2+2) * f_alpha(n,nu,p+2);

					//recursion
					v_aq[q] = (c1/c0) * v_aq[q-1];
				}
			} else if (mu == -m) {
				v_aq[0] = f_a0(m,n,mu,nu);
				v_aq[1] = f_a1norm(m,n,mu,nu)*v_aq[0];

				// backward recursion
				for (int q = 2; q <= qmaq; ++q) {
					// calculate pre-coefficient
					int p = n + nu - 2*q;

					// calculate recursion coefficients
					double c0 = f_alpha(n, nu, p+1);
					double c1 = 4*isq(m) + f_alpha(n,nu,p+2) + f_alpha(n,nu,p+3);
					double c2 = - f_alpha(n, nu, p+4);

					// recursion
					v_aq[q] = (c1/c0)*v_aq[q-1] + (c2/c0)*v_aq[q-2];
				}
				
				// forward recursion	
		                // Primo valore per la forward recursion,errore relativo e suo swap
		                double aq_fwd=f_aqmax(m,n,mu,nu,qmax)
		                double res=fabs(aq_fwd-v_aq[qmax])/fabs(aq_fwd);
		
		                //Se non ho precisione, sostituisco i valori
		                if (res>BF_PREC) {
		                        v_aq[qmax]=aq_fwd;
		                        int qi=1;
					int zeroswitch = 0; // black magic (gevero's "switch")
		                        //Entro nel ciclo della sostituzione valori
		                        for( int q=qmax-1,q>=0,--q) { // tre_f_do
						switch(qmax-q) {// FIXME switch -> if/else
							case 1: {// q==qmax-1
		                                        	//Calcolo v_aq[qmax-1]
			                                        int p=n+nu-2*(q+2);
			                                        double c1=4*isq(m)+f_alpha(n,nu,p+2)+f_alpha(n,nu,p+3);
		                                        	double c2=-f_alpha(n,nu,p+4)
			                                        double aq_fwd=-(c1/c2)*v_aq[qmax];
		
		                                        	//If a secondo che v_aq[qmax-1] sia zero o no
		                                        	if  (aq_fwd/v_aq[max]<ZERO_THRESHOLD) {
		                                          	      v_aq[q]=aq_fwd;  //Zero
		                                              	      zeroswitch=1; 
		                                              	      continue; // tre_f_do
								} else
		                                                	res=fabs(aq_fwd-v_aq[q])/fabs[aq_fwd];//Diverso da zero
								}
								break; // case 1
							default: { //Per tutti gli altri q
		                                        //Calcolo v_aq[qmax-1]
		                                        int p=n+nu-2*(q+2);
		                                        double c0=f_alpha(n,nu,p+1);
		                                        double c1=4*(m**2)+f_alpha(n,nu,p+2)+f_alpha(n,nu,p+3);
		                                        double c2=-f_alpha(n,nu,p+4);
		                                        aq_fwd=-(c1/c2)*v_aq[q+1]+(c0/c2)*v_aq[q+2];
							switch(zeroswitch}{ //FIXME switch -> if/else
		                                        	case 1: //Il valore precedente e' zero
			                                                res=fabs(aq_fwd-v_aq[q])/fabs(aq_fwd);
									break;
								case 0: //Il valore precedente e' diverso da zero
		                                                	//If a secondo che v_aq[q] sia zero o no
		                                                	if ((aq_fwd/v_aq[q+1])<ZERO_THRESHOLD) {
		                                                        	v_aq[q]=aq_fwd;  //Zero
		                                                        	zeroswitch=1;
										continue; // tre_f_do
									} else //Diverso da zero
										res=fabs(aq_fwd-v_aq[q])/fabs(aq_fwd);
									break;
								default:
									assert(0);
							}
						}
		                        	//Adesso se la precisione e' raggiunta esco dal ciclo, 
						//se no sostituisco e rimango
		                        	if (res<BF_PREC || q==0) 
							break; //tre_f_do
		                        	//Sono nel ciclo, allora sostituisco eaggiorno indice e residuo
			                        v_aq[q]=aq_fwd;
		                        	qi=q;
					}
		                // Check sul ciclo di sostituzione
				assert(q);
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
				}
			} else { // caso generale (4)
				// backward
				// Calcolo Ap4 per q=0, se e' uguale a zero chiamo cruzan ed esco dal ciclo
				int Ap4=f_Ap(m,n,mu,nu,n+nu-qmax+4);
		
				if (Ap4==0) { // cz_if1
					gaunt_cz(m,n,mu,nu,qmax,v_aq,error);
					return;
				}
		
				// Calcolo direttamente i primi tre valori della ricorsione
				v_aq[0]=f_a0(m,n,mu,nu);
				v_aq[1]=v_aq[0]*f_a1norm(m,n,mu,nu);
		
				for (int q=2; q <= qmax; ++q) {//gen_b_do: DO q=2,qmax 
							//Calcolo pre-coeff. : questi li calcoli qui per poter uscire
							int p=n+nu-2*q;
							int Ap2=f_Ap(m,n,mu,nu,p+2);
							if (Ap2==0) break; //Esco dal loop perche non posso fare la fwd recursion
		
							//Calcolo pre-coefficienti
							int p1=p-m-mu;
							int p2=p+m+mu;
							double alphap1=f_alpha(n,nu,p+1);
							double alphap2=f_alpha(n,nu,p+2);
							double alphap3=f_alpha(n,nu,p+3);
							double alphap4=f_alpha(n,nu,p+4);
							int Ap3=f_Ap(m,n,mu,nu,p+3);
							int Ap4=f_Ap(m,n,mu,nu,p+4);
		
							//Calcolo coefficienti ricorsione
							double c0=(p+2)*(p+3)*(p1+1)*(p1+2)*Ap4*alphap1;
							double c1=Ap2*Ap3*Ap4+(p+1)*(p+3)*(p1+2)*(p2+2)*Ap4*alphap2+ \
							  (p+2)*(p+4)*(p1+3)*(p2+3)*Ap2*alphap3;
							double c2=-(p+2)*(p+3)*(p2+3)*(p2+4)*Ap2*alphap4;
		
							//Ricorsione
							v_aq[q]=(c1/c0)*v_aq[q-1]+(c2/c0)*v_aq[q-2];
				} // END DO gen_b_do

				if (Ap2==0) { //cz_if2
					gaunt_cz(m,n,mu,nu,qmax,v_aq,error);
					return;
				}
		
				//!-----------------------------------
				//!FORWARD
				//!-----------------------------------
		
				//Calcolo pmin,Apmin e la mia variabile logica
				int pmin=n+nu-2*qmax;
				int Apmin=f_Ap(m,n,mu,nu,pmin);
				int test=(Apmin==0) && ((pmin==m+mu+1) || (pmin==-m-mu+1));
		
				//!........................................................
				//!Se la mia variabile logica e' vera, Faccio il mio conto
				//!........................................................
				if (test) { //Apmin_if
					//Il valore per qmax allora e' zero
					v_aq[qmax]=0;
		
					//Calcolo il secondo valore, e se la precisione e' raggiunta esco
					double aq_fwd=f_aqmax_1(m,n,mu,nu,qmax);
					double res=fabs(aq_fwd-v_aq(qmax-1))/fabs(aq_fwd);
					if (res<BF_PREC) 
						return;
		
					// Assegno il secondo valore e faccio il ciclo
					v_aq[qmax-1]=aq_fwd;
					int qi=1;

					for (int i=qmax-2; i>=0; --i) { //Apmin_do
						//Calcolo pre-coefficienti
						int p=n+nu-2*REALq+2;
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
						double c1=Ap2*Ap3*Ap4+(p+1)*(p+3)*(p1+2)*(p2+2)*Ap4*alphap2+ \
						  (p+2)*(p+4)*(p1+3)*(p2+3)*Ap2*alphap3;
						double c2=-(p+2)*(p+3)*(p2+3)*(p2+4)*Ap2*alphap4;
		
						//Ricorsione e residuo
						double aq_fwd=-(c1/c2)*v_aq(q+1)+(c0/c2)*v_aq(q+2);
						res=fabs(aq_fwd-v_aq[q])/fabs(aq_fwd);
						if (res<BF_PREC) 
							return;
		
						v_aq[q]=aq_fwd;
						qi=q;
					} // Apmin_do
		
					// Check sul ciclo di sostituzione
					if (qi == 0) //Apmin_error_if1
						assert(0); 
					/*
								WRITE(*,*) "Si e' verificato un errore nella subroutine gaunt_xu, caso generale, Apmin=0:"
								WRITE(*,*) "la precisione richiesta per i coefficienti di Gaunt nella backward"
								WRITE(*,*) "e forward recursion non e' stata raggiunta"
								WRITE(*,*)
								error=1
								RETURN
					END IF Apmin_error_if1
					*/
		
					//Esco dalla subroutine gaunt_xu
					return
		
				} // Apmin_if
		
				//!..............................................
				//!Qui lavoro se la mia variabile logica e' falsa
				//!Tutto e' identico al caso (3)
				//!..............................................
		
				//!Primo valore per la forward recursion,errore relativo e suo swap
				double aq_fwd=f_aqmax(m,n,mu,nu,qmax);
				double res = fabs(aq_fwd-v_aq[qmax])/fabs(aq_fwd);
				int qi=1;
		
				//Se non ho precisione, sostituisco i valori
				if (res > BF_PREC) { // gen_f_if
					v_aq[qmax]=aq_fwd;

					int zeroswitch = 0; // black magic
		
					//Entro nel ciclo della sostituzione valori
					for (int q = qmax-1; q >= 0; --q) { //gen_f_do
						switch(qmax-q) { //gen_q_case
						case 1: //q=qmax-1
							//Calcolo aq
							int p=n+nu-2*(q+2);
							int p1=p-m-mu;
							int p2=p+m+mu;
							double alphap1=f_alpha(n,nu,p+1);
							double alphap2=f_alpha(n,nu,p+2);
							double alphap3=f_alpha(n,nu,p+3);
							double alphap4=f_alpha(n,nu,p+4);
							int Ap2=f_Ap(m,n,mu,nu,p+2);
							int Ap3=f_Ap(m,n,mu,nu,p+3);
							int Ap4=f_Ap(m,n,mu,nu,p+4);
							double c1=Ap2*Ap3*Ap4+(p+1)*(p+3)*(p1+2)*(p2+2)*Ap4*alphap2+ \
							  (p+2)*(p+4)*(p1+3)*(p2+3)*Ap2*alphap3;
							double c2=-(p+2)*(p+3)*(p2+3)*(p2+4)*Ap2*alphap4;
							double aq_fwd=-(c1/c2)*v_aq[qmax] //E' qui che lo calcolo
		
							//If a secondo che v_aq(qmax-1) sia zero o no
							if (aq_fwd/v_aq[qmax] < ZERO_THRESHOLD) { // gen_zero1_if
									v_aq[q]=aq_fwd;	//Zero
									zeroswitch=1;
									continue; // gen_f_do
							} else 
									res=fabs(aq_fwd-v_aq[q])/fabs(aq_fwd);	//Diverso da zero
							break;
						default: // Per tutti gli altri q
							//Calcolo aq
							int p=n+nu-2*(q+2);
							int p1=p-m-mu;
							int p2=p+m+mu;
							double alphap1=f_alpha(n,nu,p+1);
							double alphap2=f_alpha(n,nu,p+2);
							double alphap3=f_alpha(n,nu,p+3);
							double alphap4=f_alpha(n,nu,p+4);
							int Ap2=f_Ap(m,n,mu,nu,p+2);
							int Ap3=f_Ap(m,n,mu,nu,p+3);
							int Ap4=f_Ap(m,n,mu,nu,p+4);
							double c0=(p+2)*(p+3)*(p1+1)*(p1+2)*Ap4*alphap1;
							double c1=Ap2*Ap3*Ap4+(p+1)*(p+3)*(p1+2)*(p2+2)*Ap4*alphap2+ \
							  (p+2)*(p+4)*(p1+3)*(p2+3)*Ap2*alphap3;
							double c2=-(p+2)*(p+3)*(p2+3)*(p2+4)*Ap2*alphap4;
							double aq_fwd=-(c1/c2)*v_aq(q+1)+(c0/c2)*v_aq(q+2);
		
							//Case a seconda che il valore appena precedente sia zero o meno
							switch (zeroswitch) { //gen_switch_case
								case 1: //Il valore precedente e' zero		
									res=fabs(aq_fwd-v_aq[q])/fabs(aq_fwd);
									break;
								case 0: //Il valore precedente e' diverso da zero
									//If a secondo che v_aq[q] sia zero o no
									if (aq_fwd/v_aq[q+1] < ZERO_THRESHOLD/10.) {// gen_zero2_if

										v_aq[q]=aq_fwd;	//Zero
										switch=1
										continue; // gen_f_do
									} else {
										res=fabs(aq_fwd-v_aq[q])/fabs(aq_fwd); // Diverso da zero
									} // gen_zero2_if
									break;
								default:
									assert(0);
							} // gen_switch_case
						} // END SELECT gen_q_case
		
					//Adesso se la precisione e' raggiunta esco dal ciclo, se no sostituisco e rimango
					if ((res<BF_PREC) || (q == 0))
						break;
		
					//Sono nel ciclo, allora sostituisco eaggiorno indice e residuo
					v_aq[q]=aq_fwd;
					qi=q;
		
					} // gen_f_do
		
					// Check sul ciclo di sostituzione
					assert(qi);
					/*
					gen_error_if1: IF (qi==0) THEN
								WRITE(*,*)
								WRITE(*,*) "Si e' verificato un errore nella subroutine gaunt_xu,caso generale:"
								WRITE(*,*) "la precisione richiesta per i coefficienti di Gaunt nella backward"
								WRITE(*,*) "e forward recursion non e' stata raggiunta"
								WRITE(*,*)
								error=1
								RETURN
					END IF gen_error_if1
					*/
		
		
				}//END IF gen_f_if
		
			} // END IF big_if
	} // END SELECT qmax_case
}	
		
		



		



/*
// gaunt_xu2 from vec_trans.f90
// btw, what is the difference from gaunt_xu?
double gaunt_xu2(int m, int n, int mu, int nu, int qmax, double *v_aq, int *error) {
	
	int v_zero[qmax] = {0};
	*error = 0;

	if(abs(m)>n || abs(mu)=nu) {
		*error = 1;
		fprintf(stderr, "invalid values for m, n, mu or nu\n")
		return NAN;
	}

	switch(qmax) {
		case 0:
			v_aq[0]  = f_a0(m,n,mu,nu);
			break;
		case 1:
			v_aq[0] = f_a0(m,n,mu,nu);
			v_aq[1] = v_aq[0] + f_a1norm(m,n,mu,nu);

			// Check for zeros
			if (fabs(v_aq[1]/v_aq[0]) < ZERO_THRESHOLD) {
				v_aq[1] = 0.;
				v_zero[1] = 0.;
			}
			break;
		case 2:
			v_aq[0] = f_a0(m,n,mu,nu);
			v_aq[1] = v_aq[0] + f_a1norm(m,n,mu,nu);
	
			// Check for zeros
			if (fabs(v_aq[1]/v_aq[0]) < ZERO_THRESHOLD) {
				v_aq[1] = 0.;
				v_zero[1] = 0.;
			}

			v_aq[2] = v_aq[0] * f_a2norm_
			break;
		

			!!!!!!!!!TODO CONTINUE HERE
*/


