VSWF conventions
================


| Source	| VSWF definition  	| VSWF norm  	| CS Phase  	|  Field expansion 	|  Radiated power | Notes |
|---	|---	|---	|---	|---	|---	|--- 	|
| Kristensson I \cite kristensson_spherical_2014 	|  \f[ \wfkcreg, \wfkcout= \dots \f] 	|   	|   	| \f[ 
	\vect E = k \sqrt{\eta_0\eta} \sum_n \left( \wckcreg_n  \wfkcreg_n + \wckcout_n \wfkcout_n  \right), 
	\\ 
	\vect H =  \frac{k \sqrt{\eta_0\eta}}{i\eta_0\eta} \sum_n \left( \wckcreg_n  \wfkcreg_n + \wckcout_n \wfkcout_n  \right)
\f] 	| \f[
	 P = \frac{1}{2} \sum_n \left( \abs{\wckcout_n}^2 +\Re \left(\wckcout_n\wckcreg_n^{*}\right)\right)
 \f]	| The \f$ \wckcreg, \wckcout \f$	coefficients have dimension \f$ \sqrt{\mathrm{W}} \f$. |
| Kristensson II \cite kristensson_scattering_2016	| \f[ \wfkrreg, \wfkrout= \dots \f] 	|   	|   	| \f[ 
	\vect E = \sum_n \left( \wckrreg_n  \wfkrreg_n + \wckrout_n \wfkrout_n  \right), 
	\\ 
	\vect H =  \frac{1}{i\eta_0\eta} \sum_n \left( \wckrreg_n  \wfkrreg_n + \wckrout_n \wfkrout_n  \right)
\f] 	| \f[
	 P = \frac{1}{2k^2\eta_0\eta} \sum_n \left( \abs{\wckrout_n}^2 +\Re \left(\wckrout_n\wckrreg_n^{*}\right)\right)
 \f]	| The \f$ \wckrreg, \wckrout \f$ coefficients have dimension \f$ \mathrm{V/m} \f$. |
| Reid \cite reid_electromagnetism_2016	|   |  	|   	|   	| 	| 	|
| Taylor \cite taylor_optical_2011	| \f[
	\wfet_{mn}^{(j)}	=	\frac{n(n+1)}{kr}\sqrt{\frac{2n+1}{4\pi}\frac{\left(n-m\right)!}{\left(n+m\right)!}}P_{n}^{m}\left(\cos\theta\right)e^{im\phi}z_{n}^{j}\left(kr\right)\uvec{r} \\
		+\left[\tilde{\tau}_{mn}\left(\cos\theta\right)\uvec{\theta}+i\tilde{\pi}_{mn}\left(\cos\theta\right)\uvec{\phi}\right]e^{im\phi}\frac{1}{kr}\frac{\ud\left(kr\,z_{n}^{j}\left(kr\right)\right)}{\ud(kr)}, \\ 
	\wfmt_{mn}^{(j)}	=	\left[i\tilde{\pi}_{mn}\left(\cos\theta\right)\uvec{\theta}-\tilde{\tau}_{mn}\left(\cos\theta\right)\uvec{\phi}\right]e^{im\phi}z_{n}^{j}\left(kr\right)
\f]  | 	|   	| \f[ 
	\vect E = \sum_{mn} \pr{-i \pr{\wcetreg_{mn}\wfetreg_{mn} + \wcmtreg_{mn}\wfmtreg{mn}} +i \pr{\wcetout_{mn}\wfetout_{mn} + \wcmtout_{mn}\wfmtout_{mn}}}, \\
	\vect H = n_{ext}\sum_{mn} \pr{- \pr{\wcmtreg_{mn}\wfetreg_{mn} + \wcetreg_{mn}\wfmtreg{mn}} + \pr{\wcmtout_{mn}\wfetout_{mn} + \wcetout_{mn}\wfmtout_{mn}}},
\f]  	| 	| Different sign for regular/scattered waves! Also WTF are the units of \f$ n_{ext} \f$? |

