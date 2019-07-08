VSWF conventions
================

In general, the (transversal) VSWFs can be defined using (some) vector spherical harmonics
as follows: \f[
	\wfm\pr{k\vect r}_{lm} = \sphbes_l(kr) \vshrot_{lm} (\uvec r),\\
	\wfe\pr{k\vect r}_{lm} = \frac{\frac{\ud}{\ud(kr)}\pr{kr\sphbes_l(kr)}}{kr} \vshgrad_{lm}(\uvec r)
				+ \sqrt{l(l+1)} \frac{\sphbes_l(kr)}{kr} \vshrad_{lm}(\uvec r),
\f]
where at this point, we don't have much expectations regarding the
normalisations and phases of the
"rotational", "gradiental" and "radial" vector spherical harmonics
\f$  \vshrot, \vshgrad, \vshrad \f$, and the waves can be of whatever "direction"
(regular, outgoing, etc.) depending on the kind of the spherical Bessel function
\f$ \sphbes \f$. 
We only require that the spherical harmonic degree \f$ l \f$
is what it is supposed to be. The meaning of the order $m$ may vary depending
on convention. Moreover, in order to \f$ \wfe \f$ be a valid "electric" multipole wave,
there is a fixed relation between radial and gradiental vector spherical harmonics
(more on that later).

Let us define the "dual" vector spherical harmonics \f$ \vshD_{\tau lm} \f$ as follows:
\f[
	\int_\Omega \vsh_{\tau lm} (\uvec r) \cdot \vshD_{\tau' l'm} (\uvec r) \, \ud \Omega 
		= \delta_{\tau', \tau}\delta_{l',l} \delta_{m',m}
\f]
where the \f$ \cdot \f$ symbol here means the bilinear form of the vector components
without complex conjugation (which is included in the "duality" mapping).

For the sake of non-ambiguity, let us define the "canonical" associated Legendre polynomials
as in \cite DLMF TODO exact refs:
\f[
	\rawLeg{l}{0}(x) = \frac{1}{2^n n!} \frac{\ud^n}{\ud x^n} \pr{x^2-1}^n , \\
	\rawLeg{l}{m}(x) = \pr{1-x^2}^{m/2} \frac{\ud^m}{\ud x^m} \rawLeg{l}{0},\quad\abs{x}\le 1, m \ge 0, \\
	\rawLeg{l}{m}(x) = (-1)^\abs{m} \frac{(l-\abs{m})!}{(l+\abs{m})!} \rawLeg{l}{\abs{m}}, 
		\quad \abs{x} \le 1, m < 0.
\f]


Literature convention table
---------------------------

| Source	| VSWF definition  	| E/M interrelations | VSWF norm  	| CS Phase  	|  Field expansion 	|  Radiated power | Notes |
|---	|---	|---	|---	|---	|---	|--- 	|--- |
| Kristensson I \cite kristensson_spherical_2014 	|  \f[ \wfkcreg, \wfkcout= \dots \f] 	| \f[
	\wfkcreg_{1lm} = \frac{1}{k}\nabla\times\wfkcreg_{2lm}, \\
	\wfkcreg_{2lm} = \frac{1}{k}\nabla\times\wfkcreg_{1lm},
\f] and analogously for outgoing waves \f$ \wfkcout \f$, eq. (2.8) onwards. 	|  	| Yes, in the spherical harmonics definition, cf. sect. D.2.  	| \f[ 
	\vect E = k \sqrt{\eta_0\eta} \sum_n \left( \wckcreg_n  \wfkcreg_n + \wckcout_n \wfkcout_n  \right), 
	\\ 
	\vect H =  \frac{k \sqrt{\eta_0\eta}}{i\eta_0\eta} \sum_n \left( \wckcreg_n  \wfkcreg_n + \wckcout_n \wfkcout_n  \right),
\f] but for plane wave expansion \cite kristensson_spherical_2014 sect. 2.5 K. uses a different definition (same as in Kristensson II).  	| \f[
	 P = \frac{1}{2} \sum_n \left( \abs{\wckcout_n}^2 +\Re \left(\wckcout_n\wckcreg_n^{*}\right)\right)
 \f]	| The \f$ \wckcreg, \wckcout \f$	coefficients have dimension \f$ \sqrt{\mathrm{W}} \f$. |
| Kristensson II \cite kristensson_scattering_2016	| \f[ \wfkrreg, \wfkrout= \dots \f] 	|  \f[
	\nabla\times\wfkrreg_{\tau n} = k\wfkrreg_{\overline{\tau} n},
\f] eq. (7.7) and analogously for outgoing waves \f$ \wfkrout \f$. 	| 	|   	| \f[ 
	\vect E = \sum_n \left( \wckrreg_n  \wfkrreg_n + \wckrout_n \wfkrout_n  \right), 
	\\ 
	\vect H =  \frac{1}{i\eta_0\eta} \sum_n \left( \wckrreg_n  \wfkrreg_n + \wckrout_n \wfkrout_n  \right)
\f] 	| \f[
	 P = \frac{1}{2k^2\eta_0\eta} \sum_n \left( \abs{\wckrout_n}^2 +\Re \left(\wckrout_n\wckrreg_n^{*}\right)\right)
 \f]	| The \f$ \wckrreg, \wckrout \f$ coefficients have dimension \f$ \mathrm{V/m} \f$. |
| Reid \cite reid_electromagnetism_2016	|   | \f[
	\nabla\times\wfr_{lmM} = -ik\wfr_{lmN}, \\ \nabla\times\wfr_{lmN} = +ik\wfr_{lmM}. 
\f] 	|	|  |  \f[
	\vect E = \sum_\alpha \pr{ \wcrreg_\alpha \wfrreg_\alpha + \wcrout_\alpha \wfrout_\alpha }, \\
	\vect H = \frac{1}{Z_0Z^r} \sum_\alpha \pr{ \wcrreg_\alpha \sigma_\alpha\wfrreg_\overline{\alpha} +
		 \wcrout_\alpha \sigma_\alpha\wfrout_\overline{\alpha}},
\f] where \f$ \sigma_{lmM} = +1, \sigma_{lmN}=-1, \overline{lmM}=lmM, \overline{lmN}=lmM, \f$  cf. eq. (6). The notation is not extremely consistent throughout Reid's memo.	| 	| 	|
| Taylor \cite taylor_optical_2011	| \f[
	\wfet_{mn}^{(j)}	=	\frac{n(n+1)}{kr}\sqrt{\frac{2n+1}{4\pi}\frac{\left(n-m\right)!}{\left(n+m\right)!}}P_{n}^{m}\left(\cos\theta\right)e^{im\phi}z_{n}^{j}\left(kr\right)\uvec{r} \\
		+\left[\tilde{\tau}_{mn}\left(\cos\theta\right)\uvec{\theta}+i\tilde{\pi}_{mn}\left(\cos\theta\right)\uvec{\phi}\right]e^{im\phi}\frac{1}{kr}\frac{\ud\left(kr\,z_{n}^{j}\left(kr\right)\right)}{\ud(kr)}, \\ 
	\wfmt_{mn}^{(j)}	=	\left[i\tilde{\pi}_{mn}\left(\cos\theta\right)\uvec{\theta}-\tilde{\tau}_{mn}\left(\cos\theta\right)\uvec{\phi}\right]e^{im\phi}z_{n}^{j}\left(kr\right)
\f]  	|	|	\f[
	\int_{S(kr)} \wfmt_{mn}^{(j)} \wfmt_{m'n'}^{(j)}\,\ud S = n(n+1) \abs{z_n^{(j)}}^2 \delta_{m,m'}\delta_{n,n'} ,\\
	\int_{S(kr)} \wfet_{mn}^{(j)} \wfet_{m'n'}^{(j)}\,\ud S =
           \pr{\pr{n(n+1)}^2 \abs{\frac{z_n^{(j)}}{kr}}^2 + n(n+1)\abs{\frac{1}{kr}\frac{\ud}{\ud(kr)}\pr{kr z_n^{(j)}}} } \delta_{m,m'}\delta_{n,n'} ,
\f] cf. \cite taylor_optical_2011, eqs. (2.40–41). I suspect that this is also wrong and \f$\delta_{m,m'}\f$ should be replaced with \f$\delta_{m,-m'}\f$. |	| \f[ 
	\vect E = \sum_{mn} \pr{-i \pr{\wcetreg_{mn}\wfetreg_{mn} + \wcmtreg_{mn}\wfmtreg{mn}} +i \pr{\wcetout_{mn}\wfetout_{mn} + \wcmtout_{mn}\wfmtout_{mn}}}, \\
	\vect H = n_{ext}\sum_{mn} \pr{- \pr{\wcmtreg_{mn}\wfetreg_{mn} + \wcetreg_{mn}\wfmtreg{mn}} + \pr{\wcmtout_{mn}\wfetout_{mn} + \wcetout_{mn}\wfmtout_{mn}}},
\f] 	| 	| Different sign for regular/scattered waves! Also WTF are the units of \f$ n_{ext} \f$?  The whole definition seems rather inconsistent. |

