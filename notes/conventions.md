VSWF conventions {#vswf_conventions}
====================================

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

The problem with conventions starts with the very definition of associated Legendre / Ferrers functions.

For the sake of non-ambiguity, let us first define the "canonical" associated Legendre/Ferrers polynomials
*without* the Condon-Shortley phase.
\f[
	\rawLeg{l}{0}(x) = \frac{1}{2^n n!} \frac{\ud^n}{\ud x^n} \pr{x^2-1}^n , \\
	\rawLeg{l}{m}(x) = \pr{1-x^2}^{m/2} \frac{\ud^m}{\ud x^m} \rawLeg{l}{0},\quad\abs{x}\le 1, m \ge 0, \\
	\rawLeg{l}{m}(x) = (-1)^\abs{m} \frac{(l-\abs{m})!}{(l+\abs{m})!} \rawLeg{l}{\abs{m}}, 
		\quad \abs{x} \le 1, m < 0.
\f]
DLMF \cite NIST:DLMF has for non-negative integer \f$m\f$ (18.5.5), (14.6.1), (14.9.3):
\f[
	\dlmfFer{\nu}{} = \dlmfLeg{\nu}{} = \frac{1}{2^n n!} \frac{\ud^n}{\ud x^n} \pr{x^2-1}^n , \\
	\dlmfFer{\nu}{m}\left(x\right)=(-1)^{m}\left(1-x^2\right)^{m/2}\frac{{
	\ud}^{m}\dlmfFer{\nu}{}\left(x\right)}{{\ud x}^{m}},\\
	%\dlmfLeg{\nu}{m}\left(x\right)=\left(-1+x^2\right)^{m/2}\frac{{
	%\ud}^{m}\dlmfLeg{\nu}{}\left(x\right)}{{\ud x}^{m}},\\
\f]
where the connection to negative orders is
\f[
	\dlmfFer{\nu}{m}(x) = (-1)^m \frac{\Gamma\pr{\nu-m+1}}{\Gamma\pr{\nu+m+1}}\dlmfFer{\nu}{m}(x),\\
	%\dlmfLeg{\nu}{m}(x) =        \frac{\Gamma\pr{\nu-m+1}}{\Gamma\pr{\nu+m+1}}\dlmfLeg{\nu}{m}(x).\\
\f]
Note that there are called "Ferrers" functions in DLMF, while the "Legendre" functions have slightly
different meaning / conventions (Ferrers functions being defined for \f$ \abs{x} \le 1 \f$, whereas
Legendre for \f$ \abs{x} \ge 1 \f$. We will not use the DLMF "Legendre" functions here.

One sees that \f$ \dlmfFer{l}{m} = (-1)^m \rawFer{l}{m} \f$, i.e. the Condon-Shortley phase is
already included in the DLMF definitions of Ferrers functions.

GSL computes \f$ \rawFer{l}{m} \f$ unless the corresponding `csphase` argument is set to 
\f$-1\f$ (then it computes \f$ \dlmfFer{l}{m} \f$). This is not explicitly obvious from the docs 
\cite GSL,
but can be tested by running `gsl_sf_legendre_array_e` for some specific arguments and comparing signs.


Convention effect on translation operators
------------------------------------------

Let us declare VSWFs in Kristensson's conventions below, 
\f$ \wfkc \f$ \cite kristensson_spherical_2014, 
\f$ \wfkr \f$ \cite kristensson_scattering_2016, as the "canonical"
spherical waves based on complex and real spherical harmonics, respectively.
They both have the property that the translation operators \f$ \tropRrr{}{},\tropSrr{}{} \f$ 
that transform
the VSWF field expansion coefficients between different origins, e.g.
\f[
	\wfkcreg(\vect{r}) = \tropRrr{\vect r}{\vect r'} \wfkcreg(\vect{r'}),
\f]
actually consist of two different submatrices $A,B$ for the same-type and different-type
(in the sense of "electric" versus "magnetic" waves) that repeat themselves once:
\f[
	\begin{bmatrix} \wfkcreg_1(\vect{r}) \\ \wfkcreg_2(\vect{r}) \end{bmatrix} 
	= \begin{bmatrix} A & B \\ B & A \end{bmatrix}(\vect{r} \leftarrow \vect{r'})
	\begin{bmatrix} \wfkcreg_1(\vect{r'}) \\ \wfkcreg_2(\vect{r'}) \end{bmatrix}.
\f]
(This symmetry holds also for singular translation operators \f$ \tropSrr{}{} \f$
and real spherical harmonics based VSWFs \f$ \wfkr \f$.)

However, the symmetry above will not hold like this in some stupider convention.
Let's suppose that one uses a different convention with some additional coefficients
compared to the canonical one,
\f[ 
	\wfm_{lm} = \alpha_{\wfm lm} \wfkc_{1lm},\\
	\wfe_{lm} = \alpha_{\wfe lm} \wfkc_{2lm}.\\
\f]
and with field expansion (WLOG assume regular fields only)
\f[ \vect E = c_{\wfe l m} \wfe_{lm} + c_{\wfm l m } \wfm_{lm}. \f]
Under translations, the coefficients then transform like
\f[
	\begin{bmatrix} \alpha_\wfe(\vect{r}) \\ \alpha_\wfm(\vect{r}) \end{bmatrix} 
	= \begin{bmatrix} R_{\wfe\wfe} & R_{\wfe\wfm} \\ 
	 	          R_{\wfm\wfe} & R_{\wfm\wfm} 
	  \end{bmatrix}(\vect{r} \leftarrow \vect{r'})
	\begin{bmatrix} \alpha_\wfe(\vect{r'}) \\ \alpha_\wfm(\vect{r'}) \end{bmatrix},
\f]
and by substituting and comparing the expressions for canonical waves above, one gets
\f[
	R_{\wfe,lm;\wfe,l'm'} = \alpha_{\wfe lm}^{-1} A \alpha_{\wfe l'm'},\\
	R_{\wfe,lm;\wfm,l'm'} = \alpha_{\wfe lm}^{-1} B \alpha_{\wfm l'm'},\\
	R_{\wfm,lm;\wfe,l'm'} = \alpha_{\wfm lm}^{-1} B \alpha_{\wfe l'm'},\\
	R_{\wfm,lm;\wfm,l'm'} = \alpha_{\wfm lm}^{-1} A \alpha_{\wfm l'm'}.
\f]
	
If the coefficients for magnetic and electric waves are the same,
\f$ \alpha_{\wfm lm} = \alpha_{\wfe lm} \f$, the translation operator 
can be written in the same symmetric form as with the canonical convention,
just the matrices \f$ A, B\f$ will be different inside. 

If the coefficients differ (as in SCUFF-EM convention, where there
is a relative \a i -factor between electric and magnetic waves),
the functions such as qpms_trans_calculator_get_AB_arrays() will 
compute \f$ R_{\wfe\wfe}, R_{\wfe\wfm} \f$ for A, B arrays.


Literature convention tables
----------------------------

### Legendre functions and spherical harmonics

| Source                 | Ferrers function      | Negative \f$m\f$   | Spherical harmonics |
|------------------------|-----------------------|--------------------|---------------------|
| DLMF \cite NIST:DLMF   | \f[
	\dlmfFer{\nu}{m}\left(x\right)=(-1)^{m}\left(1-x^2\right)^{m/2}\frac{{
	\ud}^{m}\dlmfFer{\nu}{}\left(x\right)}{{\ud x}^{m}}
                                             \f] | \f[
	\dlmfFer{\nu}{m}(x) = (-1)^m \frac{\Gamma\pr{\nu-m+1}}{\Gamma\pr{\nu+m+1}}\dlmfFer{\nu}{m}(x)
                                                                  \f] |  Complex (14.30.1): \f[
		\dlmfYc{l}{m} = \sqrt{\frac{(l-m)!(2l+1)}{4\pi(l+m)!}} e^{im\phi} \dlmfFer{l}{m}(\cos\theta).
	\f] Real, unnormalized (14.30.2): \f$ 
		\dlmfYrUnnorm{l}{m}\pr{\theta,\phi} = \cos\pr{m\phi} \dlmfFer{l}{m}\pr{\cos\theta} 
	\f$ or \f$
		\dlmfYrUnnorm{l}{m}\pr{\theta,\phi} = \sin\pr{m\phi} \dlmfFer{l}{m}\pr{\cos\theta} 
                                                                                     \f$. |
| GSL \cite GSL         |  \f[
	\Fer[GSL]{l}{m} = \csphase^m N \rawFer{l}{m}
                                  \f]  for non-negative \f$m\f$. \f$
	\csphase\f$ is one by default and can be set to \f$
	-1\f$ using the functions ending with \_e with argument `csphase = -1`. \f$
	N\f$ is a positive normalisation factor from from `gsl_sf_legendre_t`. | N/A. Must be calculated manually. | The asimuthal part must be calculated manually. Use `norm = GSL_SF_LEGENDRE_SPHARM` to get the usual normalisation factor \f$
	N= \sqrt{\frac{(l-m)!(2l+1)}{4\pi(l+m)!}} \f$. |	
| Kristensson I \cite kristensson_spherical_2014 	| \f$ \rawFer{l}{m} \f$ | As in \f$ \rawFer{l}{m} \f$. | \f[
	\spharm[Kc]{l}{m} = (-1)^m \sqrt{\frac{(l-m)!(2l+1)}{4\pi(l+m)!}} \rawFer{l}{m}(\cos\theta) e^{im\phi},
                \f] (cf. Sec. D.2), therefore it corresponds to the DLMF sph. harms.: \f[ \spharm[Kc]{l}{m} = \dlmfYc{l}{m}. \f]  |
| Kristensson II \cite kristensson_scattering_2016	| \f$ \rawFer{l}{m} \f$ | As in \f$ \rawFer{l}{m} \f$. | \f[
	\spharm[Kr]{\begin{Bmatrix}e \\ o\end{Bmatrix}}{l}{m} = 
		\sqrt{2-\delta_{m0}}\sqrt{\frac{(l-m)!(2l+1)}{4\pi(l+m)!}}
		\rawFer{l}{m}(\cos\theta) 
		\begin{Bmatrix}\cos\phi \\ \sin\phi\end{Bmatrix},
                                                                                            \f] \f$ m \ge 0 \f$. Cf. Appendix C.3.  |
| Reid \cite reid_electromagnetism_2016	|  Not described in the memos. Superficial look into the code suggests that the `GetPlm` function *does* include the Condon-Shortley phase and spherical harmonic normalisation, so \f[
	\Fer[GetPlm]{l}{m} = (-1)^m \sqrt{\frac{(l-m)!(2l+1)}{4\pi(l+m)!}} \rawFer{l}{m}
\f] for non-negative \f$ m \f$.    |  N/A. Must be calculated manually.         |   \f[
	\spharm[GetYlm]{l}{m}(\theta,\phi) = \Fer[GetPlm]{l}{m}(\cos\theta) e^{im\phi},\quad m\le 0, \\
	\spharm[GetYlm]{l}{m}(\theta,\phi) = (-1)^m\Fer[GetPlm]{l}{\abs{m}}(\cos\theta) e^{-im\phi},\quad m<0, 
                                                   \f] and the negative sign in the second line's exponent is quite concerning, because that would mean the asimuthal part is actually \f$ e^{i\abs{m}\phi} \f$. _Is this a bug in scuff-em_? Without it, it would be probably equivalent to DLMF's \f$ \dlmfYc{l}{m} \f$s for both positive and negative \f$ m\f$s. However, it seems that `GetYlmDerivArray` has it consistent, with \f[
	\spharm[GetYlmDerivArray]{l}{m} = \dlmfYc{l}{m}   
		\f] for all \f$m\f$, and this is what is actually used in `GetMNlmArray` (used by  both `SphericalWave` in `libIncField` (via `GetMNlm`) and `GetSphericalMoments` in `libscuff` (via `GetWaveMatrix`)) and `GetAngularFunctionArray` (not used).      |



### VSWF conventions

| Source	| VSWF definition  	| E/M interrelations | VSWF norm  	| CS Phase  	|  Field expansion 	|  Radiated power | Notes |
|---	|---	|---	|---	|---	|---	|--- 	|--- |
| Kristensson I \cite kristensson_spherical_2014 	|  \f[ \wfkc = \dots \f] where \f$\wfkc\f$ is either of \f$ \wfkcreg, \wfkcout, \dots \f$ based on the radial (spherical Bessel) function type.	| \f[
	\wfkcreg_{1lm} = \frac{1}{k}\nabla\times\wfkcreg_{2lm}, \\
	\wfkcreg_{2lm} = \frac{1}{k}\nabla\times\wfkcreg_{1lm},
\f] and analogously for outgoing waves \f$ \wfkcout \f$, eq. (2.8) onwards. 	|  	| Yes, in the spherical harmonics definition, cf. sect. D.2.  	| \f[ 
	\vect E = k \sqrt{\eta_0\eta} \sum_n \left( \wckcreg_n  \wfkcreg_n + \wckcout_n \wfkcout_n  \right), 
	\\ 
	\vect H =  \frac{k \sqrt{\eta_0\eta}}{i\eta_0\eta} \sum_n \left( \wckcreg_n  \wfkcreg_n + \wckcout_n \wfkcout_n  \right),
\f] but for plane wave expansion \cite kristensson_spherical_2014 sect. 2.5 K. uses a different definition (same as in Kristensson II).  	| \f[
	 P = \frac{1}{2} \sum_n \left( \abs{\wckcout_n}^2 +\Re \left(\wckcout_n\wckcreg_n^{*}\right)\right)
 \f]	| The \f$ \wckcreg, \wckcout \f$	coefficients have dimension \f$ \sqrt{\mathrm{W}} \f$. |
| Kristensson II \cite kristensson_scattering_2016	| \f[ \wfkr = \dots \f] where \f$\wfkr\f$ is either of \f$ \wfkrreg, \wfkrout, \dots \f$ based on the radial (spherical Bessel) function type. 	|  \f[
	\nabla\times\wfkrreg_{\tau n} = k\wfkrreg_{\overline{\tau} n},
\f] eq. (7.7) and analogously for outgoing waves \f$ \wfkrout \f$. 	| 	|   	| \f[ 
	\vect E = \sum_n \left( \wckrreg_n  \wfkrreg_n + \wckrout_n \wfkrout_n  \right), 
	\\ 
	\vect H =  \frac{1}{i\eta_0\eta} \sum_n \left( \wckrreg_n  \wfkrreg_n + \wckrout_n \wfkrout_n  \right)
\f] 	| \f[
	 P = \frac{1}{2k^2\eta_0\eta} \sum_n \left( \abs{\wckrout_n}^2 +\Re \left(\wckrout_n\wckrreg_n^{*}\right)\right)
 \f]	| The \f$ \wckrreg, \wckrout \f$ coefficients have dimension \f$ \mathrm{V/m} \f$. |
| Reid \cite reid_electromagnetism_2016	| By examining the code, it appears that both `GetMNlmArray()` and `GetWaveMatrix()` with argument `MaxwellWaves = true` (with `MaxwellWaves = false` it seems to calculate nonsense) return the following w.r.t. Kristensson's "complex VSWFs": \f[
	\wfr_{lmM} = i\wfkc_{1lm}, \\
	\wfr_{lmN} = -\wfkc_{2lm}.
	\f] | \f[
	\nabla\times\wfr_{lmM} = -ik\wfr_{lmN}, \\ \nabla\times\wfr_{lmN} = +ik\wfr_{lmM}. 
\f] 	|	|  |  \f[
	\vect E = \sum_\alpha \pr{ \wcrreg_\alpha \wfrreg_\alpha + \wcrout_\alpha \wfrout_\alpha }, \\
	\vect H = \frac{1}{Z_0Z^r} \sum_\alpha \pr{ \wcrreg_\alpha \sigma_\alpha\wfrreg_\overline{\alpha} +
		 \wcrout_\alpha \sigma_\alpha\wfrout_\overline{\alpha}},
\f] where \f$ \sigma_{lmM} = +1, \sigma_{lmN}=-1, \overline{lmM}=lmN, \overline{lmN}=lmM, \f$  cf. eq. (6). The notation is not extremely consistent throughout Reid's memo.	| 	| 	|
| Taylor \cite taylor_optical_2011	| \f[
	\wfet_{mn}^{(j)}	=	\frac{n(n+1)}{kr}\sqrt{\frac{2n+1}{4\pi}\frac{\left(n-m\right)!}{\left(n+m\right)!}}\Fer[Taylor]{n}{m}\left(\cos\theta\right)e^{im\phi}z_{n}^{j}\left(kr\right)\uvec{r} \\
		+\left[\tilde{\tau}_{mn}\left(\cos\theta\right)\uvec{\theta}+i\tilde{\pi}_{mn}\left(\cos\theta\right)\uvec{\phi}\right]e^{im\phi}\frac{1}{kr}\frac{\ud\left(kr\,z_{n}^{j}\left(kr\right)\right)}{\ud(kr)}, \\ 
	\wfmt_{mn}^{(j)}	=	\left[i\tilde{\pi}_{mn}\left(\cos\theta\right)\uvec{\theta}-\tilde{\tau}_{mn}\left(\cos\theta\right)\uvec{\phi}\right]e^{im\phi}z_{n}^{j}\left(kr\right).
\f] Assuming the Legendre functions \f$ \Fer[Taylor]{n}{m} \f$ here do contain the Condon-Shortley phase (AFAIK not explicitly stated in the book), i.e. \f$\Fer[Taylor]{l}{m} = \dlmfFer{l}{m} \f$, then the relation to Kristensson's waves is \f[
	\wfmt_{mn} = \sqrt{n(n+1)} \wfkc_{1nm}, \\ \wfet_{mn} = \sqrt{n(n+1)} \wfkc_{2nm}. 
		\f]	|	|	\f[
	\int_{S(kr)} \wfmt_{mn}^{(j)} \wfmt_{m'n'}^{(j)}\,\ud S = n(n+1) \abs{z_n^{(j)}}^2 \delta_{m,m'}\delta_{n,n'} ,\\
	\int_{S(kr)} \wfet_{mn}^{(j)} \wfet_{m'n'}^{(j)}\,\ud S =
           \pr{\pr{n(n+1)}^2 \abs{\frac{z_n^{(j)}}{kr}}^2 + n(n+1)\abs{\frac{1}{kr}\frac{\ud}{\ud(kr)}\pr{kr z_n^{(j)}}} } \delta_{m,m'}\delta_{n,n'} ,
\f] cf. \cite taylor_optical_2011, eqs. (2.40–41). I suspect that this is also wrong and \f$\delta_{m,m'}\f$ should be replaced with \f$\delta_{m,-m'}\f$. |	| \f[ 
	\vect E = \sum_{mn} \pr{-i \pr{\wcetreg_{mn}\wfetreg_{mn} + \wcmtreg_{mn}\wfmtreg{mn}} +i \pr{\wcetout_{mn}\wfetout_{mn} + \wcmtout_{mn}\wfmtout_{mn}}}, \\
	\vect H = n_{ext}\sum_{mn} \pr{- \pr{\wcmtreg_{mn}\wfetreg_{mn} + \wcetreg_{mn}\wfmtreg{mn}} + \pr{\wcmtout_{mn}\wfetout_{mn} + \wcetout_{mn}\wfmtout_{mn}}},
\f] 	| 	| Different sign for regular/scattered waves! Also WTF are the units of \f$ n_{ext} \f$?  The whole definition seems rather inconsistent. |

