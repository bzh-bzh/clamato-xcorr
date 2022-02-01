
; This file includes the subroutines and procedures necessary to do a
; basic PCA fit on to SDSS/BOSS quasar spectra. 
;
; Requires: 
; - Nao Suzuki's HST eigenspectra and Isabelle Paris' SDSS DR7
;   eigenspectra (xi_hst.txt and xi_sdss.txt)... XIDIR in READPCA_HST
;   and READPCA_SDSS need to point to a directory where
;   these are located.
; - David Schlegel's IDL SPECTRO utilities (specifically READSPEC)
; - Craig Markwardt's MPFIT fitting routines
; 
; Quickstart: 
; For plug-in functionality, use PCA_CHISQ to get the best-fit PCA
; parameters for a given spectrum, then PCA_FUNC to generate the
; corresponding spectrum. For more generality, use the wrapper
; function PCASPEC_FUNC in which also includes mean-flux regulation.
;
; Contact Khee-Gan Lee (lee@astro.princeton.edu) for questions/comments
;
; Revision history:
; KG Lee 08/12/2010 - Finished preliminary coding.
; KG Lee 09/14/2010 - Cleaned up for beta release to
;                     collaboration. Current version does *not* have
;                     decent metal line  masking 
; KG Lee 12/03/2010 - Modified to compute for BOSS VAC3 sightlines:
;                   1. Reads in Anze's cleaned VAC3 catalog
;                   (qsovac10101.txt) and reads spectrum from local
;                   Princeton BOSS repository  
; KG Lee 12/08/2010 - Turned off line masking for initial BOSS fits.

; KG Lee 04/17/2011 - Modified and cleaned up PCA_CHISQ extensively
;                     with line-masking iterations. Simplified
;                     ABSFLUXERR.   
; KG Lee 06/10/2011 - Modified from pca_chisq.pro, to enable PCA_FUNC
;                     to return either Suzuki+ or Paris+ fits. Only
;                     modifications are to PCA_FUNC and READPCA_SDSS 
; KG Lee 08/30/2011 - Added PCA_FUNC_DR7 to carry out fits with SDSS
;                     eigenspectra without hardcoding. PCA_CHISQ
;                     modified to accommodate. 
; KG Lee 10/20/2011 - Parameter constraints limited to 2-sigma in
;                     PARAM_CONSTRAINTS
; KG Lee 12/22/2011 - READPCA_SDSS and READPCA_HST now loads
;                     eigenspectra from specific directory rather than
;                     expecting them to be local in the calling
;                     directory 
; KG Lee 30/12/2011 - Added option to input inverse variances instead
;                     of errors. Disabled option to calculate
;                     chi-squared. 
; KG Lee  3/7/2012  - Fixed bug with ~1280A normalization region when
;                     inverse variances are used 

function	smoothflux,flux,sig
;	Smooths the flux by a window of width sig pixels.
nn = size(flux,/n_elem)
kk = lindgen(nn)
w  = where(kk GT nn/2)
kk[w]=nn-kk[w]
kk = 2*!pi*kk/nn
fk = FFT(flux,1,/double)*exp(-0.5*(kk*sig)^2)
sm = FFT(fk,-1,/double)
return,float(sm)
end

FUNCTION CALC_CHISQ, ff, sigma, lambda, pcaparams, dof=dof, $
                     masklist=masklist
; Evaluate the reduced chi-squared from the best-fit PCA parameters
; (using PCA_FUNC) and the observed spectrum. Carry this out redwards
; of Lya. 
;
; Inputs:
; FF : Observed flux (ideally already masked)
; SIGMA : Errors on the flux
; LAMBDA: Wavelength corresponding to ff and sigma.
; PCAPARAMS: Vector of best-fit PCA parameters as returned from PCA_CHISQ 
;
; Returns:
; Reduced chi-squared value between the observed spectrum and the model

lambda_red = lambda[where(lambda GT 1220. and lambda LE 1600.)]
ff_red = ff[where(lambda GT 1220. and lambda LE 1600.)]
sigma_red = sigma[where(lambda GT 1220. and lambda LE 1600.)]

cz = pcaparams[0]

ff_model = pca_func(lambda_red, pcaparams)

;print, 'dof = '+strtrim(n_elements(ff_red),2)+'-'+strtrim(n_elements(pcaparams),2)

if not keyword_set(dof) then dof = n_elements(ff_red) - $
  n_elements(pcaparams)

chisq = (ff_red - ff_model)^2 / sigma_red^2

;print, 'chisq = ', strtrim(total(chisq), 2)

chisq_dof = total(chisq) / dof

return, chisq_dof
END

FUNCTION ABSFLUXERR, ff_in, lambda, pcaparams, dr7eigen=dr7eigen
; Computes absolute flux error between reconstruction and prediction,
; \int (ff_obs - ff_model)/ff_model dlambda / \int dlambda
; c.f. Eq 11 from Suzuki et al 2005.
; *Does not take noise into account*

; Inputs:
; FF : Observed flux (ideally already masked)
; LAMBDA: Wavelength corresponding to ff and sigma.
; PCAPARAMS: Vector of best-fit PCA parameters as returned from PCA_CHISQ 
cz = pcaparams[0]

ff = ff_in

ff = smooth(ff,15)

lambda_red = lambda[where(lambda GT 1230.*cz and lambda LE 1600.*cz)]
ff_red = ff[where(lambda GT 1230.*cz and lambda LE 1600.*cz)]

ff_model = pca_func(lambda_red, pcaparams,dr7eigen=dr7eigen)

ff_dev = abs((ff_model/ff_red) - 1.)

return, avg(ff_dev)
END

;******************************************************************


PRO READPCA_HST, mmax_in=mmax_in, plotcomp=plotcomp
; Read Nao Suzuki's HST PCA eigenspectra (Suzuki et al 2005) into
; common block, up to mmax components
;
; Optional Keywords:
; mmax : Max number of components required (up to 10)
; plotcomp: Specify component to plot to x-window. 
;
; Reads to PCABLOCK:
; lambxi  - restframe wavelength
; mu_spec - mean quasar spectrum from the same, i.e. the 0th component
; sig_spec - standard deviation of the original sample
; xi_pca - mmax x npix array of eigenspectra
; m - number of components stored in the block

COMMON PCABLOCK, lambxi, mu_pca, sig_pca, xi_pca, mmax

if not keyword_set(mmax_in) then mmax_in = 10

xidir = getenv('MF_PCA_DIR')+'/data/'

xihst_file = xidir+'xi_hst.txt' 

readcol, xihst_file, lambxi, mu_pca, sig_pca, xi1, xi2, xi3, $
  xi4, xi5, xi6, xi7, xi8, xi9, xi10, skipline=15, /silent
 
npix = n_elements(lambxi)

xi_tmp = fltarr(10, npix)
xi_pca = fltarr(mmax_in, npix)
mmax = mmax_in

xi_tmp[0,*] = xi1
xi_tmp[1,*] = xi2
xi_tmp[2,*] = xi3
xi_tmp[3,*] = xi4
xi_tmp[4,*] = xi5
xi_tmp[5,*] = xi6
xi_tmp[6,*] = xi7
xi_tmp[7,*] = xi8
xi_tmp[8,*] = xi9
xi_tmp[9,*] = xi10

if keyword_set(plotcomp) then begin
    plot, lambxi, xi_tmp[(plotcomp-1),*]
endif

xi_pca = xi_tmp[0:(mmax-1),*]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO READPCA_SDSS, mmax_in=mmax_in, plotcomp=plotcomp
; Read Isabelle Paris' SDSS PCA eigenspectra (Paris et al, 2011) into
; common block, up to mmax components
;
; Optional Keywords:
; mmax : Max number of components required (up to 10)
; plotcomp: Specify component to plot to x-window. 
;
; Reads to PCABLOCK:
; lambxi  - restframe wavelength
; mu_spec - mean quasar spectrum from the same, i.e. the 0th component
; sig_spec - standard deviation of the original sample
; xi_pca - mmax x npix array of eigenspectra
; m - number of components stored in the block

COMMON PCABLOCK_SDSS, lambxi_sdss, mu_pca_sdss, sig_pca_sdss, $
  xi_pca_sdss, mmax_sdss

if not keyword_set(mmax_in) then mmax_in = 10

xidir = getenv('MF_PCA_DIR')+'/data/'

xisdss_file = xidir+'xi_sdss.txt'

readcol, xisdss_file, lambxi_sdss, mu_pca_sdss, xi1, xi2, xi3, $
  xi4, xi5, xi6, xi7, xi8, xi9, xi10, /silent

sig_pca_sdss = fltarr(n_elements(mu_pca_sdss))

npix = n_elements(lambxi_sdss)

xi_tmp = fltarr(10, npix)
xi_pca_sdss = fltarr(mmax_in, npix)
mmax_sdss = mmax_in

xi_tmp[0,*] = xi1
xi_tmp[1,*] = xi2
xi_tmp[2,*] = xi3
xi_tmp[3,*] = xi4
xi_tmp[4,*] = xi5
xi_tmp[5,*] = xi6
xi_tmp[6,*] = xi7
xi_tmp[7,*] = xi8
xi_tmp[8,*] = xi9
xi_tmp[9,*] = xi10

; For my purposes, want only 1020-1600A
wavecut = where(lambxi_sdss LE 1600.)

lambxi_sdss = lambxi_sdss[wavecut]
mu_pca_sdss = mu_pca_sdss[wavecut]
xi_tmp = xi_tmp[*,wavecut]
sig_pca_sdss = sig_pca_sdss[wavecut]


if keyword_set(plotcomp) then begin
    plot, lambxi_sdss, xi_tmp[(plotcomp-1),*]
endif

xi_pca_sdss = xi_tmp[0:(mmax_sdss-1),*]
end


;----------------------------------------------------------------------

function PARAM_CONSTRAINTS, fixz=fixz, dr7eigen=dr7eigen
; Returns a PARINFO structure to constrain the range of the parameters
; in the MPFIT routine
; 
; The first n_notpca parameters are not PCA weights
;
; Optional keyword:
; mmax - number of PCA components to set the structure for
; 
; KG Lee 08/12/2010: 
; 0th parameter - redshift correction factor, unconstrained
; 1st parameter - normalization factor, unconstrained
; 2nd parameter - powerlaw index, constrained to +-1
; The rest - PCA components, constrained to within roughly 3-4 sigma
;            of their eigenvalues as listed in Suzuki 2006 
;
; KG Lee 08/29/2010: Added fixz keyword, which fixes the redshift.
; KG Lee 06/11/2011: Added separate set of constraints for Paris+11
;                    eigenspectra 

COMMON PCABLOCK
COMMON PCABLOCK_SDSS

if keyword_set(dr7eigen) then mmax_tmp = mmax_sdss else $
  mmax_tmp = mmax

n_notpca = 3

paramconst = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]}, $
                       (mmax_tmp+n_notpca))

; Fix redshifts for use with Hewett & Wild 2010 redshifts
if keyword_set(fixz) then paramconst[0].fixed = 1

; Constrain cz
paramconst[0].limited[0] = 1
paramconst[0].limited[1] = 1
paramconst[0].limits[0] = 0.9
paramconst[0].limits[1] = 1.1

; Power law is supposed to be just a correction, so limit it to a
; relatively small value
paramconst[2].limited[0] = 1
paramconst[2].limited[1] = 1
paramconst[2].limits[0] = -1.
paramconst[2].limits[1] = 1.
; This tells MPFIT to constrain the values of the PCA weights
paramconst[n_notpca:mmax_tmp+n_notpca-1].limited[0] = 1
paramconst[n_notpca:mmax_tmp+n_notpca-1].limited[1] = 1

if not keyword_set(dr7eigen) then begin 
; And these are the limits we want for the cij, based on the
; distribution shown in Suzuki 2006
if mmax_tmp GE 1 then begin
paramconst[n_notpca].limits[0] = -15.
paramconst[n_notpca].limits[1] = 15.
endif
if mmax_tmp GE 2 then begin
paramconst[n_notpca+1:mmax_tmp+n_notpca-1].limits[0] = -7.2
paramconst[n_notpca+1:mmax_tmp+n_notpca-1].limits[1] = 7.2
endif
if mmax_tmp GE 3 then begin
paramconst[n_notpca+2:mmax_tmp+n_notpca-1].limits[0] = -5.
paramconst[n_notpca+2:mmax_tmp+n_notpca-1].limits[1] = 5.
endif
if mmax_tmp GE 5 then begin
paramconst[n_notpca+4:mmax_tmp+n_notpca-1].limits[0] = -3.
paramconst[n_notpca+4:mmax_tmp+n_notpca-1].limits[1] = 3.
endif
if mmax_tmp GE 7 then begin
paramconst[n_notpca+6:mmax_tmp+n_notpca-1].limits[0] = -2.
paramconst[n_notpca+6:mmax_tmp+n_notpca-1].limits[1] = 2.
endif
if mmax_tmp GE 8 then begin
paramconst[n_notpca+7:mmax_tmp+n_notpca-1].limits[0] = -1.5
paramconst[n_notpca+7:mmax_tmp+n_notpca-1].limits[1] = 1.5
endif

endif else begin

if mmax_tmp GE 1 then begin
paramconst[n_notpca].limits[0] = -7.
paramconst[n_notpca].limits[1] = 5.
endif
if mmax_tmp GE 2 then begin
paramconst[n_notpca+1:mmax_tmp+n_notpca-1].limits[0] = -1.6
paramconst[n_notpca+1:mmax_tmp+n_notpca-1].limits[1] = 1.6
endif
if mmax_tmp GE 5 then begin
paramconst[n_notpca+3:mmax_tmp+n_notpca-1].limits[0] = -1.
paramconst[n_notpca+3:mmax_tmp+n_notpca-1].limits[1] = 1.
endif
if mmax_tmp GE 7 then begin
paramconst[n_notpca+6:mmax_tmp+n_notpca-1].limits[0] = -0.65
paramconst[n_notpca+6:mmax_tmp+n_notpca-1].limits[1] = 0.65
endif

endelse   

return, paramconst
end

;-----------------------------------------------------------------


function PCAWEIGHTS,ff_in, lambda_in, interpflag=interpflag, $
                    dr7eigen=dr7eigen
; For input spectrum ff and sigma as a function of quasar restframe
; wavelength lambda_r, estimate the PCA weights c_ij from the xi_j PCA
; components by doing matrix integral of (f_lambda,j-mu)*xi_j over the
; portion of  spectrum redwards of Lya emission. Carry this out on the
; smoothed  spectrum. 
;
; *Does not take spectral noise into account*

; Make estimates up to mmax PCA components (default mmax=10)

; Inputs: 
; ff          = Flux spectrum of the quasar
; lambda_r    = Wavelength in restframe of quasar
; mmax        = Number of PCA components to evaluate. Defaults to 10
; interpflag  = Flag to determine interpolation method to rebin the
;               eigenspectra to lambda_r.
;               interpflag=0: Linear interpolation
;               interpflag=1: Spline interpolation

; Common block PCABLOCK needs to be initialized before calling this
; function

; Returns: n=mmax vector of PCA weights
;
; Revisions:
; Khee-Gan Lee 08/20/2010 -  Hacked version that needs to be fixed
; Khee-Gan Lee 08/24/2010 -  Amended to carry it out as a matrix
;                          multiplication 
; Khee-Gan Lee 04/15/2011 - Fixed bug which smoothed the the input
;                           spectrum 

COMMON PCABLOCK
COMMON PCABLOCK_SDSS

if keyword_set(dr7eigen) then begin
    lambxi_tmp=lambxi_sdss
    mu_pca_tmp = mu_pca_sdss
    sig_pca_tmp = sig_pca_sdss
    xi_pca_tmp=xi_pca_sdss
    mmax_tmp = mmax_sdss
endif else begin
    lambxi_tmp=lambxi
    mu_pca_tmp = mu_pca
    sig_pca_tmp = sig_pca
    xi_pca_tmp=xi_pca
    mmax_tmp = mmax
endelse 

ff = ff_in
lambda_r = lambda_in

ff = smooth(ff,10)
; Make sure it's normalized
normpix = where(lambda_r GE 1275.  AND lambda_r LE 1285.)
if normpix NE [-1] then $
  normfac = AVG(ff[where(lambda_r GE 1275. $
                         AND lambda_r LE 1285.)]) else $
  normfac = AVG(ff[where(lambda_r GE 1450. $
                         AND lambda_r LE 1470.)])

ff = ff / normfac

if not keyword_set(interpflag) then interpflag=0

; The wavelength range redwards of Lya to fit the PCA weights
lamb_red = lambda_r[where(lambda_r GE 1210. AND lambda_r LT 1600.)]
ff_red = ff[where(lambda_r GE 1210. AND lambda_r LT 1600.)]

; Interpolate the eigenspectra to the lamb_red grid
if interpflag EQ 0 then begin
    mu =  interpol(mu_pca_tmp, lambxi_tmp, lamb_red)
endif else begin
    mu = spline(lambxi_tmp, mu_pca_tmp, lamb_red)
endelse

xi = fltarr(mmax_tmp, n_elements(lamb_red))
for m=0L, mmax_tmp-1 do begin
    if interpflag EQ 0 then begin
        xi[m,*] = interpol(xi_pca_tmp[m,*], lambxi_tmp, lamb_red)
    endif else begin
        xi[m,*] =  spline(lambxi_tmp, xi_pca_tmp[m,*], lamb_red)
    endelse
endfor

; Now integrate over the deviations from the mean to get the PCA
; weights
; Now done via matrix multiplication
c_ij = fltarr(mmax_tmp)
resid = reform(ff_red - mu)

for m=0L, mmax_tmp-1 do begin
    xi_tmp = reform(xi[m,*])
    c_ij[m] = transpose(xi_tmp) # resid
    ;c_ij[m]=int_tabulated(lamb_red, (ff_red - mu)*xi[m,*]) 
endfor

return, c_ij
end

;---------------------------------------------------------------------
function PCA_FUNC, x, p, dr7eigen=dr7eigen
; For a given value of wavelength, return normalized flux for the
; corresponding set of PCA weights. To be used in conjunction with
; MPFITFUN 

; PCABLOCK needs to be initialized beforehand.
; 
; Optional Keywords:
;   DR7EIGEN  - If set, use Paris+11 DR7 eigenspectra. This only works
;               for recreating spectra... during fitting this function
;               will revert to its normal use
;
; Inputs:
; x: wavelength in the quasar restframe
; p[0]: correction factor to multiply the wavelength; used to fit for
; new redshift
; p[1]: Correction to flux normalization relative to flux at 1280A
; p[2]: Power law component
; p[3:mmax+1]: weights for PCA eigenvectors
;
; Return normalized flux corresponding to the wavelength

COMMON PCABLOCK
COMMON PCABLOCK_SDSS

if keyword_set(dr7eigen) then begin

    lambxi_tmp=lambxi_sdss
    mu_pca_tmp = mu_pca_sdss
    sig_pca_tmp = sig_pca_sdss
    xi_pca_tmp=xi_pca_sdss
    mmax_tmp = mmax_sdss
endif else begin
    lambxi_tmp=lambxi
    mu_pca_tmp = mu_pca
    sig_pca_tmp = sig_pca
    xi_pca_tmp=xi_pca
    mmax_tmp = mmax
endelse 

cij = p[3:*]
cz = p[0]
fnorm = p[1]
alpha = p[2]

weightedpcs = xi_pca_tmp

for mm = 0, mmax_tmp-1 do begin
    weightedpcs[mm,*] = cij[mm] * xi_pca_tmp[mm,*] 
endfor

sum_pcs = total(weightedpcs,1)

pcaspec =  (mu_pca_tmp + sum_pcs)

pcaspec = pcaspec / avg(pcaspec[where(lambxi_tmp GE 1275. AND $
                                      lambxi_tmp LE 1285.)])

; Now interpolate for value of normalized flux corresponding to input
; wavelength * cz
ff = fnorm * interpol(pcaspec, lambxi_tmp, cz*x) * (cz*x/1280.)^alpha

return, ff
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
function PCA_FUNC_DR7, x, p
; this version of PCA_FUNC is specifically intended to be called from
; PCA_CHISQ during fitting with Paris+11 eigenspectra. 

COMMON PCABLOCK_SDSS 

lambxi_tmp=lambxi_sdss
mu_pca_tmp = mu_pca_sdss
sig_pca_tmp = sig_pca_sdss
xi_pca_tmp=xi_pca_sdss
mmax_tmp = mmax_sdss

cij = p[3:*]
cz = p[0]
fnorm = p[1]
alpha = p[2]

weightedpcs = xi_pca_tmp

for mm = 0, mmax_tmp-1 do begin
    weightedpcs[mm,*] = cij[mm] * xi_pca_tmp[mm,*] 
endfor

sum_pcs = total(weightedpcs,1)

pcaspec =  (mu_pca_tmp + sum_pcs)

pcaspec = pcaspec / avg(pcaspec[where(lambxi_tmp GE 1275. AND $
                                      lambxi_tmp LE 1285.)])

; Now interpolate for value of normalized flux corresponding to input
; wavelength * cz
ff = fnorm * interpol(pcaspec, lambxi_tmp, cz*x) * (cz*x/1280.)^alpha

return, ff
end 

;---------------------------------------------------------------------

function PCA_CHISQ, ff, lambda_r, sigma, ivar=ivar_in, dof=dof, $
                    afluxerr=afluxerr,quiet=quiet, $
                    contiter=contiter, maskwave=maskwave, $
                    dr7eigen=dr7eigen, m_fit=m_fit1
                    
; For an input rest-frame quasar spectrum and associated errors, 
; 1. Mask out absorption lines
; 2. Estimate PCA weights by integrating over the flux
; 3. Carry out Chi-squared minization over the redshift correction, flux
;    normalization, power-law index (alpha_lambda), and PCA weights,
;    in 3 steps:
;    a) Using the first-guess PCA weights as input, do the chi-squared
;    fit. This provides an improved redshift.
;    b) With the redshift correction, select pixels near the emission
;    lines to be given extra weight. 
;    c) Now do the chi-squared fit again to get improved emission line
;    shapes 
; 4. Iterate twice, removing points more than 2 sigma below the
;    continuum to avoid absorption lines
; 
; Uses the MPFIT routine
;
; PCABLOCK should be initialized beforehand using the READ_HST procedure
;
; Returns: 
; Vector of fit parameters in the sequence: redshift
; correction factor, flux normalization, powerlaw index (alpha_lambda)
; then best-fit PCA weights. *The fitted parameters are with respect
; to an observed spectrum normalized at 1280A in the restframe of the
; pipeline redshift*

; Optional keywords:
; IVAR       - Inverse variances, entered in lieu of errors
; DOF        - Degrees of freedom in the fit, # of points fitted - #
;              fitted  parameters
; MASKWAVE - List of wavelengths indicating pixels that have been
;              masked 
; AFLUXERR - Returns absolute flux error value as evaluated by
;              ABSFLUXERR 
; CONTITER - Returns 3xnpix_red array of wavelength and 2 continuum
;            iterations during line-masking routine.
; M_FIT    - Limits number of PCA components used in the fitting. If
;            M_FIT = 0 (i.e use only mean spectrum) is required, set
;            it to a floating point number that rounds to zero,
;            e.g. M_FIT = 0.1
; 
; Revision history:
; Khee-Gan Lee (08/20/2010) - Added evaluations of chi-squared and
;                             absolute flux error in the fit. Still
;                             needs to be debugged and cleaned.  
; Khee-Gan Lee (04/06/2011) - Line masking now carried out by 2-step
;                             process, by throwing out pixels more
;                             than nu s.d. from first PCA fit
; Khee-Gan Lee (04/14/2011) - Added iteration loops to
;                             line-masking. Also fixed bug in
;                             PCAWEIGHTS which returned a smoothed
;                             spectrum. Removed weights on emission
;                             lines, and now fits from 1216 onwards.
; Khee-Gan Lee (05/03/2011) - Added hack of increasing the noise of
;                             high-SN objects to improve the fits :/
; Khee-Gan Lee (08/30/2011) - Now returns PCA fit with the same flux
;                             normalization as the input spectrum
; Khee-Gan Lee (09/05/2011) - Added M_FIT, which limits the number of
;                             PCA components used in the fits by
;                             setting components beyond that to zero,
;                             and fixing it in MPFIT. 
; Khee-Gan Lee (10/20/2011) - Loads PCA templates into common block
;                             automatically 

COMMON PCABLOCK
COMMON PCABLOCK_SDSS

; Load necessary PCA templates from common block.
if keyword_set(dr7eigen) AND not keyword_set(xi_sdss) then $
  readpca_sdss, m=8

if not keyword_set(dr7eigen) AND not keyword_set(xi_hst) then $
  readpca_hst, m=8

if keyword_set(dr7eigen) then mmax_tmp = mmax_sdss else $
  mmax_tmp = mmax
if keyword_set(m_fit1) then m_fit=m_fit1

if not keyword_set(ivar_in) AND not keyword_set(sigma) then $
   print, 'Error: Need to set either SIGMA or IVAR'

if not keyword_set(ivar_in) then ivar = $
   (sigma NE 0.)/(sigma + (sigma EQ 0.))^2 else ivar = ivar_in

; Normalize flux by the region around 1280A or 1450A restframe
normpix1 = where(lambda_r GE 1275.  AND lambda_r LE 1285.)
if normpix1 NE [-1] AND total(ivar[normpix1]) NE 0 then $
  normfac = AVG(ff[normpix1]*ivar[normpix1]) / $
  avg(ivar[normpix1]) $
else begin
    normpix2 = where(lambda_r GE 1450.  AND lambda_r LE 1470.)
    normfac = AVG(ff[normpix2]*ivar[normpix2])/ $
      avg(ivar[normpix2])
endelse

ff_norm = ff / normfac
ivar_norm = ivar * normfac^2

; Select pixels redwards of Lya
redlist = where(lambda_r GE 1216. AND lambda_r LT 1600.)
lambda_red = lambda_r[redlist]
ff_red = ff_norm[redlist]

ivar_red = ivar_norm[redlist]

if keyword_set(sigma) then begin
   sigma_norm= sigma/normfac
   sigma_red = sigma_norm[redlist]
endif

snmed = median(ff_red* sqrt(ivar_red))

; Hack: if S/N is more than 8, then increase noise since high-SN
; objects seem to have crappier fits
if snmed GT 8. then begin
   ivar_red = ivar_red * (8./snmed)^2
   if keyword_set(sigma_red) then sigma_red = sigma_red * snmed / 8.
end

; Now compute PCA weights for first guess. Avoid the actual Lya peak
; which might be absorbed
ff_weights = ff_red
cij = PCAWEIGHTS(ff_weights, lambda_red,interpflag=1, $
                 dr7eigen=dr7eigen) 

;print, 'initial guess: ', cij

; Use MPFIT to fit for best-fit PCA weights cij, redshift correction
; factor cz and flux normalization
n_notpca = 3

; Set up constraints on the parameters, to be passed to PARINFO
; keyword in MPFIT
paramconst = param_constraints(dr7eigen=dr7eigen);,/fixz)

; In case the parameters from the initial guess are beyond the
; parameter bounds, set the guesses to the upper/lower limits
cij_constraint = paramconst[n_notpca:*].limits[1]

for mm=0, mmax_tmp-1 do begin
    if abs(cij[mm]) GT cij_constraint[mm] then begin
        signconst = 1.
        if cij[mm] LT 0 then signconst = -1.
        cij[mm] = cij_constraint[mm] * signconst*0.98
    endif
endfor

; If M_FIT is set, that means only want to fit up to M_FIT PCA
; components  
if keyword_set(m_fit) then begin
    m_fit = round(m_fit)
    if m_fit LT mmax then begin
        if m_fit EQ 0 then cij=replicate(0., mmax) else $
          cij[m_fit-1:*] = 0.
        paramconst[m_fit+3:*].fixed = 1
    endif
endif

; Parameters for PCA_FUNC: 0=redshift correction, 1=normalization
; factor, 2=powerlaw, 3 to (3+mmax_tmp-1): PCA weights
inparam= fltarr(mmax_tmp+n_notpca)
inparam[0] = 1.
inparam[1] = 1.
inparam[2] = 0.
inparam[n_notpca:mmax_tmp+n_notpca-1] = cij

;print, paramconst.fixed

; FIRST LINE MASKING ITERATION

; Was the inverse variance specified by user? If so then use inverse
; variances in MPFIT. Also define a dummy SIGMA_RED vector for
; brevity, which will be ignored by MPFIT
if keyword_set(ivar_in) then begin
   ivar_mpfit = ivar_red
   sigma_red = fltarr(n_elements(ivar_red))
endif
if keyword_set(dr7eigen) then pcafunc_str = 'PCA_FUNC_DR7' $
  else pcafunc_str = 'PCA_FUNC'

pcafit = mpfitfun(pcafunc_str, lambda_red, ff_red, sigma_red, $
                  inparam, PARINFO=paramconst,/quiet,WEIGHTS=ivar_mpfit) 

;print, '# of params: ', n_elements(pcafit)
pcaspec_first = pca_func(lambda_red, pcafit, dr7eigen=dr7eigen)

; Now should have decent redshift correction, so now choose the pixels
; in the vicinity of metal lines for extra weight 
cz = pcafit[0]

; Remove absorption lines by throwing out pixels which are more than
; nu_abs standard deviations from the first continuum-fit.
; For higher-S/N spectra, rescale nu_abs

; Criterion for throwing out pixels
nu_abs = 2.5                     
; For high-sn pixels, impose more stringent criteria to avoid throwing
; pixels due to imperfect fitting
 if snmed GT 6. then nu_abs = nu_abs * snmed/ 6.

abslist = where(ff_red - pcaspec_first LT -nu_abs/sqrt(ivar_red))
; Remove masked lines from copy of vectors
    lambda2 = lambda_red
    ff2 = ff_red
    ivar2 = ivar_red
    sigma2 = sigma_red
    
if abslist NE [-1] then remove, abslist, lambda2,ff2, ivar2, sigma2

; Criteria for loop iterations
crit_delta_fits = 0.02
n_iter_max = 5
delta_fits = 1.
n_iter = 0

pcaspec_init = pcaspec_first

; LINE MASKING ITERATIONS************************
while n_iter LT n_iter_max AND delta_fits GT crit_delta_fits do begin
; Redo fit with absorption from previous iteration removed
    pcafit2 = mpfitfun(pcafunc_str, lambda2, ff2, sigma2, $
                       pcafit, PARINFO=paramconst,/quiet, $
                       dof=dof, weights=ivar2,bestnorm=bestnorm)
    pcaspec_second = pca_func(lambda_red, pcafit2, $
                             dr7eigen=dr7eigen)
    
; Absolute residual between this fit and previous iterations
    delta_fits = avg(abs(pcaspec_second/pcaspec_first - 1.))
;    print, 'delta_fits = ', delta_fits
; Redo the metal line masking with the improved continuum
    abslist2 = where(ff_red - pcaspec_second LT -nu_abs/sqrt(ivar_red))

; Remove masked lines from copy of vectors
    lambda2 = lambda_red
    ff2 = ff_red
    sigma2 = sigma_red
    ivar2 = ivar_red
    if abslist2 NE [-1] then maskwave = lambda_red[abslist2]
    if abslist2 NE [-1] then remove, abslist2, lambda2,ff2, sigma2, ivar2

    n_iter++
    pcaspec_first = pcaspec_second
endwhile
; DONE WITH LINE MASKING ********************************

pcafit_fin = pcafit2
; Rescale flux normalization to that of the input spectrum
pcafit_fin[1] = pcafit_fin[1] * normfac

;chisq_dof = calc_chisq(ff_red, sigma_red, lambda_red, pcafit_fin)
afluxerr = absfluxerr(ff_norm[redlist], lambda_r[redlist], pcafit_fin,dr7eigen=dr7eigen)

;if not keyword_set(quiet) then begin
;    print, 'Chi-squared/d.o.f. = ', strtrim(chisq_dof,2)
;endif 


contiter = fltarr(3, n_elements(redlist))
contiter[0,*] = lambda_r[redlist]
contiter[1,*] = pcaspec_init
contiter[2,*] = pcaspec_second

; maskwave isn't set by now, then it's null
if not keyword_set(maskwave) then maskwave = [-1]

return, pcafit_fin 
end  

; ----------------------------------------------------------------
