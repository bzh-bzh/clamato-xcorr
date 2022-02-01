


;
; Contact Khee-Gan Lee at lee@astro.princeton.edu for questions

; This file includes:
; RESAMPLE.PRO      - Code to support bootstrap resampling, poached
;                     off internet
; MFLUXCORR         - Functional form of fitting function for
;                     mean-flux regulation
; MFLUX_TAUEVO_FUNC - Wrapper to invert MFLUXCORR during mean-flux
;                     regulation 
; CALC_MEANFLUX     - Evaluate mean-flux in rest-frame bins, and
;                     returns error estimates by bootstrapping pixels
;                     within each bin
; MF_PCA            - Carries out mean-flux regulation on an input
;                     quasar restframe spectrum. Uses PCA_CHISQ to
;                     carry out PCA fit.
;
; Revision history: V
; KG Lee 09/14/10 -  Cleaned-up beta version for release to collaboration
; KG Lee 12/04/10 - Hacked version to try stuff out for BOSS
;                   scripts. plot_mffit now plots in observed frame. 
; KG Lee 12/08/10 -  Fixed mistake in multiplying by redshift
;                   correction cz when converting back to observed
;                   frame. Also fixed bug in pcaspec_func so that
;                   mean-flux correction is applied only bluewards of
;                   1280A.
; KG LEe 04/19/11 - Cleaned up version, added quadratic fitting func.
; KG Lee 02/11/11 - Fitting function changed to linear function
;                   pivoting at 1113A, still with 2 free parameters 
; KG Lee 03/07/12 - Fixed normalization bug when all pix in 1280A are
;                   masked and ivar's are used.


;------------------------------------------------------------------------


;-------------------------------------------------------------------------

function MF_PCA, ff_in, lambda_r_in, sigma_in, ivar=ivar_in, $
                 delta_mf=delta_mf, pcaparams=pcaparams, $
                 afluxerr=afluxerr, dr7eigen=dr7eigen, $
                 m_fit=m_fit, forest_range=forest_range
;+
; NAME:
;   MF_PCA
;
; PURPOSE: 
; For input QSO restframe spectrum, carry out mean-flux regulated PCA
; fit. Formerly known as PCA_FOREST
; 
; Doesn't actually return anything, but use keywords to returns
; the corresponding delta_mf (exp correction factor) and
; pcaparams, the best-fit PCA parameters. 
;
; See PCASPEC_FUNC for a wrapper for this. 
;
; CALLING SEQUENCE:
;  mf_pca = mf_pca(ff_in, lambda_r_in, sigma_in, ivar=ivar_in, $
;                  delta_mf=delta_mf, pcaparams=pcaparams, $
;                  afluxerr=afluxerr, dr7eigen=dr7eigen, $
;                  m_fit=m_fit, forest_range=forest_range)
;
; Inputs:
; ff : Quasar flux
; lambda_r: Restframe Wavelength
;
; Optional inputs:
; sigma: Flux Errors
; 
; Optional Keywords:
; ivar     - Inverse variances, input in lieu of SIGMA_IN
; delta_mf: if set, returns the meanflux correction parameter, where
;           Exp(-delta_mf*(lambda/1280  -1))  
; pcaparams: returns, the best-fit parameters from PCA_CHISQ
; AFLUXERR - Returns absolute flux error as provided by PCA_CHISQ
; FOREST_RANGE - 2-element vector defining the low- and upper-limit of
;                the Lya forest range 

; KG Lee 08/12/2010 - Adapted from PLAW_FOREST
; KG Lee 08/22/2010 - Added redshift correction to the forest
;                    extraction 
; KG Lee 02/19/2011 - Added limits to constrain dekta_mf to +-10. 
; KG Lee 04/19/2011 - Changed forest range to 1041-1185A. Removed
;                     rebinning.  
; KG Lee 05/25/2011 - Now does only linear fit if just 1 LyaF flux bin is
;                     covered by spectrum. Assume CALC_MEANFLUX
;                     computes nbins=3 bins.
; KG Lee 01/27/2012 - Added IVAR keyword
; KG Lee 02/12/2012 - Fitting is now carried out per-pixel with the
;                     intrinsic forest variance for weighting supplied
;                     by FORESTVAR_FUNC. Everything is converted to
;                     inverse-variance weighting before fitting.
; KG Lee 07/06/2013 - Added FOREST_RANGE keyword to allow modification
;                     of forest range
;-

common mfblock

ff = ff_in
lambda_r = lambda_r_in
if keyword_set(ivar_in) then ivar=ivar_in else sigma = sigma_in

; First fit the input spectrum, then generate the forest region using
; the best-fit parameters. If PCAPARAMS is input, then just recreate
; it 
if not keyword_set(pcaparams) then begin
   fitpca = PCA_CHISQ(ff, lambda_r, sigma,ivar=ivar, $
                      /quiet, afluxerr=afluxerr, $
                      dr7eigen=dr7eigen, m_fit=m_fit)
    ;print, fitpca
endif else begin
    fitpca = pcaparams
endelse 

pcaflux= PCA_FUNC(lambda_r, fitpca, dr7eigen=dr7eigen)
cz = fitpca[0]
;zqso = cz * (1. + zqso) -1.

; Lower and upper wavelength limit of Lya forest
if keyword_set(forest_range) then begin
   for_low = (forest_range[0])[0]
   for_hi  = (forest_range[1])[0]
endif else begin
   for_low = 1041.
   for_hi = 1185.
endelse 

; This is to optionally return the best-fit paramaters to the calling
; procedure 
pcaparams=fitpca  

; Everything will be carried out through inverse-variances from here
; onwards. 
if not keyword_set(ivar_in) then ivar = $
   (sigma NE 0.) / (sigma^2 + (sigma EQ 0)) else ivar=ivar_in

; Normalize flux and errors. Here we don't correct for cz because the
; parameters output from PCA_CHISQ operate on spectra normalized at
; 1280A in the restframe of the original pipeline redshift 
normpix1 = where(lambda_r GE 1275.  AND lambda_r LE 1285.)
if normpix1 NE [-1] AND total(ivar[normpix1]) NE 0 then $
  normfac = AVG(ff[normpix1]*ivar[normpix1])/avg(ivar[normpix1]) $
else begin
    normpix2 = where(lambda_r GE 1450. AND lambda_r LE 1470.)
    normfac = AVG(ff[normpix2]*ivar[normpix2])/avg(ivar[normpix2])
endelse

ivar = ivar * normfac^2 
ff = ff / normfac
pcaflux = pcaflux / normfac

; Now divide the forest region by the fitted PCA continuum
forestrange = where(lambda_r GE for_low AND lambda_r LE for_hi AND $
                    lambda_r*(1.+zqso) GT wavemin)
lamb_forest = lambda_r[forestrange]
z_for = (lamb_forest/1216.7)*(1.+zqso) -1.  ; Redshift at each pixel
fforest = ff[forestrange]/pcaflux[forestrange]
ivarforest = ivar[forestrange] * pcaflux[forestrange]^2

; Estimate weights for each pixel
var_F = forestvar_func(z_for) * (exp(-taueff_evo(z_for)))^2
var_noise = (ivarforest NE 0) / (ivarforest + (ivarforest EQ 0))
var_total = var_F + var_noise
weights_forest = (var_total NE 0) / (var_total + (var_total EQ 0))

; But need to make sure that masked pixels remain masked...
maskedpix = where(ivarforest EQ 0)
if maskedpix NE [-1] then weights_forest[maskedpix] = 0
;-------------------------------------------------------------

pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},2)
delta_mf = [0.,0.]
; Note: we are using errors here instead of inverse variances because
; the errors have been directly estimated from bootstrap in
; CALC_MEANFLUX. Inverse variances, if set, have already been used to
; veto pixels
delta_mf = mpfitfun('MFLUX_TAUEVO_FUNC', lamb_forest, $
                    fforest, replicate(0., n_elements(z_for)), $
                    delta_mf, bestnorm=chisq, $
                    weights=weights_forest,dof=dof,/quiet, $
                    parinfo=pi)
end
