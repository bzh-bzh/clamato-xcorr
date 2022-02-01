function PCASPEC_FUNC, lambda_r_in, ff_in, sigma_in, ivar=ivar_in, $
                       pca_only=pca_only, pcacont=pcacont,$
                       dr7eigen=dr7eigen, m_fit=m_fit, $
                       zqso=zqso1,  mfflag=mfflag1, wavemin=wavemin1, $
                       delta_mf=delta_mf,pcaparams=pcaparams, $
                       afluxerr=afluxerr, dla_in=dla_in , ff_dla=ff_dla, $
                       forest_range=forest_range
;
;+
; NAME:
;   PCASPEC_FUNC 
;
;
; PURPOSE:
;
;   Returns best-fit continuum corresponding to the input spectrum,
;   using mean-flux regulated PCA technique (Lee et al 2011,
;   arXiv:1108:6080). This is just a wrapper for MF_PCA.PRO and
;   PCA_CHISQ.PRO  
;
;   The code can be set to mask DLAs and correct its damping
;   wings. Note that since the errors are also rescaled by the damping
;   corrections.
;
;
; CATEGORY:
;   Wrapper function
;
;
; CALLING SEQUENCE:
;   mfpca_spec = pcaspec_func(wave, flux, [sigma, ivar=ivar_in,
;   /pca_only, /dr7eigen, pcacont=pcacont,
;   zqso=zqso1, mfflag=mfflag1, wavemin=wavemin1, m_fit=m_fit,
;   afluxerr=afluxerr,pcaparams=pcaparams, delta_mf=delta_mf,
;   dla_in=dla_in, forest_range=forest_range]
;
;
; INPUTS:
;   WAVE        - Wavelength grid, *IN QUASAR RESTFRAME*
;   FLUX        - Quasar spectrum
;
; OPTIONAL INPUTS:
;   SIGMA_IN     - Errors on spectrum. Either this or IVAR must be
;                  input. SIGMA=0 is equivalent to IVAR=0. This gets
;                  converted into IVAR in the code. It is recommended
;                  to stick to IVAR since masks can more easily be
;                  incorporated there.  
;
;
; KEYWORD PARAMETERS:
;   IVAR         - Inverse variances of the spectrum, instead of
;                  errors. Either this or SIGMA_IN must be input. 
;   PCA_ONLY     - If set, will do continuum fitting using only PCA
;                  (without mean-flux regulation)
;   PCACONT      - Returns pure PCA continuum, even when mean-flux
;                  regulation is being used.
;   DR7EIGEN     - If set, will use Paris et al 2011 eigenspectra
;                  templates for fitting. Otherwise, defaults to using
;                  Suzuki et al 2005 templates
;   M_FIT        - Number of principal components with which to do the
;                  fit. Defaults to 8.
;   ZQSO         - Quasar redshift. Alternatively, the calling program
;                  can define the MFBLOCK common block with this
;                  parameter
;   MFFLAG       - Integer specifying which mean-flux measurement to
;                  use. Defaults to 2, which is Faucher-Giguere 2008's
;                  power-law with metals. Alternatively, the calling
;                  program can define the MFBLOCK common block with
;                  this parameter 
;   WAVEMIN      - Minimum *observed* wavelength value - anything
;                  below this value will be excised. Defaults to
;                  3700. Alternatively, the calling
;                  program can define the MFBLOCK common block with
;                  this parameter 
;   DLA_IN       - List of the absorbers redshift and log10(N_HI)of
;                  known DLAs wihtin the spectrum. The 
;                  central EW of the DLA will be masked, while a
;                  correction will be applied to the wings [2, NDLA]
;   PCAPARAMS    - Can be set on input to supply PCA parameters if
;                  only mean-flux regulation is required. Will be
;                  ignored if PCA_ONLY keyword is set.
;   FOREST_RANGE - 2-element vector defining the lower- and
;                  upper-limit of restframe wavelengths of the Lya
;                  forest region
; 
;
; OUTPUTS:
;   PCASPEC      - Best-fit continuum. Corresponds to input wavelength
;                  grid, but only pixels in rest wavelengths 1030-1600
;                  are evaluated - all else are zeroed
;
;
; OPTIONAL OUTPUTS:
;   PCAPARAMS    - Best-fit PCA parameters (inclusive of things like
;                  normalization and redshift correction) for the
;                  fit. See PCA_CHISQ.PRO for full list of parameters
;   DELTA_MF     - 2-parameter vector with the linear and quadratic
;                  mean-flux regulation parameters 
;   PCACONT      - Returns pure PCA continuum without
;                  MF-regulation. This is returned even when
;                  MF is being enforced. Somewhat degenerate with
;                  PCA_ONLY keyword
;   AFLUXERR     - Returns absolute flux error in the range
;                  1230-1600A, i.e. ABS(cont/flux - 1)
;   FF_DLA       - If DLAs are specified, then return the flux vector
;                  after EW has been masked (FF_DLA = 0) and wings
;                  corrected.
;
;
; COMMON BLOCKS:
;   MFBLOCK      - Contains ZQSO, MFFLAG and WAVEMIN. Can also be
;                  called from the calling function.
;
;
; PROCEDURES/FILES:
;   pca_chisq, mf_pca, readpca_hst, readpca_sdss, mfluxcorr,
;   xi_hst.txt, xi_sdss.txt
;
; COMMENTS:
;   Need to set $MF_PCA_DIR environment variable. $MF_PCA_DIR/pro to
;   contain routines, and $MF_PCA_DIR/data contains eigenspectra templates. 
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
; Khee-Gan Lee 9/1/2010 - Created to clean up misc lines of code
; Khee-Gan Lee 10/21/2011 - Cleaned up to work with DR9 versions of
;                           code. Should now function as a blackbox
;                           wrapper for MF-PCA
; Khee-Gan Lee 01/27/2012 - Added option to input inverse variances
;                           instead of errors
; Khee-Gan Lee 02/14/2012 - Added option to mask/correct known
;                           DLAs 
; Khee-Gan Lee 04/05/2012 - Use only inverse variances
;                           internally. Input errors get converted to
;                           IVAR. 
; Khee-Gan Lee 07/06/2012 - Added FOREST_RANGE keyword
; Khee-Gan Lee 09/15/2015 - Fixed bug where PCA_ONLY fits were fed the
;                           raw ivar_in instead of with various cuts   
;-


  COMMON MFBLOCK, zqso, mfflag, wavemin, lamb_piv


  if keyword_set(zqso1) then zqso=zqso1 else $
        print, 'ERROR: ZQSO needs to be set'

  if not keyword_set(mfflag1) then mfflag = 2 else mfflag = mfflag1
  if not keyword_set(wavemin1) then wavemin = 3600. else wavemin=wavemin1

  if not keyword_set(ivar_in) AND not keyword_set(sigma_in) then begin
     print, 'ERROR: Either flux errors or inverse variances must be input' 
     err_out = 1
  endif 

  if not keyword_set(err_out) then begin
     ff = ff_in
     lambda_r = lambda_r_in
     
     if keyword_set(ivar_in) then ivar = ivar_in else $
        ivar = (sigma_in NE 0)/(sigma_in^2 + (sigma_in EQ 0))

 ; Now, deal with DLAs
     if keyword_set(dla_in) then begin
        if (size(dla_in))[0] EQ 2 then n_dla = (size(dla_in))[2] else $
           n_dla = 1
                                          
        ff_dla = ff
         maskedpix = where(ivar EQ 0) 
        ff_dla[maskedpix] = 0.

        for rr=0L, n_dla -1L do begin
           
           zabs_tmp = (reform(dla_in[0,rr]))[0]
           nhi = (10.^reform(dla_in[1,rr]))[0]

           ;; Mask pixels which are within the equivalent width
           ;; of DLA
           dlamask = dla_mask_kg(lambda_r*(1.+zqso), nhi, zabs_tmp, ew_dla=ew_dla)
           cutpix = where(dlamask LT 1.e-6, complement=dla_notmask)    
           ff_dla[cutpix] = 0.
           
           ivar[cutpix] = 0. 

           exptau = dla_correctpix_kg(lambda_r[dla_notmask]*(1.+zqso), $
                                      nhi, z_abs=zabs_tmp)
           ff[dla_notmask] = ff[dla_notmask] * exptau
           ff_dla[dla_notmask] = ff_dla[dla_notmask] * exptau
           ivar[dla_notmask] = ivar[dla_notmask] / exptau^2
           
        endfor
     endif 

; Deal with case in which input wavelength range includes stuff
; outside 1030-1600A
     waveexcess = where(lambda_r LT 1030. OR lambda_r GT 1600., $
                        complement=wavecut)
     ;waveexcess = where(lambda_r LT 0. OR lambda_r GT 3000., $
                        ;complement=wavecut)
     if waveexcess NE [-1] then begin
        lambda_r_full = lambda_r_in
        pcaspec_full = replicate(0.,n_elements(ff))
        lambda_r = lambda_r[wavecut]
        ff = ff[wavecut]
        ivar = ivar[wavecut]
     endif 

     if keyword_set(forest_range) then $
        for_hi = (forest_range[1])[0] else for_hi=1185.
     
     if not keyword_set(pca_only) then begin
        
        forest_tmp = mf_pca(ff, lambda_r, sigma,ivar=ivar, $
                            delta_mf=delta_mf, $
                            pcaparams=pcaparams, afluxerr=afluxerr, $
                            dr7eigen=dr7eigen, m_fit=m_fit)
        cz = pcaparams[0]
        pcaspec = pca_func(lambda_r, pcaparams, dr7eigen=dr7eigen) 
        pcacont=pcaspec
        pcaspec[where(lambda_r LE for_hi)] = $
           pcaspec[where(lambda_r LE for_hi)] * $
           mfluxcorr(lambda_r[where(lambda_r LE for_hi)], delta_mf)
     endif else begin 
        pcaparams = pca_chisq(ff, lambda_r, sigma, ivar=ivar, $
                              afluxerr=afluxerr,/quiet, $
                              dr7eigen=dr7eigen, m_fit=m_fit)
        cz = pcaparams[0]

        pcaspec = pca_func(lambda_r, pcaparams, dr7eigen=dr7eigen)
        pcacont = pcaspec
     endelse

        if waveexcess NE [-1] then begin
        pcaspec_full[wavecut] = pcaspec
        pcaspec = pcaspec_full
        pcacont_tmp = replicate(0., n_elements(pcaspec))
        pcacont_tmp[wavecut] = pcacont
        pcacont = pcacont_tmp
     end

  endif else pcaspec = [-1]

  return, pcaspec
end
