; This version has the LBG normalization set at 1245A instead of 1350A!
;+
;
; Estimate Ly-a forest SNR by extrapolating a power law from >1216A. We mask
;absorption lines. Incorporates <F>(z).
;
; For quasars, use PCA continuum estimate.
;
; snr = cnr_forest(wave, flux_sig, zem=zem, /qso)
;
;-


function PLAWFUNC, x, p

  RETURN, p[0] * (x/1245.)^p[1]

end 

function cnr_forest, wave, flux, sig, zem=zem, qso=qso, plot=plot, $
                     cont_out=cont, snr_range=snr_range, dr7eigen=dr7eigen, $
                     minpix =minpix

  ;; Minimum number of pixels within selected range to evaluate
  ;; S/N. Function will return -9. if less than this number of pixels
  ;; are selected
  if not keyword_set(minpix) then minpix= 20.
  
  waverest = wave/ (1.+zem)

  ;; qsos and gals should have different wavelength ranges for
  ;; estimating SNR
  if keyword_set(qso) then $
     wrange = where(waverest GE 1000 AND waverest LT 1600) $
  else wrange = where(waverest GE 1250 AND waverest LT 1500)

  ivar = (sig NE 0) / (sig^2 + (sig EQ 0))

  ;; Mask absorption lines in the case of galaxies
  wmask_min = [1292., 1327., 1384.]
  wmask_max = [1312., 1341., 1411.]
  if not keyword_set(qso) then begin
     maskpix = where( (waverest GE wmask_min[0] AND waverest LT wmask_max[0]) OR $
                      (waverest GE wmask_min[1] AND waverest LT wmask_max[1]) OR $
                      (waverest GE wmask_min[2] AND waverest LT wmask_max[2]) )
     ivar[maskpix] = 0.
  endif

  if not keyword_set(qso) then begin
     fnorm = median(flux[where(waverest GE 1240. AND waverest LT 1250.)])
     slope = -1.
     
     plawpar = mpfitfun('PLAWFUNC',waverest[wrange], flux[wrange], $
                        sig[wrange], [fnorm,slope], weights=ivar[wrange],/quiet)
     cont = plawfunc(waverest, plawpar)
  endif else begin
                                ;stop
     cont = fltarr(n_elements(wave))
     conttmp = pcaspec_func(waverest[wrange], flux[wrange], sig[wrange],ivar=ivar[wrange], $
                         /pca_only,zqso=zem, dr7eigen=dr7eigen)
     cont[wrange] = conttmp
  endelse

     
  if keyword_set(plot) then begin
     set_plot, 'x'
     plot, wave, smooth(flux,3), xran=[3800, 5100], xsty=1
     oplot, wave, sig, color=djs_icolor('red')
     oplot, wave, cont, color=djs_icolor('green')
  endif

  if not keyword_set(snr_range) then $
     snr_range = (1.+zem) * [1041., 1195.]
     
  snrcut = where(wave GE snr_range[0] AND wave LT snr_range[1],npixcut)

  if npixcut GE minpix then begin
     cnr = median(cont[snrcut]* $
                  exp(-taueff_evo(wave[snrcut]/1216.-1.))*(sig[snrcut] GT 0)/ $
                  (sig[snrcut] + (sig[snrcut] EQ 0)) )
  endif else cnr = -9.
     
return, cnr
end
