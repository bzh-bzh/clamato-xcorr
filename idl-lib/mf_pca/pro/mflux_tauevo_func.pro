;------------------------------------------------------------------------

function MFLUX_TAUEVO_FUNC, x, p
; This function is the mean flux evolution * exp(delta*(lambda/1280
; -1)), meant for use with MPFITFUNC to correct the fitted Lya forest
; continuum
;
; Fit to preliminary Lya transmission spectrum with some binning
;
; Abscissa parameter x is the restframe wavelength, and free parameter
; p is delta. This function NEEDS the quasar redshift zqso to be in
; the common block 

common mfblock

zfor = (x/1216.)* (1. + zqso) - 1.

tau = taueff_evo(zfor)

fmean = exp(-tau)

mfluxtauevo = fmean * MFLUXCORR(x,p, LAMB_PIV=lamb_piv)

return, mfluxtauevo
end
