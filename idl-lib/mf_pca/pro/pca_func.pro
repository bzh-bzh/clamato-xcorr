;---------------------------------------------------------------------
function PCA_FUNC, x, p, dr7eigen=dr7eigen
; For a given value of wavelength, return normalized flux for the
; corresponding set of PCA weights. To be used in conjunction with
; MPFITFUN 

; PCABLOCK needs to be initialized beforehand.
; 
; Optional Keywords:
;   DR7EIGEN  - If set, use Paris+11 DR7 eigenspectra.
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

pcaspec = pcaspec / avg(pcaspec[where(lambxi GE 1275. AND $
                                      lambxi LE 1285.)])

; Now interpolate for value of normalized flux corresponding to input
; wavelength * cz
ff = fnorm * interpol(pcaspec, lambxi_tmp, cz*x) * (cz*x/1280.)^alpha

return, ff
end
