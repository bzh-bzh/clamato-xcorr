function MFLUXCORR, x, p, LAMB_PIV=lamb_piv
; Correction factor to power-law fit
if not keyword_set(LAMB_PIV) then $
   lamb_piv = 1113.             ; This is the pivot point in the restframe spectrum

return, p[0] + p[1]*(x/lamb_piv - 1.)
end
