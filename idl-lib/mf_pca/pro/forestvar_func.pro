function FORESTVAR_FUNC, z_in
; Return intrinsic variance of LyaF variance for weighting. This
; estimate is roughly from McDonald et al 2006

return, 0.065 * (( 1.+z_in)/(1.+2.25))^3.8
end
