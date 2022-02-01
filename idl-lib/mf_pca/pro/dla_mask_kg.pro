;+
; NAME:
;   dla_mask_kg
;;
; PURPOSE:
;   For an input wavelength vector, mask DLA region based on
;   equivalent width (using input HI column density) 
;
; CALLING SEQUENCE:
;   dla_mask_kg, wave, n_hi, z_abs, [ew_dla=ew, /log_column]
;
;
; INPUTS:
;
;   WAVE     - Wavelength vector in observed frame
;   N_HI     - Neutral hydrogen column density of DLA
;   Z_ABS    - Absorption redshift of DLA
;
; KEYWORD PARAMETERS:
;
;   LOG_COLUMN   - If set, then N_HI will be assumed to be log10(N_HI)
;   EW_DLA       - Outputs EW of the DLA, in restframe
;
; OUTPUTS:
;
;   MASKFLAG   - A vector corresponding to input wavelength vector,
;                where masked pixels are set to 0, and 1 otherwise  
;
; INTERNAL SUPPORT ROUTINES:
;
;   EW_DLA
;

; MODIFICATION HISTORY:
;
;       Thu Jul 28 15:03:45 2011, Khee-Gan Lee
;
;-
FUNCTION EW_DLA, N_HI
; Calculate equivalent width of damped hydrogen Lya
; absorber, using formula from Eq 9.22 of Drain's 'Physics of the ISM' 

e = 4.80320d-10
m_e = 9.109389d-28
c=2.99792d10
gamma_lambda = 7616.d
lambda_cm = 1.21567d-5 
f_lu = 0.4164d

ew = sqrt( (e^2/m_e/c^2) * N_HI * f_lu * lambda_cm * $
           gamma_lambda / c)

return, ew * 1215.67
end


FUNCTION DLA_MASK_KG, WAVE, N_HI, Z_ABS, LOG_COLUMN=LOG_COLUMN, $
                      EW_DLA=ew

if keyword_set(log_column) then nhi_tmp=10^n_hi else $
  nhi_tmp = n_hi

; Translate observed wavelength scale to DLA restframe
wave_dla = wave / (1.+ z_abs)

ew = ew_dla(nhi_tmp)

maskflag = replicate(1., n_elements(wave))

masklist = where(wave_dla GE (1215.67-ew/2.) AND $
                 wave_dla LT (1215.67+ew/2.))

maskflag[masklist] = 0.

RETURN, MASKFLAG
END
