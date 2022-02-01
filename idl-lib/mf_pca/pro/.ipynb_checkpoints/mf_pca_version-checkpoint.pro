;+
; NAME:
;   mf_pca_version()
; PURPOSE:
;   Return the version name for the mf_pca product
; CALLING SEQUENCE:
;   vers = mf_pca_version()
; OUTPUTS:
;   vers       - Version name for the product mf_pca
; COMMENTS:
;   Requires that the MF_PCA_DIR environment variable be set
; VERSION:
;   $Id: mf_pca_version.pro 132130 2012-03-30 08:13:20Z kheegan $
;-
;------------------------------------------------------------------------------
FUNCTION mf_pca_version
    RETURN, FILE_BASENAME(GETENV('MF_PCA_DIR'))
END
;------------------------------------------------------------------------------
