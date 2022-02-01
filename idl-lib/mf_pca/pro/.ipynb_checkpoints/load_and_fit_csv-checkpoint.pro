harris_lya = READ_CSV('../../notebooks/harris_lya.csv')
mfpca_spec = PCASPEC_FUNC(harris_lya.FIELD1, harris_lya.FIELD2, harris_lya.FIELD3, zqso=2.5, pca_only=0, dr7eigen=1)
WRITE_CSV, 'harris_mfpca.csv', mfpca_spec
