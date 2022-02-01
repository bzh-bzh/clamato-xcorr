#!/bin/bash
module load idl/8.5
#export IDL_DIR=$IDL_DIR:/global/u1/b/bzh/quasar-covar/mf_pca/idlutils-v5_5_32/external
export IDLUTILS_DIR=/global/u1/b/bzh/idl-lib/idlutils-v5_5_32
export PATH=$IDLUTILS_DIR/bin:$PATH
export IDL_PATH=$IDL_PATH:$IDLUTILS_DIR/pro
export IDL_PATH=$IDL_PATH:$IDLUTILS_DIR/goddard/pro
export IDL_PATH=$IDL_PATH:/global/u1/b/bzh/idl-lib/mf_pca/pro
export IDL_PATH=$IDL_PATH:/global/u1/b/bzh/idl-lib/astrolib/pro
export IDL_PATH=$IDL_PATH:/global/u1/b/bzh/idl-lib/astrolib/pro/coyote
export IDL_PATH=$IDL_PATH:/global/u1/b/bzh/idl-lib/mpfit
export MF_PCA_DIR=/global/u1/b/bzh/idl-lib/mf_pca
