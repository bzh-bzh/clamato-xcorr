import papermill
import os
os.chdir('..')

import constants


if __name__ == '__main__':
    papermill.execute_notebook('mocksurvey_lyagen.ipynb',
                              os.path.join(constants.BIAS_DIR_BASE, 'papermill', 'mocksurvey_lyagen_out.ipynb'),
                              kernel_name='clamato-xcorr',
                              parameters={
                                  'outdir': os.path.join(constants.BIAS_DIR_BASE, 'xcorr', 'mock', 'skewers'),
                                  'nmock_x': 1,
                                  'nmock_y': 1,
                                  'nabs': 1,
                                  'width_x': 250,
                                  'width_y': 250
                              })