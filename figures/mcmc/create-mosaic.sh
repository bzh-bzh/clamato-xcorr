#!/bin/bash
surveys=(3dhst clamato mosdef vuds zDeep)
montage -mode concatenate -tile 3x0 -background white -geometry -15-0 corner_3dhst.png corner_clamato.png corner_mosdef.png corner_vuds.png corner_zDeep.png corner_mosaic.png
montage -mode concatenate -tile 3x0 -background white bestmodel_3dhst.png bestmodel_clamato.png bestmodel_mosdef.png bestmodel_vuds.png bestmodel_zDeep.png bestmodel_mosaic.png

mkdir -p combined
for survey in "${surveys[@]}" ; do
    montage -mode concatenate -gravity Center -tile 2x1 -background white corner_$survey.png bestmodel_$survey.png combined/combined_$survey.png
done
