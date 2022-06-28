#!/bin/bash
surveys=(3dhst clamato mosdef vuds zDeep)
montage -mode concatenate -tile 3x0 -background white -geometry -15-0 titledcorner_3dhst.png titledcorner_clamato.png titledcorner_mosdef.png titledcorner_vuds.png titledcorner_zDeep.png corner_mosaic.png
montage -mode concatenate -tile 3x0 -background white titledbestmodel_3dhst.png titledbestmodel_clamato.png titledbestmodel_mosdef.png titledbestmodel_vuds.png titledbestmodel_zDeep.png bestmodel_colorbar.png bestmodel_mosaic.png

mkdir -p combined
for survey in "${surveys[@]}" ; do
    montage -mode concatenate -gravity Center -tile 0x1 -background white titledcorner_$survey.png bestmodel_$survey.png bestmodel_colorbar.png combined/combined_$survey.png
done
