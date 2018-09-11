#!/usr/bin/env zsh



NLAYS=5
WL=0.870

TAUA=0.4
#SZA=(10  30  50  70)
SZA=(25)
GALBEDO=(0.0 0.3)
#GALBEDO=(0.0 0.3 0.6  0.9)

for sza_i in $SZA; 
do
    if [ ! -d sza_$sza_i ];
    then
        mkdir sza_$sza_i
    fi;

    for galbedo_i in $GALBEDO; 
    do
    	echo $sza_i $galbedo_i
        if [ -f sza_$sza_i/conparizon_g$galbedo_i\_t=$TAUA.png ];
        then
            rm sza_$sza_i/conparizon_g$galbedo_i\_t=$TAUA.png
        fi;


        ./rtmain.py prepare_files --layfile=atmoslay.lay --hpbl=3 --taua=$TAUA --wl=$WL --nlays $NLAYS
    	./rtmain.py ssrt --layfile=atmoslay.lay --sza=$sza_i
    	./rtmain.py rt3  --layfile=atmoslay.lay --sza=$sza_i --galbedo=$galbedo_i --wavelen=$WL
    	./rtmain.py vizualize --ssrt ssrt.npz --rt3 rt3.npz --prepare prepare.npz --savef sza_$sza_i/conparizon_g$galbedo_i\_t=$TAUA

    done;

done;