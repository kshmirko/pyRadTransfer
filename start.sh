#!/usr/bin/env zsh


# 15 km atmosphere layer
NLAYS=15
WL=0.75

TAUA=(0.0 0.1 0.3 0.5)
SZA=(10 30 60 80)
GALBEDO=0.0

for sza_i in $SZA; 
do
    if [ -d sza_$sza_i ];
    then
        rm -rf sza_$sza_i
    fi;

    mkdir sza_$sza_i

    for taua_i in $TAUA; 
    do
    	echo $sza_i $taua_i


        ./rtmain.py prepare_files --layfile=atmoslay.lay --hpbl=3 --taua=$taua_i --wl=$WL --nlays $NLAYS
    	./rtmain.py ssrt --layfile=atmoslay.lay --sza=$sza_i
    	./rtmain.py rt3  --layfile=atmoslay.lay --sza=$sza_i --galbedo=$GALBEDO --wavelen=$WL
    	./rtmain.py vizualize --ssrt ssrt.npz --rt3 rt3.npz --prepare prepare.npz --savef sza_$sza_i/comparison_t$taua_i_g$GALBEDO.png

    done;

done;