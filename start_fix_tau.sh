#!/usr/bin/env zsh



NLAYS=15
WL=0.75

TAUA=0.5
SZA=(10  30  50  70)
GALBEDO=(0.0 0.3 0.6  0.9)

for sza_i in $SZA; 
do
    #if [ -d sza_$sza_i ];
    #then
    #    rm -rf sza_$sza_i
    #fi;

    #mkdir sza_$sza_i

    for galbedo_i in $GALBEDO; 
    do
    	echo $sza_i $galbedo_i


        ./rtmain.py prepare_files --layfile=atmoslay.lay --hpbl=3 --taua=$TAUA --wl=$WL --nlays $NLAYS
    	./rtmain.py ssrt --layfile=atmoslay.lay --sza=$sza_i
    	./rtmain.py rt3  --layfile=atmoslay.lay --sza=$sza_i --galbedo=$galbedo_i --wavelen=$WL
    	./rtmain.py vizualize --ssrt ssrt.npz --rt3 rt3.npz --prepare prepare.npz --savef sza_$sza_i/conparizon_g$galbedo_i\_t=$TAUA.png

    done;

done;