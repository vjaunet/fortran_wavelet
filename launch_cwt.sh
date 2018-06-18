#!/bin/bash
#
#
#=====================================================

. /home/jaunetv/DATA_BASE/UTILITIES/nrwait.sh

# trap ctrl-c and call ctrl_c()
trap ctrl_c INT

function ctrl_c() {
    killall cwvttr
    echo "** Trapped CTRL-C"
    exit 0
}


for imic in 27
do
    for ifile in ../../DATA/PLAQUE_R050_4D/MACH_080/M*_T0.bin \
		 ../../DATA/PLAQUE_R050_4D/MACH_078/M*_T0.bin \
		 ../../DATA/PLAQUE_R050_2D/MACH_060/M*_T0.bin \
		 ../../DATA/PLAQUE_R050_3D/MACH_080/M*_T0.bin
    do

	echo "Processing "$ifile

	#build the output file name and folder
	ofile=$(basename $ifile)
	ofolder=$(dirname $ifile)/CWT
	mkdir -p $ofolder
	ofile=${ofolder}/${ofile%.*}_mic${imic}_bump_cwt.dat

	#launch the computation
	(./cwvttr -i $ifile -o $ofile \
		  -imic $imic \
		  -type 'Bump' -nsamp 4000000 ) &
	#-type 'Morl' -nsamp 1000000 ) &
	#-type 'Paul' -order 4 -nsamp 1000000 ) &



	nrwait 4

    done
done
wait
