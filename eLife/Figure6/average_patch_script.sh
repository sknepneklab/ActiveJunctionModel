#!/bin/bash

topdir_in='/home/sh18581/Documents/AJM_Quarantine/Active_patch/'
topdir_out='/media/sh18581/BristolBackup/AJM_data/Patch/'

inputfile='random_40x40_patch.json'


# Values for Panel A in Figure 7
seeds='1'

betaval='0.5'
fpullval='0.2'


for beta in ${betaval}
do
	for fpull in ${fpullval}
	do
		for s in ${seeds}
		do

            echo ${s}
			confdir=${topdir_out}/beta${beta}/fpull${fpull}/seed${s}/
			mkdir -p ${confdir}
			#infile=${topdir_in}${inputfile}
			#cp ${infile} ${confdir}
			cp random_set_params.py ${confdir}
			cp make_random_passive.py ${confdir}
			cp ../AJM/analysis/runscripts/Analyze_AJM_patch.py ${confdir}
			cd ${confdir}
						
			echo python3 make_random_passive.py
			#python3 make_random_passive.py
			infile=${topdir_in}${inputfile}
						
			echo python3 random_set_params.py --beta ${beta} --forcex ${fpull} --seed ${s} --graner
			python3 random_set_params.py --beta ${beta} --forcex ${fpull} --seed ${s} --graner
			
			echo python3 Analyze_AJM_patch.py
# 			python3 Analyze_AJM_patch.py --filemax 400 --start 300
 			python3 Analyze_AJM_patch.py
			
			cp tensor_timeseries.p ${topdir_in}/picklefiles/beta${beta}/fpull${fpull}/tensor_timeseries_${s}.p

			
			cd ${topdir_in}
		done
	done
done
