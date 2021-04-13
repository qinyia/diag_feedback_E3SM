
module load nco

regrid=False

run_id1=$1 #cori-haswell.20190513.F2010C5-CMIP6-LR.ne30_oECv3
run_id2=$2 #cori-haswell.20190513.F2010C5-CMIP6-LR.plus4K.ne30_oECv3
rgr_map=$3 #~zender/data/maps/map_ne30np4_to_cmip6_72x144_aave.20181001.nc
int_year=$4
end_year=$5
datadir_in1=$6
datadir_in2=$7
outdir_out=$8

run_id=(${run_id1} ${run_id2})
echo ${run_id}
exp_id=(FC5 FC5_4K)

ncase=2 # two cases: one is control simulation, the other is plus4K simulation
echo $ncase

#rgr_map=~zender/data/maps/map_ne30np4_to_cmip6_72x144_aave.20181001.nc

#vars=(FISCCP1_COSP,FSDSC,FSNSC,TREFHT,FSNT,FSNTC,FLUT,FLUTC,TS,T,CLOUD,CLDLIQ,CLDICE,FSDS,FSNS,Q)
#var_list=(FISCCP1_COSP FSDSC FSNSC TREFHT FSNT FSNTC FLUT FLUTC TS T CLOUD CLDLIQ CLDICE FSDS FSNS Q)
#var_new_list=(clisccp rsdscs rsnsc tas rsnt rsntcs rlut rlutcs ts ta CLOUD CLDLIQ CLDICE rsds rsns hus)

#vars=(FISCCP1_COSP,FSDSC,FSNSC,TREFHT,FSNT,FSNTC,FLUT,FLUTC,TS,FSDS,FSNS,TGCLDLWP,TGCLDIWP,T,Q)
#var_list=(FISCCP1_COSP FSDSC FSNSC TREFHT FSNT FSNTC FLUT FLUTC TS FSDS FSNS TGCLDLWP TGCLDIWP T Q)
#var_new_list=(clisccp rsdscs rsnsc tas rsnt rsntcs rlut rlutcs ts rsds rsns TGCLDLWP TGCLDIWP ta hus)

#vars=(FISCCP1_COSP,FSDSC,FSNSC,TREFHT,FSNT,FSNTC,FLUT,FLUTC,TS,T,CLOUD,CLDLIQ,CLDICE,FSDS,FSNS,Q,SOLIN,FSUTOA,FSUTOAC,PSL,PS)
#var_list=(FISCCP1_COSP FSDSC FSNSC TREFHT FSNT FSNTC FLUT FLUTC TS T CLOUD CLDLIQ CLDICE FSDS FSNS Q SOLIN FSUTOA FSUTOAC PSL PS)
#var_new_list=(clisccp rsdscs rsnsc tas rsnt rsntcs rlut rlutcs ts ta CLOUD CLDLIQ CLDICE rsds rsns hus rsdt FSUTOA FSUTOAC psl ps)

vars=(FISCCP1_COSP,FSDSC,FSNSC,TREFHT,FSNT,FSNTC,FLUT,FLUTC,TS,T,FSDS,FSNS,Q,SOLIN,FSUTOA,FSUTOAC,PSL,PS,OMEGA)
var_list=(FISCCP1_COSP FSDSC FSNSC TREFHT FSNT FSNTC FLUT FLUTC TS T FSDS FSNS Q SOLIN FSUTOA FSUTOAC PSL OMEGA)
var_new_list=(clisccp rsdscs rsnsc tas rsnt rsntcs rlut rlutcs ts ta rsds rsns hus rsdt FSUTOA FSUTOAC psl ps wap)


nvar=${#var_list[@]}
echo $nvar

int_mon=1
end_mon=12

int_year_4d=`printf %04d $int_year`
end_year_4d=`printf %04d $end_year`
int_mon_2d=`printf %02d $int_mon`
end_mon_2d=`printf %02d $end_mon`

for ii in `seq 0 $[ncase-1]`
do
#	datadir_in=/global/cscratch1/sd/qinyi/E3SM_predata/${run_id[ii]}/archive/atm/hist/
#	outdir=/global/cscratch1/sd/qinyi/E3SM_middata/${run_id[ii]}/

    if [ "${run_id[ii]}" == "20200428.DECKv1b_amip1-CFMIP.ne30_oEC.cori-knl-L" ] ; then
        datadir=${datadir_in1}/${run_id[ii]}/${datadir_in2}
    else
        datadir=${datadir_in1}/${run_id[ii]}/${datadir_in2}
    fi

    echo $datadir

    outdir=${outdir_out}/${run_id[ii]}/

	if [ ! -d "${outdir}" ] ; then
		mkdir -p ${outdir}
    #else
	#	# if the directory is not empty, pls don't continue the following codes
	#	if [ -z "$(ls -A ${outdir})" ]; then
	#		echo "outdir is Empy"
	#	else
	#		echo "outdir Not Empty"
	#		continue
	#	fi
	fi
	
	if [ "$regrid" == "True" ] ; then
		# use ncremap to regrid monthly data to lat-lon grid by myself, rather than using the output from e3sm_diags
		# because e3sm_diags always needs the December data from the previous year. This leads to the waste of the previous year's data.
		all_file_list=''
		for iyr in `seq ${int_year} ${end_year}`
		do
			for imon in `seq ${int_mon} ${end_mon}`
			do
				iyr_4d=`printf %04d $iyr`
				imon_2d=`printf %02d $imon`
				echo "yr=" ${iyr_4d},"month=" ${imon_2d}
				file_tmp=`ls $datadir/${run_id[ii]}.cam.h0.${iyr_4d}-${imon_2d}.nc`
				out_file=${run_id[ii]}.cam.h0.${iyr_4d}-${imon_2d}_regrid.nc
		
				if [ ! -f "$outdir/${out_file}" ] ; then
					echo "regrid file to", $out_file
					ncremap -v ${vars} -m ${rgr_map} ${file_tmp} $outdir/${out_file}
                else
                    echo "The output file is ready:" $out_file
				fi

				all_file_list="${all_file_list} ${out_file}"

			done
		done
	fi

	# use 'ncrcat' to get monthly-series variables
    cd $outdir
    checkfile=rsuscs_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc

    tmpfile=${run_id[ii]}.cam.h0.${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc

    if [ ! -f "$outdir/$checkfile" ] ; then

        echo "ncrcat all necessary monthly files..."
		ncrcat -O ${all_file_list} $outdir/${tmpfile}

	    in_file=$outdir/${run_id[ii]}.cam.h0.${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc
	    cd $outdir
	    
	    for ivar in `seq 0 $[$nvar-1]`
	    do
	    	echo ${var_new_list[ivar]}
	    	out_file=${var_new_list[ivar]}_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc
        	# extract each variable from the input file
        	ncks -O -v ${var_list[ivar]} ${in_file} ${out_file}
        	# change variable name
            if [ "${var_list[ivar]}" != "${var_new_list[ivar]}" ] ; then
        		ncrename -v ${var_list[ivar]},${var_new_list[ivar]} ${out_file}
            fi
	    done
	    
	    # get surface up SW at clear-sky
        ncks -O -v FSDSC,FSNSC,FSDS,FSNS,SOLIN,FSNT,FSNTC ${in_file} tmp.nc
        ncap2 -O -s 'rsuscs=FSDSC-FSNSC' tmp.nc rsuscs_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc
        ncks -O -v rsuscs rsuscs_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc rsuscs_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc

        # get surface up SW at all-sky
	    ncap2 -O -s 'rsus=FSDS-FSNS' tmp.nc rsus_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc
        ncks -O -v rsus rsus_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc rsus_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc

        # get TOA up SW at all-sky
        ncap2 -O -s 'rsut=SOLIN-FSNT' tmp.nc rsut_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc
        ncks -O -v rsut rsut_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc rsut_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc

        # get TOA up SW at clear-sky
        ncap2 -O -s 'rsutcs=SOLIN-FSNTC' tmp.nc rsutcs_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc
        ncks -O -v rsutcs rsutcs_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc rsutcs_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc

        rm tmp.nc
        rm $outdir/${run_id[ii]}.cam.h0.${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc

    else
        echo "All necessary output files are ready."
        echo `ls *${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc`
    fi

done 

echo "The final outdir =" $outdir
echo "Well Done!"

