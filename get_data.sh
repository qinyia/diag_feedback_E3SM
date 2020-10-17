
module load nco

regrid=True

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

vars=(FISCCP1_COSP,FSDSC,FSNSC,TREFHT,FSNT,FSNTC,FLUT,FLUTC,TS,FSDS,FSNS,TGCLDLWP,TGCLDIWP,T,Q)
var_list=(FISCCP1_COSP FSDSC FSNSC TREFHT FSNT FSNTC FLUT FLUTC TS FSDS FSNS TGCLDLWP TGCLDIWP T Q)
var_new_list=(clisccp rsdscs rsnsc tas rsnt rsntcs rlut rlutcs ts rsds rsns TGCLDLWP TGCLDIWP ta hus)

nvar=${#var_list[@]}
echo $nvar

#int_year=1
#end_year=5
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

    datadir=${datadir_in1}/${run_id[ii]}/${datadir_in2}/
    outdir=${outdir_out}/${run_id[ii]}/

	if [ ! -d "${outdir}" ] ; then
		mkdir -p ${outdir}
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
		#		all_file_list="${all_file_list} ${list_tmp}"
				out_file=${run_id[ii]}.cam.h0.${iyr_4d}-${imon_2d}_regrid.nc
		
				#if [ ! -f "$outdir/${out_file}" ] ; then
					echo "regrid file to", $out_file
					ncremap -v ${vars} -m ${rgr_map} ${file_tmp} $outdir/${out_file}
				#fi
			done
		done

		# use 'ncrcat' to get monthly-series variables
		ncrcat -O $outdir/*.cam.h0.????-??_regrid.nc $outdir/${run_id[ii]}.cam.h0.${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc
	fi
	
	echo $outdir
	in_file=$outdir/${run_id[ii]}.cam.h0.${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc
	
	cd $outdir
	
	for ivar in `seq 0 $[$nvar-1]`
	do
		echo ${var_new_list[ivar]}
		out_file=${var_new_list[ivar]}_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc
		# extract each variable from the input file
		ncks -O -v ${var_list[ivar]} ${in_file} ${out_file}
		# change variable name
		ncrename -v ${var_list[ivar]},${var_new_list[ivar]} ${out_file}
	done
	
	# get surface up SW at clear-sky
	ncap2 -O -s 'rsuscs=FSDSC-FSNSC' ${in_file} rsuscs_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc
    ncks -O -v rsuscs rsuscs_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc rsuscs_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc

    # get surface up SW at all-sky
	ncap2 -O -s 'rsus=FSDS-FSNS' ${in_file} rsus_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc
    ncks -O -v rsus rsus_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc rsus_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc

	echo $outdir
done 



