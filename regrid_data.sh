#****************************************************************
#
#    Filename: get_data.sh
#
#    Author: Yi Qin - qin4@llnl.gov
#    Description: only used to regrid all monthly data to FV grids, including all variables.
#    Input: 
#    Output: 
#    Create: 2020-??-??
#    Last Modified: 2021-07-26 22:59:58
#    Aug 11, 2021: only used to regrid all monthly data to FV grids, including all variables.
#****************************************************************

#module load nco

run_id1=$1 #cori-haswell.20190513.F2010C5-CMIP6-LR.ne30_oECv3
run_id2=$2 #cori-haswell.20190513.F2010C5-CMIP6-LR.plus4K.ne30_oECv3
rgr_map=$3 #~zender/data/maps/map_ne30np4_to_cmip6_72x144_aave.20181001.nc
int_year1=$4
end_year1=$5
int_year2=$6
end_year2=$7
datadir_in1=$8
datadir_in2=$9
outdir_out=${10}
comp=${11}

echo ${run_id1}
echo ${run_id2}
echo ${rgr_map}
echo ${int_year1}
echo ${end_year1}
echo ${int_year2}
echo ${end_year2}
echo ${datadir_in1}
echo ${datadir_in2}
echo ${outdir_out}
echo ${comp}

int_year1_4d=`printf %04d $int_year1`
end_year1_4d=`printf %04d $end_year1`

run_id=(${run_id1} ${run_id2})
int_years=(${int_year1} ${int_year2})
end_years=(${end_year1} ${end_year2})

echo ${run_id}
exp_id=(FC5 FC5_4K)

ncase=2 # two cases: one is control simulation, the other is plus4K simulation
echo $ncase

vars=(FISCCP1_COSP,FSDSC,FSNSC,TREFHT,FSNT,FSNTC,FLUT,FLUTC,TS,T,FSDS,FSNS,Q,SOLIN,FSUTOA,FSUTOAC,PSL,PS,T,Q,CLOUD,CLDLIQ,CLDICE,PRECC,PRECL,U,V,OMEGA,TGCLDCWP,TGCLDIWP,TGCLDLWP,CLDLOW,CLDLOW_CAL)
var_list=(FISCCP1_COSP FSDSC FSNSC TREFHT FSNT FSNTC FLUT FLUTC TS T FSDS FSNS Q SOLIN FSUTOA FSUTOAC PSL PS T Q CLOUD CLDLIQ CLDICE PRECC PRECL U V OMEGA TGCLDCWP TGCLDIWP TGCLDLWP CLDLOW CLDLOW_CAL)
var_new_list=(clisccp rsdscs rsnsc tas rsnt rsntcs rlut rlutcs ts ta rsds rsns hus rsdt FSUTOA FSUTOAC psl ps T Q CLOUD CLDLIQ CLDICE PRECC PRECL U V OMEGA TGCLDCWP TGCLDIWP TGCLDLWP CLDLOW CLDLOW_CAL)

nvar=${#var_list[@]}
echo $nvar

int_mon=1
end_mon=12

int_mon_2d=`printf %02d $int_mon`
end_mon_2d=`printf %02d $end_mon`

for ii in `seq 0 $[ncase-1]`
do

    if [ "${run_id[ii]}" == "20200428.DECKv1b_amip1-CFMIP.ne30_oEC.cori-knl-L" ] ; then
        datadir=${datadir_in1}/${run_id[ii]}/${datadir_in2}
    else
        datadir=${datadir_in1}/${run_id[ii]}/${datadir_in2}
    fi

    echo $datadir

    int_year=${int_years[ii]}
    end_year=${end_years[ii]}

    int_year_4d=`printf %04d $int_year`
    end_year_4d=`printf %04d $end_year`

    #<qinyi 2021-08-12 #------------------
    # if the destinated grid is 72x144. save them in the file name with this tag. 
    echo ${rgr_map} | grep "72x144" 
    if [ $? != 0 ] ; then
        outdir=${outdir_out}/${run_id[ii]}/regrid/
    else
        outdir=${outdir_out}/${run_id[ii]}_72x144/regrid/
    fi

	if [ ! -d "${outdir}" ] ; then
		mkdir -p ${outdir}
    fi


    #=============================================================
    # use ncremap to regrid monthly data to lat-lon grid by myself, rather than using the output from e3sm_diags
    # because e3sm_diags always needs the December data from the previous year. This leads to the waste of the previous year's data.
    #=============================================================

	for iyr in `seq ${int_year} ${end_year}`
	do
		for imon in `seq ${int_mon} ${end_mon}`
		do
			iyr_4d=`printf %04d $iyr`
			imon_2d=`printf %02d $imon`
			echo "yr=" ${iyr_4d},"month=" ${imon_2d}
			file_tmp=`ls $datadir/${run_id[ii]}.${comp}.${iyr_4d}-${imon_2d}.nc`
			out_file=${run_id[ii]}.${comp}.${iyr_4d}-${imon_2d}_regrid.nc
	
            #<qinyi 2021-07-26 #------------------
            # check whether the last variable is in the file. If not, still need run the regrid process to generate the new xxx_regrid.nc
            ncks -O -v ${var_list[-1]} $outdir/${out_file} $outdir/tmp.nc >& /dev/null

			if [ $? != 0 -o ! -f "$outdir/${out_file}" ] ; then
            #>qinyi 2021-07-26 #------------------
				echo "regrid file to", $out_file
				ncremap -m ${rgr_map} ${file_tmp} $outdir/${out_file}
            else
                echo "The output file is ready:" $out_file
			fi


		done
	done

done 

echo "Regridding is successfully."
echo "Well Done!"
