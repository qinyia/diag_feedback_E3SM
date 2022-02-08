#****************************************************************
#
#    Filename: get_data.sh
#
#    Author: Yi Qin - qin4@llnl.gov
#    Description: 
#    Input: 
#    Output: 
#    Create: 2020-??-??
#    Last Modified: 2021-07-26 22:59:58
#    Jul 26, 2021: add TGCLDCWP,TGCLDIWP,TGCLDLWP,CLDLOW,CLDLOW_CAL in output
#                  change the condition to decide whether executing regriding process
#****************************************************************


#module load nco

regrid=False

run_id1=$1 #cori-haswell.20190513.F2010C5-CMIP6-LR.ne30_oECv3
run_id2=$2 #cori-haswell.20190513.F2010C5-CMIP6-LR.plus4K.ne30_oECv3
rgr_map=$3 #~zender/data/maps/map_ne30np4_to_cmip6_72x144_aave.20181001.nc
int_year1=$4 # control start year
end_year1=$5 # control end year
int_year2=$6 # warming start year
end_year2=$7 # warming end year
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

int_year2_4d=`printf %04d $int_year2`
end_year2_4d=`printf %04d $end_year2`

run_id=(${run_id1} ${run_id2})
int_years=(${int_year1} ${int_year2})
end_years=(${end_year1} ${end_year2})

echo ${run_id}
exp_id=(FC5 FC5_4K)

ncase=2 # two cases: one is control simulation, the other is plus4K simulation
echo $ncase

vars=(FISCCP1_COSP,FSDSC,FSNSC,TREFHT,FSNT,FSNTC,FLUT,FLUTC,TS,T,FSDS,FSNS,Q,SOLIN,FSUTOA,FSUTOAC,PSL,PS,CLOUD,CLDLIQ,CLDICE,PRECC,PRECL,U,V,OMEGA,TGCLDCWP,TGCLDIWP,TGCLDLWP,CLDLOW,CLDLOW_CAL,Z3)
var_list=(FISCCP1_COSP FSDSC FSNSC TREFHT FSNT FSNTC FLUT FLUTC TS T FSDS FSNS Q SOLIN FSUTOA FSUTOAC PSL PS CLOUD CLDLIQ CLDICE PRECC PRECL U V OMEGA TGCLDCWP TGCLDIWP TGCLDLWP CLDLOW CLDLOW_CAL Z3)
#<qinyi 2021-08-12 #------------------
# to reduce the time consuming for ncrename, I change the output variable names to its default model name, rather than the CMIP conventional name.
#var_new_list=(clisccp rsdscs rsnsc tas rsnt rsntcs rlut rlutcs ts ta rsds rsns hus rsdt FSUTOA FSUTOAC psl ps CLOUD CLDLIQ CLDICE PRECC PRECL U V OMEGA TGCLDCWP TGCLDIWP TGCLDLWP CLDLOW CLDLOW_CAL)
var_new_list=(FISCCP1_COSP rsdscs rsnsc tas rsnt rsntcs rlut rlutcs ts T rsds rsns Q rsdt FSUTOA FSUTOAC psl ps CLOUD CLDLIQ CLDICE PRECC PRECL U V OMEGA TGCLDCWP TGCLDIWP TGCLDLWP CLDLOW CLDLOW_CAL Z3)
#>qinyi 2021-08-12 #------------------


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
        outdir=${outdir_out}/${run_id[ii]}/
    else
        outdir=${outdir_out}/${run_id[ii]}_72x144/
    fi

	if [ ! -d "${outdir}" ] ; then
		mkdir -p ${outdir}
    fi

    #=============================================================
    # get all necessary files for extacting variables        
    #=============================================================
	all_file_list=''
	for iyr in `seq ${int_year} ${end_year}`
	do
		for imon in `seq ${int_mon} ${end_mon}`
		do
			iyr_4d=`printf %04d $iyr`
			imon_2d=`printf %02d $imon`
			echo "yr=" ${iyr_4d},"month=" ${imon_2d}
			out_file=${outdir}/regrid/${run_id[ii]}.${comp}.${iyr_4d}-${imon_2d}_regrid.nc
			#out_file=${outdir}/${run_id[ii]}.${comp}.${iyr_4d}-${imon_2d}_regrid.nc

			all_file_list="${all_file_list} ${out_file}"

		done
	done

    #=============================================================
    # extract variables     
    #=============================================================

    cd $outdir
    
    #<qinyi 2021-07-26 #------------------
    # change checkfile variable from 'rsuscs' to the last variable in var_new_list to ensure all variables are processed.
    checkfile=${var_new_list[-1]}_${exp_id[ii]}_${int_year2_4d}${int_mon_2d}-${end_year2_4d}${end_mon_2d}.nc
    echo "checkfile=" $checkfile
    #>qinyi 2021-07-26 #------------------

    if [ ! -f "$outdir/$checkfile" ] ; then

        in_file=${all_file_list}

	    for ivar in `seq 0 $[$nvar-1]`
	    do
	    	echo ${var_new_list[ivar]}
	    	out_file=${var_new_list[ivar]}_${exp_id[ii]}_${int_year2_4d}${int_mon_2d}-${end_year2_4d}${end_mon_2d}.nc

            if [ ! -f "$outdir/${out_file}" ] ; then

                ncrcat -O -v ${var_list[ivar]} ${in_file} tmp.nc
        	    # change variable name
                if [ "${var_list[ivar]}" != "${var_new_list[ivar]}" ] ; then
                    #<qinyi 2021-08-12 #------------------ CDO rename seems faster than NCO. Change it.
                    # too bad. cdo chname cannot process 5-D variables like clisccp.
        	    	ncrename -v ${var_list[ivar]},${var_new_list[ivar]} tmp.nc ${out_file}
        	    	#cdo chname,${var_list[ivar]},${var_new_list[ivar]} tmp.nc ${out_file}
                    rm tmp.nc
                else
                    mv tmp.nc ${out_file}
                fi
            fi

	    done
	    
        ncrcat -O -v FSDSC,FSNSC,FSDS,FSNS,SOLIN,FSNT,FSNTC ${in_file} tmp.nc

	    # get surface up SW at clear-sky
        out_file=rsuscs_${exp_id[ii]}_${int_year2_4d}${int_mon_2d}-${end_year2_4d}${end_mon_2d}.nc
        if [ ! -f "$outdir/${out_file}" ] ; then
            ncap2 -O -s 'rsuscs=FSDSC-FSNSC' tmp.nc ${out_file}
            ncks -O -v rsuscs ${out_file} ${out_file}
        fi

        # get surface up SW at all-sky
        out_file=rsus_${exp_id[ii]}_${int_year2_4d}${int_mon_2d}-${end_year2_4d}${end_mon_2d}.nc
        if [ ! -f "$outdir/${out_file}" ] ; then
	        ncap2 -O -s 'rsus=FSDS-FSNS' tmp.nc ${out_file}
            ncks -O -v rsus ${out_file} ${out_file}
        fi

        # get TOA up SW at all-sky
        out_file=rsut_${exp_id[ii]}_${int_year2_4d}${int_mon_2d}-${end_year2_4d}${end_mon_2d}.nc
        if [ ! -f "$outdir/${out_file}" ] ; then
            ncap2 -O -s 'rsut=SOLIN-FSNT' tmp.nc ${out_file}
            ncks -O -v rsut ${out_file} ${out_file}
        fi

        # get TOA up SW at clear-sky
        out_file=rsutcs_${exp_id[ii]}_${int_year2_4d}${int_mon_2d}-${end_year2_4d}${end_mon_2d}.nc
        if [ ! -f "$outdir/${out_file}" ] ; then
            ncap2 -O -s 'rsutcs=SOLIN-FSNTC' tmp.nc ${out_file}
            ncks -O -v rsutcs ${out_file} ${out_file}
        fi

        rm tmp.nc

    else
        echo "All necessary output files are ready."
        echo `ls *${exp_id[ii]}_${int_year2_4d}${int_mon_2d}-${end_year2_4d}${end_mon_2d}.nc`
    fi

    ##=============================================================
    ## regrid data to coarse grid to run cal_RadKernel_xxx.py
    ##=============================================================
    ##<qinyi 2021-08-12 #------------------
    ## if the destinated grid is not 72x144, but 180x360. regrid it to 72x144.
    #echo ${rgr_map} | grep "72x144" 
    #if [ $? != 0 ] ; then

    #    all_file_list=`ls *${exp_id[ii]}_${int_year2_4d}${int_mon_2d}-${end_year2_4d}${end_mon_2d}.nc`
    #    for ii in ${all_file_list}
    #    do
    #        echo $ii
    #        echo 'regrid data with 180x360 to 72x144 grids'
    #        # this map file is copied from cori: /global/homes/z/zender/data/maps/map_cmip6_180x360_to_cmip6_72x144_aave.20181001.nc
    #        ncremap -m /g/g90/qin4/scripts/diag_feedback_E3SM/map_cmip6_180x360_to_cmip6_72x144_aave.20181001.nc $ii ${ii}_new
    #    done
    #fi

done 

echo "The final outdir =" $outdir
echo "Well Done!"

