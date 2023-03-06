#!/bin/bash

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
regrid_SE2FV=${12}

run_id=(${run_id1} ${run_id2})
int_years=(${int_year1} ${int_year2})
end_years=(${end_year1} ${end_year2})
exp_id=(FC5 FC5_4K)

ncase=2 # two cases: one is control simulation, the other is plus4K simulation
echo $ncase

# regular variables 
reg_var_list=(FISCCP1_COSP FSDSC FSNSC TREFHT FSNT FSNTC FLUT FLUTC TS T FSDS FSNS Q SOLIN FSUTOA FSUTOAC PSL PS CLOUD CLDLIQ CLDICE PRECC PRECL U V OMEGA TGCLDCWP TGCLDIWP TGCLDLWP CLDLOW CLDLOW_CAL Z3 CLDMED CLDHGH CLDTOT RELHUM LHFLX SHFLX SWCF LWCF PBLH FREQZM CAPE QFLX QREFHT)
reg_var_new_list=(FISCCP1_COSP rsdscs rsnsc tas rsnt rsntcs rlut rlutcs ts T rsds rsns Q rsdt FSUTOA FSUTOAC psl ps CLOUD CLDLIQ CLDICE PRECC PRECL U V OMEGA TGCLDCWP TGCLDIWP TGCLDLWP CLDLOW CLDLOW_CAL Z3 CLDMED CLDHGH CLDTOT RELHUM LHFLX SHFLX SWCF LWCF PBLH FREQZM CAPE QFLX QREFHT)

#var_list=(DPDLFLIQ FREQZM dlf_out du_out ql_out mb_out dltaa_out dadt_out cape_out dcape_out deltat_out deltaq_out t_out q_out t_star_out q_star_out TUQ TVQ TUH TVH)
#var_new_list=(DPDLFLIQ FREQZM dlf_out du_out ql_out mb_out dltaa_out dadt_out cape_out dcape_out deltat_out deltaq_out t_out q_out t_star_out q_star_out TUQ TVQ TUH TVH)

# tendency variables 
ten_var_list=(DPDLFLIQ DPDLFICE DPDLFT TTEND_CLUBB RCMTEND_CLUBB RIMTEND_CLUBB RVMTEND_CLUBB TTGWORO DTCORE EVAPTZM EVAPQZM ZMDT ZMDQ ZMDLIQ ZMDICE ZMMTT FZSNTZM EVSNTZM PTTEND PTEQ PTECLDLIQ PTECLDICE DTCOND DCQ DCCLDLIQ DCCLDICE TTEND_TOT MPDT MPDQ MPDLIQ MPDICE EVAPSNOW EVAPPREC QISEVAP QCSEVAP QIRESO QCRESO PSACWSO PRACSO MPDW2I MNUCCRO CMEIOUT BERGSO MELTSDT FRZRDT QRL QRS MPDW2P MPDW2I MPDW2V QCSEDTEN PRAO PRCO PSACWSO BERGSO)
ten_var_new_list=(DPDLFLIQ DPDLFICE DPDLFT TTEND_CLUBB RCMTEND_CLUBB RIMTEND_CLUBB RVMTEND_CLUBB TTGWORO DTCORE EVAPTZM EVAPQZM ZMDT ZMDQ ZMDLIQ ZMDICE ZMMTT FZSNTZM EVSNTZM PTTEND PTEQ PTECLDLIQ PTECLDICE DTCOND DCQ DCCLDLIQ DCCLDICE TTEND_TOT MPDT MPDQ MPDLIQ MPDICE EVAPSNOW EVAPPREC QISEVAP QCSEVAP QIRESO QCRESO PSACWSO PRACSO MPDW2I MNUCCRO CMEIOUT BERGSO MELTSDT FRZRDT QRL QRS MPDW2P MPDW2I MPDW2V QCSEDTEN PRAO PRCO PSACWSO BERGSO)

# CLUBB related variables 
    clb1_var_list=(WP2_CLUBB WP3_CLUBB WPRCP_CLUBB WPRTP_CLUBB WPTHLP_CLUBB WPTHVP_CLUBB CLOUDCOVER_CLUBB CLOUDFRAC_CLUBB RTP2_CLUBB RTPTHLP_CLUBB THLP2_CLUBB SKW_ZM_af SKW_ZM_bf)
clb1_var_new_list=(WP2_CLUBB WP3_CLUBB WPRCP_CLUBB WPRTP_CLUBB WPTHLP_CLUBB WPTHVP_CLUBB CLOUDCOVER_CLUBB CLOUDFRAC_CLUBB RTP2_CLUBB RTPTHLP_CLUBB THLP2_CLUBB SKW_ZM_af SKW_ZM_bf)

# additional diagnostic variables from CLUBB
clb2_var_list=(DPBLH PBLHI PBLHOLD FRQPBL ZLCL DCPPBL DTHVINV DTINV DZINV ENTEFF ENTEFFN ENTEFFP ENTEFF_ALL ENTEFF_ALLN ENTEFF_ALLP ENTRAT ENTRATN ENTRATP ENTRAT_ALL ENTRAT_ALLN ENTRAT_ALLP)
clb2_var_new_list=(DPBLH PBLHI PBLHOLD FRQPBL ZLCL DCPPBL DTHVINV DTINV DZINV ENTEFF ENTEFFN ENTEFFP ENTEFF_ALL ENTEFF_ALLN ENTEFF_ALLP ENTRAT ENTRATN ENTRATP ENTRAT_ALL ENTRAT_ALLN ENTRAT_ALLP)

# merge as the final list
#var_list=( "${reg_var_list[@]}" "${ten_var_list[@]}" "${clb1_var_list[@]}" "${clb2_var_list[@]}" )
#var_new_list=( "${reg_var_new_list[@]}" "${ten_var_new_list[@]}" "${clb1_var_new_list[@]}" "${clb2_var_new_list[@]}" )

var_list=( "${reg_var_list[@]}" )
var_new_list=( "${reg_var_new_list[@]}" )

echo ${var_list[@]}
echo ${var_new_list[@]}


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

    outdir=${outdir_out}/${run_id[ii]}

	if [ ! -d "${outdir}/tmp" ] ; then
		mkdir -p ${outdir}/tmp
    fi

    echo ${outdir}
    cd ${outdir}

    # create symbolic links for input files 
    for (( year=${int_year}; year<=${end_year}; year++ ))
    do
      YYYY=`printf "%04d" ${year}`
      for file in ${datadir}/${run_id[ii]}.${comp}.${YYYY}-*.nc
      do
        if [ -f "$file" ] ; then
            ln -sf ${file} $outdir/tmp/
        fi
      done
    done
    
    # Generate time series files for each variable 
    for ivar in `seq 0 $[$nvar-1]`
    do
        var=${var_list[ivar]}
        varout=${var_new_list[ivar]}
        echo ${var}, ${varout}

        outfile=${varout}_${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc

        if [ ! -f "${outfile}" ] ; then
            ncclimo \
            -c ${run_id[ii]} \
            -v ${var} \
            --yr_srt=${int_year} \
            --yr_end=${end_year} \
            --map=${rgr_map} \
            -o trash \
            -O output \
            ${outdir}/tmp/${run_id[ii]}.${comp}.????-*.nc
    
            # rename the output file to the desired format
            mv ${outdir}/output/${var}* ${outdir}/${outfile}

            # rename from var_list to var_new_list if varout != var
            if [ "${var}" != "${varout}" ] ; then
                ncrename -v ${var},${varout} ${outdir}/${outfile}
            fi

        else
            echo ${outfile} "is ready"
        fi
    done # for ivar 

    # Process some specific variables through calculation
    ## prepare unit file for necessary variables 
    append=${exp_id[ii]}_${int_year_4d}${int_mon_2d}-${end_year_4d}${end_mon_2d}.nc

    cd ${outdir}

    ## get surface up SW at clear-sky
    out_file=rsuscs_${append}
    if [ ! -f "$outdir/${out_file}" ] ; then
        ncks -A rsdscs_${append} tmp.nc
        ncks -A rsnsc_${append} tmp.nc 
        ncap2 -O -s 'rsuscs=rsdscs-rsnsc' tmp.nc ${out_file}
        ncks -O -v rsuscs ${out_file} ${out_file}
        rm -f tmp.nc 
    fi

    # get surface up SW at all-sky
    out_file=rsus_${append}
    if [ ! -f "$outdir/${out_file}" ] ; then
        ncks -A rsds_${append} tmp.nc
        ncks -A rsns_${append} tmp.nc 
        ncap2 -O -s 'rsus=rsds-rsns' tmp.nc ${out_file}
        ncks -O -v rsus ${out_file} ${out_file}
        rm -f tmp.nc
    fi

    # get TOA up SW at all-sky
    out_file=rsut_${append}
    if [ ! -f "$outdir/${out_file}" ] ; then
        ncks -A rsdt_${append} tmp.nc
        ncks -A rsnt_${append} tmp.nc 
        ncap2 -O -s 'rsut=rsdt-rsnt' tmp.nc ${out_file}
        ncks -O -v rsut ${out_file} ${out_file}
        rm -f tmp.nc
    fi

    # get TOA up SW at clear-sky
    out_file=rsutcs_${append}
    if [ ! -f "$outdir/${out_file}" ] ; then
        ncks -A rsdt_${append} tmp.nc
        ncks -A rsntcs_${append} tmp.nc 
        ncap2 -O -s 'rsutcs=rsdt-rsntcs' tmp.nc ${out_file}
        ncks -O -v rsutcs ${out_file} ${out_file}
        rm -f tmp.nc
    fi

done

