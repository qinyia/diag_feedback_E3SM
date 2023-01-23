#****************************************************************
#
#    Filename: archive_data.sh
#
#    Author: Yi Qin - qin4@llnl.gov
#    Description: run case.st_archive to archive data to archive directory
#    Input: 
#    Output: 
#    Create: 2022-01-19 12:21:41
#    Last Modified: 2022-01-19 12:21:41
#****************************************************************


run_id1=$1 #cori-haswell.20190513.F2010C5-CMIP6-LR.ne30_oECv3
run_id2=$2 #cori-haswell.20190513.F2010C5-CMIP6-LR.plus4K.ne30_oECv3
datadir_in1=$3
datadir_in2=$4

echo ${run_id1}
echo ${run_id2}
echo ${datadir_in1}
echo ${datadir_in2}

run_id=(${run_id1} ${run_id2})
echo ${run_id}

ncase=2 # two cases: one is control simulation, the other is plus4K simulation
echo $ncase

for ii in `seq 0 $[ncase-1]`
do

    if [ "${run_id[ii]}" == "20200428.DECKv1b_amip1-CFMIP.ne30_oEC.cori-knl-L" ] ; then
        datadir=${datadir_in1}/${run_id[ii]}/${datadir_in2}
    else
        datadir=${datadir_in1}/${run_id[ii]}/${datadir_in2}
    fi

    echo $datadir
    cd $datadir 

    #./case.st_archive

    if [ ! -d "${datadir_in1}/${run_id[ii]}/archive" ] ; then
        echo ">>>>>> Start archiving."
        ./case.st_archive
    else
        echo ">>>>>> Already archived."
    fi

done 

echo "Archiving data is successfully."
echo "Well Done!"
