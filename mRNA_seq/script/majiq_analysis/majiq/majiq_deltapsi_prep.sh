current_dir=`pwd`
echo $current_dir

cd majiq_build_res
all_pairs=`ls`
for pair in $all_pairs
do
    tissue1=${pair%_*}
    tissue2=${pair#*_}

    cd $pair
    #echo 'In '$pair'.'

    files_1="`ls -d "$PWD/"*| grep $tissue1'_ENC' | grep '\.majiq'`"
    files_2="`ls -d "$PWD/"*| grep $tissue2'_ENC' | grep '\.majiq'`"
    echo 'python3 $MAJIQ deltapsi -j 8 -o '$current_dir'/majiq_deltapsi_res/'$pair' -n '$tissue1' '$tissue2' -grp1 '$files_1' -grp2 '$files_2 >> $current_dir/majiq_deltapsi.sh
    cd ..
done
cd ..


