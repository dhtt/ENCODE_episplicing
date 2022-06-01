for f in majiq_to_be_mapped/*
do
    NAME=${f##*/}
    NAME=${NAME%%.*}
    echo $f $NAME
    echo $NAME >> majiq_all_AS.txt
    awk -F';' '{print $5}' $f | grep -c True >> majiq_all_AS.txt
    awk -F';' '{print $6}' $f | grep -c True >> majiq_all_AS.txt
    awk -F';' '{print $7}' $f | grep -c True >> majiq_all_AS.txt
done
