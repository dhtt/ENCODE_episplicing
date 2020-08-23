INPUT_PATH=$1

#=== GET SUBMETADATA ===
echo "Get submetadata from metadata.tsv"
for file in $INPUT_PATH/*
do (
    ID=${file%%.*}
    echo $ID
    grep "^$ID" $INPUT_PATH/metadata.tsv >> $INPUT_PATH/submetadata.tsv
)
done
