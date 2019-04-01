folder=~/Dropbox/PhD/scripts/msm

cp $folder/EMMA/* .

mkdir tpt its matrix pcca discretized chapman

for i in $(seq 6)
do
    $folder/extractCentroid.sh i
done
