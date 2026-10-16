for file in *.$3
do
echo $file
replace $1 $2 < $file > junk.txt
mv junk.txt $file
done
