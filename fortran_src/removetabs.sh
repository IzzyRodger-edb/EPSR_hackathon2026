for file in *.$1
do
expand -t 6 $file > junk.txt
mv junk.txt $file
done
