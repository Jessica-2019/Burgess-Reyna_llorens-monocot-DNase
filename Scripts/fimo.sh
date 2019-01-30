#!bin/bash

#run fimo on DGFs

fasta=$1

cd ./meme

for d in */; 
do
	samp=`basename ${d::-17}`
	cd ${d};
	echo $samp
	~/meme2/bin/fimo --o ${samp}.fimo meme.txt ${fasta}
	cd ..;

done



### Run Fimo in random DGFs generated outside the DHS

#mkdir fimo

for a in */
do
	samp=`basename ${a::-5}`
	
	cd ${a};
	echo $samp
	mv */*.fimo ../fimo
	cd ..;
done



