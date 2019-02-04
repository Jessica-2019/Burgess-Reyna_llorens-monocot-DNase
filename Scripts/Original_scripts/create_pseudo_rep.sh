# Shell script to create pseduo reps
# The input needs to be a TAGALIGN file
# Adapted from https://sites.google.com/site/anshulkundaje/projects/idr

# Define outfiles
echo "Running Pseudo Rep generator"
outputStub="${1::-12}.IDR."
pr1ALIGN="${outputStub}pr1.tagAlign.gz"
pr2ALIGN="${outputStub}pr2.tagAlign.gz"

nlines=$( zcat $1 | wc -l ) # Number of reads in the tagAlign file
nlines=$(( (nlines + 1) / 2 )) # half that number
echo `zcat $1 | shuf | split -d -l ${nlines} - ${outputStub}` # This will shuffle the lines in the file and split it into two parts
echo `gzip ${outputStub}00`
echo `gzip ${outputStub}01`
echo `mv ${outputStub}00.gz ${pr1ALIGN}`
echo `mv ${outputStub}01.gz ${pr2ALIGN}`
