#!/bin/bash
#1 Working Directory
#2 Command Options
#3 Castor File
#4 Output Directory
REMOTEDIR=$PWD
cd $1
eval `scramv1 runtime -sh`
STAGE_SVCCLASS=default

#echo "#### Begin Debug Output ####"
#date
#id
#env | grep RFIO
#env | grep STAGE
#env | grep PATH
#which rfdir
#echo "#### End Debug Output ####"

cd -
rfdir $3
rfcp $3 $REMOTEDIR
FILENAME=`basename $3`
/bin/ls *.root
echo "VertexError $2 $FILENAME"
VertexError "$2" $FILENAME
if [ ! -d "$1/$4" ]; then
	mkdir $1/$4
fi
cp "Vertex_$FILENAME" $1/$4
