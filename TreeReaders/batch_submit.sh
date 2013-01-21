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
rfdir $3 >& /dev/null
if [ $? -eq 0 ]; then
	echo "File is CASTOR"
	rfdir $3
	rfcp $3 $REMOTEDIR
fi

/afs/cern.ch/project/eos/installation/0.1.0-22d/bin/eos.select ls $3 >& /dev/null
if [ $? -eq 0 ]; then
	echo "File is EOS"
	/afs/cern.ch/project/eos/installation/0.1.0-22d/bin/eos.select find -f $3
	xrdcp -np -v root://eoscms/$3 .
fi

FILENAME=`basename $3`
/bin/ls *.root
echo "VertexError $2 $FILENAME"
VertexError "$2" $FILENAME
if [ ! -d "$1/$4" ]; then
	mkdir $1/$4
fi
cp "Vertex_$FILENAME" $1/$4
