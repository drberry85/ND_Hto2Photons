#!/bin/bash
#1 Working Directory
#2 Command Options
#3 Castor File
#4 Output Directory
REMOTEDIR=$PWD
cd $1
rfcp $3 $REMOTEDIR
eval `scramv1 runtime -sh`
FILENAME=`basename $3`
VertexError "$2" $REMOTEDIR/$FILENAME
if [ ! -d "$4" ]; then
	mkdir $4
fi
mv "Vertex_$FILENAME" $4
