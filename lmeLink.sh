#! /bin/sh
# $Id: lmeLink.sh,v 1.1 1999/10/13 00:50:10 saikat Exp $
# Create the symbolic links in the lme subdirectory
if [ ! -d ../NLMEDATA ]
then
   echo "both src/NLME and src/NLMEDATA must present to create the lme package"
   exit 1
fi
#  for i in *
#  do                               # remove old symbolic links
#    if [ -L $i ]
#    then 
#      rm $i
#    fi
#  done 
#  for i in ../../*.q               # link the .q files in the parent
#  do
#    ln -s $i `basename $i .q`.R
#  done
mkdir -p ./lme/data; cd ./lme/data
for i in *; do if [ -L $i ]; then rm $i; fi; done  # remove old symbolic links
for i in ../../../NLMEDATA/*.q; do ln -s $i `basename $i .q`.R; done
mkdir -p ../src; cd ../src
ln -sf ../../nlmefit.c .
cd ..
if [ -d ../../SAS_Mixed ]
then
  cd data
  for i in ../../../SAS_Mixed/data/*.q; do ln -s $i `basename $i .q`.R; done
  cd ..
  mkdir -p SAS_Mixed; cd SAS_Mixed;
  for i in *; do if [ -L $i ]; then rm $i; fi; done
  for i in ../../../SAS_Mixed/*.{tex,bib,ps}; do ln -s $i `basename $i`; done
fi
exit 0
