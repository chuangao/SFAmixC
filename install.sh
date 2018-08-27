

mkdir -p ./include
mkdir -p ./lib

tar -xzvf gsl-2.2.tar.gz
cd gsl-2.2

./configure --prefix=$PWD
make
make install
cp -r include/gsl ../include/
cp -r lib/* ../lib/

#cd ../

#mkdir -p ./include/R
#echo "loc = R.home('include');write.table(loc,'loc_include.txt',row.names=F,col.names=F,sep=',',quote=F)"  | R --vanilla
#locinclude=$(cat loc_include.txt)
#cp -r $locinclude/* ./include/R/

#mkdir -p ./lib
#echo "loc = R.home('lib');write.table(loc,'loc_lib.txt',row.names=F,col.names=F,sep=',',quote=F)"  | R --vanilla
#loclib=$(cat loc_lib.txt)
#cp -r $loclib/* ./lib/

#ls
