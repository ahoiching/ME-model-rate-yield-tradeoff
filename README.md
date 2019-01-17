mkdir rate_yield
virtualenv ./rate_yield/
cd rate_yield/
source ./bin/activate

git clone https://github.com/SBRG/cobrame.git
cd cobrame
python setup.py develop


cd qminos1114/
cp Makefile.defs ./minos56/
cp Makefile.defs ./qminos56/

cd minos56/
make clean
make
 
cd ..
cd qminos56/
make clean
make

cd ..
cd ..
git clone https://github.com/SBRG/solvemepy.git 
cd solvemepy/
cp ../qminos1114/qminos56/lib/libquadminos.a ./
cp ../qminos1114/minos56/lib/libminos.a ./
sudo apt-get install gfortran

git clone https://github.com/SBRG/ecolime.git
cd ecolime/
python setup.py develop
