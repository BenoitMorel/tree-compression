
deps="$PWD/deps/"
mkdir $deps

git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
./install.sh
cp $HOME/include/div* $deps
cp -r $HOME/include/sdsl $deps
cp $HOME/lib/*.a $deps
cd ..

git clone --recursive https://github.com/ddarriba/pll-modules.git
cd pll-modules
mkdir build
cd build
cmake ..
make -j 4 
cp src/*/*.a $deps
cd ../..


cd src
libpll=$deps/libpll
mkdir libpll
cp ../pll-modules/src/tree/pll_tree.h $libpll
cp ../pll-modules/libs/libpll/src/pll.h  $libpll
cp ../pll-modules/build/libs/libpll/src/libpll.a $deps
cp ../pll-modules/build/src/*/libpllmodtree.a $deps
make
 

