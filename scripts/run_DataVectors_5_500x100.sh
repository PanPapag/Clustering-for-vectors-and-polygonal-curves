cd ../app
make
cd build/
./cluster -i ../../datasets/vectors_clustering/DataVectors_5_500x100.csv -c ../../config/cluster.conf -o ../../results/DataVectors_5_500x100 --complete
cd ..
make clean
