cd ../app
make
cd build/
./cluster -i ../../datasets/vectors_clustering/DataVectors_5_500x100.csv -c ../../config/cluster.conf -o ../../results/DataVectors_5_500x100 --complete --init k-means++ --assign range-lsh --update pam
cd ..
make clean
