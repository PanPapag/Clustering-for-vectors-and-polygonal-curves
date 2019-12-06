cd ../app
make
cd build/
<<<<<<< HEAD
./cluster -i ../../datasets/vectors_clustering/DataVectors_5_500x100.csv -c ../../config/cluster.conf -o ../../results/DataVectors_5_500x100 --complete --init random --assign lloyd --update pam
=======
./cluster -i ../../datasets/vectors_clustering/DataVectors_5_500x100.csv -c ../../config/cluster.conf -o ../../results/DataVectors_5_500x100 --complete --init k-means++ --assign lloyd --update pam
>>>>>>> c1183df3b3529bb1e567df35e761454a8a6c634c
cd ..
make clean
