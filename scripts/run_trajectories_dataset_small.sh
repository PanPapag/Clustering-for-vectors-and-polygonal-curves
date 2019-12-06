cd ../app
make
cd build/
<<<<<<< HEAD
./cluster -i ../../datasets/curves_clustering/input_projection6.csv -c ../../config/cluster.conf -o ../../results/input_projection6.csv --complete --init random --assign lloyd --update pam
=======
./cluster -i ../../datasets/curves_clustering/input_projection10.csv -c ../../config/cluster.conf -o ../../results/trajectories_dataset_small.csv --complete --init k-means++ --assign range-lsh --update pam
>>>>>>> c1183df3b3529bb1e567df35e761454a8a6c634c
cd ..
make clean
