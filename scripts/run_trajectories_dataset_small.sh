cd ../app
make
cd build/
./cluster -i ../../datasets/curves_clustering/input_projection10.csv -c ../../config/cluster.conf -o ../../results/trajectories_dataset_small.csv --complete --init k-means++ --assign lloyd --update pam
cd ..
make clean
