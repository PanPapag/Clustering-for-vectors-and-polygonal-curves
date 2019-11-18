cd ../app
make
cd build/
./cluster -i ../../datasets/curves_clustering/trajectories_dataset_small.csv -c ../../config/cluster.conf -o ../../results/trajectories_dataset_small.csv --complete
cd ..
make clean
