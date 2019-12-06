cd ../app
make
cd build/
./cluster -i ../../datasets/curves_clustering/input_projection6.csv -c ../../config/cluster.conf -o ../../results/input_projection6.csv --complete --init random --assign lloyd --update pam
cd ..
make clean
