The c code is used to run standard algorithm for K-means clustering. 

Compile Command:
gcc -fopenmp Kmeans.c -o Kmeans
-------------------------------
Run Command:
./Kmeans k iterations path 
i.e. 
./Kmeans 10 250 ./digits_train.csv
--------------------------------
In the code, I have an export function to rewrite data.csv, meanPoints.csv, clustering.csv to label the data. In the Jupyter notebook Visualization.ipynb, I basically used python to visualize the graph for K-means clustering for random dimension 2 vectors and the meanPoints for MNIST data.
