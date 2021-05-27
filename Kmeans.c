#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include <omp.h>
#define ROWS 3000
#define COLS 784
#define MIN 0
#define MAX 257

/*Preprocessing and Utils*/
int* IntVecStore(int size){
	int *array = (int *) malloc(size*sizeof(int));
	return array;
}
int ** IntMatStore(int rows, int cols){
	int i;
	int **mat = (int **) malloc(rows*sizeof(int *));
	for(i=0;i<rows;i++){
		mat[i] = IntVecStore(cols);
	}
	return mat;
}
double* DoubleVecStore(int size){
	double *array = (double *) malloc(size*sizeof(double));
	return array;
}
double ** DoubleMatStore(int rows, int cols){
	int i;
	double **mat = (double **) malloc(rows*sizeof(double *));
	for(i=0;i<rows;i++){
		mat[i] = DoubleVecStore(cols);
	}
	return mat;
}
double ** doubleMatConvert(double* data,int rows, int cols){
	int i,j;
    double **array= (double **) malloc(rows*sizeof(double*));

    for (i=0; i<rows; i++)
		for (j=0; j<cols; j++) 
        	array[i] = &(data[cols*i]);

    return array;
}
int ** IntMatConvert(int* data,int rows, int cols){
	int i,j;
    int **array= (int **)malloc(rows*sizeof(double*));

    for (i=0; i<rows; i++)
		for (j=0; j<cols; j++) 
        	array[i] = &(data[cols*i]);



    return array;
}

double ** IntMatToDoubleMat(int **data, int rows, int cols){
	int i,j;
	double **newdata = DoubleMatStore(rows,cols);
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++) 
        	newdata[i][j] = (double) data[i][j];
    return newdata;
}
void DisplayDoubleMat(char name[],double **data,int rows, int cols){
    int i,j;
    	printf("------------------------%s-----------------------\n",name);
    	for(i=0;i<rows;i++){
	        for(j=0;j<cols;j++){
	            printf("%f ",data[i][j]);
	        }
        	printf("\n");
    	}
    printf("\n\n");
    
}

void DisplaySparseIntMat(char name[],int **data,int rows, int cols){
    int i,j;

    printf("--------------------------------%s------------------------\n",name);
    for(i=0;i<rows;i++){
	    for(j=0;j<cols;j++){
	        if(data[i][j]){
	        	printf("(%d,%d):%d ",i,j,data[i][j]);
	        	
	        }       
	    }
	    printf("\n");
    }
    printf("\n\n");
}	        
void DisplayIntVec(char name[],int *assignment,int size){
    int i;
    printf("------------------------%s-----------------------\n",name);
    for(i=0;i<size;i++){
        printf("%d ",assignment[i]);
        printf("\n");
    }
    printf("\n\n");
}
void DisplayDoubleVec(char name[],double *assignment,int size){
    int i;
    printf("------------------------%s-----------------------\n",name);
    for(i=0;i<size;i++){
        printf("%f ",assignment[i]);
        printf("\n");
    }
    printf("\n\n");
}
int ** ReadInCSV(char location[],int rows, int cols){
	int i,j;
	i=0;
	int N = rows*cols;
	//int counts=0;
	int **data = IntMatStore(rows,cols);
	char *record, *line;
	FILE* fstream = fopen(location, "r");
	int size = cols*5*sizeof(char);
	char* buffer = (char *) malloc(size);
	if (fstream != NULL){
		printf("Read Success\n");
	}else{
		printf("Read Fail, fstream=null\n");
	}
	printf("size:%d vectors with size %d\n",rows,cols);
	while((line=fgets(buffer,size,fstream))!=NULL && i<rows){
		record = strtok(line,",");
		j=0;
		while(record != NULL && j<cols){
			//counts+=1;
			//printf("%d ",counts);
			data[i][j]=atoi(record);
			j++;
			record = strtok(NULL,",");
		}
		++i;
	}
	return data;

}
/*Export Data to Plot in Python*/
void export(double **data, double ** meanPoints, int *clustering, int rows, int cols, int k){
    FILE *datafile = fopen("data.csv", "w");
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            fprintf(datafile, "%g, ", data[i][j]);
        }
        fprintf(datafile, "\n");
    }
    datafile = fopen("meanPoints.csv", "w");
    for(int i = 0; i < k; i++) {
        for(int j = 0; j < cols; j++) {
            fprintf(datafile, "%g, ", meanPoints[i][j]);
        }
        fprintf(datafile, "\n");
    }
    datafile = fopen("clustering.csv", "w");
    for(int i = 0; i < rows; i++) {
  
        fprintf(datafile, "%d, ", clustering[i]);
    }
}

double* zero(int size){
    double *zero = (double *) malloc(size*sizeof(double));
    memset(zero,0,size);
    return zero;
}
double* VecSum(double* point1, double* point2, int size){
    int i;
    double *sum = (double *) malloc(size*sizeof(double));
    for(i=0;i<size;i++){
        sum[i] = point1[i]+point2[i];
    }
    return sum;
}

double* para_VecSum(double* point1, double* point2, int size){
	int i;
    double *sum = (double *) malloc(size*sizeof(double));
    #pragma omp parallel for
    for(i=0;i<size;i++){
        sum[i] = point1[i]+point2[i];
    }
    return sum;
}
/*Random Data*/
double randomNum(double min, double max){
    return (max - min) * ( (double)rand() / (double)RAND_MAX ) + min;
}

double* randomDoublePoints(int size, double min, double max){
    int i;
    double *data = (double *) malloc(size*sizeof(double));

    for(i=0;i<size;i++){
        data[i]= randomNum(min,max);
    }

    return data;
}

double** randomDoubleMat(int rows,int cols, double min, double max){
    int i;
    double **data = (double **) malloc(rows*sizeof(double*));
    for(i=0;i<rows;i++){
        data[i] = randomDoublePoints(cols,min,max);
    }
    return data;
}

int * randomIntPoints(int size, int min, int max){
    int i;
    int *data = (int *) malloc(size*sizeof(int));
    for(i=0;i<size;i++){
        data[i] = rand()%(max-min)+min;
    }
    return data;
}
/*K means Clustering Non-parallel version*/

double square(double a){
    return a*a;
}
double squared_l2_distance(double *point1, double* point2, int size){
    int i;
    double sum;
    sum=0;
    for(i=0;i<size;i++){
        sum += square(point1[i]-point2[i]);
    }
    return sum;
}
double variance(int *assignment, double** data, double **meanPoints,int rows, int cols){
    int i;
    double sum=0;
    for(i=0;i<rows;i++){
        sum += squared_l2_distance(data[i],meanPoints[assignment[i]],cols);
    }
    return sum;
}
int* assignment(double **meanPoints, double **data, int rows, int cols, int k){
    int i,j,min_index;
    double distance,min;
    int *clustering = (int *) malloc(rows*sizeof(int));
    for(i=0;i<rows;i++){
        min = squared_l2_distance(data[i],meanPoints[0],cols);
 
        min_index=0;
        for(j=1;j<k;j++){
            distance = squared_l2_distance(data[i],meanPoints[j],cols);
            if(min>distance){
                min=distance;
                min_index=j;
            }
        }
        clustering[i] = min_index;
    }

    return clustering;
}
double** update(int *assignment, double **data, int rows, int cols, int k){
    int i,j;
    double **meanPoints = (double **) malloc(k*sizeof(double*));
    int *counts = (int *) zero(k);

    for(i=0;i<k;i++){

        meanPoints[i]=zero(cols);
    }

    for(i=0;i<rows;i++){
        meanPoints[assignment[i]]=VecSum(meanPoints[assignment[i]],data[i],cols);
        // if(assignment[i]==2){
        // 	DisplayDoubleVec("meanPoints",meanPoints[3],cols);
        // }
        counts[assignment[i]]+=1;
    }

    for(i=0;i<k;i++){
        for(j=0;j<cols;j++){
            meanPoints[i][j] /= (double) counts[i];
        }
    }
    
    return meanPoints;
}

double para_variance(int *assignment, double** data, double **meanPoints,int rows, int cols){
    int i;
    double sum=0;
    #pragma omp parallel for reduction(+:sum)
    for(i=0;i<rows;i++){
        sum += squared_l2_distance(data[i],meanPoints[assignment[i]],cols);
    }
    return sum;
}

int* para_assignment(double **meanPoints, double **data, int rows, int cols, int k){
    int i,j,min_index;
    long double distance,min;
    int *clustering = (int *) malloc(rows*sizeof(int));
    #pragma omp parallel for private(distance,min,min_index,j)
    for(i=0;i<rows;i++){
        min = squared_l2_distance(data[i],meanPoints[0],cols);
 		
        min_index=0;
        for(j=1;j<k;j++){
            distance = squared_l2_distance(data[i],meanPoints[j],cols);
            if(min>distance){
                min=distance;
                min_index=j;
            }
        }
        clustering[i] = min_index;
    }
    return clustering;
}
double** para_update(int *assignment, double **data, int rows, int cols, int k){
    int i,j;
    double **meanPoints = (double **) malloc(k*sizeof(double*));
    int *counts = (int *) zero(k);

    #pragma omp parallel for shared(meanPoints)
    for(i=0;i<k;i++){

        meanPoints[i]=zero(cols);
    }

    for(i=0;i<rows;i++){
        meanPoints[assignment[i]]=para_VecSum(meanPoints[assignment[i]],data[i],cols);
        counts[assignment[i]]+=1;
    }

    #pragma omp parallel for shared(meanPoints,counts) collapse(2)
    for(i=0;i<k;i++){
        for(j=0;j<cols;j++){
            meanPoints[i][j] /= counts[i];
        }
    }

    return meanPoints;
}


void NormalKMC(double **data, int rows, int cols, int k, int min, int max,int iters){
    int i;
    //Initializing Clustering Randomly
    double var;
    clock_t tic = clock();
    srand(time(NULL));
    double **meanPoints = (double **) malloc(k*sizeof(double*));
    int *clustering = (int *) malloc(rows*sizeof(int));
    
    clustering = randomIntPoints(rows,0,k);
    for(i=0;i<k;i++){
    	clustering[i]=i;
    }
    
    //algorithm iteration part
    for(i=0;i<iters;i++){
    	//DisplayIntVec("clustering",clustering,rows);
    	free(meanPoints);
    	meanPoints = update(clustering,data,rows,cols,k);
    	//DisplayDoubleMat("meanPoints",meanPoints,k,cols);
    	free(clustering);
        clustering = assignment(meanPoints,data,rows,cols,k);
        
        
        var = variance(clustering,data,meanPoints,rows,cols); 
          
    }

    clock_t toc = clock();
    
    printf("Normal %d means Clustering:%fs\n",k,(double) (toc-tic)/CLOCKS_PER_SEC);
    printf("Variance:%f\n\n",var);
    /*the final step is to save the data in the .csv files.*/
    //export(data,meanPoints,clustering,rows,cols,k);
    
    
}

void ParallelKMC(double **data, int rows, int cols, int k, int min, int max, int iters){
    int i;
    //Initializing Clustering Randomly
    double var;
    double tt = omp_get_wtime();
    srand(time(NULL));
    double **meanPoints = (double **) malloc(k*sizeof(double*));
    int *clustering = (int *) malloc(rows*sizeof(int));

    clustering = randomIntPoints(rows,0,k);
    for(i=0;i<k;i++){
    	clustering[i]=i;
    }

    //algorithm iteration part
    for(i=0;i<iters;i++){
    	//DisplayIntVec("clustering",clustering,rows);
    	free(meanPoints);
    	meanPoints = para_update(clustering,data,rows,cols,k);
    	//DisplayDoubleMat("meanPoints",meanPoints,k,cols);
    	free(clustering);
        clustering = para_assignment(meanPoints,data,rows,cols,k);
        
        var = variance(clustering,data,meanPoints,rows,cols);
         
    }
    
    printf("Parallel %d means Clustering:%fs\n",k,omp_get_wtime()-tt);
    printf("Variance:%f\n\n",var);  
    /*the final step is to save the data in the .csv files.*/
    export(data,meanPoints,clustering,rows,cols,k);
    
}

int main(int argc, char** argv) {
    srand(time(NULL));
    int k,iters;
    char *filename;
    k = atoi(argv[1]);
    iters = atoi(argv[2]);
    filename = argv[3];
    int **mat = (int **) malloc(ROWS*sizeof(int*));
    double **dmat = DoubleMatStore(ROWS,COLS);
    mat = ReadInCSV(filename,ROWS,COLS); 
    dmat = IntMatToDoubleMat(mat,ROWS,COLS); 
 
    //DisplaySparseIntMat("Data",mat,ROWS,COLS); 
    //for(k=1;k<=10;k++){
    	//NormalKMC(dmat,ROWS,COLS,k,MIN,MAX,iters);
    	ParallelKMC(dmat,ROWS,COLS,k,MIN,MAX,iters);
    //}
    

    free(dmat);
    free(mat);
    return 0;

}