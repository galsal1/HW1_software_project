//
// Created by galha on 27/01/2024.
//
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int K,N,d,iter=200;
char* path;
double* data;
int* places;
int* Group_size;
double* K_centroid;
double* Mean_vector;

//getline func delete in nova server
size_t getline(char **lineptr, size_t *n, FILE *stream) {
    char *bufptr = NULL;
    char *p = bufptr;
    size_t size;
    int c;

    if (lineptr == NULL) {
        return -1;
    }
    if (stream == NULL) {
        return -1;
    }
    if (n == NULL) {
        return -1;
    }
    bufptr = *lineptr;
    size = *n;

    c = fgetc(stream);
    if (c == EOF) {
        return -1;
    }
    if (bufptr == NULL) {
        bufptr = malloc(128);
        if (bufptr == NULL) {
            return -1;
        }
        size = 128;
    }
    p = bufptr;
    while(c != EOF) {
        if ((p - bufptr) > (size - 1)) {
            size = size + 128;
            bufptr = realloc(bufptr, size);
            if (bufptr == NULL) {
                return -1;
            }
        }
        *p++ = c;
        if (c == '\n') {
            break;
        }
        c = fgetc(stream);
    }

    *p++ = '\0';
    *lineptr = bufptr;
    *n = size;

    return p - bufptr - 1;
}

void initparametrs(int argc, char *argv[]){
    char *p;
    //init N
    long Nlong = strtol(argv[2], &p, 10);
    if (errno != 0 || *p != '\0' || Nlong < 2) {
        fprintf(stderr, "Invalid number of points!");
        exit(-1);
    } else {
        // No error
        N = Nlong;
    }

    //init K
    long Klong = strtol(argv[1], &p, 10);
    if (errno != 0 || *p != '\0' || Klong > N-1 || Klong < 2) {
        fprintf(stderr, "Invalid number of clusters!");
        exit(-1);
    } else {
        // No error
        K = Klong;
    }


    //init d
    long dlong = strtol(argv[3], &p, 10);
    if (errno != 0 || *p != '\0' || dlong < 1) {
        fprintf(stderr, "Invalid dimension of point!");
        exit(-1);
    } else {
        // No error
        d = dlong;
    }


    //init iter and file path
    if(argc<6){
        path= argv[4];
    } else {
        long iterlong = strtol(argv[4], &p, 10);
        if (errno != 0 || *p != '\0' || iterlong > 1000 || iterlong < 1) {
            fprintf(stderr, "Invalid maximum iteration!");
            exit(-1);
        } else {
            // No error
            iter = iterlong;
            path= argv[5];
        }
    }
}

void initK_centroid(){
    FILE * fp;
    // Open a file in read mode
    fp = fopen(path, "r");
    char * line = NULL;
    size_t len = 0;
    size_t read;
    if (fp == NULL)
        exit(EXIT_FAILURE);
    Group_size = calloc(K,sizeof(double));
    K_centroid = calloc(K*d,sizeof(double));
    // Read data from file
    for (int i = 0; i < K; i++) {
        if (getline(&line, &len, fp) != -1) {
            // Parse doubles from the line
            char *start = line;
            char *end;
            for (int j = 0; j < d; j++) {
                K_centroid[(i*d)+j] = strtod(start, &end);
                if (start == end) {
                    fprintf(stderr, "Error parsing double at line %d, column %d\n", i+1, j+1);
                    break;
                }
                start = end;
                while (*start == ',') {
                    start++;
                }
            }
        } else {
            fprintf(stderr, "File contains fewer lines than expected.\n");
            break;
        }
        Group_size[i]=0;
    }

    // Close the file
    fclose(fp);
    if (line)
        free(line);
}

void init_data(){
    FILE * fp;
    // Open a file in read mode
    fp = fopen(path, "r");
    char * line = NULL;
    size_t len = 0;
    size_t read;
    if (fp == NULL)
        exit(EXIT_FAILURE);
    data = calloc(N*d,sizeof(double));
    places = calloc(N,sizeof(int ));
    // Read data from file
    int cnt=0;
    for (int i = 0; i < N; i++) {
        places[i]=-1;
        if (getline(&line, &len, fp) != -1) {
            // Parse doubles from the line
            char *start = line;
            char *end;
            for (int j = 0; j < d; j++) {
                data[(i * d) + j] = strtod(start, &end);
                if (start == end) {
                    fprintf(stderr, "Error parsing double at line %d, column %d\n", i + 1, j + 1);
                    break;
                }
                start = end;
                while (*start == ',') {
                    start++;
                }
            }
            cnt++;
        }
    }
    //update the actual number of line
    N=cnt;
}

void groupsMean(){
    Mean_vector = calloc(K*d,sizeof(double));
    for(int i=0;i<K*d;i++) {
        Mean_vector[i]=0;
    }
    for(int i=0;i<N;i++) {
        for(int j=0;j<d;j++){
            Mean_vector[places[i]*d+j] += data[i*d+j];
        }
    }
    for(int i=0;i<K;i++) {
        for(int j=0;j<d;j++){
            Mean_vector[i*d+j] /= Group_size[i];
        }
    }
}

double All_euclidean_distance() {
    double res = 0;
    for (int i = 0; i < K; i++) {
        double temp = 0;
        for (int j = 0; j < d; j++){
            temp+= pow(K_centroid[i*d+j]-Mean_vector[i*d+j],2);
        }
        res += sqrt(temp);
    }
    return res;
}

double euclidean_distance(int vector1,int Kvector2) {
    double res = 0;
    for (int j = 0; j < d; j++){
        res+= pow(data[vector1*d+j]-K_centroid[Kvector2*d+j],2);
    }
    res = sqrt(res);
    return res;
}

int indexofSmallestElement(double array[], int size)
{
    int index = 0;
    if (size != 1)
    {

        double n = array[0];
        for (int i = 1; i < size; i++)
        {
            if (array[i] < n)
            {
                n = array[i];
                index = i;
            }
        }
    }
    return index;
}

int moveVector(int vector_index,int cnt){
    int current_group = places[vector_index];
    double distances[K] ={};
    for(int i=0;i<K;i++){
        distances[i] = euclidean_distance(vector_index,i);
    }
    int min_index = indexofSmallestElement(distances,K);
    if (min_index == current_group)
        return 0;
    else {
        places[vector_index] = min_index;
        Group_size[min_index] += 1;
        if(cnt) {
            Group_size[current_group] -= 1;
        }
        return 1;
    }
}

int main(int argc, char *argv[]) {

}
