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
    long Nlong;
    long Klong;
    long dlong;

    if(argc<5){
        printf( "An Error Has Occurred");
        exit(1);
    }

    Nlong = strtol(argv[2], &p, 10);
    if (*p != '\0' || Nlong < 2) {
        printf("Invalid number of points!");
        exit(1);
    } else {
        N = Nlong;
    }

    Klong = strtol(argv[1], &p, 10);
    if ( *p != '\0' || Klong > N-1 || Klong < 2) {
        printf( "Invalid number of clusters!");
        exit(1);
    } else {
        K = Klong;
    }

    dlong = strtol(argv[3], &p, 10);
    if ( *p != '\0' || dlong < 1) {
        printf( "Invalid dimension of point!");
        exit(1);
    } else {
        d = dlong;
    }

    if(argc<6){
        path= argv[4];
    } else {
        long iterlong = strtol(argv[4], &p, 10);
        if ( *p != '\0' || iterlong > 1000 || iterlong < 1) {
            printf( "Invalid maximum iteration!");
            exit(1);
        } else {
            iter = iterlong;
            path= argv[5];
        }
    }
}

void initK_centroid(){
    FILE * fp;
    char * line;
    size_t len;
    int i,j;

    fp = fopen(path, "r");
    line = NULL;
    len = 0;
    if (fp == NULL)
        exit(EXIT_FAILURE);
    Group_size = calloc(K,sizeof(double));
    K_centroid = calloc(K*d,sizeof(double));
    for (i= 0; i < K; i++) {
        if (getline(&line, &len, fp) != -1) {
            char *start = line;
            char *end;
            for (j = 0; j < d; j++) {
                K_centroid[(i*d)+j] = strtod(start, &end);
                if (start == end) {
                    fprintf(stderr, "An Error Has Occurred");
                    exit(1);
                }
                start = end;
                while (*start == ',') {
                    start++;
                }
            }
        } else {
            fprintf(stderr, "An Error Has Occurred");
            exit(1);
        }
        Group_size[i]=0;
    }

    fclose(fp);
    if (line)
        free(line);
}

void init_data(){
    FILE * fp;
    char *line,*start,*end;
    size_t len;
    int i,j,cnt;

    fp = fopen(path, "r");
    line = NULL;
    len = 0;
    if (fp == NULL)
        exit(EXIT_FAILURE);
    data = calloc(N*d,sizeof(double));
    places = calloc(N,sizeof(int ));

    cnt=0;
    for (i = 0; i < N; i++) {
        places[i]=-1;
        if (getline(&line, &len, fp) != -1) {
            start = line;
            for (j = 0; j < d; j++) {
                data[(i * d) + j] = strtod(start, &end);
                if (start == end) {
                    fprintf(stderr, "An Error Has Occurred");
                    exit(1);
                }
                start = end;
                while (*start == ',') {
                    start++;
                }
            }
            cnt++;
        }
    }
    N=cnt;

    fclose(fp);
    if (line)
        free(line);
}

void groupsMean(){
    int i,j;
    for(i=0;i<K*d;i++) {
        Mean_vector[i]=0;
    }
    for(i=0;i<N;i++) {
        for(j=0;j<d;j++){
            Mean_vector[places[i]*d+j] += data[i*d+j];
        }
    }
    for(i=0;i<K;i++) {
        for(j=0;j<d;j++){
            Mean_vector[i*d+j] /= Group_size[i];
        }
    }
}

double All_euclidean_distance() {
    int i,j;
    double res = 0;
    for (i = 0; i < K; i++) {
        double temp = 0;
        for (j = 0; j < d; j++){
            temp+= pow(K_centroid[i*d+j]-Mean_vector[i*d+j],2);
        }
        res += sqrt(temp);
    }
    return res;
}

double euclidean_distance(int vector1,int Kvector2) {
    int j;
    double res = 0;
    for (j = 0; j < d; j++){
        res+= pow(data[vector1*d+j]-K_centroid[Kvector2*d+j],2);
    }
    res = sqrt(res);
    return res;
}

int indexofSmallestElement(double array[], int size)
{
    int i;
    int index = 0;
    if (size != 1)
    {

        double n = array[0];
        for (i = 1; i < size; i++)
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
    int i,min_index;
    int current_group = places[vector_index];
    double *distances= calloc(K,sizeof(double ));
    for(i=0;i<K;i++){
        distances[i] = euclidean_distance(vector_index,i);
    }
    min_index = indexofSmallestElement(distances,K);
    if (min_index == current_group) {
        free(distances);
        return 0;
    }
    else {
        places[vector_index] = min_index;
        Group_size[min_index] += 1;
        if(cnt) {
            Group_size[current_group] -= 1;
        }
        free(distances);
        return 1;
    }

}

int main(int argc, char *argv[]) {
    int cnt,i,j,moveflag;
    double deltasum;

    initparametrs(argc,argv);
    initK_centroid();
    init_data();
    Mean_vector = calloc(K*d,sizeof(double));

    for(cnt=0;cnt<iter;cnt++){
        moveflag=0,deltasum=0;
        for(i=0;i<N;i++) {
            moveflag = moveVector(i,cnt) || moveflag ;
        }
        groupsMean();
        deltasum += All_euclidean_distance();
        K_centroid = Mean_vector;
        if(cnt!=0 && (deltasum<0.001 || !moveflag)){
            break;
        }
    }
    for(i=0;i<K;i++){
        for(j=0;j<d;j++){
            if(j<d-1) {
                printf("%.4f ,", K_centroid[i * d + j]);
            }
            else
                printf("%.4f",K_centroid[i*d+j]);
        }
        printf("\n");
    }
    return 0;
}