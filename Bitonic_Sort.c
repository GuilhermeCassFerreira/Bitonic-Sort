#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <pthread.h>
#include <omp.h>

#define PRINTTIME
//#define SAVETOFILE

typedef int val_t;                          //alias para alteracao de tipo do array
struct timeval t1, t2;

struct sarg {                               //struct para passar argumentos para pthreads
    val_t* arr;
    int dir, len, num_layers;
};

/* funcao que inicia o array */
void populate_array(val_t* array, int len) {

    for (int i = 0; i < len; i++) {
            *(array + i) = (val_t)(rand() % 1000);
    }
    return;
}
/* funcao que checa se array esta ordenado */
int debug_check(val_t* arr, int len, int dir) {

    for (val_t* p = arr; p < (arr + len - 1); p++) {
        val_t* q =  p + 1;
		if ((*p!=*q) && (dir == (*p > *q))) {
            return 0;
		}
	}
    return 1;
}

/* funcao de swap usada em todos os algoritmos */
void swap(val_t* p, val_t* q) {

    val_t tmp;
    tmp = *p;
    *p = *q;
    *q = tmp;
    return;

}

/* bubble sort tirado de:
https://www.geeksforgeeks.org/bubble-sort/
*/
void bubble_sort(val_t* arr, int len, int dir) {

    int i, j, swapped;

    for (i = 0; i < len - 1; i++) {
        swapped = 0;
        for (j = 0; j < len - i - 1; j++) {
            if (dir == (*(arr + j) > *(arr + j + 1))) {
                swap(arr + j, arr + j + 1);
                swapped = 1;
            }
        }

        if (!swapped)
            break;
    }
    return;
}

/*=====quick sort de particao simples======*/

/* tirado de:
https://www.programiz.com/dsa/quick-sort
*/

int partition(val_t* arr, int low, int high, int dir) {

    int pivot = *(arr + high);

    int i = (low - 1);

    for (int j = low; j < high; j++) {
        if (dir == (*(arr + j) <= pivot)) {
            i++;
            swap(arr + i, arr + j);
        }
    }
    swap(arr + i + 1, arr + high);
    return (i + 1);
}

void quick_sort(val_t* arr, int low, int high, int dir) {

    if (low < high) {

        int pi = partition(arr, low, high, dir);
        quick_sort(arr, low, pi - 1, dir);
        quick_sort(arr, pi + 1, high, dir);
    }
    return;
}

/*======bitonic sort recursivo em camadas com uso de pthreads======*/

/* adaptado de:
https://github.com/RidipDe/Bitonic-Sorting-using-PThreads/blob/master/ParallelBitonicSort.c
https://github.com/Theodosis/parallel-programming/blob/master/pthreads/bitonic/pbitonic.c
*/

void seq_merge(val_t* arr, int dir, int len) {                  //funcao de merge recursiva sequencial

        if (len > 1) {
                for (val_t *p = arr; p < (arr + len/2); p++) {  //etapa de comparacao
                        if ((dir) == (*p > *(p + (len/2)))) {
                            swap(p, (p + len/2));
                        }
                }
                seq_merge(arr, dir, len/2);
                seq_merge(arr + len/2, dir, len/2);
        }
        return;
}

void* p_merge(void* arg) {                                       //funcao de merge recursiva paralela, usa pthreads relativas ao numero de camadas

        val_t* arr = ((struct sarg*) arg) -> arr;
        int dir = ((struct sarg*) arg) -> dir;
        int len = ((struct sarg*) arg) -> len;
        int cur_layer = ((struct sarg*) arg) -> num_layers;

        if (len > 1) {
            for (val_t *p = arr; p < (arr + len/2); p++) {      //etapa de comparacao
                if ((dir) == (*p > *(p + (len/2)))) {
                    swap(p, (p + len/2));
                }
            }
            if (cur_layer) {                                    //enquanto possui camadas disponiveis, continua dividindo em paralelo
                pthread_t thread1, thread2;
                struct sarg arg1 = {arr, dir, len/2, cur_layer-1};
                struct sarg arg2 = {arr + len/2, dir, len/2, cur_layer-1};
                pthread_create(&thread1, NULL, p_merge, (void *) &arg1);
                pthread_create(&thread2, NULL, p_merge, (void *) &arg2);

                pthread_join(thread1, NULL);
                pthread_join(thread2, NULL);

            } else {                                            //caso contrario, segue sequencialmente
                seq_merge(arr, dir, len/2);
                seq_merge(arr + len/2, dir, len/2);
            }
        }
        return NULL;
}

void seq_rec_bitonic_sort(val_t* arr, int dir, int len) {       //funcao recursiva de divisao do array, sequencial

        if (len > 1) {
                seq_rec_bitonic_sort(arr, 1, len/2);
                seq_rec_bitonic_sort(arr + (len/2), 0, len/2);
                seq_merge(arr, dir, len);
        }
        return;
}
void* p_rec_bitonic_sort(void* arg) {                            //funcao recursiva paralela de divisao, usa pthreads em camadas

        val_t* arr = ((struct sarg*) arg) -> arr;
        int dir = ((struct sarg*) arg) -> dir;
        int len = ((struct sarg*) arg) -> len;
        int cur_layer = ((struct sarg*) arg) -> num_layers;

        if (len > 1) {                                          //enquanto possui camadas disponiveis, continua dividindo em paralelo, e ordena bitonicamente as duas metades
            if (cur_layer) {
                pthread_t thread1, thread2;
                struct sarg arg1 = {arr, 1, len/2, cur_layer-1};
                struct sarg arg2 = {arr + len/2, 0, len/2, cur_layer-1};
                pthread_create(&thread1, NULL, p_rec_bitonic_sort, (void *) &arg1);
                pthread_create(&thread2, NULL, p_rec_bitonic_sort, (void *) &arg2);

                pthread_join(thread1, NULL);
                pthread_join(thread2, NULL);
                struct sarg arg3 = {arr, dir, len, cur_layer-1};
                p_merge(&arg3);
            } else {                                            //caso contrario, segue sequencialmente, incluindo a ordenacao bitonica
                seq_rec_bitonic_sort(arr, 1, len/2);
                seq_rec_bitonic_sort(arr + (len/2), 0, len/2);
                seq_merge(arr, dir, len);
            }
        }
        return NULL;
}

/*=====Bitonic Sort Paralelo com uso de OpenMP=====*/

/* adaptado de:
https://people.cs.rutgers.edu/~venugopa/parallel_summer2012/bitonic_openmp.html#graph
*/
void omp_merge (val_t* arr, int dir, int len, int m) {          //funcao de merge recursiva paralela com openmp

    int i;
    int half = len/2;

    if (len>1) {
#pragma omp parallel for shared(arr, half, dir) private(i)      //na etapa de comparacao e swap, cada setor é paralelizado internamente

        for (i = 0; i < (half); i++) {
            if ((dir) == (arr[i] > arr[i + half])) {
                swap(arr + i, arr + i + half);
            }
        }
        if (half > m) {                                         //enquanto a metade da secao bitonica for maior que o tamanho do array/numero de threads, corre recursivamente e em paralelo
            omp_merge(arr, dir, half, m);
            omp_merge(arr + half, dir, half, m);
        }
    }
    return;
}

void omp_bitonic_sort(val_t* arr, int dir, int len, int num_threads) {  //funcao que torna a sequencia bitonica

    omp_set_num_threads(num_threads);
    int m = len / num_threads;
    int i, j, sub_dir;

    for (i = 2; i <= m; i = 2 * i) {                                    //divisao em sequencias de 2, 4, 8.. ate chegar ao tamanho maximo de divisao igual entre threads
#pragma omp parallel for schedule(dynamic, m) shared(i, arr, dir) private(j, sub_dir)   //para entao cada thread fazer as comparacoes e swaps de um chunk, recursiva e sequencialmente
        for (j = 0; j < len; j += i) {
            sub_dir = ((dir)==((j/i) % 2 == 0));
            seq_merge(arr + j, sub_dir, i);
        }
    }
    for (i = 2; i <= num_threads; i = 2 * i) {                          //segunda etapa da transformacao em sequencia bitonica, comeca onde a parte de cima parou. dessa vez a paralelizacao
        for (j = 0; j < num_threads; j += i) {                          //é feita dentro da funcao de merge (ver acima)

            sub_dir = ((dir)==((j/i) % 2 == 0));
            omp_merge(arr+(j*m), sub_dir, i*m, m);

        }
#pragma omp parallel for private(j)
        for (j = 0; j < num_threads; j++)
        {
            sub_dir = ((dir)==((j/i) % 2 == 0));
            seq_merge(arr+(j*m), sub_dir, m);

        }
    }
    return;
}


int main (int argc, char *argv[]) {

/* Argumento possui quatro valores:
    - Exponencial de 2 para o tamanho do array
    - Número de threads/"camadas" de threads
    - Orientação do sorting (1 para crescente, 0 para decrescente)
    - Número de testes
*/
    if( argc != 5 ) {
		printf("Fornecer quatro valores de entrada\n");
		return 0;
	}

    int len = 1 << atoi(argv[1]);                                   //tamanho do array, len = 2^argv[1]
    int num_threads = atoi(argv[2]);                                //numero de threads/camadas
	int direction = atoi(argv[3]);                                  //1 = crescente, 0 = decrescente
	int num_runs = atoi(argv[4]);                                   //numero de testes

    val_t *array = (val_t *)malloc(len*sizeof(val_t));              //alocacao de memoria para o array

    double *t_time = (double *)malloc(4*num_runs*sizeof(double));   //alocacao para guardar tempos

    printf("Rodando algoritmos de sorting em arrays de tamanho %d, %d vezes, com %d threads/camadas\n\n", len, num_runs, num_threads);

    for (int t = 0; t < num_runs; t++) {

        printf("teste numero %d\n", t);

        populate_array(array, len);

        gettimeofday(&t1, NULL);
        //bubble_sort(array, len, direction);
        gettimeofday(&t2, NULL);
        t_time[4*t] = (t2.tv_sec - t1.tv_sec) + ((t2.tv_usec - t1.tv_usec)/1000000.0);
        printf(debug_check(array, len, direction) ? "ORDENACAO OK\n" : "ARRAY NAO ORDENADO\n");
        #ifdef PRINTTIME
        printf("Tempo do bubble sort: %f\n", t_time[4*t]);
        #endif // PRINTTIME

        populate_array(array, len);

        gettimeofday(&t1, NULL);
        quick_sort(array, 0, len - 1, direction);
        gettimeofday(&t2, NULL);
        t_time[4*t+1] = (t2.tv_sec - t1.tv_sec) + ((t2.tv_usec - t1.tv_usec)/1000000.0);
        printf(debug_check(array, len, direction) ? "ORDENACAO OK\n" : "ARRAY NAO ORDENADO\n");
        #ifdef PRINTTIME
        printf("Tempo do quick sort: %f\n", t_time[4*t+1]);
        #endif // PRINTTIME

        populate_array(array, len);

        gettimeofday(&t1, NULL);
        struct sarg arg0 = {array, direction, len, num_threads};
        p_rec_bitonic_sort(&arg0);
        gettimeofday(&t2, NULL);
        t_time[4*t+2] = (t2.tv_sec - t1.tv_sec) + ((t2.tv_usec - t1.tv_usec)/1000000.0);
        printf(debug_check(array, len, direction) ? "ORDENACAO OK\n" : "ARRAY NAO ORDENADO\n");
        #ifdef PRINTTIME
        printf("Tempo do bitonic sort com uso de pthreads: %f\n", t_time[4*t+2]);
        #endif // PRINTTIME

        populate_array(array, len);

        gettimeofday(&t1, NULL);
        omp_bitonic_sort(array, direction, len, num_threads);
        gettimeofday(&t2, NULL);
        t_time[4*t+3] = (t2.tv_sec - t1.tv_sec) + ((t2.tv_usec - t1.tv_usec)/1000000.0);
        printf(debug_check(array, len, direction) ? "ORDENACAO OK\n" : "ARRAY NAO ORDENADO\n");
        #ifdef PRINTTIME
        printf("Tempo do bitonic sort com uso do openmp: %f\n", t_time[4*t+3]);
        #endif // PRINTTIME

    }

    free(array);
    free(t_time);

    #ifdef SAVETOFILE
    FILE *myfile= fopen ("results.csv","w");
    for (int i=0; i < num_runs; i++) {
        for (int j=0; j< 4; j++)
            fprintf(myfile, "%f\n",t_time[4*j + i]);
    }
    fclose(myfile);
    #endif // SAVETOFILE

    return 0;
}
