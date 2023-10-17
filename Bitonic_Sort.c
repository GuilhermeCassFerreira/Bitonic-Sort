#include <stdio.h>

void compareAndSwap(int a[], int i, int j, int dir) {
    if ((a[i] > a[j] && dir) || (a[i] < a[j] && !dir)) {
        int temp = a[i];
        a[i] = a[j];
        a[j] = temp;
    }
}

void bitonicMerge(int a[], int low, int cnt, int dir) {
    if (cnt > 1) {
        int k = cnt / 2;
        for (int i = low; i < low + k; i++)
            compareAndSwap(a, i, i + k, dir);
        bitonicMerge(a, low, k, dir);
        bitonicMerge(a, low + k, k, dir);
    }
}

void bitonicSort(int a[], int low, int cnt, int dir) {
    if (cnt > 1) {
        int k = cnt / 2;
        bitonicSort(a, low, k, 1);
        bitonicSort(a, low + k, k, 0);
        bitonicMerge(a, low, cnt, dir);
    }
}

void sort(int arr[], int n) {
    bitonicSort(arr, 0, n, 1);
}

int main() {
    int arr[] = {30, 11, 23, 4, 20, 10};
    int n = sizeof(arr) / sizeof(arr[0);

    printf("Vetor original: \n");
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");

    sort(arr, n);

    printf("Vetor ordenado em ordem ascendente: \n");
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");

    return 0;
}
