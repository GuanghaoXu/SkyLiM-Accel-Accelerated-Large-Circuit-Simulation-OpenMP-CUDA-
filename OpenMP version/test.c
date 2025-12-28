#include <stdio.h>
#include <omp.h>

int main() {
    double A[1000];
    omp_set_num_threads(4);
    #pragma omp parallel
    {
        int ID = omp_get_thread_num();
        printf("Hello, World from C! (Thread %d)\n", ID);
        printf("GCC Version: %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
        printf("This is a C program test!\n");
        printf("Thread ID: %d\n", ID);
    }
    return 0;
}
