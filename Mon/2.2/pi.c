#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

int main(int argc, char *argv[])
{
  struct timeval my_timeval;
  int iterations = 0;
  int inside_circle = 0;
  int i;
  double x, y, pi_estimate;

  gettimeofday(&my_timeval, NULL);
  srand48(my_timeval.tv_sec ^ my_timeval.tv_usec);

  if (argc == 2) {
    iterations = atoi(argv[1]);
  } else {
    printf("usage: circlepi ITERATIONS\n");
    exit(1);
  }

  for (i = 0; i < iterations; i++) {
    x = (drand48() - 0.5) * 2.0;
    y = (drand48() - 0.5) * 2.0;
    if (((x * x) + (y * y)) <= 1.0) {
      inside_circle++;
    }
  }
  pi_estimate = 4.0 * ((double) inside_circle / (double) iterations);
  printf("%d iterations, %d inside; pi = %f\n", iterations, inside_circle, pi_estimate);
  return 0;
}
