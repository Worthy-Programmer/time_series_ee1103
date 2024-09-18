#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.141592653589793

// Function to generate normal-distributed noise using the Box-Muller transform
double generate_noise()
{
  double u1 = ((double)rand() / RAND_MAX);
  double u2 = ((double)rand() / RAND_MAX);
  return sqrt(-2 * log(u1)) * cos(2 * PI * u2);
}

// Lorentzian function
double lorentzian(double t, double a, double mT)
{
  return (a * a) / ((t - mT) * (t - mT) + a * a);
}

// Structure to define a Lorentzian peak with location, amplitude, and width
typedef struct
{
  double location;
  double amplitude;
  double width;
} LorentzianPeak;

// Generate Lorentzian peaks with added noise
void generate_lorentzian_peaks(LorentzianPeak *peaks, int M, double T, double a, double eps_T, double eps_a)
{
  for (int m = 1; m <= M; m++)
  {
    peaks[m - 1].location = m * T + generate_noise() * eps_T;
    peaks[m - 1].amplitude = 1.0; // Let's assume amplitude is always 1.0 for simplicity
    peaks[m - 1].width = a + generate_noise() * eps_a;
  }
}

// Simulate the time series with noisy Lorentzian peaks and noise added to amplitude
void simulate_time_series(double *series, int num_points, double time_span, LorentzianPeak *peaks, int M, double noise_level)
{
  double t;
  for (int i = 0; i < num_points; i++)
  {
    t = i * (time_span / num_points);
    series[i] = 0.0;
    for (int m = 0; m < M; m++)
    {
      series[i] += lorentzian(t, peaks[m].width, peaks[m].location);
    }
    series[i] += generate_noise() * noise_level; // Adding noise to amplitude
  }
}

// Estimate the average T and a, as well as their standard deviations
void estimate_parameters(LorentzianPeak *peaks, int M, double *avg_T, double *avg_a, double *std_T, double *std_a)
{
  double sum_T = 0, sum_a = 0;
  double sum_T2 = 0, sum_a2 = 0;
  for (int m = 1; m < M; m++)
  {
    double diff_T = peaks[m].location - peaks[m - 1].location;
    sum_T += diff_T;
    sum_T2 += diff_T * diff_T;

    sum_a += peaks[m].width;
    sum_a2 += peaks[m].width * peaks[m].width;
  }

  *avg_T = sum_T / (M - 1);
  *avg_a = sum_a / M;

  *std_T = sqrt((sum_T2 / (M - 1)) - (*avg_T * *avg_T));
  *std_a = sqrt((sum_a2 / M) - (*avg_a * *avg_a));
}

void display_array(double *array, int size)
{
  // Loop through each item in the array
  for (int i = 0; i < size; i++)
  {
    // Print each item to the console
    printf("Item %d: %f\n", i, array[i]);
  }
}

void plot_series_gnuplot(double *series, int size)
{
  // Open a pipe to Gnuplot
  FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

  if (gnuplotPipe == NULL)
  {
    printf("Error opening pipe to Gnuplot.\n");
    return;
  }

  // Send Gnuplot commands to setup the plot
  fprintf(gnuplotPipe, "set title 'Time Series Data'\n");
  fprintf(gnuplotPipe, "set xlabel 'Index'\n");
  fprintf(gnuplotPipe, "set ylabel 'Value'\n");
  fprintf(gnuplotPipe, "plot '-' with lines title 'Series'\n");

  // Send the array data to Gnuplot
  for (int i = 0; i < size; i++)
  {
    fprintf(gnuplotPipe, "%d %f\n", i, series[i]);
  }

  // Signal Gnuplot to end the data series
  fprintf(gnuplotPipe, "e\n");

  // Close the pipe
  pclose(gnuplotPipe);
}

// Main function
int main(int argc, char *argv[])
{
  srand(time(NULL));

  int M = 10;                // Default number of peaks
  double T = 10.0;           // Default period
  double a = 1.0;            // Default width
  double eps_T = 0.1;        // Noise level for location
  double eps_a = 0.1;        // Noise level for width
  double noise_level = 0.05; // Amplitude noise level

  if (argc == 4)
  {
    M = atoi(argv[1]);
    T = atof(argv[2]);
    a = atof(argv[3]);
  }

  int num_points = 1000;
  double time_span = M * T;
  double *series = (double *)malloc(num_points * sizeof(double));
  LorentzianPeak *peaks = (LorentzianPeak *)malloc(M * sizeof(LorentzianPeak));

  // Generate Lorentzian peaks
  generate_lorentzian_peaks(peaks, M, T, a, eps_T, eps_a);

  // Simulate time series with noise
  simulate_time_series(series, num_points, time_span, peaks, M, noise_level);

  // Estimate parameters from the simulated data
  double avg_T, avg_a, std_T, std_a;
  estimate_parameters(peaks, M, &avg_T, &avg_a, &std_T, &std_a);

  plot_series_gnuplot(series, num_points);

  // Output results
  printf("Estimated Avg T: %f, Std T: %f\n", avg_T, std_T);
  printf("Estimated Avg a: %f, Std a: %f\n", avg_a, std_a);

  // Free memory
  free(series);
  free(peaks);

  return 0;
}