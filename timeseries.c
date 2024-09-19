#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.141592653589793
#define WINDOW_SIZE 10

#define MIN_PEAK_HEIGHT 0.5 // Threshold to ignore noise
#define MAX_PEAKS 100       // Maximum number of peaks to store

// Function to compute standard deviation
double compute_std(double *data, int n, double mean)
{
  if (n <= 1)
    return 0.0; // Standard deviation is not defined for n <= 1
  double sum = 0.0;
  for (int i = 0; i < n; i++)
  {
    sum += (data[i] - mean) * (data[i] - mean);
  }
  return sqrt(sum / (n - 1));
}

// Function to calculate the moving average
void moving_average(double series[], int length, double smoothed[])
{
  int i, j;
  double sum;

  // Loop through the series
  for (i = 0; i < length; i++)
  {
    sum = 0.0;
    int count = 0;

    // Sum elements within the window, handling boundaries
    for (j = i - WINDOW_SIZE / 2; j <= i + WINDOW_SIZE / 2; j++)
    {
      if (j >= 0 && j < length)
      {
        sum += series[j];
        count++;
      }
    }

    // Calculate the average for the window and assign it to smoothed array
    smoothed[i] = sum / count;
  }
}

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

double change_index_to_time(int index, int num_points, double time_span)
{
  return index * (time_span / num_points);
}

// Simulate the time series with noisy Lorentzian peaks and noise added to amplitude
void simulate_time_series(double *series, int num_points, double time_span, LorentzianPeak *peaks, int M, double noise_level)
{
  double t;
  for (int i = 0; i < num_points; i++)
  {
    t = change_index_to_time(i, num_points, time_span);
    series[i] = 0.0;
    for (int m = 0; m < M; m++)
    {
      series[i] += lorentzian(t, peaks[m].width, peaks[m].location);
    }
    series[i] += generate_noise() * noise_level; // Adding noise to amplitude
  }
}

void plot_series_gnuplot(double *series, int size, char * title)
{
  // Open a pipe to Gnuplot
  FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

  if (gnuplotPipe == NULL)
  {
    printf("Error opening pipe to Gnuplot.\n");
    return;
  }

  // Send Gnuplot commands to setup the plot
  fprintf(gnuplotPipe, "set title '%s'\n", title);
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

// Function to find Half Width at Half Maximum (HFHM) for a peak
double find_hwhm(double smoothed[], int length, int peak_index, double peak_height)
{
  double half_height = peak_height / 2.0;
  int left = peak_index, right = peak_index;

  // Find the left point where signal drops below half-height
  while (left > 0 && smoothed[left] > half_height)
  {
    left--;
  }

  // Ensure we don't go out of bounds
  if (smoothed[left] <= half_height && left > 0)
    left++;

  // Find the right point where signal drops below half-height
  while (right < length - 1 && smoothed[right] > half_height)
  {
    right++;
  }

  // Ensure we don't go out of bounds
  if (smoothed[right] <= half_height && right < length - 1)
    right--;

  // Return the width between the two points (FWHM) and divide it by 2 (HWHM)
  return (right - left) / 2.0;
}

// Function to find peaks in the smoothed data and store them in LorentzianPeak structs
int find_peaks(double smoothed[], int length, LorentzianPeak peaks[], int max_peaks, int num_points, double time_span)
{
  int peak_count = 0;

  for (int i = 1; i < length - 1; i++)
  {
    // Check if the current value is greater than its neighbors (local maxima)
    if (smoothed[i] > smoothed[i - 1] && smoothed[i] > smoothed[i + 1])
    {
      // Check if the peak height exceeds a minimum threshold (to avoid noise)
      if (smoothed[i] > MIN_PEAK_HEIGHT)
      {
        if (peak_count < max_peaks)
        {
          // Store the peak information in the LorentzianPeak struct
          peaks[peak_count].location = change_index_to_time(i, num_points, time_span) ;
          peaks[peak_count].amplitude = smoothed[i];
          peaks[peak_count].width = change_index_to_time(find_hwhm(smoothed, length, i, smoothed[i]), num_points, time_span);
          peak_count++;
        }
      }
    }
  }
  return peak_count; // Return the number of peaks found
}

// Function to estimate the average T and a, and their standard deviations
void estimate_avg_and_std(LorentzianPeak *peaks, int num_peaks, double *avg_T, double *std_T, double *avg_a, double *std_a)
{
  double sum_T = 0.0, sum_a = 0.0;
  double *time_intervals = (double *)malloc(num_peaks - 1); // Store the intervals between peaks for calculating standard deviation

  if (time_intervals == NULL)
  {
    fprintf(stderr, "Memory allocation for time_intervals failed\n");
    exit(1);
  }

  // Calculate average T (time between peaks)
  for (int i = 1; i < num_peaks; i++)
  {
    double interval = peaks[i].location - peaks[i - 1].location;
    sum_T += interval;
    time_intervals[i - 1] = interval; // Store the interval
  }
  *avg_T = sum_T / (num_peaks - 1);
  *std_T = compute_std(time_intervals, num_peaks - 1, *avg_T);

  // Calculate average a (width of peaks)
  for (int i = 0; i < num_peaks; i++)
  {
    sum_a += peaks[i].width;
  }
  *avg_a = sum_a / num_peaks;

  double * widthList = (double * ) malloc(num_peaks*sizeof(double));

  if (widthList == NULL)
  {
    fprintf(stderr, "Memory allocation for widthList failed\n");
    exit(1);
  }

  for (int i = 0; i < num_peaks; i++)
  {
    widthList[i] = peaks[i].width;
  }

  

  *std_a = compute_std( widthList, num_peaks, *avg_a);

  free(widthList);

  free(time_intervals);
}

// Function to process time series data, find peaks, and estimate parameters
void process_time_series(double *smoothed, int length, int max_peaks, double actual_T, double actual_a, int num_points, double time_span)
{
  LorentzianPeak *detected_peaks = (LorentzianPeak *)malloc(max_peaks * sizeof(LorentzianPeak));
  int peak_count = find_peaks(smoothed, length, detected_peaks, max_peaks, num_points, time_span);

  double avg_T, std_T, avg_a, std_a;

  // Ensure we have enough peaks to calculate averages
  if (peak_count < 2)
  {
    fprintf(stderr, "Not enough peaks found to estimate parameters.\n");
    return;
  }

  // Estimate average T and a
  estimate_avg_and_std(detected_peaks, peak_count, &avg_T, &std_T, &avg_a, &std_a);

  // Output results
  printf("%f,%f,%f,%f\n", avg_T, avg_a, std_T, std_a);
  free(detected_peaks);
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


  // Array to store the smoothed data
  double smoothed[num_points];

  // Perform the moving average
  moving_average(series, num_points, smoothed);

  double avg_T, avg_a, std_T, std_a;

  // Ensure valid memory allocation
  if (series == NULL || peaks == NULL)
  {
    fprintf(stderr, "Memory allocation failed\n");
    return 1;
  }

  process_time_series(smoothed, num_points, M*2, T, a, num_points, time_span);


  // Free memory
  free(series);
  free(peaks);

  return 0;
}