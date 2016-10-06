#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "cpgplot.h"

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

#define MAX_SWEEPS 1100
#define MEASUREMENT_SWEEP_SEPARATION 10

#define MIN_LATTICE_SIDE 15
#define MAX_LATTICE_SIDE 25
#define LATTICE_SIDE_INCR 5

#define LOWER_LIMIT_OF_BETA_VALUES 0.1
#define UPPER_LIMIT_OF_BETA_VALUES 1.0
#define BETA_INCR 0.01

#define ARRAY_SIZE 10000

double boltzmann_const = 1.38064852e-23;
// critical temperature for phase change
double beta_sub_c = 0.4407;

//
// create 1 D array of char in memory
//
char * create_1d_array_c(int cols)
{
    char * p1dArray = malloc(cols*sizeof(char));
    return p1dArray;
}

//
// create 1 D array of int in memory
//
int * create_1d_array_i(int cols)
{
    int * p1dArray = malloc(cols*sizeof(int));
    return p1dArray;
}

//
// create 1 D array of double in memory
//
double * create_1d_array_d(int cols)
{
    double * p1dArray = malloc(cols*sizeof(double));
    return p1dArray;
}

//
// Destroy a created 1D array of int
//
void destroy_1d_array_i(int p1dArray[])
{
    free(p1dArray);
}

//
// Destroy a created 1D array of char
//
void destroy_1d_array_c(char p1dArray[])
{
    free(p1dArray);
}

//
// Destroy a created 1D array of double
//
void destroy_1d_array_d(double p1dArray[])
{
    free(p1dArray);
}

//
// Create an array of dimensions (rows,cols) of char in memory
//
char ** create_2d_array_c(int rows, int cols)
{
    char ** p2dArray;
    p2dArray = (char **) malloc(rows*sizeof(char *));
    int i = 0;
    for(i=0; i< rows; i++)
    {
        p2dArray[i] = (char *) malloc(cols*sizeof(char));
    }
    return p2dArray;
}


//
// Destroy a 2D array of dimensions (rows,[]) of char in memory
//
void destroy_2d_array_c(char *array[], int rows)
{
    int i=0;
    for(i=0;i<rows;i++)
    {
        free(array[i]);
    }
    free(array);
}

//
// Create an array of dimensions (rows,cols) of int in memory
//
int ** create_2d_array_i(int rows, int cols)
{
    int ** p2dArray;
    p2dArray = (int **) malloc(rows*sizeof(char *));
    int i = 0;
    for(i=0; i< rows; i++)
    {
        p2dArray[i] = (int *) malloc(cols*sizeof(int));
    }
    return p2dArray;
}


//
// Destroy a 2D array of dimensions (rows,[]) of int in memory
//
void destroy_2d_array_i(int *array[], int rows)
{
    int i=0;
    for(i=0;i<rows;i++)
    {
        free(array[i]);
    }
    free(array);
}

void print_array_c(char ** array, int rows, int cols)
{
    int i,j;
    for(i=0;i<rows;i++)
    {
        for(j=0;j<cols;j++)
        {
//            printf("x: %d, y: %d\n", i, j);
            printf(" %c", array[i][j]);
        }
        printf("\n");
    }
}

void display_scatter_arr_char(char **array, int rows, int cols, int pos_incr, char *heading, char *x_label,
                              char *y_label)
{
    static float fxp[1];
    static float fyp[1];
    int j, k, char_code;

    // set starting position
    int starting_pos = pos_incr * rows;
    int x_pos;
    int y_pos = starting_pos;

    float fxmin = 0.0;
    float fxmax = (float) pos_incr * (cols + 1);
    float fymin = 0.0;
    float fymax = (float) pos_incr * (rows + 1);

    //
    // Set up the display region
    //
    cpgenv(fxmin, fxmax, fymin, fymax, 0, 1);
    cpgbbuf();

    for(j=0; j < rows; j++) {
        x_pos = pos_incr;
        fyp[0] = (float)y_pos;

        for(k=0; k < cols; k++) {
            if (array[j][k] == 'u') {
                char_code = 30;
            } else {
                char_code = 31;
            }
            fxp[0] = (float)x_pos;

            /*
             * Now make scatter plot with correct symbol
             */
            cpgpt(1, fxp, fyp, char_code);
            x_pos = x_pos + pos_incr;
        }
        y_pos = y_pos - pos_incr;
    }


    // Label the plot

    cpglab(x_label, y_label,heading);
    // cpgsave saves the current graphics
    cpgsave();

}

char query_array(char **array, int x_pos, int y_pos, int x_size, int y_size)
{
    // first deal with out-of-bounds indices
    if (x_pos > x_size - 1) {
        x_pos = 0;
    }

    if (x_pos < 0) {
        x_pos = x_size - 1;
    }

    if (y_pos > y_size - 1) {
        y_pos = 0;
    }

    if (y_pos < 0) {
        y_pos = y_size - 1;
    }

    return array[x_pos][y_pos];
}

int get_spin_val_from_char(char this_spin_char)
{
    int val = this_spin_char == 'u' ? 1 : -1;
//    printf("Val: %d\n", val);
    return val;
}

double calc_delta_e(char **array, int x_pos, int y_pos, int x_size, int y_size)
{
    double val;

    char this_spin_char = query_array(array, x_pos, y_pos, x_size, y_size);
    char x_minus_one_neigh = query_array(array, x_pos - 1, y_pos, x_size, y_size);
    char x_plus_one_neigh = query_array(array, x_pos + 1, y_pos, x_size, y_size);
    char y_minus_one_neigh = query_array(array, x_pos, y_pos - 1, x_size, y_size);
    char y_plus_one_neigh = query_array(array, x_pos, y_pos + 1, x_size, y_size);
    val = -2.0 * get_spin_val_from_char(this_spin_char) * (get_spin_val_from_char(x_minus_one_neigh) + get_spin_val_from_char(x_plus_one_neigh) + get_spin_val_from_char(y_minus_one_neigh) + get_spin_val_from_char(y_plus_one_neigh));
    return val;
}

void sweep_2d_array(char **array, int x_size, int y_size, int *knuth_arr, int knuth_size, double beta)
{
    int x_pos, y_pos, k_count, this_k_val;
    double this_delta_e;

    for (k_count = 0; k_count < knuth_size; k_count++) {
        this_k_val = knuth_arr[k_count];
        x_pos = this_k_val % x_size;
        y_pos = (int)floor((double)(this_k_val / x_size));

        // flip spin
        if (array[x_pos][y_pos] == 'u') {
            array[x_pos][y_pos] = 'd';
        } else {
            array[x_pos][y_pos] = 'u';
        }

        this_delta_e = calc_delta_e(array, x_pos, y_pos, x_size, y_size);

//        printf("This delta_e: %lf\n", this_delta_e);

        if (!(this_delta_e < 0.0 || exp(-beta * this_delta_e) > drand48())) {
            // flip spin back
            if (array[x_pos][y_pos] == 'u') {
                array[x_pos][y_pos] = 'd';
            } else {
                array[x_pos][y_pos] = 'u';
            }
        }
    }
}

void fill_2d_array_c(char ** array, int x_size, int y_size, char * start)
{
    int j, k, boundary;
    // set default value for cold start
    char spin = 'u';

    for (j = 0; j < x_size; j++)
    {

        for(k = 0; k < y_size; k++)
        {
            boundary = 0;

            if (j == x_size - 1)
            {
                spin = array[0][k];
                boundary = 1;
            }

            if (j == x_size - 2)
            {
                spin = array[1][k];
                boundary = 1;
            }

            if (k == y_size - 1)
            {
                spin = array[j][0];
                boundary = 1;
            }

            if (k == y_size - 2)
            {
                spin = array[j][1];
                boundary = 1;
            }

            if (boundary == 0) {
                // do this if we're not at a boundary
                if (start = "hot") {
                    // randomise
                    if (drand48() > 0.5) {
                        spin = 'u';
                    } else {
                        spin = 'd';
                    }
                }
            }

            array[j][k] = spin;
        }
    }
}

double calc_mag_site(char **array, int x_pos, int y_pos)
{
    if (array[x_pos][y_pos] == 'u') {
        return 1.0;
    } else {
        return 0.0;
    }
}

double calc_mag_of_lattice(char **array, int x_size, int y_size)
{
    double mag_sum = 0.0;
    int j, k;

    for (j = 0; j < x_size; j++) {
        for (k = 0; k < y_size; k++) {
            mag_sum = mag_sum + calc_mag_site(array, j, k);
        }
    }

    return mag_sum / (x_size * y_size);
}

double calc_variance_of_lattice(double *mag_arr_across_sweeps, int no_of_measurements, double magnetisation_mean,
                                int lattice_size)
{
    double var_sum = 0.0;
    int i;

    for (i = 0; i < no_of_measurements; i++) {
        var_sum = var_sum + pow(mag_arr_across_sweeps[i] - magnetisation_mean, 2.0);
    }

    return var_sum / (lattice_size - 1);
}

void fill_knuth_1d(int *knuth_arr, int size, int *counter)
{
    int i;

    for (i = 0; i < size; i++)
    {
        knuth_arr[i] = *counter;
        * counter = * counter + 1;
    }
}

void print_array_knuth(int * array, int size)
{
    int i;
    for(i=0; i < size; i++)
    {
        printf("%d ", array[i]);
    }

    printf("\n");
}

void swap_elements_one_d_arr(int one_d_array[], int pos_1, int pos_2)
{
    int temp = one_d_array[pos_1];
    one_d_array[pos_1] = one_d_array[pos_2];
    one_d_array[pos_2] = temp;
}

void knuth(int one_d_array[], int size)
{
    int i, j;
    double boundary;
    int found;

    for (i = 0; i < size; i++) {
        boundary = (double)i + drand48() * (double)(size - i - 1);
        found = 0;
        for (j = i + 1; j < size && 0 == found; j++) {
            if (boundary < (double)j) {
//                printf("knuth i: %d, j: %d\n", i, j);
                swap_elements_one_d_arr(one_d_array, i, j);
                found = 1;
            }
        }
    }
}

double calc_energy_of_lattice(char ** Array, int x_size, int y_size)
{
    double energy = 0.0;
    int j, k;

    for (j = 0; j < x_size; j++) {
        for (k = 0; k < y_size; k++) {
            energy = energy + calc_delta_e(Array, j, k, x_size, y_size);
        }
    }

    return energy;
}

void disp_line(int num, double x_vals[], double y_vals[], char * heading, char * x_label, char * y_label)
{
    int j = 0;
    static float f_x_vals[ARRAY_SIZE];
    static float f_y_vals[ARRAY_SIZE];
    float x_min = 1.0e30;
    float x_max = -1.0e30;
    float y_min = 1.0e30;
    float y_max = -1.0e30;
    for(j = 0; j < num; j++)
    {
        f_x_vals[j] = x_vals[j];
        f_y_vals[j] = y_vals[j];
        if (f_x_vals[j] < x_min)
        {
            x_min = f_x_vals[j];
        }
        if (f_x_vals[j] > x_max)
        {
            x_max = f_x_vals[j];
        }
        if (f_y_vals[j] < y_min)
        {
            y_min = f_y_vals[j];
        }
        if (f_y_vals[j] > y_max)
        {
            y_max = f_y_vals[j];
        }
    }

    // start buffering
    cpgbbuf();

    cpgenv(x_min, x_max, y_min, y_max, 0, 2);
    /*
     * Now plot a histogram
     */

    // for outputting to ps file
    cpgsci(1);

    cpgline(num, f_x_vals, f_y_vals);

//    printf("v_min: %f, v_max: %f", v_min, v_max);
    cpglab(x_label, y_label, heading);

//    cpgsave();

//    cpgpage();

    // end buffering
    cpgebuf();
}

void disp_line_spec_axis(int num, double x_vals[], double y_vals[], double err_vals[], float x_min, float x_max, float y_min, float y_max, char *heading, char *x_label, char *y_label)
{
    int j = 0;
    static float f_x_vals[ARRAY_SIZE];
    static float f_y_vals[ARRAY_SIZE];
    static float f_err_vals[ARRAY_SIZE];
    for(j = 0; j < num; j++)
    {
        f_x_vals[j] = x_vals[j];
        f_y_vals[j] = y_vals[j];
        f_err_vals[j] = err_vals[j];
    }
    cpgbbuf();
    /*
     * Now plot a histogram
     */

    cpgenv(x_min, x_max, y_min, y_max, 0, 2);

    // for outputting to ps file
    cpgsci(1);

    cpgline(num, f_x_vals, f_y_vals);

    cpgerrb(6, num, f_x_vals, f_y_vals, f_err_vals, 1.0);

//    printf("v_min: %f, v_max: %f", v_min, v_max);
    cpglab(x_label, y_label, heading);

    cpgebuf();
}

void disp_err_bar(int num, double x_vals[], double y_vals[], float x_min, float x_max, float y_min, float y_max, char *heading, char *x_label, char *y_label)
{
    int j = 0;
    static float f_x_vals[ARRAY_SIZE];
    static float f_y_vals[ARRAY_SIZE];
    for(j = 0; j < num; j++)
    {
        f_x_vals[j] = x_vals[j];
        f_y_vals[j] = y_vals[j];
    }
    cpgbbuf();
    /*
     * Now plot a histogram
     */

    cpgenv(x_min, x_max, y_min, y_max, 0, 2);

    // for outputting to ps file
    cpgsci(1);

    cpgline(num, f_x_vals, f_y_vals);

//    printf("v_min: %f, v_max: %f", v_min, v_max);
    cpglab(x_label, y_label, heading);

    cpgebuf();
}

int main(void)
{
    int j, k, curr_latt_side;
    char spin;
    srand48((int) time(NULL));

    int i;
    int no_of_beta_vals = (int) ((UPPER_LIMIT_OF_BETA_VALUES - LOWER_LIMIT_OF_BETA_VALUES) / BETA_INCR);
    int lattice_size;
    int knuth_arr_size, sweep_counter, beta_int;
    double beta, magnetisation_mean, magnetisation_within_sweep, energy, variance;
    char plot_heading[256];

    //    if (cpgbeg(0, "?", 1, 1) != 1) {
    if (cpgbeg(0, "/XWINDOW", 1, 1) != 1) {
//        if (cpgbeg(0, "proj_4_plot.ps/CPS", 1, 1) != 1) {
        exit(EXIT_FAILURE);
    }
    cpgask(1);

    for (curr_latt_side = MIN_LATTICE_SIDE; curr_latt_side <= MAX_LATTICE_SIDE; curr_latt_side = curr_latt_side + LATTICE_SIDE_INCR) {
        // Allocate a region of memory of size curr_latt_side*curr_latt_side*sizeof(double) and return a pointer to it
        char **Array = create_2d_array_c(curr_latt_side, curr_latt_side);
        fill_2d_array_c(Array, curr_latt_side, curr_latt_side, "cold");

//        print_array_c(Array, curr_latt_side, curr_latt_side);

        lattice_size = curr_latt_side * curr_latt_side;
        int *knuth_arr = create_1d_array_i(lattice_size);
        double *beta_arr_for_beta_val = create_1d_array_d(no_of_beta_vals);
        double *mag_arr_for_beta_val = create_1d_array_d(no_of_beta_vals);
        double *variance_arr_for_beta_val = create_1d_array_d(no_of_beta_vals);
        int number_of_sweeps = (MAX_SWEEPS - 100) / MEASUREMENT_SWEEP_SEPARATION;
        double *mag_arr_across_sweeps = create_1d_array_d(number_of_sweeps);
        int count_mag_across_sweeps;
        double mag_median_sum_across_sweeps = 0.0;

        //    display_scatter_arr_char(Array, curr_latt_side, curr_latt_side, 1, "Plotting a 2-D array of char representing spins", "x-axis", "y-axis");

        knuth_arr_size = 0;
        fill_knuth_1d(knuth_arr, lattice_size, &knuth_arr_size);

        for (beta_int = 1; beta_int < no_of_beta_vals; beta_int++) {
            beta = BETA_INCR * beta_int;
            beta_arr_for_beta_val[beta_int] = beta;
            for (sweep_counter = 0; sweep_counter < MAX_SWEEPS; sweep_counter++) {
                knuth(knuth_arr, knuth_arr_size);
                //            print_array_knuth(knuth_arr, knuth_arr_size);

                sweep_2d_array(Array, curr_latt_side, curr_latt_side, knuth_arr, knuth_arr_size, beta);
                //        printf("Flipped element of array at x = %d, y = %d: %c\n", x_rand_i, y_rand_i, query_array(Array, x_rand_i, y_rand_i, curr_latt_side, curr_latt_side));

                count_mag_across_sweeps = 0;

                if (sweep_counter > 100 && sweep_counter % MEASUREMENT_SWEEP_SEPARATION == 0) {
                    // take a measurement
//                    energy = calc_energy_of_lattice(Array, curr_latt_side, curr_latt_side);
                    magnetisation_within_sweep = calc_mag_of_lattice(Array, curr_latt_side, curr_latt_side);
                    mag_arr_across_sweeps[count_mag_across_sweeps] = magnetisation_within_sweep;
                    count_mag_across_sweeps++;

                    mag_median_sum_across_sweeps = mag_median_sum_across_sweeps + magnetisation_within_sweep;
                }
            }

            magnetisation_mean = mag_median_sum_across_sweeps / count_mag_across_sweeps;
            variance = calc_variance_of_lattice(mag_arr_across_sweeps, count_mag_across_sweeps, magnetisation_mean, lattice_size);
            variance_arr_for_beta_val[beta_int] = variance;
            mag_arr_for_beta_val[beta_int] = magnetisation_mean;
        }

        //    display_scatter_arr_char(Array, curr_latt_side, curr_latt_side, 1, "Plotting a 2-D array of char representing spins", "x-axis", "y-axis");

        snprintf(plot_heading, sizeof(plot_heading), "Magnetisation versus beta - lattice side: %d", curr_latt_side);
//        disp_line(beta_int, beta_arr_for_beta_val, mag_arr_for_beta_val, plot_heading, "Beta", "Magnetisation");


        int test_counter;
        for (test_counter = 0; test_counter < beta_int; test_counter++) {
            printf("count_beta_vals: %d, beta: %lf, mag: %lf\n", test_counter, beta_arr_for_beta_val[test_counter], mag_arr_for_beta_val[test_counter]);
        }

//        disp_line_spec_axis(beta_int, beta_arr_for_beta_val, mag_arr_for_beta_val, variance_arr_for_beta_val, 0.0, 1.0, 0.0, 1.0, plot_heading, "Beta", "Magnetisation");

        destroy_1d_array_i(knuth_arr);
        destroy_2d_array_c(Array, curr_latt_side); // deallocate the memory
        destroy_1d_array_d(beta_arr_for_beta_val);
        destroy_1d_array_d(mag_arr_for_beta_val);
    }

    cpgend();
//    printf("Critical temperature for phase change, T_sub_c: %.6E\n", 1.0 / boltzmann_const * beta_sub_c);
}