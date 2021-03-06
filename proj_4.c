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

#define LATTICE_SIDE_X 64
#define LATTICE_SIDE_Y 64

#define BETA_LOWER 0.1
#define BETA_UPPER 1.0
#define BETA_INCR 0.01

#define FERRO_MAG_DOMAINS_BETA 0.4407

#define TYPE_OF_START_HOT "hot"
#define TYPE_OF_START "cold"

#define ARRAY_SIZE 90

double boltzmann_const = 1.38064852e-23;
// critical temperature for phase change
double beta_sub_c = 0.4407;

//
// create 1 D array of int in memory
//
static int * create_1d_array_i(int cols)
{
    int * p1dArray = malloc(cols*sizeof(int));
    return p1dArray;
}

//
// create 1 D array of double in memory
//
static double * create_1d_array_d(int cols)
{
    double * p1dArray = malloc(cols*sizeof(double));
    return p1dArray;
}

//
// Destroy a created 1D array of int
//
static void destroy_1d_array_i(int p1dArray[])
{
    free(p1dArray);
}

//
// Destroy a created 1D array of double
//
static void destroy_1d_array_d(double p1dArray[])
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
static void destroy_2d_array_c(char *array[], int rows)
{
    int i;
    for(i = 0; i < rows; i++)
    {
        free(array[i]);
    }
    free(array);
}

static void print_array_c(char ** array, int rows, int cols)
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

static void display_scatter_arr_char(char **array, int rows, int cols, int pos_incr, char *heading, char *x_label, char *y_label)
{
    static float fxp[1];
    static float fyp[1];
    int j, k, char_code, char_colour;

    // set starting position
    int starting_pos = pos_incr * rows;
    int x_pos;
    int y_pos = starting_pos;

    float fxmin = 0.0;
    float fxmax = (float) pos_incr * (cols + 1);
    float fymin = 0.0;
    float fymax = (float) pos_incr * (rows + 1);

    cpgbbuf();
    cpgsci(14);
    cpgenv(fxmin, fxmax, fymin, fymax, 0, 1);

    for(j=0; j < rows; j++) {
        x_pos = pos_incr;
        fyp[0] = (float)y_pos;

        for(k=0; k < cols; k++) {
            if (array[j][k] == 'u') {
                char_code = 30;
                char_colour = 8;
            } else {
                char_code = 31;
                char_colour = 10;
            }
            fxp[0] = (float)x_pos;

            /*
             * Now make scatter plot with correct symbol and colour
             */
            cpgsci(char_colour);
            cpgpt(1, fxp, fyp, char_code);
            x_pos = x_pos + pos_incr;
        }
        y_pos = y_pos - pos_incr;
    }

    cpgsci(14);
    cpglab(x_label, y_label,heading);

    cpgebuf();

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

static int get_spin_val_from_char(char this_spin_char)
{
    int val = this_spin_char == 'u' ? 1 : -1;
//    printf("Val: %d\n", val);
    return val;
}

static double calc_delta_e(char **array, int x_pos, int y_pos, int x_size, int y_size)
{
    double val;

    char this_spin_char = query_array(array, x_pos, y_pos, x_size, y_size);

    // flip spin
    if (this_spin_char == 'u') {
        this_spin_char = 'd';
    } else {
        this_spin_char = 'u';
    }

    char x_minus_one_neigh = query_array(array, x_pos - 1, y_pos, x_size, y_size);
    char x_plus_one_neigh = query_array(array, x_pos + 1, y_pos, x_size, y_size);
    char y_minus_one_neigh = query_array(array, x_pos, y_pos - 1, x_size, y_size);
    char y_plus_one_neigh = query_array(array, x_pos, y_pos + 1, x_size, y_size);
    val = -2.0 * get_spin_val_from_char(this_spin_char) * (get_spin_val_from_char(x_minus_one_neigh) + get_spin_val_from_char(x_plus_one_neigh) + get_spin_val_from_char(y_minus_one_neigh) + get_spin_val_from_char(y_plus_one_neigh));
    return val;
}

static void sweep_2d_array(char **array, int x_size, int y_size, int *knuth_arr, int knuth_size, double beta)
{
    int x_pos, y_pos, k_count, this_k_val;
    double this_delta_e;

    for (k_count = 0; k_count < knuth_size; k_count++) {
        this_k_val = knuth_arr[k_count];
        x_pos = this_k_val % x_size;
        y_pos = (int)floor((double)(this_k_val / x_size));

        this_delta_e = calc_delta_e(array, x_pos, y_pos, x_size, y_size);

//        printf("This delta_e: %lf\n", this_delta_e);

        if (this_delta_e < 0.0 || exp(-beta * this_delta_e) > drand48()) {
            // flip spin
            if (array[x_pos][y_pos] == 'u') {
                array[x_pos][y_pos] = 'd';
            } else {
                array[x_pos][y_pos] = 'u';
            }
        }
    }
}

static void fill_2d_array_c(char ** array, int x_size, int y_size, char * start)
{
    int j, k, boundary;
    // set default value for cold start
    char spin;

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

//            if (j == x_size - 2)
//            {
//                spin = array[1][k];
//                boundary = 1;
//            }

            if (k == y_size - 1)
            {
                spin = array[j][0];
                boundary = 1;
            }

//            if (k == y_size - 2)
//            {
//                spin = array[j][1];
//                boundary = 1;
//            }

            if (boundary == 0) {
                // do this if we're not at a boundary

                if (start == TYPE_OF_START_HOT) {
                    // randomise
                    if (drand48() > 0.5) {
                        spin = 'u';
                    } else {
                        spin = 'd';
                    }
                } else {
                    spin = 'u';
                }
            }

            array[j][k] = spin;
        }
    }
}

static double calc_mag_site(char **array, int x_pos, int y_pos)
{
    if (array[x_pos][y_pos] == 'u') {
        return 1.0;
    } else {
        return -1.0;
    }
}

static double calc_mag_of_lattice(char **array, int x_size, int y_size)
{
    double mag_sum = 0.0;
    int j, k;

    for (j = 0; j < x_size; j++) {
        for (k = 0; k < y_size; k++) {
            mag_sum = mag_sum + calc_mag_site(array, j, k);
        }
    }

    return fabs(mag_sum) / (x_size * y_size);
}

static double calc_variance_of_lattice(double *mag_arr_across_sweeps, int no_of_measurements, double magnetisation_mean, int lattice_size)
{
    double var_sum = 0.0;
    int i;

    for (i = 0; i < no_of_measurements; i++) {
        var_sum = var_sum + pow(mag_arr_across_sweeps[i] - magnetisation_mean, 2.0);
    }

    return sqrt(var_sum / (lattice_size - 1));
}

static void fill_knuth_1d(int *knuth_arr, int size)
{
    int i;

    for (i = 0; i < size; i++)
    {
        knuth_arr[i] = i + 1;
    }
}

static void print_array_knuth(int * array, int size)
{
    int i;
    for(i=0; i < size; i++)
    {
        printf("%d ", array[i]);
    }

    printf("\n");
}

static void swap_elements_one_d_arr(int one_d_array[], int pos_1, int pos_2)
{
    int temp = one_d_array[pos_1];
    one_d_array[pos_1] = one_d_array[pos_2];
    one_d_array[pos_2] = temp;
}

static void knuth(int one_d_array[], int size)
{
    int i, spot, temp;

    for (i = 0; i < size; i++) {
        spot = (int)(drand48() * (double)size);

        temp = one_d_array[i];
        one_d_array[i] = one_d_array[spot];
        one_d_array[spot] = temp;
    }
}

static void disp_line_spec_axis(int num, double x_vals[], double y_vals[], double std_dev_vals[], float x_min, float x_max, float y_min, float y_max, char *heading, char *x_label, char *y_label)
{
    int i, j;
    double largest_std_dev = 0.0;
    double beta_for_lgst_std_dev = 0.0;
    double mag_for_lgst_std_dev = 0.0;
    static float f_x_vals[ARRAY_SIZE];
    static float f_y_vals[ARRAY_SIZE];
    static float f_err_vals[ARRAY_SIZE];
    static float f_crit_temp_beta_x_vals[ARRAY_SIZE];
    static float f_crit_temp_beta_y_vals[ARRAY_SIZE];
    static float f_crit_temp_mag_x_vals[ARRAY_SIZE];
    static float f_crit_temp_mag_y_vals[ARRAY_SIZE];

    static float f_crit_temp_beta_pt_x_val[1];
    static float f_crit_temp_beta_pt_y_val[1];

    for (i = 0; i < num; i++)
    {
        f_x_vals[i] = (float) x_vals[i];
        f_y_vals[i] = (float) y_vals[i];
        f_err_vals[i] = (float) std_dev_vals[i];

        if (std_dev_vals[i] > largest_std_dev) {
            largest_std_dev = std_dev_vals[i];
            beta_for_lgst_std_dev = x_vals[i];
            mag_for_lgst_std_dev = y_vals[i];
        }
    }

    // create array for critical temp
    for (j = 0; j < num; j++) {
        f_crit_temp_beta_x_vals[j] = (float) beta_for_lgst_std_dev;
        f_crit_temp_beta_y_vals[j] = j;
        f_crit_temp_mag_x_vals[j] = j;
        f_crit_temp_mag_y_vals[j] = (float) mag_for_lgst_std_dev;
    }

    f_crit_temp_beta_pt_x_val[1] = (float) beta_for_lgst_std_dev;
    f_crit_temp_beta_pt_y_val[1] = 0.0;

    cpgbbuf();

//    printf("beta_half: %lf\n", largest_std_dev);

    cpgsci(15);
    cpgenv(x_min, x_max, y_min, y_max, 0, 2);

    cpgsci(8);
    // main line of plot
    cpgline(num, f_x_vals, f_y_vals);

    // print line for crit temp beta
    cpgsci(4);
    cpgline(num, f_crit_temp_beta_x_vals, f_crit_temp_beta_y_vals);

    // print line for crit temp mag
    cpgsci(4);
    cpgline(num, f_crit_temp_mag_x_vals, f_crit_temp_mag_y_vals);

    // print point for crit temp
    cpgsci(5);
    cpgpt(1, f_crit_temp_beta_pt_x_val, f_crit_temp_beta_pt_y_val, -8);

    // print label for crit temp beta point
    char pt_label[256];
    snprintf(pt_label, sizeof(pt_label), "Critical temp: %lf", beta_for_lgst_std_dev);
    cpgsci(6);
    cpgtext((float) (beta_for_lgst_std_dev + 0.05), 0.05, pt_label);

    cpgsci(10);
    cpgerrb(6, num, f_x_vals, f_y_vals, f_err_vals, 1.0);

    cpgsci(1);
    cpglab(x_label, y_label, heading);

    cpgebuf();
}

static void disp_ferro_mag_domains_scttr_plt(char ** array, int x_size, int y_size, char * heading, char *x_label, char * y_label)
{
    static float f_xp[LATTICE_SIDE_X * LATTICE_SIDE_Y];
    static float f_yp[LATTICE_SIDE_X * LATTICE_SIDE_Y];

    int j, k;
    int array_counter = 0;

    for (j = 0; j < x_size; j++) {
        for (k = 0; k < y_size; k++) {
            if (array[j][k] == 'u') {
                f_xp[array_counter] = (float) j;
                f_yp[array_counter] = (float) k;
                array_counter++;
            }
        }
    }

    float fxmin, fxmax, fymin, fymax;
    fxmin = 0.0;
    fxmax = (float) x_size;
    fymin = 0.0;
    fymax = (float) y_size;

    cpgbbuf();

    cpgsci(14);
    cpgenv(fxmin, fxmax, fymin, fymax, 0, 1);

    cpgsci(1);
    cpgpt(array_counter, f_xp, f_yp, -1);

    cpgsci(14);
    cpglab(x_label, y_label, heading);

    cpgebuf();
}

int main(void)
{
    int j, k;
    char spin;
    srand48((int) time(NULL));

    int i;
    int no_of_beta_vals = (int) ((BETA_UPPER + BETA_INCR - BETA_LOWER) / BETA_INCR);
    int lattice_size, sweep_counter, number_of_measurements_per_beta, mag_measurement_index_per_beta;
    int count_beta_vals;
    double beta, magnetisation_mean, this_magnetisation, energy, variance, mag_mean_sum_per_beta;
    char plot_heading_line_plot[256];
    char x_label_line_plot[256];
    char y_label_line_plot[256];

//    if (cpgbeg(0, "?", 1, 1) != 1) {
    if (cpgbeg(0, "/XWINDOW", 1, 1) != 1) {
//    if (cpgbeg(0, "/media/sf_Comp/code/proj_4/proj_4_plot.ps/CPS", 1, 1) != 1) {
//    if (cpgbeg(0, "proj_4_plot.ps/CPS", 1, 1) != 1) {
//    if (cpgbeg(0, "/PS", 1, 1) != 1) {
        exit(EXIT_FAILURE);
    }
    cpgask(1);

    double *beta_arr_for_beta_val = create_1d_array_d(no_of_beta_vals);
    double *mag_arr_for_beta_val = create_1d_array_d(no_of_beta_vals);
    double *variance_arr_for_beta_val = create_1d_array_d(no_of_beta_vals);
    number_of_measurements_per_beta = (MAX_SWEEPS - 100) / MEASUREMENT_SWEEP_SEPARATION;
    double *mag_arr_per_beta = create_1d_array_d(number_of_measurements_per_beta);

    // Allocate a region of memory of size LATTICE_SIDE_X*LATTICE_SIDE_Y*sizeof(double) and return a pointer to it
    char **Array = create_2d_array_c(LATTICE_SIDE_X, LATTICE_SIDE_Y);
    fill_2d_array_c(Array, LATTICE_SIDE_X, LATTICE_SIDE_Y, TYPE_OF_START);

//        print_array_c(Array, LATTICE_SIDE_X, LATTICE_SIDE_Y);

    lattice_size = LATTICE_SIDE_X * LATTICE_SIDE_Y;
    int *knuth_arr = create_1d_array_i(lattice_size);

    //    display_scatter_arr_char(Array, LATTICE_SIDE_X, LATTICE_SIDE_Y, 1, "Plotting a 2-D array of char representing spins", "x-axis", "y-axis");

    fill_knuth_1d(knuth_arr, lattice_size);

    count_beta_vals = 0;

    for (beta = BETA_LOWER; beta <= BETA_UPPER; beta = beta + BETA_INCR) {
        beta_arr_for_beta_val[count_beta_vals] = beta;
        mag_measurement_index_per_beta = 0;
        mag_mean_sum_per_beta = 0.0;

        for (sweep_counter = 0; sweep_counter < MAX_SWEEPS; sweep_counter++) {
            knuth(knuth_arr, lattice_size);
            //            print_array_knuth(knuth_arr, knuth_arr_size);

            sweep_2d_array(Array, LATTICE_SIDE_X, LATTICE_SIDE_Y, knuth_arr, lattice_size, beta);

            // allow 100 sweeps before measuring and only measure magnetisation every so often
            if (sweep_counter > 100 && sweep_counter % MEASUREMENT_SWEEP_SEPARATION == 0) {
                // take a measurement
                this_magnetisation = calc_mag_of_lattice(Array, LATTICE_SIDE_X, LATTICE_SIDE_Y);
                mag_arr_per_beta[mag_measurement_index_per_beta] = this_magnetisation;
                mag_mean_sum_per_beta = mag_mean_sum_per_beta + this_magnetisation;

                mag_measurement_index_per_beta++;
            }
        }

        magnetisation_mean = mag_mean_sum_per_beta / (mag_measurement_index_per_beta + 1);
        variance = calc_variance_of_lattice(mag_arr_per_beta, mag_measurement_index_per_beta, magnetisation_mean, lattice_size);
        variance_arr_for_beta_val[count_beta_vals] = variance;
        mag_arr_for_beta_val[count_beta_vals] = magnetisation_mean;

        count_beta_vals++;

        if (beta >= FERRO_MAG_DOMAINS_BETA && beta < FERRO_MAG_DOMAINS_BETA + BETA_INCR) {
            char plot_heading_mag_dmns[256];
            snprintf(plot_heading_mag_dmns, sizeof(plot_heading_mag_dmns), "Ferro-magnetic domains of spin-up lattice sites - side x: %d, y: %d, beta: %.4lf, start: %s", LATTICE_SIDE_X, LATTICE_SIDE_Y, beta, TYPE_OF_START);
            disp_ferro_mag_domains_scttr_plt(Array, LATTICE_SIDE_X, LATTICE_SIDE_Y, plot_heading_mag_dmns, "Lattice side - x direction", "Lattice side - y direction");
        }

        if (LATTICE_SIDE_Y <= 50 && LATTICE_SIDE_Y <= 50 && beta == BETA_LOWER) {
            // show initial plot of spins
            char plot_heading_spins_plot[256];
            snprintf(plot_heading_spins_plot, sizeof(plot_heading_spins_plot), "Plot of 2D character array representing lattice spins, lattice sides x: %d, y: %d, beta: %.4lf", LATTICE_SIDE_X, LATTICE_SIDE_Y, beta);
            display_scatter_arr_char(Array, LATTICE_SIDE_X, LATTICE_SIDE_Y, 1, plot_heading_spins_plot, "Lattice side - x direction", "Lattice side - y direction");
        }
    }

    snprintf(plot_heading_line_plot, sizeof(plot_heading_line_plot), "Magnetisation versus beta - lattice side x: %d, y: %d, type of start: %s", LATTICE_SIDE_X, LATTICE_SIDE_Y, TYPE_OF_START);
    snprintf(x_label_line_plot, sizeof(x_label_line_plot), "Beta (1/J) - from %.1lf to %.1lf with %.2lf increments", BETA_LOWER, BETA_UPPER, BETA_INCR);
    snprintf(y_label_line_plot, sizeof(y_label_line_plot), "Magnetisation (dimensionless)");

    disp_line_spec_axis(count_beta_vals, beta_arr_for_beta_val, mag_arr_for_beta_val, variance_arr_for_beta_val, 0.1, 1.0, 0.0, 1.0, plot_heading_line_plot, x_label_line_plot, y_label_line_plot);


    destroy_1d_array_i(knuth_arr);
    destroy_2d_array_c(Array, LATTICE_SIDE_X); // deallocate the memory

    destroy_1d_array_d(beta_arr_for_beta_val);
    destroy_1d_array_d(mag_arr_for_beta_val);
    destroy_1d_array_d(variance_arr_for_beta_val);
    destroy_1d_array_d(mag_arr_per_beta);

    cpgend();

    printf("Critical temperature for phase change, T_sub_c: %.6E\n", 1.0 / boltzmann_const * beta_sub_c);
}