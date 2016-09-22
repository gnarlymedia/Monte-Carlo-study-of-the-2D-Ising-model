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

#define MAX_SIZE 10000000
#define ARRAY_SIZE 10000000

//
// create 1 D array of char in memory
//
char * create_1d_array_c(int cols)
{
    char * p1dArray = malloc(cols*sizeof(char));
    return p1dArray;
}
//
// Destroy a created 1D array of char
//
void destroy_1d_array_c(char p1dArray[])
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

void print_array_c(char ** array, int rows, int cols)
{
    int i,j;
    for(i=0;i<rows;i++)
    {
        for(j=0;j<cols;j++)
        {
            printf(" %c", array[i][j]);
        }
        printf("\n");
    }
}

void display_scatter_char(char ** array, int rows, int cols, int pos_incr, char * heading, char *x_label, char * y_label)
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
            if ('u' == array[j][k]) {
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

char interr_array(char ** array, int x_pos, int y_pos)
{
    return array[x_pos][y_pos];
}

int main(void)
{
    int j, k;
    int xsize = 25;
    int ysize = 50;
    char spin;
    srand48((int) time(NULL));

    // Allocate a region of memory of size xsize*ysize*sizeof(double) and return a pointer to it
    char ** Array = create_2d_array_c(xsize, ysize);
    for(j=0; j < xsize; j++)
    {
        for(k=0; k < ysize; k++)
        {
            if (drand48() > 0.5) {
                spin = 'u';
            } else {
                spin = 'd';
            }
            Array[j][k] = spin;
        }
    }

//    print_array_c(Array, xsize, ysize);

    // Initiate cpgplot and open output device

//    if (1 != cpgbeg(0, "?", 1, 1))
    if (1 != cpgbeg(0, "/XWINDOW", 1, 1))
//    if (1 != cpgbeg(0, "proj_3_plot.ps/CPS", 1, 1))
    {
        exit(EXIT_FAILURE);
    }
    cpgask(1);

    display_scatter_char(Array, xsize, ysize, 1, "Plotting a 2-D array of char representing spins", "x-axis", "y-axis");

    cpgend();

    // find 4,7th element of array
    printf("The 4, 7th element of array: %c", interr_array(Array, 4, 7));

    destroy_2d_array_c(Array, xsize); // deallocate the memory
}