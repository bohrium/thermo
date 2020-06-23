#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*=============================================================================
======  0. GLOBAL VARIABLES  ==================================================
=============================================================================*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~  0.0. Hyperparameters  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------  0.0.0. thermodynamic parameters  --------------------*/ 

const float avg_energy = 0.1;
//const float density = 256.0;

//const int nb_molecules = (int)density;
#define nb_molecules 256

/*--------------------  0.0.1. interaction parameters  ----------------------*/ 

const float delta = 0.01;  
const float epsilon = 0.1;  

const float spring = 2.0 * epsilon / (delta * delta); 
float force_law(float dist)
{
    return (dist < delta) ? spring * dist : 0.0; 
}

/*--------------------  0.0.2. simulation parameters  -----------------------*/ 

const float dt = 0.001;
float t = 0.0;
const float T = 10.0;
const float print_t = 0.1;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~  0.1. Simulation State  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

float q[2][nb_molecules]; 
float p[2][nb_molecules]; 
float f[2][nb_molecules]; 
float f_new[2][nb_molecules]; 

/*=============================================================================
======  1. MAIN LOOP  =========================================================
=============================================================================*/

void init_state();
void reset_forces();
void update_positions(); 
void accumulate_forces();
void update_momenta(); 
void print_state(); 
void print_histogram();

void main()
{
    srand(time(0));

    printf("hi!\n");

    {
        init_state();
        reset_forces();
        reset_forces();
        accumulate_forces();
        //print_state();
        print_histogram();
    }
    for ( t = 0.0; t < T; t += dt ) {
        reset_forces();
        update_positions(); 
        accumulate_forces();
        update_momenta(); 

        if ( floor((t + dt)/print_t) == floor(t/print_t) ) { continue; }
        //print_state();
        print_histogram();
    }

    printf("bye!\n");
}

/*=============================================================================
======  2. IMPLEMENTATION  ====================================================
=============================================================================*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~  2.0. Randomness  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

float unif() { return ((float)rand())/RAND_MAX; } 
float coin() { return (rand()%2) * 2 - 1; }

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~  2.1. Physical Simulation  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void init_state()
{
    float p_scale = sqrt(2.0 * avg_energy); 
    for ( int i=0; i!=2; ++i ) {
        for ( int j=0; j!=nb_molecules; ++j ) {
            q[i][j] = unif() * 0.5 + 0.25;
            p[i][j] = (i==0) ? 0 : (p_scale*coin());
        }
    }
}

void reset_forces()
{
    for ( int i=0; i!=2; ++i ) {
        for ( int j=0; j!=nb_molecules; ++j ) {
            f[i][j] = f[i][j];
            f_new[i][j] = 0.0;
        }
    }
}

void accumulate_forces() 
{
    for ( int j=0; j!=nb_molecules; ++j ) {
        for ( int k=j+1; k!=nb_molecules; ++k ) {
            float diff[2];
            diff[0] = q[0][k] - q[0][j];  if ( 0.5 < diff[0] ) { diff[0] -= 1.0; } else if ( diff[0] < -0.5 ) { diff[0] += 1.0; }
            diff[1] = q[1][k] - q[1][j];  if ( 0.5 < diff[1] ) { diff[1] -= 1.0; } else if ( diff[1] < -0.5 ) { diff[1] += 1.0; }

            float dist = sqrt(diff[0]*diff[0] + diff[1]*diff[1]);
            float mag = force_law(dist) / dist;

            f_new[0][j] -= diff[0] * mag; 
            f_new[1][j] -= diff[1] * mag;

            f_new[0][k] += diff[0] * mag; 
            f_new[1][k] += diff[1] * mag;
        }
    }
}

void update_positions()
{
    for ( int i=0; i!=2; ++i ) {
        for ( int j=0; j!=nb_molecules; ++j ) {
            q[i][j] += (p[i][j] + 0.5*f[i][j]*dt) * dt; 
            q[i][j] -= floor(q[i][j]);
        }
    }
}

void update_momenta()
{
    for ( int i=0; i!=2; ++i ) {
        for ( int j=0; j!=nb_molecules; ++j ) {
            p[i][j] += 0.5 * (f[i][j] + f_new[i][j]) * dt; 
        }
    }
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~  2.2. Display  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#define GRIDSIZE 20
void print_state()
{
    int grid[GRIDSIZE][GRIDSIZE]; 
    for ( int r=0; r!=GRIDSIZE; ++r ) {
        for ( int c=0; c!=GRIDSIZE; ++c ) {
            grid[r][c] = 0;
        }
    }

    for ( int j = 0; j!=nb_molecules; ++j ) {
        int r = (int)(q[0][j] * ((float)GRIDSIZE-0.01));
        int c = (int)(q[1][j] * ((float)GRIDSIZE-0.01));
        grid[r][c] += 1;
    }

    printf(" ");
    for ( int c=0; c!=2*GRIDSIZE+1; ++c ) {
        printf("_");
    }
    printf(" \n");
    for ( int r=0; r!=GRIDSIZE; ++r ) {
        printf("| ");
        for ( int c=0; c!=GRIDSIZE; ++c ) {
            printf("%c ",
                grid[r][c] == 0 ? ' ' :
                grid[r][c] == 1 ? 'o' :     
                grid[r][c] == 2 ? '8' : '%' );
        }
        printf("|\n");
    }
}

#define HIST_MAX    1.0  
#define NB_BINS     20  
void print_histogram()
{
    int grid[NB_BINS];
    for ( int r=0; r!=NB_BINS; ++r ) {
        grid[r] = 0;
    }

    for ( int j = 0; j!=nb_molecules; ++j ) {
        float speed = sqrt(p[0][j]*p[0][j] + p[1][j]*p[1][j]);
        int bin = (speed/HIST_MAX) * NB_BINS;  if ( NB_BINS <= bin ) { bin = NB_BINS-1; }
        grid[bin] += 1;
    }

    for ( int c=0; c!=30; ++c ) { printf("="); }
    printf("  SPEED  at time %.4f  ", t);
    for ( int c=0; c!=30; ++c ) { printf("="); }
    printf("\n");
    for ( int r=0; r!=NB_BINS; ++r ) {
        printf("|%.2f|", (HIST_MAX * r)/NB_BINS);
        for ( int c=0; c!=grid[r]; ++c ) {
            printf("-");
        }
        printf("\n");
    }
}


