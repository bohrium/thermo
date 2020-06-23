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
#define nb_molecules 256

/*--------------------  0.0.1. interaction parameters  ----------------------*/ 

const float delta = 0.005;  
const float epsilon = 0.1;  

const float spring = 2.0 * epsilon / (delta * delta); 
float force_law(float dist)
{
    return (dist < delta) ? spring * dist : 0.0; 
}

/*--------------------  0.0.2. simulation parameters  -----------------------*/ 

const float dt = delta / 100;
const float T = 50.0;
const float print_t = 0.03;   
float t = 0.0;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~  0.1. Simulation State  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#define KDSIZE 16 
int lengths[KDSIZE][KDSIZE];
int items[KDSIZE][KDSIZE][nb_molecules];

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
void print_histogram(int offset);

void main()
{
    srand(time(0));

    printf("hi!\n");

    {
        init_state();
        reset_forces();
        reset_forces();
        accumulate_forces();
    }
    clock_t start = clock();
    for ( t = 0.0; t < T; t += dt ) {
        if ( floor((t + dt)/print_t) != floor(t/print_t) ) {
            print_state();
            print_histogram(210);
        }
        reset_forces();
        update_positions(); 
        accumulate_forces();
        update_momenta(); 

        while ( clock() < start + (int)(CLOCKS_PER_SEC * t) ) {}
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
~~~~~~~~~~~~~~  2.1. Leapfrog Integration Helpers  ~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void init_state()
{
    float p_scale = sqrt(2.0 * avg_energy); 
    for ( int i=0; i!=2; ++i ) {
        for ( int j=0; j!=nb_molecules; ++j ) {
            q[i][j] = unif()*0.8 + 0.1;
            p[i][j] = (i==0) ? 0.0 : p_scale * coin();
            //p[i][j] = p_scale * coin();
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
~~~~~~~~~~~~~~  2.2. Smart Force Computations  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void assign_to_kd() 
{
    for ( int r=0; r!=KDSIZE; ++r ) {
        for ( int c=0; c!=KDSIZE; ++c ) {
            lengths[r][c] = 0;
        }
    }

    for ( int j = 0; j!=nb_molecules; ++j ) {
        int r = (int)(q[0][j] * ((float)KDSIZE-0.0001));
        int c = (int)(q[1][j] * ((float)KDSIZE-0.0001));
        items[r][c][lengths[r][c]] = j;
        lengths[r][c] += 1;
    }
}

void accumulate_forces_inner(int* selves, int nb_selves, int* others, int nb_others); 

void accumulate_forces()
{
    assign_to_kd();
    for ( int r=0; r!=KDSIZE; ++r ) {
        for ( int c=0; c!=KDSIZE; ++c ) {
            int nb_others = 0;
            int others[nb_molecules];

            for ( int dr=0; dr!=2; ++dr) {
                for ( int dc=0; dc!=2; ++dc) {
                    if ( dr==0 && dc==0 ) { continue; }
                    int rr = (r+dr+KDSIZE) % KDSIZE;
                    int cc = (c+dc+KDSIZE) % KDSIZE;
                    for ( int j=0; j!=lengths[rr][cc]; ++j ) {
                        others[nb_others] = items[rr][cc][j];  
                        nb_others += 1;
                    }
                }
            }
            accumulate_forces_inner(items[r][c], lengths[r][c], others, nb_others);
        }
    }
}

void accumulate_forces_inner(int* selves, int nb_selves, int* others, int nb_others) 
{
    for ( int j_idx=0; j_idx!=nb_selves; ++j_idx ) {
        for ( int k_idx=j_idx+1; k_idx!=nb_selves + nb_others; ++k_idx ) {
            int j = selves[j_idx];
            int k = (k_idx < nb_selves) ? selves[k_idx] : others[k_idx-nb_selves];

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

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~  2.3. Display  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#define GRIDSIZE 100
void print_state()
{
    int grid[GRIDSIZE][2*GRIDSIZE]; 
    for ( int r=0; r!=GRIDSIZE; ++r ) {
        for ( int c=0; c!=2*GRIDSIZE; ++c ) {
            grid[r][c] = 0;
        }
    }

    for ( int j = 0; j!=nb_molecules; ++j ) {
        int r = (int)(q[0][j] * ((float)GRIDSIZE-0.01));
        int c = (int)(q[1][j] * (2*(float)GRIDSIZE-0.01));
        grid[r][c] += 1;
    }

    printf(" ");
    for ( int c=0; c!=2*GRIDSIZE+1; ++c ) {
        printf("_");
    }
    printf(" \n");
    for ( int r=0; r!=GRIDSIZE; ++r ) {
        printf("| ");
        for ( int c=0; c!=2*GRIDSIZE; ++c ) {
            printf("%c",
                grid[r][c] == 0 ? ' ' :
                grid[r][c] == 1 ? 'o' :     
                grid[r][c] == 2 ? '8' : '%' );
        }
        printf("|\n");
    }
    for ( int r=0; r!=GRIDSIZE+1; ++r ) {
        printf("\033[1A");
    }
}

#define HIST_MAX    1.0  
#define NB_BINS     50  
void print_histogram(int offset)
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

    for ( int c=0; c!=offset; ++c ) { printf("\033[1C"); }
    for ( int c=0; c!=30; ++c ) { printf("="); }
    printf("  SPEED  at time %.4f  ", t);
    for ( int c=0; c!=30; ++c ) { printf("="); }
    printf("\n"); 

    for ( int r=0; r!=NB_BINS; ++r ) {
        for ( int c=0; c!=offset; ++c ) { printf("\033[1C"); }
        printf("|%.2f|", (HIST_MAX * r)/NB_BINS);
        int c = 0;
        for ( ; c!=grid[r]; ++c ) {
            printf("-");
            if ( c == 75 ) { printf("[-]"); break; }
        }
        for ( ; c!= 80; ++c ) { printf(" "); }
        printf("\n");
    }
    for ( int r=0; r!=NB_BINS+1; ++r ) {
        printf("\033[1A");
    }
}
