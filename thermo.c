/*  author: samtenka
 *  change: 2020-06-23
 *  create: 2020-06-22
 *  descrp: Simulation of a gas with very short-range interactions.
 *  to use: Run `make all`, then run `./out.o`.  Two displays will appear in
 *          the terminal: to the left, a large square with moving dots
 *          representing the Particle locations; to the right, a histogram of
 *          speeds, showing both Measurement and equilibrium-theory.  Observe
 *          how the system equilibrates With time.  One may change the
 *          hardcoded hyperparameters in the Hyperparameters section below.
 */

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
#define nb_molecules 512  

/*--------------------  0.0.1. interaction parameters  ----------------------*/ 

const float delta = 0.001;   
const float epsilon = 0.5;  

const float spring = 2.0 * epsilon / (delta * delta); 
float force_law(float dist)
{
    return (dist < delta) ? spring * (delta-dist) : 0.0; 
}

/*--------------------  0.0.2. simulation parameters  -----------------------*/ 

const float dt = delta / 25; 
const float T = 5000.0;
const float print_t = 0.04;   
float t = 0.0;
float animation_rate = 1.0 ;

/*--------------------  0.0.3. display parameters  --------------------------*/ 

#define GRIDSIZE (2*50)

#define HIST_MAX    1.25 
#define NB_BINS     25
int hist_meas[NB_BINS];
float hist_pred[NB_BINS];

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~  0.1. Simulation State  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#define KDSIZE 32
int indices[nb_molecules];
int lengths[KDSIZE][KDSIZE];
int offsets[KDSIZE][KDSIZE];

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
void measure_histogram();
void predict_histogram();
void print_histogram(int offset, float scale);

void main()
{
    srand(time(0));

    printf("hi!\n");
    if ( 1.0/KDSIZE < delta ) {
        printf(" UH OH! \n");
        exit(0);
    }

    predict_histogram();

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
            measure_histogram();
            print_histogram(GRIDSIZE+5, nb_molecules/256.0);
        }
        reset_forces();
        update_positions(); 
        accumulate_forces();
        update_momenta(); 

        while ( clock() < start + (int)(CLOCKS_PER_SEC * t / animation_rate) ) {}
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
    for ( int j=0; j!=nb_molecules; ++j ) {
        q[0][j] = unif()*0.4 + 0.3;
        q[1][j] = unif();
        float angle = 2*M_PI*unif();
        p[0][j] = p_scale * sin(angle); 
        p[1][j] = p_scale * cos(angle); 
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
    int visited[KDSIZE][KDSIZE];

    for ( int r=0; r!=KDSIZE; ++r ) {
        for ( int c=0; c!=KDSIZE; ++c ) {
            lengths[r][c] = 0;
            visited[r][c] = 0;
        }
    }

    for ( int j = 0; j!=nb_molecules; ++j ) {
        int r = (int)(q[0][j] * ((float)KDSIZE-0.0001));
        int c = (int)(q[1][j] * ((float)KDSIZE-0.0001));
        lengths[r][c] += 1;
    }

    int accum = 0;
    for ( int r=0; r!=KDSIZE; ++r ) {
        for ( int c=0; c!=KDSIZE; ++c ) {
            offsets[r][c] = accum;
            accum += lengths[r][c]; 
        }
    }

    for ( int j = 0; j!=nb_molecules; ++j ) {
        int r = (int)(q[0][j] * ((float)KDSIZE-0.0001));
        int c = (int)(q[1][j] * ((float)KDSIZE-0.0001));
        indices[offsets[r][c] + visited[r][c]] = j;
        visited[r][c] += 1;
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
                        others[nb_others] = indices[offsets[rr][cc] + j];
                        nb_others += 1;
                    }
                }
            }
            accumulate_forces_inner(indices+offsets[r][c], lengths[r][c], others, nb_others);
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

void print_state()
{
    int grid_special[GRIDSIZE][GRIDSIZE]; 
    int grid_normal[GRIDSIZE][GRIDSIZE]; 
    for ( int r=0; r!=GRIDSIZE; ++r ) {
        for ( int c=0; c!=GRIDSIZE; ++c ) {
            grid_special[r][c] = 0;
            grid_normal[r][c] = 0;
        }
    }

    for ( int j = 0; j!=nb_molecules; ++j ) {
        int r = (int)(q[0][j] * ((float)GRIDSIZE-0.0001));
        int c = (int)(q[1][j] * ((float)GRIDSIZE-0.0001));
        if ( j < 8 ) {
            grid_special[r][c] += 1;
        } else {
            grid_normal[r][c] += 1;
        }
    }

    printf(" ");
    for ( int c=0; c!=GRIDSIZE+1; ++c ) {
        printf("_");
    }
    printf(" \n");
    for ( int r=0; r!=GRIDSIZE; r+=2 ) {
        printf("| ");
        printf("\033[34m");
        for ( int c=0; c!=GRIDSIZE; ++c ) {
            printf("%s",
                grid_special[r][c] ? ( grid_special[r+1][c] ? "\033[31m8\033[34m" : "\033[31m\u1d3c\033[34m")
                                   : ( grid_special[r+1][c] ? "\033[31mo\033[34m" :
                grid_normal[r][c] ? ( grid_normal[r+1][c] ? ":" : "\u02d9")
                                  : ( grid_normal[r+1][c] ? "." : " " )
                )
            );
        }
        printf("\033[36m");
        printf("|\n");
    }
    for ( int r=0; r!=GRIDSIZE/2+1; ++r ) {
        printf("\033[1A");
    }
}

float maxwell(float speed)
{
    float norm = avg_energy;
    return speed * exp(- speed*speed/(2.0*avg_energy)) / norm;
}

void predict_histogram()
{
    for ( int r=0; r!=NB_BINS; ++r ) {
        float speed_a = (HIST_MAX * (r+0.0))/NB_BINS;
        float speed_b = (HIST_MAX * (r+0.5))/NB_BINS;
        float speed_c = (HIST_MAX * (r+1.0))/NB_BINS;
        float prob_density = (maxwell(speed_a)+2*maxwell(speed_b)+maxwell(speed_c))/4.0;
        float prob_mass = prob_density * nb_molecules * ((float)HIST_MAX)/NB_BINS;
        hist_pred[r] = prob_mass;
    }
}

void measure_histogram()
{
    for ( int r=0; r!=NB_BINS; ++r ) { hist_meas[r] = 0; }

    for ( int j = 0; j!=nb_molecules; ++j ) {
        float speed = sqrt(p[0][j]*p[0][j] + p[1][j]*p[1][j]);
        int bin = (speed/HIST_MAX) * NB_BINS;
        if ( NB_BINS <= bin ) { bin = NB_BINS-1; }
        hist_meas[bin] += 1;
    }
}

void print_many(char const* str, int times) {
    for ( int i=0; i!=times; ++i ) { printf("%s", str); }
} 

void print_histogram(int offset, float scale)
{
    print_many("\033[1C", offset);
    print_many("+", 30); printf("  SPEED  at time %.4f  ", t); print_many("+", 30); printf("\n"); 

    float chi_squared = 0.0; 

    for ( int r=0; r!=NB_BINS; ++r ) {
        float speed_a = (HIST_MAX * (r+0.0))/NB_BINS;
        print_many("\033[1C", offset);
        printf("|%.2f|", speed_a);

        {
            int c = 0;
            for ( ; c!=75; ++c ) {
                if      ( c==(int)(hist_pred[r]/scale) ) { printf("|"); }
                else if ( c==(int)(hist_meas[r]/scale) ) { printf("<"); }
                else if ( c<=(int)(hist_meas[r]/scale) ) { printf("-"); }
                else                                     { printf(" "); }
            }
        }

        printf("\n");

        float diff = ((float)(hist_meas[r]) - hist_pred[r]);
        chi_squared += diff*diff/hist_pred[r];
    }

    print_many("\033[1C", offset);
    printf("Chi squared is {%8.2f} (vs {%5.2f} at equilibrium)\n",
           chi_squared, (float)NB_BINS
    );
    if ( chi_squared < NB_BINS ) { exit(0); }

    print_many("\033[1A", NB_BINS+2);
}
