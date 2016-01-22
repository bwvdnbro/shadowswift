/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 *               2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/******************************************************************************
 * The Riemann solver in this file was written by Bert Vandenbroucke as part of
 * the moving mesh code Shadowfax and adapted for use with SWIFT. It consists
 * of an exact Riemann solver as described in
 *  Toro, Eleuterio F., Riemann Solvers and Numerical Methods for Fluid
 *  Dynamics, Springer (2009, 3rd edition)
 *
 ******************************************************************************/

#ifndef RIEMANN_EXACT_H
#define RIEMANN_EXACT_H
/* gives us const_hydro_gamma and tells us which floating point type to use */
#include "const.h"
#include "math.h"
#include "stdio.h"
#include "float.h"
#include "stdlib.h"
#include "error.h"

/* frequently used combinations of const_hydro_gamma */
#define const_riemann_gp1d2g    (0.5f*(const_hydro_gamma+1.0f)/const_hydro_gamma)
#define const_riemann_gm1d2g    (0.5f*(const_hydro_gamma-1.0f)/const_hydro_gamma)
#define const_riemann_gm1dgp1   ((const_hydro_gamma-1.0f)/(const_hydro_gamma+1.0f))
#define const_riemann_tdgp1     (2.0f/(const_hydro_gamma+1.0f))
#define const_riemann_tdgm1     (2.0f/(const_hydro_gamma-1.0f))
#define const_riemann_gm1d2     (0.5f*(const_hydro_gamma-1.0f))
#define const_riemann_tgdgm1    (2.0f*const_hydro_gamma/(const_hydro_gamma-1.0f))
#define const_riemann_ginv      (1.0f/const_hydro_gamma)

/**
 * @brief Functions (4.6) and (4.7) in Toro.
 *
 * @param p The current guess for the pressure
 * @param W The left or right state vector
 * @param a The left or right sound speed
 */
static float riemann_fb(float p, float* W, float a){
    float fval = 0.;
    float A, B;
    if(p > W[4]){
        A = const_riemann_tdgp1 / W[0];
        B = const_riemann_gm1dgp1 * W[4];
        fval = (p-W[4])*sqrtf(A/(p+B));
    } else {
        fval = const_riemann_tdgm1 * a * ( powf( p / W[4], const_riemann_gm1d2g ) - 1.0f );
    }
    return fval;
}

/**
 * @brief Function (4.5) in Toro
 *
 * @param p The current guess for the pressure
 * @param WL The left state vector
 * @param WR The right state vector
 * @param vL The left velocity along the interface normal
 * @param vR The right velocity along the interface normal
 * @param aL The left sound speed
 * @param aR The right sound speed
 */
static float riemann_f(float p, float* WL, float* WR, float vL, float vR, float aL, float aR){
    return riemann_fb(p, WL, aL) + riemann_fb(p, WR, aR) + ( vR - vL );
}

/**
 * @brief Function (4.37) in Toro
 *
 * @param p The current guess for the pressure
 * @param W The left or right state vector
 * @param a The left or right sound speed
 */
static float riemann_fprimeb(float p, float* W, float a){
    float fval = 0.;
    float A, B;
    if(p > W[4]){
        A = const_riemann_tdgp1 / W[0];
        B = const_riemann_gm1dgp1 * W[4];
        fval = (1.0f-0.5f*(p-W[4])/(B+p))*sqrtf(A/(p+B));
    } else {
        fval = 1.0f / ( W[0] * a ) * powf( p / W[4] , -const_riemann_gp1d2g );
    }
    return fval;
}

/**
 * @brief The derivative of riemann_f w.r.t. p
 *
 * @param p The current guess for the pressure
 * @param WL The left state vector
 * @param WR The right state vector
 * @param aL The left sound speed
 * @param aR The right sound speed
 */
static float riemann_fprime(float p, float* WL, float* WR, float aL, float aR){
    return riemann_fprimeb(p, WL, aL) + riemann_fprimeb(p, WR, aR);
}

/**
 * @brief Bottom function of (4.48) in Toro
 *
 * @param p The current guess for the pressure
 * @param W The left or right state vector
 */
static float riemann_gb(float p, float* W){
    float A, B;
    A = const_riemann_tdgp1 / W[0];
    B = const_riemann_gm1dgp1 * W[4];
    return sqrtf(A/(p+B));
}

/**
 * @brief Get a good first guess for the pressure in the iterative scheme
 *
 * This function is based on (4.47) and (4.48) in Toro and on the
 * FORTRAN code provided in Toro p.156-157
 *
 * @param WL The left state vector
 * @param WR The right state vector
 * @param vL The left velocity along the interface normal
 * @param vR The right velocity along the interface normal
 * @param aL The left sound speed
 * @param aR The right sound speed
 */
static float riemann_guess_p(float* WL, float* WR, float vL, float vR, float aL, float aR){
    float pguess, pmin, pmax, qmax;
    float ppv;

    pmin = fminf(WL[4], WR[4]);
    pmax = fmaxf(WL[4], WR[4]);
    qmax = pmax/pmin;
    ppv = 0.5f * ( WL[4] + WR[4] ) - 0.125f * ( vR - vL ) * ( WL[0] + WR[0] ) * ( aL + aR );
    ppv = fmaxf(1.e-8f, ppv);
    if(qmax <= 2.0f && pmin <= ppv && ppv <= pmax){
        pguess = ppv;
    } else {
        if(ppv < pmin){
            /* two rarefactions */
            pguess = powf( ( aL + aR - const_riemann_gm1d2 * ( vR - vL ) ) /
                           ( aL / powf( WL[4], const_riemann_gm1d2g ) + aR / powf( WR[4], const_riemann_gm1d2g ) )
                           , const_riemann_tgdgm1 );
        } else {
            /* two shocks */
            pguess = ( riemann_gb(ppv, WL) * WL[4] + riemann_gb(ppv, WR) * WR[4] - vR + vL ) /
                     ( riemann_gb(ppv, WL) + riemann_gb(ppv, WR) );
        }
    }
    /* Toro: "Not that approximate solutions may predict, incorrectly, a negative value for pressure (...).
       Thus in order to avoid negative guess values we introduce the small positive constant _tolerance" */
    pguess = fmaxf(1.e-8f, pguess);
    return pguess;
}

/**
 * @brief Find the zeropoint of riemann_f(p) using Brent's method
 *
 * @param lower_limit Lower limit for the method (riemann_f(lower_limit) < 0)
 * @param upper_limit Upper limit for the method (riemann_f(upper_limit) > 0)
 * @param error_tol Tolerance used to decide if the solution is converged
 * @param WL Left state vector
 * @param WR Right state vector
 * @param vL The left velocity along the interface normal
 * @param vR The right velocity along the interface normal
 * @param aL The left sound speed
 * @param aR The right sound speed
 */
static float riemann_solve_brent(float lower_limit, float upper_limit, float lowf, float upf, float error_tol, float* WL, float* WR, float vL, float vR, float aL, float aR) {
    float a, b, c, d, s;
    float fa, fb, fc, fs;
    float tmp, tmp2;
    int mflag;
    int i;

    a = lower_limit;
    b = upper_limit;
    c = 0.0f;
    d = FLT_MAX;

    fa = lowf;
    fb = upf;

    fc = 0.0f;
    s = 0.0f;
    fs = 0.0f;

    /* if f(a) f(b) >= 0 then error-exit */
    if ( fa * fb >= 0.0f ) {
        error("Brent's method called with equal sign function values!\n"
              "f(%g) = %g, f(%g) = %g\n", a, fa, b, fb);
        /* return NaN */
        return 0.0f / 0.0f;
    }

    /* if |f(a)| < |f(b)| then swap (a,b) */
    if ( fabs(fa) < fabs(fb) ){
        tmp = a;
        a = b;
        b = tmp;
        tmp = fa;
        fa = fb;
        fb = tmp;
    }

    c = a;
    fc = fa;
    mflag = 1;
    i = 0;

    while ( !( fb == 0.0f ) && ( fabs( a - b ) > error_tol * 0.5f * ( a + b) ) ){
        if ( ( fa != fc ) && ( fb != fc ) )
            /* Inverse quadratic interpolation */
            s = a * fb * fc / ( fa - fb ) / ( fa - fc )
                    + b * fa * fc / ( fb - fa ) / ( fb - fc )
                    + c * fa * fb / ( fc - fa ) / ( fc - fb );
        else
            /* Secant Rule */
            s = b - fb * ( b - a ) / ( fb - fa );

        tmp2 = 0.25f * (3.0f * a + b);
        if ( !( ( ( s > tmp2 ) && ( s < b ) ) || ( ( s < tmp2 ) && ( s > b ) ) )
             || ( mflag && ( fabs( s - b ) >= ( 0.5f * fabs( b - c ) ) ) )
             || ( !mflag && ( fabs( s - b ) >= ( 0.5f * fabs( c - d ) ) ) )
             || ( mflag && ( fabs( b - c ) < error_tol ) )
             || ( !mflag && ( fabs( c - d ) < error_tol ) ) ){
            s = 0.5f * (a + b);
            mflag = 1;
        } else {
            mflag = 0;
        }
        fs = riemann_f(s, WL, WR, vL, vR, aL, aR);
        d = c;
        c = b;
        fc = fb;
        if ( fa * fs < 0. ){
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }

        /* if |f(a)| < |f(b)| then swap (a,b) */
        if( fabs(fa) < fabs(fb) ){
            tmp = a;
            a = b;
            b = tmp;
            tmp = fa;
            fa = fb;
            fb = tmp;
        }
        i++;
    }
    return b;
}

/**
 * @brief Vacuum Riemann solver, based on section 4.6 in Toro
 *
 * @param WL The left state vector
 * @param WR The right state vector
 * @param vL The left velocity along the interface normal
 * @param vR The right velocity along the interface normal
 * @param aL The left sound speed
 * @param aR The right sound speed
 * @param Whalf Empty state vector to store the solution in
 * @param n_unit Normal vector of the interface
 */
static void riemann_solve_vacuum(float* WL, float* WR, float vL, float vR, float aL, float aR, float* Whalf, float* n_unit){
    float SL, SR;
    float vhalf;

//    printf("Vacuum solver...\n");

    if( !WR[0] && !WL[0] ){
        /* if both states are vacuum, the solution is also vacuum */
        Whalf[0] = 0.0f;
        Whalf[1] = 0.0f;
        Whalf[2] = 0.0f;
        Whalf[3] = 0.0f;
        Whalf[4] = 0.0f;
        return;
    }
    if( !WR[0] ){
        Whalf[1] = WL[1];
        Whalf[2] = WL[2];
        Whalf[3] = WL[3];
        /* vacuum right state */
        if( vL < aL ){
            SL = vL + const_riemann_tdgm1 * aL;
            if( SL > 0.0f ){
                Whalf[0] = WL[0] * powf( const_riemann_tdgp1 + const_riemann_gm1dgp1 / aL * vL, const_riemann_tdgm1 );
                vhalf = const_riemann_tdgp1 * ( aL + const_riemann_gm1d2 * vL ) - vL;
                Whalf[4] = WL[4] * powf( const_riemann_tdgp1 + const_riemann_gm1dgp1 / aL * vL, const_riemann_tgdgm1 );
            } else {
                Whalf[0] = 0.0f;
                Whalf[1] = 0.0f;
                Whalf[2] = 0.0f;
                Whalf[3] = 0.0f;
                Whalf[4] = 0.0f;
                return;
            }
        } else {
            Whalf[0] = WL[0];
            vhalf = 0.0f;
            Whalf[4] = WL[4];
        }
    } else {
        if( !WL[0] ){
            Whalf[1] = WR[1];
            Whalf[2] = WR[2];
            Whalf[3] = WR[3];
            /* vacuum left state */
            if( -vR < aR ){
                SR = vR - const_riemann_tdgm1 * aR;
                if( SR >= 0.0f ){
                    Whalf[0] = 0.0f;
                    Whalf[1] = 0.0f;
                    Whalf[2] = 0.0f;
                    Whalf[3] = 0.0f;
                    Whalf[4] = 0.0f;
                    return;
                } else {
                    Whalf[0] = WR[0] * powf( const_riemann_tdgp1 - const_riemann_gm1dgp1 / aR * vR, const_riemann_tdgm1 );
                    vhalf = const_riemann_tdgp1 * ( -aR + const_riemann_gm1d2 * vR ) - vR;
                    Whalf[4] = WR[4] * powf( const_riemann_tdgp1 - const_riemann_gm1dgp1 / aR * vR, const_riemann_tgdgm1 );
                }
            } else {
                Whalf[0] = WR[0];
                vhalf = 0.0f;
                Whalf[4] = WR[4];
            }
        } else {
            /* vacuum generation */
            SR = vR - const_riemann_tdgm1 * aR;
            SL = vL + const_riemann_tdgm1 * aL;
            if( SR > 0.0f && SL < 0.0f ){
                Whalf[0] = 0.0f;
                Whalf[1] = 0.0f;
                Whalf[2] = 0.0f;
                Whalf[3] = 0.0f;
                Whalf[4] = 0.0f;
                return;
            } else {
                if( SL >= 0.0f ){
                    Whalf[1] = WL[1];
                    Whalf[2] = WL[2];
                    Whalf[3] = WL[3];
                    if( aL > vL ){
                        Whalf[0] = WL[0] * powf( const_riemann_tdgp1 + const_riemann_gm1dgp1 / aL * vL, const_riemann_tdgm1 );
                        vhalf = const_riemann_tdgp1 * ( aL + const_riemann_gm1d2 * vL ) - vL;
                        Whalf[4] = WL[4] * powf( const_riemann_tdgp1 + const_riemann_gm1dgp1 / aL * vL, const_riemann_tgdgm1 );
                    } else {
                        Whalf[0] = WL[0];
                        vhalf = 0.0f;
                        Whalf[4] = WL[4];
                    }
                } else {
                    Whalf[1] = WR[1];
                    Whalf[2] = WR[2];
                    Whalf[3] = WR[3];
                    if( -vR < aR ){
                        Whalf[0] = WR[0] * powf( const_riemann_tdgp1 - const_riemann_gm1dgp1 / aR * vR, const_riemann_tdgm1 );
                        vhalf = const_riemann_tdgp1 * ( -aR + const_riemann_gm1d2 * vR ) - vR;
                        Whalf[4] = WR[4] * powf( const_riemann_tdgp1 - const_riemann_gm1dgp1 / aR * vR, const_riemann_tgdgm1 );
                    } else {
                        Whalf[0] = WR[0];
                        vhalf = 0.0f;
                        Whalf[4] = WR[4];
                    }
                }
            }
        }
    }

    /* Add the velocity solution along the interface normal to the velocities */
    Whalf[1] += vhalf * n_unit[0];
    Whalf[2] += vhalf * n_unit[1];
    Whalf[3] += vhalf * n_unit[2];
}

/* Solve the Riemann problem between the states WL and WR and store the result in Whalf
 * The Riemann problem is solved in the x-direction; the velocities in the y- and z-direction
 * are simply advected.
 */
/**
 * @brief Solve the Riemann problem between the given left and right state and along the given interface normal
 *
 * Based on chapter 4 in Toro
 *
 * @param WL The left state vector
 * @param WR The right state vector
 * @param Whalf Empty state vector in which the result will be stored
 * @param n_unit Normal vector of the interface
 */
static void riemann_solver_solve(float* WL, float* WR, float* Whalf, float *n_unit){
    /* velocity of the left and right state in a frame aligned with n_unit */
    float vL, vR, vhalf;
    /* sound speeds */
    float aL, aR;
    /* variables used for finding pstar */
    float p, pguess, fp, fpguess;
    /* variables used for sampling the solution */
    float u;
    float pdpR, SR;
    float SHR, STR;
    float pdpL, SL;
    float SHL, STL;
    int errorFlag = 0;

    /* sanity checks */
    if(WL[0] != WL[0] || WL[4] != WL[4]){
        printf("NaN WL!\n");
        errorFlag = 1;
    }
    if(WR[0] != WR[0] || WR[4] != WR[4]){
        printf("NaN WR!\n");
        errorFlag = 1;
    }
    if(WL[0] < 0.0f || WL[4] < 0.0f){
        printf("Negative WL!\n");
        errorFlag = 1;
    }
    if(WR[0] < 0.0f || WR[4] < 0.0f){
        printf("Negative WR!\n");
        errorFlag = 1;
    }
    if(errorFlag){
        printf("WL: %g %g %g %g %g\n", WL[0], WL[1], WL[2], WL[3], WL[4]);
        printf("WR: %g %g %g %g %g\n", WR[0], WR[1], WR[2], WR[3], WR[4]);
        error("Riemman solver input error!\n");
    }

    /* calculate velocities in interface frame */
    vL = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
    vR = WR[1] * n_unit[0] + WR[2] * n_unit[1] + WR[3] * n_unit[2];

    /* calculate sound speeds */
    aL = sqrtf( const_hydro_gamma * WL[4] / WL[0] );
    aR = sqrtf( const_hydro_gamma * WR[4] / WR[0] );

    if(!WL[0] || !WR[0]){
        /* vacuum: we need a vacuum riemann solver */
        riemann_solve_vacuum(WL, WR, vL, vR, aL, aR, Whalf, n_unit);
        return;
    }

    /* check vacuum generation condition */
    if( 2.0f * aL / ( const_hydro_gamma - 1.0f ) + 2.0f * aR / ( const_hydro_gamma - 1.0f ) < fabs( vL - vR ) ){
        /* vacuum generation: need a vacuum riemann solver */
        riemann_solve_vacuum(WL, WR, vL, vR, aL, aR, Whalf, n_unit);
        return;
    } else {
        /* values are ok: let's find pstar (riemann_f(pstar) = 0)! */
        /* We normally use a Newton-Raphson iteration to find the zeropoint
           of riemann_f(p), but if pstar is close to 0, we risk negative p values.
           Since riemann_f(p) is undefined for negative pressures, we don't
           want this to happen.
           We therefore use Brent's method if riemann_f(0) is larger than some
           value. -5 makes the iteration fail safe while almost never invoking
           the expensive Brent solver. */
        p = 0.;
        /* obtain a first guess for p */
        pguess = riemann_guess_p(WL, WR, vL, vR, aL, aR);
        fp = riemann_f(p, WL, WR, vL, vR, aL, aR);
        fpguess = riemann_f(pguess, WL, WR, vL, vR, aL, aR);
        /* if fp is negative enough, we don't want to use Brent's method ever */
//        if( fp < -5.0f ){
//            unsigned int counter = 0;
//            while( fabs(p-pguess) > 1.e-6f * 0.5f * ( p + pguess ) ){
//                p = pguess;
//                pguess = pguess - fpguess / riemann_fprime(pguess, WL, WR, aL, aR);
//                if( pguess < 0.0f ){
//                    pguess = 0.0f;
//                }
//                fpguess = riemann_f(pguess, WL, WR, vL, vR, aL, aR);
//                counter++;
//                if(counter > 1000){
//                    message("WL: %g %g %g\n", WL[0], vL, WL[4]);
//                    message("WR: %g %g %g\n", WR[0], vR, WR[4]);
//                    error("Stuck in Newton-Raphson!\n");
//                }
//            }
//            p = pguess;
//        } else {
            /* ok, pstar is close to 0, better use Brent's method... */
            /* we use Newton-Raphson until we find a suitable interval */
            if( fp * fpguess >= 0.0f ){
                /* Newton-Raphson until convergence or until suitable interval is found to use Brent's method */
                unsigned int counter = 0;
                while( fabs( p - pguess ) > 1.e-6f * 0.5f * ( p + pguess ) && fpguess < 0.0f ){
                    p = pguess;
                    pguess = pguess - fpguess / riemann_fprime(pguess, WL, WR, aL, aR);
                    fpguess = riemann_f(pguess, WL, WR, vL, vR, aL, aR);
                    counter++;
                    if(counter > 1000){
                        error("Stuck in Newton-Raphson!\n");
                    }
                }
            }
            /* As soon as there is a suitable interval: use Brent's method */
            if( 1.e6 * fabs( p - pguess ) > 0.5f * ( p + pguess ) && fpguess > 0.0f ){
                p = 0.0f;
                fp = riemann_f(p, WL, WR, vL, vR, aL, aR);
                /* use Brent's method to find the zeropoint */
                p = riemann_solve_brent(p, pguess, fp, fpguess, 1.e-6, WL, WR, vL, vR, aL, aR);
            } else {
                p = pguess;
            }
//        }

        /* calculate the velocity in the intermediate state */
        u = 0.5f * ( vL + vR ) + 0.5f * ( riemann_fb(p, WR, aR) - riemann_fb(p, WL, aL) );

        /* sample the solution */
        /* This corresponds to the flow chart in Fig. 4.14 in Toro */
        if( u < 0.0f ){
            /* advect velocity components */
            Whalf[1] = WR[1];
            Whalf[2] = WR[2];
            Whalf[3] = WR[3];
            pdpR = p / WR[4];
            if( p > WR[4] ){
                /* shockwave */
                SR = vR + aR * sqrtf( const_riemann_gp1d2g * pdpR + const_riemann_gm1d2g );
                if( SR > 0.0f ){
                    Whalf[0] = WR[0] * ( pdpR + const_riemann_gm1dgp1 ) / (const_riemann_gm1dgp1 * pdpR + 1.0f );
                    vhalf = u - vR;
                    Whalf[4] = p;
                } else {
                    Whalf[0] = WR[0];
                    vhalf = 0.0f;
                    Whalf[4] = WR[4];
                }
            } else {
                /* rarefaction wave */
                SHR = vR + aR;
                if( SHR > 0.0f ){
                    STR = u + aR * powf( pdpR, const_riemann_gm1d2g );
                    if( STR <= 0.0f ){
                        Whalf[0] = WR[0] * powf( const_riemann_tdgp1 - const_riemann_gm1dgp1 / aR * vR, const_riemann_tdgm1 );
                        vhalf = const_riemann_tdgp1 * ( -aR + const_riemann_gm1d2 * vR ) - vR;
                        Whalf[4] = WR[4] * powf( const_riemann_tdgp1 - const_riemann_gm1dgp1 / aR * vR, const_riemann_tgdgm1 );
                    } else {
                        Whalf[0] = WR[0] * powf( pdpR, const_riemann_ginv );
                        vhalf = u - vR;
                        Whalf[4] = p;
                    }
                } else {
                    Whalf[0] = WR[0];
                    vhalf = 0.0f;
                    Whalf[4] = WR[4];
                }
            }
        } else {
            Whalf[1] = WL[1];
            Whalf[2] = WL[2];
            Whalf[3] = WL[3];
            pdpL = p / WL[4];
            if( p > WL[4] ){
                /* shockwave */
                SL = vL - aL * sqrtf( const_riemann_gp1d2g * pdpL + const_riemann_gm1d2g );
                if( SL < 0.0f ){
                    Whalf[0] = WL[0] * ( pdpL + const_riemann_gm1dgp1 ) / ( const_riemann_gm1dgp1 * pdpL + 1.0f );
                    vhalf = u - vL;
                    Whalf[4] = p;
                } else {
                    Whalf[0] = WL[0];
                    vhalf = 0.0f;
                    Whalf[4] = WL[4];
                }
            } else {
                /* rarefaction wave */
                SHL = vL - aL;
                if( SHL < 0.0f ){
                    STL = u - aL * powf( pdpL, const_riemann_gm1d2g );
                    if( STL > 0.0f ){
                        Whalf[0] = WL[0] * powf( const_riemann_tdgp1 + const_riemann_gm1dgp1 / aL * vL, const_riemann_tdgm1 );
                        vhalf = const_riemann_tdgp1 * ( aL + const_riemann_gm1d2 * vL ) - vL;
                        Whalf[4] = WL[4] * powf( const_riemann_tdgp1 + const_riemann_gm1dgp1 / aL * vL, const_riemann_tgdgm1);
                    } else {
                        Whalf[0] = WL[0] * powf( pdpL, const_riemann_ginv );
                        vhalf = u - vL;
                        Whalf[4] = p;
                    }
                } else {
                    Whalf[0] = WL[0];
                    vhalf = 0.0f;
                    Whalf[4] = WL[4];
                }
            }
        }
    }

    /* add the velocity solution along the interface normal to the velocities */
    Whalf[1] += vhalf*n_unit[0];
    Whalf[2] += vhalf*n_unit[1];
    Whalf[3] += vhalf*n_unit[2];
}

#endif // RIEMANN_EXACT_H
