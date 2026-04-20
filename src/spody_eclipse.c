#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "spody_eclipse.h"


double spody_get_suneclipse( double pos[3], double unit_sat2sun_pos[3], double body_app_r, double sun_app_r){
    // Montenbruck and Gill 
    // GMAT like 
    double r = sqrt( pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2] ); 
    double dot_c = -(pos[0]*unit_sat2sun_pos[0] + pos[1]*unit_sat2sun_pos[1] + pos[2]*unit_sat2sun_pos[2]) / r;
    if (dot_c > 1.0) dot_c = 1.0;
    if (dot_c < -1.0) dot_c = -1.0;
    double c = acos(dot_c);

    double aa = sun_app_r *sun_app_r;
    double bb = body_app_r *body_app_r;
    double x = (c * c + aa - bb) / (2 * c);
    double y = sqrt(fmax(0.0, aa - x * x)); 

    double arg1 = x / sun_app_r;
    double acos_safe1 = acos((arg1 < -1.0) ? -1.0 : (arg1 > 1.0 ? 1.0 : arg1));
    double arg2 = (c - x) / body_app_r;
    double acos_safe2 = acos((arg2 < -1.0) ? -1.0 : (arg2 > 1.0 ? 1.0 : arg2));

    double unmbra_area = (aa * acos_safe1) + (bb * acos_safe2) - (c * y);

    #if DEBUG_ECLIPSE == 1
    printf("[eclipse_08] pos: (%.12f, %.12f, %.12f) | unit_sat2sun_pos: (%.12f, %.12f, %.12f) | body_app_r: %.12f | sun_app_r: %.12f \n", pos[0], pos[1], pos[2], unit_sat2sun_pos[0], unit_sat2sun_pos[1], unit_sat2sun_pos[2], body_app_r, sun_app_r);
    printf("[eclipse_08] r: %.12f\n", r);
    printf("[eclipse_08] c: %.12f\n", c);
    printf("[eclipse_08] aa: %.12f\n", aa);
    printf("[eclipse_08] bb: %.12f\n", bb);
    printf("[eclipse_08] x: %.12f\n", x);
    printf("[eclipse_08] y: %.12f\n", y);
    printf("acos(x / sun_app_r): %.12f\n", acos(x / sun_app_r));
    printf("acos((c - x) / body_app_r): %.12f\n", acos((c - x) / body_app_r));
    printf("c * y: %.12f\n", c * y);
    printf("[eclipse_08] r: %.12f | c: %.12f | aa: %.12f | bb: %.12f | x: %.12f | y: %.12f | unmbra_area: %.12f\n", r, c, aa, bb, x, y, unmbra_area);
    printf("sun area: %.12f\n", (PI * aa));
    #endif

    return 1.0 - unmbra_area / (PI * aa); //percent of sunlight
}

double spody_get_sateclipsestatus( double occulting2sat_pos[3], double occulting2sun_pos[3], double sat2sun_pos[3], double sun_rad, double body_rad){
    // Montenbruck and Gill     
    // GMAT like    

    double occulting2sun_r = sqrt( occulting2sun_pos[0]*occulting2sun_pos[0] + occulting2sun_pos[1]*occulting2sun_pos[1] + occulting2sun_pos[2]*occulting2sun_pos[2] );
    double unit_occulting2sun_pos[3];
    unit_occulting2sun_pos[0] = occulting2sun_pos[0] / occulting2sun_r;
    unit_occulting2sun_pos[1] = occulting2sun_pos[1] / occulting2sun_r;
    unit_occulting2sun_pos[2] = occulting2sun_pos[2] / occulting2sun_r;
    double pos_dot_occulting2sun = occulting2sat_pos[0]*unit_occulting2sun_pos[0] + occulting2sat_pos[1]*unit_occulting2sun_pos[1] + occulting2sat_pos[2]*unit_occulting2sun_pos[2];
    
    if (pos_dot_occulting2sun > 0){
        
        #if DEBUG_ECLIPSE == 1
        printf("[eclipse_00][FAST] sunlight case, no other evaluation because -> pos dot occulting2sun : %.6f\n",pos_dot_occulting2sun);
        #endif
        
        return 1.0; //sunlight
    }else{
        double sat2sun_r = sqrt( sat2sun_pos[0]*sat2sun_pos[0] + sat2sun_pos[1]*sat2sun_pos[1] + sat2sun_pos[2]*sat2sun_pos[2] );
        double r = sqrt( occulting2sat_pos[0]*occulting2sat_pos[0] + occulting2sat_pos[1]*occulting2sat_pos[1] + occulting2sat_pos[2]*occulting2sat_pos[2] );
        
        #if DEBUG_ECLIPSE == 1
        printf("[eclipse_01] occulting2sun_r: %.12f | sat2sun_r: %.12f | r: %.12f \n", occulting2sun_r, sat2sun_r, r);
        printf("[eclipse_01] sun_rad: %.12f | body_rad: %.12f \n", sun_rad, body_rad);
        #endif
 
    
        double a = asin( sun_rad / sat2sun_r );

        
        if(body_rad >= r){
            return 0; // the sat is inside the body
        }
        
        
        double b = asin( body_rad / r );
        double unit_sat2sun_pos[3];
        unit_sat2sun_pos[0] = sat2sun_pos[0] / sat2sun_r;
        unit_sat2sun_pos[1] = sat2sun_pos[1] / sat2sun_r;
        unit_sat2sun_pos[2] = sat2sun_pos[2] / sat2sun_r;
        double unit_pos[3];
        unit_pos[0] = -occulting2sat_pos[0] / r; //sat2occulting x
        unit_pos[1] = -occulting2sat_pos[1] / r; //sat2occulting y
        unit_pos[2] = -occulting2sat_pos[2] / r; //sat2occulting z
        
        // a good way to avoid numerical errors
        double dot = unit_pos[0]*unit_sat2sun_pos[0] + unit_pos[1]*unit_sat2sun_pos[1] + unit_pos[2]*unit_sat2sun_pos[2];
        if (dot > 1.0) dot = 1.0;
        if (dot < -1.0) dot = -1.0;
        double c = acos(dot);
        
        double aplusb = a+b;
        double aminusb = fabs(a - b);

        #if DEBUG_ECLIPSE == 1
        printf("[eclipse_02] a (SUN) : %.12f | b (occulting body) : %.12f \n", a, b);
        printf("[eclipse_03] unit_pos: (%.12f, %.12f, %.12f) \n", unit_pos[0], unit_pos[1], unit_pos[2]);
        printf("[eclipse_04] unit_sat2sun_pos: (%.12f, %.12f, %.12f) \n", unit_sat2sun_pos[0], unit_sat2sun_pos[1], unit_sat2sun_pos[2]);
        printf("[eclipse_05] aplusb: %.12f | aminusb (abs): %.12f | c: %.12f \n", aplusb, aminusb, c);
        #endif

        if ( aplusb <= c ){

            #if DEBUG_ECLIPSE == 1
            printf("[eclipse_06] sunlight case \n");
            #endif

            return 1.0; //sunlight
        }else if ( c <= aminusb && b >= a ){

            #if DEBUG_ECLIPSE == 1
            printf("[eclipse_07] full eclipse case \n");
            #endif

            return 0; //full eclipse
        }else if ( aminusb < c && c < aplusb ){

            #if DEBUG_ECLIPSE == 1
            printf("[eclipse_08] partial eclipse case \n");
            printf("[eclipse_08] occulting2sat_pos: (%.12f, %.12f, %.12f) \n", occulting2sat_pos[0], occulting2sat_pos[1], occulting2sat_pos[2]);

            #endif

            double se = spody_get_suneclipse(occulting2sat_pos, unit_sat2sun_pos, b, a);

            #if DEBUG_ECLIPSE == 1
            printf("[eclipse_08] partial eclipse case -> suneclipse: %.12f \n", se);
            #endif 

            return se; //penumbra
        }else{

            #if DEBUG_ECLIPSE == 1
            printf("[eclipse_09] anteumbra case \n");
            #endif

            return 1 - (b*b)/(a*a); //anteumbra
        }
    }

}