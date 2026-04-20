#ifndef SPODY_ECLIPSE_H
#define SPODY_ECLIPSE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

//debug
#define DEBUG_ECLIPSE 0 // 0 = no debug | 1 = debug |---> CODE TESTING

#define BUFFER_SIZE_ECL 256
#define PI 3.14159265358979323846264338327950288419716939937511

double spody_get_suneclipse( double pos[3], double unit_sat2sun_pos[3], double body_app_r, double sun_app_r);
double spody_get_sateclipsestatus( double occulting2sat_pos[3], double occulting2sun_pos[3], double sat2sun_pos[3], double sun_rad, double body_rad);

#ifdef __cplusplus
}
#endif

#endif // SPODY_ECLIPSE_H