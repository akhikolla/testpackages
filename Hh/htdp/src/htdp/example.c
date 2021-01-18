/**
 * Sample displacement and velocity calculation.
 *
 * gcc example.c -o example.exe -L. -lhtdp
 * 
 * Output:
 *
 * DN:   0.074013 meters
 * DE:  -0.001450 meters
 * DU:  -0.004262 meters
 * VN:  37.151775 mm/yr
 * VE: -25.826116 mm/yr
 * VU:  -1.329391 mm/yr
 */

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>

#include "htdp.h"

int main()
{   
    // Obtain parameters defining crustal motion model
    model_();
    
    // Initialize transformation parameters from ITRF94 to other reference frames
    settp_();
    
    // Initialize conversion table between reference frame identifiers
    setrf_();

    // Values for sample point 'beta'
    double xlat = 36.6698;
    double xlon = 121.7722;
    int m0 = 10;
    int d0 = 16;
    int y0 = 1989;
    int m1 = 10;
    int d1 = 18;
    int y1 = 1989;
    int iopt = 1;
    
    // local variables
    int min0 = 0, min1 = 0;
    double date0 = 0, date1 = 0;
    double lat0 = 0, lon0 = 0, eht0 = 0;
    double lat1 = 0, lon1 = 0, eht1 = 0;
    int jregn = 0;
    double vn = 0, ve = 0, vu = 0;
    double dn = 0, de = 0, du = 0;
    
    // lat/lon in positive N/W, radians
    lat0 = xlat * (M_PI / 180.0);
    lon0 = xlon * (M_PI / 180.0);
    
    // Calculate modified Julian date
    int rc = 0;
    rc = c_getmdy(m0, d0, y0, &date0, &min0);
    if (rc == 1) return 1;
    rc = c_getmdy(m1, d1, y1, &date1, &min1);
    if (rc == 1) return 1;
    
    // Predict velocity in iopt reference frame
    predv_(&lat0, &lon0, &eht0, &date0, &iopt, &jregn, &vn, &ve, &vu);
    
    // Predict coordinates and displacements from time MIN1 to time MIN2.
    newcor_(&lat0, &lon0, &eht0, &min0, &min1, &lat1, &lon1, &eht1, &dn, &de, &du, &vn, &ve, &vu);
    
    printf("DN: %10f meters\n", dn);
    printf("DE: %10f meters\n", de);
    printf("DU: %10f meters\n", du);
    printf("VN: %10f mm/yr\n", vn);
    printf("VE: %10f mm/yr\n", ve);
    printf("VU: %10f mm/yr\n", vu);
    
    return 0;
}

// Reimplementation of interactive subroutine GETMDY in htdp.for
int c_getmdy(const int month, const int iday, const int iyear, double *date, int *mins)
{   
    // The model is not valid for dates prior to 1906.
    if (iyear <= 1906) return 1;
    
    // Convert MDY to modified Julian date
    int mjd = 0, mjd0 = 0, mjd1 = 0, i = 1;
    iymdmj_(&iyear, &month, &iday, &mjd);

    // calculate time in decimal years (date)
    iymdmj_(&iyear, &i, &i, &mjd0);
    int iyear1 = iyear + 1;
    iymdmj_(&iyear1, &i, &i, &mjd1);
    int day = mjd - mjd0;
    int denom = mjd1 - mjd0;
    *date = iyear + (day / denom);
    
    // calculate Julian time in minutes (mins)
    *mins = mjd * 24 * 60;
    
    return 0;
}
