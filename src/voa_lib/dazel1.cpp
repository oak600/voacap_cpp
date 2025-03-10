/** Computes end point given distance and bearing.

C#  SUB DAZEL1                 Same as CALL DAZEL(1)
       IMPLICIT DOUBLE PRECISION(A-L,N-Y)
C
C     TWO MODES--   0   INPUT LAT AND LON OF END POINT
C                       RETURN DISTANCE AND AZIMUTH TO END PT WITH ELEVATIONS
C                   1   INPUT BEARING (AZIMUTH) OF END POINT
C                       RETURN LAT AND LON OF END POINT WITH ELEVATIONS
C
C   MODE 0
C   INPUT PARAMETERS (THESE DEFINE LOCATION OF POINTS T (TRANSMITTER)
C     AND R (RECEIVER) RELATIVE TO A SPHERICAL EARTH.
C     ZTLAT - 
C     ZTLON - 
C     ZTHT  - 
C     ZRLAT - 
C     ZRLON - 
C     ZRHT  - 
C
C   OUTPUT PARAMETERS
C     ZTAZ  - 
C     ZRAZ  - 
C     ZTELV - 
C     ZRELV - 
C     ZTAKOF - 
C     ZRAKOF - 
C     ZD    - 
C     ZDGC  - 
C
C   MODE 1
C   INPUT PARAMETERS                    OUTPUT PARAMETERS
C     ZTLAT                                ZRLAT
C     ZTLON                                ZRLON
C     ZTAZ                                 RELEV,ZRAKOF
C     ZDGC                                 TELEV,ZTAKOF
C
C
C     ALL OF THE ABOVE PARAMETERS START WITH THE LETTER Z AND ARE SINGLE
C     PRECISION.  ALL PROGRAM VARIABLES ARE DOUBLE PRECISION.
C     PROGRAM IS UNPREDICTABLE FOR SEPARATIONS LESS THAN 0.00005 DEGREES,
C     ABOUT 5 METERS.
C
C   WRITTEN BY KEN SPIES 5/79
C   REFRACTION AND ST. LINE ELEVATIONS BY EJH
C
      COMMON/AZEL/ ZTLAT,ZTLON,ZTHT,ZRLAT,ZRLON,ZRHT,ZTAZ,ZRAZ,
     * ZTELV,ZRELV,ZD,ZDGC,ZTAKOF,ZRAKOF
      DATA PI/3.141592653589793238462643D0/,RERTH/6370.0D0/
      DATA DEGREES_TO_RADIANS/0.01745329252D0/,RTOD/57.29577951D0/
C
C     COMPUTE END POINT GIVEN DISTANCE AND BEARING
C
*/


#include <cmath>
#include <vector>

// Define constants
const double RADIUS_EARTH_KM = 6370.0; //radius in earth in km, todo can we update this to 6378?
constexpr double RADIUS_4_3_KM = RADIUS_EARTH_KM * 4.0 / 3.0;
constexpr double DEGREES_TO_RADIANS = M_PI / 180.0; //degrees to radians
constexpr double RADIANS_TO_DEGREES = 180.0 / M_PI;

/** Structure for AZEL common block variables
 * 
 * @note The original Fortran used a mixture of floats and doubles. Moving to all doubles for this iteration.
 */ 
typedef struct DAZEL1_Params {
    double tLatDegrees; // LATITUDE (DECIMAL DEGREES NORTH OF EQUATOR) OF POINT T
    double tLongDegrees; // LONGITUDE (DECIMAL DEGREES EAST OF PRIME (GREENWICH) MERIDIAN OF POINT T
    double tElevSeaLvlMeters;  // HEIGHT (METERS ABOVE MEAN SEA LEVEL) OF POINT T
    double rLatDegrees; // LATITUDE (DECIMAL DEGREES NORTH OF EQUATOR) OF POINT R
    double rLongDegrees; // LONGITUDE (DECIMAL DEGREES EAST OF PRIME MERIDIAN OF POINT R
    double rElevSeaLvlMeters;  // HEIGHT (METERS ABOVE MEAN SEA LEVEL) OF POINT R
    double tAzimuthOfRDegrees;  // AZUMUTH (DECIMAL DEGREES CLOCKWISE FROM NORTH) AT T OF R
    double rAzimuthOfTDegrees;  // AZIMUTH (DECIMAL DEGREES CLOCKWISE FROM NORTH) AT R OF T
    double tElevationAngleDegrees; // ELEVATION ANGLE (DECIMAL DEGREES ABOVE HORIZONTAL AT T OF STRAIGHT LINE BETWEEN T AND R
    double rElevationAngleDegrees; // ELEVATION ANGLE (DECIMAL DEGREES ABOVE HORIZONTAL AT R OF STRAIGHT LINE BETWEEN T AND R
    double trDistanceKilometers;    // STRAIGHT LINE DISTANCE (KILOMETERS) BETWEEN T AND R
    double trGreatCircleDistanceKilometers;  // GREAT CIRCLE DISTANCE (KILOMETERS) BETWEEN T AND R
    double tTakeoffAngleDegrees; // TAKE-OFF ANGLE (DECIMAL DEGREES ABOVE HORIZONTAL AT T OF REFRACTED RAY BETWEEN T AND R (ASSUMED 4/3 EARTH RADIUS)
    double rTakeoffAngleDegrees; // TAKE-OFF ANGLE (DECIMAL DEGREES ABOVE HORIZONTAL AT R OF REFRACTED RAY BETWEEN T AND R (ASSUMED 4/3 EARTH RADIUS)
} DAZEL1_Params;

void DAZEL1(DAZEL1_Params &azel) {
    const double tLongRadians = azel.tLongDegrees * DEGREES_TO_RADIANS; //todo, unused?
    const double tAzimuthRadians = azel.tAzimuthOfRDegrees * DEGREES_TO_RADIANS;
    double greatCircleFraction = azel.trGreatCircleDistanceKilometers / RADIUS_EARTH_KM; //Great Circle fraction
    const double COLAT = M_PI_2 - (azel.tLatDegrees * DEGREES_TO_RADIANS);
    const double COSCO = cos(COLAT);
    const double SINCO = sin(COLAT);
    const double COSGC = cos(greatCircleFraction);
    const double SINGC = sin(greatCircleFraction);

    const double COSB = COSCO * COSGC + SINCO * SINGC * cos(azel.tAzimuthOfRDegrees * DEGREES_TO_RADIANS);

    const double B = atan2(sqrt(std::max(0.0, (1.0 - COSB * COSB))), COSB);
    const double ARC = (COSGC - COSCO * COSB) / (SINCO * sin(B));

    const double RDLON = atan2(sqrt(std::max(0.0, (1.0 - ARC * ARC))), ARC);
    azel.rLatDegrees = (M_PI_2 - fabs(B)) * RADIANS_TO_DEGREES;

    const double DRLAT = azel.rLatDegrees;
    azel.rLatDegrees = copysign(DRLAT, COSB);
    azel.rLongDegrees = azel.tLongDegrees + (fabs(RDLON) * RADIANS_TO_DEGREES);
    if (azel.tAzimuthOfRDegrees > 180) 
    {
       azel.rLongDegrees = azel.tLongDegrees - (fabs(RDLON) * RADIANS_TO_DEGREES);
    }
    
    double tElevSeaLvlKilometers = azel.tElevSeaLvlMeters * 1.0E-3; /// Height above sea level in km of point T
    double rElevationSeaLvlKilometers = azel.rElevSeaLvlMeters * 1.0E-3; /// Height above sea level in km of point R
    double deltaElevationKilometers = rElevationSeaLvlKilometers - tElevSeaLvlKilometers; /// Difference in height in km between R and T
    double sinGreatCircle = sin(0.5 * greatCircleFraction);
    double D = sqrt((deltaElevationKilometers * deltaElevationKilometers) + 
                   4.0 * (RADIUS_EARTH_KM + tElevSeaLvlKilometers) * (RADIUS_EARTH_KM + rElevationSeaLvlKilometers) * sinGreatCircle * sinGreatCircle);
    
    double aElevSeaLvlKilometers, bElevSeaLvlKilometers;
    if (deltaElevationKilometers >= 0) {
        aElevSeaLvlKilometers = tElevSeaLvlKilometers;
        bElevSeaLvlKilometers = rElevationSeaLvlKilometers;
    } else {
        aElevSeaLvlKilometers = rElevationSeaLvlKilometers;
        bElevSeaLvlKilometers = tElevSeaLvlKilometers;
    }
    
    //Compute take off angles assuming 4/3 earth radius
    double SAELV = 0.5 * (D * D + fabs(deltaElevationKilometers) * (RADIUS_EARTH_KM + aElevSeaLvlKilometers + RADIUS_EARTH_KM + bElevSeaLvlKilometers)) / (D * (RADIUS_EARTH_KM + aElevSeaLvlKilometers));
    double AELV = atan2(SAELV, sqrt(std::max(0.0001, (1.0 - SAELV * SAELV))));
    double BELV = (AELV - greatCircleFraction) * RADIANS_TO_DEGREES;
    AELV = -AELV * RADIANS_TO_DEGREES;
    
    greatCircleFraction = 0.75 * greatCircleFraction;
    sinGreatCircle = sin(0.5 * greatCircleFraction);
    double P = 2.0 * sinGreatCircle * sinGreatCircle;
    double AALT = RADIUS_4_3_KM + aElevSeaLvlKilometers;
    double BALT = RADIUS_4_3_KM + bElevSeaLvlKilometers;
    double DA = sqrt(deltaElevationKilometers * deltaElevationKilometers + 2.0 * AALT * BALT * P);
    SAELV = 0.5 * (DA * DA + fabs(deltaElevationKilometers) * (AALT + BALT)) / (DA * AALT);

    double ATAKOF = atan(SAELV / sqrt(std::max(0.0001, (1.0 - SAELV * SAELV))));
    double BTAKOF = (ATAKOF - greatCircleFraction) * RADIANS_TO_DEGREES;
    ATAKOF = -ATAKOF * RADIANS_TO_DEGREES;
    
    if (deltaElevationKilometers < 0) {
        azel.tElevationAngleDegrees = AELV;
        azel.rElevationAngleDegrees = BELV;
        azel.tTakeoffAngleDegrees = ATAKOF;
        azel.rTakeoffAngleDegrees = BTAKOF;
    } else {
        azel.tElevationAngleDegrees = BELV;
        azel.rElevationAngleDegrees = AELV;
        azel.tTakeoffAngleDegrees = BTAKOF;
        azel.rTakeoffAngleDegrees = ATAKOF;
    }
}

