#ifndef MATH_CONVERSIONS_H
#define	MATH_CONVERSIONS_H

namespace liver
{
    static double deg_to_rad(double deg)
    {
        return deg * PI / 180.;
    }

    static double rad_to_deg(double rad)
    {
        return rad / PI * 180.;
    }
}

#endif
