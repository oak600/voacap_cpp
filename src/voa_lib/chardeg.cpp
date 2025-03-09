#include <string>

/** Converts a lat/long string to degrees.
 *
 * @param cardCStr A C string representing the lat/long in the forms
 *             dddXmm'ss", dddXmm, dddX, ddd.dddX, ddd.ddd, or -ddd.ddd
 *             where  X=N or S for latitude
 *                     =E or W for longitude
 *             W longitude is negative
 *             E latitude  is negative
 *
 * @param latLon 0 = latitude, 1 = longitude, anything else = error return
 * @param deg A reference to a double to return the data to
 * @param ierr A reference to an int32_teger to return function success.
 *
 * @returns None, uses legacy return architecture...see params deg and ierr
 *
 * @note This is a direct conversion of CHARDEG in Fortran to C++. This will need to be reworked in the future
 *       to optimize and use full C++ feature set.
 */
void charDeg(const char *cardCStr, int32_t latlon, double &deg, int32_t &ierr)
{
      std::string card;
      try
      {
            card = std::string(cardCStr);
      } catch(...)
      {
            ierr = -9; //BAD VALUE
      }

      const double DMAX[2] = {90.0, 360.0}; // Max values for LAT/LON
      ierr = 0;
      deg = 0.0;
      int32_t min = 0;
      int32_t isec = 0;
      int32_t isign = 1;
      int32_t idec = 0;
      int32_t idel = 0;
      bool ionce = false;
      size_t ii = 0;

      for (size_t i = 0; i < card.length(); ++i)
      {
            char ich = card[i];

            if ((!ionce) && (ich == ' '))
            {
                  continue; // Skip leading spaces
            }
            ionce = true;

            if ((ich == ' ') || (ich == ','))
            {
                  break; // End of string
            }

            if (idel == 3)
            {
                  ierr = -2; // Something after seconds
                  return;
            } 

            if (ich == '-')
            {
                  if (ii != 0)
                  {
                        ierr = -3; // '-' must be first
                        return;
                  } 
                  isign = -1;
                  continue;
            }

            if (ich == '.')
            {
                  if (idec != 0)
                  {
                        ierr = -4; // Only one decimal point allowed
                        return;
                  } 

                  if (idel != 0)
                  {
                        ierr = -5; // Must be before delimiters
                        return;
                  } 

                  idec = 1;
                  continue;
            }

            if (std::isdigit(ich))
            {
                  int32_t n = ich - '0';
                  if (idel == 0)
                  {
                        if (idec == 0)
                        {
                              deg = deg * 10.0 + n;
                        }
                        else
                        {
                              deg += n / pow(10.0, idec);
                              idec++;
                        }
                  }
                  else if (idel == 1)
                  {
                        min = min * 10 + n;
                        if (min >= 60)
                        {
                              ierr = -6;
                              return;
                        }
                  }
                  else if (idel == 2)
                  {
                        isec = isec * 10 + n;
                        if (isec >= 60)
                        {
                              ierr = -7;
                              return;
                        }
                  }
                  continue;
            }

            idel++;
            if (idel == 1)
            {
                  // Direction check
                  if (latlon == 0)
                  {
                        if ((ich == 'N') || (ich == 'n'))
                        {
                              continue;
                        }  
                        if ((ich == 'S') || (ich == 's'))
                        {
                              isign = -isign;
                              continue;
                        }
                        ierr = -10;
                        return;
                  }
                  else
                  {
                        if ((ich == 'E') || (ich == 'e'))
                        {
                              continue;
                        }  
                        if ((ich == 'W') || (ich == 'w'))
                        {
                              isign = -isign;
                              continue;
                        }
                        ierr = -11;
                        return;
                  }
            }
      }

      deg += (static_cast<double>(min) / 60.0) + (static_cast<double>(isec) / 3600.0);
      if (deg > DMAX[latlon])
      {
            ierr = -9; //Over size
            return;
      }
      if (latlon == 1 && deg > 180.0)
      {
            deg -= 360.0;
      }
            
      deg *= isign;
}