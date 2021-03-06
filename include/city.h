/*    This file is part of geoloc_celine
      Copyright (C) 2017  Julien Thevenon ( julien_thevenon at yahoo.fr )

      This program is free software: you can redistribute it and/or modify
      it under the terms of the GNU General Public License as published by
      the Free Software Foundation, either version 3 of the License, or
      (at your option) any later version.

      This program is distributed in the hope that it will be useful,
      but WITHOUT ANY WARRANTY; without even the implied warranty of
      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
      GNU General Public License for more details.

      You should have received a copy of the GNU General Public License
      along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

#ifndef GEOLOC_CELINE_CITY_H
#define GEOLOC_CELINE_CITY_H

#include <string>

class city
{
  public:
    city(unsigned int p_line,
         const std::string & p_name,
         const double & p_lon,
         const double & p_lat,
         const std::string & p_postal_code,
         const std::string & p_city_code,
         const std::string & p_department
        );

    unsigned int
    get_line() const;

    const std::string &
    get_name() const;

    const double &
    get_lon() const;

    const double &
    get_lat() const;

    const std::string &
    get_postal_code() const;

    const std::string &
    get_city_code() const;

    const std::string &
    get_department() const;

  private:
    unsigned int m_line;
    std::string m_name;
    double m_lon;
    double m_lat;
    std::string m_postal_code;
    std::string m_city_code;
    std::string m_department;
};

//-----------------------------------------------------------------------------
city::city(unsigned int p_line,
           const std::string & p_name,
           const double & p_lon,
           const double & p_lat,
           const std::string & p_postal_code,
           const std::string & p_city_code,
           const std::string & p_department
          ):
    m_line(p_line),
    m_name(p_name),
    m_lon(p_lon),
    m_lat(p_lat),
    m_postal_code(p_postal_code),
    m_city_code(p_city_code),
    m_department(p_department)
{

}

//-----------------------------------------------------------------------------
unsigned int
city::get_line() const
{
    return m_line;
}

//-----------------------------------------------------------------------------
const std::string &
city::get_name() const
{
    return m_name;
}

//-----------------------------------------------------------------------------
const double &
city::get_lon() const
{
    return m_lon;
}

//-----------------------------------------------------------------------------
const double &
city::get_lat() const
{
    return m_lat;
}

//-----------------------------------------------------------------------------
const std::string &
city::get_postal_code() const
{
    return m_postal_code;
}

//-----------------------------------------------------------------------------
const std::string &
city::get_city_code() const
{
    return m_city_code;
}

//-----------------------------------------------------------------------------
const std::string &
city::get_department() const
{
    return m_department;
}

#endif //GEOLOC_CELINE_CITY_H
