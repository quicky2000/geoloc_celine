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

#ifndef GEOLOC_CELINE_GENERATE_CSV_H
#define GEOLOC_CELINE_GENERATE_CSV_H

#include "quicky_exception.h"
#include "city.h"
#include <string>
#include <fstream>
#include <iomanip>

class generate_csv
{
  public:
    inline generate_csv(const std::string & p_file_name);
    inline void add_city(const city & p_city,
                         const unsigned int & p_red,
                         const unsigned int & p_green,
                         const unsigned int & p_blue,
                         const std::string & p_legend
                        );
    inline ~generate_csv(void);
  private:
    inline void
    print_color_componant(unsigned int p_componant);
    std::ofstream m_stream;
};

//-----------------------------------------------------------------------------
generate_csv::generate_csv(const std::string & p_file_name):
    m_stream(nullptr)
{
    m_stream.open(p_file_name);
    if(!m_stream.is_open())
    {
        throw quicky_exception::quicky_runtime_exception("Unable to create file \"" + p_file_name + "\"", __LINE__, __FILE__);
    }
    m_stream << "code commune;";
    m_stream << "code postale;";
    m_stream << "Ville;";
    m_stream << "latitude;";
    m_stream << "longitude;";
    m_stream << "departement;";
    m_stream << "Couleur;"; // Red
    m_stream << "Legende" << std::endl;
}

//-----------------------------------------------------------------------------
generate_csv::~generate_csv(void)
{
    m_stream.close();
}

//-----------------------------------------------------------------------------
void
generate_csv::add_city(const city & p_city,
                       const unsigned int & p_red,
                       const unsigned int & p_green,
                       const unsigned int & p_blue,
                       const std::string & p_legend
                      )
{
    m_stream << p_city.get_city_code() << ";";
    m_stream << p_city.get_postal_code() << ";";
    m_stream << p_city.get_name() << ";";
    m_stream << p_city.get_lat() << ";";
    m_stream << p_city.get_lon() << ";";
    m_stream << p_city.get_department() << ";";
    m_stream << "#" ;
    print_color_componant(p_red);
    print_color_componant(p_green);
    print_color_componant(p_blue);
    m_stream << ";" ;
    m_stream << "\"" << p_legend << "\"" << std::endl;
}

//-----------------------------------------------------------------------------
void
generate_csv::print_color_componant(unsigned int p_componant)
{
    size_t l_width = m_stream.width();
    m_stream.width(2);

    char l_fill = m_stream.fill();
    m_stream.fill('0');

    m_stream << std::hex << p_componant;
    // Restablish stream characterisitcs
    m_stream << std::dec << std::setfill(l_fill) << std::setw(l_width);
}

#endif //GEOLOC_CELINE_GENERATE_CSV_H
// EOF
