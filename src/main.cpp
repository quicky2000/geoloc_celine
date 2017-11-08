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

#include "my_bmp.h"
#include "quicky_exception.h"
#include "simple_gui.h"
#include "projections.h"
#include "city.h"
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <string>
#include <set>
#include <fstream>
#include <map>
#include <vector>
#include <limits>

//#define LAMBERT_PROJECTION
#define LAMBERT2_PROJECTION

//------------------------------------------------------------------------------
double get_pixel_coordinate(const double & p_input1,
                            const double & p_input2,
                            const double & p_output_1,
                            const double & p_output_2,
                            const double & p_input
                                 )
{
    double l_delta_den = p_input2 - p_input1;
    double l_delta_num = (double)p_output_2 - (double)p_output_1;
    double l_const = ((double)p_output_1) - (l_delta_num / l_delta_den) * p_input1;
    return p_input * (l_delta_num / l_delta_den) + l_const;
}

//------------------------------------------------------------------------------
void get_x_y(const double & p_lat,
             const double & p_lon,
             unsigned int & p_x,
             unsigned int & p_y,
             simple_gui::simple_gui & l_gui
            )
{
    // Brest :     48.4004997828    -4.5027907853
    // Lochrist  : 48.6198454528    -4.1987464419  (48,106)
    // Bray-Dunes; 51.0730436867     2.5276487046  (210,16)
    // Dunkerque : 51.0307229078     2.3375241410
    double l_lochrist_lon = -4.1987464419;
    unsigned int l_lochrist_x = 48;
    double l_bray_dunes_lon = 2.5276487046;
    unsigned int l_bray_dunes_x = 210;

    double l_x_north = get_pixel_coordinate(l_lochrist_lon, l_bray_dunes_lon, l_lochrist_x, l_bray_dunes_x,p_lon);

    double l_lochrist_lat = 48.6198454528;
    unsigned int l_lochrist_y = 106;
    double l_bray_dunes_lat = 51.0730436867;
    unsigned int l_bray_dunes_y = 16;

    double l_y_north = get_pixel_coordinate(l_lochrist_lat, l_bray_dunes_lat, l_lochrist_y, l_bray_dunes_y,p_lat);

    // Argeles   : 42.5352193463     3.0242986289
    // Cerbère   : 42.4439076617     3.1483550073  (239,313)
    // Antibes   : 43.5874651460     7.1063541826  (338,266)
    double l_cerbere_lon = 3.1483550073;
    unsigned int l_cerbere_x = 239;
    double l_antibes_lon = 7.1063541826;
    unsigned int l_antibes_x = 338;

    double l_x_south = get_pixel_coordinate(l_cerbere_lon, l_antibes_lon, l_cerbere_x, l_antibes_x,p_lon);

    double l_cerbere_lat = 42.4439076617;
    unsigned int l_cerbere_y = 313;
    double l_antibes_lat = 43.5874651460;
    unsigned int l_antibes_y = 266;

    double l_y_south = get_pixel_coordinate(l_cerbere_lat, l_antibes_lat, l_cerbere_y, l_antibes_y,p_lat);

    double l_x_west = get_pixel_coordinate(l_lochrist_lon, l_cerbere_lon, l_lochrist_x, l_cerbere_x,p_lon);
    double l_y_west = get_pixel_coordinate(l_lochrist_lat, l_cerbere_lat, l_lochrist_y, l_cerbere_y,p_lat);

    double l_x_east = get_pixel_coordinate(l_bray_dunes_lon, l_antibes_lon, l_bray_dunes_x, l_antibes_x,p_lon);
    double l_y_east = get_pixel_coordinate(l_bray_dunes_lat, l_antibes_lat, l_bray_dunes_y, l_antibes_y,p_lat);


    unsigned int l_blue_code = l_gui.get_color_code(0, 0, 255);
    unsigned int l_black_code = l_gui.get_color_code(128, 128, 128);
    unsigned int l_violet_code = l_gui.get_color_code(255, 0, 255);
    unsigned int l_white_code = l_gui.get_color_code(255, 255, 255);
    int l_size = 2;
    for(int l_offset = -l_size;
        l_offset <= l_size;
        ++l_offset)
    {
        l_gui.set_pixel_without_lock(l_x_east + l_offset,l_y_east,l_blue_code);
        l_gui.set_pixel_without_lock(l_x_east,l_y_east + l_offset,l_blue_code);
        l_gui.refresh();
    }
    for(int l_offset = -l_size;
        l_offset <= l_size;
        ++l_offset)
    {
        l_gui.set_pixel_without_lock(l_x_west + l_offset,l_y_west,l_violet_code);
        l_gui.set_pixel_without_lock(l_x_west,l_y_west + l_offset,l_violet_code);
        l_gui.refresh();
    }
    std::cout << l_x_north << "," << l_y_north << std::endl;
    std::cout << l_x_south << "," << l_y_south << std::endl;
    for(int l_offset = -l_size;
        l_offset <= l_size;
        ++l_offset)
    {
        l_gui.set_pixel_without_lock(l_x_north + l_offset,l_y_north,l_white_code);
        l_gui.set_pixel_without_lock(l_x_north,l_y_north + l_offset,l_white_code);
        l_gui.refresh();
    }
    for(int l_offset = -l_size;
        l_offset <= l_size;
        ++l_offset)
    {
        l_gui.set_pixel_without_lock(l_x_south + l_offset,l_y_south,l_black_code);
        l_gui.set_pixel_without_lock(l_x_south,l_y_south + l_offset,l_black_code);
        l_gui.refresh();
    }
}

//------------------------------------------------------------------------------
void draw_cross(unsigned int p_x,
                unsigned int p_y,
                int p_size,
                simple_gui::simple_gui & p_gui,
                unsigned int p_color
               )
{
    for(int l_offset = -p_size;
        l_offset <= p_size;
        ++l_offset)
    {
        p_gui.set_pixel_without_lock(p_x + l_offset,p_y,p_color);
        p_gui.set_pixel_without_lock(p_x,p_y + l_offset,p_color);
        p_gui.refresh();
    }

}

//------------------------------------------------------------------------------
int main(int argc,char ** argv)
{
  try
  {
      std::string l_map = "map_without_legend.bmp";
      std::string l_simple_map = "simple_map.bmp";
      lib_bmp::my_bmp l_bmp(l_map);
      lib_bmp::my_bmp l_simple_bmp(l_bmp.get_width(),l_bmp.get_height(),24);
      simple_gui::simple_gui l_gui;
      l_gui.create_window(400,400);
      lib_bmp::my_color l_white(255,255,255);
      std::vector<lib_bmp::my_color> l_simple_colors;
      l_simple_colors.push_back(lib_bmp::my_color(255,255,255));
      l_simple_colors.push_back(lib_bmp::my_color(0,0,0));
      l_simple_colors.push_back(lib_bmp::my_color(255,0,0));
      l_simple_colors.push_back(lib_bmp::my_color(255,153,51));
      l_simple_colors.push_back(lib_bmp::my_color(255,255,0));
      l_simple_colors.push_back(lib_bmp::my_color(204,255,204));
      l_simple_colors.push_back(lib_bmp::my_color(102,153,0));
      std::set<lib_bmp::my_color> l_colors;

      // Store min max values of pixel x,y not white
      unsigned int l_min_pix_x = std::numeric_limits<unsigned int>::max();
      unsigned int l_max_pix_x = std::numeric_limits<unsigned int>::min();
      unsigned int l_min_pix_y = std::numeric_limits<unsigned int>::max();
      unsigned int l_max_pix_y = std::numeric_limits<unsigned int>::min();

      std::map<unsigned int, std::pair<unsigned int, unsigned int> > l_line_x_min_max;
      for(unsigned int l_y = 0;
                l_y < l_bmp.get_height();
                ++l_y
         )
      {
          unsigned int l_line_x_min = std::numeric_limits<unsigned int>::max();
          unsigned int l_line_x_max = std::numeric_limits<unsigned int>::min();
            for(unsigned int l_x = 0;
                    l_x < l_bmp.get_width();
                    ++l_x
               )
            {
                lib_bmp::my_color l_color = l_bmp.get_pixel_color(l_x,l_y);
                uint32_t l_min = std::numeric_limits<uint32_t>::max();
                lib_bmp::my_color l_simple_color;
                for(auto l_iter:l_simple_colors)
                {
                    double l_distance = pow(((int32_t)l_color.get_red()) - ((int32_t)l_iter.get_red()),2);
                    l_distance += pow(((int32_t)l_color.get_green()) - ((int32_t)l_iter.get_green()),2);
                    l_distance += pow(((int32_t)l_color.get_blue()) - ((int32_t)l_iter.get_blue()),2);
                    if(((uint32_t)l_distance) < l_min)
                    {
                        l_min = (uint32_t)l_distance;
                        l_simple_color = l_iter;
                    }
                }
                l_simple_bmp.set_pixel_color(l_x,l_y,lib_bmp::my_color_alpha(l_simple_color.get_red(),l_simple_color.get_green(),l_simple_color.get_blue(),0));
                if(l_color != l_white)
                {
                    if(l_x < l_min_pix_x)
                    {
                        l_min_pix_x = l_x;
                    }
                    if(l_x < l_line_x_min)
                    {
                        l_line_x_min = l_x;
                    }
                    if(l_y < l_min_pix_y)
                    {
                        l_min_pix_y = l_y;
                    }
                    if(l_x > l_max_pix_x)
                    {
                        l_max_pix_x = l_x;
                    }
                    if(l_x > l_line_x_max)
                    {
                        l_line_x_max = l_x;
                    }
                    if(l_y > l_max_pix_y)
                    {
                        l_max_pix_y = l_y;
                    }
                    std::cout << "(" << l_x << "," << l_y << ") = " << l_color << std::endl;
                    l_colors.insert(l_color);
                    uint8_t l_red = l_color.get_red();
                    uint8_t l_green = l_color.get_green();
                    uint8_t l_blue = l_color.get_blue();
                    unsigned int l_color = l_gui.get_color_code(l_red, l_green, l_blue);
                    l_gui.set_pixel_without_lock(l_x,l_y,l_color);
                }
            }
          l_line_x_min_max.insert(std::map<unsigned int,std::pair<unsigned int,unsigned int> >::value_type(l_y,std::pair<unsigned int,unsigned int>(l_line_x_min,l_line_x_max)));
        }
      l_gui.refresh();
      l_simple_bmp.save(l_simple_map);
      std::cout << "Nb colors : " << l_colors.size() << std::endl;
      std::cout << "Min x = " << l_min_pix_x << std::endl;
      std::cout << "Max x = " << l_max_pix_x << std::endl;
      std::cout << "Min y = " << l_min_pix_y << std::endl;
      std::cout << "Max y = " << l_max_pix_y << std::endl;

//      unsigned int l_x= get_pixel_coordinate(l_lochrist_lon, l_bray_dunes_lon, l_lochrist_x, l_bray_dunes_x,l_bray_dunes_lon);
//      std::cout << "l_x = " << l_x << std::endl;
      unsigned int l_x;
      unsigned int l_y;
      //Feurget_x_y(45.7323040855,4.2232415379,l_x,l_y,l_gui);
      // Grenoble get_x_y(45.1821215167,5.7213305175,l_x,l_y,l_gui);
      // Marseille get_x_y(43.2999009436,5.3822786980,l_x,l_y,l_gui);
      // Calais get_x_y(50.9502072754,1.8757556613,l_x,l_y,l_gui);
      // Biarritz get_x_y(43.4695847227,-1.5530985752,l_x,l_y,l_gui);
      // Perpignan get_x_y(42.6965954131,2.8993695398,l_x,l_y,l_gui);
      // Toulon
      get_x_y(43.1361589728,5.9323963425,l_x,l_y,l_gui);
      double l_dx_lochrist;
      double l_dy_lochrist;
      WGS84ToLambert93(48.6198454528,-4.1987464419,l_dx_lochrist,l_dy_lochrist);
      std::cout << l_dx_lochrist << "," << l_dy_lochrist << std::endl;
      double l_dx_bray_dunes;
      double l_dy_bray_dunes;
      // Antibes   : 43.5874651460     7.1063541826  (338,266)
      WGS84ToLambert93(43.5874651460,7.1063541826,l_dx_bray_dunes,l_dy_bray_dunes);
      std::cout << l_dx_bray_dunes << "," << l_dy_bray_dunes << std::endl;
      double l_dx_feurs;
      double l_dy_feurs;
      WGS84ToLambert93(42.6965954131,2.8993695398,l_dx_feurs,l_dy_feurs);
      std::cout << l_dx_feurs << "," << l_dy_feurs << std::endl;
      unsigned int l_x_feurs = get_pixel_coordinate(l_dx_lochrist,l_dx_bray_dunes,48,338,l_dx_feurs);
      unsigned int l_y_feurs = get_pixel_coordinate(l_dy_lochrist,l_dy_bray_dunes,106,266,l_dy_feurs);
      std::cout << l_x_feurs << std::endl;
      std::cout << l_y_feurs << std::endl;
      int l_size = 10;

      unsigned int l_red_code = l_gui.get_color_code(255, 0, 0);
      for(int l_offset = -l_size;
          l_offset <= l_size;
          ++l_offset)
      {
          l_gui.set_pixel_without_lock(l_x_feurs + l_offset,l_y_feurs,l_red_code);
          l_gui.set_pixel_without_lock(l_x_feurs,l_y_feurs + l_offset,l_red_code);
          l_gui.refresh();
      }

      std::string l_csv_file("liste_villes_2017.csv");
      std::ifstream l_file;
      l_file.open(l_csv_file);
      if(!l_file.is_open())
      {
          throw quicky_exception::quicky_runtime_exception("Unable to open "+l_csv_file,__LINE__,__FILE__);
      }
      std::string l_line;
      std::map<unsigned int,city> l_cities;
      double l_lat_min = std::numeric_limits<double>::max();
      double l_lat_max = std::numeric_limits<double>::min();
      double l_lon_min = std::numeric_limits<double>::max();
      double l_lon_max = std::numeric_limits<double>::min();
      unsigned int l_line_number = 1;
      std::set<std::string> l_ignored_cities =
              {
                      "",
                      "Ville",
                      "Ouessant",
                      "Île-de-Sein",
                      "Île-Molène",
                      "Locmaria",
                      "Sainte-Marie-de-Ré",
                      "Couarde-sur-Mer",
                      "Groix",
                      "Île-d'Yeu",
                      "Grand-Village-Plage",
                      "Saint-Martin-de-Ré",
                      "Saint-Georges-d'Oléron",
                      "Saint-Clément-des-Baleines",
                      "Palais",
                      "Bois-Plage-en-Ré",
                      "Sauzon",
                      "Noirmoutier-en-l'Île",
                      "Bangor",
                      "Île-de-Batz",
                      "Ars-en-Ré",
                      "Flotte",
                      "Île-d'Houat",
                      "Saint-Pierre-d'Oléron",
                      "Saint-Denis-d'Oléron",
                      "Loix",
                      "Brée-les-Bains",
                      "Dolus-d'Oléron",
                      "Guérinière",
                      "Hoedic",
                      "Portes-en-Ré"
              };
      std::set<unsigned int> l_ignored_lines =
              {
                      16587 // Épine
              };
      while(!l_file.eof())
      {
          std::getline(l_file,l_line);
          if(l_line.size())
          {
              size_t l_pos = 0;
              size_t l_previous_pos = 0;
              std::string l_name;
              double l_lat;
              double l_lon;
              unsigned int l_index = 0;
              do
              {
                  l_previous_pos = l_pos;
                  l_pos = l_line.find(';',l_previous_pos);
                  std::string l_substring(l_line.substr(l_previous_pos, l_pos - l_previous_pos));
                  switch(l_index)
                  {
                      case 2:
                          l_name = l_substring;
                          break;
                      case 4:
                          if("" != l_substring)
                          {
                              l_lat = std::stod(l_substring);
                          }
                          break;
                      case 5:
                          if("" != l_substring)
                          {
                              l_lon = std::stod(l_substring);
                          }
                          break;
                  }
                  ++l_index;
                  ++l_pos;
              }while(l_index < 7 && std::string::npos != l_pos);
              // Corsica cities have no name
              if(l_ignored_cities.end() == l_ignored_cities.find(l_name) && !l_ignored_lines.count(l_line_number))
              {
                  l_cities.insert(std::map<unsigned int,city>::value_type(l_line_number,city(l_line_number,l_name,l_lon,l_lat)));
                  if(l_lat < l_lat_min)
                  {
                      l_lat_min = l_lat;
                  }
                  if(l_lat > l_lat_max)
                  {
                      l_lat_max = l_lat;
                  }
                  if(l_lon < l_lon_min)
                  {
                      l_lon_min = l_lon;
                  }
                  if(l_lon > l_lon_max)
                  {
                      l_lon_max = l_lon;
                  }
              }
              ++l_line_number;
          }
      }
      l_file.close();

      std::cout << "Min lon = " << l_lon_min << std::endl;
      std::cout << "Max lon = " << l_lon_max << std::endl;
      std::cout << "Min lat = " << l_lat_min << std::endl;
      std::cout << "Max lat = " << l_lat_max << std::endl;

      double l_lambert_lat_min = std::numeric_limits<double>::max();
      double l_lambert_lat_max = std::numeric_limits<double>::min();
      double l_lambert_lon_min = std::numeric_limits<double>::max();
      double l_lambert_lon_max = std::numeric_limits<double>::min();

      for(auto l_iter:l_cities)
      {
          double l_lon = l_iter.second.get_lon();
          double l_lat = l_iter.second.get_lat();
          double l_lambert_lon;
          double l_lambert_lat;
#ifdef LAMBERT2_PROJECTION
          WGS84ToLambert2e(l_lon,l_lat,l_lambert_lon,l_lambert_lat);
#else
          WGS84ToLambert93(l_lat,l_lon,l_lambert_lon,l_lambert_lat);
#endif
          if(l_lambert_lat < l_lambert_lat_min)
          {
              l_lambert_lat_min = l_lambert_lat;
          }
          if(l_lambert_lat > l_lambert_lat_max)
          {
              l_lambert_lat_max = l_lambert_lat;
          }
          if(l_lambert_lon < l_lambert_lon_min)
          {
              l_lambert_lon_min = l_lambert_lon;
          }
          if(l_lambert_lon > l_lambert_lon_max)
          {
              l_lambert_lon_max = l_lambert_lon;
          }
      }

      std::cout << "Lambert Min lon = " << l_lambert_lon_min << std::endl;
      std::cout << "Lambert Max lon = " << l_lambert_lon_max << std::endl;
      std::cout << "Lambert Min lat = " << l_lambert_lat_min << std::endl;
      std::cout << "Lambert Max lat = " << l_lambert_lat_max << std::endl;

      for(auto l_iter:l_cities)
      {
          double l_lon = l_iter.second.get_lon();
          double l_lat = l_iter.second.get_lat();
          double l_dx;
          double l_dy;
#ifdef LAMBERT_PROJECTION
#ifdef LAMBERT2_PROJECTION
          WGS84ToLambert2e(l_lon,l_lat,l_dx,l_dy);
#else
          WGS84ToLambert93(l_lat,l_lon,l_dx,l_dy);
#endif
          unsigned int l_y = get_pixel_coordinate(l_lambert_lat_max,l_lambert_lat_min,l_min_pix_y,l_max_pix_y,l_dy);
//          auto l_iter_min_max = l_line_x_min_max.find(l_y);
//          assert(l_line_x_min_max.end() != l_iter_min_max);
//          l_min_pix_x = l_iter_min_max->second.first;
//          l_max_pix_x = l_iter_min_max->second.second;
          unsigned int l_x = get_pixel_coordinate(l_lambert_lon_min,l_lambert_lon_max,l_min_pix_x,l_max_pix_x,l_dx);
#else
          unsigned int l_y = get_pixel_coordinate(l_lat_max,l_lat_min,l_min_pix_y,l_max_pix_y,l_lat);
//          auto l_iter_min_max = l_line_x_min_max.find(l_y);
//          assert(l_line_x_min_max.end() != l_iter_min_max);
//          l_min_pix_x = l_iter_min_max->second.first;
//          l_max_pix_x = l_iter_min_max->second.second;
          unsigned int l_x = get_pixel_coordinate(l_lon_min,l_lon_max,l_min_pix_x,l_max_pix_x,l_lon);
#endif
//          uint32_t l_color_code = l_gui.get_pixel(l_x,l_y);
//          assert(l_color_code != l_gui.get_color_code(0,0,0));

          unsigned int l_blue = l_gui.get_color_code(0, 0, 255);
          l_gui.set_pixel_without_lock(l_x,l_y,l_blue);
      }
      l_gui.refresh();
      sleep(100);

  }
  catch(quicky_exception::quicky_runtime_exception & e)
    {
      std::cout << "ERROR : " << e.what() << " " << e.get_file() << ":" << e.get_line() << std::endl ;
      return(-1);
    }
  catch(quicky_exception::quicky_logic_exception & e)
    {
      std::cout << "ERROR : " << e.what() << e.get_file() << ":" << e.get_line() <<  std::endl ;
      return(-1);
    }
  return 0;
  
}

//EOF
