/*    This file is part of quicky_utils
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
#include <iostream>
#include <string>
#include <set>
#include "simple_gui.h"
#include "unistd.h"
#include <math.h>
#include <fstream>
#include <map>
#include <vector>
#include <gmpxx.h>

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

    void WGS84toLambert2e(double d_long,
                          double m_long,
                          double s_long,
                          char orientation_long,
                          double d_lat,
                          double m_lat,
                          double s_lat,
                          char orientation_lat,
                          double & p_x,
                          double & p_y
                         )
    {
        const double pi = acos(-1);
        double lambda_w, phi_w;

        /**************************************************************************************************************/
        /**        0) degres-minutes-secondes + orientation (d,m,s,o) -> radian                                           **/
        /**************************************************************************************************************/

        // Pour la longitude
        lambda_w = d_long + m_long / 60 + s_long / 3600;
        if (orientation_long == 'W') lambda_w = -1 * lambda_w; // Le système de coordonnées géographiques utilisé est postif vers le Nord et vers l'Est

        lambda_w = lambda_w * pi / 180;

        // Pour la latitude
        phi_w = d_lat + m_lat / 60 + s_lat / 3600;
        if (orientation_lat == 'S') phi_w = -1 * phi_w;          // Le système de coordonnées géographiques utilisé est postif vers le Nord et vers l'Est

        phi_w = phi_w * pi / 180;

        /**************************************************************************************************************/
        /**        1) coordonnées géographiques WGS84 (phi_w,lambda_w) -> coordonnées cartésiennes WGS84 (X_w,Y_w,Z_w)  **/
        /**************************************************************************************************************/

        // J'ai utilisé 2 formules données par l'IGN dans un de leur document ainsi que deux constantes de
        // l'ellipsoide de référence du système WGS84 (les deux demi-axes) :

        double a_w = 6378137.0;
        double b_w = 6356752.314;

        // d'où
        double e2_w = (a_w * a_w - b_w * b_w) / (a_w * a_w);

        // et on a la grande normale de l'ellipsoide WGS84
        double N = a_w / sqrt(1 - e2_w * pow(sin(phi_w), 2));

        // ainsi on obtient
        double X_w = N * cos(phi_w) * cos(lambda_w);
        double Y_w = N * cos(phi_w) * sin(lambda_w);
        double Z_w = N * (1 - e2_w) * sin(phi_w);

        /**************************************************************************************************************/
        /**        2) coordonnées cartésiennes WGS84 (X_w,Y_w,Z_w) -> coordonnées cartésiennes NTF (X_n,Y_n,Z_n)          **/
        /**************************************************************************************************************/

        // Là il n'y a qu'un translation à effectuer :

        // on a donc
        double dX = 168.0;
        double dY = 60.0;
        double dZ = -320.0;

        // et...
        double X_n = X_w + dX;
        double Y_n = Y_w + dY;
        double Z_n = Z_w + dZ;

        /**************************************************************************************************************/
        /**        3) coordonnées cartésiennes NTF (X_n,Y_n,Z_n) -> coordonnées géographiques NTF (phi_n,lambda_n)      **/
        /**************************************************************************************************************/

        // J'ai utilisé 1 formule donnée par l'IGN toujours dans le même document ainsi que deux constantes de l'ellipsoide
        // de référence du système NTF, Clarke 1880 :

        double a_n = 6378249.2;
        double b_n = 6356515.0;

        // d'où
        double e2_n = (a_n * a_n - b_n * b_n) / (a_n * a_n);

        // on définit une tolérance de convergence
        double epsilon = pow(10, -10);

        // puis on amorce une boucle de calcul
        double p0 = atan(Z_n / sqrt(X_n * X_n + Y_n * Y_n) * (1 - (a_n * e2_n) / (sqrt(X_n * X_n + Y_n * Y_n + Z_n * Z_n))));
        double p1 = atan((Z_n / sqrt(X_n * X_n + Y_n * Y_n)) / (1 - (a_n * e2_n * cos(p0)) / (sqrt((X_n * X_n + Y_n * Y_n) * (1 - e2_n * pow(sin(p0), 2))))));

        while (!(abs(p1 - p0) < epsilon))
        {

            p0 = p1; p1 = atan((Z_n / sqrt(X_n * X_n + Y_n * Y_n)) / (1 - (a_n * e2_n * cos(p0)) / (sqrt((X_n * X_n + Y_n * Y_n) * (1 - e2_n * pow(sin(p0), 2))))));

        }

        double phi_n = p1;
        double lambda_n = atan(Y_n / X_n);

        /**********************************************************************************************************************/
        /**        4) coordonnées géographiques NTF (phi_n,lambda_n)  coordonnées projetées en Lambert II étendu (X_l2e, Y_l2e) **/
        /**********************************************************************************************************************/

        // J'utilise les formules de projection et les constantes fournies par l'IGN dans un autre document :

        // avant tout on définit quelques constantes
        double n = 0.7289686274;
        double c = 11745793.39;
        double Xs = 600000.0;
        double Ys = 8199695.768;

        double e_n = sqrt(e2_n);
        double lambda0 = 0.04079234433198;   //correspond à la longitude en radian de Paris (2°20'14.025" E) par rapport à Greenwich
                                                // puis on calcule la latitude isométrique
        double L = log(tan(pi / 4 + phi_n / 2) * pow(((1 - e_n * sin(phi_n)) / (1 + e_n * sin(phi_n))), (e_n / 2)));

        // enfin on projette

        double X_l2e = Xs + c * exp((-n * L)) * sin(n * (lambda_n - lambda0));
        double Y_l2e = Ys - c * exp((-n * L)) * cos(n * (lambda_n - lambda0));

        p_x = X_l2e;
        p_y = Y_l2e;
    }

    void WGS84ToLambert2e(double longitude, double latitude, double & p_x, double & p_y)
    {
        // Le système de coordonnées géographiques utilisé est postif vers le Nord et vers l'Est
        char orientation_lat = 'N';
        if (latitude < 0)
        {
            orientation_lat = 'S';
            latitude = -1 * latitude;
        }
        double degree_lat = trunc(latitude);
        double decPart = latitude - degree_lat;

        double multipliedBy3600 = decPart * 3600;
        double dividedBy60 = multipliedBy3600 / 60.0;
        double minute_lat = trunc(dividedBy60);
        decPart = dividedBy60 - minute_lat;
        double seconde_lat = decPart * 60;

        // Le système de coordonnées géographiques utilisé est postif vers le Nord et vers l'Est
        char orientation_long = 'E';
        if (longitude < 0)
        {
            orientation_long = 'W';
            longitude = -1 * longitude;
        }
        double degree_long = trunc(longitude);
        decPart = longitude - degree_long;

        multipliedBy3600 = decPart * 3600;
        dividedBy60 = multipliedBy3600 / 60.0;
        double minute_long = trunc(dividedBy60);
        decPart = dividedBy60 - minute_long;
        double seconde_long = decPart * 60;

        return WGS84toLambert2e(degree_long, minute_long, seconde_long, orientation_long, degree_lat, minute_lat, seconde_lat, orientation_lat, p_x, p_y);
    }

void WGS84ToLambert93(double latitude, double longitude, double & p_x, double & p_y)
{
    const double pi = acos(-1);
    /**** Conversion latitude,longitude en coordonée lambert 93 ****/
    // Projection conforme sécante, algo détailler dans NTG_71.pdf : http://www.ign.fr/affiche_rubrique.asp?rbr_id=1700&lng_id=FR
    //  > ACCUEIL > L'offre IGN Pro > Géodésie > RGF93 > Outils

    //variables:
    double a = 6378137; //demi grand axe de l'ellipsoide (m)
    double e = 0.08181919106; //première excentricité de l'ellipsoide
    double l0 = (pi / 180) * 3;
    double lc = l0;
    double phi0 = (pi / 180) * 46.5; //latitude d'origine en radian
    double phi1 = (pi / 180) * 44; //1er parallele automécoïque
    double phi2 = (pi / 180) * 49; //2eme parallele automécoïque

    double x0 = 700000; //coordonnées à l'origine
    double y0 = 6600000; //coordonnées à l'origine

    double phi = (pi / 180) * latitude;
    double l = (pi / 180) * longitude;

    //calcul des grandes normales
    double gN1 = a / sqrt(1 - e * e * sin(phi1) * sin(phi1));
    double gN2 = a / sqrt(1 - e * e * sin(phi2) * sin(phi2));

    //calculs des latitudes isométriques
    double gl1 = log(tan(pi / 4 + phi1 / 2) * pow((1 - e * sin(phi1)) / (1 + e * sin(phi1)), e / 2));
    double gl2 = log(tan(pi / 4 + phi2 / 2) * pow((1 - e * sin(phi2)) / (1 + e * sin(phi2)), e / 2));
    double gl0 = log(tan(pi / 4 + phi0 / 2) * pow((1 - e * sin(phi0)) / (1 + e * sin(phi0)), e / 2));
    double gl = log(tan(pi / 4 + phi / 2) * pow((1 - e * sin(phi)) / (1 + e * sin(phi)), e / 2));

    //calcul de l'exposant de la projection
    double n = (log((gN2 * cos(phi2)) / (gN1 * cos(phi1)))) / (gl1 - gl2);//ok

    //calcul de la constante de projection
    double c = ((gN1 * cos(phi1)) / n) * exp(n * gl1);//ok

    //calcul des coordonnées
    double ys = y0 + c * exp(-1 * n * gl0);

    double x93 = x0 + c * exp(-1 * n * gl) * sin(n * (l - lc));
    double y93 = ys - c * exp(-1 * n * gl) * cos(n * (l - lc));

    p_x = x93;
    p_y = y93;
}

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
          l_gui.refresh();
        }
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
      std::map<std::string,std::pair<double,double> > l_cities;
      double l_lat_min = std::numeric_limits<double>::max();
      double l_lat_max = std::numeric_limits<double>::min();
      double l_lon_min = std::numeric_limits<double>::max();
      double l_lon_max = std::numeric_limits<double>::min();
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
              if("" != l_name && "Ville" != l_name && "Ouessant" != l_name && "Île-de-Sein" != l_name && "Île-Molène" != l_name)
              {
                  l_cities.insert(std::map<std::string, std::pair<double, double> >::value_type(l_name,
                                                                                                std::pair<double,
                                                                                                          double
                                                                                                         >(l_lon,
                                                                                                           l_lat
                                                                                                          )
                                                                                               )
                                 );
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
          double l_lon = l_iter.second.first;
          double l_lat = l_iter.second.second;
          double l_lambert_lon;
          double l_lambert_lat;
#if 1
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
          double l_lon = l_iter.second.first;
          double l_lat = l_iter.second.second;
          double l_dx;
          double l_dy;
#if 0
#if 1
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
