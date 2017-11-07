//
// Created by quickux on 07/11/17.
//
// Code has been adapted from:
// https://www.developpez.net/forums/d1562113/dotnet/langages/csharp/conversion-projections-wgs84-lambert-ii-etendu-lambert-93-a/

#ifndef GEOLOC_CELINE_PROJECTIONS_H
#define GEOLOC_CELINE_PROJECTIONS_H

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

#endif //GEOLOC_CELINE_PROJECTIONS_H
// EOF
