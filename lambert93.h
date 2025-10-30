// MIT License
// Lambert-93 (EPSG:2154) <-> WGS84, plus generic LCC 2SP

void lambert93_latlon2xy(double lat_deg, double lon_deg, double* out_x, double* out_y);
int lambert93_xy2latlon(double x, double y, double* out_lat_deg, double* out_lon_deg);
