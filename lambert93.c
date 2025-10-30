#include <math.h>
#include <float.h>

#include "lambert93.h"

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif


typedef struct {
    // Ellipsoid
    double a;      // semi-major
    double f;      // flattening
    double e;      // eccentricity
    // Projection params (radians unless noted)
    double phi1;   // 1st standard parallel
    double phi2;   // 2nd standard parallel
    double phi0;   // latitude of origin
    double lam0;   // central meridian
    double X0;     // false easting (m)
    double Y0;     // false northing (m)
    // Derived
    double n;
    double F;
    double rho0;
} lcc2sp_t;

/*
// Initialize generic LCC 2SP (angles in degrees).
void lcc2sp_init(lcc2sp_t* p,
                 double a, double f,
                 double phi1_deg, double phi2_deg,
                 double phi0_deg, double lam0_deg,
                 double X0, double Y0);

// Forward/inverse using a configured projection.
void lcc2sp_forward(const lcc2sp_t* p,
                    double lon_deg, double lat_deg,
                    double* out_x, double* out_y);

int lcc2sp_inverse(const lcc2sp_t* p,
                    double x, double y,
                    double* out_lon_deg, double* out_lat_deg);

// Ready-to-use Lambert-93 wrappers (RGF93/GRS80).
void lambert93_forward(double lon_deg, double lat_deg, double* out_x, double* out_y);
int lambert93_inverse(double x, double y, double* out_lon_deg, double* out_lat_deg);
*/


static inline double deg2rad(double d) { return d * (M_PI / 180.0); }
static inline double rad2deg(double r) { return r * (180.0 / M_PI); }

static inline double m_func(double phi, double e2) {
    double s = sin(phi);
    return cos(phi) / sqrt(1.0 - e2 * s * s);
}
static inline double t_func(double phi, double e) {
    double es = e * sin(phi);
    double num = tan(M_PI/4.0 - phi/2.0);
    double pow_e = pow((1.0 - es) / (1.0 + es), e/2.0);
    return num / pow_e;
}

static void lcc2sp_init(lcc2sp_t* p,
                 double a, double f,
                 double phi1_deg, double phi2_deg,
                 double phi0_deg, double lam0_deg,
                 double X0, double Y0)
{
    p->a = a;
    p->f = f;
    p->e = sqrt(f * (2.0 - f));
    p->phi1 = deg2rad(phi1_deg);
    p->phi2 = deg2rad(phi2_deg);
    p->phi0 = deg2rad(phi0_deg);
    p->lam0 = deg2rad(lam0_deg);
    p->X0 = X0;
    p->Y0 = Y0;

    double e2 = p->e * p->e;

    double m1 = m_func(p->phi1, e2);
    double m2 = m_func(p->phi2, e2);
    double t1 = t_func(p->phi1, p->e);
    double t2 = t_func(p->phi2, p->e);
    double t0 = t_func(p->phi0, p->e);

    p->n = (log(m1) - log(m2)) / (log(t1) - log(t2));
    p->F = m1 / (p->n * pow(t1, p->n));
    p->rho0 = p->a * p->F * pow(t0, p->n);
}


static lcc2sp_t L93;
static int L93_inited = 0;

__attribute__ ((constructor)) static void local_ctr(void)
{
    lcc2sp_init(&L93,
        6378137.0, 1.0/298.257222101,
        49.0, 44.0, 46.5, 3.0,
        700000.0, 6600000.0
    );
    L93_inited = 1;
}



void lcc2sp_forward(const lcc2sp_t* p,
                    double lat_deg, double lon_deg,
                    double* out_x, double* out_y)
{
    // Clamp latitude to avoid tan(π/4 - φ/2) singularities
    double lat = deg2rad(fmax(fmin(lat_deg, 89.999999), -89.999999));
    double lon = deg2rad(lon_deg);

    double t  = t_func(lat, p->e);
    double rho = p->a * p->F * pow(t, p->n);
    double theta = p->n * (lon - p->lam0);

    double x = p->X0 + rho * sin(theta);
    double y = p->Y0 + p->rho0 - rho * cos(theta);

    *out_x = x;
    *out_y = y;
}

int lcc2sp_inverse(const lcc2sp_t* p,
                    double x, double y,
                    double* out_lat_deg, double* out_lon_deg)
{
    double dx = x - p->X0;
    double dy = p->rho0 - (y - p->Y0);
    double rho = hypot(dx, dy);
    if (!(rho > 0.0)) return -1;

    double theta = atan2(dx, dy);
    double t = pow(rho / (p->a * p->F), 1.0 / p->n);

    // Iterative solve for φ from t (EPSG:9802 recommendation).
    // Start with spherical approximation:
    double phi = M_PI/2.0 - 2.0 * atan(t);
    for (int i = 0; i < 10; ++i) {
        double es = p->e * sin(phi);
        double phi_next = M_PI/2.0 - 2.0 * atan( t * pow((1.0 - es)/(1.0 + es), p->e/2.0) );
        if (fabs(phi_next - phi) < 1e-12) { phi = phi_next; break; }
        phi = phi_next;
    }
    double lam = p->lam0 + theta / p->n;

    *out_lon_deg = rad2deg(lam);
    *out_lat_deg = rad2deg(phi);
    return 0;
}

// -------- Lambert-93 ready-to-use (RGF93 / GRS80) --------

void lambert93_latlon2xy(double lat_deg, double lon_deg, double* out_x, double* out_y) {
    lcc2sp_forward(&L93, lat_deg, lon_deg, out_x, out_y);
}

int lambert93_xy2latlon(double x, double y, double* out_lat_deg, double* out_lon_deg) {
    return lcc2sp_inverse(&L93, x, y, out_lat_deg, out_lon_deg);
}
