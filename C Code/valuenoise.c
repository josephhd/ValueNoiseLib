//Created by Joseph Dunbar -- June 6, 2018

#include "valuenoise.h"
#include "vectors.h"

#define EASE(x) ((x)*(x)*(x)*((x)*((x)*6.0-15.0)+10))   //6x^5-15x^4+10x^3
#define DERIV(x) (30*(x)*(x)*((x)*((x)-2)+1))           //30x^4-60x^3+30x^2

static const unsigned char hash[] = {
        151,160,137, 91, 90, 15,131, 13,201, 95, 96, 53,194,233,  7,225,
		140, 36,103, 30, 69,142,  8, 99, 37,240, 21, 10, 23,190,  6,148,
		247,120,234, 75,  0, 26,197, 62, 94,252,219,203,117, 35, 11, 32,
		 57,177, 33, 88,237,149, 56, 87,174, 20,125,136,171,168, 68,175,
		 74,165, 71,134,139, 48, 27,166, 77,146,158,231, 83,111,229,122,
		 60,211,133,230,220,105, 92, 41, 55, 46,245, 40,244,102,143, 54,
		 65, 25, 63,161,  1,216, 80, 73,209, 76,132,187,208, 89, 18,169,
		200,196,135,130,116,188,159, 86,164,100,109,198,173,186,  3, 64,
		 52,217,226,250,124,123,  5,202, 38,147,118,126,255, 82, 85,212,
		207,206, 59,227, 47, 16, 58, 17,182,189, 28, 42,223,183,170,213,
		119,248,152,  2, 44,154,163, 70,221,153,101,155,167, 43,172,  9,
		129, 22, 39,253, 19, 98,108,110, 79,113,224,232,178,185,112,104,
		218,246, 97,228,251, 34,242,193,238,210,144, 12,191,179,162,241,
		 81, 51,145,235,249, 14,239,107, 49,192,214, 31,181,199,106,157,
		184, 84,204,176,115,121, 50, 45,127,  4,150,254,138,236,205, 93,
		222,114, 67, 29, 24, 72,243,141,128,195, 78, 66,215, 61,156,180,
        151,160,137, 91, 90, 15,131, 13,201, 95, 96, 53,194,233,  7,225,
		140, 36,103, 30, 69,142,  8, 99, 37,240, 21, 10, 23,190,  6,148,
		247,120,234, 75,  0, 26,197, 62, 94,252,219,203,117, 35, 11, 32,
		 57,177, 33, 88,237,149, 56, 87,174, 20,125,136,171,168, 68,175,
		 74,165, 71,134,139, 48, 27,166, 77,146,158,231, 83,111,229,122,
		 60,211,133,230,220,105, 92, 41, 55, 46,245, 40,244,102,143, 54,
		 65, 25, 63,161,  1,216, 80, 73,209, 76,132,187,208, 89, 18,169,
		200,196,135,130,116,188,159, 86,164,100,109,198,173,186,  3, 64,
		 52,217,226,250,124,123,  5,202, 38,147,118,126,255, 82, 85,212,
		207,206, 59,227, 47, 16, 58, 17,182,189, 28, 42,223,183,170,213,
		119,248,152,  2, 44,154,163, 70,221,153,101,155,167, 43,172,  9,
		129, 22, 39,253, 19, 98,108,110, 79,113,224,232,178,185,112,104,
		218,246, 97,228,251, 34,242,193,238,210,144, 12,191,179,162,241,
		 81, 51,145,235,249, 14,239,107, 49,192,214, 31,181,199,106,157,
		184, 84,204,176,115,121, 50, 45,127,  4,150,254,138,236,205, 93,
		222,114, 67, 29, 24, 72,243,141,128,195, 78, 66,215, 61,156,180
    };




static inline double HashValue1D (int x) {
    x &= 255;

    return (double)hash[x] / 255;
}

static inline double HashValue2D (int x, int y) {
    x &= 255;
    y &= 255;

    return (double)hash[hash[x] + y]  / 255;
}

static inline double HashValue3D (int x, int y, int z) {
    x &= 255;
    y &= 255;
    z &= 255;

    return (double)hash[hash[hash[x] + y] + z] / 255;
}

static inline double HashValue4D (int x, int y, int z, int w) {
    x &= 255;
    y &= 255;
    z &= 255;
    w &= 255;

    return (double)hash[hash[hash[hash[x] + y] + z] + w] / 255;
}

struct double2 ValueNoise1D (double x, double frequency, int octaves, double lacunarity, double persistence) {
    if (x < 0) x *= -1;

    x *= frequency;

    double range = 0;
    double l = 1;
    double p = 1;
    struct double2 result = {0, 0};

    for (int i = 0; i < octaves; i++) {     
        //find lattice points
        int x0 = (int)(x);
        int x1 = x0 + 1;

        //find distance from lattice point
        double xd = EASE (x - x0);

        //get hash values at lattice points
        double h0 = HashValue1D (x0);
        double h1 = HashValue1D (x1);

        //interpolate in x-dimension
        double h = h0 * (1 - xd) + h1 * xd;

        //store derivative for xd
        double dxd = DERIV (x - x0);

        //calculate the gradient, multiply by the lacunarity because of chainrule
        double dhx = l * dxd * (h1 - h0);

        //add results
        result.x += p * h;
        result.y += p * dhx;
        
        range += p;
        l *= lacunarity;
        p *= persistence;
        x *= lacunarity;
        range += p;
    }

    result.x /= range;
    result.y /= range;

    return result;
}

struct double3 ValueNoise2D (double x, double y, double frequency, int octaves, double lacunarity, double persistence) {
    if (x < 0) x *= -1;
    if (y < 0) y *= -1;
    
    x *= frequency;
    y *= frequency;

    double range = 0;
    double l = 1;
    double p = 1;
    struct double3 result = {0,0,0};

    for (int i = 0; i < octaves; i++) {
        //find lattice points
        int x0 = (int)(x);
        int x1 = x0 + 1;
        int y0 = (int)(y);
        int y1 = y0 + 1;

        //find distance from lattice point
        double xd = EASE(x - x0);
        double yd = EASE(y - y0);

        //get hash values at lattice points (name format is hxy)
        double h00 = HashValue2D(x0, y0);
        double h10 = HashValue2D(x1, y0);
        double h01 = HashValue2D(x0, y1);
        double h11 = HashValue2D(x1, y1);

        //interpolate in x-dimension
        double h0 = h00 * (1 - xd) + h10 * xd;
        double h1 = h01 * (1 - xd) + h11 * xd;

        //interpolate in y-dimension
        double h = h0 * (1 - yd) + h1 * yd;

        //store the derivative for xd, and yd
        double dxd = DERIV(x - x0);
        double dyd = DERIV(y - y0);

        //calculate the gradient, multiply by the lacunarity because of chainrule
        double dhx = l * dxd * ((1 - yd) * (h10 - h00) +
                                     yd * (h11 - h01));

        double dhy = l * dyd * (h1 - h0);

        //add results
        result.x += p * h;
        result.y += p * dhx;
        result.z += p * dhy;

        range += p;
        l *= lacunarity;
        p *= persistence;
        x *= lacunarity;
        y *= lacunarity;
    }

    result.x /= range;
    result.y /= range;
    result.z /= range;

    return result;
}

struct double4 ValueNoise3D (double x, double y, double z, double frequency, int octaves, double lacunarity, double persistence) {
    if (x < 0) x *= -1;
    if (y < 0) y *= -1;
    if (z < 0) z *= -1;

    x *= frequency;
    y *= frequency;
    z *= frequency;

    double range = 0;
    double l = 1;
    double p = 1;
    struct double4 result = {0, 0, 0, 0};

    for (int i = 0; i < octaves; i++) {
        //get lattice points
        int x0 = (int)(x);
        int x1 = x0+1;
        int y0 = (int)(y);
        int y1 = y0 + 1;
        int z0 = (int)(z);
        int z1 = z0 + 1;

        //find distance to lattice point
        double xd = EASE (x - x0);
        double yd = EASE (y - y0);
        double zd = EASE (z - z0);

        //get hash values at lattice points
        double h000 = HashValue3D (x0, y0, z0);
        double h100 = HashValue3D (x1, y0, z0);
        double h010 = HashValue3D (x0, y1, z0);
        double h110 = HashValue3D (x1, y1, z0);
        double h001 = HashValue3D (x0, y0, z1);
        double h101 = HashValue3D (x1, y0, z1);
        double h011 = HashValue3D (x0, y1, z1);
        double h111 = HashValue3D (x1, y1, z1);

        //interpolate in x-dimension
        double h00 = h000 * (1 - xd) + xd * h100;
        double h10 = h010 * (1 - xd) + xd * h110;
        double h01 = h001 * (1 - xd) + xd * h101;
        double h11 = h011 * (1 - xd) + xd * h111;

        //interpolate in y-dimension
        double h0 = h00 * (1 - yd) + yd * h10;
        double h1 = h01 * (1 - yd) + yd * h11;

        //interpolate in z-dimension
        double h = h0 * (1 - zd) + zd * h1;
        
        //derivatives of distance values
        double dxd = DERIV (x - x0);
        double dyd = DERIV (y - y0);
        double dzd = DERIV (z - z0);

        //calculate the gradient, multiply by the lacunarity because of chainrule
        double dhx = l * dxd * ((1 - zd) * ((1 - yd) * (h100 - h000) + yd * (h110 - h010)) + 
                                      zd * ((1 - yd) * (h101 - h001) + yd * (h111 - h011)));
        
        double dhy = l * dyd * ((1 - zd) * (h10 - h00) + 
                                      zd * (h11 - h01));
        
        double dhz = l * dzd * (h1 - h0);

        //add final results        
        result.x += p * h;
        result.y += p * dhx;
        result.z += p * dhy;
        result.w += p * dhz;

        range += p;
        l *= lacunarity;
        p *= persistence;
        x *= lacunarity;
        y *= lacunarity;
        z *= lacunarity;
    }

    result.x /= range;
    result.y /= range;
    result.z /= range;
    result.w /= range;

    return result;
}

struct double5 ValueNoise4D (double x, double y, double z, double w, double frequency, int octaves, double lacunarity, double persistence) {
    if (x < 0) x *= -1;
    if (y < 0) y *= -1;
    if (z < 0) z *= -1;
    if (w < 0) w *= -1;

    x *= frequency;
    y *= frequency;
    z *= frequency;
    w *= frequency;

    double range = 0;
    double l = 1;
    double p = 1;
    struct double5 result = {0, 0, 0, 0, 0};

    for (int i = 0; i < octaves; i++) {
        //get lattice points
        int x0 = (int)(x);
        int x1 = x0 + 1;
        int y0 = (int)(y);
        int y1 = y0 + 1;
        int z0 = (int)(z);
        int z1 = z0 + 1;
        int w0 = (int)(w);
        int w1 = w0 + 1;

        //find distance to lattice point
        double xd = EASE (x - x0);
        double yd = EASE (y - y0);
        double zd = EASE (z - z0);
        double wd = EASE (w - w0);

        //get hash values at lattice points
        double h0000 = HashValue4D (x0, y0, z0, w0);
        double h1000 = HashValue4D (x1, y0, z0, w0);
        double h0100 = HashValue4D (x0, y1, z0, w0);
        double h1100 = HashValue4D (x1, y1, z0, w0);
        double h0010 = HashValue4D (x0, y0, z1, w0);
        double h1010 = HashValue4D (x1, y0, z1, w0);
        double h0110 = HashValue4D (x0, y1, z1, w0);
        double h1110 = HashValue4D (x1, y1, z1, w0);
        double h0001 = HashValue4D (x0, y0, z0, w1);
        double h1001 = HashValue4D (x1, y0, z0, w1);
        double h0101 = HashValue4D (x0, y1, z0, w1);
        double h1101 = HashValue4D (x1, y1, z0, w1);
        double h0011 = HashValue4D (x0, y0, z1, w1);
        double h1011 = HashValue4D (x1, y0, z1, w1);
        double h0111 = HashValue4D (x0, y1, z1, w1);
        double h1111 = HashValue4D (x1, y1, z1, w1);

        //interpolate in x-dimension
        double h000 = h0000 * (1 - xd) + h1000 * xd;
        double h100 = h0100 * (1 - xd) + h1100 * xd;
        double h010 = h0010 * (1 - xd) + h1010 * xd;
        double h110 = h0110 * (1 - xd) + h1110 * xd;
        double h001 = h0001 * (1 - xd) + h1001 * xd;
        double h101 = h0101 * (1 - xd) + h1101 * xd;
        double h011 = h0011 * (1 - xd) + h1011 * xd;
        double h111 = h0111 * (1 - xd) + h1111 * xd;

        //interpolate in y-dimension
        double h00 = h000 * (1 - yd) + h100 * yd;
        double h10 = h010 * (1 - yd) + h110 * yd;
        double h01 = h001 * (1 - yd) + h101 * yd;
        double h11 = h011 * (1 - yd) + h111 * yd;

        //interpolate in z-dimension
        double h0 = h00 * (1 - zd) + h10 * zd;
        double h1 = h01 * (1 - zd) + h11 * zd;

        //interpolate in w-dimension
        double h = h0 * (1 - wd) + h1 * wd;

        //derivatives of distance values
        double dxd = DERIV (x - x0);
        double dyd = DERIV (y - y0);
        double dzd = DERIV (z - z0);
        double dwd = DERIV (w - w0);

        //calculate the gradient, multiply by the lacunarity because of chainrule
        //if you're trying to make sense of this part, stop now
        double dhx = l * dxd * (1 - wd) * ((1 - zd) * ((1 - yd) * (h1000 - h0000) + yd * (h1100 - h0100)) + zd * ((1 - yd) * (h1010 - h0010) + yd * (h1110 - h0110))) + 
                                     wd * ((1 - zd) * ((1 - yd) * (h1001 - h0001) + yd * (h1101 - h0101)) + zd * ((1 - yd) * (h1011 - h0011) + yd * (h1111 - h0111)));

        double dhy = l * (1 - wd) * ((1 - zd) * (h100 * dyd - h000 * dyd) + zd * (h110 * dyd - h010 * dyd))
                         + wd * ((1 - zd) * (h101 * dyd - h001 * dyd) + zd * (h111 * dyd - h011 * dyd));
        
        double dhz = l * (1 - wd) * (h10 * dzd - h00 * dzd)
                         + wd * (h11 * dzd - h01 * dzd);
        
        double dhw = l * (h1 * dwd - h0 * dwd);

        //add final results        
        result.x += p * h;
        result.y += p * dhx;
        result.z += p * dhy;
        result.w += p * dhz;
        result.t += p * dhw;

        range += p;
        l *= lacunarity;
        p *= persistence;
        x *= lacunarity;
        y *= lacunarity;
        z *= lacunarity;
        w *= lacunarity;
    }

    result.x /= range;
    result.y /= range;
    result.z /= range;
    result.w /= range;
    result.t /= range;

    return result;
}
