#include <stdio.h>
#include "valuenoise.h"
#include "vectors.h"

#define width 128

int main () {
    struct double3 value;
    
    FILE* fp = fopen ("valuenoise.pgm", "w");
    fprintf (fp, "P2\n%d %d\n255\n", width, width);

    for (int y = 0; y < width; y++) {
        for (int x = 0; x < width; x++) {
            value = ValueNoise2D (x, y, 0.1, 5, 2, 0.5);
            fprintf (fp, "%u ", (unsigned)(value.x * 255));
        }
    }

    fclose (fp);
    return 0;
}