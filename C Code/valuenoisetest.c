#include <stdio.h>
#include "valuenoise.h"
#include "vectors.h"

#define width 256
#define height 256

int main () {
    struct double3 value;
    
    FILE* fp = fopen ("valuenoise.pgm", "w");
    fprintf (fp, "P2\n%d %d\n255\n", width, height);

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            value = ValueNoise2D (x, y, 0.05, 5, 2, 0.5);
            fprintf (fp, "%u ", (unsigned)(value.x * 255));
        }

        fprintf (fp, "\n");
    }

    fclose (fp);
    return 0;
}