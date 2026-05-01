#pragma once
#ifndef FREE_ARRAY_H
#define FREE_ARRAY_H

extern double*** Ex, *** Ey, *** Ez, *** Hx, *** Hy, *** Hz;
extern double*** CBx, *** CBy, *** CBz;
extern double*** CQx, *** CQy, *** CQz;
extern int MAX_X, MAX_Y, MAX_Z;

void free_array() {

    free3DArray(Ex, MAX_X, MAX_Y + 1);
    free3DArray(Ey, MAX_X + 1, MAX_Y);
    free3DArray(Ez, MAX_X + 1, MAX_Y + 1);  // ÐÞÕý£ºEz µÄÎ¬¶È
    free3DArray(Hx, MAX_X + 1, MAX_Y);
    free3DArray(Hy, MAX_X, MAX_Y + 1);
    free3DArray(Hz, MAX_X, MAX_Y);


    free3DArray(CBx, MAX_X, MAX_Y + 1);
    free3DArray(CBy, MAX_X + 1, MAX_Y);
    free3DArray(CBz, MAX_X + 1, MAX_Y + 1);
    free3DArray(CQx, MAX_X + 1, MAX_Y);
    free3DArray(CQy, MAX_X, MAX_Y + 1);
    free3DArray(CQz, MAX_X, MAX_Y);

    free1DArray(source);

}
#endif

