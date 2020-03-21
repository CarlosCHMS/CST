#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

typedef double FTYPE;

struct gridStruct{
    int Nm;
    int Np;
    int Ne;
    int maxElemEntry;
    int type1, type2;

    int* marker;

    int* pointMarker;

    FTYPE Lref;
    FTYPE* x;
    FTYPE* y;
    FTYPE* dualArea;

    int** elem;
};

void gridInit(struct gridStruct* grid, char* fileName, FTYPE Lref);

void gridPrintNumbers(struct gridStruct* grid);

void gridPrintMarkers(struct gridStruct* grid);

void gridPrintCoords(struct gridStruct* grid);

void gridPrintElem(struct gridStruct* grid);

void gridPrintArea(struct gridStruct* grid);

float gridCalcArea(struct gridStruct* grid, int ii);

void gridCalcDualArea(struct gridStruct* grid);

void gridCalcGradCoef(struct gridStruct* grid, int ie, int ip, FTYPE* a0, FTYPE* a1);

void gridPrintGradCoef(struct gridStruct* grid);

void gridCheckGradCoef(struct gridStruct* grid);

void gridPrintDualArea(struct gridStruct* grid);

void gridCheckOrient(struct gridStruct* grid);

void gridDelete(struct gridStruct* grid);

void gridSetPointMarkers(struct gridStruct* grid);

void gridPrintPointMarkers(struct gridStruct* grid);

void gridGetElemPoints(struct gridStruct* grid, int ie, int* p0, int* p1, int* p2);

int gridElementForPoint(struct gridStruct* grid, FTYPE x, FTYPE y);

void gridTest();

#endif // GRID_H_INCLUDED
