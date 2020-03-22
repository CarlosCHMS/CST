#ifndef SOLVER_H_INCLUDED
#define SOLVER_H_INCLUDED

struct solverStruct{

    //Number of iteractions
    int N;

    FTYPE c;
    FTYPE k;
    FTYPE* T;
    FTYPE* Tnew;
    FTYPE* Tnew2;
    FTYPE* flux;
    FTYPE* fluxnew;

    FTYPE tol;
    int NmaxImplicit;

    int Nsave;

    int* TPointShare;

    char* outputName;
    FILE outFile;

    struct gridStruct* grid;
    struct boundaryStruct* bound;

};

void solverAloc(struct solverStruct* solver);

void solverInit(struct solverStruct* solver, struct gridStruct* grid, struct inputStruct* input);

void solverInitTemp(struct solverStruct* solver);

void solverPrintT(struct solverStruct* solver);

void solverCalcTPointShare(struct solverStruct* solver);

void solverPrintParameters(struct solverStruct* solver);

void solverPropagates(struct solverStruct* solver, FTYPE* T, FTYPE* flux);

void solverSimulateExplicity(struct solverStruct* solver);

void solverSimulateImplicity(struct solverStruct* solver);

void solverSimulateImplicity2(struct solverStruct* solver);

FTYPE solverImplicitError(struct solverStruct* solver);

void solverSalveVTK(struct solverStruct* solver);

#endif // SOLVER_H_INCLUDED
