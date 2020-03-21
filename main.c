#include<stdio.h>
#include<stdlib.h>
#include "grid.h"
#include "solver.h"
#include "input.h"

int main(){

    printf("CSU Core - Code of Simulation for Unstructured meshes Core");

    struct inputStruct* input = malloc(sizeof(struct inputStruct));
    struct gridStruct* grid = malloc(sizeof(struct gridStruct));
    struct solverStruct* solver = malloc(sizeof(struct solverStruct));

    inputInit(input);

    inputPrintParameters(input);

    gridInit(grid, input->meshName, input->Lref);

    gridPrintNumbers(grid);

    gridPrintMarkers(grid);

    //gridPrintCoords(grid);

    //gridPrintElem(grid);

    //gridCheckOrient(grid);

    //gridPrintArea(grid);

    //gridPrintDualArea(grid);

    //gridPrintGradCoef(grid);

    //gridCheckGradCoef(grid);

    //gridPrintPointMarkers(grid);

    //printf("\n%i", gridElementForPoint(grid, -0.5, -0.2));

    solverInit(solver, grid, input);

    //solverPrintParameters(solver);

    //solver->saveStep = input->N;

    //solverPrintT(solver);

    //solverSimulateImplicity(solver);

    //solverSimulateImplicity2(solver);

    solverSimulateExplicity(solver);

    //solverSalveVTK(solver);

    //solverPrintT(solver);

    printf("\n");

    return 0;
}
