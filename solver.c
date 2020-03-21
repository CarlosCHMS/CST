#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "grid.h"
#include "boundary.h"
#include "input.h"
#include "solver.h"

void solverAloc(struct solverStruct* solver){

    solver->T = (FTYPE*)malloc(solver->grid->Np*sizeof(FTYPE));
    solver->Tnew = (FTYPE*)malloc(solver->grid->Np*sizeof(FTYPE));
    solver->flux = (FTYPE*)malloc(solver->grid->Np*sizeof(FTYPE));

}

void solverInit(struct solverStruct* solver, struct gridStruct* grid, struct inputStruct* input){

    //inputPrintParameters(input);

    solver->c = input->c;
    //printf("\n%f", input->c);
    solver->N = input->N;

    solver->grid = grid;
    solver->bound = (struct boundaryStruct*)malloc(sizeof(struct boundaryStruct));

    boundaryInit(solver->bound, solver->grid, input->markers, input->fInputs, input->types);

    //boundaryPrintTypes(solver->bound);

    solver->tol = 0.000001;
    solver->NmaxImplicit = 100;

    solver->saveStep = input->saveStep;
    solver->outputName = input->outName;

    solverAloc(solver);

    solverInitTemp(solver);

}

void solverInitTemp(struct solverStruct* solver){

    int ii, p0, p1, p2;

    for(ii=0; ii<solver->grid->Ne; ii++){

        if(boundaryElemIsDomain(solver->bound, ii)){

            gridGetElemPoints(solver->grid, ii, &p0, &p1, &p2);

            solver->T[p0] = solver->bound->fInputs[solver->bound->elemIndex[ii]];
            solver->T[p1] = solver->bound->fInputs[solver->bound->elemIndex[ii]];
            solver->T[p2] = solver->bound->fInputs[solver->bound->elemIndex[ii]];

        };

    };

    for(ii=0; ii<solver->grid->Ne; ii++){

        if(boundaryElemIsTemp(solver->bound, ii)){

            gridGetElemPoints(solver->grid, ii, &p0, &p1, &p2);

            solver->T[p0] = solver->bound->fInputs[solver->bound->elemIndex[ii]];
            solver->T[p1] = solver->bound->fInputs[solver->bound->elemIndex[ii]];

        };

    };

    for(ii=0; ii<solver->grid->Np; ii++){

        // Tnew must be initializaded beacuse of the points that are not propagate.
        solver->Tnew[ii] = solver->T[ii];

    };

}

void solverPrintT(struct solverStruct* solver){

    int ii;

    printf("\nTemperature:");

    for(ii=0; ii<solver->grid->Np; ii++){

        printf("\n%i: %f", ii, solver->T[ii]);

    };

}

void solverPrintParameters(struct solverStruct* solver){

    int ii;

    printf("\nParameters:");
    printf("\nc: %f", solver->c);
    printf("\nN: %i", solver->N);

    printf("\nBoundary markers: ");
    for(ii=0; ii<solver->grid->Nm; ii++){
        printf("%i, ", solver->bound->markers[ii]);
    };

    printf("\nBoundary inputs: ");
    for(ii=0; ii<solver->grid->Nm; ii++){
        printf("%f, ", solver->bound->fInputs[ii]);
    };

    printf("\nBoundary types: ");
    for(ii=0; ii<solver->grid->Nm; ii++){
        printf("%c, ", solver->bound->types[ii]);
    };


};

void solverPropagates(struct solverStruct* solver, FTYPE* T, FTYPE* flux){

    int ii, p0, p1, p2;
    FTYPE a0, a1, Tm, q, s;

    for(ii=0; ii<solver->grid->Np; ii++){

        flux[ii] = 0.0;

    };

    for(ii=0; ii<solver->grid->Ne; ii++){

        if(boundaryElemIsDomain(solver->bound, ii)){

            gridGetElemPoints(solver->grid, ii, &p0, &p1, &p2);

            Tm = (T[p0] + T[p1] + T[p2])/3;

            gridCalcGradCoef(solver->grid, ii, 0, &a0, &a1);
            flux[p0] += a0*(T[p1] - Tm) + a1*(T[p2] - Tm);

            gridCalcGradCoef(solver->grid, ii, 1, &a0, &a1);
            flux[p1] += a0*(T[p2] - Tm) + a1*(T[p0] - Tm);

            gridCalcGradCoef(solver->grid, ii, 2, &a0, &a1);
            flux[p2] += a0*(T[p0] - Tm) + a1*(T[p1] - Tm);

        }else if(boundaryElemIsHeat(solver->bound, ii)){

            //Constant heat flux boundary condition

            gridGetElemPoints(solver->grid, ii, &p0, &p1, &p2);
            s = gridCalcArea(solver->grid, ii);
            q = boundaryGetElemInput(solver->bound, ii);

            flux[p0] -= s*q/2;
            flux[p1] -= s*q/2;

        };

    };

}

void solverSimulateExplicity(struct solverStruct* solver){

    int ii, jj;
    FTYPE* aux;
    FILE* outFile;

    printf("\nEuler explicit solver");

    outFile = fopen(solver->outputName, "w");

    fprintf(outFile," %i, %i\n", solver->grid->Np, 0);
    for(jj=0; jj<solver->grid->Np; jj++){
        fprintf(outFile," %f,\n", solver->T[jj]);
    };

    for(ii=0; ii<solver->N; ii++){

        solverPropagates(solver, solver->T, solver->flux);

            //printf("\nFlux:");
        for(jj=0; jj<solver->grid->Np; jj++){

            if(solver->bound->pointPropag[jj]){

                solver->Tnew[jj] = solver->T[jj] - solver->c*solver->flux[jj]/solver->grid->dualArea[jj];
                //printf("\n%f: %f", solver->T[ii], solver->Tnew[ii]);

            };

        };

        aux = solver->T;
        solver->T = solver->Tnew;
        solver->Tnew = aux;

        if((ii+1)%solver->saveStep == 0){
            fprintf(outFile," %i, %i\n", solver->grid->Np, ii);
            for(jj=0; jj<solver->grid->Np; jj++){
                fprintf(outFile," %f,\n", solver->T[jj]);
            };
            printf("\nCalculating interaction %i", ii);
        };

    };

    fclose(outFile);

}

void solverSimulateImplicity(struct solverStruct* solver){

    int ii, jj, kk, flag;
    FTYPE* aux;
    FILE* outFile;
    FTYPE error;

    printf("\nEuler Implicit solver 1st order");

    solver->Tnew2 = (FTYPE*)malloc(solver->grid->Np*sizeof(FTYPE));
    //solver->fluxnew = (FTYPE*)malloc(solver->grid->Np*sizeof(FTYPE));

    outFile = fopen(solver->outputName, "w");

    for(ii=0; ii<solver->N; ii++){

        for(jj=0; jj<solver->grid->Np; jj++){

            solver->Tnew[jj] = solver->T[jj];
            solver->Tnew2[jj] = solver->T[jj];

        };

        //printf("\nFlux:");
        flag = 1;
        kk = 0;
        do{

            if(flag){
                flag = 0;
            }else{
                aux = solver->Tnew;
                solver->Tnew = solver->Tnew2;
                solver->Tnew2 = aux;
            };

            solverPropagates(solver, solver->Tnew, solver->flux);

            for(jj=0; jj<solver->grid->Np; jj++){

                if(solver->bound->pointPropag[jj]){

                    solver->Tnew2[jj] = solver->T[jj] - solver->c*solver->flux[jj]/solver->grid->dualArea[jj];
                    //printf("\n%f: %f", solver->Tnew[ii], solver->Tnew2[ii]);

                };

            };

            kk++;
            error = solverImplicitError(solver);
            printf("\nStep %i, Iteraction %i, Error: %f", ii, kk, error);

        }while( ( (error>solver->tol) & (kk<solver->NmaxImplicit) ) );

        aux = solver->T;
        solver->T = solver->Tnew;
        solver->Tnew = aux;

        if((ii+1)%solver->saveStep == 0){
            fprintf(outFile," %i, %i\n", solver->grid->Np, ii);
            for(jj=0; jj<solver->grid->Np; jj++){
                fprintf(outFile," %f,\n", solver->T[jj]);
            };
        };

    };

    fclose(outFile);

}

void solverSimulateImplicity2(struct solverStruct* solver){

    // Implicity second order in time
    printf("\nEuler Implicit solver 2nd order in time");

    int ii, jj, kk, flag;
    FTYPE* aux;
    FILE* outFile;
    FTYPE error;

    solver->Tnew2 = (FTYPE*)malloc(solver->grid->Np*sizeof(FTYPE));
    solver->fluxnew = (FTYPE*)malloc(solver->grid->Np*sizeof(FTYPE));

    outFile = fopen(solver->outputName, "w");

    for(ii=0; ii<solver->N; ii++){

        solverPropagates(solver, solver->T, solver->flux);

        for(jj=0; jj<solver->grid->Np; jj++){

            solver->Tnew[jj] = solver->T[jj];
            solver->Tnew2[jj] = solver->T[jj];

        };

        //printf("\nFlux:");
        flag = 1;
        kk = 0;
        do{

            if(flag){
                flag = 0;
            }else{
                aux = solver->Tnew;
                solver->Tnew = solver->Tnew2;
                solver->Tnew2 = aux;
            };

            solverPropagates(solver, solver->Tnew, solver->fluxnew);

            for(jj=0; jj<solver->grid->Np; jj++){

                if(solver->bound->pointPropag[jj]){

                    solver->Tnew2[jj] = solver->T[jj] - 0.5*solver->c*(solver->flux[jj] + solver->fluxnew[jj])/solver->grid->dualArea[jj];
                    //printf("\n%f: %f", solver->Tnew[ii], solver->Tnew2[ii]);

                };

            };

            kk++;
            error = solverImplicitError(solver);
            printf("\nStep %i, Iteraction %i, Error: %f", ii, kk, error);

        }while( ( (error>solver->tol) & (kk<solver->NmaxImplicit) ));

        aux = solver->T;
        solver->T = solver->Tnew;
        solver->Tnew = aux;

        if((ii+1)%solver->saveStep == 0){
            fprintf(outFile," %i, %i\n", solver->grid->Np, ii);
            for(jj=0; jj<solver->grid->Np; jj++){
                fprintf(outFile," %f,\n", solver->T[jj]);
            };
        };

    };

    fclose(outFile);

}

FTYPE solverImplicitError(struct solverStruct* solver){

    int jj;
    FTYPE error = 0.0;

    for(jj=0; jj<solver->grid->Np; jj++){

        error += (solver->Tnew2[jj] - solver->Tnew[jj])*(solver->Tnew2[jj] - solver->Tnew[jj]);

    };

    error /= solver->grid->Np;
    error = sqrt(error);

    return error;

};

void solverSalveVTK(struct solverStruct* solver){

    int Ne5, ii, p0, p1, p2;
    FILE* f1;

    f1 = fopen("output.vtk", "w");
    fprintf(f1, "# vtk DataFile Version 2.0\n");
    fprintf(f1, "CST output\n");
    fprintf(f1, "ASCII\n");
    fprintf(f1, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(f1, "POINTS %i float\n", solver->grid->Np);

    for(ii=0; ii<solver->grid->Np; ii++){
        fprintf(f1, "%f %f %f\n", solver->grid->x[ii], solver->grid->y[ii], 0.0);
    };
    fprintf(f1, "\n");

    Ne5 = 0;

    for(ii=0; ii<solver->grid->Ne; ii++){
        if(solver->grid->elem[ii][0] == 5){
            Ne5++;
        };
    };

    fprintf(f1, "CELLS %i %i\n", Ne5, Ne5*4);
    for(ii=0; ii<solver->grid->Ne; ii++){
        if(solver->grid->elem[ii][0] == 5){
            gridGetElemPoints(solver->grid, ii, &p0, &p1, &p2);
            fprintf(f1, "%i %i %i %i\n", 3, p0, p1, p2);
        };
    };
    fprintf(f1, "\n");

    fprintf(f1, "CELL_TYPES %i\n", Ne5);
    for(ii=0; ii<Ne5; ii++){
        fprintf(f1, "5\n");
    };
    fprintf(f1, "\n");

    fprintf(f1, "POINT_DATA %i\n", solver->grid->Np);

    //Scalar data
    fprintf(f1, "SCALARS temperature float 1\n");
    fprintf(f1, "LOOKUP_TABLE default\n");

    for(ii=0; ii<solver->grid->Np; ii++){
        fprintf(f1, "%f\n", solver->T[ii]);
    };

    fclose(f1);

}
