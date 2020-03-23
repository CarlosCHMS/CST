#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "grid.h"
#include "input.h"
#include "boundary.h"


void boundaryAlloc(struct boundaryStruct* bound){

    bound->elemIndex = (int*)malloc(bound->grid->Ne*sizeof(int));
    bound->pointPropag = (int*)malloc(bound->grid->Np*sizeof(int));

}

void boundaryInit(struct boundaryStruct* bound, struct gridStruct* grid, struct inputStruct* input){

    bound->markers = input->markers;

    bound->C = input->C;
    bound->K = input->K;

    bound->fInputs = input->fInputs;
    bound->fInputs1 = input->fInputs1;
    bound->fInputs4 = input->fInputs4;
    bound->types = input->types;
    bound->grid = grid;

    boundaryAlloc(bound);
    boundaryInitIndex(bound);
    boundarySetPointPropag(bound);

}

void boundaryInitIndex(struct boundaryStruct* bound){

    int ii, jj;

    for(ii=0; ii<bound->grid->Ne; ii++){

        for(jj=0; jj<bound->grid->Nm; jj++){

            if(bound->grid->elem[ii][1] == bound->markers[jj]){

                bound->elemIndex[ii] = jj;

            };

        };

    };

}

void boundarySetPointPropag(struct boundaryStruct* bound){

    int ii, p0, p1, p2;

    for(ii=0; ii<bound->grid->Np; ii++){

        bound->pointPropag[ii] = 1;

    };

    for(ii=0; ii<bound->grid->Ne; ii++){

        if(boundaryElemIsTemp(bound, ii)){

            gridGetElemPoints(bound->grid, ii, &p0, &p1, &p2);

            bound->pointPropag[p0] = 0;
            bound->pointPropag[p1] = 0;

        };

    };

}

FTYPE boundaryGetElemInput(struct boundaryStruct* bound, int ii){

    return bound->fInputs[bound->elemIndex[ii]];

}

void boundaryPrintTypes(struct boundaryStruct* bound){

    int ii, jj;

    printf("\nElement types:");
    for(ii=0; ii<bound->grid->Ne; ii++){
        jj = bound->elemIndex[ii];
        printf("\n%i, %i, %f", ii, bound->types[jj], bound->fInputs[jj]);
        //printf("\n%i,", bound->elemTypes[ii]);

    };

}

int boundaryElemIsDomain(struct boundaryStruct* bound, int ie){

    int flag;

    if(bound->types[bound->elemIndex[ie]] == 'D'){

        flag = 1;

    }else{

        flag = 0;

    };

    return flag;

}

int boundaryElemIsTemp(struct boundaryStruct* bound, int ie){

    int flag;

    if(bound->types[bound->elemIndex[ie]] == 'T'){

        flag = 1;

    }else{

        flag = 0;

    };

    return flag;

}

int boundaryElemIsHeat(struct boundaryStruct* bound, int ie){

    int flag;

    if(bound->types[bound->elemIndex[ie]] == 'Q'){

        flag = 1;

    }else{

        flag = 0;

    };

    return flag;

}

