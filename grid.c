#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "input.h"
#include "grid.h"

void gridAloc(struct gridStruct* grid){

    int ii;

    grid->marker = (int*)malloc(grid->Nm*sizeof(int));

    grid->x = (FTYPE*)malloc(grid->Np*sizeof(FTYPE));
    grid->y = (FTYPE*)malloc(grid->Np*sizeof(FTYPE));
    grid->dualArea = (FTYPE*)malloc(grid->Np*sizeof(FTYPE));

    grid->elem = (int**)malloc(grid->Ne*sizeof(int*));
    grid->pointMarker = (int*)malloc(grid->Np*sizeof(int));

    for(ii=0; ii<grid->Ne; ii++){
        grid->elem[ii] = (int*)malloc(grid->maxElemEntry*sizeof(int));
    };

}

void gridInit(struct gridStruct* grid, struct inputStruct* input){

    FILE* f1;
    int ii, jj, kk, alocFlag;
    char c;
    char s[100];

    grid->type1 = 4;
    grid->type2 = 5;
    grid->maxElemEntry = 5;
    f1 = fopen(input->meshName, "r");
    grid->geo = input->geo;

    ii = 0;
    jj = 0;
    kk = 0;
    alocFlag = 1;
    //Read mesh file
    while(c != EOF){
        c = getc(f1);

        if(c == ','){
            //When the string gets a coma:
            s[jj] = '\0';
            jj = 0;
            //printf("%s, ", s);
            //printf("\ni: %i, k: %i, s: %s", ii, kk, s);
            if(ii==0){
                //This if gets the number of points and elments
                if(kk==0){
                    grid->Nm = strtod(s, NULL);
                    //N1 = grid->Nm + 1;
                }else if(kk==1){
                    grid->Np = strtod(s, NULL);
                    //N2 = grid->Np + N1;
                }else{
                    grid->Ne = strtod(s, NULL);

                };

            }else if(ii<grid->Nm + 1){

                if(alocFlag){
                    //This if works for allocation
                    gridAloc(grid);
                    alocFlag=0;
                };

                if(kk==0){
                    grid->marker[ii-1] = strtod(s, NULL);
                };

            }else if(ii<grid->Nm + grid->Np + 1){

                //This if get coordinates properties
                if(kk==0){
                    grid->x[ii-(grid->Nm + 1)] = strtod(s, NULL);
                }else{
                    grid->y[ii-(grid->Nm + 1)] = strtod(s, NULL);
                };

            }else if(ii<grid->Nm + grid->Np + grid->Ne + 1){
                //This if get element properties
                //printf("\nelem: %i", strtol(s, NULL, 10));//, NULL));
                grid->elem[ii-(grid->Nm + grid->Np + 1)][kk] = strtod(s, NULL);
                //printf("\nelem: %i", grid->elem[ii-(grid->Np + 1)]);

            };

            kk ++;
        }else if(c == '\n'){
            ii++;
            kk = 0;
        }else{
            s[jj] = c;
            jj++;
        };
    };

    gridCalcDualArea(grid);

    gridSetPointMarkers(grid);

}

void gridDelete(struct gridStruct* grid){

    int ii;

    free(grid->x);
    free(grid->y);

    for(ii=0; ii<grid->Ne; ii++){
        free(grid->elem[ii]);
    };

    free(grid->elem);
    free(grid);

}

void gridPrintNumbers(struct gridStruct* grid){

    printf("\nNumbers:");
    printf("\nNm:%i, Np:%i, Ne:%i", grid->Nm, grid->Np, grid->Ne);

}

void gridPrintMarkers(struct gridStruct* grid){

    int ii;

    printf("\nMarkers:");

    for(ii=0;ii<grid->Nm;ii++){
        printf("\n%i,", grid->marker[ii]);
    };

}

void gridPrintCoords(struct gridStruct* grid){

    int ii;

    printf("\nCoordinates:");
    printf("\ni, x, y,");

    for(ii=0;ii<grid->Np;ii++){
        printf("\n%i, %f, %f,", ii, grid->x[ii], grid->y[ii]);
    };

}

void gridPrintElem(struct gridStruct* grid){

    int ii, jj;

    printf("\nElements:");
    //printf("\ni, p1, p2, p3,", grid->Np, grid->Ne);

    //printf("%i, ", grid->elem[0][0]);

    for(ii=0;ii<grid->Ne;ii++){
        printf("\n");
        //printf("%i, ", grid->elem[ii][0]);
        for(jj=0;jj<grid->elem[ii][0];jj++){
            printf("%i, ", grid->elem[ii][jj]);
        };

    };

}

float gridCalcArea(struct gridStruct* grid, int ii){

    /*
    For 2d elements it returns its area
    For 1d elements it returns its length
    */

    int p1, p2, p3;
    FTYPE dx21, dx31, dy21, dy31, A;

    gridGetElemPoints(grid, ii, &p1, &p2, &p3);

    if(grid->elem[ii][0]==grid->type2){

        dx21 = grid->x[p2] - grid->x[p1];
        dx31 = grid->x[p3] - grid->x[p1];

        dy21 = grid->y[p2] - grid->y[p1];
        dy31 = grid->y[p3] - grid->y[p1];

        A = (dx21*dy31 - dy21*dx31)/2;

    }else{

        dx21 = grid->x[p2] - grid->x[p1];
        dy21 = grid->y[p2] - grid->y[p1];

        A = sqrt(dx21*dx21 + dy21*dy21);

    };

    return A;

}

void gridCalcDualArea(struct gridStruct* grid){

    int p1, p2, p3, ii;
    FTYPE A, y1, y2, y3;

    for(ii=0; ii<grid->Np; ii++){
        grid->dualArea[ii] = 0.0;
    };

    for(ii=0; ii<grid->Ne; ii++){
        if(grid->elem[ii][0]==grid->type2){

            gridGetElemPoints(grid, ii, &p1, &p2, &p3);

            A = gridCalcArea(grid, ii)/3;

            if(grid->geo == 1){

                y1 = grid->y[p1];
                y2 = grid->y[p2];
                y3 = grid->y[p3];

                grid->dualArea[p1] += A*(22*y1 + 7*y2 + 7*y3)/36;
                grid->dualArea[p2] += A*(22*y2 + 7*y3 + 7*y1)/36;
                grid->dualArea[p3] += A*(22*y3 + 7*y1 + 7*y2)/36;

            }else{

                grid->dualArea[p1] += A;
                grid->dualArea[p2] += A;
                grid->dualArea[p3] += A;

            };


        };

    };

}

void gridCalcGradCoef(struct gridStruct* grid, int ie, int ip, FTYPE* a0, FTYPE* a1){

    int p0, p1, p2;
    FTYPE dx1, dx2, dy1, dy2, xm, ym, dr1xdr2, dr1dr2;

    if(ip == 0){

        gridGetElemPoints(grid, ie, &p0, &p1, &p2);

    }else if(ip == 1){

        gridGetElemPoints(grid, ie, &p2, &p0, &p1);

    }else if(ip == 2){

        gridGetElemPoints(grid, ie, &p1, &p2, &p0);

    }else{

        printf("Erro 1 in: gridCalcGradCoef");

    };

    xm = (grid->x[p0] + grid->x[p1] + grid->x[p2])/3;
    ym = (grid->y[p0] + grid->y[p1] + grid->y[p2])/3;

    dx1 = grid->x[p1] - xm;
    dx2 = grid->x[p2] - xm;

    dy1 = grid->y[p1] - ym;
    dy2 = grid->y[p2] - ym;

    dr1xdr2 = dx1*dy2 - dx2*dy1;
    dr1dr2 = dx1*dx2 + dy1*dy2;

    *a0 = 0.5*(dr1dr2 - (dx2*dx2 + dy2*dy2))/dr1xdr2;
    *a1 = 0.5*(dr1dr2 - (dx1*dx1 + dy1*dy1))/dr1xdr2;

}

void gridCalcGradCoef2(struct gridStruct* grid, int ie, int ip, FTYPE* a1, FTYPE* a2){

    int p1, p2, p3;
    FTYPE dx1, dx2, dy1, dy2, xm, ym, dr1xdr2, dr1dr2, yb;

    if(ip == 0){

        gridGetElemPoints(grid, ie, &p1, &p2, &p3);

    }else if(ip == 1){

        gridGetElemPoints(grid, ie, &p3, &p1, &p2);

    }else if(ip == 2){

        gridGetElemPoints(grid, ie, &p2, &p3, &p1);

    }else{

        printf("Erro 1 in: gridCalcGradCoef");

    };

    xm = (grid->x[p1] + grid->x[p2] + grid->x[p3])/3;
    ym = (grid->y[p1] + grid->y[p2] + grid->y[p3])/3;

    dx1 = grid->x[p1] - xm;
    dx2 = grid->x[p2] - xm;

    dy1 = grid->y[p1] - ym;
    dy2 = grid->y[p2] - ym;

    dr1xdr2 = dx1*dy2 - dx2*dy1;
    dr1dr2 = dx1*dx2 + dy1*dy2;

    *a1 = 0.5*(dr1dr2 + (dx2*dx2 + dy2*dy2))/dr1xdr2;
    *a2 = -0.5*(dr1dr2 + (dx1*dx1 + dy1*dy1))/dr1xdr2;

    if(grid->geo == 1){

        yb = (5*grid->y[p1] + 5*grid->y[p2] + 2*grid->y[p3])/12;
        *a1 *= yb;
        *a2 *= yb;

    };

}

void gridCheckGradCoef(struct gridStruct* grid){

    int ii;
    FTYPE a0, a1;
    FTYPE T0 = 100.0, T1 = 200.0, T2 = 400.0, Tm, f;

    Tm = (T0 + T1 + T2)/3;

    printf("\nCheck gradient coeficients:");
    for(ii=0; ii<grid->Ne; ii++){

        if(grid->elem[ii][0]==grid->type2){

            f = 0.0;

            gridCalcGradCoef(grid, ii, 0, &a0, &a1);
            f += a0*(T1 - Tm) + a1*(T2 - Tm);
            //printf("\n%f, %f, %f", f, a0, a1);

            gridCalcGradCoef(grid, ii, 1, &a0, &a1);
            f += a0*(T2 - Tm) + a1*(T0 - Tm);
            //printf("\n%f, %f, %f", f, a0, a1);

            gridCalcGradCoef(grid, ii, 2, &a0, &a1);
            f += a0*(T0 - Tm) + a1*(T1 - Tm);
            //printf("\n%f, %f, %f", f, a0, a1);

            printf("\n%i, %e", ii, f);
        };

    };

}

void gridPrintGradCoef(struct gridStruct* grid){

    int ii;
    FTYPE a0, a1;

    printf("\nGradient coeficients:");
    for(ii=0; ii<grid->Ne; ii++){
        if(grid->elem[ii][0]==grid->type2){
            gridCalcGradCoef(grid, ii, 0, &a0, &a1);
            printf("\n%i, %f, %f,", ii, a0, a1);
            gridCalcGradCoef(grid, ii, 1, &a0, &a1);
            printf(" %f, %f,", a0, a1);
            gridCalcGradCoef(grid, ii, 2, &a0, &a1);
            printf(" %f, %f,", a0, a1);
        };

    };

}

void gridPrintArea(struct gridStruct* grid){

    int ii;

    printf("\nArea:");

    for(ii=0;ii<grid->Ne;ii++){
        printf("\n%i, %f,", ii, gridCalcArea(grid, ii));
    };

}

void gridPrintDualArea(struct gridStruct* grid){

    int ii;

    printf("\nDual area:");

    for(ii=0;ii<grid->Np;ii++){
        printf("\n%i, %f,", ii, grid->dualArea[ii]);
    };

}

void gridCheckOrient(struct gridStruct* grid){

    int ii;
    FTYPE A;

    for(ii=0; ii<grid->Ne; ii++){
        if(grid->elem[ii][0]==grid->type2){
            A = gridCalcArea(grid, ii);
            printf("\n%i, %f", ii, A);
        };
    };

}

void gridSetPointMarkers(struct gridStruct* grid){

    int ii, p0, p1, p2;

    for(ii=0; ii<grid->Ne; ii++){

        if(grid->elem[ii][0]==grid->type2){

            gridGetElemPoints(grid, ii, &p0, &p1, &p2);

            grid->pointMarker[p0] = grid->elem[ii][1];
            grid->pointMarker[p1] = grid->elem[ii][1];
            grid->pointMarker[p2] = grid->elem[ii][1];

        };

    };

    for(ii=0; ii<grid->Ne; ii++){

        if(grid->elem[ii][0]==grid->type1){

            gridGetElemPoints(grid, ii, &p0, &p1, &p2);

            grid->pointMarker[p0] = grid->elem[ii][1];
            grid->pointMarker[p1] = grid->elem[ii][1];

        };

    };

}

void gridPrintPointMarkers(struct gridStruct* grid){

    int ii;

    for(ii=0; ii<grid->Np; ii++){

        printf("\n%i: %i", ii, grid->pointMarker[ii]);

    };

}

void gridGetElemPoints(struct gridStruct* grid, int ie, int* p0, int* p1, int* p2){

    *p0 = grid->elem[ie][2]-1;
    *p1 = grid->elem[ie][3]-1;

    if(grid->elem[ie][0]==grid->type2){

        *p2 = grid->elem[ie][4]-1;

    };

}

int gridElementForPoint(struct gridStruct* grid, FTYPE x, FTYPE y){

    // Precisa  melhorar, ainda não esta completa.

    int ii, inside, p0, p1, p2;

    for(ii=0; ii<grid->Np; ii++){

        if(grid->elem[ii][0]==grid->type2){

            gridGetElemPoints(grid, ii, &p0, &p1, &p2);

            inside = 1;

            if( (x<grid->x[p0]) & (x<grid->x[p1]) & (x<grid->x[p2])){
                inside = 0;

            }else if( (x>grid->x[p0]) & (x>grid->x[p1]) & (x>grid->x[p2])){
                inside = 0;

            }else if( (y<grid->y[p0]) & (y<grid->y[p1]) & (y<grid->y[p2])){
                inside = 0;

            }else if( (y>grid->y[p0]) & (y>grid->y[p1]) & (y>grid->y[p2])){
                inside = 0;

            };

            if(inside){
                break;
            };

        };

    };

    // preciso acrescentar mais checagens para garantir que esteja dentro do triangulo

    if(1-inside){
        ii++;
    };

    return ii;

}

void gridTest(){
    printf("ok");
}
