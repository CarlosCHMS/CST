#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "grid.h"
#include "input.h"

void inputInit(struct inputStruct* input){

    FILE* f1;
    int ii, jj, kk;
    char c;
    char s[100];

    input->meshName = "default.meshsim";

    f1 = fopen("aux.csv", "r");

    ii = 0;
    jj = 0;
    kk = 0;

    //Read input file
    while(c!=EOF){
        c = getc(f1);

        if(c == ','){
            //When the string gets a coma:
            s[jj] = '\0';
            jj = 0;
            //printf("%s, ", s);
            //printf("\ni: %i, k: %i, s: %s", ii, kk, s);
            if(ii==0){

                if(kk==1){
                    input->c = strtod(s, NULL);
                };

            }else if(ii==1){

                if(kk==1){
                    input->N = strtod(s, NULL);
                };

            }else if(ii==2){

                if(kk==1){
                    input->Nmarkers = strtod(s, NULL);
                };

                input->markers = (int*)malloc(input->Nmarkers*sizeof(int));
                input->fInputs = (FTYPE*)malloc(input->Nmarkers*sizeof(FTYPE));
                input->fInputs1 = (FTYPE*)malloc(input->Nmarkers*sizeof(FTYPE));
                input->fInputs4 = (FTYPE*)malloc(input->Nmarkers*sizeof(FTYPE));
                input->types = (char*)malloc(input->Nmarkers*sizeof(char));

            }else if(ii==3){

                if(kk>0){
                    input->markers[kk-1] = strtod(s, NULL);
                };

            }else if(ii==4){

                if(kk>0){
                    input->fInputs[kk-1] = strtod(s, NULL);
                };

            }else if(ii==5){

                if(kk>0){
                    input->types[kk-1] = s[0];
                };

            }else if(ii==6){

                if(kk>0){

                    //sprintf(a, "%s", s);
                    //printf("\noi:%s.", a);
                    input->meshName = strdup(s);
                    //printf("\noi:%s.", input->meshName);
                };

            }else if(ii==7){

                if(kk>0){

                    input->Nsave = strtod(s, NULL);

                };

            }else if(ii==8){

                if(kk>0){

                    input->outName = strdup(s);

                };

            }else if(ii==9){

                if(kk>0){
                    //input->Lref = strtod(s, NULL);
                };

            }else if(ii==10){

                if(kk>0){
                    input->fInputs1[kk-1] = strtod(s, NULL);
                };

            }else if(ii==11){

                if(kk>0){
                    input->k = strtod(s, NULL);
                };

            }else if(ii==10){

                if(kk>0){
                    input->fInputs4[kk-1] = strtod(s, NULL);
                };

            };

            kk ++;
        }else if(c == '\n'){
            ii++;
            kk = 0;
        }else{
            if(c != ' '){
                s[jj] = c;
                jj++;
            };
        };
        //printf("\n%i", ii);
    };
    //printf("\noiii");

    fclose(f1);

}

void inputPrintParameters(struct inputStruct* input){

    int ii;

    printf("\nParameters:");
    printf("\nc: %f", input->c);
    printf("\nN: %i", input->N);
    printf("\nNmarkers: %i", input->Nmarkers);

    printf("\nBoundary markers: ");
    for(ii=0; ii<input->Nmarkers; ii++){
        printf("%i, ", input->markers[ii]);
    };

    printf("\nBoundary inputs: ");
    for(ii=0; ii<input->Nmarkers; ii++){
        printf("%f, ", input->fInputs[ii]);
    };

    printf("\nBoundary inputs 1: ");
    for(ii=0; ii<input->Nmarkers; ii++){
        printf("%f, ", input->fInputs1[ii]);
    };

    printf("\nBoundary types: ");
    for(ii=0; ii<input->Nmarkers; ii++){
        printf("%c, ", input->types[ii]);
    };

    printf("\nmeshName: %s", input->meshName);
    printf("\nNsave: %i", input->Nsave);
    printf("\noutName: %s", input->outName);
    //printf("\nLref: %f", input->Lref);

    printf("\nBoundary inputs 4: ");
    for(ii=0; ii<input->Nmarkers; ii++){
        printf("%f, ", input->fInputs4[ii]);
    };



};
