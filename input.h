#ifndef INPUT_H_INCLUDED
#define INPUT_H_INCLUDED

    struct inputStruct{

        FTYPE c;
        FTYPE k;
        int N;
        int Nmarkers;
        int* markers;
        FTYPE* fInputs;
        FTYPE* fInputs1;
        FTYPE* fInputs4;
        char* types;
        char* meshName;
        int Nsave;
        char* outName;

    };

    void inputInit(struct inputStruct* input);

    void inputPrintParameters(struct inputStruct* input);

#endif // INPUT_H_INCLUDED
