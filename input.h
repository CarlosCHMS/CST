#ifndef INPUT_H_INCLUDED
#define INPUT_H_INCLUDED

    typedef double FTYPE;

    struct inputStruct{

        //FTYPE c;
        //FTYPE k;
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
        FTYPE* C;
        FTYPE* K;
        int geo;

    };

    void inputInit(struct inputStruct* input);

    void inputPrintParameters(struct inputStruct* input);

#endif // INPUT_H_INCLUDED
