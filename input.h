#ifndef INPUT_H_INCLUDED
#define INPUT_H_INCLUDED

    struct inputStruct{

        FTYPE c;
        int N;
        int Nmarkers;
        int* markers;
        FTYPE* fInputs;
        char* types;
        char* meshName;
        int saveStep;
        char* outName;
        FTYPE Lref;

    };

    void inputInit(struct inputStruct* input);

    void inputPrintParameters(struct inputStruct* input);

#endif // INPUT_H_INCLUDED
