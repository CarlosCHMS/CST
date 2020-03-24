#ifndef BOUDARY_H_INCLUDED
#define BOUDARY_H_INCLUDED

    struct boundaryStruct{

        int* markers;
        FTYPE* fInputs;
        FTYPE* fInputs1;
        FTYPE* fInputs4;
        FTYPE* C;
        FTYPE* K;
        char* types;

        int* pointPropag;
        int* elemIndex;

        struct gridStruct* grid;

    };

    void boundaryAlloc(struct boundaryStruct* bound);

    void boundaryInit(struct boundaryStruct* bound, struct gridStruct* grid, struct inputStruct* input);

    void boundarySetPointPropag(struct boundaryStruct* bound);

    void boundaryInitIndex(struct boundaryStruct* bound);

    void boundaryPrintTypes(struct boundaryStruct* bound);

    int boundaryElemIsDomain(struct boundaryStruct* bound, int ie);

    int boundaryElemIsTemp(struct boundaryStruct* bound, int ie);

    int boundaryElemIsHeat(struct boundaryStruct* bound, int ie);

    //int boundaryElemIsConvec(struct boundaryStruct* bound, int ie);

    FTYPE boundaryGetElemInput(struct boundaryStruct* bound, int ii);

    FTYPE boundaryGetElemInput2(struct boundaryStruct* bound, int ii);

#endif // BOUDARY_H_INCLUDED
