#ifndef I2PHASE_H
#define I2PHASE_H

typedef long long ll;

class Util;
class CubieCube;
class CoordCube;
class Search;

class Util {
public:
    static const unsigned char Ux1 = 0;
    static const unsigned char Ux2 = 1;
    static const unsigned char Ux3 = 2;
    static const unsigned char Rx1 = 3;
    static const unsigned char Rx2 = 4;
    static const unsigned char Rx3 = 5;
    static const unsigned char Fx1 = 6;
    static const unsigned char Fx2 = 7;
    static const unsigned char Fx3 = 8;
    static const unsigned char Dx1 = 9;
    static const unsigned char Dx2 = 10;
    static const unsigned char Dx3 = 11;
    static const unsigned char Lx1 = 12;
    static const unsigned char Lx2 = 13;
    static const unsigned char Lx3 = 14;
    static const unsigned char Bx1 = 15;
    static const unsigned char Bx2 = 16;
    static const unsigned char Bx3 = 17;

    static const unsigned char U1 = 0;
    static const unsigned char U2 = 1;
    static const unsigned char U3 = 2;
    static const unsigned char U4 = 3;
    static const unsigned char U5 = 4;
    static const unsigned char U6 = 5;
    static const unsigned char U7 = 6;
    static const unsigned char U8 = 7;
    static const unsigned char U9 = 8;
    static const unsigned char R1 = 9;
    static const unsigned char R2 = 10;
    static const unsigned char R3 = 11;
    static const unsigned char R4 = 12;
    static const unsigned char R5 = 13;
    static const unsigned char R6 = 14;
    static const unsigned char R7 = 15;
    static const unsigned char R8 = 16;
    static const unsigned char R9 = 17;
    static const unsigned char F1 = 18;
    static const unsigned char F2 = 19;
    static const unsigned char F3 = 20;
    static const unsigned char F4 = 21;
    static const unsigned char F5 = 22;
    static const unsigned char F6 = 23;
    static const unsigned char F7 = 24;
    static const unsigned char F8 = 25;
    static const unsigned char F9 = 26;
    static const unsigned char D1 = 27;
    static const unsigned char D2 = 28;
    static const unsigned char D3 = 29;
    static const unsigned char D4 = 30;
    static const unsigned char D5 = 31;
    static const unsigned char D6 = 32;
    static const unsigned char D7 = 33;
    static const unsigned char D8 = 34;
    static const unsigned char D9 = 35;
    static const unsigned char L1 = 36;
    static const unsigned char L2 = 37;
    static const unsigned char L3 = 38;
    static const unsigned char L4 = 39;
    static const unsigned char L5 = 40;
    static const unsigned char L6 = 41;
    static const unsigned char L7 = 42;
    static const unsigned char L8 = 43;
    static const unsigned char L9 = 44;
    static const unsigned char B1 = 45;
    static const unsigned char B2 = 46;
    static const unsigned char B3 = 47;
    static const unsigned char B4 = 48;
    static const unsigned char B5 = 49;
    static const unsigned char B6 = 50;
    static const unsigned char B7 = 51;
    static const unsigned char B8 = 52;
    static const unsigned char B9 = 53;

    static const unsigned char U = 0;
    static const unsigned char R = 1;
    static const unsigned char F = 2;
    static const unsigned char D = 3;
    static const unsigned char L = 4;
    static const unsigned char B = 5;

    static const unsigned char cornerFacelet[8][3];
    static const unsigned char edgeFacelet[12][2];

    static int Cnk[13][13];
    static int fact[14];
    static int permMult[24][24];
    static const char* move2str[18];
    static int preMove[9];
    static int ud2std[10];
    static int std2ud[18];
    static bool ckmv2[11][10];

    static void toCubieCube(char *f, CubieCube &ccRet);
    static char* toFaceCube(CubieCube cc);
    static int getNParity(int idx, int n);
    static unsigned char setVal(int val0, int val, bool isEdge);
    static int getVal(int val0, bool isEdge);
    static void set8Perm(unsigned char *arr, int idx, bool isEdge);
    static int get8Perm(unsigned char *arr, bool isEdge);
    static void setNPerm(unsigned char *arr, int idx, int n, bool isEdge);
    static int getNPerm(unsigned char *arr, int n, bool isEdge);
    static int getComb(unsigned char *arr, int mask, bool isEdge, int end);
    static void setComb(unsigned char *arr, int idx, int mask, bool isEdge, int end);

    static void init();
    static bool inited;
};

class CubieCube {
public:
    static CubieCube CubeSym[16];
    static CubieCube moveCube[18];
    static ll moveCubeSym[18];
    static int firstMoveSym[48];

    static int preMove[9];

    static int SymInv[16];
    static int SymMult[16][16];
    static int SymMove[16][18];
    static int SymMultInv[16][16];
    static int Sym8Mult[8 * 8];
    static int Sym8Move[8 * 18];
    static int Sym8MultInv[8 * 8];
    static int SymMoveUD[16][10];

    static int FlipS2R[336];
    static int TwistS2R[324];
    static int EPermS2R[2768];
    static int UDSliceFlipS2R[64430];

    static unsigned char e2c[16];

    static int MtoEPerm[40320];

    static int FlipSlice2UDSliceFlip[336 * 495];

    static int FlipR2S[2048];// = new char[2048];
    static int TwistR2S[2187];// = new char[2187];
    static int EPermR2S[40320];// = new char[40320];
    static int FlipS2RF[336 * 8];
    static int TwistS2RF[324 * 8];

    static int SymStateTwist[324];
    static int SymStateFlip[336];
    static int SymStatePerm[2768];
    static int SymStateUDSliceFlip[64430];

    static CubieCube urf1;
    static CubieCube urf2;
    static unsigned char urfMove[6][18];

    static unsigned char ica[8];
    static unsigned char iea[12];
    unsigned char ca[8];
    unsigned char ea[12];

    CubieCube();
    CubieCube(int cperm, int twist, int eperm, int flip);
    CubieCube(const CubieCube &c);

    bool equalsCorn(CubieCube c);
    bool equalsEdge(CubieCube c);
    void copy(const CubieCube &c);
    void construct(int cperm, int twist, int eperm, int flip);
    void invCubieCube();

    static void CornMult(CubieCube a, CubieCube b, CubieCube &prod);
    static void EdgeMult(CubieCube a, CubieCube b, CubieCube &prod);
    static void CornConjugate(CubieCube a, int idx, CubieCube &b);
    static void EdgeConjugate(CubieCube a, int idx, CubieCube &b);

    void URFConjugate();
    int getFlip();
    void setFlip(int idx);
    int getFlipSym();
    int getTwist();
    void setTwist(int idx);
    int getTwistSym();
    int getUDSlice();
    void setUDSlice(int idx);
    int getU4Comb();
    int getD4Comb();
    int getCPerm();
    void setCPerm(int idx);
    int getCPermSym();
    int getEPerm();
    void setEPerm(int idx);
    int getEPermSym();
    int getMPerm();
    void setMPerm(int idx);
    int getCComb();
    void setCComb(int idx);
    int verify();
    ll selfSymmetry();
    void setUDSliceFlip(int idx);
    int getUDSliceFlip();
    int getUDSliceFlipSym();

    static void initMove();
    char* toString();
    static void initSym();
    static void initFlipSym2Raw();
    static void initTwistSym2Raw();
    static int Perm2Comb[2768];
    static void initPermSym2Raw();
    static void initUDSliceFlipSym2Raw();
};

class CoordCube {
public:
    static const int N_MOVES = 18;
    static const int N_MOVES2 = 10;

    static const int N_SLICE = 495;
    static const int N_TWIST = 2187;
    static const int N_TWIST_SYM = 324;
    static const int N_FLIP = 2048;
    static const int N_FLIP_SYM = 336;
    static const int N_PERM = 40320;
    static const int N_PERM_SYM = 2768;
    static const int N_MPERM = 24;
    static const int N_COMB = 70;

    static int UDSliceMove[N_SLICE][N_MOVES];
    static int TwistMove[N_TWIST_SYM][N_MOVES];
    static int FlipMove[N_FLIP_SYM][N_MOVES];
    static int UDSliceConj[N_SLICE][8];
    static int UDSliceTwistPrun[N_SLICE * N_TWIST_SYM / 8 + 1];
    static int UDSliceFlipPrun[N_SLICE * N_FLIP_SYM / 8];
    static int TwistFlipPrun[N_FLIP * N_TWIST_SYM / 8];

    static int CPermMove[N_PERM_SYM][N_MOVES];
    static int EPermMove[N_PERM_SYM][N_MOVES2];
    static int MPermMove[N_MPERM][N_MOVES2];
    static int MPermConj[N_MPERM][16];
    static int CCombMove[N_COMB][N_MOVES];
    static int CCombConj[N_COMB][16];
    static int MCPermPrun[N_MPERM * N_PERM_SYM / 8];
    static int MEPermPrun[N_MPERM * N_PERM_SYM / 8];
    static int EPermCCombPrun[N_COMB * N_PERM_SYM / 8];

    static void init();
    static void setPruning(int *table, int index, int value);
    static int getPruning(int *table, int index);
    static void initUDSliceMoveConj();
    static void initFlipMove();
    static void initTwistMove();
    static void initCPermMove();
    static void initEPermMove();
    static void initMPermMoveConj();
    static void initCombMoveConj();
    static void initTwistFlipPrun();
    static void initRawSymPrun(int *PrunTable, const int INV_DEPTH, int (&RawMove)[N_SLICE][N_MOVES], int (&RawConj)[N_SLICE][8], int (&SymMove)[N_TWIST_SYM][N_MOVES], int *SymState, const int PrunFlag, const int N_RAW, const int N_SYM, const int NN_MOVES);
    static void initRawSymPrun(int *PrunTable, const int INV_DEPTH, int (&RawMove)[N_SLICE][N_MOVES], int (&RawConj)[N_SLICE][8], int (&SymMove)[N_FLIP_SYM][N_MOVES], int *SymState, const int PrunFlag, const int N_RAW, const int N_SYM, const int NN_MOVES);
    static void initRawSymPrun(int *PrunTable, const int INV_DEPTH, int (&RawMove)[N_MPERM][N_MOVES2], int (&RawConj)[N_MPERM][16], int (&SymMove)[N_PERM_SYM][N_MOVES2], int *SymState, const int PrunFlag, const int N_RAW, const int N_SYM, const int NN_MOVES);
    static void initRawSymPrun(int *PrunTable, const int INV_DEPTH, int (&RawMove)[N_MPERM][N_MOVES2], int (&RawConj)[N_MPERM][16], int (&SymMove)[N_PERM_SYM][N_MOVES], int *SymState, const int PrunFlag, const int N_RAW, const int N_SYM, const int NN_MOVES);
    static void initRawSymPrun(int *PrunTable, const int INV_DEPTH, int (&RawMove)[N_COMB][N_MOVES], int (&RawConj)[N_COMB][16], int (&SymMove)[N_PERM_SYM][N_MOVES2], int *SymState, const int PrunFlag, const int N_RAW, const int N_SYM, const int NN_MOVES);
    static void initSliceTwistPrun();
    static void initSliceFlipPrun();
    static void initMEPermPrun();
    static void initMCPermPrun();
    static void initPermCombPrun();

    int twist;
    int tsym;
    int flip;
    int fsym;
    int slice;
    int prun;

    CoordCube();
    void set(CoordCube node);
    void calcPruning(bool isPhase1);
    void set(CubieCube cc);
    int doMovePrun(CoordCube cc, int m, bool isPhase1);
};

class Search {
public:
    static const bool USE_TWIST_FLIP_PRUN = true;
    static const int EXTRA_PRUN_LEVEL = 0;

    static const bool USF_FULL_PRUN = EXTRA_PRUN_LEVEL > 0;
    static const bool USF_HUGE_PRUN = EXTRA_PRUN_LEVEL > 1;

    static const bool TRY_PRE_MOVE = true;
    static const bool TRY_INVERSE = true;
    static const bool TRY_THREE_AXES = true;

    static const int MAX_DEPTH2 = 13;

    static const int PRE_IDX_MAX = TRY_PRE_MOVE ? 9 : 1;

    static bool inited;

    int move[31];

    int corn0[6][PRE_IDX_MAX];
    int ud8e0[6][PRE_IDX_MAX];

    CoordCube nodeUD[21];
    CoordCube nodeRL[21];
    CoordCube nodeFB[21];

    CoordCube node0[6][PRE_IDX_MAX];

    char f[54];

    ll selfSym;
    int preIdxMax;
    int conjMask;
    int urfIdx;
    int preIdx;
    int length1;
    int depth1;
    int maxDep2;
    int sol;
    char *solution;
    ll probe;
    ll probeMax;
    ll probeMin;
    int verbose;
    CubieCube cc;

    bool isRec;

    static const int USE_SEPARATOR = 0x1;
    static const int INVERSE_SOLUTION = 0x2;
    static const int APPEND_LENGTH = 0x4;
    static const int OPTIMAL_SOLUTION = 0x8;

    Search();

    char *Solution(char *facelets, int maxDepth, ll probeMax, ll probeMin, int verbose);
    void initSearch();
    char *next(ll probeMax, ll probeMin, int verbose);
    bool isInited();
    ll numberOfProbes();
    int length();
    static void init();
    int verify(char *facelets);
    char *search();
    int phase1(CoordCube node, ll ssym, int maxl, int lm);
    int initPhase2();
    int phase2(int eidx, int esym, int cidx, int csym, int mid, int maxl, int depth, int lm);
    char *solutionToString();
};
#endif // I2PHASE_H
