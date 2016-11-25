#include "I2Phase.h"
#include <cstring>

const unsigned char Util::cornerFacelet[8][3] = {
    { U9, R1, F3 }, { U7, F1, L3 }, { U1, L1, B3 }, { U3, B1, R3 },
    { D3, F9, R7 }, { D1, L9, F7 }, { D7, B9, L7 }, { D9, R9, B7 }
};
const unsigned char Util::edgeFacelet[12][2] = {
    { U6, R2 }, { U8, F2 }, { U4, L2 }, { U2, B2 }, { D6, R8 }, { D2, F8 },
    { D4, L8 }, { D8, B8 }, { F6, R4 }, { F4, L6 }, { B6, L4 }, { B4, R6 }
};
const char* Util::move2str[18] = {
    "U ", "U2", "U'", "R ", "R2", "R'", "F ", "F2", "F'",
    "D ", "D2", "D'", "L ", "L2", "L'", "B ", "B2", "B'"
};
int Util::preMove[9] = { -1, Rx1, Rx3, Fx1, Fx3, Lx1, Lx3, Bx1, Bx3};
int Util::ud2std[10] = {Ux1, Ux2, Ux3, Rx2, Fx2, Dx1, Dx2, Dx3, Lx2, Bx2};
int Util::Cnk[13][13];
int Util::fact[14];
int Util::permMult[24][24];
int Util::std2ud[18];
bool Util::ckmv2[11][10];
bool Util::inited = false;

void Util::toCubieCube(char *f, CubieCube &ccRet) {
    unsigned char ori;
    for (int i = 0; i < 8; i++)
        ccRet.ca[i] = 0;
    for (int i = 0; i < 12; i++)
        ccRet.ea[i] = 0;
    unsigned char col1, col2;
    for (unsigned char i = 0; i < 8; i++) {
        for (ori = 0; ori < 3; ori++)
            if (f[cornerFacelet[i][ori]] == U || f[cornerFacelet[i][ori]] == D)
                break;
        col1 = f[cornerFacelet[i][(ori + 1) % 3]];
        col2 = f[cornerFacelet[i][(ori + 2) % 3]];

        for (unsigned char j = 0; j < 8; j++) {
            if (col1 == cornerFacelet[j][1] / 9 && col2 == cornerFacelet[j][2] / 9) {
                ccRet.ca[i] = (unsigned char) (ori % 3 << 3 | j);
                break;
            }
        }
    }
    for (unsigned char i = 0; i < 12; i++) {
        for (unsigned char j = 0; j < 12; j++) {
            if (f[edgeFacelet[i][0]] == edgeFacelet[j][0] / 9
                    && f[edgeFacelet[i][1]] == edgeFacelet[j][1] / 9) {
                ccRet.ea[i] = (unsigned char) (j << 1);
                break;
            }
            if (f[edgeFacelet[i][0]] == edgeFacelet[j][1] / 9
                    && f[edgeFacelet[i][1]] == edgeFacelet[j][0] / 9) {
                ccRet.ea[i] = (unsigned char) (j << 1 | 1);
                break;
            }
        }
    }
}

char *Util::toFaceCube(CubieCube cc) {
    static char f[54];
    char ts[6] = {'U', 'R', 'F', 'D', 'L', 'B'};
    for (int i = 0; i < 54; i++) {
        f[i] = ts[i / 9];
    }
    for (unsigned char c = 0; c < 8; c++) {
        int j = cc.ca[c] & 0x7;
        int ori = cc.ca[c] >> 3;
        for (unsigned char n = 0; n < 3; n++)
            f[cornerFacelet[c][(n + ori) % 3]] = ts[cornerFacelet[j][n] / 9];
    }
    for (unsigned char e = 0; e < 12; e++) {
        int j = cc.ea[e] >> 1;
        int ori = cc.ea[e] & 1;
        for (unsigned char n = 0; n < 2; n++)
            f[edgeFacelet[e][(n + ori) % 2]] = ts[edgeFacelet[j][n] / 9];
    }
    return f;
}

int Util::getNParity(int idx, int n) {
    int p = 0;
    for (int i = n - 2; i >= 0; i--) {
        p ^= idx % (n - i);
        idx /= (n - i);
    }
    return p & 1;
}

unsigned char Util::setVal(int val0, int val, bool isEdge) {
    return (unsigned char) (isEdge ? (val << 1 | val0 & 1) : (val | val0 & 0xf8));
}

int Util::getVal(int val0, bool isEdge) {
    return isEdge ? val0 >> 1 : val0 & 7;
}

void Util::set8Perm(unsigned char *arr, int idx, bool isEdge) {
    int val = 0x76543210;
    for (int i = 0; i < 7; i++) {
        int p = fact[7 - i];
        int v = idx / p;
        idx -= v * p;
        v <<= 2;
        arr[i] = setVal(arr[i], (val >> v & 0x7), isEdge);
        int m = (1 << v) - 1;
        val = val & m | val >> 4 & ~m;
    }
    arr[7] = setVal(arr[7], val, isEdge);
}

int Util::get8Perm(unsigned char *arr, bool isEdge) {
    int idx = 0;
    int val = 0x76543210;
    for (int i = 0; i < 7; i++) {
        int v = getVal(arr[i], isEdge) << 2;
        idx = (8 - i) * idx + (val >> v & 0x7);
        val -= 0x11111110 << v;
    }
    return idx;
}

void Util::setNPerm(unsigned char *arr, int idx, int n, bool isEdge) {
    arr[n - 1] = setVal(arr[n - 1], 0, isEdge);
    for (int i = n - 2; i >= 0; i--) {
        int arri = idx % (n - i);
        arr[i] = setVal(arr[i], arri, isEdge);
        idx /= (n - i);
        for (int j = i + 1; j < n; j++) {
            int arrj = getVal(arr[j], isEdge);
            if (arrj >= arri) {
                arr[j] = setVal(arr[j], ++arrj, isEdge);
            }
        }
    }
}

int Util::getNPerm(unsigned char *arr, int n, bool isEdge) {
    int idx = 0;
    for (int i = 0; i < n; i++) {
        idx *= (n - i);
        int arri = getVal(arr[i], isEdge);
        for (int j = i + 1; j < n; j++) {
            int arrj = getVal(arr[j], isEdge);
            if (arrj < arri) {
                idx++;
            }
        }
    }
    return idx;
}

int Util::getComb(unsigned char *arr, int mask, bool isEdge, int end) {
    int idxC = 0, idxP = 0, r = 4, val = 0x0123;
    for (int i = end; i >= 0; i--) {
        int perm = getVal(arr[i], isEdge);
        if ((perm & 0xc) == mask) {
            int v = (perm & 3) << 2;
            idxP = r * idxP + (val >> v & 0xf);
            val -= 0x0111 >> (12 - v);
            idxC += Cnk[i][r--];
        }
    }
    return idxP << 9 | Cnk[end + 1][4] - 1 - idxC;
}

void Util::setComb(unsigned char *arr, int idx, int mask, bool isEdge, int end) {
    int r = 4, fill = end, val = 0x0123;
    int idxC = Cnk[end + 1][4] - 1 - (idx & 0x1ff);
    int idxP = idx >> 9;
    for (int i = end; i >= 0; i--) {
        if (idxC >= Cnk[i][r]) {
            idxC -= Cnk[i][r--];
            int p = fact[r];
            int v = idxP / p << 2;
            idxP %= p;
            arr[i] = setVal(arr[i], val >> v & 3 | mask, isEdge);
            int m = (1 << v) - 1;
            val = val & m | val >> 4 & ~m;
        } else {
            if ((fill & 0xc) == mask) {
                fill -= 4;
            }
            arr[i] = setVal(arr[i], fill--, isEdge);
        }
    }
}

void Util::init() {
    if (inited) return; else inited = true;
    for (int i = 0; i < 10; i++) {
        std2ud[ud2std[i]] = i;
    }
    for (int i = 0; i < 10; i++) {
        int ix = ud2std[i];
        for (int j = 0; j < 10; j++) {
            int jx = ud2std[j];
            ckmv2[i][j] = (ix / 3 == jx / 3) || ((ix / 3 % 3 == jx / 3 % 3) && (ix >= jx));
        }
        ckmv2[10][i] = false;
    }
    fact[0] = 1;
    for (int i = 0; i < 13; i++) {
        Cnk[i][0] = Cnk[i][i] = 1;
        fact[i + 1] = fact[i] * (i + 1);
        for (int j = 1; j < i; j++) {
            Cnk[i][j] = Cnk[i - 1][j - 1] + Cnk[i - 1][j];
        }
    }
    unsigned char arr1[4];
    unsigned char arr2[4];
    unsigned char arr3[4];
    memset(arr1, 0, sizeof(arr1));
    memset(arr2, 0, sizeof(arr2));
    for (int i = 0; i < 24; i++) {
        setNPerm(arr1, i, 4, false);
        for (int j = 0; j < 24; j++) {
            setNPerm(arr2, j, 4, false);
            for (int k = 0; k < 4; k++) {
                arr3[k] = arr1[arr2[k]];
            }
            permMult[i][j] = getNPerm(arr3, 4, false);
        }
    }
}
