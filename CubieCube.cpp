#include "I2Phase.h"
#include <cstring>
#include <cassert>
#include <algorithm>
#include <cstdio>
#include <cstdlib>

int CubieCube::preMove[9] = { -1, Util::Rx1, Util::Rx3, Util::Fx1, Util::Fx3, Util::Lx1, Util::Lx3, Util::Bx1, Util::Bx3};
unsigned char CubieCube::e2c[16] = {0, 0, 0, 0, 1, 3, 1, 3, 1, 3, 1, 3, 0, 0, 0, 0};
CubieCube CubieCube::urf1(2531, 1373, 67026819, 1367);
CubieCube CubieCube::urf2(2089, 1906, 322752913, 2040);
unsigned char CubieCube::urfMove[6][18] = {
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17},
    {6, 7, 8, 0, 1, 2, 3, 4, 5, 15, 16, 17, 9, 10, 11, 12, 13, 14},
    {3, 4, 5, 6, 7, 8, 0, 1, 2, 12, 13, 14, 15, 16, 17, 9, 10, 11},
    {2, 1, 0, 5, 4, 3, 8, 7, 6, 11, 10, 9, 14, 13, 12, 17, 16, 15},
    {8, 7, 6, 2, 1, 0, 5, 4, 3, 17, 16, 15, 11, 10, 9, 14, 13, 12},
    {5, 4, 3, 8, 7, 6, 2, 1, 0, 14, 13, 12, 17, 16, 15, 11, 10, 9}
};
unsigned char CubieCube::ica[8] = {0, 1, 2, 3, 4, 5, 6, 7};
unsigned char CubieCube::iea[12] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22};

CubieCube CubieCube::CubeSym[16];
CubieCube CubieCube::moveCube[18];
ll CubieCube::moveCubeSym[18];
int CubieCube::firstMoveSym[48];
int CubieCube::SymInv[16];
int CubieCube::SymMult[16][16];
int CubieCube::SymMove[16][18];
int CubieCube::SymMultInv[16][16];
int CubieCube::Sym8Mult[8 * 8];
int CubieCube::Sym8Move[8 * 18];
int CubieCube::Sym8MultInv[8 * 8];
int CubieCube::SymMoveUD[16][10];
int CubieCube::FlipS2R[336];
int CubieCube::TwistS2R[324];
int CubieCube::EPermS2R[2768];
int CubieCube::UDSliceFlipS2R[64430];
int CubieCube::MtoEPerm[40320];

int CubieCube::FlipSlice2UDSliceFlip[336 * 495];

int CubieCube::FlipR2S[2048];// = new unsigned char[2048];
int CubieCube::TwistR2S[2187];// = new unsigned char[2187];
int CubieCube::EPermR2S[40320];// = new unsigned char[40320];
int CubieCube::FlipS2RF[336 * 8];
int CubieCube::TwistS2RF[324 * 8];

int CubieCube::SymStateTwist[324];
int CubieCube::SymStateFlip[336];
int CubieCube::SymStatePerm[2768];
int CubieCube::SymStateUDSliceFlip[64430];
int CubieCube::Perm2Comb[2768];

CubieCube::CubieCube() {
    Util::init();
    memcpy((*this).ca, CubieCube::ica, sizeof((*this).ca));
    memcpy((*this).ea, CubieCube::iea, sizeof((*this).ea));
}

CubieCube::CubieCube(int cperm, int twist, int eperm, int flip) {
    Util::init();
    memcpy((*this).ca, CubieCube::ica, sizeof((*this).ca));
    memcpy((*this).ea, CubieCube::iea, sizeof((*this).ea));
    (*this).setCPerm(cperm);
    (*this).setTwist(twist);
    Util::setNPerm(ea, eperm, 12, true);
    (*this).setFlip(flip);
}

CubieCube::CubieCube(const CubieCube &c) {
    copy(c);
}

bool CubieCube::equalsCorn(CubieCube c) {
    for (int i = 0; i < 8; i++) {
        if (ca[i] != c.ca[i]) {
            return false;
        }
    }
    return true;
}

bool CubieCube::equalsEdge(CubieCube c) {
    for (int i = 0; i < 12; i++) {
        if (ea[i] != c.ea[i]) {
            return false;
        }
    }
    return true;
}

void CubieCube::copy(const CubieCube &c) {
    for (int i = 0; i < 8; i++) {
        (*this).ca[i] = c.ca[i];
    }
    for (int i = 0; i < 12; i++) {
        (*this).ea[i] = c.ea[i];
    }
}

void CubieCube::construct(int cperm, int twist, int eperm, int flip) {
    CubieCube c(cperm, twist, eperm, flip);
    (*this).copy(c);
}

void CubieCube::invCubieCube() {
    CubieCube temps;
    for (int edge = 0; edge < 12; edge++) {
        temps.ea[ea[edge] >> 1] = (int) (edge << 1 | ea[edge] & 1);
    }
    for (int corn = 0; corn < 8; corn++) {
        int ori = ca[corn] >> 3;
        ori = 4 >> ori & 3; //0->0, 1->2, 2->1
        temps.ca[ca[corn] & 0x7] = (int) (corn | ori << 3);
    }
    copy(temps);
}

void CubieCube::CornMult(CubieCube a, CubieCube b, CubieCube &prod) {
    for (int corn = 0; corn < 8; corn++) {
        int oriA = a.ca[b.ca[corn] & 7] >> 3;
        int oriB = b.ca[corn] >> 3;
        int ori = oriA;
        ori += (oriA < 3) ? oriB : 6 - oriB;
        ori %= 3;
        if ((oriA >= 3) ^ (oriB >= 3)) {
            ori += 3;
        }
        prod.ca[corn] = (int) (a.ca[b.ca[corn] & 7] & 7 | ori << 3);
    }
}

void CubieCube::EdgeMult(CubieCube a, CubieCube b, CubieCube &prod) {
    for (int ed = 0; ed < 12; ed++) {
        prod.ea[ed] = (int) (a.ea[b.ea[ed] >> 1] ^ (b.ea[ed] & 1));
    }
}

void CubieCube::CornConjugate(CubieCube a, int idx, CubieCube &b) {
    CubieCube sinv = CubeSym[SymInv[idx]];
    CubieCube s = CubeSym[idx];
    for (int corn = 0; corn < 8; corn++) {
        int oriA = sinv.ca[a.ca[s.ca[corn] & 7] & 7] >> 3;
        int oriB = a.ca[s.ca[corn] & 7] >> 3;
        int ori = (oriA < 3) ? oriB : (3 - oriB) % 3;
        b.ca[corn] = (int) (sinv.ca[a.ca[s.ca[corn] & 7] & 7] & 7 | ori << 3);
    }
}

void CubieCube::EdgeConjugate(CubieCube a, int idx, CubieCube &b) {
    CubieCube sinv = CubeSym[SymInv[idx]];
    CubieCube s = CubeSym[idx];
    for (int ed = 0; ed < 12; ed++) {
        b.ea[ed] = (int) (sinv.ea[a.ea[s.ea[ed] >> 1] >> 1] ^ (a.ea[s.ea[ed] >> 1] & 1) ^ (s.ea[ed] & 1));
    }
}

void CubieCube::URFConjugate() {
    CubieCube temps;
    CornMult(urf2, (*this), temps);
    CornMult(temps, urf1, (*this));
    EdgeMult(urf2, (*this), temps);
    EdgeMult(temps, urf1, (*this));
}

int CubieCube::getFlip() {
    int idx = 0;
    for (int i = 0; i < 11; i++) {
        idx = idx << 1 | ea[i] & 1;
    }
    return idx;
}

void CubieCube::setFlip(int idx) {
    int parity = 0;
    for (int i = 10; i >= 0; i--) {
        int val = idx & 1;
        ea[i] = (int) (ea[i] & 0xfe | val);
        parity ^= val;
        idx >>= 1;
    }
    ea[11] = (int) (ea[11] & 0xfe | parity);
}

int CubieCube::getFlipSym() {
    /*
    if (FlipR2S != null) {
        return FlipR2S[getFlip()];
    }
    */
    CubieCube temps;
    for (int k = 0; k < 16; k += 2) {
        EdgeConjugate((*this), SymInv[k], temps);
        int c = temps.getFlip();
        int idx = std::lower_bound(FlipS2R, FlipS2R + 336, c) - FlipS2R;
        if (idx < 336)
            if (FlipS2R[idx] == c) {
                return idx << 3 | k >> 1;
            }
    }
    assert(0);
    return 0;
}

int CubieCube::getTwist() {
    int idx = 0;
    for (int i = 0; i < 7; i++) {
        idx += (idx << 1) + (ca[i] >> 3);
    }
    return idx;
}

void CubieCube::setTwist(int idx) {
    int twst = 0;
    for (int i = 6; i >= 0; i--) {
        int val = idx % 3;
        ca[i] = (int) (ca[i] & 0x7 | val << 3);
        twst += val;
        idx /= 3;
    }
    ca[7] = (int) (ca[7] & 0x7 | ((15 - twst) % 3) << 3);
}

int CubieCube::getTwistSym() {
    /*
    if (TwistR2S != null) {
        return TwistR2S[getTwist()];
    }
    */
    CubieCube temps;
    for (int k = 0; k < 16; k += 2) {
        CornConjugate((*this), SymInv[k], temps);
        int c = temps.getTwist();
        int idx = std::lower_bound(TwistS2R, TwistS2R + 324, c) - TwistS2R;
        if (idx < 324)
            if (TwistS2R[idx] == c) {
                return idx << 3 | k >> 1;
            }
    }
    assert(0);
    return 0;
}

int CubieCube::getUDSlice() {
    return Util::getComb(ea, 8, true, 11);
}

void CubieCube::setUDSlice(int idx) {
    Util::setComb(ea, idx, 8, true, 11);
}

int CubieCube::getU4Comb() {
    return Util::getComb(ea, 0, true, 11);
}

int CubieCube::getD4Comb() {
    return Util::getComb(ea, 4, true, 11);
}

int CubieCube::getCPerm() {
    return Util::get8Perm(ca, false);
}

void CubieCube::setCPerm(int idx) {
    Util::set8Perm(ca, idx, false);
}

//changed
int CubieCube::getCPermSym() {
    /*
    if (EPermR2S != null) {
        int idx = EPermR2S[getCPerm()];
        return idx ^ e2c[idx & 0xf];
    }
    */
    CubieCube temps;
    for (int k = 0; k < 16; k++) {
        CornConjugate((*this), SymInv[k], temps);
        int idx = std::lower_bound(EPermS2R, EPermS2R + 2768, (int) temps.getCPerm()) - EPermS2R;
        if (idx < 2768)
            if (EPermS2R[idx] == (int) temps.getCPerm()) {
                return idx << 4 | k;
            }
    }
    assert(0);
    return 0;
}

int CubieCube::getEPerm() {
    return Util::get8Perm(ea, true);
}

void CubieCube::setEPerm(int idx) {
    Util::set8Perm(ea, idx, true);
}

//changed
int CubieCube::getEPermSym() {
    /*
    if (EPermR2S != null) {
        return EPermR2S[getEPerm()];
    }
    */
    CubieCube temps;
    for (int k = 0; k < 16; k++) {
        EdgeConjugate((*this), SymInv[k], temps);
        int idx = std::lower_bound(EPermS2R, EPermS2R + 2768, (int) temps.getEPerm()) - EPermS2R;
        if (idx < 2768)
            if (EPermS2R[idx] == (int)temps.getEPerm()) {
                return idx << 4 | k;
            }
    }
    assert(0);
    return 0;
}

int CubieCube::getMPerm() {
    return Util::getComb(ea, 8, true, 11) >> 9;
}

void CubieCube::setMPerm(int idx) {
    Util::setComb(ea, idx << 9, 8, true, 11);
}

int CubieCube::getCComb() {
    return 69 - (Util::getComb(ca, 0, false, 7) & 0x1ff);
}

void CubieCube::setCComb(int idx) {
    Util::setComb(ca, 69 - idx, 0, false, 7);
}

int CubieCube::verify() {
    int sum = 0;
    int edgeMask = 0;
    for (int e = 0; e < 12; e++) {
        edgeMask |= 1 << (ea[e] >> 1);
        sum ^= ea[e] & 1;
    }
    if (edgeMask != 0xfff) {
        return -2;// missing edges
    }
    if (sum != 0) {
        return -3;
    }
    int cornMask = 0;
    sum = 0;
    for (int c = 0; c < 8; c++) {
        cornMask |= 1 << (ca[c] & 7);
        sum += ca[c] >> 3;
    }
    if (cornMask != 0xff) {
        return -4;// missing corners
    }
    if (sum % 3 != 0) {
        return -5;// twisted corner
    }
    if ((Util::getNParity(Util::getNPerm(ea, 12, true), 12) ^ Util::getNParity(getCPerm(), 8)) != 0) {
        return -6;// parity error
    }
    return 0;// cube ok
}

ll CubieCube::selfSymmetry() {
    CubieCube c((*this));
    CubieCube d;
    ll sym = 0LL;
    for (int i = 0; i < 48; i++) {
        CornConjugate(c, SymInv[i % 16], d);
        if (d.equalsCorn((*this))) {
            EdgeConjugate(c, SymInv[i % 16], d);
            if (d.equalsEdge((*this))) {
                sym |= 1LL << i;
            }
        }
        if (i % 16 == 15) {
            c.URFConjugate();
        }
    }
    c.invCubieCube();
    for (int i = 0; i < 48; i++) {
        CornConjugate(c, SymInv[i % 16], d);
        if (d.equalsCorn((*this))) {
            EdgeConjugate(c, SymInv[i % 16], d);
            if (d.equalsEdge((*this))) {
                sym |= 1LL << 48;
                break;
            }
        }
        if (i % 16 == 15) {
            c.URFConjugate();
        }
    }
    return sym;
}

void CubieCube::setUDSliceFlip(int idx) {
    setFlip(idx & 0x7ff);
    setUDSlice(idx >> 11);
}

int CubieCube::getUDSliceFlip() {
    return (getUDSlice() & 0x1ff) << 11 | getFlip();
}

int CubieCube::getUDSliceFlipSym() {
    int flip = getFlipSym();
    int fsym = flip & 0x7;
    flip >>= 3;
    int udslice = getUDSlice() & 0x1ff;
    int udsliceflip = FlipSlice2UDSliceFlip[flip * 495 + CoordCube::UDSliceConj[udslice][fsym]];
    return udsliceflip & 0xfffffff0 | SymMult[udsliceflip & 0xf][fsym << 1];
}

// ********************************************* Initialization functions *********************************************

void CubieCube::initMove() {
    moveCube[0].construct(15120, 0, 119750400, 0);
    moveCube[3].construct(21021, 1494, 323403417, 0);
    moveCube[6].construct(8064, 1236, 29441808, 550);
    moveCube[9].construct(9, 0, 5880, 0);
    moveCube[12].construct(1230, 412, 2949660, 0);
    moveCube[15].construct(224, 137, 328552, 137);
    for (int a = 0; a < 18; a += 3) {
        for (int p = 0; p < 2; p++) {
            EdgeMult(moveCube[a + p], moveCube[a], moveCube[a + p + 1]);
            CornMult(moveCube[a + p], moveCube[a], moveCube[a + p + 1]);
        }
    }
}

//changed
char *CubieCube::toString() {
    static char sb[101];
    memset(sb, 0, sizeof(sb));
    for (int i = 0; i < 8; i++) {
        sprintf(sb, "|%d %d", ca[i] & 7, ca[i] >> 3);
    }
    sprintf(sb, "\n");
    for (int i = 0; i < 12; i++) {
        sprintf(sb, "|%d %d", ea[i] >> 1, ea[i] & 1);
    }
    return sb;
}

void CubieCube::initSym() {
    CubieCube c;
    CubieCube d;
    CubieCube t;

    CubieCube f2(28783, 0, 259268407, 0);
    CubieCube u4(15138, 0, 119765538, 7);
    CubieCube lr2(5167, 0, 83473207, 0);
    for (int i = 0; i < 8; i++) {
        lr2.ca[i] |= 3 << 3;
    }

    for (int i = 0; i < 16; i++) {
        CubeSym[i].copy(c);
        CornMult(c, u4, d);
        EdgeMult(c, u4, d);
        t = d;  d = c;  c = t;
        if (i % 4 == 3) {
            CornMult(c, lr2, d);
            EdgeMult(c, lr2, d);
            t = d;  d = c;  c = t;
        }
        if (i % 8 == 7) {
            CornMult(c, f2, d);
            EdgeMult(c, f2, d);
            t = d;  d = c;  c = t;
        }
    }
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
            CornMult(CubeSym[i], CubeSym[j], c);
            for (int k = 0; k < 16; k++) {
                if (CubeSym[k].equalsCorn(c)) {
                    SymMult[i][j] = k;
                    if (k == 0) {
                        SymInv[i] = j;
                    }
                    break;
                }
            }
        }
    }
    for (int j = 0; j < 18; j++) {
        for (int s = 0; s < 16; s++) {
            CornConjugate(moveCube[j], SymInv[s], c);
            for (int m = 0; m < 18; m++) {
                if (c.equalsCorn(moveCube[m])) {
                    SymMove[s][j] = m;
                    break;
                }
            }
        }
    }
    for (int s = 0; s < 16; s++) {
        for (int j = 0; j < 10; j++) {
            SymMoveUD[s][j] = Util::std2ud[SymMove[s][Util::ud2std[j]]];
        }
        for (int j = 0; j < 16; j++) {
            SymMultInv[j][s] = SymMult[j][SymInv[s]];
        }
    }
    for (int s = 0; s < 8; s++) {
        for (int j = 0; j < 8; j++) {
            Sym8Mult[s << 3 | j] = SymMult[j << 1][s << 1] >> 1;
            Sym8MultInv[j << 3 | s] = SymMult[j << 1][SymInv[s << 1]]>>1;
        }
        for (int j = 0; j < 18; j++) {
            Sym8Move[j << 3 | s] = SymMove[s << 1][j];
        }
    }
    for (int i = 0; i < 18; i++) {
        moveCubeSym[i] = moveCube[i].selfSymmetry();
    }
    for (int i = 0; i < 18; i++) {
        int j = i;
        for (int s = 0; s < 48; s++) {
            if (SymMove[s % 16][j] < i) {
                firstMoveSym[s] |= 1 << i;
            }
            if (s % 16 == 15) {
                j = urfMove[2][j];
            }
        }
    }
}


void CubieCube::initFlipSym2Raw() {
    CubieCube c;
    CubieCube d;
    memset(FlipR2S, 0, sizeof(FlipR2S));
    int count = 0;
    for (int i = 0; i < 2048; i++) {
        if (FlipR2S[i] != 0) {
            continue;
        }
        c.setFlip(i);
        for (int s = 0; s < 16; s += 2) {
            EdgeConjugate(c, s, d);
            int idx = d.getFlip();
            if (idx == i) {
                SymStateFlip[count] |= 1 << (s >> 1);
            }
            FlipR2S[idx] = (int) (count << 3 | s >> 1);
            if (Search::USE_TWIST_FLIP_PRUN) {
                FlipS2RF[count << 3 | s >> 1] = (int) idx;
            }
        }
        FlipS2R[count++] = (int) i;
    }
    assert(count == 336);
}

void CubieCube::initTwistSym2Raw() {
	CubieCube c;
	CubieCube d;
	int count = 0;
	memset(TwistR2S, 0, sizeof(TwistR2S));
	for (int i = 0; i < 2187; i++) {
		if (TwistR2S[i] != 0) {
			continue;
		}
		c.setTwist(i);
		for (int s = 0; s < 16; s += 2) {
			CornConjugate(c, s, d);
			int idx = d.getTwist();
			if (idx == i) {
				SymStateTwist[count] |= 1 << (s >> 1);
			}
			TwistR2S[idx] = (int) (count << 3 | s >> 1);
			if (Search::EXTRA_PRUN_LEVEL > 0) {
				TwistS2RF[count << 3 | s >> 1] = (int) idx;
			}
		}
		TwistS2R[count++] = (int) i;
	}
	assert(count == 324);
}

void CubieCube::initPermSym2Raw() {
    CubieCube c;
    CubieCube d;
    int count = 0;
    memset(EPermR2S, 0, sizeof(EPermR2S));
    for (int i = 0; i < 40320; i++) {
        if (EPermR2S[i] != 0) {
            continue;
        }
        c.setEPerm(i);
        for (int s = 0; s < 16; s++) {
            EdgeConjugate(c, s, d);
            int idx = d.getEPerm();
            if (idx == i) {
                SymStatePerm[count] |= 1 << s;
            }
            int a = d.getU4Comb();
            int b = d.getD4Comb() >> 9;
            int m = 494 - (a & 0x1ff) + (a >> 9) * 70 + b * 1680;
            MtoEPerm[m] = EPermR2S[idx] = (int) (count << 4 | s);
            if (s == 0) {
                Perm2Comb[count] = (int) (494 - (a & 0x1ff));
            }
        }
        EPermS2R[count++] = (int) i;
    }
    assert(count == 2768);
}

 void CubieCube::initUDSliceFlipSym2Raw() {
    CubieCube c;
    CubieCube d;
    int occ[2048 * 495 >> 5];
    memset(occ, 0, sizeof(occ));
    int count = 0;
    for (int i = 0; i < 2048 * 495; i++) {
        if ((occ[i >> 5] & 1 << (i & 0x1f)) != 0) {
            continue;
        }
        c.setUDSliceFlip(i);
        for (int s = 0; s < 16; s++) {
            EdgeConjugate(c, s, d);
            int idx = d.getUDSliceFlip();
            if (idx == i) {
                SymStateUDSliceFlip[count] |= 1 << s;
            }
            occ[idx >> 5] |= 1 << (idx & 0x1f);
            int fidx = std::lower_bound(FlipS2R, FlipS2R + 336, (int) (idx & 0x7ff)) - FlipS2R;
            if (fidx < 336)
                if (FlipS2R[fidx] == (int) (idx & 0x7ff)) {
                    FlipSlice2UDSliceFlip[fidx * CoordCube::N_SLICE + (idx >> 11)] = count << 4 | s;
                }
        }
        UDSliceFlipS2R[count++] = i;
    }
    assert(count == 64430);
}
