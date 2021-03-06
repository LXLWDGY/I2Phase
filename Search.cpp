#include "I2Phase.h"
#include <cmath>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>

bool Search::inited = false;
Search::Search() {
    (*this).isRec = false;
}
char *Search::Solution(char *facelets, int maxDepth, ll probeMax, ll probeMin, int verbose) {
	int check = verify(facelets);
	if (check != 0) {
        static char s[10] = "Error";
		return s;
	}
	(*this).sol = maxDepth + 1;
	(*this).probe = 0;
	(*this).probeMax = probeMax;
	(*this).probeMin = std::min(probeMin, probeMax);
	(*this).verbose = verbose;
	(*this).solution = NULL;
	(*this).isRec = false;

	init();

	initSearch();

	return search();
}

void Search::initSearch() {
	conjMask = (TRY_INVERSE ? 0 : 0x38) | (TRY_THREE_AXES ? 0 : 0x36);
	CubieCube pc;
	selfSym = cc.selfSymmetry();
	if ((selfSym >> 48) != 0) {
		conjMask |= 0x38;
	}
	if ((selfSym >> 16 & 0xffff) != 0) {
		conjMask |= 0x12;
	}
	if ((selfSym >> 32 & 0xffff) != 0) {
		conjMask |= 0x24;
	}
	preIdxMax = conjMask > 7 ? 1 : PRE_IDX_MAX;
	for (int i = 0; i < 6; i++) {
		node0[i][0].set(cc);
		corn0[i][0] = cc.getCPermSym();
		ud8e0[i][0] = cc.getU4Comb() << 16 | cc.getD4Comb();
		if ((conjMask & 1 << i) == 0) {
			for (int j = 1; j < preIdxMax; j++) {
				CubieCube::CornMult(CubieCube::moveCube[CubieCube::preMove[j]], cc, pc);
				CubieCube::EdgeMult(CubieCube::moveCube[CubieCube::preMove[j]], cc, pc);
				node0[i][j].set(pc);
				corn0[i][j] = pc.getCPermSym();
				ud8e0[i][j] = pc.getU4Comb() << 16 | pc.getD4Comb();
			}
		}
		cc.URFConjugate();
		if (i % 3 == 2) {
			cc.invCubieCube();
		}
	}
	selfSym = selfSym & 0xffffffffffffLL;
}

char *Search::next(ll probeMax, ll probeMin, int verbose) {
	(*this).probe = 0;
	(*this).probeMax = probeMax;
	(*this).probeMin = std::min(probeMin, probeMax);
	(*this).solution = NULL;
	(*this).isRec = ((*this).verbose & OPTIMAL_SOLUTION) == (verbose & OPTIMAL_SOLUTION);
	(*this).verbose = verbose;
	return search();
}

bool Search::isInited() {
	return inited;
}

ll Search::numberOfProbes() {
	return probe;
}

int Search::length() {
	return sol;
}

void Search::init() {
	if (inited) {
		return;
	}
	CubieCube::initMove();
	CubieCube::initSym();

    CoordCube::init();

	inited = true;
}

int Search::verify(char *facelets) {
	int count = 0x000000;
    char center[6] = {
        facelets[Util::U5],
        facelets[Util::R5],
        facelets[Util::F5],
        facelets[Util::D5],
        facelets[Util::L5],
        facelets[Util::B5]
    };
    for (int i = 0; i < 54; i++) {
        f[i] = (char) (strchr(center, facelets[i]) - center);
        if (f[i] == -1) {
            return -1;
        }
        count += 1 << (f[i] << 2);
    }
	if (count != 0x999999) {
		return -1;
	}
	Util::toCubieCube(f, cc);
	return cc.verify();
}

char *Search::search() {
	for (length1 = isRec ? length1 : 0; length1 < sol; length1++) {
		maxDep2 = std::min(13, sol - length1);
		for (urfIdx = isRec ? urfIdx : 0; urfIdx < 6; urfIdx++) {
			if ((conjMask & 1 << urfIdx) != 0) {
				continue;
			}
			for (preIdx = isRec ? preIdx : 0; preIdx < preIdxMax; preIdx++) {
				if (preIdx != 0 && preIdx % 2 == 0) {
					continue;
				}
				node0[urfIdx][preIdx].calcPruning(true);
				int ssym = (int) (0xffff & selfSym);
				if (preIdx != 0) {
					ssym &= CubieCube::moveCubeSym[CubieCube::preMove[preIdx]];
				}
				depth1 = length1 - (preIdx == 0 ? 0 : 1);
				if (node0[urfIdx][preIdx].prun <= depth1) {
						if (phase1(node0[urfIdx][preIdx], ssym, depth1, -1) == 0) {
                            return solution;
                        }
				}
			}
		}
	}
	return solution;
}

int Search::phase1(CoordCube node, ll ssym, int maxl, int lm) {
	if (node.prun == 0 && maxl < 5) {
		if (maxl == 0) {
			int ret = initPhase2();
			if (ret == 0 || preIdx == 0) {
				return ret;
			}
			preIdx++;
			ret = std::min(initPhase2(), ret);
			preIdx--;
			return ret;
		} else {
			return 1;
		}
	}

	int skipMoves = 0;
	int i = 1;
	for (ll s = ssym; (s >>= 1) != 0; i++) {
		if ((s & 1) == 1) {
			skipMoves |= CubieCube::firstMoveSym[i];
		}
	}

	for (int axis = 0; axis < 18; axis += 3) {
		if (axis == lm || axis == lm - 9
				|| (isRec && axis < move[depth1 - maxl] - 2)) {
			continue;
		}
		for (int power = 0; power < 3; power++) {
			int m = axis + power;

			if (isRec && m != move[depth1 - maxl]
					|| ssym != 1 && (skipMoves & 1 << m) != 0) {
				continue;
			}

			int prun = nodeUD[maxl].doMovePrun(node, m, true);
			if (prun > maxl) {
				break;
			} else if (prun == maxl) {
				continue;
			}

			move[depth1 - maxl] = m;
			int ret = phase1(nodeUD[maxl], ssym & CubieCube::moveCubeSym[m], maxl - 1, axis);
			if (ret == 0) {
				return 0;
			} else if (ret == 2) {
				break;
			}
		}
	}
	return 1;
}

int Search::initPhase2() {
	isRec = false;
	if (probe >= (solution == NULL ? probeMax : probeMin)) {
		return 0;
	}
	++probe;
	int cidx = corn0[urfIdx][preIdx] >> 4;
	int csym = corn0[urfIdx][preIdx] & 0xf;
	int mid = node0[urfIdx][preIdx].slice;
	for (int i = 0; i < depth1; i++) {
		int m = move[i];
		cidx = CoordCube::CPermMove[cidx][CubieCube::SymMove[csym][m]];
		csym = CubieCube::SymMult[cidx & 0xf][csym];
		cidx >>= 4;

		int cx = CoordCube::UDSliceMove[mid & 0x1ff][m];
		mid = Util::permMult[mid >> 9][cx >> 9] << 9 | cx & 0x1ff;
	}
	mid >>= 9;
	int prun = CoordCube::getPruning(CoordCube::MCPermPrun, cidx * 24 + CoordCube::MPermConj[mid][csym]);
	if (prun >= maxDep2) {
		return prun > maxDep2 ? 2 : 1;
	}

	int u4e = ud8e0[urfIdx][preIdx] >> 16;
	int d4e = ud8e0[urfIdx][preIdx] & 0xffff;
	for (int i = 0; i < depth1; i++) {
		int m = move[i];

		int cx = CoordCube::UDSliceMove[u4e & 0x1ff][m];
		u4e = Util::permMult[u4e >> 9][cx >> 9] << 9 | cx & 0x1ff;

		cx = CoordCube::UDSliceMove[d4e & 0x1ff][m];
		d4e = Util::permMult[d4e >> 9][cx >> 9] << 9 | cx & 0x1ff;
	}

	int edge = CubieCube::MtoEPerm[494 - (u4e & 0x1ff) + (u4e >> 9) * 70 + (d4e >> 9) * 1680];
	int esym = edge & 0xf;
	edge >>= 4;

	prun = std::max(prun, std::max(
						CoordCube::getPruning(CoordCube::MEPermPrun,
								edge * 24 + CoordCube::MPermConj[mid][esym]),
						CoordCube::getPruning(CoordCube::EPermCCombPrun,
								edge * 70 + CoordCube::CCombConj[CubieCube::Perm2Comb[cidx]][CubieCube::SymMultInv[esym][csym]])));

	if (prun >= maxDep2) {
		return prun > maxDep2 ? 2 : 1;
	}

	int lm = 10;
	if (depth1 >= 2 && move[depth1 - 1] / 3 % 3 == move[depth1 - 2] / 3 % 3) {
		lm = Util::std2ud[std::max(move[depth1 - 1], move[depth1 - 2]) / 3 * 3 + 1];
	} else if (depth1 >= 1) {
		lm = Util::std2ud[move[depth1 - 1] / 3 * 3 + 1];
		if (move[depth1 - 1] > Util::Fx3) {
			lm = -lm;
		}
	}

	int depth2;
	for (depth2 = maxDep2 - 1; depth2 >= prun; depth2--) {
		int ret = phase2(edge, esym, cidx, csym, mid, depth2, depth1, lm);
		if (ret < 0) {
			break;
		}
		depth2 = depth2 - ret;
		sol = depth1 + depth2;
		if (preIdx != 0) {
			assert(depth2 > 0); //If depth2 == 0, the solution is optimal. In (*this) case, we won't try preScramble to find inter solutions.
			int axisPre = Util::preMove[preIdx] / 3;
			int axisLast = move[sol - 1] / 3;
			if (axisPre == axisLast) {
				int pow = (Util::preMove[preIdx] % 3 + move[sol - 1] % 3 + 1) % 4;
				move[sol - 1] = axisPre * 3 + pow;
			} else if (depth2 > 1
					   && axisPre % 3 == axisLast % 3
					   && move[sol - 2] / 3 == axisPre) {
				int pow = (Util::preMove[preIdx] % 3 + move[sol - 2] % 3 + 1) % 4;
				move[sol - 2] = axisPre * 3 + pow;
			} else {
				move[sol++] = Util::preMove[preIdx];
			}
		}
		solution = solutionToString();
	}

	if (depth2 != maxDep2 - 1) { //At least one solution has been found.
		maxDep2 = std::min(13, sol - length1);
		return probe >= probeMin ? 0 : 1;
	} else {
		return 1;
	}
}

int Search::phase2(int eidx, int esym, int cidx, int csym, int mid, int maxl, int depth, int lm) {
	if (eidx == 0 && cidx == 0 && mid == 0) {
		return maxl;
	}
	for (int m = 0; m < 10; m++) {
		if (lm < 0 ? (m == -lm) : Util::ckmv2[lm][m]) {
			continue;
		}
		int midx = CoordCube::MPermMove[mid][m];
		int cidxx = CoordCube::CPermMove[cidx][CubieCube::SymMove[csym][Util::ud2std[m]]];
		int csymx = CubieCube::SymMult[cidxx & 0xf][csym];
		cidxx >>= 4;
		if (CoordCube::getPruning(CoordCube::MCPermPrun,
								 cidxx * 24 + CoordCube::MPermConj[midx][csymx]) >= maxl) {
			continue;
		}
		int eidxx = CoordCube::EPermMove[eidx][CubieCube::SymMoveUD[esym][m]];
		int esymx = CubieCube::SymMult[eidxx & 0xf][esym];
		eidxx >>= 4;
		if (CoordCube::getPruning(CoordCube::EPermCCombPrun,
								 eidxx * 70 + CoordCube::CCombConj[CubieCube::Perm2Comb[cidxx]][CubieCube::SymMultInv[esymx][csymx]]) >= maxl) {
			continue;
		}
		if (CoordCube::getPruning(CoordCube::MEPermPrun,
								 eidxx * 24 + CoordCube::MPermConj[midx][esymx]) >= maxl) {
			continue;
		}
		int ret = phase2(eidxx, esymx, cidxx, csymx, midx, maxl - 1, depth + 1, (lm < 0 && m + lm == -5) ? -lm : m);
		if (ret >= 0) {
			move[depth] = Util::ud2std[m];
			return ret;
		}
	}
	return -1;
}

char *Search::solutionToString() {
    static char sb[101];
    int ls = 0;
    memset(sb, 0, sizeof(sb));
	int urf = (verbose & INVERSE_SOLUTION) != 0 ? (urfIdx + 3) % 6 : urfIdx;
	if (urf < 3) {
		for (int s = 0; s < sol; s++) {
			if ((verbose & USE_SEPARATOR) != 0 && s == depth1) {
				ls += sprintf(sb + ls, ".  ");
			}
			ls += sprintf(sb + ls, "%s ", Util::move2str[CubieCube::urfMove[urf][move[s]]]);
		}
	} else {
		for (int s = sol - 1; s >= 0; s--) {
            ls += sprintf(sb + ls, "%s ", Util::move2str[CubieCube::urfMove[urf][move[s]]]);
			if ((verbose & USE_SEPARATOR) != 0 && s == depth1) {
				ls += sprintf(sb + ls, ".  ");
			}
		}
	}
	if ((verbose & APPEND_LENGTH) != 0) {
        ls += sprintf(sb + ls, "(%df)", sol);
	}
	return sb;
}

