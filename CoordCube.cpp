#include "I2Phase.h"
#include <algorithm>
#include <cstring>

int CoordCube::UDSliceMove[N_SLICE][N_MOVES];
int CoordCube::TwistMove[N_TWIST_SYM][N_MOVES];
int CoordCube::FlipMove[N_FLIP_SYM][N_MOVES];
int CoordCube::UDSliceConj[N_SLICE][8];
int CoordCube::UDSliceTwistPrun[N_SLICE * N_TWIST_SYM / 8 + 1];
int CoordCube::UDSliceFlipPrun[N_SLICE * N_FLIP_SYM / 8];
int CoordCube::TwistFlipPrun[N_FLIP * N_TWIST_SYM / 8];

int CoordCube::CPermMove[N_PERM_SYM][N_MOVES];
int CoordCube::EPermMove[N_PERM_SYM][N_MOVES2];
int CoordCube::MPermMove[N_MPERM][N_MOVES2];
int CoordCube::MPermConj[N_MPERM][16];
int CoordCube::CCombMove[N_COMB][N_MOVES];
int CoordCube::CCombConj[N_COMB][16];
int CoordCube::MCPermPrun[N_MPERM * N_PERM_SYM / 8];
int CoordCube::MEPermPrun[N_MPERM * N_PERM_SYM / 8];
int CoordCube::EPermCCombPrun[N_COMB * N_PERM_SYM / 8];

void CoordCube::init() {
	CubieCube::initPermSym2Raw();

	initCPermMove();
	initEPermMove();
	initMPermMoveConj();
	initCombMoveConj();

	initMEPermPrun();
	initMCPermPrun();
	initPermCombPrun();

	CubieCube::initFlipSym2Raw();
	CubieCube::initTwistSym2Raw();
	initFlipMove();
	initTwistMove();
	initUDSliceMoveConj();

	if (Search::USE_TWIST_FLIP_PRUN) {
		initTwistFlipPrun();
	}
	initSliceTwistPrun();
	initSliceFlipPrun();
}

void CoordCube::setPruning(int *table, int index, int value) {
	table[index >> 3] ^= (0xf ^ value) << ((index & 7) << 2);
}

int CoordCube::getPruning(int *table, int index) {
	return table[index >> 3] >> ((index & 7) << 2) & 0xf;
}

void CoordCube::initUDSliceMoveConj() {
	CubieCube c;
	CubieCube d;
	for (int i = 0; i < N_SLICE; i++) {
		c.setUDSlice(i);
		for (int j = 0; j < N_MOVES; j += 3) {
			CubieCube::EdgeMult(c, CubieCube::moveCube[j], d);
			UDSliceMove[i][j] = (int) d.getUDSlice();
		}
		for (int j = 0; j < 16; j += 2) {
			CubieCube::EdgeConjugate(c, CubieCube::SymInv[j], d);
			UDSliceConj[i][j >> 1] = (int) (d.getUDSlice() & 0x1ff);
		}
	}
	for (int i = 0; i < N_SLICE; i++) {
		for (int j = 0; j < N_MOVES; j += 3) {
			int udslice = UDSliceMove[i][j];
			for (int k = 1; k < 3; k++) {
				int cx = UDSliceMove[udslice & 0x1ff][j];
				udslice = Util::permMult[udslice >> 9][cx >> 9] << 9 | cx & 0x1ff;
				UDSliceMove[i][j + k] = (int) udslice;
			}
		}
	}
}

void CoordCube::initFlipMove() {
	CubieCube c;
	CubieCube d;
	for (int i = 0; i < N_FLIP_SYM; i++) {
		c.setFlip(CubieCube::FlipS2R[i]);
		for (int j = 0; j < N_MOVES; j++) {
			CubieCube::EdgeMult(c, CubieCube::moveCube[j], d);
			FlipMove[i][j] = (int) d.getFlipSym();
		}
	}
}

void CoordCube::initTwistMove() {
	CubieCube c;
	CubieCube d;
	for (int i = 0; i < N_TWIST_SYM; i++) {
		c.setTwist(CubieCube::TwistS2R[i]);
		for (int j = 0; j < N_MOVES; j++) {
			CubieCube::CornMult(c, CubieCube::moveCube[j], d);
			TwistMove[i][j] = (int) d.getTwistSym();
		}
	}
}

void CoordCube::initCPermMove() {
	CubieCube c;
	CubieCube d;
	for (int i = 0; i < N_PERM_SYM; i++) {
		c.setCPerm(CubieCube::EPermS2R[i]);
		for (int j = 0; j < N_MOVES; j++) {
			CubieCube::CornMult(c, CubieCube::moveCube[j], d);
			CPermMove[i][j] = (int) d.getCPermSym();
		}
	}
}

void CoordCube::initEPermMove() {
	CubieCube c;
	CubieCube d;
	for (int i = 0; i < N_PERM_SYM; i++) {
		c.setEPerm(CubieCube::EPermS2R[i]);
		for (int j = 0; j < N_MOVES2; j++) {
			CubieCube::EdgeMult(c, CubieCube::moveCube[Util::ud2std[j]], d);
			EPermMove[i][j] = (int) d.getEPermSym();
		}
	}
}

void CoordCube::initMPermMoveConj() {
	CubieCube c;
	CubieCube d;
	for (int i = 0; i < N_MPERM; i++) {
		c.setMPerm(i);
		for (int j = 0; j < N_MOVES2; j++) {
			CubieCube::EdgeMult(c, CubieCube::moveCube[Util::ud2std[j]], d);
			MPermMove[i][j] = (int) d.getMPerm();
		}
		for (int j = 0; j < 16; j++) {
			CubieCube::EdgeConjugate(c, CubieCube::SymInv[j], d);
			MPermConj[i][j] = (int) d.getMPerm();
		}
	}
}

void CoordCube::initCombMoveConj() {
	CubieCube c;
	CubieCube d;
	for (int i = 0; i < N_COMB; i++) {
		c.setCComb(i);
		for (int j = 0; j < N_MOVES; j++) {
			CubieCube::CornMult(c, CubieCube::moveCube[j], d);
			CCombMove[i][j] = (int) d.getCComb();
		}
		for (int j = 0; j < 16; j++) {
			CubieCube::CornConjugate(c, CubieCube::SymInv[j], d);
			CCombConj[i][j] = (int) d.getCComb();
		}
	}
}

void CoordCube::initTwistFlipPrun() {
	int depth = 0;
	int done = 1;
	bool inv;
	int select;
	int check;
	const int N_SIZE = N_FLIP * N_TWIST_SYM;
	for (int i = 0; i < N_SIZE / 8; i++) {
		TwistFlipPrun[i] = -1;
	}
	setPruning(TwistFlipPrun, 0, 0);

	while (done < N_SIZE) {
		inv = depth > 6;
		select = inv ? 0xf : depth;
		check = inv ? depth : 0xf;
		depth++;
		int val = 0;
		for (int i = 0; i < N_SIZE; i++, val >>= 4) {
			if ((i & 7) == 0) {
				val = TwistFlipPrun[i >> 3];
				if (!inv && val == -1) {
					i += 7;
					continue;
				}
			}
			if ((val & 0xf) != select) {
				continue;
			}
			int twist = i >> 11;
			int flip = CubieCube::FlipR2S[i & 0x7ff];
			int fsym = flip & 7;
			flip >>= 3;
			for (int m = 0; m < N_MOVES; m++) {
				int twistx = TwistMove[twist][m];
				int tsymx = twistx & 7;
				twistx >>= 3;
				int flipx = FlipMove[flip][CubieCube::Sym8Move[m << 3 | fsym]];
				int fsymx = CubieCube::Sym8MultInv[CubieCube::Sym8Mult[flipx & 7 | fsym << 3] << 3 | tsymx];
				flipx >>= 3;
				int idx = twistx << 11 | CubieCube::FlipS2RF[flipx << 3 | fsymx];
				if (getPruning(TwistFlipPrun, idx) != check) {
					continue;
				}
				done++;
				if (inv) {
					setPruning(TwistFlipPrun, i, depth);
					break;
				}
				setPruning(TwistFlipPrun, idx, depth);
				int sym = CubieCube::SymStateTwist[twistx];
				if (sym == 1) {
					continue;
				}
				for (int k = 0; k < 8; k++) {
					if ((sym & 1 << k) == 0) {
						continue;
					}
					int idxx = twistx << 11 | CubieCube::FlipS2RF[flipx << 3 | CubieCube::Sym8MultInv[fsymx << 3 | k]];
					if (getPruning(TwistFlipPrun, idxx) == 0xf) {
						setPruning(TwistFlipPrun, idxx, depth);
						done++;
					}
				}
			}
		}
		// System.out.println(String.format("%2d%10d", depth, done));
	}
}

void CoordCube::initRawSymPrun(int *PrunTable, const int INV_DEPTH, int (&RawMove)[N_SLICE][N_MOVES], int (&RawConj)[N_SLICE][8], int (&SymMove)[N_TWIST_SYM][N_MOVES], int *SymState, const int PrunFlag, const int N_RAW, const int N_SYM, const int NN_MOVES) {

	const int SYM_SHIFT = PrunFlag & 0xf;
	const bool SymSwitch = ((PrunFlag >> 4) & 1) == 1;
	const bool MoveMapSym = ((PrunFlag >> 5) & 1) == 1;
	const bool MoveMapRaw = ((PrunFlag >> 6) & 1) == 1;

	const int SYM_MASK = (1 << SYM_SHIFT) - 1;
	const int N_SIZE = N_RAW * N_SYM;
	const int N_MOVES = MoveMapRaw ? 10 : NN_MOVES;

	for (int i = 0; i < (N_RAW * N_SYM + 7) / 8; i++) {
		PrunTable[i] = -1;
	}
	setPruning(PrunTable, 0, 0);

	int depth = 0;
	int done = 1;

	while (done < N_SIZE) {
		bool inv = depth > INV_DEPTH;
		int select = inv ? 0xf : depth;
		int check = inv ? depth : 0xf;
		depth++;
		int val = 0;
		for (int i = 0; i < N_SIZE; i++, val >>= 4) {
			if ((i & 7) == 0) {
				val = PrunTable[i >> 3];
				if (!inv && val == -1) {
					i += 7;
					continue;
				}
			}
			if ((val & 0xf) != select) {
				continue;
			}
			int raw = i % N_RAW;
			int sym = i / N_RAW;
			for (int m = 0; m < N_MOVES; m++) {
				int symx = SymMove[sym][MoveMapSym ? Util::ud2std[m] : m];
				int rawx = RawConj[RawMove[raw][MoveMapRaw ? Util::ud2std[m] : m] & 0x1ff][symx & SYM_MASK];
				symx >>= SYM_SHIFT;
				int idx = symx * N_RAW + rawx;
				if (getPruning(PrunTable, idx) != check) {
					continue;
				}
				done++;
				if (inv) {
					setPruning(PrunTable, i, depth);
					break;
				}
				setPruning(PrunTable, idx, depth);
				for (int j = 1, symState = SymState[symx]; (symState >>= 1) != 0; j++) {
					if ((symState & 1) != 1) {
						continue;
					}
					int idxx = symx * N_RAW + RawConj[rawx][j ^ (SymSwitch ? CubieCube::e2c[j] : 0)];
					if (getPruning(PrunTable, idxx) == 0xf) {
						setPruning(PrunTable, idxx, depth);
						done++;
					}
				}
			}
		}
		// System.out.println(String.format("%2d%10d", depth, done));
	}
}

void CoordCube::initRawSymPrun(int *PrunTable, const int INV_DEPTH, int (&RawMove)[N_SLICE][N_MOVES], int (&RawConj)[N_SLICE][8], int (&SymMove)[N_FLIP_SYM][N_MOVES], int *SymState, const int PrunFlag, const int N_RAW, const int N_SYM, const int NN_MOVES) {

	const int SYM_SHIFT = PrunFlag & 0xf;
	const bool SymSwitch = ((PrunFlag >> 4) & 1) == 1;
	const bool MoveMapSym = ((PrunFlag >> 5) & 1) == 1;
	const bool MoveMapRaw = ((PrunFlag >> 6) & 1) == 1;

	const int SYM_MASK = (1 << SYM_SHIFT) - 1;
	const int N_SIZE = N_RAW * N_SYM;
	const int N_MOVES = MoveMapRaw ? 10 : NN_MOVES;

	for (int i = 0; i < (N_RAW * N_SYM + 7) / 8; i++) {
		PrunTable[i] = -1;
	}
	setPruning(PrunTable, 0, 0);

	int depth = 0;
	int done = 1;

	while (done < N_SIZE) {
		bool inv = depth > INV_DEPTH;
		int select = inv ? 0xf : depth;
		int check = inv ? depth : 0xf;
		depth++;
		int val = 0;
		for (int i = 0; i < N_SIZE; i++, val >>= 4) {
			if ((i & 7) == 0) {
				val = PrunTable[i >> 3];
				if (!inv && val == -1) {
					i += 7;
					continue;
				}
			}
			if ((val & 0xf) != select) {
				continue;
			}
			int raw = i % N_RAW;
			int sym = i / N_RAW;
			for (int m = 0; m < N_MOVES; m++) {
				int symx = SymMove[sym][MoveMapSym ? Util::ud2std[m] : m];
				int rawx = RawConj[RawMove[raw][MoveMapRaw ? Util::ud2std[m] : m] & 0x1ff][symx & SYM_MASK];
				symx >>= SYM_SHIFT;
				int idx = symx * N_RAW + rawx;
				if (getPruning(PrunTable, idx) != check) {
					continue;
				}
				done++;
				if (inv) {
					setPruning(PrunTable, i, depth);
					break;
				}
				setPruning(PrunTable, idx, depth);
				for (int j = 1, symState = SymState[symx]; (symState >>= 1) != 0; j++) {
					if ((symState & 1) != 1) {
						continue;
					}
					int idxx = symx * N_RAW + RawConj[rawx][j ^ (SymSwitch ? CubieCube::e2c[j] : 0)];
					if (getPruning(PrunTable, idxx) == 0xf) {
						setPruning(PrunTable, idxx, depth);
						done++;
					}
				}
			}
		}
		// System.out.println(String.format("%2d%10d", depth, done));
	}
}

void CoordCube::initRawSymPrun(int *PrunTable, const int INV_DEPTH, int (&RawMove)[N_MPERM][N_MOVES2], int (&RawConj)[N_MPERM][16], int (&SymMove)[N_PERM_SYM][N_MOVES2], int *SymState, const int PrunFlag, const int N_RAW, const int N_SYM, const int NN_MOVES) {

	const int SYM_SHIFT = PrunFlag & 0xf;
	const bool SymSwitch = ((PrunFlag >> 4) & 1) == 1;
	const bool MoveMapSym = ((PrunFlag >> 5) & 1) == 1;
	const bool MoveMapRaw = ((PrunFlag >> 6) & 1) == 1;

	const int SYM_MASK = (1 << SYM_SHIFT) - 1;
	const int N_SIZE = N_RAW * N_SYM;
	const int N_MOVES = MoveMapRaw ? 10 : NN_MOVES;

	for (int i = 0; i < (N_RAW * N_SYM + 7) / 8; i++) {
		PrunTable[i] = -1;
	}
	setPruning(PrunTable, 0, 0);

	int depth = 0;
	int done = 1;

	while (done < N_SIZE) {
		bool inv = depth > INV_DEPTH;
		int select = inv ? 0xf : depth;
		int check = inv ? depth : 0xf;
		depth++;
		int val = 0;
		for (int i = 0; i < N_SIZE; i++, val >>= 4) {
			if ((i & 7) == 0) {
				val = PrunTable[i >> 3];
				if (!inv && val == -1) {
					i += 7;
					continue;
				}
			}
			if ((val & 0xf) != select) {
				continue;
			}
			int raw = i % N_RAW;
			int sym = i / N_RAW;
			for (int m = 0; m < N_MOVES; m++) {
				int symx = SymMove[sym][MoveMapSym ? Util::ud2std[m] : m];
				int rawx = RawConj[RawMove[raw][MoveMapRaw ? Util::ud2std[m] : m] & 0x1ff][symx & SYM_MASK];
				symx >>= SYM_SHIFT;
				int idx = symx * N_RAW + rawx;
				if (getPruning(PrunTable, idx) != check) {
					continue;
				}
				done++;
				if (inv) {
					setPruning(PrunTable, i, depth);
					break;
				}
				setPruning(PrunTable, idx, depth);
				for (int j = 1, symState = SymState[symx]; (symState >>= 1) != 0; j++) {
					if ((symState & 1) != 1) {
						continue;
					}
					int idxx = symx * N_RAW + RawConj[rawx][j ^ (SymSwitch ? CubieCube::e2c[j] : 0)];
					if (getPruning(PrunTable, idxx) == 0xf) {
						setPruning(PrunTable, idxx, depth);
						done++;
					}
				}
			}
		}
		// System.out.println(String.format("%2d%10d", depth, done));
	}
}

void CoordCube::initRawSymPrun(int *PrunTable, const int INV_DEPTH, int (&RawMove)[N_MPERM][N_MOVES2], int (&RawConj)[N_MPERM][16], int (&SymMove)[N_PERM_SYM][N_MOVES], int *SymState, const int PrunFlag, const int N_RAW, const int N_SYM, const int NN_MOVES) {

	const int SYM_SHIFT = PrunFlag & 0xf;
	const bool SymSwitch = ((PrunFlag >> 4) & 1) == 1;
	const bool MoveMapSym = ((PrunFlag >> 5) & 1) == 1;
	const bool MoveMapRaw = ((PrunFlag >> 6) & 1) == 1;

	const int SYM_MASK = (1 << SYM_SHIFT) - 1;
	const int N_SIZE = N_RAW * N_SYM;
	const int N_MOVES = MoveMapRaw ? 10 : NN_MOVES;

	for (int i = 0; i < (N_RAW * N_SYM + 7) / 8; i++) {
		PrunTable[i] = -1;
	}
	setPruning(PrunTable, 0, 0);

	int depth = 0;
	int done = 1;

	while (done < N_SIZE) {
		bool inv = depth > INV_DEPTH;
		int select = inv ? 0xf : depth;
		int check = inv ? depth : 0xf;
		depth++;
		int val = 0;
		for (int i = 0; i < N_SIZE; i++, val >>= 4) {
			if ((i & 7) == 0) {
				val = PrunTable[i >> 3];
				if (!inv && val == -1) {
					i += 7;
					continue;
				}
			}
			if ((val & 0xf) != select) {
				continue;
			}
			int raw = i % N_RAW;
			int sym = i / N_RAW;
			for (int m = 0; m < N_MOVES; m++) {
				int symx = SymMove[sym][MoveMapSym ? Util::ud2std[m] : m];
				int rawx = RawConj[RawMove[raw][MoveMapRaw ? Util::ud2std[m] : m] & 0x1ff][symx & SYM_MASK];
				symx >>= SYM_SHIFT;
				int idx = symx * N_RAW + rawx;
				if (getPruning(PrunTable, idx) != check) {
					continue;
				}
				done++;
				if (inv) {
					setPruning(PrunTable, i, depth);
					break;
				}
				setPruning(PrunTable, idx, depth);
				for (int j = 1, symState = SymState[symx]; (symState >>= 1) != 0; j++) {
					if ((symState & 1) != 1) {
						continue;
					}
					int idxx = symx * N_RAW + RawConj[rawx][j ^ (SymSwitch ? CubieCube::e2c[j] : 0)];
					if (getPruning(PrunTable, idxx) == 0xf) {
						setPruning(PrunTable, idxx, depth);
						done++;
					}
				}
			}
		}
		// System.out.println(String.format("%2d%10d", depth, done));
	}
}

void CoordCube::initRawSymPrun(int *PrunTable, const int INV_DEPTH, int (&RawMove)[N_COMB][N_MOVES], int (&RawConj)[N_COMB][16], int (&SymMove)[N_PERM_SYM][N_MOVES2], int *SymState, const int PrunFlag, const int N_RAW, const int N_SYM, const int NN_MOVES) {

	const int SYM_SHIFT = PrunFlag & 0xf;
	const bool SymSwitch = ((PrunFlag >> 4) & 1) == 1;
	const bool MoveMapSym = ((PrunFlag >> 5) & 1) == 1;
	const bool MoveMapRaw = ((PrunFlag >> 6) & 1) == 1;

	const int SYM_MASK = (1 << SYM_SHIFT) - 1;
	const int N_SIZE = N_RAW * N_SYM;
	const int N_MOVES = MoveMapRaw ? 10 : NN_MOVES;

	for (int i = 0; i < (N_RAW * N_SYM + 7) / 8; i++) {
		PrunTable[i] = -1;
	}
	setPruning(PrunTable, 0, 0);

	int depth = 0;
	int done = 1;

	while (done < N_SIZE) {
		bool inv = depth > INV_DEPTH;
		int select = inv ? 0xf : depth;
		int check = inv ? depth : 0xf;
		depth++;
		int val = 0;
		for (int i = 0; i < N_SIZE; i++, val >>= 4) {
			if ((i & 7) == 0) {
				val = PrunTable[i >> 3];
				if (!inv && val == -1) {
					i += 7;
					continue;
				}
			}
			if ((val & 0xf) != select) {
				continue;
			}
			int raw = i % N_RAW;
			int sym = i / N_RAW;
			for (int m = 0; m < N_MOVES; m++) {
				int symx = SymMove[sym][MoveMapSym ? Util::ud2std[m] : m];
				int rawx = RawConj[RawMove[raw][MoveMapRaw ? Util::ud2std[m] : m] & 0x1ff][symx & SYM_MASK];
				symx >>= SYM_SHIFT;
				int idx = symx * N_RAW + rawx;
				if (getPruning(PrunTable, idx) != check) {
					continue;
				}
				done++;
				if (inv) {
					setPruning(PrunTable, i, depth);
					break;
				}
				setPruning(PrunTable, idx, depth);
				for (int j = 1, symState = SymState[symx]; (symState >>= 1) != 0; j++) {
					if ((symState & 1) != 1) {
						continue;
					}
					int idxx = symx * N_RAW + RawConj[rawx][j ^ (SymSwitch ? CubieCube::e2c[j] : 0)];
					if (getPruning(PrunTable, idxx) == 0xf) {
						setPruning(PrunTable, idxx, depth);
						done++;
					}
				}
			}
		}
		// System.out.println(String.format("%2d%10d", depth, done));
	}
}

void CoordCube::initSliceTwistPrun() {
	initRawSymPrun(UDSliceTwistPrun, 6,	UDSliceMove, UDSliceConj, TwistMove, CubieCube::SymStateTwist, 0x3, N_SLICE, N_TWIST_SYM, N_MOVES);
}

void CoordCube::initSliceFlipPrun() {
	initRawSymPrun(UDSliceFlipPrun, 6, UDSliceMove, UDSliceConj, FlipMove, CubieCube::SymStateFlip, 0x3, N_SLICE, N_FLIP_SYM, N_MOVES);
}

void CoordCube::initMEPermPrun() {
	initRawSymPrun(MEPermPrun, 7, MPermMove, MPermConj, EPermMove, CubieCube::SymStatePerm, 0x4, N_MPERM, N_PERM_SYM, N_MOVES2);
}

void CoordCube::initMCPermPrun() {
	initRawSymPrun(MCPermPrun, 10, MPermMove, MPermConj, CPermMove, CubieCube::SymStatePerm, 0x34, N_MPERM, N_PERM_SYM, N_MOVES2);
}

void CoordCube::initPermCombPrun() {
	initRawSymPrun(EPermCCombPrun, 8, CCombMove, CCombConj, EPermMove, CubieCube::SymStatePerm, 0x44, N_COMB, N_PERM_SYM, N_MOVES);
}

CoordCube::CoordCube() {
    twist = 0;
    tsym = 0;
    flip = 0;
    fsym = 0;
    slice = 0;
    prun = 0;
}


void CoordCube::set(CoordCube node) {
	this->twist = node.twist;
	this->tsym = node.tsym;
	this->flip = node.flip;
	this->fsym = node.fsym;
	this->slice = node.slice;
	this->prun = node.prun;
}

void CoordCube::calcPruning(bool isPhase1) {
	prun = std::max(
			   std::max(
				   getPruning(UDSliceTwistPrun,
							  twist * N_SLICE + UDSliceConj[slice & 0x1ff][tsym]),
				   getPruning(UDSliceFlipPrun,
							  flip * N_SLICE + UDSliceConj[slice & 0x1ff][fsym])),
			   Search::USE_TWIST_FLIP_PRUN ? getPruning(TwistFlipPrun,
					   twist << 11 | CubieCube::FlipS2RF[flip << 3 | CubieCube::Sym8MultInv[fsym << 3 | tsym]]) : 0);
}

void CoordCube::set(CubieCube cc) {
    twist = cc.getTwistSym();
    flip = cc.getFlipSym();
    slice = cc.getUDSlice();
    tsym = twist & 7;
    twist = twist >> 3;
    fsym = flip & 7;
    flip = flip >> 3;
}

int CoordCube::doMovePrun(CoordCube cc, int m, bool isPhase1) {
	slice = UDSliceMove[cc.slice & 0x1ff][m] & 0x1ff;

	flip = FlipMove[cc.flip][CubieCube::Sym8Move[m << 3 | cc.fsym]];
	fsym = CubieCube::Sym8Mult[flip & 7 | cc.fsym << 3];
	flip >>= 3;

	twist = TwistMove[cc.twist][CubieCube::Sym8Move[m << 3 | cc.tsym]];
	tsym = CubieCube::Sym8Mult[twist & 7 | cc.tsym << 3];
	twist >>= 3;

	prun = std::max(
			   std::max(
				   getPruning(UDSliceTwistPrun,
							  twist * N_SLICE + UDSliceConj[slice][tsym]),
				   getPruning(UDSliceFlipPrun,
							  flip * N_SLICE + UDSliceConj[slice][fsym])),
			   Search::USE_TWIST_FLIP_PRUN ? getPruning(TwistFlipPrun,
					   twist << 11 | CubieCube::FlipS2RF[flip << 3 | CubieCube::Sym8MultInv[fsym << 3 | tsym]]) : 0);
	return prun;
}

