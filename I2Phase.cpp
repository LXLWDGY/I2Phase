#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include "I2Phase.h"
using namespace std;

char s[101];
Search S;

int main() {
	scanf("%s", s);
	int last = 21;
	S.Solution(s, 21, 100, 0, 0x4);
	while (S.solution == NULL) S.next(100, 0, 0x4);
	do {
        last = S.sol;
        printf("%s\n", S.solution);
        S.next(100, 0, 0x4);
	} while (S.sol != last);
    return 0;
}
