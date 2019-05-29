#include "grid.h"
#include "solver.h"

int main() {
    SmokeSolver solver(32,64,32);
	for (int i = 0; i < 200; ++i) {
		solver.run_dt(i);
		solver.dump_density(i+1);
	}
	system("pause");
}