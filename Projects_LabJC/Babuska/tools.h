#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzinterpolationspace.h"
#include "TPZMatLaplacian.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"

#include "pzgengrid.h"
#include "pzfunction.h"

#include "pzlog.h"

#include <iostream>
#include <math.h>

#include "TPZVTKGeoMesh.h"
#include <iostream>
using namespace std;
static void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du) {

	const REAL n = 0;

	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL r = sqrt(x*x + y * y);
	const REAL t = atan2(y, x);
	const REAL sol = pow((REAL)2, 0.25 + n / 2.)*pow(r, 0.5 + n)*cos((0.5 + n)*t);
	u[0] = sol;

	du(0, 0) = pow((REAL)2, -0.75 + n / 2.)*(1 + 2 * n)*pow(pow(x, 2) + pow(y, 2), -0.75 + n / 2.)*(x*cos((0.5 + n)*atan2(y, x)) + y * sin((0.5 + n)*atan2(y, x)));
	du(1, 0) = pow((REAL)2, -0.75 + n / 2.)*(1 + 2 * n)*pow(pow(x, 2) + pow(y, 2), -0.75 + n / 2.)*(y*cos((0.5 + n)*atan2(y, x)) - x * sin((0.5 + n)*atan2(y, x)));

}

static void NeumannEsquerda(const TPZVec<REAL> &loc, TPZVec<STATE> &result) {
	REAL normal[2] = { -1,0 };

	TPZManVector<STATE> u(1);
	TPZFNMatrix<10, STATE> du(2, 1);
	SolExataSteklov(loc, u, du);

	result.Resize(1);
	result[0] = du(0, 0)*normal[0] + du(1, 0)*normal[1];
}

static void NeumannDireita(const TPZVec<REAL> &loc, TPZVec<STATE> &result) {
	REAL normal[2] = { +1,0 };

	TPZManVector<STATE> u(1);
	TPZFNMatrix<10, STATE> du(2, 1);
	SolExataSteklov(loc, u, du);

	result.Resize(1);
	result[0] = du(0, 0)*normal[0] + du(1, 0)*normal[1];
}

static void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result) {
	REAL normal[2] = { 0,+1 };

	TPZManVector<STATE> u(1);
	TPZFNMatrix<10, STATE> du(2, 1);
	SolExataSteklov(loc, u, du);

	result.Resize(1);
	result[0] = du(0, 0)*normal[0] + du(1, 0)*normal[1];
}

static void NeumannAbaixo(const TPZVec<REAL> &loc, TPZVec<STATE> &result) {
	REAL normal[2] = { 0,+1 };

	TPZManVector<STATE> u(1);
	TPZFNMatrix<10, STATE> du(2, 1);
	SolExataSteklov(loc, u, du);

	result.Resize(1);
	result[0] = du(0, 0)*normal[0] + du(1, 0)*normal[1];
}