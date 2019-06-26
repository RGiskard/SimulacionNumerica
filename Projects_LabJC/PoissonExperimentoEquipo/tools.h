#pragma once
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
/*Funcion para pruebas de elMEnto computacional
posiblemente calculo de gradiente, solucion, otros
*/
void fun(TPZCompEl *cEl)
{
	TPZVec< REAL > qsi(3,0.);
	//calculo del centro de un elemento en coordenadas master
	(cEl->Reference())->CenterPoint(6, qsi);
	TPZVec<REAL> punto(3, 0.);
	//Transformamos a coordenadas cartesianas
	(cEl->Reference())->X(qsi, punto);
	cout << punto[0] << endl;
	cout << punto[1] << endl;
	TPZSolVec sol;
	TPZGradSolVec dsol;
	TPZFMatrix<REAL> axes;
	cEl->ComputeSolution(qsi, sol, dsol,axes);
	cout<<"Tamaño sol:" << sol.NElements() << endl;
	cout << "sol []" << sol[0] << endl;
	cout << "Tamaño dsol:" << dsol.NElements() << endl;
	cout << "dsol []" << dsol[0] << endl;
}

/*Vecinos dimensión 2*/
void getVecinos(TPZGeoEl *geoEl, TPZVec<TPZGeoElSide> &vecinos)
{
	set<TPZGeoElSide> conjuntoVecinos;
	set<int> conjuntoIndex;
	int dim = 2;
	for (int i = geoEl->NCornerNodes(); i < geoEl->NSides(); i++)
	{
		TPZStack<TPZGeoElSide> pila;
		TPZGeoElSide temp(geoEl, i);
		temp.ComputeNeighbours(pila);
		while (pila.NElements() > 0)
		{
			TPZGeoElSide tempVecino = pila.Pop();
			int dimension = (tempVecino.Element())->Dimension();
			int index = (tempVecino.Element())->Index();
			if (dimension == dim && conjuntoIndex.find(index) == conjuntoIndex.end()) {
				conjuntoVecinos.insert(tempVecino);
				conjuntoIndex.insert(index);
			}
		}
	}
	vecinos.resize(conjuntoVecinos.size());
	set<TPZGeoElSide>::iterator iter;
	int j = 0;
	for (iter = conjuntoVecinos.begin(); iter != conjuntoVecinos.end(); iter++, j++)
		vecinos[j] = *iter;
}

template<typename T>
void printVector(TPZVec<T> vect)
{
	cout << "[";
	for (int i = 0; i < vect.NElements(); i++)
		cout << vect[i] << ",";
	cout << "]" << endl;
}

REAL computeEtha(TPZGeoEl *geoEl, TPZVec<TPZGeoElSide> vecinos)
{
	TPZVec <int> lados;//indice cero corresponde al elemento geoelemento
	REAL etha = 0.;
	TPZVec< REAL > qsi(3, 0.);
	//calculo del centro de un elemento en coordenadas master
	//geoEl->CenterPoint(cEl->Reference()->NSides() - 1, qsi);
	TPZSolVec sol;
	TPZGradSolVec dsol;
	TPZFNMatrix<3, REAL> axes;
	TPZMaterialData material;
	
	for (int i = 0; i < vecinos.NElements(); i++);
	return etha;
}

void UniformHRefinementSelection(const int nDiv, TPZGeoMesh* gmesh, TPZVec<int64_t>& indexsubelements) {
	if (!indexsubelements.NElements() || !gmesh) {
		return;
	}
	if (nDiv != 1) {
		std::cout << "\nNeed implementation.\n"; return;
	}

	//	TPZStack<int64_t> subelements;
	TPZManVector<TPZGeoEl*> filhos;
	// Loop over (number of refinements) nDiv - It refines all geometric elements without subelements
//	for (int D = 0; D < nDiv; D++)
//	{
	int64_t nels = indexsubelements.NElements();
	for (int64_t elem = 0; elem < nels; elem++)
	{
		int64_t index = indexsubelements[elem];
		if (index < 0) continue;
		TPZCompEl* cel = gmesh->Reference()->ElementVec()[indexsubelements[elem]];
		TPZGeoEl* gel = cel->Reference();
		if (!gel || gel->HasSubElement() || !gel->Dimension())
			continue;
		else
			gel->Divide(filhos);
	}
	//	}
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}

void printDsol(TPZGradSolVec dsol)
{
	cout << "____Imprimiendo Dsol_____" << endl;
	for (int j = 0; j < dsol.NElements(); j++)
	{
		cout << "Numero Filas:" << dsol[j].Rows() << endl;
		cout << "Numero Columnas:" << dsol[j].Cols() << endl;
		for (int k = 0; k < dsol[j].Rows(); k++) {
			for (int l = 0; l < dsol[j].Cols(); l++)
				cout << dsol[j].GetVal(k, l) << " ";
			cout<<endl;
		}
	}
	cout << "_______________________" << endl;
}


void printAxes(TPZFNMatrix<3, REAL> axes)
{
	cout << "...Imprimiendo Axes..." << endl;
	cout << "Numero Filas:" << axes.Rows() << endl;
	cout << "Numero Columnas:" << axes.Cols() << endl;
	for (int k = 0; k < axes.Rows(); k++) {
		for (int l = 0; l < axes.Cols(); l++)
			cout << axes.GetVal(k, l) << " ";
		cout << endl;
	}
	cout << "..................." << endl;
}

void IdentifyingElementsByEthakBig(TPZAnalysis &an, REAL Tolerance, TPZManVector<int64_t>& elfathers)
{
	int64_t i, nels = an.Mesh()->NElements();
	elfathers.Resize(nels);
	int64_t nelstorefine = 0;
	TPZVec<TPZGeoElSide> vecinos;
	for (i = 0; i < nels; i++) {
		TPZInterpolatedElement * intEl;
		TPZCompEl * cEl = an.Mesh()->ElementVec()[i];
		if (!cEl || cEl->Dimension() != an.Mesh()->Dimension() || cEl->IsInterface())
			continue;
		TPZGeoEl *geoEl = an.Mesh()->Reference()->ElementVec()[i];
		TPZCompEl *comEl = an.Mesh()->ElementVec()[i];
		if (geoEl->Dimension() == 2) {
			getVecinos(geoEl, vecinos);
			//TPZInterpolatedElement *intel = (TPZInterpolatedElement*)comEl;
			
		}
	}
	elfathers.Resize(nelstorefine);
	//elfathers[0] = 1;
}

