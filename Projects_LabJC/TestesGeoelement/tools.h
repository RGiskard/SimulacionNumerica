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
//test neopz
using namespace std;
//Se cambió de long a Int64_t
typedef int64_t tipoEntrada;
void CreateGeoMesh(TPZGeoMesh *gmesh) {
	///Crear nodos

	const int nnodes = 8;
	/*crear un arreglo de nodos (ternas) y define geometria*/
	/* */
	double coord[nnodes][3] = { { 1.,0.,0. },{ 3.,0.,0. },{ 5.,0.,0. },{ 8.,0.,0. },{ 8.,1.,0. },{ 5.,1.,0. },{ 2.,1.,0. },{ 1.,1.,0. } };
	for (int i = 0; i < nnodes; i++) {
		int nodind = gmesh->NodeVec().AllocateNewElement();//establece elementos posibles en la malla
		//cout << nodind << endl;
		/*Adapta los elementos para indexarlos a un
		tpzManvector
		*/
		TPZManVector<REAL, 3> nodeCoord(3);
		nodeCoord[0] = coord[i][0];
		nodeCoord[1] = coord[i][1];
		nodeCoord[2] = coord[i][2];

		// asocia cordenadas del nodo a la malla
		gmesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gmesh);
	}

	///Crear elementos
	const int nel = 6;
	int els[nel][3] = { { 0,6,7 },{ 0,1,6 },{ 1,5,6 },{ 1,2,5 },{ 2,4,5 },{ 2,3,4 } };//similar a una matriz
	//acceso en el for 
	for (int iel = 0; iel < nel; iel++) {
		TPZManVector<tipoEntrada, 3> nodind(3);
		tipoEntrada index;
		nodind[0] = els[iel][0];
		nodind[1] = els[iel][1];
		nodind[2] = els[iel][2];

		//E index se modificará dado que pasa por referencia
		gmesh->CreateGeoElement(ETriangle, nodind, 4, index);
	}

	///Crear elementos de contorno
	const int nelbc = 8;
	int bcels[nelbc][3] = { { 0,1,-20 },{ 1,2,-30 },{ 2,3,-20 },{ 3,4,-32 },{ 4,5,-33 },{ 5,6,-33 },{ 6,7,-30 },{ 7,0,-31 } };
	for (int iel = 0; iel < nelbc; iel++) {
		TPZManVector<tipoEntrada, 2> nodind(2);
		tipoEntrada index;
		nodind[0] = bcels[iel][0];
		nodind[1] = bcels[iel][1];
		int matid = bcels[iel][2];
		gmesh->CreateGeoElement(EOned, nodind, matid, index);
	}

	///Construyendo la conectividad de la malla
	gmesh->BuildConnectivity();


	//cout << "Se ejecutó la función" << endl;
}

//Refinamiento de la malla geométrica
void RefinaGeoMesh(TPZGeoMesh *gmesh, int ndiv) // función refina malla
{
	///Inicializando padrões de refinamento uniforme
	gRefDBase.InitializeUniformRefPattern(EOned);
	gRefDBase.InitializeUniformRefPattern(ETriangle);
	//const int ndiv = 2;

	for (int i = 0; i < ndiv; i++) {
		int nel = gmesh->NElements();
		for (int iel = 0; iel < nel; iel++) {
			TPZGeoEl * gel = gmesh->ElementVec()[iel];
			if (!gel) continue;
			if (gel->HasSubElement()) continue;//para no dividir elementos ya divididos
			TPZVec<TPZGeoEl*> filhos;
			gel->Divide(filhos);
			//para publicar información de cada elemento de la malla refinada
			string Ruta = "F:\\ResultadosTestAQP\\";
			string filename = "geoElement" + std::to_string(i) + " " + std::to_string(iel) + ".txt";
			string exit = Ruta + filename;
			ofstream os(exit.c_str());
			gel->Print(os);


		}///iel
	}///i
	ofstream myfile("geoRefinado.txt");
	gmesh->Print(myfile);
}

void testGeoEL(TPZGeoEl * geoEl)
{
	cout << "Testando geoElemento:" << endl;
	cout<<"Indice en el vector:"<<geoEl->Index() << endl;
	cout << "Identificador:" << geoEl->Id() << endl;
	int lados = geoEl->NSides();
	int total = 0;
	cout << "---------------------------" << endl;
	TPZVec<REAL> qsi(3, 0.);//dar puntos
	int lado = 3;
	if (geoEl != NULL)
		geoEl->CenterPoint(lado, qsi);
	cout << "El lado: " << lado << endl;
	cout << qsi[0] << endl;
	cout << qsi[1] << endl;
	cout << "***Cambio de cordenadas***" << endl;
	TPZVec<REAL> coor(2, 0.);
	geoEl->X(qsi, coor);
	cout << coor[0] << endl;
	cout << coor[1] << endl;
	cout << "***Calculo Gradiente de la transformacio\'on X en qsi***" << endl;
	TPZFMatrix<double> gradx(2, 2);

	geoEl->GradX(qsi, gradx);
	for (int i = 0; i < gradx.Rows(); i++) {
		for (int j = 0; j < gradx.Cols(); j++)
			cout << gradx.GetVal(0, 0) << " ";
		cout << endl;
	}
	cout << "---------------------------" << endl;	
}

void CreateCompMesh(TPZCompMesh *cmesh)
{
	cout << "Direccion de la malla computacional " << cmesh << endl;

}
