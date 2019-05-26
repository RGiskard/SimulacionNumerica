/*Solución del ejemplo notas de Feng*/
#include "tools.h"

class PoissonVariable :public TPZMatLaplacian
{
	public:
		PoissonVariable(int nummat, int dim) :
			TPZMatLaplacian(nummat,dim) {}
};

void crearMallaGeometrica(TPZGeoMesh *& mallageometrica)
{
	const int nodos = 4;
	REAL cordenadasNodo[nodos][3]=
	{ 
	{-1.,-1.,0.},
	{1.,-1.,0.}, 
	{1.,1.,0.},
	{-1.0,1.,0.}
	};

	for (int i=0;i<4;i++)
	{
		for (int j = 0; j < 3; j++)
			cout << cordenadasNodo[i][j] << " ";
		cout << endl;
	}
	
	for (int i = 0; i < nodos; i++)
	{
		int nodid = mallageometrica->NodeVec().AllocateNewElement();
		TPZManVector<REAL, 3> coordenadasVec(3);
		coordenadasVec[0] = cordenadasNodo[i][0];
		coordenadasVec[1] = cordenadasNodo[i][1];
		coordenadasVec[2] = cordenadasNodo[i][2];
		mallageometrica->NodeVec()[nodid].Initialize(i,coordenadasVec,*mallageometrica);
	}
	const int numeroElementos = 2;
	int elementosMatrix[numeroElementos][3] =
	{
	{0,1,3},
	{1,2,3},
	};
	for (int iel = 0; iel < numeroElementos; iel++) {
		TPZManVector<int64_t, 3> nodosDelemento(3);
		int64_t index;
		nodosDelemento[0] = elementosMatrix[iel][0];
		nodosDelemento[1] = elementosMatrix[iel][1];
		nodosDelemento[2] = elementosMatrix[iel][2];

		mallageometrica->CreateGeoElement(ETriangle, nodosDelemento, 4, index);
	}
	const int totalElsFront = 4;
	int elementosFrontera[totalElsFront][3] =
	{
		{0,1,-1},
		{1,2,-2},
		{2,3,-3},
		{3,0,-4}
	};
	for (int iel = 0; iel < totalElsFront; iel++) {
		TPZManVector<int64_t, 2> nodosElFront(2);
		int64_t index;
		nodosElFront[0] = elementosFrontera[iel][0];
		nodosElFront[1] = elementosFrontera[iel][1];
		int idMaterial = elementosFrontera[iel][2];
		mallageometrica->CreateGeoElement(EOned, nodosElFront, idMaterial, index);
	}
	mallageometrica->BuildConnectivity();
}

//Refinamiento de la malla geométrica
void RefinaGeoMesh(TPZGeoMesh *gmesh, int ndiv) // función refina malla
{
	///Inicializando padrões de refinamento uniforme
	gRefDBase.InitializeUniformRefPattern(EOned);
	gRefDBase.InitializeUniformRefPattern(ETriangle);
	//const int ndiv = 3;

	for (int i = 0; i < ndiv; i++) {
		int nel = gmesh->NElements();
		for (int iel = 0; iel < nel; iel++) {
			TPZGeoEl * gel = gmesh->ElementVec()[iel];
			if (!gel) continue;
			if (gel->HasSubElement()) continue;//para no dividir elementos ya divididos
			TPZVec<TPZGeoEl*> filhos;
			gel->Divide(filhos);
		}///iel
	}///i
	{
		ofstream myfile("geoRefinado.txt");
		gmesh->Print(myfile);
		ofstream salidaVTK("mallaGeometricaRefinada.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh, salidaVTK, true);
	}
}

/*
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
*/

void funcionFuerza(const TPZVec<REAL> &loc, TPZVec<STATE> &u) {
	const REAL x1 = loc[0];
	const REAL y1 = loc[1];
	const REAL value = -x1 -y1;
	u[0] = value;
	//cout << "El valor de result:" << u[0] << endl;
}

static void funDirichletAbajo(const TPZVec<REAL> &loc, TPZVec<STATE> &result) {  
	const REAL x1 = loc[0];
	const REAL y1 = loc[1];
	result.Resize(1);
	result[0] =(1.0/6.0)*(x1*x1*x1-1); 
}



static void funDirichletDerecha(const TPZVec<REAL> &loc, TPZVec<STATE> &result) {
	const REAL x1 = loc[0];
	const REAL y1 = loc[1];
	result.Resize(1);
	result[0] = (1.0 / 6.0)*(y1*y1*y1 +1);
}




static void funNewmannArriba(const TPZVec<REAL> &loc, TPZVec<STATE> &result) {
	const REAL x1 = loc[0];
	const REAL y1 = loc[1];
	result.Resize(1);
	result[0] = 1.0/2.0;
}


static void funNewmannIzquierda(const TPZVec<REAL> &loc, TPZVec<STATE> &result) {
	const REAL x1 = loc[0];
	const REAL y1 = loc[1];
	result.Resize(1);
	result[0] =- 1.0 / 2.0;
}

void crearMallaComputacional(TPZCompMesh *& mallaComputacional)
{
	const int dimension = 2;
	const int idMaterial = 4;
	TPZMatLaplacian *material = new TPZMatLaplacian(idMaterial, dimension);
	const STATE K = 1, F = 1;
	const int pOrden = 1;

	material->SetParameters(K, F);
	material->SetForcingFunction(new TPZDummyFunction<STATE>(funcionFuerza, pOrden));
	mallaComputacional->InsertMaterialObject(material);
	/*hipotesis: val2 esta asociado a dirichelet, val1´para newman*/
	TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);//antes val 1 val1(1, 1, 3.)

	/*Condiciones Dirichelet*/
	TPZMaterial * condDirichletAbajo = material->CreateBC(material, -1, 0, val1, val2);
	condDirichletAbajo->SetForcingFunction(funDirichletAbajo,pOrden);
	mallaComputacional->InsertMaterialObject(condDirichletAbajo);

	TPZMaterial * condDirichletDerecha = material->CreateBC(material, -2, 0, val1, val2);
	condDirichletDerecha->SetForcingFunction(funDirichletDerecha, pOrden);
	mallaComputacional->InsertMaterialObject(condDirichletDerecha);

	/*Condiciones newmann*/
	TPZMaterial * condNewmannArriba = material->CreateBC(material, -3, 1, val1, val2);
	condNewmannArriba->SetForcingFunction(funNewmannArriba, pOrden);
	mallaComputacional->InsertMaterialObject(condNewmannArriba);

	TPZMaterial * condNewmannIzquierda = material->CreateBC(material, -4, 1, val1, val2);
	condNewmannIzquierda->SetForcingFunction(funNewmannIzquierda, pOrden);
	mallaComputacional->InsertMaterialObject(condNewmannIzquierda);

	mallaComputacional->SetDefaultOrder(pOrden);
	mallaComputacional->SetDimModel(dimension);
	mallaComputacional->SetAllCreateFunctionsContinuous();
	mallaComputacional->AutoBuild();

	{
		std::ofstream myfile("geoPosAutoBuild.txt");
		mallaComputacional->Reference()->Print(myfile);
	}

	{
		std::ofstream myfile("cmesh.txt");
		mallaComputacional->Print(myfile);
	}

}


int main()
{
	TPZGeoMesh * mallaGeo =new TPZGeoMesh() ;
	crearMallaGeometrica(mallaGeo);
	ofstream salidaVTK("mallaGeometrica.vtk");
	ofstream mallaInfo("mallaInfo.txt");
	TPZVTKGeoMesh::PrintGMeshVTK(mallaGeo, salidaVTK, true);
	mallaGeo->Print(mallaInfo);
	RefinaGeoMesh(mallaGeo, 4);
	TPZCompMesh *mallaComp = new TPZCompMesh(mallaGeo);
	crearMallaComputacional(mallaComp);

	TPZAnalysis an(mallaComp);
	TPZSkylineNSymStructMatrix matriz(mallaComp);
	an.SetStructuralMatrix(matriz);
	///Decomposio LU
	TPZStepSolver<STATE> step;
	step.SetDirect(ELU);
	an.SetSolver(step);

	///Assemble da matriz de rigidez e vetor de carga
	an.Assemble();


	///Resoluo do sistema devuelve una matriz
	an.Solve();
	TPZFMatrix< STATE > vectorF = an.Rhs();
	vectorF.Print(cout);
	TPZFMatrix<REAL> solution = an.Solution();
	ofstream solucionDisk("sol4refin.txt");
	solution.Print(solucionDisk);

	TPZVec<std::string> scalarVars(1), vectorVars(0);
	scalarVars[0] = "Solution";
	an.DefineGraphMesh(2, scalarVars, vectorVars, "solucao.vtk");

	
	delete mallaComp;
	delete mallaGeo;
	return 0;
}