/*Solución del ejemplo 3.5 de babuska*/
#include "tools.h"


void crearMallaGeometrica(TPZGeoMesh *& mallageometrica)
{
	const int nodos = 16;
	REAL cordenadasNodo[nodos][3]=
	{ 
	{0.,0.,0.},
	{1/3.,0.,0.}, 
	{2.0/3,0.,0.},
	{1.0,0.,0.},
	{0.,1./3,0.},
	{1./3,1./3,0.},
	{2./3,1./3,0.},
	{1.0,1./3,0.},
	{0.,2./3,0.},
	{1./3,2./3,0.},
	{2./3,2./3,0.},
	{1.,2./3,0.},
	{0.,1.,0.},
	{1./3,1,0.},
	{2./3,1.,0.},
	{1.,1.,0.},
	};
	/*
	for (int i=0;i<16;i++)
	{
		for (int j = 0; j < 3; j++)
			cout << cordenadasNodo[i][j] << " ";
		cout << endl;
	}
	*/
	for (int i = 0; i < nodos; i++)
	{
		int nodid = mallageometrica->NodeVec().AllocateNewElement();
		TPZManVector<REAL, 3> coordenadasVec(3);
		coordenadasVec[0] = cordenadasNodo[i][0];
		coordenadasVec[1] = cordenadasNodo[i][1];
		coordenadasVec[2] = cordenadasNodo[i][2];
		mallageometrica->NodeVec()[nodid].Initialize(i,coordenadasVec,*mallageometrica);
	}
	const int numeroElementos = 18;
	int elementosMatrix[numeroElementos][3] =
	{
	{0,1,4},
	{1,5,4},
	{1,6,5},
	{1,2,6},
	{2,7,6},
	{2,3,7},
	{4,5,9},
	{4,9,8},
	{5,6,10},
	{5,10,9},
	{6,7,11},
	{6,11,10},
	{8,9,13},
	{8,13,12},
	{9,10,13},
	{10,14,13},
	{10,11,15},
	{10,15,14},
	};
	for (int iel = 0; iel < numeroElementos; iel++) {
		TPZManVector<int64_t, 3> nodosDelemento(3);
		int64_t index;
		nodosDelemento[0] = elementosMatrix[iel][0];
		nodosDelemento[1] = elementosMatrix[iel][1];
		nodosDelemento[2] = elementosMatrix[iel][2];

		mallageometrica->CreateGeoElement(ETriangle, nodosDelemento, 4, index);
	}
	const int totalElsFront = 12;
	int elementosFrontera[totalElsFront][3] =
	{
		{0,1,-1},
		{1,2,-2},
		{2,3,-3},
		{3,7,-4},
		{7,11,-5},
		{11,15,-6},
		{15,14,-7},
		{14,13,-8},
		{13,12,-9},
		{12,8,-10},
		{8,4,-11},
		{4,0,-12},
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



/*Funciones para condiciones de frontera*/
static void condicionNeumann(const TPZVec<REAL> &loc, TPZVec<STATE> &result) {  
	result.Resize(1);
	result[0] = 0;
}

static void fun(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du) {
	const REAL x1 = loc[0];
	const REAL y1 = loc[1];
	const REAL sol = x1 * x1 + y1 * y1;
	u[0] = sol;
}

template<typename T>
class function:TPZFunction<T>
{
	
};



void crearMallaComputacional(TPZCompMesh *& mallaComputacional)
{
	const int dimension = 2;
	const int idMaterial = 4;
	TPZMatLaplacian *material = new TPZMatLaplacian(idMaterial,dimension);
	const STATE K = 2,F=300 ;
	material->SetParameters(K,F);
	const int pOrden = 1;
	//TPZAutoPointer<TPZFunction<STATE>> fun2;
	//material->SetForcingFunction(fun, pOrden);
	mallaComputacional->InsertMaterialObject(material);
	
	/*hipotesis: val2 esta asociado a dirichelet, val1´para newman*/
	TPZFMatrix<STATE> val1(1, 1, -1000), val2(1, 1, 30.);//antes val 1 val1(1, 1, 3.)
	TPZBndCond * condDirichlet3_3a7 = material->CreateBC(material, -4,0, val1, val2);
	mallaComputacional->InsertMaterialObject(condDirichlet3_3a7);

	TPZBndCond * condDirichlet3_7a11 = material->CreateBC(material, -5, 0, val1, val2);
	mallaComputacional->InsertMaterialObject(condDirichlet3_7a11);

	TPZBndCond * condDirichlet3_11a15 = material->CreateBC(material, -6, 0, val1, val2);
	mallaComputacional->InsertMaterialObject(condDirichlet3_11a15);

	TPZBndCond * condDirichlet3_15a14 = material->CreateBC(material, -7, 0, val1, val2);
	mallaComputacional->InsertMaterialObject(condDirichlet3_15a14);

	TPZBndCond * condDirichlet3_14a13 = material->CreateBC(material, -8, 0, val1, val2);
	mallaComputacional->InsertMaterialObject(condDirichlet3_14a13);

	TPZBndCond * condDirichlet3_13a12 = material->CreateBC(material, -9, 0, val1, val2);
	mallaComputacional->InsertMaterialObject(condDirichlet3_13a12);
	////-------

	TPZMaterial * condNeumannCero_0a1 = material->CreateBC(material, -1, 1, val1, val2);//1 = Neumann
	condNeumannCero_0a1->SetForcingFunction(condicionNeumann,pOrden);
	mallaComputacional->InsertMaterialObject(condNeumannCero_0a1);

	TPZMaterial * condNeumannCero_1a2 = material->CreateBC(material, -2, 1, val1, val2);//1 = Neumann
	condNeumannCero_0a1->SetForcingFunction(condicionNeumann, pOrden);
	mallaComputacional->InsertMaterialObject(condNeumannCero_1a2);

	TPZMaterial * condNeumannCero_2a3 = material->CreateBC(material, -3, 1, val1, val2);//1 = Neumann
	condNeumannCero_0a1->SetForcingFunction(condicionNeumann, pOrden);
	mallaComputacional->InsertMaterialObject(condNeumannCero_2a3);

	TPZMaterial * condNeumannCero_0a4 = material->CreateBC(material, -12, 1, val1, val2);//1 = Neumann
	condNeumannCero_0a1->SetForcingFunction(condicionNeumann, pOrden);
	mallaComputacional->InsertMaterialObject(condNeumannCero_0a4);

	TPZMaterial * condNeumannCero_4a8 = material->CreateBC(material, -11, 1, val1, val2);//1 = Neumann
	condNeumannCero_0a1->SetForcingFunction(condicionNeumann, pOrden);
	mallaComputacional->InsertMaterialObject(condNeumannCero_4a8);

	TPZMaterial * condNeumannCero_8a12 = material->CreateBC(material, -10, 1, val1, val2);//1 = Neumann
	condNeumannCero_0a1->SetForcingFunction(condicionNeumann, pOrden);
	mallaComputacional->InsertMaterialObject(condNeumannCero_8a12);
	
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
	mallaGeo->Print(mallaInfo);
	TPZVTKGeoMesh::PrintGMeshVTK(mallaGeo,salidaVTK,true);
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
	TPZAutoPointer< TPZStructMatrix > MatrizEstructurada = an.StructMatrix();
	cout << "Ecuaciones despues delfiltro: " << MatrizEstructurada->NReducedEquations() << endl;
	
	TPZFMatrix< REAL > 	fillin;
	mallaComp->ComputeFillIn(1,fillin);

	fillin.Print(cout);

	TPZFMatrix<REAL> solution= an.Solution();
	TPZMatrixSolver< STATE > *state = &an.Solver();
	ofstream posible("MatrizPosible.txt");
	(state->Matrix())->Print(posible);
	/*Impresion de la matriz solución*/
	solution.Print("Solución",cout);

	

	TPZVec<string> names;
	TPZVec<string> names2;
	//an.AnimateRun(10, names, names2, "salida");
	const string FileNameSolver = "MatrizK.txt";
	ofstream ossolver("MatrizK.txt");
	an.Print(FileNameSolver, ossolver);

	cout << "Numero de threads: " << matriz.GetNumThreads() << endl;
	set<int> materiales = matriz.MaterialIds();
	set<int>::iterator iter;
	for (iter = materiales.begin(); iter != materiales.end(); iter++)
		cout << *iter << endl;


	std::ofstream anPostProcessFile("postprocess.txt");

		///Exportando para Paraview
	TPZVec<std::string> scalarVars(1), vectorVars(0);
	scalarVars[0] = "Solution";
	an.DefineGraphMesh(2, scalarVars, vectorVars, "solucao.vtk");
	int resolution = 0;
	an.PostProcess(resolution);

	ofstream salidaSolucion("Solucao.txt");
	an.Print("solucion.txt", salidaSolucion);

	delete mallaComp;
	delete mallaGeo;
	return 0;
}