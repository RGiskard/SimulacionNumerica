//test neopz
#include "tools.h"
//Se cambió de long a Int64_t
typedef int64_t tipoEntrada;
void CreateGeoMesh(TPZGeoMesh *gmesh) {
	///Criando nodos

	const int nnodes = 8;
	/*crear un arreglo de nodos (ternas) y define geometria*/
	double coord[nnodes][3] = { { 1.,0.,0. },{ 3.,0.,0. },{ 5.,0.,0. },{ 8.,0.,0. },{ 8.,1.,0. },{ 5.,1.,0. },{ 2.,1.,0. },{ 1.,1.,0. } };
	for (int i = 0; i < nnodes; i++) {
		int nodind = gmesh->NodeVec().AllocateNewElement();//establece elementos posibles en la malla
		TPZManVector<REAL, 3> nodeCoord(3);
		nodeCoord[0] = coord[i][0];
		nodeCoord[1] = coord[i][1];
		nodeCoord[2] = coord[i][2];

		// asocia cordenadas del nodo a la malla
		gmesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gmesh);
	}

	///Criando elementos
	const int nel = 6;
	int els[nel][3] = { { 0,6,7 },{ 0,1,6 },{ 1,5,6 },{ 1,2,5 },{ 2,4,5 },{ 2,3,4 } };//similar a una matriz
	//acceso en el for 
	for (int iel = 0; iel < nel; iel++) {
		TPZManVector<int64_t, 3> nodind(3);
		int64_t index;
		nodind[0] = els[iel][0];
		nodind[1] = els[iel][1];
		nodind[2] = els[iel][2];

		//E index se modificara dado que pasa por referencia
		gmesh->CreateGeoElement(ETriangle, nodind, 4, index);
	}

	///Criando elementos de contorno
	const int nelbc = 8;
	int bcels[nelbc][3] = { { 0,1,-20 },{ 1,2,-30 },{ 2,3,-20 },{ 3,4,-32 },{ 4,5,-33 },{ 5,6,-33 },{ 6,7,-30 },{ 7,0,-31 } };
	for (int iel = 0; iel < nelbc; iel++) {
		TPZManVector<int64_t, 2> nodind(2);
		int64_t index;
		nodind[0] = bcels[iel][0];
		nodind[1] = bcels[iel][1];
		int matid = bcels[iel][2];
		gmesh->CreateGeoElement(EOned, nodind, matid, index);
	}

	///Construindo conectividade da malha
	gmesh->BuildConnectivity();
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
			string Ruta="D:\\salidaRTest\\";
			string FileName = "GeoElemento" + std::to_string(i) + " " + std::to_string(iel) + ".txt";
			string exit = Ruta + FileName;
			ofstream os(exit.c_str());
			string FileNameTopological = "GeoElementopological" + std::to_string(i) + " " + std::to_string(iel) + ".txt";
			ofstream osTopologicaInfo(Ruta+FileNameTopological);
			gel->Print(os);
			gel->PrintTopologicalInfo(osTopologicaInfo);

		}///iel
	}///i
	ofstream myfile("geoRefinado.txt");
	gmesh->Print(myfile);
}

void CreateCompMesh(TPZCompMesh *cmesh) // funcion crea malla computacional
{
	/*Indicacao de malha DGFem. Se false, vamos criar malha H1
	Espacio de aproximacion
	*/
	const bool DGFEM = false;//Metodo de galerkim discontinuo en caso sea true

	/// crear materiales
	const int dim = 2; 
	/*identificador de material Para frontera es negativo, */
	// del primer material
	const int matid1 = 4;
	TPZMatLaplacian * mat1 = new TPZMatLaplacian(matid1, dim);
	//utilizar ecuacion de poisson, (ver contribute)
	const STATE K1 = 1., F1 = 0.;
	mat1->SetParameters(K1, F1);//como se da de parámetro 0, se estable que es ecuación de laplace

	// el material creado es insertado en la ecuacion diferencial
	cmesh->InsertMaterialObject(mat1);
	const int pOrder = 2;

	///Condicioneses de frontera
	/*Inicializando la matriz solución val1->solucion U. val2->a su gradiente*/
	TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);

	TPZBndCond * BCondDirichletNulo = mat1->CreateBC(mat1, -3, 0, val1, val2);//0 = Dirichlet
	cmesh->InsertMaterialObject(BCondDirichletNulo);

	TPZBndCond * BCondNeumannZero = mat1->CreateBC(mat1, -2, 1, val1, val2);//1 = Neumann
	cmesh->InsertMaterialObject(BCondNeumannZero);

	TPZMaterial * BCondNeumannEsq = mat1->CreateBC(mat1, -5, 1, val1, val2);//1 = Neumann
	BCondNeumannEsq->SetForcingFunction(NeumannEsquerda, pOrder);
	cmesh->InsertMaterialObject(BCondNeumannEsq);

	TPZMaterial * BCondNeumannDir = mat1->CreateBC(mat1, -4, 1, val1, val2);//1 = Neumann
	BCondNeumannDir->SetForcingFunction(NeumannDireita, pOrder);
	cmesh->InsertMaterialObject(BCondNeumannDir);

	TPZMaterial * BCondNeumannAbaixo = mat1->CreateBC(mat1, -6, 1, val1, val2);//1 = Neumann
	BCondNeumannAbaixo->SetForcingFunction(NeumannAbaixo, pOrder);
	cmesh->InsertMaterialObject(BCondNeumannAbaixo);
	

	cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);//dim = 2
	cmesh->SetAllCreateFunctionsContinuous();
	cmesh->AutoBuild();
	{
		std::ofstream myfile("geoPosAutoBuild.txt");
		cmesh->Reference()->Print(myfile);
	}

	{
		std::ofstream myfile("cmesh.txt");
		cmesh->Print(myfile);
	}
}



int main()
{
	TPZGeoMesh * gmesh=new TPZGeoMesh();
	CreateGeoMesh(gmesh);
	string Ruta = "D:\\salidaRTest\\salida\\";
	ofstream salidatemp(Ruta+"Mallainfo.txt");
	gmesh->Print(salidatemp);
	cout << "Nro de nodos " << gmesh -> NNodes() << endl;
	cout << "Nro de elementos " << (*gmesh).NElements() << endl;
	//const int ndiv = 3;
	//RefinaGeoMesh(gmesh, ndiv);// refinar malla
	ofstream salidaVTK(Ruta+"MallaGeometricaOriginal.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, salidaVTK, true);

	const int ndiv = 1;
	RefinaGeoMesh(gmesh, ndiv);// refinar malla
	ofstream salidaRefin(Ruta+"MallaGeometricaRefinada.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, salidaRefin, true);

	ofstream infoRefin(Ruta+"infoRefin.txt");
	gmesh->Print(infoRefin);
	return 0;
}