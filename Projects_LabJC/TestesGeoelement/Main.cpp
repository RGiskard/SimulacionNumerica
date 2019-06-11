#include"tools.h"


int main()
{
	TPZGeoMesh * gmesh=new TPZGeoMesh();
	CreateGeoMesh(gmesh);
	string Ruta = "D:\\ResultadosTestAQP\\salida\\";
	ofstream salidatemp(Ruta+"MallaInfo.txt");//Para archivo informac�on malla
	gmesh->Print(salidatemp);
	
	/*Esta nos muestra informaci�n de la malla*/
	
	ofstream salidaVTK(Ruta+"MallaGeometrica.vtk");//archivo salida vtk
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, salidaVTK, true);
	const int ndiv = 1;
	
	//RefinaGeoMesh(gmesh, ndiv);// refinar malla
	
	/*Para probar informaci�n de un geoelemento*/

	testGeoEL(gmesh->ElementVec()[0]);//analizaremos informaci�n del primer geoelemento

	return 0;
}