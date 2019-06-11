# Jorge  -  Common functions
	add_subdirectory(Commom_JC)
	include_directories(${PZSOURCEDIR}/Projects/ABelinha/Projects_LabJC/Commom_JC)

#Jorge - Biblioteca para adaptatividade
	add_subdirectory(AdaptivityLib)
	include_directories(${PZSOURCEDIR}/Projects/ABelinha/Projects_LabJC/AdaptivityLib)


#Jorge 2017 - Implementing error analysis to MHM Method
	add_subdirectory(MHM_Jorge)			

#Jorge 2018 - Implementing error analysis from Paredes Thesis
	add_subdirectory(HDiv_Jorge)							
#Testing generating ellipse from a vector of points by least square method
	add_subdirectory(EllipseGenerator)

# Tutorial about hp-adaptive techniques
	add_subdirectory(Poisson3D_Fichera)
	add_subdirectory(Poisson3D_Shock)
	add_subdirectory(Poisson3D_AutoHP)

	add_subdirectory(Elasticity2D_Wav_Gid)

	add_subdirectory(Laplacian_3D)
	add_subdirectory(Laplacian_3D_Wav)
	add_subdirectory(Laplacian_3D_ErrorReal)
	add_subdirectory(Laplacian_3D_ErrorEstConMedia)

#Project developed by Olivier
	#add_subdirectory(OlivierRoussel)

#Projects to Conservation Laws
	add_subdirectory(ThesisPhD)
	add_subdirectory(EquationTransport)

#Proyectos grupo simulacion numerica
	add_subdirectory(TestAQP)
	add_subdirectory(Babuska)
	add_subdirectory(FengPoisson)
	add_subdirectory(TestesGeoelement)