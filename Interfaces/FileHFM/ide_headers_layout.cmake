set(Project_Headers 
 	${Base_Headers}
 	${Base_Implementation_Headers}
 	${Specializations_Headers}
 	${Specializations_Implementation_Headers}
 	${Experimental_Headers}
 	${Experimental_Implementation_Headers}
 	${ExtraAlgorithms_Headers}
 	${ExtraAlgorithms_Implementation_Headers}
 	${Root_Headers}

 	${DataStructures_Headers}
 	${LinearAlgebra_Headers}
 	${Output_Headers}
 	${Macros_Headers}
)


 # ---- Headers IDE layout -----
source_group("DataStructures" FILES ${DataStructures_Headers})
source_group("Macros" FILES ${Macros_Headers})
source_group("LinearAlgebra" FILES ${LinearAlgebra_Headers})
source_group("Output" FILES ${Output_Headers})

source_group("Base" FILES ${Base_Headers})
source_group("Base\\Implementation" FILES ${Base_Implementation_Headers})
source_group("Specializations" FILES ${Specializations_Headers})
source_group("Specializations\\Implementation" FILES ${Specializations_Implementation_Headers})
source_group("Experimental" FILES ${Experimental_Headers})
source_group("Experimental\\Implementation" FILES ${Experimental_Implementation_Headers})
source_group("ExtraAlgorithms" FILES ${ExtraAlgorithms_Headers})
source_group("ExtraAlgorithms\\Implementation" FILES ${ExtraAlgorithms_Implementation_Headers})
