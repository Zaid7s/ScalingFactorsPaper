F90      = gfortran
F77      = gfortran

F90FLAGS = -g 													\
		   -ffpe-trap='zero','overflow','underflow','invalid'	\
		   -static												\
		   -c 													\
		   -std=f2008 											\
		   -O 													\
		   -W 													\
		   -fcheck=all 											\
		   -Wline-truncation 									\
		   -Wuninitialized 										\
		   -fbounds-check 										\
		   -Wconversion 										\
		   -Wsurprising 										\
		   -Wpedantic 											\
		   -Warray-temporaries 									\
		   -Wconversion-extra 									\
		   -Wextra 												\
		   -Wunused-parameter 									\
		   -fbacktrace 											\
		   -static

F77FLAGS  = -c -O -fbounds-check -fcheck-array-temporaries -Wconversion -Wsurprising -Wall
LDFLAGS  = 
LIBS     =

MOD_FILES = ModuleFiles
OBJ_FILES = ObjectFiles

COMPILE_OBJ = $(F90) $(F90FLAGS) $(MODS)

LINK_AND_LOAD = $(F90) $(LDFLAGS) $(LIBS)

CMD = run.exe
#working

MV_MOD = mv *.mod $(MOD_FILES)
MV_OBJ = mv *.o $(OBJ_FILES)

OBJS = $(OBJ_FILES)/Main.o \
	   $(OBJ_FILES)/GetGlobalParameter.o \
	   $(OBJ_FILES)/GetNormalizeMean.o \
	   $(OBJ_FILES)/GetMeanValues.o \
	   $(OBJ_FILES)/GetMean_And_Pet.o \
	   $(OBJ_FILES)/akima433.o \
	   $(OBJ_FILES)/Akima1D.o \
	   $(OBJ_FILES)/GetAreaDerivative.o \
	   $(OBJ_FILES)/GetdQdT_BDF.o \
	   $(OBJ_FILES)/GetJacobians.o \
	   $(OBJ_FILES)/GetLUSGS.o \
	   $(OBJ_FILES)/GetSolveLUSGS.o \
	   $(OBJ_FILES)/GetDX.o \
	   $(OBJ_FILES)/GetiMaxValue.o \
       $(OBJ_FILES)/Compressible1DFunctions.o \
	   $(OBJ_FILES)/GetExactSol.o \
	   $(OBJ_FILES)/GetLMatrix.o \
	   $(OBJ_FILES)/GetUMatrix.o \
	   $(OBJ_FILES)/GetDMatrix.o \
	   $(OBJ_FILES)/GetRHSVector.o \
	   $(OBJ_FILES)/GetDissipation.o \
	   $(OBJ_FILES)/GetdEdXVector.o \
	   $(OBJ_FILES)/GetInflowBC.o \
 	   $(OBJ_FILES)/GetInflowBC_Mean_And_Pet.o \
	   $(OBJ_FILES)/GetLDU_Vectors.o \
	   $(OBJ_FILES)/GetDeltaA.o \
	   $(OBJ_FILES)/GetSpatialNOdQdT.o \
	   $(OBJ_FILES)/GetSpatialDerivative.o \
	   $(OBJ_FILES)/GetCompactScheme.o \
	   $(OBJ_FILES)/GetSolvePreFac.o \
	   $(OBJ_FILES)/GetOutflowBC.o \
	   $(OBJ_FILES)/GetCFL_Ramping.o \
 	   $(OBJ_FILES)/GetOutflowBC_Mean_And_Pet.o \
	   $(OBJ_FILES)/ForwardSweep.o \
	   $(OBJ_FILES)/BackwardSweep.o \
	   $(OBJ_FILES)/GetRungeKutta.o \
	   $(OBJ_FILES)/GetTriDi.o \
	   $(OBJ_FILES)/GetFixBoundaries.o \
	   $(OBJ_FILES)/GetInverse.o \
	   $(OBJ_FILES)/GetDiag.o \
	   $(OBJ_FILES)/fft1d.o \
	   $(OBJ_FILES)/Precision_Def.o \
	   $(OBJ_FILES)/BlockTriDiSolver1D.o

MODS = \
	   -I$(MOD_FILES)

$(CMD): $(OBJS)
	($(LINK_AND_LOAD) -o  $(@) $(OBJS));
	(./compile.sh)

$(OBJ_FILES)/Main.o: \
 Main.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o \
 $(OBJ_FILES)/GetNormalizeMean.o \
 $(OBJ_FILES)/GetMean_And_Pet.o \
 $(OBJ_FILES)/GetMeanValues.o
	$(COMPILE_OBJ) Main.f90
	$(MV_OBJ)

$(OBJ_FILES)/GetGlobalParameter.o: \
 GetGlobalParameter.f90
	$(COMPILE_OBJ) GetGlobalParameter.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetMeanValues.o: \
 GetMeanValues.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o \
 $(OBJ_FILES)/GetLUSGS.o	\
 $(OBJ_FILES)/GetJacobians.o	\
 $(OBJ_FILES)/GetTriDi.o	\
 $(OBJ_FILES)/GetDX.o	\
 $(OBJ_FILES)/GetCFL_Ramping.o	\
 $(OBJ_FILES)/GetiMaxValue.o	\
 $(OBJ_FILES)/GetExactSol.o	\
 $(OBJ_FILES)/GetAreaDerivative.o	
	$(COMPILE_OBJ) GetMeanValues.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/Akima1D.o: \
 Akima1D.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o \
 $(OBJ_FILES)/akima433.o
	$(COMPILE_OBJ) Akima1D.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/akima433.o: \
 akima433.f
	$(F77) $(F77FLAGS) akima433.f
	$(MV_OBJ)

$(OBJ_FILES)/GetiMaxValue.o: \
 GetiMaxValue.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetiMaxValue.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetDX.o: \
 GetDX.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetDX.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetExactSol.o: \
 GetExactSol.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o \
 $(OBJ_FILES)/Compressible1DFunctions.o
	$(COMPILE_OBJ) GetExactSol.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/Compressible1DFunctions.o: \
 Compressible1DFunctions.f90 
	$(COMPILE_OBJ) Compressible1DFunctions.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetAreaDerivative.o: \
 GetAreaDerivative.f90 \
 $(OBJ_FILES)/BlockTriDiSolver1D.o \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetAreaDerivative.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetJacobians.o: \
 GetJacobians.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetJacobians.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetLUSGS.o: \
 GetLUSGS.f90 \
 $(OBJ_FILES)/GetRHSVector.o \
 $(OBJ_FILES)/GetGlobalParameter.o \
 $(OBJ_FILES)/GetRungeKutta.o \
 $(OBJ_FILES)/GetSolveLUSGS.o \
 $(OBJ_FILES)/GetLDU_Vectors.o \
 $(OBJ_FILES)/GetLMatrix.o \
 $(OBJ_FILES)/GetDMatrix.o \
 $(OBJ_FILES)/GetUMatrix.o \
 $(OBJ_FILES)/ForwardSweep.o \
 $(OBJ_FILES)/BackwardSweep.o 
	$(COMPILE_OBJ) GetLUSGS.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetSolveLUSGS.o: \
 $(OBJ_FILES)/GetGlobalParameter.o \
 GetSolveLUSGS.f90
	$(COMPILE_OBJ) GetSolveLUSGS.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetLMatrix.o: \
 GetLMatrix.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetLMatrix.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetUMatrix.o: \
 GetUMatrix.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetUMatrix.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetDMatrix.o: \
 GetDMatrix.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetDMatrix.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetRHSVector.o: \
 GetRHSVector.f90 \
 $(OBJ_FILES)/GetdEdXVector.o \
 $(OBJ_FILES)/GetdQdT_BDF.o \
 $(OBJ_FILES)/GetGlobalParameter.o \
 $(OBJ_FILES)/GetDissipation.o
	$(COMPILE_OBJ) GetRHSVector.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetdQdT_BDF.o: \
 GetdQdT_BDF.f90	\
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetdQdT_BDF.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetSolvePreFac.o: \
 GetSolvePreFac.f90
	$(COMPILE_OBJ) GetSolvePreFac.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetSpatialNOdQdT.o: \
 GetSpatialNOdQdT.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetSpatialNOdQdT.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetDissipation.o: \
 GetDissipation.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetDissipation.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetLDU_Vectors.o: \
 $(OBJ_FILES)/GetGlobalParameter.o \
 GetLDU_Vectors.f90
	$(COMPILE_OBJ) GetLDU_Vectors.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetdEdXVector.o: \
 GetdEdXVector.f90 \
 $(OBJ_FILES)/GetSpatialDerivative.o	\
 $(OBJ_FILES)/GetCompactScheme.o	\
 $(OBJ_FILES)/GetInflowBC.o \
 $(OBJ_FILES)/GetGlobalParameter.o \
 $(OBJ_FILES)/GetOutflowBC.o
	$(COMPILE_OBJ) GetdEdXVector.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetInflowBC.o: \
 GetInflowBC.f90 \
 $(OBJ_FILES)/BlockTriDiSolver1D.o \
 $(OBJ_FILES)/GetGlobalParameter.o \
 $(OBJ_FILES)/GetDeltaA.o
	$(COMPILE_OBJ) GetInflowBC.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetDeltaA.o: \
 GetDeltaA.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetDeltaA.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetInflowBC_Mean_And_Pet.o: \
 GetInflowBC_Mean_And_Pet.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o \
 $(OBJ_FILES)/GetDeltaA.o
	$(COMPILE_OBJ) GetInflowBC_Mean_And_Pet.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetSpatialDerivative.o: \
 GetSpatialDerivative.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetSpatialDerivative.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetCompactScheme.o: \
 GetCompactScheme.f90 \
 $(OBJ_FILES)/GetSolvePreFac.o \
 $(OBJ_FILES)/BlockTriDiSolver1D.o  
	$(COMPILE_OBJ) GetCompactScheme.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetOutflowBC.o: \
 GetOutflowBC.f90 \
 $(OBJ_FILES)/BlockTriDiSolver1D.o \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetOutflowBC.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetOutflowBC_Mean_And_Pet.o: \
 GetOutflowBC_Mean_And_Pet.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetOutflowBC_Mean_And_Pet.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/ForwardSweep.o: \
 ForwardSweep.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) ForwardSweep.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/BackwardSweep.o: \
 BackwardSweep.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) BackwardSweep.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetRungeKutta.o: \
 GetRungeKutta.f90 \
 $(OBJ_FILES)/GetdEdXVector.o \
 $(OBJ_FILES)/GetdQdT_BDF.o \
 $(OBJ_FILES)/GetGlobalParameter.o \
 $(OBJ_FILES)/GetSpatialNOdQdT.o \
 $(OBJ_FILES)/GetDissipation.o
	$(COMPILE_OBJ) GetRungeKutta.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/fft1d.o:   \
 fft1d.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o \
 $(OBJ_FILES)/Precision_Def.o
	$(COMPILE_OBJ) fft1d.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/Precision_Def.o: \
 Precision_Def.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ)  Precision_Def.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetCFL_Ramping.o: \
 GetCFL_Ramping.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetCFL_Ramping.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetTriDi.o: \
 GetTriDi.f90 \
 $(OBJ_FILES)/GetRHSVector.o	\
 $(OBJ_FILES)/GetRungeKutta.o \
 $(OBJ_FILES)/GetFixBoundaries.o \
 $(OBJ_FILES)/GetDiag.o \
 $(OBJ_FILES)/GetGlobalParameter.o \
 $(OBJ_FILES)/BlockTriDiSolver1D.o
	$(COMPILE_OBJ) GetTriDi.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetFixBoundaries.o: \
 GetFixBoundaries.f90			\
 $(OBJ_FILES)/GetGlobalParameter.o \
 $(OBJ_FILES)/GetInverse.o
	$(COMPILE_OBJ) GetFixBoundaries.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetInverse.o: \
 GetInverse.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetInverse.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetDiag.o: \
 GetDiag.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) GetDiag.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/BlockTriDiSolver1D.o: \
 BlockTriDiSolver1D.f90 \
 $(OBJ_FILES)/GetGlobalParameter.o 
	$(COMPILE_OBJ) BlockTriDiSolver1D.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetMean_And_Pet.o: \
 GetMean_And_Pet.f90 \
 $(OBJ_FILES)/GetLUSGS.o	\
 $(OBJ_FILES)/GetJacobians.o	\
 $(OBJ_FILES)/GetTriDi.o	\
 $(OBJ_FILES)/GetDX.o	\
 $(OBJ_FILES)/GetiMaxValue.o	\
 $(OBJ_FILES)/GetExactSol.o	\
 $(OBJ_FILES)/fft1d.o \
 $(OBJ_FILES)/Precision_Def.o \
 $(OBJ_FILES)/GetGlobalParameter.o \
 $(OBJ_FILES)/GetAreaDerivative.o	
	$(COMPILE_OBJ) GetMean_And_Pet.f90
	$(MV_OBJ); $(MV_MOD)

$(OBJ_FILES)/GetNormalizeMean.o: \
 GetNormalizeMean.f90 \
 $(OBJ_FILES)/GetLUSGS.o	\
 $(OBJ_FILES)/GetJacobians.o	\
 $(OBJ_FILES)/GetTriDi.o	\
 $(OBJ_FILES)/GetDX.o	\
 $(OBJ_FILES)/GetiMaxValue.o	\
 $(OBJ_FILES)/GetExactSol.o	\
 $(OBJ_FILES)/fft1d.o \
 $(OBJ_FILES)/Precision_Def.o \
 $(OBJ_FILES)/GetGlobalParameter.o \
 $(OBJ_FILES)/GetAreaDerivative.o	
	$(COMPILE_OBJ) GetNormalizeMean.f90
	$(MV_OBJ); $(MV_MOD)

z:
	rm -f *.o 
	rm -f *.dat 
	rm -f *.mod 
	rm -f *.exe 
	rm -f *.eps 
	rm -f TimeStepSol/*.dat
	rm -f ROC/*.dat

clean:
	rm -f *.o 
	rm -f *.dat 
	rm -f *.mod 
	rm -f *.exe 
	rm -f *.eps 
	rm -f $(MOD_FILES)/*.mod
	rm -f $(OBJ_FILES)/*.o
	rm -f TimeStepSol/Normalize/*.dat
	rm -f ROC/*.dat
	rm -f FinalSol/*.dat
	rm -f FFT/*.dat
