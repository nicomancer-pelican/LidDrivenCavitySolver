##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=LidDrivenCavitySolver
ConfigurationName      :=Debug
WorkspacePath          :=/home/np3217/LidDrivenCavitySolver
ProjectPath            :=/home/np3217/LidDrivenCavitySolver/LidDrivenCavitySolver
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Nicole Pellizzon
Date                   :=22/03/20
CodeLitePath           :=/home/np3217/.codelite
LinkerName             :=g++
SharedObjectLinkerName :=g++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.o.i
DebugSwitch            :=-gstab
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=$(IntermediateDirectory)/$(ProjectName)
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E 
ObjectsFileList        :="LidDrivenCavitySolver.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  -llapack  -lblas
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := ar rcus
CXX      := g++
CC       := gcc
CXXFLAGS :=  -g -O0 -std=c++14 -Wall $(Preprocessors)
CFLAGS   :=  -g -O0 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IntermediateDirectory)/LidDrivenCavity.cpp$(ObjectSuffix) $(IntermediateDirectory)/PoissonSolver.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


$(IntermediateDirectory)/.d:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/main.cpp$(ObjectSuffix): main.cpp $(IntermediateDirectory)/main.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/np3217/LidDrivenCavitySolver/LidDrivenCavitySolver/main.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/main.cpp$(DependSuffix): main.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/main.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/main.cpp$(DependSuffix) -MM main.cpp

$(IntermediateDirectory)/main.cpp$(PreprocessSuffix): main.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/main.cpp$(PreprocessSuffix) main.cpp

$(IntermediateDirectory)/LidDrivenCavity.cpp$(ObjectSuffix): LidDrivenCavity.cpp $(IntermediateDirectory)/LidDrivenCavity.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/np3217/LidDrivenCavitySolver/LidDrivenCavitySolver/LidDrivenCavity.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/LidDrivenCavity.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/LidDrivenCavity.cpp$(DependSuffix): LidDrivenCavity.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/LidDrivenCavity.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/LidDrivenCavity.cpp$(DependSuffix) -MM LidDrivenCavity.cpp

$(IntermediateDirectory)/LidDrivenCavity.cpp$(PreprocessSuffix): LidDrivenCavity.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/LidDrivenCavity.cpp$(PreprocessSuffix) LidDrivenCavity.cpp

$(IntermediateDirectory)/PoissonSolver.cpp$(ObjectSuffix): PoissonSolver.cpp $(IntermediateDirectory)/PoissonSolver.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/np3217/LidDrivenCavitySolver/LidDrivenCavitySolver/PoissonSolver.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/PoissonSolver.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/PoissonSolver.cpp$(DependSuffix): PoissonSolver.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/PoissonSolver.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/PoissonSolver.cpp$(DependSuffix) -MM PoissonSolver.cpp

$(IntermediateDirectory)/PoissonSolver.cpp$(PreprocessSuffix): PoissonSolver.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/PoissonSolver.cpp$(PreprocessSuffix) PoissonSolver.cpp


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


