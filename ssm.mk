##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Release
ProjectName            :=ssm
ConfigurationName      :=Release
WorkspacePath          :=/home/taylor/Documents/ssmworkspace
ProjectPath            :=/home/taylor/ssm
IntermediateDirectory  :=./Release
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=taylor
Date                   :=09/01/18
CodeLitePath           :=/home/taylor/.codelite
LinkerName             :=/usr/bin/x86_64-linux-gnu-g++
SharedObjectLinkerName :=/usr/bin/x86_64-linux-gnu-g++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.i
DebugSwitch            :=-g 
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=$(IntermediateDirectory)/lib$(ProjectName).a
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="ssm.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch)/home/taylor/ssm/include $(IncludeSwitch)/home/taylor/ssm/include2 $(IncludeSwitch)/usr/include/eigen3 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/x86_64-linux-gnu-ar rcu
CXX      := /usr/bin/x86_64-linux-gnu-g++
CC       := /usr/bin/x86_64-linux-gnu-gcc
CXXFLAGS :=  -pg -O3 -std=c++11 $(Preprocessors)
CFLAGS   :=   $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/x86_64-linux-gnu-as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/src_utilities_transformations.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_utilities_convenience_funcs.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_distributions_multinomial_resampler.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_distributions_densities.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_distributions_pmfs.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_mcmc_bases2_ada_pmmh.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_filter_bases_apf_smooth.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_filter_bases_apf_filter.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_filter_bases_hmm_rbpf.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_filter_bases_fshmm.cpp$(ObjectSuffix) \
	$(IntermediateDirectory)/src_filter_bases_sisr_smooth.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_filter_bases_kalman_rbpf.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_filter_bases_bootstrap_filter.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_filter_bases_lgssm.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_filter_bases_sisr_filter.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_mcmc_bases_pmmh.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(IntermediateDirectory) $(OutputFile)

$(OutputFile): $(Objects)
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(AR) $(ArchiveOutputSwitch)$(OutputFile) @$(ObjectsFileList) $(ArLibs)
	@$(MakeDirCommand) "/home/taylor/Documents/ssmworkspace/.build-release"
	@echo rebuilt > "/home/taylor/Documents/ssmworkspace/.build-release/ssm"

MakeIntermediateDirs:
	@test -d ./Release || $(MakeDirCommand) ./Release


./Release:
	@test -d ./Release || $(MakeDirCommand) ./Release

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/src_utilities_transformations.cpp$(ObjectSuffix): src/utilities/transformations.cpp $(IntermediateDirectory)/src_utilities_transformations.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/taylor/ssm/src/utilities/transformations.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_utilities_transformations.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_utilities_transformations.cpp$(DependSuffix): src/utilities/transformations.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_utilities_transformations.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_utilities_transformations.cpp$(DependSuffix) -MM src/utilities/transformations.cpp

$(IntermediateDirectory)/src_utilities_transformations.cpp$(PreprocessSuffix): src/utilities/transformations.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_utilities_transformations.cpp$(PreprocessSuffix) src/utilities/transformations.cpp

$(IntermediateDirectory)/src_utilities_convenience_funcs.cpp$(ObjectSuffix): src/utilities/convenience_funcs.cpp $(IntermediateDirectory)/src_utilities_convenience_funcs.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/taylor/ssm/src/utilities/convenience_funcs.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_utilities_convenience_funcs.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_utilities_convenience_funcs.cpp$(DependSuffix): src/utilities/convenience_funcs.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_utilities_convenience_funcs.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_utilities_convenience_funcs.cpp$(DependSuffix) -MM src/utilities/convenience_funcs.cpp

$(IntermediateDirectory)/src_utilities_convenience_funcs.cpp$(PreprocessSuffix): src/utilities/convenience_funcs.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_utilities_convenience_funcs.cpp$(PreprocessSuffix) src/utilities/convenience_funcs.cpp

$(IntermediateDirectory)/src_distributions_multinomial_resampler.cpp$(ObjectSuffix): src/distributions/multinomial_resampler.cpp $(IntermediateDirectory)/src_distributions_multinomial_resampler.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/taylor/ssm/src/distributions/multinomial_resampler.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_distributions_multinomial_resampler.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_distributions_multinomial_resampler.cpp$(DependSuffix): src/distributions/multinomial_resampler.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_distributions_multinomial_resampler.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_distributions_multinomial_resampler.cpp$(DependSuffix) -MM src/distributions/multinomial_resampler.cpp

$(IntermediateDirectory)/src_distributions_multinomial_resampler.cpp$(PreprocessSuffix): src/distributions/multinomial_resampler.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_distributions_multinomial_resampler.cpp$(PreprocessSuffix) src/distributions/multinomial_resampler.cpp

$(IntermediateDirectory)/src_distributions_densities.cpp$(ObjectSuffix): src/distributions/densities.cpp $(IntermediateDirectory)/src_distributions_densities.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/taylor/ssm/src/distributions/densities.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_distributions_densities.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_distributions_densities.cpp$(DependSuffix): src/distributions/densities.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_distributions_densities.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_distributions_densities.cpp$(DependSuffix) -MM src/distributions/densities.cpp

$(IntermediateDirectory)/src_distributions_densities.cpp$(PreprocessSuffix): src/distributions/densities.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_distributions_densities.cpp$(PreprocessSuffix) src/distributions/densities.cpp

$(IntermediateDirectory)/src_distributions_pmfs.cpp$(ObjectSuffix): src/distributions/pmfs.cpp $(IntermediateDirectory)/src_distributions_pmfs.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/taylor/ssm/src/distributions/pmfs.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_distributions_pmfs.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_distributions_pmfs.cpp$(DependSuffix): src/distributions/pmfs.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_distributions_pmfs.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_distributions_pmfs.cpp$(DependSuffix) -MM src/distributions/pmfs.cpp

$(IntermediateDirectory)/src_distributions_pmfs.cpp$(PreprocessSuffix): src/distributions/pmfs.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_distributions_pmfs.cpp$(PreprocessSuffix) src/distributions/pmfs.cpp

$(IntermediateDirectory)/src_mcmc_bases2_ada_pmmh.cpp$(ObjectSuffix): src/mcmc_bases2/ada_pmmh.cpp $(IntermediateDirectory)/src_mcmc_bases2_ada_pmmh.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/taylor/ssm/src/mcmc_bases2/ada_pmmh.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_mcmc_bases2_ada_pmmh.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_mcmc_bases2_ada_pmmh.cpp$(DependSuffix): src/mcmc_bases2/ada_pmmh.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_mcmc_bases2_ada_pmmh.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_mcmc_bases2_ada_pmmh.cpp$(DependSuffix) -MM src/mcmc_bases2/ada_pmmh.cpp

$(IntermediateDirectory)/src_mcmc_bases2_ada_pmmh.cpp$(PreprocessSuffix): src/mcmc_bases2/ada_pmmh.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_mcmc_bases2_ada_pmmh.cpp$(PreprocessSuffix) src/mcmc_bases2/ada_pmmh.cpp

$(IntermediateDirectory)/src_filter_bases_apf_smooth.cpp$(ObjectSuffix): src/filter_bases/apf_smooth.cpp $(IntermediateDirectory)/src_filter_bases_apf_smooth.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/taylor/ssm/src/filter_bases/apf_smooth.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_filter_bases_apf_smooth.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_filter_bases_apf_smooth.cpp$(DependSuffix): src/filter_bases/apf_smooth.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_filter_bases_apf_smooth.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_filter_bases_apf_smooth.cpp$(DependSuffix) -MM src/filter_bases/apf_smooth.cpp

$(IntermediateDirectory)/src_filter_bases_apf_smooth.cpp$(PreprocessSuffix): src/filter_bases/apf_smooth.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_filter_bases_apf_smooth.cpp$(PreprocessSuffix) src/filter_bases/apf_smooth.cpp

$(IntermediateDirectory)/src_filter_bases_apf_filter.cpp$(ObjectSuffix): src/filter_bases/apf_filter.cpp $(IntermediateDirectory)/src_filter_bases_apf_filter.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/taylor/ssm/src/filter_bases/apf_filter.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_filter_bases_apf_filter.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_filter_bases_apf_filter.cpp$(DependSuffix): src/filter_bases/apf_filter.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_filter_bases_apf_filter.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_filter_bases_apf_filter.cpp$(DependSuffix) -MM src/filter_bases/apf_filter.cpp

$(IntermediateDirectory)/src_filter_bases_apf_filter.cpp$(PreprocessSuffix): src/filter_bases/apf_filter.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_filter_bases_apf_filter.cpp$(PreprocessSuffix) src/filter_bases/apf_filter.cpp

$(IntermediateDirectory)/src_filter_bases_hmm_rbpf.cpp$(ObjectSuffix): src/filter_bases/hmm_rbpf.cpp $(IntermediateDirectory)/src_filter_bases_hmm_rbpf.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/taylor/ssm/src/filter_bases/hmm_rbpf.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_filter_bases_hmm_rbpf.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_filter_bases_hmm_rbpf.cpp$(DependSuffix): src/filter_bases/hmm_rbpf.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_filter_bases_hmm_rbpf.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_filter_bases_hmm_rbpf.cpp$(DependSuffix) -MM src/filter_bases/hmm_rbpf.cpp

$(IntermediateDirectory)/src_filter_bases_hmm_rbpf.cpp$(PreprocessSuffix): src/filter_bases/hmm_rbpf.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_filter_bases_hmm_rbpf.cpp$(PreprocessSuffix) src/filter_bases/hmm_rbpf.cpp

$(IntermediateDirectory)/src_filter_bases_fshmm.cpp$(ObjectSuffix): src/filter_bases/fshmm.cpp $(IntermediateDirectory)/src_filter_bases_fshmm.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/taylor/ssm/src/filter_bases/fshmm.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_filter_bases_fshmm.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_filter_bases_fshmm.cpp$(DependSuffix): src/filter_bases/fshmm.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_filter_bases_fshmm.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_filter_bases_fshmm.cpp$(DependSuffix) -MM src/filter_bases/fshmm.cpp

$(IntermediateDirectory)/src_filter_bases_fshmm.cpp$(PreprocessSuffix): src/filter_bases/fshmm.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_filter_bases_fshmm.cpp$(PreprocessSuffix) src/filter_bases/fshmm.cpp

$(IntermediateDirectory)/src_filter_bases_sisr_smooth.cpp$(ObjectSuffix): src/filter_bases/sisr_smooth.cpp $(IntermediateDirectory)/src_filter_bases_sisr_smooth.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/taylor/ssm/src/filter_bases/sisr_smooth.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_filter_bases_sisr_smooth.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_filter_bases_sisr_smooth.cpp$(DependSuffix): src/filter_bases/sisr_smooth.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_filter_bases_sisr_smooth.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_filter_bases_sisr_smooth.cpp$(DependSuffix) -MM src/filter_bases/sisr_smooth.cpp

$(IntermediateDirectory)/src_filter_bases_sisr_smooth.cpp$(PreprocessSuffix): src/filter_bases/sisr_smooth.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_filter_bases_sisr_smooth.cpp$(PreprocessSuffix) src/filter_bases/sisr_smooth.cpp

$(IntermediateDirectory)/src_filter_bases_kalman_rbpf.cpp$(ObjectSuffix): src/filter_bases/kalman_rbpf.cpp $(IntermediateDirectory)/src_filter_bases_kalman_rbpf.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/taylor/ssm/src/filter_bases/kalman_rbpf.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_filter_bases_kalman_rbpf.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_filter_bases_kalman_rbpf.cpp$(DependSuffix): src/filter_bases/kalman_rbpf.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_filter_bases_kalman_rbpf.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_filter_bases_kalman_rbpf.cpp$(DependSuffix) -MM src/filter_bases/kalman_rbpf.cpp

$(IntermediateDirectory)/src_filter_bases_kalman_rbpf.cpp$(PreprocessSuffix): src/filter_bases/kalman_rbpf.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_filter_bases_kalman_rbpf.cpp$(PreprocessSuffix) src/filter_bases/kalman_rbpf.cpp

$(IntermediateDirectory)/src_filter_bases_bootstrap_filter.cpp$(ObjectSuffix): src/filter_bases/bootstrap_filter.cpp $(IntermediateDirectory)/src_filter_bases_bootstrap_filter.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/taylor/ssm/src/filter_bases/bootstrap_filter.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_filter_bases_bootstrap_filter.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_filter_bases_bootstrap_filter.cpp$(DependSuffix): src/filter_bases/bootstrap_filter.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_filter_bases_bootstrap_filter.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_filter_bases_bootstrap_filter.cpp$(DependSuffix) -MM src/filter_bases/bootstrap_filter.cpp

$(IntermediateDirectory)/src_filter_bases_bootstrap_filter.cpp$(PreprocessSuffix): src/filter_bases/bootstrap_filter.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_filter_bases_bootstrap_filter.cpp$(PreprocessSuffix) src/filter_bases/bootstrap_filter.cpp

$(IntermediateDirectory)/src_filter_bases_lgssm.cpp$(ObjectSuffix): src/filter_bases/lgssm.cpp $(IntermediateDirectory)/src_filter_bases_lgssm.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/taylor/ssm/src/filter_bases/lgssm.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_filter_bases_lgssm.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_filter_bases_lgssm.cpp$(DependSuffix): src/filter_bases/lgssm.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_filter_bases_lgssm.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_filter_bases_lgssm.cpp$(DependSuffix) -MM src/filter_bases/lgssm.cpp

$(IntermediateDirectory)/src_filter_bases_lgssm.cpp$(PreprocessSuffix): src/filter_bases/lgssm.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_filter_bases_lgssm.cpp$(PreprocessSuffix) src/filter_bases/lgssm.cpp

$(IntermediateDirectory)/src_filter_bases_sisr_filter.cpp$(ObjectSuffix): src/filter_bases/sisr_filter.cpp $(IntermediateDirectory)/src_filter_bases_sisr_filter.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/taylor/ssm/src/filter_bases/sisr_filter.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_filter_bases_sisr_filter.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_filter_bases_sisr_filter.cpp$(DependSuffix): src/filter_bases/sisr_filter.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_filter_bases_sisr_filter.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_filter_bases_sisr_filter.cpp$(DependSuffix) -MM src/filter_bases/sisr_filter.cpp

$(IntermediateDirectory)/src_filter_bases_sisr_filter.cpp$(PreprocessSuffix): src/filter_bases/sisr_filter.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_filter_bases_sisr_filter.cpp$(PreprocessSuffix) src/filter_bases/sisr_filter.cpp

$(IntermediateDirectory)/src_mcmc_bases_pmmh.cpp$(ObjectSuffix): src/mcmc_bases/pmmh.cpp $(IntermediateDirectory)/src_mcmc_bases_pmmh.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/taylor/ssm/src/mcmc_bases/pmmh.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_mcmc_bases_pmmh.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_mcmc_bases_pmmh.cpp$(DependSuffix): src/mcmc_bases/pmmh.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_mcmc_bases_pmmh.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_mcmc_bases_pmmh.cpp$(DependSuffix) -MM src/mcmc_bases/pmmh.cpp

$(IntermediateDirectory)/src_mcmc_bases_pmmh.cpp$(PreprocessSuffix): src/mcmc_bases/pmmh.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_mcmc_bases_pmmh.cpp$(PreprocessSuffix) src/mcmc_bases/pmmh.cpp


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Release/


