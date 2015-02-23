
PROGRAM1=extractData
PROGSRC1=extractData.cpp


PROGRAM2=submitJobs
PROGSRC2=submitJobs.cpp

PROGRAM3=simulation
PROGSRC3=lib/simulation.cpp


PROGRAM4=submitSpecificJobs
PROGSRC4=submitSpecificJobs.cpp

CXXFLAGS+=-g -O3 -w



all:$(PROGRAM1) $(PROGRAM2) $(PROGRAM3) #$(PROGRAM4)

$(PROGRAM1) : $(PROGSRC1)
	$(CXX) $(CXXFLAGS) $(PROGSRC1) -o $(PROGRAM1)

$(PROGRAM2) : $(PROGSRC2)
	$(CXX) $(CXXFLAGS) $(PROGSRC2) -o $(PROGRAM2)

$(PROGRAM3) : $(PROGSRC3)
	$(CXX) $(CXXFLAGS) $(PROGSRC3) -o $(PROGRAM3)

$(PROGRAM4) : $(PROGSRC4)
	$(CXX) $(CXXFLAGS) $(PROGSRC4) -o $(PROGRAM4)

clean :
	rm -f $(PROGRAM1)
	rm -f $(PROGRAM2)
	rm -f $(PROGRAM3)
	rm -f $(PROGRAM4)

cleanCondor :
	rm -f condor/log/*.condor
	rm -f condor/error/*.condor
	rm -f condor/output/*.condor
	rm -f condor/runJob/*.condor