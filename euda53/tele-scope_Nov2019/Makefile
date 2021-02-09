
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)

ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS   = $(shell $(ROOTSYS)/bin/root-config --glibs)

# -g for gdb
# -pg for gprof
# -std=c++11

CXXFLAGS = -O2 -Wall -Wextra $(ROOTCFLAGS) -I../main/include/
#CXXFLAGS = -g -Wall -Wextra $(ROOTCFLAGS) -I../main/include/
#-I/home/pitzl/GBL/V01-17-00/cpp/include

#exportLD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/user/g/gennai/Phase2/Desy/euda53/lib/
#export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH


tele: tele.cc
	g++ tele.cc $(CXXFLAGS) -o tele \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@install_name_tool -add_rpath @executable_path ./tele
   	@echo 'done: tele'

scope53m_25x100: scope53m_25x100.cc
	g++ $(CXXFLAGS) scope53m_25x100.cc -o scope53m_25x100 \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@install_name_tool -add_rpath @executable_path ./scope53m_25x100
	@echo 'done: scope53m_25x100'

scope53m_25x100_test: scope53m_25x100_test.cc
	g++ $(CXXFLAGS) scope53m_25x100_test.cc -o scope53m_25x100_test \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@install_name_tool -add_rpath @executable_path ./scope53m_25x100_test
	@echo 'done: scope53m_25x100_test'

scope53m_50x50: scope53m_50x50.cc
	g++ $(CXXFLAGS) scope53m_50x50.cc -o scope53m_50x50 \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@install_name_tool -add_rpath @executable_path ./scope53m_50x50
	@echo 'done: scope53m_50x50'

scope53m: scope53m.cc
	g++ $(CXXFLAGS) scope53m.cc -o scope53m \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@install_name_tool -add_rpath @executable_path ./scope53m
	@echo 'done: scope53m'

scope53m_test: scope53m_test.cc
	g++ $(CXXFLAGS) scope53m_test.cc -o scope53m_test \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@install_name_tool -add_rpath @executable_path ./scope53m_test
	@echo 'done: scope53m_test'


scope53m_deltaRays: scope53m_deltaRays.cc
	g++ $(CXXFLAGS) scope53m_deltaRays.cc -o scope53m_deltaRays \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@install_name_tool -add_rpath @executable_path ./scope53m_deltaRays
	@echo 'done: scope53m_deltaRays'

ed53: ed53.cc
	g++ $(CXXFLAGS) ed53.cc -o ed53 \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@echo 'done: ed53'

edg53: edg53.cc
	g++ $(CXXFLAGS) edg53.cc -o edg53 \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@install_name_tool -add_rpath @executable_path ./edg53
	@echo 'done: edg53'

edg53_original: edg53_original.cc
	g++ $(CXXFLAGS) edg53_original.cc -o edg53_original \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@install_name_tool -add_rpath @executable_path ./edg53_original
	@echo 'done: edg53_original'
