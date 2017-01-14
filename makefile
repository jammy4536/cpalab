CXX=g++
CFLAGS=-g -fbounds-check

CpAjm13368 : CpAjm13368_Master.cpp
	$(CXX) $(CFLAGS) -o  $@ $<

test: CpAjm13368 input.txt
	< input.txt ./CpAjm13368
