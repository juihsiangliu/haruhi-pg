HEADGSL = ../lib/gsl-1.9/
LIBGSL = ../lib/gsl-1.9/.libs/libgsl.a ../lib/gsl-1.9/cblas/.libs/libgslcblas.a
HEADGDSL = ../lib/gdsl-1.4/src/
LIBGDSL = ../lib/gdsl-1.4/src/.libs/libgdsl.a
HEADAMD = ../lib/AMD/Include/
LIBAMD = ../lib/AMD/Lib/libamd.a
HEADUFCONFIG = ../lib/UFconfig/

# HEADMETIS = ../lib/metis-4.0/Lib/
# LIBMETIS = ../lib/metis-4.0/libmetis.a

HEADMETIS = ../lib/metis-5.0pre2/include/
LIBMETIS = ../lib/metis-5.0pre2/libmetis.a

#simulator.out: main.o quadelement.o quadmatrix.o dlmread.o parser.o gsl_extern.o quadlinear.o fgetl.o quadnonlinear.o plot.o montesample.o montelinear.o montenonlinear.o parallel.o sparsequadmatrix.o
#	gcc -o simulator.out -O2 main.o quadelement.o quadmatrix.o dlmread.o parser.o gsl_extern.o quadlinear.o fgetl.o quadnonlinear.o plot.o montesample.o montelinear.o montenonlinear.o parallel.o sparsequadmatrix.o sparsedoublematrix.o ${LIBGDSL} ${LIBGSL} -lm -lpthread -pg
test.out: test.o quadelement.o quadmatrix.o dlmread.o parser.o gsl_extern.o fgetl.o plot.o sparsequadmatrix.o quadlinear.o sparsegvarient.o quadnonlinear.o mempool.o mymatrix.o parser_util.o montesample.o transfernetlist.o sparsedoublematrix.o sparsegvarient_monte.o montelinear.o montenonlinear.o parallel_lu_double.o partition_double.o partition_quad.o solvesparsedoublematrix.o solvesparsequadmatrix.o parallel_lu_quad.o parallel_lu_common.o parallel.o min_max_heap.o dqueue.o postorder.o
	gcc -o ../bin/test.out -O2 test.o quadelement.o quadmatrix.o dlmread.o parser.o gsl_extern.o fgetl.o plot.o sparsequadmatrix.o quadlinear.o sparsegvarient.o quadnonlinear.o mempool.o mymatrix.o parser_util.o montesample.o transfernetlist.o sparsedoublematrix.o sparsegvarient_monte.o montelinear.o montenonlinear.o parallel_lu_double.o partition_double.o partition_quad.o solvesparsedoublematrix.o solvesparsequadmatrix.o parallel_lu_quad.o parallel_lu_common.o parallel.o min_max_heap.o dqueue.o postorder.o ${LIBGDSL} ${LIBGSL} ${LIBAMD} ${LIBMETIS} -lm -lpthread -pg
haruhi_lu: haruhi_lu.o mempool.o sparsedoublematrix.o solvesparsedoublematrix.o dqueue.o partition_double.o mymatrix.o parallel_common.o parallel_lu_double.o postorder.o parse_spice.o union_find.o parallel_pcg_double.o parallel.o
	gcc -static -o ../bin/haruhi_lu -O2 haruhi_lu.o mempool.o sparsedoublematrix.o solvesparsedoublematrix.o dqueue.o partition_double.o mymatrix.o parallel_lu_double.o parallel_common.o postorder.o parse_spice.o union_find.o parallel_pcg_double.o parallel.o ${LIBGDSL} ${LIBAMD} ${LIBMETIS} ${LIBGSL} -lm -lpthread 
test_io: test_io.o mempool.o sparsedoublematrix.o solvesparsedoublematrix.o dqueue.o partition_double.o mymatrix.o parallel_lu_common.o parallel_lu_double.o postorder.o parse_spice.o
	gcc -o ../bin/test_io -O2 test_io.o mempool.o sparsedoublematrix.o solvesparsedoublematrix.o dqueue.o partition_double.o mymatrix.o parallel_lu_double.o parallel_lu_common.o postorder.o parse_spice.o ${LIBGDSL} ${LIBAMD} ${LIBMETIS} ${LIBGSL} -lm -lpthread 
gen2Dr: gen2Dr.o mempool.o sparsedoublematrix.o dqueue.o mymatrix.o 
	gcc -o ../bin/gen2Dr -O2 gen2Dr.o mempool.o sparsedoublematrix.o dqueue.o mymatrix.o -lm -lpthread
parse_spice: parse_spice.o sparsedoublematrix.o mempool.o dqueue.o mymatrix.o
	gcc -o ../bin/parse_spice -O2 parse_spice.o sparsedoublematrix.o mempool.o dqueue.o mymatrix.o -lm -lpthread
parse_spice.o: parse_spice.c sparsedoublematrix.h union_find.h
	gcc -c -O2 parse_spice.c
haruhi_lu.o: haruhi_lu.c mempool.h sparsedoublematrix.h solvesparsedoublematrix.h parallel_lu_double.h parse_spice.h parallel_pcg_double.h partition_double.h haruhi_lu.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADMETIS} haruhi_lu.c
gen2Dr.o: gen2Dr.c mempool.h sparsedoublematrix.h
	gcc -c -O2 gen2Dr.c
test.o: test.c quadelement.h quadmatrix.h dlmread.h parser.h gsl_extern.h quadlinear.h sparsequadmatrix.h mempool.h quadlinear.h quadnonlinear.h sparsegvarient.h montesample.h transfernetlist.h sparsedoublematrix.h sparsegvarient_monte.h montelinear.h plot.h montenonlinear.h sparsedoublematrix.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADMETIS} test.c
main.o: main.c quadelement.h quadmatrix.h dlmread.h parser.h gsl_extern.h quadlinear.h montesample.h sparsequadmatrix.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} main.c
quadelement.o: quadelement.c quadelement.h gsl_extern.h mempool.h mymatrix.h
	gcc -c -O2 -I ${HEADGSL} quadelement.c 
quadmatrix.o: quadmatrix.c quadmatrix.h quadelement.h
	gcc -c -O2 -I ${HEADGSL} quadmatrix.c
dlmread.o: dlmread.c dlmread.h fgetl.h
	gcc -c -O2 -I ${HEADGSL} dlmread.c
parser.o: parser.h parser.c fgetl.h sparsegvarient.h sparsequadmatrix.h quadmatrix.h dlmread.h quadelement.h parser_util.h mempool.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} parser.c
gsl_extern.o: gsl_extern.h gsl_extern.c
	gcc -c -O2 -I ${HEADGSL} gsl_extern.c
quadlinear.o: quadlinear.h quadlinear.c quadmatrix.h parser.h sparsequadmatrix.h partition_double.h partition_quad.h solvesparsequadmatrix.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADMETIS} quadlinear.c 
fgetl.o: fgetl.c fgetl.h
	gcc -c -O2 -c fgetl.c
quadnonlinear.o: quadnonlinear.h quadnonlinear.c quadmatrix.h parser.h mempool.h sparsequadmatrix.h mymatrix.h sparsegvarient.h solvesparsequadmatrix.h parallel.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADMETIS} quadnonlinear.c
plot.o: plot.c quadmatrix.h plot.h parser.h transfernetlist.h mymatrix.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} plot.c 
montesample.o: montesample.c montesample.h parser.h mempool.h gsl_extern.h quadelement.h mymatrix.h parser_util.h sparsedoublematrix.h 
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} montesample.c
montelinear.o: montelinear.c montelinear.h montesample.h transfernetlist.h sparsedoublematrix.h mempool.h parallel_lu_double.h partition_double.h
	gcc -c -g -O2 -I ${HEADGDSL} -I ${HEADGSL} -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADMETIS} montelinear.c -pg
montenonlinear.o: montenonlinear.h montenonlinear.c transfernetlist.h mymatrix.h mempool.h sparsegvarient_monte.h sparsedoublematrix.h solvesparsedoublematrix.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADMETIS} montenonlinear.c
sparsequadmatrix.o: sparsequadmatrix.c sparsequadmatrix.h quadmatrix.h quadelement.h mempool.h sparsedoublematrix.h
	gcc -c -O2 -I ${HEADGSL} sparsequadmatrix.c 
sparsegvarient.o: sparsegvarient.h sparsegvarient.c quadmatrix.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} sparsegvarient.c
mempool.o: mempool.c mempool.h dqueue.h
	gcc -c -O2 mempool.c -std=c99
mymatrix.o: mymatrix.c mymatrix.h
	gcc -c -O2 mymatrix.c
parser_util.o: parser_util.h parser_util.c quadelement.h quadmatrix.h
	gcc -c -O2 -I ${HEADGSL} parser_util.c
sparsedoublematrix.o: sparsedoublematrix.c sparsedoublematrix.h mempool.h mymatrix.h parallel.h
	gcc -c -O2 sparsedoublematrix.c
sparsegvarient_monte.o: sparsegvarient_monte.c sparsegvarient_monte.h quadmatrix.h mempool.h sparsegvarient.h montesample.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADGSL} sparsegvarient_monte.c
transfernetlist.o: transfernetlist.c transfernetlist.h parser.h montesample.h sparsegvarient_monte.h sparsedoublematrix.h sparsequadmatrix.h mempool.h
	gcc -c -O2 -I ${HEADGDSL} -I ${HEADGSL} transfernetlist.c
parallel_lu_double.o: parallel_lu_double.c parallel_lu_double.h mempool.h sparsedoublematrix.h partition_double.h parallel_common.h postorder.h solvesparsedoublematrix.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADMETIS} parallel_lu_double.c
parallel_pcg_double.o: parallel_pcg_double.c parallel_pcg_double.h parallel_lu_double.h mempool.h sparsedoublematrix.h partition_double.h parallel_common.h postorder.h solvesparsedoublematrix.h 
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADMETIS} parallel_pcg_double.c
partition_double.o: partition_double.h partition_double.c sparsedoublematrix.h mempool.h mymatrix.h
	gcc -c -O2 -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADGDSL} -I ${HEADMETIS} partition_double.c
partition_quad.o: partition_quad.h partition_quad.c partition_double.h sparsedoublematrix.h sparsequadmatrix.h mempool.h 
	gcc -c -O2 -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADGDSL} -I ${HEADGSL} -I ${HEADMETIS} partition_quad.c 
solvesparsedoublematrix.o: solvesparsedoublematrix.h solvesparsedoublematrix.c mempool.h sparsedoublematrix.h partition_double.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADGDSL} -I ${HEADMETIS} solvesparsedoublematrix.c
solvesparsequadmatrix.o: solvesparsequadmatrix.h solvesparsequadmatrix.c quadelement.h quadmatrix.h mempool.h sparsedoublematrix.h sparsequadmatrix.h partition_quad.h parallel_lu_quad.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADGDSL} -I ${HEADMETIS} solvesparsequadmatrix.c
parallel_lu_quad.o: parallel_lu_quad.c parallel_lu_quad.h mempool.h sparsequadmatrix.h partition_double.h parallel_lu_common.h solvesparsequadmatrix.h
	gcc -c -O2 -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADGDSL} -I ${HEADGSL} -I ${HEADMETIS} parallel_lu_quad.c
parallel_common.o: parallel_common.c parallel_common.h mempool.h partition_double.h postorder.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} -I ${HEADMETIS} parallel_common.c
parallel.o: parallel.h parallel.c
	gcc -c -O2 parallel.c
gdsl_extern.o: gdsl_extern.h gdsl_extern.c
	gcc -c -O2 -I ${HEADGDSL} gdsl_extern.c
min_max_heap.o: min_max_heap.h min_max_heap.c
	gcc -c -O2 min_max_heap.c -std=c99
dqueue.o: dqueue.h dqueue.c
	gcc -c -O2 dqueue.c
postorder.o: postorder.h postorder.c
	gcc -c -O2 postorder.c
test_io.o: test_io.c mempool.h sparsedoublematrix.h
	gcc -c -O2 test_io.c
union_find.o: union_find.c union_find.h
	gcc -c -O2 union_find.c
clean:
	rm *.o
