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

haruhi-pg: haruhi_pg.o mempool.o sparsedoublematrix.o solvesparsedoublematrix.o dqueue.o partition_double.o mymatrix.o parallel_common.o parallel_lu_double.o postorder.o parse_spice.o union_find.o parallel_pcg_double.o parallel.o
	gcc -static -o ../bin/haruhi-pg -O2 haruhi_pg.o mempool.o sparsedoublematrix.o solvesparsedoublematrix.o dqueue.o partition_double.o mymatrix.o parallel_lu_double.o parallel_common.o postorder.o parse_spice.o union_find.o parallel_pcg_double.o parallel.o ${LIBGDSL} ${LIBAMD} ${LIBMETIS} ${LIBGSL} -lm -lpthread 
parse_spice.o: parse_spice.c sparsedoublematrix.h union_find.h
	gcc -c -O2 parse_spice.c
haruhi_pg.o: haruhi_pg.c mempool.h sparsedoublematrix.h solvesparsedoublematrix.h parallel_lu_double.h parse_spice.h parallel_pcg_double.h partition_double.h haruhi_pg.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADMETIS} haruhi_pg.c
mempool.o: mempool.c mempool.h dqueue.h
	gcc -c -O2 mempool.c -std=c99
mymatrix.o: mymatrix.c mymatrix.h
	gcc -c -O2 mymatrix.c
sparsedoublematrix.o: sparsedoublematrix.c sparsedoublematrix.h mempool.h mymatrix.h parallel.h
	gcc -c -O2 sparsedoublematrix.c
parallel_lu_double.o: parallel_lu_double.c parallel_lu_double.h mempool.h sparsedoublematrix.h partition_double.h parallel_common.h postorder.h solvesparsedoublematrix.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADMETIS} parallel_lu_double.c
parallel_pcg_double.o: parallel_pcg_double.c parallel_pcg_double.h parallel_lu_double.h mempool.h sparsedoublematrix.h partition_double.h parallel_common.h postorder.h solvesparsedoublematrix.h 
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADMETIS} parallel_pcg_double.c
partition_double.o: partition_double.h partition_double.c sparsedoublematrix.h mempool.h mymatrix.h
	gcc -c -O2 -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADGDSL} -I ${HEADMETIS} partition_double.c
solvesparsedoublematrix.o: solvesparsedoublematrix.h solvesparsedoublematrix.c mempool.h sparsedoublematrix.h partition_double.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADAMD} -I ${HEADUFCONFIG} -I ${HEADGDSL} -I ${HEADMETIS} solvesparsedoublematrix.c
parallel_common.o: parallel_common.c parallel_common.h mempool.h partition_double.h postorder.h
	gcc -c -O2 -I ${HEADGSL} -I ${HEADGDSL} -I ${HEADMETIS} parallel_common.c
parallel.o: parallel.h parallel.c
	gcc -c -O2 parallel.c
dqueue.o: dqueue.h dqueue.c
	gcc -c -O2 dqueue.c
postorder.o: postorder.h postorder.c
	gcc -c -O2 postorder.c
union_find.o: union_find.c union_find.h
	gcc -c -O2 union_find.c
clean:
	rm *.o
