cd collector/helperFunctions
mex getRaster.c
cd ../../

cd prediction/predVariationalSubFunctions/
mex optQCC.c
cd ../../

cd training/helperFunctions
mex cglasso.c
mex getCondTransMatrixC.c

cd ../../
cd prediction/helperFunctions
mex CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" sumProductSparseC.c
cd ../../
