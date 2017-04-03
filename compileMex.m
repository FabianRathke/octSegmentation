cd collector/helperFunctions
mex getRaster.c
cd ../../

cd prediction/predVariationalSubFunctions/CCode
%mex optQCC.c
mex optQCC.c CFLAGS="\$CFLAGS -march=native -Wall -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
mex calcFuncValsC.c CFLAGS="\$CFLAGS -march=native -Wall -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
mex calcCondQB.c CFLAGS="\$CFLAGS -march=native -Wall -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
mex calcMuAB2.c CFLAGS="\$CFLAGS -march=native -Wall -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
cd ../../../

cd training/helperFunctions
mex cglasso.c
mex getCondTransMatrixC.c

cd ../../
cd prediction/helperFunctions
mex sumProductSparseC.c CFLAGS="\$CFLAGS -march=native -Wall -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" 
%mex sumProductSparseC.c
cd ../../
