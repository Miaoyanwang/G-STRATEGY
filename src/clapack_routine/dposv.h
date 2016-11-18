#ifndef __CLAPACK_DGESVD_H
#define __CLAPACK_DGESVD_H

#include "f2c.h"
#include "blaswrap.h"

#ifdef __cplusplus 	
extern "C" {	
#endif	


/* Subroutine */ int dposv_(char *uplo, integer *n, integer *nrhs, doublereal 
                                *a, integer *lda, doublereal *b, integer *ldb, integer *info);
#ifdef __cplusplus
}
#endif


#endif /* __CLAPACK_DGESVD_H */
