/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   Copyright (C) 2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 *
 *   Copyright (C) 2007 Oak Ridge National Laboratory
 *
 *   Copyright (C) 2008 Sun Microsystems, Lustre group
 */

#include "ad_lustre.h"
#include <shmem.h>

struct ADIOI_Fns_struct ADIO_LUSTRE_operations = {
    ADIOI_LUSTRE_Open, /* Open */
    ADIOI_GEN_OpenColl, /* OpenColl */
    ADIOI_LUSTRE_ReadContig, /* ReadContig */
    ADIOI_LUSTRE_WriteContig, /* WriteContig */
    ADIOI_GEN_ReadStridedColl, /* ReadStridedColl */
    ADIOI_LUSTRE_WriteStridedColl, /* WriteStridedColl */
    ADIOI_GEN_SeekIndividual, /* SeekIndividual */
    ADIOI_GEN_Fcntl, /* Fcntl */
    ADIOI_LUSTRE_SetInfo, /* SetInfo */
    ADIOI_GEN_ReadStrided, /* ReadStrided */
    ADIOI_LUSTRE_WriteStrided, /* WriteStrided */
    ADIOI_GEN_Close, /* Close */
#if defined(ROMIO_HAVE_WORKING_AIO) && !defined(CRAY_XT_LUSTRE)
    ADIOI_LUSTRE_IreadContig, /* IreadContig */
    ADIOI_LUSTRE_IwriteContig, /* IwriteContig */
#else
    ADIOI_FAKE_IreadContig, /* IreadContig */
    ADIOI_FAKE_IwriteContig, /* IwriteContig */
#endif
    ADIOI_GEN_IODone, /* ReadDone */
    ADIOI_GEN_IODone, /* WriteDone */
    ADIOI_GEN_IOComplete, /* ReadComplete */
    ADIOI_GEN_IOComplete, /* WriteComplete */
    ADIOI_GEN_IreadStrided, /* IreadStrided */
    ADIOI_GEN_IwriteStrided, /* IwriteStrided */
    ADIOI_GEN_Flush, /* Flush */
    ADIOI_GEN_Resize, /* Resize */
    ADIOI_GEN_Delete, /* Delete */
    ADIOI_GEN_Feature, /* Features */
};

void *LUSTRE_shmalloc_fn(size_t size, int lineno, const char *fname)
{
    void *new;
    new = shmalloc(size);
    if (!new) {
	    fprintf(stderr, "Out of symmetric heap memory in file %s, line %d\n",
                fname, lineno);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    DBG_FPRINTF(stderr, "LUSTRE_shmalloc %s:<%d> %p (%#zX)\n", fname, lineno, new, size);
    return new;
}

void LUSTRE_shfree_fn(void *ptr, int lineno, const char *fname)
{
    DBG_FPRINTF(stderr, "LUSTRE_shfree %s:<%d> %p\n", fname, lineno, ptr);
    if (!ptr) {
        fprintf(stderr, "Attempt to free symmetric heap null pointer in file %s, line %d\n", fname, lineno);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    shfree(ptr);
}
