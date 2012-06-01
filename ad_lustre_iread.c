/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *
 *   Copyright (C) 2004 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "adio.h"

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_SIGNAL_H
#include <signal.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#include <time.h>

#include "../../mpi-io/mpioimpl.h"
#include "../../mpi-io/mpioprof.h"

/* ADIOI_LUSTRE_IreadContig
 *
 * When direct I/O is in effect, this code reverts to a blocking read
 * by calling ADIOI_LUSTRE_ReadContig; otherwise,
 * ADIOI_GEN_IreadContig is called to perform the non-blocking
 * operation.
 */
void ADIOI_LUSTRE_IreadContig(ADIO_File fd, void *buf, int count,
			    MPI_Datatype datatype, int file_ptr_type,
			    ADIO_Offset offset, ADIO_Request *request,
			    int *error_code)
{
    MPI_Offset nbytes=0;
    int typesize;

    if (!(fd->direct_read || fd->direct_write)) {
        ADIOI_GEN_IreadContig(fd, buf, count, datatype, file_ptr_type,
                              offset, request, error_code);
    } else {
        ADIOI_LUSTRE_ReadContig(fd, buf, count, datatype, file_ptr_type,
                                offset, MPI_STATUS_IGNORE, error_code);
        if (*error_code == MPI_SUCCESS) {
            MPI_Type_size(datatype, &typesize);
            nbytes = (MPI_Offset)count*(MPI_Offset)typesize;
        }
        MPIO_Completed_request_create(&fd, nbytes, error_code, request);
    }
}
