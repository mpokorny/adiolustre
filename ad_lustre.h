/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   Copyright (C) 1997 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 *
 *   Copyright (C) 2007 Oak Ridge National Laboratory
 *
 *   Copyright (C) 2008 Sun Microsystems, Lustre group
 */

#ifndef AD_UNIX_INCLUDE
#define AD_UNIX_INCLUDE

/* temp*/
#define HAVE_ASM_TYPES_H 1

#include <unistd.h>
#include <linux/types.h>

#ifdef __linux__
#  include <sys/ioctl.h>                            /* necessary for: */
#  include <time.h>
#  define __USE_GNU                                 /* O_DIRECT and */
#  include <fcntl.h>                                /* IO operations */
#  undef __USE_GNU
#endif /* __linux__ */

/*#include <fcntl.h>*/
#include <sys/ioctl.h>
#include <lustre/lustre_user.h>
#include "adio.h"
/*#include "adioi.h"*/

#ifdef HAVE_SIGNAL_H
#include <signal.h>
#endif

#ifdef HAVE_AIO_H
#include <aio.h>
#ifdef HAVE_SYS_AIO_H
#include <sys/aio.h>
#endif
#endif /* End of HAVE_SYS_AIO_H */

/* TODO: does the following type definition belong here? */
struct req_info {
	void *source;
	ADIO_Offset length;
	int aggregator_index;
	int stripe_num;
	int relative_offset;
};

void ADIOI_LUSTRE_Open(ADIO_File fd, int *error_code);
void ADIOI_LUSTRE_Close(ADIO_File fd, int *error_code);
void ADIOI_LUSTRE_ReadContig(ADIO_File fd, void *buf, int count,
                             MPI_Datatype datatype, int file_ptr_type,
                             ADIO_Offset offset, ADIO_Status *status,
                             int *error_code);
void ADIOI_LUSTRE_WriteContig(ADIO_File fd, void *buf, int count,
                              MPI_Datatype datatype, int file_ptr_type,
                              ADIO_Offset offset, ADIO_Status *status,
                              int *error_code);
void ADIOI_LUSTRE_WriteStrided(ADIO_File fd, void *buf, int count,
			       MPI_Datatype datatype, int file_ptr_type,
			       ADIO_Offset offset, ADIO_Status *status,
			       int *error_code);
void ADIOI_LUSTRE_WriteStridedColl(ADIO_File fd, void *buf, int count,
		                   MPI_Datatype datatype, int file_ptr_type,
		                   ADIO_Offset offset, ADIO_Status *status,
                                   int *error_code);
void ADIOI_LUSTRE_ReadStridedColl(ADIO_File fd, void *buf, int count,
		                  MPI_Datatype datatype, int file_ptr_type,
		                  ADIO_Offset offset, ADIO_Status *status,
                                  int *error_code);
void ADIOI_LUSTRE_ReadStrided(ADIO_File fd, void *buf, int count,
			      MPI_Datatype datatype, int file_ptr_type,
			      ADIO_Offset offset, ADIO_Status *status,
                              int *error_code);
void ADIOI_LUSTRE_Fcntl(ADIO_File fd, int flag, ADIO_Fcntl_t *fcntl_struct,
	               int *error_code);
void ADIOI_LUSTRE_SetInfo(ADIO_File fd, MPI_Info users_info, int *error_code);

/* the lustre utilities: */
int ADIOI_LUSTRE_Docollect(ADIO_File fd, int contig_access_count,
                           ADIO_Offset *len_list);

void ADIOI_LUSTRE_Get_striping_info(ADIO_File fd,
                                    ADIO_Offset **striping_info_ptr,
                                    int mode, ADIO_Offset min_offset,
                                    ADIO_Offset max_offset);
void ADIOI_LUSTRE_Calc_my_req(ADIO_File fd, ADIO_Offset *offset_list,
                              ADIO_Offset *len_list, int contig_access_count,
                              ADIO_Offset *striping_info,
                              void *buf, int buftype_is_contig,
                              MPI_Datatype datatype, int count,
                              struct req_info **req_info_ptr,
                              int *req_info_count);

int ADIOI_LUSTRE_Calc_aggregator(ADIO_File fd, ADIO_Offset off,
                                 ADIO_Offset *len, ADIO_Offset *striping_info);
void ADIOI_LUSTRE_Calc_stripe(ADIO_Offset off, ADIO_Offset *striping_info,
                              int *stripe_num, int *rel_offset);
void ADIOI_LUSTRE_IwriteContig(ADIO_File fd, void *buf, int count,
                               MPI_Datatype datatype, int file_ptr_type,
                               ADIO_Offset offset, ADIO_Request *request,
                               int *error_code);
void ADIOI_LUSTRE_IreadContig(ADIO_File fd, void *buf, int count,
                              MPI_Datatype datatype, int file_ptr_type,
                              ADIO_Offset offset, ADIO_Request *request,
                              int *error_code);

#endif /* End of AD_UNIX_INCLUDE */
