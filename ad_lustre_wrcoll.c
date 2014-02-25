/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   Copyright (C) 1997 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 *
 *   Copyright (C) 2007 Oak Ridge National Laboratory
 *
 *   Copyright (C) 2008 Sun Microsystems, Lustre group
 */

#include "ad_lustre.h"
#include "adio_extern.h"
#include <shmem.h>
#include <stdio.h>

static void ADIOI_LUSTRE_Exch_and_write(ADIO_File fd, int mype,
                                        struct req_info *req_info,
                                        int req_info_count,
                                        ADIO_Offset *striping_info,
                                        int num_buffered_stripes,
                                        void **collective_buffers,
                                        long **locks, int **mod_flags,
                                        int *error_code);

void ADIOI_LUSTRE_WriteStridedColl(ADIO_File fd, void *buf, int count,
                                   MPI_Datatype datatype,
                                   int file_ptr_type, ADIO_Offset offset,
                                   ADIO_Status *status, int *error_code)
{
	/* Uses a generalized version of the extended two-phase method described
	 * in "An Extended Two-Phase Method for Accessing Sections of
	 * Out-of-Core Arrays", Rajeev Thakur and Alok Choudhary,
	 * Scientific Programming, (5)4:301--317, Winter 1996.
	 * http://www.mcs.anl.gov/home/thakur/ext2ph.ps
	 */
	printf("in my code!!!");
	int i, filetype_is_contig, nprocs, myrank, do_collect = 0;
	int contig_access_count = 0, buftype_is_contig, interleave_count = 0;
	ADIO_Offset orig_fp, start_offset, end_offset, off, min_offset, max_offset;
	ADIO_Offset *offset_list = NULL, *st_offsets = NULL, *end_offsets = NULL;
	ADIO_Offset *len_list = NULL;
	ADIO_Offset *striping_info = NULL;
	int old_error, tmp_error;
	char *value;
	int info_flag, coll_bufsize, num_buffered_stripes;
	int avail_cb_nodes;
	void **collective_buffers;
	long **locks = NULL;
	int **mod_flags;
	struct req_info *req_info;
	int req_info_count;

	MPI_Comm_size(fd->comm, &nprocs);
	MPI_Comm_rank(fd->comm, &myrank);

	orig_fp = fd->fp_ind;

	/* IO patten identification if cb_write isn't disabled */
	if (fd->hints->cb_write != ADIOI_HINT_DISABLE) {
		/* For this process's request, calculate the list of offsets and
		   lengths in the file and determine the start and end offsets. */

		/* Note: end_offset points to the last byte-offset that will be accessed.
		 * e.g., if start_offset=0 and 100 bytes to be read, end_offset=99
		 */

		ADIOI_Calc_my_off_len(fd, count, datatype, file_ptr_type, offset,
		                      &offset_list, &len_list, &start_offset,
		                      &end_offset, &contig_access_count);

		/* each process communicates its start and end offsets to other
		 * processes. The result is an array each of start and end offsets
		 * stored in order of process rank.
		 */
		st_offsets = (ADIO_Offset *) ADIOI_Malloc(nprocs * sizeof(ADIO_Offset));
		end_offsets = (ADIO_Offset *) ADIOI_Malloc(nprocs * sizeof(ADIO_Offset));
		MPI_Allgather(&start_offset, 1, ADIO_OFFSET, st_offsets, 1,
		              ADIO_OFFSET, fd->comm);
		MPI_Allgather(&end_offset, 1, ADIO_OFFSET, end_offsets, 1,
		              ADIO_OFFSET, fd->comm);
		min_offset = st_offsets[0];
		max_offset = end_offsets[0];
		/* are the accesses of different processes interleaved? */
		for (i = 1; i < nprocs; i++) {
			if (st_offsets[i] < min_offset)
				min_offset = st_offsets[i];
			if (end_offsets[i] > max_offset)
				max_offset = end_offsets[i];
			if ((st_offsets[i] < end_offsets[i-1]) &&
			    (st_offsets[i] <= end_offsets[i]))
				interleave_count++;
		}
		/* This is a rudimentary check for interleaving, but should suffice
		   for the moment. */

		/* Two typical access patterns can benefit from collective write.
		 *   1) the processes are interleaved, and
		 *   2) the req size is small.
		 */
		if (interleave_count > 0)
			do_collect = 1;
		else
			do_collect =
				ADIOI_LUSTRE_Docollect(fd, contig_access_count, len_list);
	}
	ADIOI_Datatype_iscontig(datatype, &buftype_is_contig);

	/* Decide if collective I/O should be done */
	if ((!do_collect && fd->hints->cb_write == ADIOI_HINT_AUTO) ||
	    fd->hints->cb_write == ADIOI_HINT_DISABLE) {

		/* use independent accesses */
		if (fd->hints->cb_write != ADIOI_HINT_DISABLE) {
			ADIOI_Free(offset_list);
			ADIOI_Free(len_list);
			ADIOI_Free(st_offsets);
			ADIOI_Free(end_offsets);
		}

		fd->fp_ind = orig_fp;
		ADIOI_Datatype_iscontig(fd->filetype, &filetype_is_contig);
		if (buftype_is_contig && filetype_is_contig) {
			if (file_ptr_type == ADIO_EXPLICIT_OFFSET) {
				off = fd->disp + (ADIO_Offset)(fd->etype_size) * offset;
				ADIO_WriteContig(fd, buf, count, datatype,
				                 ADIO_EXPLICIT_OFFSET,
				                 off, status, error_code);
			} else
				ADIO_WriteContig(fd, buf, count, datatype, ADIO_INDIVIDUAL,
				                 0, status, error_code);
		} else {
			ADIO_WriteStrided(fd, buf, count, datatype, file_ptr_type,
			                  offset, status, error_code);
		}
		return;
	}

	/* get collective buffer size */
	value = (char *) ADIOI_Malloc((MPI_MAX_INFO_VAL+1)*sizeof(char));
	ADIOI_Info_get(fd->info, "cb_buffer_size", MPI_MAX_INFO_VAL, value,
	               &info_flag);
	coll_bufsize = atoi(value);
	ADIOI_Free(value);

	/* Get Lustre hints information */
	ADIOI_LUSTRE_Get_striping_info(fd, &striping_info, 1, min_offset,
	                               max_offset);
	avail_cb_nodes = striping_info[2];

	num_buffered_stripes = coll_bufsize / striping_info[0];
	if (num_buffered_stripes == 0)
		num_buffered_stripes = 1;
	
	/* Allocate memory in the symmetric heap for all collective buffers and
	 * locks (if needed) */
	collective_buffers = ADIOI_Malloc(avail_cb_nodes * sizeof(void *));
	for (i = 0; i < avail_cb_nodes; ++i)
		collective_buffers[i] =
			LUSTRE_shmalloc(num_buffered_stripes * striping_info[0]);
	if (fd->atomicity) {
		locks = ADIOI_Malloc(avail_cb_nodes * sizeof(long *));
		for (i = 0; i < avail_cb_nodes; ++i) {
			locks[i] = LUSTRE_shmalloc(num_buffered_stripes * sizeof(long));
			memset(locks[i], 0, num_buffered_stripes * sizeof(long));
		}
	}
	mod_flags = ADIOI_Malloc(avail_cb_nodes * sizeof(int *));
	for (i = 0; i < avail_cb_nodes; ++i)
		mod_flags = LUSTRE_shmalloc(num_buffered_stripes * sizeof(int));

	/* calculate what portions of the access requests of this process are
	 * located in which process
	 */
	ADIOI_LUSTRE_Calc_my_req(fd, offset_list, len_list, contig_access_count,
	                         striping_info, buf, buftype_is_contig, datatype,
	                         count, &req_info, &req_info_count);

	/* exchange data and write in sizes of no more than stripe_size. */
	ADIOI_LUSTRE_Exch_and_write(fd, myrank, req_info, req_info_count,
	                            striping_info, num_buffered_stripes,
	                            collective_buffers, locks, mod_flags,
	                            error_code);

	/* If this collective write is followed by an independent write,
	 * it's possible to have those subsequent writes on other processes
	 * race ahead and sneak in before the read-modify-write completes.
	 * We carry out a collective communication at the end here so no one
	 * can start independent i/o before collective I/O completes.
	 *
	 * need to do some gymnastics with the error codes so that if something
	 * went wrong, all processes report error, but if a process has a more
	 * specific error code, we can still have that process report the
	 * additional information */

	old_error = *error_code;
	if (*error_code != MPI_SUCCESS)
		*error_code = MPI_ERR_IO;

	/* optimization: if only one process performing i/o, we can perform
	 * a less-expensive Bcast  */
#ifdef ADIOI_MPE_LOGGING
	MPE_Log_event(ADIOI_MPE_postwrite_a, 0, NULL);
#endif
	if (fd->hints->cb_nodes == 1)
		MPI_Bcast(error_code, 1, MPI_INT,
		          fd->hints->ranklist[0], fd->comm);
	else {
		tmp_error = *error_code;
		MPI_Allreduce(&tmp_error, error_code, 1, MPI_INT,
		              MPI_MAX, fd->comm);
	}
#ifdef ADIOI_MPE_LOGGING
	MPE_Log_event(ADIOI_MPE_postwrite_b, 0, NULL);
#endif

	if ((old_error != MPI_SUCCESS) && (old_error != MPI_ERR_IO))
		*error_code = old_error;

	/* free all memory allocated for collective I/O */

	/* Free memory in the symmetric heap for all collective buffers and
	 * locks (if allocated) */
	if (fd->atomicity) {
		for (i = 0; i < avail_cb_nodes; ++i)
			LUSTRE_shfree(locks[i]);
		ADIOI_Free(locks);
	}
	for (i = 0; i < avail_cb_nodes; ++i)
		LUSTRE_shfree(collective_buffers[i]);
	ADIOI_Free(collective_buffers);
	for (i = 0; i < avail_cb_nodes; ++i)
		LUSTRE_shfree(mod_flags[i]);
	ADIOI_Free(mod_flags);

	ADIOI_Free(req_info);

	ADIOI_Free(offset_list);
	ADIOI_Free(len_list);
	ADIOI_Free(st_offsets);
	ADIOI_Free(end_offsets);
	ADIOI_Free(striping_info);

#ifdef HAVE_STATUS_SET_BYTES
	if (status) {
		int bufsize, size;
		/* Don't set status if it isn't needed */
		MPI_Type_size(datatype, &size);
		bufsize = size * count;
		MPIR_Status_set_bytes(status, datatype, bufsize);
	}
	/* This is a temporary way of filling in status. The right way is to
	 * keep track of how much data was actually written during collective I/O.
	 */
#endif

	fd->fp_sys_posn = -1;	/* set it to null. */
}

/* If successful, error_code is set to MPI_SUCCESS.  Otherwise an error
 * code is created and returned in error_code.
 */
static void ADIOI_LUSTRE_Exch_and_write(ADIO_File fd, int mype,
                                        struct req_info *req_info,
                                        int req_info_count,
                                        ADIO_Offset *striping_info,
                                        int num_buffered_stripes,
                                        void **collective_buffers,
                                        long **locks, int **mod_flags,
                                        int *error_code)
{
	int i, j, err;
	int stripe_count = striping_info[1], avail_cb_nodes = striping_info[2];
	int min_offset_rank = striping_info[5];
	ADIO_Offset stripe_size = striping_info[0], region_size = striping_info[4];
	ADIO_Offset min_aligned_offset = striping_info[3], *pes = &striping_info[7];
	ADIO_Offset group_step_size, block_offset, off;
	int group_size, max_stripe_num, stripe_num, num_groups, group;
	int req_stripe, req_agg, current_stripe, current_agg, *group_req_count;
	struct req_info **group_req_info, *ri;
	void *coll_buff;
#ifdef WCBE_AUTOZERO_WORKAROUND
	MPI_Status status;
	int count;
#endif

	*error_code = MPI_SUCCESS;	/* changed below if error */
	/* only I/O errors are currently reported */

	group_size = ADIOI_MIN(avail_cb_nodes, stripe_count);
	num_groups = avail_cb_nodes / group_size;
	group_step_size = (ADIO_Offset) group_size * stripe_size;

	/* Within each group, stripe numbers are in non-decreasing order. We can use
	 * this property to order the traversal through stripe numbers in the
	 * request list.
	 */
	max_stripe_num = 0;
	group_req_info = ADIOI_Calloc(num_groups, sizeof(struct req_info *));
	group_req_count = ADIOI_Calloc(num_groups, sizeof(int));
	for (i = 0; i < req_info_count; ++i) {
		max_stripe_num = ADIOI_MAX(max_stripe_num, req_info[i].stripe_num);
		group =
			((req_info[i].aggregator_index + avail_cb_nodes - min_offset_rank)
			 % avail_cb_nodes) / num_groups;
		if (!group_req_info[group])
			group_req_info[group] = &req_info[i];
		group_req_count[group]++;
	}
	MPI_Allreduce(MPI_IN_PLACE, &max_stripe_num, 1, MPI_INT, MPI_MAX, fd->comm);

	block_offset = -1;
	for (stripe_num = 0; stripe_num <= max_stripe_num; stripe_num += num_buffered_stripes) {
		// Note that, without knowledge of where the buffer is written, we can
		// never avoid doing read-modify-write cycles.
		for (i = 0; i < avail_cb_nodes; i++) {
			if (mype == pes[i]) {
				if (block_offset == -1)
					block_offset = min_aligned_offset
						+ (i / group_size) * region_size
						+ (i % group_size) * stripe_size;
				off = block_offset;
				coll_buff = collective_buffers[i];
				for (j = 0; j < num_buffered_stripes; j++) {
#ifndef WCBE_AUTOZERO_WORKAROUND
					ADIO_ReadContig(fd, coll_buff, stripe_size,
					                MPI_BYTE, ADIO_EXPLICIT_OFFSET,
					                off, MPI_STATUS_IGNORE, &err);
#else
					ADIO_ReadContig(fd, coll_buff, stripe_size,
					                MPI_BYTE, ADIO_EXPLICIT_OFFSET,
					                off, &status, &err);
					MPI_Get_count(&status, MPI_BYTE, &count);
					if (count < stripe_size)
						memset(coll_buff + count, 0, stripe_size - count);
#endif
					if (err != MPI_SUCCESS) {
						MPI_Abort(MPI_COMM_WORLD, 1); // FIXME
						/* *error_code = */
						/* 	MPIO_Err_create_code(err, */
						/* 	                     MPIR_ERR_RECOVERABLE, */
						/* 	                     myname, __LINE__, */
						/* 	                     MPI_ERR_IO, */
						/* 	                     "**ioRMWrdwr", 0); */
						/* return; */
					}
					coll_buff += stripe_size;
					off += group_step_size;
				}
				memset(mod_flags[i], 0, num_buffered_stripes * sizeof(int));
			}
		}
		/* Leave the following barrier outside the previous block since it's
		 * also needed to ensure that any preceding write has been completed. */
		shmem_barrier_all();
		for (group = 0; group < num_groups; group++) {
			current_stripe = -1;
			current_agg = -1;
			while (group_req_count[group] > 0) {
				ri = group_req_info[group];
				req_stripe = ri->stripe_num;
				if (req_stripe >= stripe_num + num_buffered_stripes)
					break;
				req_stripe %= num_buffered_stripes;
				req_agg = ri->aggregator_index;
				if (current_stripe != req_stripe || current_agg != req_agg) {
					if (fd->atomicity) {
						if (current_stripe != -1)
							shmem_clear_lock(&(locks[current_agg][current_stripe]));
						shmem_set_lock(&(locks[req_agg][req_stripe]));
					}
					shmem_int_inc(&(mod_flags[req_agg][req_stripe]), pes[req_agg]);
					coll_buff =
						collective_buffers[req_agg] + req_stripe * stripe_size;
					current_stripe = req_stripe;
					current_agg = req_agg;
				}
				// copy data for ri to "target" on PE pes[agg]
				shmem_putmem(coll_buff + ri->relative_offset, ri->source,
				               ri->length, pes[req_agg]);
				group_req_info[group]++;
				group_req_count[group]--;
			}
			if (fd->atomicity && current_stripe != -1)
				shmem_clear_lock(&(locks[current_agg][current_stripe]));
		}
		shmem_barrier_all();
		// write collective buffers to file
		for (i = 0; i < avail_cb_nodes; i++) {
			if (mype == pes[i]) {
				off = block_offset;
				coll_buff = collective_buffers[i];
				for (j = 0; j < num_buffered_stripes; j++) {
					if (mod_flags[i][j] > 0) {
						ADIO_WriteContig(fd, coll_buff, stripe_size,
						                 MPI_BYTE, ADIO_EXPLICIT_OFFSET,
						                 off, MPI_STATUS_IGNORE, &err);
						if (err != MPI_SUCCESS) {
							MPI_Abort(MPI_COMM_WORLD, 1); // FIXME
							/* *error_code = */
							/* 	MPIO_Err_create_code(err, */
							/* 	                     MPIR_ERR_RECOVERABLE, */
							/* 	                     myname, __LINE__, */
							/* 	                     MPI_ERR_IO, */
							/* 	                     "**ioRMWrdwr", 0); */
							/* return; */
						}
					}
					coll_buff += stripe_size;
					off += group_step_size;
				}
				block_offset = off;
			}
		}
	}
	ADIOI_Free(group_req_info);
	ADIOI_Free(group_req_count);
}
