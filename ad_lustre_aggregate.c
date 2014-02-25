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

#undef AGG_DEBUG

void ADIOI_LUSTRE_Get_striping_info(ADIO_File fd,
                                    ADIO_Offset **striping_info_ptr,
                                    int mode, ADIO_Offset min_offset,
                                    ADIO_Offset max_offset)
{
	ADIO_Offset *striping_info = NULL;
	/* get striping information:
	 *  striping_info[0]: stripe_size
	 *  striping_info[1]: stripe_count
	 *  striping_info[2]: avail_cb_nodes
	 *  striping_info[3]: min_aligned_offset
	 *  striping_info[4]: region_size - size of aar for one group
	 *  striping_info[5]: min_offset_rank
	 *  striping_info[6]: max_offset
	 *  striping_info[7..]: cb node rank
	 */
	ADIO_Offset aligned_aar_size, aar_stripes, min_stripe_offset;
	int stripe_size, stripe_count, CO = 1, num_groups, i;
	int avail_cb_nodes, divisor, nprocs_for_coll = fd->hints->cb_nodes;

	/* Get hints value */
	/* stripe size */
	stripe_size = fd->hints->striping_unit;
	/* stripe count */
	/* stripe_size and stripe_count have been validated in ADIOI_LUSTRE_Open() */
	stripe_count = fd->hints->striping_factor;

	/* Calculate the available number of I/O clients */
	if (!mode) {
		/* for collective read,
		 * if "CO" clients access the same OST simultaneously,
		 * the OST disk seek time would be much. So, to avoid this,
		 * it might be better if 1 client only accesses 1 OST.
		 * So, we set CO = 1 to meet the above requirement.
		 */
		CO = 1;
		/*XXX: maybe there are other better way for collective read */
	} else {
		/* CO also has been validated in ADIOI_LUSTRE_Open(), >0 */
		CO = fd->hints->fs_hints.lustre.co_ratio;
	}
	/* Calculate how many IO clients we need */
	/* Algorithm courtesy Pascal Deveze (pascal.deveze@bull.net) */
	/* To avoid extent lock conflicts,
	 * avail_cb_nodes should either
	 * 	- be a multiple of stripe_count,
	 *  - or divide stripe_count exactly
	 * so that each OST is accessed by a maximum of CO constant clients. */
	if (nprocs_for_coll >= stripe_count)
		/* avail_cb_nodes should be a multiple of stripe_count and the number
		 * of procs per OST should be limited to the minimum between
		 * nprocs_for_coll/stripe_count and CO
		 *
		 * e.g. if stripe_count=20, nprocs_for_coll=42 and CO=3 then
		 * avail_cb_nodes should be equal to 40 */
		avail_cb_nodes =
			stripe_count * ADIOI_MIN(nprocs_for_coll/stripe_count, CO);
	else {
		/* nprocs_for_coll is less than stripe_count */
		/* avail_cb_nodes should divide stripe_count */
		/* e.g. if stripe_count=60 and nprocs_for_coll=8 then
		 * avail_cb_nodes should be egal to 6 */
		/* This could be done with :
		   while (stripe_count % avail_cb_nodes != 0) avail_cb_nodes--;
		   but this can be optimized for large values of nprocs_for_coll and
		   stripe_count */
		divisor = 2;
		avail_cb_nodes = 1;
		/* try to divise */
		while (stripe_count >= divisor*divisor) {
			if ((stripe_count % divisor) == 0) {
				if (stripe_count/divisor <= nprocs_for_coll) {
					/* The value is found ! */
					avail_cb_nodes = stripe_count/divisor;
					break;
				}
				/* if divisor is less than nprocs_for_coll, divisor is a
				 * solution, but it is not sure that it is the best one */
				else if (divisor <= nprocs_for_coll)
					avail_cb_nodes = divisor;
			}
			divisor++;
		}
	}

	*striping_info_ptr = (ADIO_Offset *)
		ADIOI_Malloc((7 + avail_cb_nodes) * sizeof(ADIO_Offset));
	striping_info = *striping_info_ptr;
	striping_info[0] = stripe_size;
	striping_info[1] = stripe_count;
	striping_info[2] = avail_cb_nodes;
	min_stripe_offset = min_offset / stripe_size;
	striping_info[3] = min_stripe_offset * stripe_size;
	if (avail_cb_nodes % stripe_count == 0)
		num_groups = avail_cb_nodes / stripe_count;
	else
		num_groups = avail_cb_nodes / stripe_count + 1;
	if (max_offset % stripe_size == 0)
		aligned_aar_size = max_offset;
	else
		aligned_aar_size = (max_offset / stripe_size + 1) * stripe_size;
	aligned_aar_size -= striping_info[3];
	aar_stripes = aligned_aar_size / stripe_size;
	if (aar_stripes % num_groups == 0)
		striping_info[4] = aar_stripes / num_groups;
	else
		striping_info[4] = aar_stripes / num_groups + 1;
	striping_info[4] *= stripe_size;
	striping_info[5] = min_stripe_offset % avail_cb_nodes;
	striping_info[6] = max_offset;
	for (i = 0; i < avail_cb_nodes; ++i)
		striping_info[7 + i] = fd->hints->ranklist[i];
}

int ADIOI_LUSTRE_Calc_aggregator(ADIO_File fd, ADIO_Offset off,
                                 ADIO_Offset *len, ADIO_Offset *striping_info)
{
	int rank_index;
	ADIO_Offset avail_bytes;
	int stripe_size = striping_info[0];
	int stripe_count = striping_info[1];
	int avail_cb_nodes = striping_info[2];
	ADIO_Offset min_aligned_offset = striping_info[3];
	ADIO_Offset region_size = striping_info[4];
	int min_offset_rank = striping_info[5];

	/* Produce the group-cyclic pattern for Lustre */
	ADIO_Offset rel_aligned_offset = off - min_aligned_offset;
	rank_index =
		(min_offset_rank
		 + (rel_aligned_offset / region_size) * stripe_count
		 + (rel_aligned_offset / stripe_size) % stripe_count)
		% avail_cb_nodes;

	/* we index into fd_end with rank_index, and fd_end was allocated to be no
	 * bigger than fd->hints->cb_nodes.   If we ever violate that, we're
	 * overrunning arrays.  Obviously, we should never ever hit this abort
	 */
	if (rank_index >= fd->hints->cb_nodes)
		MPI_Abort(MPI_COMM_WORLD, 1);

	avail_bytes = (off / (ADIO_Offset)stripe_size + 1) *
		(ADIO_Offset)stripe_size - off;
	if (avail_bytes < *len) {
		/* this proc only has part of the requested contig. region */
		*len = avail_bytes;
	}
	/* map our index to a rank */
	/* NOTE: FOR NOW WE DON'T HAVE A MAPPING...JUST DO 0..NPROCS_FOR_COLL */
	/* return fd->hints->ranklist[rank_index]; */
	return rank_index;
}

void ADIOI_LUSTRE_Calc_stripe(ADIO_Offset off, ADIO_Offset *striping_info,
                              int *stripe_num, int *rel_offset)
{
	/* If the definition of ADIOI_Access could change to accommodate the
	 * stripe_num and rel_offset parameters calculated here, it might make sense
	 * to move this functionality into ADIOI_LUSTRE_Calc_aggregator.
	 */
	int stripe_size = striping_info[0];
	int stripe_count = striping_info[1];
	ADIO_Offset min_aligned_offset = striping_info[3];
	ADIO_Offset region_size = striping_info[4];

	/* Produce the group-cyclic pattern for Lustre */
	ADIO_Offset rel_aligned_offset = off - min_aligned_offset;

	/* compute the stripe number */
	*stripe_num =
		(rel_aligned_offset % region_size) / (stripe_size * stripe_count);

	/* compute the offset relative to the beginning of the stripe */
	*rel_offset = rel_aligned_offset % stripe_size;
}

/* ADIOI_LUSTRE_Calc_my_req() - calculate what portions of the access requests
 * of this process are located in the file domains of various processes
 * (including this one)
 */

static void init_user_buf(int buftype_is_contig, MPI_Aint buftype_extent,
                          int count, void *buf,
                          ADIOI_Flatlist_node *flat_buf, int *flat_buf_idx,
                          int *n_buftypes, void **user_buf,
                          ADIO_Offset *user_buf_sz)
{
	if (!buftype_is_contig) {
		*flat_buf_idx = 0;
		*n_buftypes = 0;
		*user_buf = buf + flat_buf->indices[0];
		*user_buf_sz = flat_buf->blocklens[0];
	} else {
		*user_buf = buf;
		*user_buf_sz = count * buftype_extent;
	}
}

#define NEXT_SEGMENT(len) do {                                          \
		user_buf_sz -= (len); \
		user_buf += (len); \
		if (!buftype_is_contig && !user_buf_sz) { \
			if (flat_buf_idx < (flat_buf->count - 1)) \
				flat_buf_idx++; \
			else { \
				flat_buf_idx = 0; \
				n_buftypes++; \
			} \
			user_buf = buf + flat_buf->indices[flat_buf_idx] + \
				(ADIO_Offset)n_buftypes * (ADIO_Offset)buftype_extent; \
			user_buf_sz = flat_buf->blocklens[flat_buf_idx]; \
		} \
		off += (len); \
		rem_len -= (len); \
	} while (0)

void ADIOI_LUSTRE_Calc_my_req(ADIO_File fd, ADIO_Offset *offset_list,
                              ADIO_Offset *len_list, int contig_access_count,
                              ADIO_Offset *striping_info,
                              void *buf, int buftype_is_contig,
                              MPI_Datatype datatype, int count,
                              struct req_info **req_info_ptr,
                              int *req_info_count)
{
	int i;
	ADIO_Offset avail_len, rem_len, off;
	MPI_Aint buftype_extent;
	ADIOI_Flatlist_node *flat_buf = NULL;
	int flat_buf_idx, n_buftypes;
	ADIO_Offset user_buf_sz;
	void *user_buf;
	struct req_info *req_info;

	MPI_Type_extent(datatype, &buftype_extent);
	if (!buftype_is_contig) {
		ADIOI_Flatten_datatype(datatype);
		flat_buf = ADIOI_Flatlist;
		while (flat_buf->type != datatype)
			flat_buf = flat_buf->next;
	}
	init_user_buf(buftype_is_contig, buftype_extent, count, buf,
	              flat_buf, &flat_buf_idx, &n_buftypes,
	              &user_buf, &user_buf_sz);

	/* one pass just to calculate how much space to allocate for my_req;
	 * contig_access_count was calculated way back in ADIOI_Calc_my_off_len()
	 */
	*req_info_count = 0;
	for (i = 0; i < contig_access_count; i++) {
		off = offset_list[i];
		rem_len = len_list[i];
		while (rem_len > 0) {
			avail_len = ADIOI_MIN(rem_len, user_buf_sz);
			ADIOI_LUSTRE_Calc_aggregator(fd, off, &avail_len, striping_info);
			*req_info_count += 1;
			NEXT_SEGMENT(avail_len);
		}
	}

	init_user_buf(buftype_is_contig, buftype_extent, count, buf,
	              flat_buf, &flat_buf_idx, &n_buftypes,
	              &user_buf, &user_buf_sz);

	*req_info_ptr = ADIOI_Malloc(*req_info_count * sizeof(struct req_info));
	req_info = *req_info_ptr;
	for (i = 0; i < contig_access_count; ++i) {
		off = offset_list[i];
		rem_len = len_list[i];
		while (rem_len > 0) {
			avail_len = ADIOI_MIN(rem_len, user_buf_sz);
			req_info->aggregator_index =
				ADIOI_LUSTRE_Calc_aggregator(fd, off, &avail_len, striping_info);
			req_info->length = avail_len;
			ADIOI_LUSTRE_Calc_stripe(off, striping_info,
			                         &(req_info->stripe_num),
			                         &(req_info->relative_offset));
			req_info->source = user_buf;
			NEXT_SEGMENT(req_info->length);
			req_info++;
		}
	}

	if (!buftype_is_contig)
		ADIOI_Delete_flattened(datatype);
}

int ADIOI_LUSTRE_Docollect(ADIO_File fd, int contig_access_count,
                           ADIO_Offset *len_list)
{
	/* If the processes are non-interleaved, we will check the req_size.
	 *   if (avg_req_size > big_req_size) {
	 *       docollect = 0;
	 *   }
	 */

	int i, docollect = 1, big_req_size = 0;
	ADIO_Offset req_size = 0, total_req_size;
	int avg_req_size, total_access_count;

	/* calculate total_req_size and total_access_count */
	for (i = 0; i < contig_access_count; i++)
		req_size += len_list[i];
	MPI_Allreduce(&req_size, &total_req_size, 1, MPI_LONG_LONG_INT, MPI_SUM,
	              fd->comm);
	MPI_Allreduce(&contig_access_count, &total_access_count, 1, MPI_INT, MPI_SUM,
	              fd->comm);
	/* estimate average req_size */
	avg_req_size = (int)(total_req_size / total_access_count);
	/* get hint of big_req_size */
	big_req_size = fd->hints->fs_hints.lustre.coll_threshold;
	/* Don't perform collective I/O if there are big requests */
	if ((big_req_size > 0) && (avg_req_size > big_req_size))
		docollect = 0;

	return docollect;
}
