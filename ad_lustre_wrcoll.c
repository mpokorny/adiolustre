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

#define USE_IWRITE_FUNCTIONS // TODO: remove this def
#undef USE_WAITSOME // TODO: remove this def

struct ExchAndWrite_RequestSeq {
  /* The following values are provided upon instantiation */
  ADIO_File fd;
  void *buf;
  ADIO_Offset *offset_list, *len_list, *striping_info, *stripe_off_list;
  ADIOI_Access *others_req, *my_req;
  ADIOI_Flatlist_node *flat_buf;
  int **buf_idx;
  int nprocs, myrank;
  int buftype_is_contig, contig_access_count;
  MPI_Aint buftype_extent;
  ADIO_Offset end_loc;
  ADIO_Offset group_step_size;
  int iter_step;

  /* The following belong to this "instance", and may be freely
   * modified by ExchAndWrite_RequestSeq functions. */
  char *write_buf;
  /* recv_count: to store count of how many off-len pairs per proc are
     satisfied in an iteration;

     send_size: total size of data to be sent to each proc. in an
     iteration. Of size nprocs so that I can use MPI_Alltoall later;

     recv_size: total size of data to be recd. from each proc. in an
     iteration;

     sent_to_proc: amount of data sent to each proc so far. Used in
     ADIOI_Fill_send_buffer;

     send_buf_idx, curr_to_proc, done_to_proc: are used in
     ADIOI_Fill_send_buffer;

     recv_start_pos: used to store the starting value of
     recv_curr_offlen_ptr[i] in an iteration
 */
  int *recv_count, *recv_size, *send_size;
  int *recv_start_pos;
  int *this_buf_idx, *send_buf_idx, *curr_to_proc, *done_to_proc;
  ADIO_Offset *srt_off;
  int *srt_len, srt_num;
  int num_sub_requests;
  MPI_Request *sub_requests;
  char **send_buf;

  int iter;
  int real_size;
  int phase;
  ADIO_Offset off;
  int hole;
  int error_code;
};
#define ITER_N(offset, si, gs) ((((offset) - (si)[3]) % (si)[4]) / (gs))

/* prototypes of functions used for collective writes only. */
static void ADIOI_LUSTRE_Exch_and_write(ADIO_File fd, void *buf,
					MPI_Datatype datatype, int nprocs,
					int myrank,
					ADIOI_Access *others_req,
					ADIOI_Access *my_req,
					ADIO_Offset *offset_list,
					ADIO_Offset *len_list,
					int contig_access_count,
					ADIO_Offset *striping_info,
                                        int **buf_idx, int *error_code);
static void ADIOI_LUSTRE_Fill_send_buffer(ADIO_File fd, void *buf,
					  ADIOI_Flatlist_node *flat_buf,
					  char **send_buf,
					  ADIO_Offset *offset_list,
					  ADIO_Offset *len_list, int *send_size,
					  MPI_Request *requests,
					  int *sent_to_proc, int nprocs,
					  int myrank, int contig_access_count,
					  ADIO_Offset *striping_info,
					  int *send_buf_idx, int *curr_to_proc,
					  int *done_to_proc, int iter,
					  MPI_Aint buftype_extent);
static void ADIOI_LUSTRE_W_Exchange_data(ADIO_File fd, void *buf,
					 char *write_buf,
					 ADIOI_Flatlist_node *flat_buf,
					 ADIO_Offset *offset_list,
					 ADIO_Offset *len_list, int *send_size,
					 int *recv_size, ADIO_Offset off,
					 int size, int *count,
					 int *start_pos,
					 int *sent_to_proc, int nprocs,
					 int myrank, int buftype_is_contig,
					 int contig_access_count,
					 ADIO_Offset *striping_info,
					 ADIOI_Access *others_req,
					 int *send_buf_idx, int *curr_to_proc,
					 int *done_to_proc, int *hole, int iter,
                          MPI_Aint buftype_extent, int *buf_idx,
					 ADIO_Offset **srt_off, int **srt_len, int *srt_num,
                          MPI_Request *requests, int *num_requests,
                          char **send_buf, int *error_code);
void ADIOI_Heap_merge(ADIOI_Access *others_req, int *count,
                      ADIO_Offset *srt_off, int *srt_len, int *start_pos,
                      int nprocs, int nprocs_recv, int total_elements);

static int ExchAndWrite_RequestSeq_test(struct ExchAndWrite_RequestSeq *rseq,
                                        int *flag);
static int ExchAndWrite_RequestSeq_poll(struct ExchAndWrite_RequestSeq *rseq);
static int ExchAndWrite_RequestSeq_wait(int count,
                                        struct ExchAndWrite_RequestSeq **rseqs);
#ifdef USE_WAITSOME
static int ExchAndWrite_RequestSeq_waitsome(int count,
                                            struct ExchAndWrite_RequestSeq **rseqs, int *outcount);
#endif
static void init_exch_and_write(ADIO_File fd, void *buf,
                                ADIO_Offset *offset_list,
                                ADIO_Offset *len_list,
                                ADIO_Offset *striping_info,
                                ADIO_Offset *stripe_off_list,
                                ADIO_Offset end_loc,
                                ADIOI_Access *others_req, ADIOI_Access *my_req,
                                ADIOI_Flatlist_node *flat_buf, int **buf_idx,
                                int nprocs, int myrank,
                                int buftype_is_contig, int contig_access_count,
                                MPI_Aint buftype_extent,
                                int iter_offset, int iter_step,
                                struct ExchAndWrite_RequestSeq *rseq);
static void finalize_exch_and_write(struct ExchAndWrite_RequestSeq *rseq);
static int start_next_exch_and_write(struct ExchAndWrite_RequestSeq *rseq,
                                     int *sent_to_proc,
                                     int *send_curr_offlen_ptr,
                                     int *recv_curr_offlen_ptr);
static int poll_exch(struct ExchAndWrite_RequestSeq *rseq);
static int poll_write(struct ExchAndWrite_RequestSeq *rseq);

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

    ADIOI_Access *my_req;
    /* array of nprocs access structures, one for each other process has
       this process's request */

    ADIOI_Access *others_req;
    /* array of nprocs access structures, one for each other process
       whose request is written by this process. */

    int i, filetype_is_contig, nprocs, myrank, do_collect = 0;
    int contig_access_count = 0, buftype_is_contig, interleave_count = 0;
    int *count_my_req_per_proc, count_my_req_procs, count_others_req_procs;
    ADIO_Offset orig_fp, start_offset, end_offset, off, min_offset, max_offset;
    ADIO_Offset *offset_list = NULL, *st_offsets = NULL, *end_offsets = NULL;
    ADIO_Offset *len_list = NULL;
    int **buf_idx = NULL;
    ADIO_Offset *striping_info = NULL;
    int old_error, tmp_error;

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
        if (interleave_count > 0) {
	    do_collect = 1;
        } else {
            do_collect = ADIOI_LUSTRE_Docollect(fd, contig_access_count,
			                        len_list, nprocs);
        }
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

    /* Get Lustre hints information */
    ADIOI_LUSTRE_Get_striping_info(fd, &striping_info, 1, min_offset,
                                   max_offset);

    /* calculate what portions of the access requests of this process are
     * located in which process
     */
    ADIOI_LUSTRE_Calc_my_req(fd, offset_list, len_list, contig_access_count,
                             striping_info, nprocs,
                             &count_my_req_procs,
                             &count_my_req_per_proc, &my_req,
                             &buf_idx);

    /* based on everyone's my_req, calculate what requests of other processes
     * will be accessed by this process.
     * count_others_req_procs = number of processes whose requests (including
     * this process itself) will be accessed by this process
     * count_others_req_per_proc[i] indicates how many separate contiguous
     * requests of proc. i will be accessed by this process.
     */

    ADIOI_Calc_others_req(fd, count_my_req_procs, count_my_req_per_proc,
                          my_req, nprocs, myrank, &count_others_req_procs,
                          &others_req);
    ADIOI_Free(count_my_req_per_proc);

    /* exchange data and write in sizes of no more than stripe_size. */
    ADIOI_LUSTRE_Exch_and_write(fd, buf, datatype, nprocs, myrank,
                                others_req, my_req, offset_list, len_list,
                                contig_access_count, striping_info,
                                buf_idx, error_code);

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


    if (!buftype_is_contig)
	ADIOI_Delete_flattened(datatype);

    /* free all memory allocated for collective I/O */
    /* free others_req */
    for (i = 0; i < nprocs; i++) {
	if (others_req[i].count) {
	    ADIOI_Free(others_req[i].offsets);
	    ADIOI_Free(others_req[i].lens);
	    ADIOI_Free(others_req[i].mem_ptrs);
	}
    }
    ADIOI_Free(others_req);
    /* free my_req here */
    for (i = 0; i < nprocs; i++) {
	if (my_req[i].count) {
	    ADIOI_Free(my_req[i].offsets);
	    ADIOI_Free(my_req[i].lens);
	}
    }
    ADIOI_Free(my_req);
    for (i = 0; i < nprocs; i++) {
        ADIOI_Free(buf_idx[i]);
    }
    ADIOI_Free(buf_idx);
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
static void ADIOI_LUSTRE_Exch_and_write(ADIO_File fd, void *buf,
					MPI_Datatype datatype, int nprocs,
					int myrank, ADIOI_Access *others_req,
                                        ADIOI_Access *my_req,
					ADIO_Offset *offset_list,
                                        ADIO_Offset *len_list,
					int contig_access_count,
                         ADIO_Offset *striping_info, int **buf_idx,
                                        int *error_code)
{
    /* Send data to appropriate processes and write in sizes of no more
     * than lustre stripe_size.
     * The idea is to reduce the amount of extra memory required for
     * collective I/O. If all data were written all at once, which is much
     * easier, it would require temp space more than the size of user_buf,
     * which is often unacceptable. For example, to write a distributed
     * array to a file, where each local array is 8Mbytes, requiring
     * at least another 8Mbytes of temp space is unacceptable.
     */

    int hole, i, j, m, flag, max_ntimes, buftype_is_contig;
    ADIO_Offset st_loc = -1, end_loc = -1, min_st_loc, max_end_loc;
    ADIO_Offset off, req_off, send_off, *off_list;
    ADIO_Offset max_size, step_size = 0, group_step_size;
    int group_size;
    int req_len, send_len;
    int *send_curr_offlen_ptr, *recv_curr_offlen_ptr, *sent_to_proc;
    ADIOI_Flatlist_node *flat_buf = NULL;
    MPI_Aint buftype_extent;
    int stripe_count = striping_info[1], avail_cb_nodes = striping_info[2];
    ADIO_Offset stripe_size = striping_info[0], region_size = striping_info[4];
    int data_sieving = 0;
    ADIO_Offset block_offset;
    int block_len;
    char *value;
    int info_flag, coll_bufsize, num_rseqs, rseq_idx;

    struct ExchAndWrite_RequestSeq *rseqs, **active_rseqs;
    int num_active;

    *error_code = MPI_SUCCESS;	/* changed below if error */
    /* only I/O errors are currently reported */

    value = (char *) ADIOI_Malloc((MPI_MAX_INFO_VAL+1)*sizeof(char));
    ADIOI_Info_get(fd->info, "cb_buffer_size", MPI_MAX_INFO_VAL, value,
                   &info_flag);
    coll_bufsize = atoi(value);
    ADIOI_Free(value);

    /* calculate the number of writes of stripe size to be done. */

    for (i = 0; i < nprocs; i++) {
      if (others_req[i].count) {
        st_loc = others_req[i].offsets[0];
        end_loc = others_req[i].offsets[0];
        break;
      }
    }
    for (i = 0; i < nprocs; i++) {
      for (j = 0; j < others_req[i].count; j++) {
        st_loc = ADIOI_MIN(st_loc, others_req[i].offsets[j]);
        end_loc = ADIOI_MAX(end_loc, (others_req[i].offsets[j] +
                                      others_req[i].lens[j] - 1));
      }
    }
    max_end_loc = striping_info[6];
    min_st_loc = striping_info[3];

    /* Each time, only avail_cb_nodes number of IO clients perform IO,
     * so, step_size=avail_cb_nodes*stripe_size IO will be performed at most,
     * and ntimes=whole_file_portion/step_size
     */
    step_size = (ADIO_Offset) avail_cb_nodes * stripe_size;
    group_size = ADIOI_MIN(avail_cb_nodes, stripe_count);
    group_step_size = (ADIO_Offset) group_size * stripe_size;
    max_ntimes = (max_end_loc - min_st_loc + 1) / step_size
      + (((max_end_loc - min_st_loc + 1) % step_size) ? 1 : 0);
    num_rseqs = ADIOI_MAX(coll_bufsize / stripe_size, 1);
    rseqs = (struct ExchAndWrite_RequestSeq *) ADIOI_Malloc(num_rseqs * sizeof(struct ExchAndWrite_RequestSeq));
    active_rseqs = (struct ExchAndWrite_RequestSeq **) ADIOI_Malloc(num_rseqs * sizeof(struct ExchAndWrite_RequestSeq *));
    recv_curr_offlen_ptr = (int *) ADIOI_Calloc(nprocs, sizeof(int));
    send_curr_offlen_ptr = (int *) ADIOI_Calloc(nprocs, sizeof(int));
    sent_to_proc = (int *) ADIOI_Malloc(nprocs * sizeof(int));

    /* calculate the start offset for each iteration */
    off_list = (ADIO_Offset *) ADIOI_Malloc(max_ntimes * sizeof(ADIO_Offset));
    for (m = 0; m < max_ntimes; m ++)
      off_list[m] = max_end_loc;
    for (i = 0; i < nprocs; i++) {
      for (j = 0; j < others_req[i].count; j ++) {
        req_off = others_req[i].offsets[j];
        m = ((req_off - min_st_loc) % region_size) / group_step_size;
        off_list[m] = ADIOI_MIN(off_list[m], req_off);
      }
    }

    ADIOI_Datatype_iscontig(datatype, &buftype_is_contig);
    if (!buftype_is_contig) {
      ADIOI_Flatten_datatype(datatype);
      flat_buf = ADIOI_Flatlist;
      while (flat_buf->type != datatype)
        flat_buf = flat_buf->next;
    }
    MPI_Type_extent(datatype, &buftype_extent);

    for (i = 0; i < num_rseqs; ++i)
      init_exch_and_write(fd, buf, offset_list, len_list, striping_info,
                          off_list, end_loc, others_req, my_req, flat_buf,
                          buf_idx, nprocs, myrank, buftype_is_contig,
                          contig_access_count, buftype_extent,
                          i, num_rseqs, &(rseqs[i]));

    /* I need to check if there are any outstanding nonblocking writes to
     * the file, which could potentially interfere with the writes taking
     * place in this collective write call. Since this is not likely to be
     * common, let me do the simplest thing possible here: Each process
     * completes all pending nonblocking operations before completing.
     */
    /*ADIOI_Complete_async(error_code);
    if (*error_code != MPI_SUCCESS) return;
    MPI_Barrier(fd->comm);
    */

    for (i = 0; i < nprocs; ++i)
      sent_to_proc[i] = 0;
    for (m = 0; m < max_ntimes; m++) {
      rseq_idx = m % num_rseqs;
#ifdef USE_WAITSOME
      num_active = 0;
      for (i = 0; i < num_rseqs; ++i) {
        j = (i + rseq_idx) % num_rseqs;
        flag = 0;
        *error_code = ExchAndWrite_RequestSeq_test(&(rseqs[j]), &flag);
        if (*error_code)
          goto over;
        if (!flag)
          active_rseqs[num_active++] = &(rseqs[j]);
      }
      while (!*error_code && num_active > 0 &&
             active_rseqs[0] == &(rseqs[rseq_idx]))
          *error_code = ExchAndWrite_RequestSeq_waitsome(num_active,
                                                         active_rseqs,
                                                         &num_active);
#else /* !USE_WAITSOME */
      for (i = 0; i < num_rseqs; ++i) {
        j = (i + rseq_idx) % num_rseqs;
        *error_code = ExchAndWrite_RequestSeq_test(&(rseqs[j]), &flag);
        if (*error_code)
          goto over;
      }
      active_rseqs[0] = &(rseqs[rseq_idx]);
      *error_code = ExchAndWrite_RequestSeq_wait(1, active_rseqs);
#endif /* USE_WAITSOME */
      if (*error_code)
        goto over;
      *error_code = start_next_exch_and_write(&(rseqs[rseq_idx]),
                                              sent_to_proc,
                                              send_curr_offlen_ptr,
                                              recv_curr_offlen_ptr);
      if (*error_code)
        goto over;
    }
    for (i = 0; i < num_rseqs; ++i) {
      j = (i + max_ntimes % num_rseqs) % num_rseqs;
      active_rseqs[i] = &(rseqs[j]);
    }
    *error_code = ExchAndWrite_RequestSeq_wait(num_rseqs, active_rseqs);

over:
    for (i = 0; i < num_rseqs; ++i)
      finalize_exch_and_write(&(rseqs[i]));
    ADIOI_Free(sent_to_proc);
    ADIOI_Free(recv_curr_offlen_ptr);
    ADIOI_Free(send_curr_offlen_ptr);
    ADIOI_Free(rseqs);
    ADIOI_Free(active_rseqs);
    ADIOI_Free(off_list);
}

/* Sets error_code to MPI_SUCCESS if successful, or creates an error code
 * in the case of error.
 */
static void ADIOI_LUSTRE_W_Exchange_data(ADIO_File fd, void *buf,
					 char *write_buf,
					 ADIOI_Flatlist_node *flat_buf,
					 ADIO_Offset *offset_list,
					 ADIO_Offset *len_list, int *send_size,
					 int *recv_size, ADIO_Offset off,
					 int size, int *count,
					 int *start_pos,
					 int *sent_to_proc, int nprocs,
					 int myrank, int buftype_is_contig,
					 int contig_access_count,
					 ADIO_Offset *striping_info,
					 ADIOI_Access *others_req,
					 int *send_buf_idx,
					 int *curr_to_proc, int *done_to_proc,
                          int *hole, int iter,
                          MPI_Aint buftype_extent,
					 int *buf_idx,
                          ADIO_Offset **srt_off, int **srt_len, int *srt_num,
                          MPI_Request *requests, int *num_requests,
                          char **send_buf, int *error_code)
{
    int i, j, nprocs_recv, nprocs_send, err;
    MPI_Datatype *recv_types;
    int sum_recv;
    int data_sieving = *hole;
    static char myname[] = "ADIOI_W_EXCHANGE_DATA";

    /* create derived datatypes for recv */
    nprocs_recv = 0;
    for (i = 0; i < nprocs; i++)
	if (recv_size[i])
	    nprocs_recv++;

    recv_types = (MPI_Datatype *) ADIOI_Malloc((nprocs_recv + 1) *
					       sizeof(MPI_Datatype));
    /* +1 to avoid a 0-size malloc */

    j = 0;
    for (i = 0; i < nprocs; i++) {
	if (recv_size[i]) {
	    MPI_Type_hindexed(count[i],
			      &(others_req[i].lens[start_pos[i]]),
			      &(others_req[i].mem_ptrs[start_pos[i]]),
			      MPI_BYTE, recv_types + j);
	    /* absolute displacements; use MPI_BOTTOM in recv */
	    MPI_Type_commit(recv_types + j);
	    j++;
	}
    }

    /* To avoid a read-modify-write,
     * check if there are holes in the data to be written.
     * For this, merge the (sorted) offset lists others_req using a heap-merge.
     */

    *srt_num = 0;
    for (i = 0; i < nprocs; i++)
        *srt_num += count[i];
    if (*srt_off)
        *srt_off = (ADIO_Offset *) ADIOI_Realloc(*srt_off, (*srt_num + 1) * sizeof(ADIO_Offset));
    else
        *srt_off = (ADIO_Offset *) ADIOI_Malloc((*srt_num + 1) * sizeof(ADIO_Offset));
    if (*srt_len)
        *srt_len = (int *) ADIOI_Realloc(*srt_len, (*srt_num + 1) * sizeof(int));
    else
        *srt_len = (int *) ADIOI_Malloc((*srt_num + 1) * sizeof(int));
    /* +1 to avoid a 0-size malloc */

    ADIOI_Heap_merge(others_req, count, *srt_off, *srt_len, start_pos,
		     nprocs, nprocs_recv, *srt_num);

    /* check if there are any holes */
    *hole = 0;
    for (i = 0; i < *srt_num - 1; i++) {
        if ((*srt_off)[i] + (*srt_len)[i] < (*srt_off)[i + 1]) {
            *hole = 1;
	    break;
	}
    }
    /* In some cases (see John Bent ROMIO REQ # 835), an odd interaction
     * between aggregation, nominally contiguous regions, and cb_buffer_size
     * should be handled with a read-modify-write (otherwise we will write out
     * more data than we receive from everyone else (inclusive), so override
     * hole detection
     */
    if (*hole == 0) {
        sum_recv = 0;
        for (i = 0; i < nprocs; i++)
            sum_recv += recv_size[i];
	if (size > sum_recv)
	    *hole = 1;
    }
    /* check the hint for data sieving */
    if (data_sieving == ADIOI_HINT_ENABLE && nprocs_recv && *hole) {
        ADIO_ReadContig(fd, write_buf, size, MPI_BYTE,
                        ADIO_EXPLICIT_OFFSET, off, MPI_STATUS_IGNORE, &err);
        // --BEGIN ERROR HANDLING--
        if (err != MPI_SUCCESS) {
            *error_code = MPIO_Err_create_code(err,
                                               MPIR_ERR_RECOVERABLE,
                                               myname, __LINE__,
                                               MPI_ERR_IO,
                                               "**ioRMWrdwr", 0);
            ADIOI_Free(recv_types);
            return;
        }
        // --END ERROR HANDLING--
    }

    nprocs_send = 0;
    for (i = 0; i < nprocs; i++)
	if (send_size[i])
	    nprocs_send++;

    if (fd->atomicity) {
	/* bug fix from Wei-keng Liao and Kenin Coloma */
      *num_requests = nprocs_send;
    } else {
      *num_requests = nprocs_send + nprocs_recv;

	/* post receives */
	j = 0;
	for (i = 0; i < nprocs; i++) {
	    if (recv_size[i]) {
		MPI_Irecv(MPI_BOTTOM, 1, recv_types[j], i,
                    myrank + i + 100 * iter, fd->comm,
                    requests + nprocs_send + j);
		j++;
	    }
	}
    }

    /* post sends.
     * if buftype_is_contig, data can be directly sent from
     * user buf at location given by buf_idx. else use send_buf.
     */
    if (buftype_is_contig) {
	j = 0;
	for (i = 0; i < nprocs; i++)
	    if (send_size[i]) {
                ADIOI_Assert(buf_idx[i] != -1);
		MPI_Isend(((char *) buf) + buf_idx[i], send_size[i],
                    MPI_BYTE, i, myrank + i + 100 * iter, fd->comm,
                    requests + j);
		j++;
	    }
    } else
        if (nprocs_send) {
	/* buftype is not contig */
	for (i = 0; i < nprocs; i++)
	    if (send_size[i])
		send_buf[i] = (char *) ADIOI_Malloc(send_size[i]);

	ADIOI_LUSTRE_Fill_send_buffer(fd, buf, flat_buf, send_buf, offset_list,
                                   len_list, send_size, requests,
                                   sent_to_proc, nprocs, myrank,
                                   contig_access_count, striping_info,
                                   send_buf_idx, curr_to_proc, done_to_proc,
                                   iter, buftype_extent);
	/* the send is done in ADIOI_Fill_send_buffer */
    }

	/* bug fix from Wei-keng Liao and Kenin Coloma */
    if (fd->atomicity) {
	j = 0;
	for (i = 0; i < nprocs; i++) {
	    if (recv_size[i]) {
		MPI_Recv(MPI_BOTTOM, 1, recv_types[j], i,
                   myrank + i + 100 * iter, fd->comm, MPI_STATUS_IGNORE);
		j++;
	    }
	}
    }

    for (i = 0; i < nprocs_recv; i++)
        MPI_Type_free(recv_types + i);
    ADIOI_Free(recv_types);
}

#define ADIOI_BUF_INCR \
{ \
    while (buf_incr) { \
        size_in_buf = ADIOI_MIN(buf_incr, flat_buf_sz); \
        user_buf_idx += size_in_buf; \
        flat_buf_sz -= size_in_buf; \
        if (!flat_buf_sz) { \
            if (flat_buf_idx < (flat_buf->count - 1)) flat_buf_idx++; \
            else { \
                flat_buf_idx = 0; \
                n_buftypes++; \
            } \
            user_buf_idx = flat_buf->indices[flat_buf_idx] + \
                (ADIO_Offset)n_buftypes*(ADIO_Offset)buftype_extent;  \
            flat_buf_sz = flat_buf->blocklens[flat_buf_idx]; \
        } \
        buf_incr -= size_in_buf; \
    } \
}


#define ADIOI_BUF_COPY \
{ \
    while (size) { \
        size_in_buf = ADIOI_MIN(size, flat_buf_sz); \
        ADIOI_Assert((((ADIO_Offset)(MPIR_Upint)buf) + user_buf_idx) == (ADIO_Offset)(MPIR_Upint)((MPIR_Upint)buf + user_buf_idx)); \
        ADIOI_Assert(size_in_buf == (size_t)size_in_buf);               \
        memcpy(&(send_buf[p][send_buf_idx[p]]), \
               ((char *) buf) + user_buf_idx, size_in_buf); \
        send_buf_idx[p] += size_in_buf; \
        user_buf_idx += size_in_buf; \
        flat_buf_sz -= size_in_buf; \
        if (!flat_buf_sz) { \
            if (flat_buf_idx < (flat_buf->count - 1)) flat_buf_idx++; \
            else { \
                flat_buf_idx = 0; \
                n_buftypes++; \
            } \
            user_buf_idx = flat_buf->indices[flat_buf_idx] + \
                (ADIO_Offset)n_buftypes*(ADIO_Offset)buftype_extent;    \
            flat_buf_sz = flat_buf->blocklens[flat_buf_idx]; \
        } \
        size -= size_in_buf; \
        buf_incr -= size_in_buf; \
    } \
    ADIOI_BUF_INCR \
}

static void ADIOI_LUSTRE_Fill_send_buffer(ADIO_File fd, void *buf,
					  ADIOI_Flatlist_node *flat_buf,
					  char **send_buf,
					  ADIO_Offset *offset_list,
					  ADIO_Offset *len_list, int *send_size,
					  MPI_Request *requests,
					  int *sent_to_proc, int nprocs,
					  int myrank,
					  int contig_access_count,
					  ADIO_Offset *striping_info,
					  int *send_buf_idx,
					  int *curr_to_proc,
					  int *done_to_proc, int iter,
					  MPI_Aint buftype_extent)
{
    /* this function is only called if buftype is not contig */
    int i, p, flat_buf_idx, size;
    int flat_buf_sz, buf_incr, size_in_buf, jj, n_buftypes;
    ADIO_Offset off, len, rem_len, user_buf_idx;

    /* curr_to_proc[p] = amount of data sent to proc. p that has already
     * been accounted for so far
     * done_to_proc[p] = amount of data already sent to proc. p in
     * previous iterations
     * user_buf_idx = current location in user buffer
     * send_buf_idx[p] = current location in send_buf of proc. p
     */

    for (i = 0; i < nprocs; i++) {
	send_buf_idx[i] = curr_to_proc[i] = 0;
	done_to_proc[i] = sent_to_proc[i];
    }
    jj = 0;

    user_buf_idx = flat_buf->indices[0];
    flat_buf_idx = 0;
    n_buftypes = 0;
    flat_buf_sz = flat_buf->blocklens[0];

    /* flat_buf_idx = current index into flattened buftype
     * flat_buf_sz = size of current contiguous component in flattened buf
     */
    for (i = 0; i < contig_access_count; i++) {
	off = offset_list[i];
	rem_len = (ADIO_Offset) len_list[i];

	/*this request may span to more than one process */
	while (rem_len != 0) {
	    len = rem_len;
	    /* NOTE: len value is modified by ADIOI_Calc_aggregator() to be no
	     * longer than the single region that processor "p" is responsible
	     * for.
	     */
	    p = ADIOI_LUSTRE_Calc_aggregator(fd, off, &len, striping_info);

         if (send_buf_idx[p] < send_size[p]) {
           if (curr_to_proc[p] + len > done_to_proc[p]) {
             if (done_to_proc[p] > curr_to_proc[p]) {
               size = (int) ADIOI_MIN(curr_to_proc[p] + len -
                                      done_to_proc[p],
                                      send_size[p] -
                                      send_buf_idx[p]);
               buf_incr = done_to_proc[p] - curr_to_proc[p];
               ADIOI_BUF_INCR
                 ADIOI_Assert((curr_to_proc[p] + len - done_to_proc[p]) == (unsigned)(curr_to_proc[p] + len - done_to_proc[p]));
               buf_incr = (int) (curr_to_proc[p] + len -
                                 done_to_proc[p]);
               ADIOI_Assert((done_to_proc[p] + size) == (unsigned)(done_to_proc[p] + size));
               curr_to_proc[p] = done_to_proc[p] + size;
               ADIOI_BUF_COPY
                 } else {
               size = (int) ADIOI_MIN(len, send_size[p] -
                                      send_buf_idx[p]);
               buf_incr = (int) len;
               ADIOI_Assert((curr_to_proc[p] + size) == (unsigned)((ADIO_Offset)curr_to_proc[p] + size));
               curr_to_proc[p] += size;
               ADIOI_BUF_COPY
                 }
             if (send_buf_idx[p] == send_size[p]) {
               MPI_Isend(send_buf[p], send_size[p], MPI_BYTE, p,
                         myrank + p + 100 * iter, fd->comm,
                         requests + jj);
               jj++;
             }
           } else {
             ADIOI_Assert((curr_to_proc[p] + len) == (unsigned)((ADIO_Offset)curr_to_proc[p] + len));
             curr_to_proc[p] += (int) len;
             buf_incr = (int) len;
             ADIOI_BUF_INCR
               }
         } else {
           buf_incr = (int) len;
           ADIOI_BUF_INCR
             }
	    off += len;
	    rem_len -= len;
	}
    }
    for (i = 0; i < nprocs; i++)
	if (send_size[i])
	    sent_to_proc[i] = curr_to_proc[i];
}

static int ExchAndWrite_RequestSeq_test(struct ExchAndWrite_RequestSeq *rseq,
                                        int *flag)
{
  int errcode;
  errcode = ExchAndWrite_RequestSeq_poll(rseq);
  *flag = (rseq->phase == 0 || rseq->phase == 3 || rseq->error_code);
  return errcode;
}

static int ExchAndWrite_RequestSeq_poll(struct ExchAndWrite_RequestSeq *rseq)
{
  int errcode = MPI_SUCCESS;
  if (!rseq->error_code) {
    switch (rseq->phase) {
     case 1:
      errcode = poll_exch(rseq);
      break;
     case 2:
      errcode = poll_write(rseq);
      break;
     default:
      break;
    }
  }
  return errcode;
}

static int ExchAndWrite_RequestSeq_wait(int count,
                                        struct ExchAndWrite_RequestSeq **rseqs)
{
  int i, j;
  int flag;
  int errcode = MPI_SUCCESS;
  MPI_Request *sub_requests = NULL;
  int num_sub_requests;
  int sub_requests_size = 0;
  for (i = 0; i < count; ++i)
    sub_requests_size += rseqs[i]->num_sub_requests;
  if (sub_requests_size)
    sub_requests = ADIOI_Malloc(sub_requests_size * sizeof(MPI_Request));
  while (!errcode && count > 0) {
    num_sub_requests = 0;
    for (i = 0; i < count; ++i)
      num_sub_requests += rseqs[i]->num_sub_requests;
    if (num_sub_requests > sub_requests_size) {
      sub_requests = ADIOI_Realloc(sub_requests,
                                   num_sub_requests * sizeof(MPI_Request));
      sub_requests_size = num_sub_requests;
    }
    num_sub_requests = 0;
    for (i = 0; i < count; ++i)
      for (j = 0; j < rseqs[i]->num_sub_requests; ++j)
        sub_requests[num_sub_requests++] = rseqs[i]->sub_requests[j];
    errcode = MPI_Waitall(num_sub_requests, sub_requests, MPI_STATUSES_IGNORE);
    for (i = 0; i < count; ++i)
      rseqs[i]->num_sub_requests = 0;
    j = 0;
    for (i = 0; !errcode && i < count; ++i) {
      errcode = ExchAndWrite_RequestSeq_test(rseqs[i], &flag);
      if (!(errcode || flag))
        rseqs[j++] = rseqs[i];
    }
    count = j;
  }
  if (sub_requests)
    ADIOI_Free(sub_requests);
  return errcode;
}

#ifdef USE_WAITSOME
static int ExchAndWrite_RequestSeq_waitsome(int count,
                                            struct ExchAndWrite_RequestSeq **rseqs,
  int *outcount)
{
  int i, j, *indices = NULL;
  int flag, sub_outcount, errcode = MPI_SUCCESS;
  MPI_Request *sub_requests = NULL;
  int num_sub_requests;
  int sub_requests_size = 0;
  *outcount = count;
  for (i = 0; i < count; ++i)
    sub_requests_size += rseqs[i]->num_sub_requests;
  if (sub_requests_size) {
    sub_requests = ADIOI_Malloc(sub_requests_size * sizeof(MPI_Request));
    indices = ADIOI_Malloc(3 * sub_requests_size * sizeof(int));
  }
  else
    return MPI_SUCCESS;
  num_sub_requests = sub_requests_size;
  while (!errcode && count == *outcount) {
    if (num_sub_requests > sub_requests_size) {
      sub_requests = ADIOI_Realloc(sub_requests,
                                   num_sub_requests * sizeof(MPI_Request));
      indices = ADIOI_Realloc(indices, 3 * num_sub_requests * sizeof(int));
      sub_requests_size = num_sub_requests;
    }
    num_sub_requests = 0;
    for (i = 0; i < count; ++i)
      for (j = 0; j < rseqs[i]->num_sub_requests; ++j)
        if (rseqs[i]->sub_requests[j] != MPI_REQUEST_NULL) {
          sub_requests[num_sub_requests] = rseqs[i]->sub_requests[j];
          indices[num_sub_requests + sub_requests_size] = i;
          indices[num_sub_requests++ + 2 * sub_requests_size] = j;
        }
    errcode = MPI_Waitsome(num_sub_requests, sub_requests, &sub_outcount,
                           indices, MPI_STATUSES_IGNORE);
    if (errcode == MPI_UNDEFINED)
      *outcount = 0;
    else {
      for (i = 0; i < sub_outcount; ++i)
        rseqs[indices[indices[i] + sub_requests_size]]->sub_requests[indices[indices[i] + 2 * sub_requests_size]] = MPI_REQUEST_NULL;
      j = 0;
      for (i = 0; !errcode && i < count; ++i) {
        errcode = ExchAndWrite_RequestSeq_test(rseqs[i], &flag);
        if (!(errcode || flag))
          rseqs[j++] = rseqs[i];
      }
      *outcount = j;
    }
    num_sub_requests = 0;
    for (i = 0; i < *outcount; ++i)
      for (j = 0; j < rseqs[i]->num_sub_requests; ++j)
        if (rseqs[i]->sub_requests[j] != MPI_REQUEST_NULL)
          num_sub_requests++;
  }
  ADIOI_Free(sub_requests);
  ADIOI_Free(indices);
  return errcode;
}
#endif /* USE_WAITSOME */

static void init_exch_and_write(ADIO_File fd, void *buf,
                                ADIO_Offset *offset_list,
                                ADIO_Offset *len_list,
                                ADIO_Offset *striping_info,
                                ADIO_Offset *stripe_off_list,
                                ADIO_Offset end_loc,
                                ADIOI_Access *others_req, ADIOI_Access *my_req,
                                ADIOI_Flatlist_node *flat_buf, int **buf_idx,
                                int nprocs, int myrank,
                                int buftype_is_contig, int contig_access_count,
                                MPI_Aint buftype_extent,
                                int iter_offset, int iter_step,
                                struct ExchAndWrite_RequestSeq *rseq)
{
  rseq->fd = fd;
  rseq->buf = buf;
  rseq->offset_list = offset_list;
  rseq->len_list = len_list;
  rseq->striping_info = striping_info;
  rseq->stripe_off_list = stripe_off_list;
  rseq->others_req = others_req;
  rseq->my_req = my_req;
  rseq->flat_buf = flat_buf;
  rseq->buf_idx = buf_idx;

  rseq->nprocs = nprocs;
  rseq->myrank = myrank;
  rseq->buftype_is_contig = buftype_is_contig;
  rseq->contig_access_count = contig_access_count;
  rseq->buftype_extent = buftype_extent;
  rseq->end_loc = end_loc;
  rseq->group_step_size =
    ADIOI_MIN(striping_info[2], striping_info[1]) * striping_info[0];
  rseq->iter_step = iter_step;

  rseq->write_buf = NULL;
  rseq->recv_count = (int *) ADIOI_Malloc(nprocs * sizeof(int));
  rseq->recv_size = (int *) ADIOI_Malloc(nprocs * sizeof(int));
  rseq->send_size = (int *) ADIOI_Malloc(nprocs * sizeof(int));
  rseq->recv_start_pos = (int *) ADIOI_Malloc(nprocs * sizeof (int));
  rseq->this_buf_idx = (int *) ADIOI_Malloc(nprocs * sizeof(int));
  rseq->send_buf_idx = (int *) ADIOI_Malloc(nprocs * sizeof(int));
  rseq->curr_to_proc = (int *) ADIOI_Malloc(nprocs * sizeof(int));
  rseq->done_to_proc = (int *) ADIOI_Malloc(nprocs * sizeof(int));
  rseq->srt_off = NULL;
  rseq->srt_len = NULL;
  rseq->srt_num = 0;
  rseq->num_sub_requests = 0;
  rseq->sub_requests =
    (MPI_Request *) ADIOI_Malloc(2 * nprocs * sizeof(MPI_Request));
  rseq->send_buf = (char **) ADIOI_Calloc(nprocs, sizeof(char *));

  rseq->iter = -iter_step + iter_offset;
  rseq->off = 0;
  rseq->real_size = 0;
  rseq->phase = 0;
  rseq->hole = 0;
}

static void finalize_exch_and_write(struct ExchAndWrite_RequestSeq *rseq)
{
  int i;
  if (rseq->write_buf)
    ADIOI_Free(rseq->write_buf);
  ADIOI_Free(rseq->recv_count);
  ADIOI_Free(rseq->recv_size);
  ADIOI_Free(rseq->send_size);
  ADIOI_Free(rseq->recv_start_pos);
  ADIOI_Free(rseq->this_buf_idx);
  ADIOI_Free(rseq->send_buf_idx);
  ADIOI_Free(rseq->curr_to_proc);
  ADIOI_Free(rseq->done_to_proc);
  if (rseq->srt_len)
    ADIOI_Free(rseq->srt_len);
  if (rseq->srt_off)
    ADIOI_Free(rseq->srt_off);
  ADIOI_Free(rseq->sub_requests);
  for (i = 0; i < rseq->nprocs; ++i)
    if (rseq->send_buf[i])
      ADIOI_Free(rseq->send_buf[i]);
  ADIOI_Free(rseq->send_buf);
}

static int start_next_exch_and_write(struct ExchAndWrite_RequestSeq *rseq,
                                     int *sent_to_proc,
                                     int *send_curr_offlen_ptr,
                                     int *recv_curr_offlen_ptr)
{
  int i, j;
  ADIO_Offset stripe_size = rseq->striping_info[0];
  ADIO_Offset send_off, req_off;

  ADIOI_Assert(rseq->phase == 0 || rseq->phase == 3);

  rseq->iter += rseq->iter_step;
  rseq->phase = 1;
  rseq->off = rseq->stripe_off_list[rseq->iter];
  rseq->real_size =
    (int) ADIOI_MIN((rseq->off / stripe_size + 1) * stripe_size - rseq->off,
                    rseq->end_loc - rseq->off + 1);
  rseq->num_sub_requests = 0;
  rseq->error_code = MPI_SUCCESS;

  /* go through all others_req and my_req to check which will be received
   * and sent in this iteration.
   */

  /* Note that MPI guarantees that displacements in filetypes are in
     monotonically nondecreasing order and that, for writes, the
     filetypes cannot specify overlapping regions in the file. This
     simplifies implementation a bit compared to reads. */

  /*
    off         = start offset in the file for the data to be written in
    this iteration
    real_size   = size of data written (bytes) corresponding to off
    req_off     = offset in the file for a particular contiguous request minus
                  what was satisfied in previous iteration
    send_off    = offset the request needed by other processes in this
                  iteration
    req_len     = size corresponding to req_off
    send_len    = size corresponding to send_off
  */

  /* first calculate what should be communicated */
  for (i = 0; i < rseq->nprocs; i++)
    rseq->recv_count[i] = rseq->recv_size[i] = rseq->send_size[i] = 0;

  ADIOI_Assert(rseq->off == rseq->striping_info[6] ||
               ITER_N(rseq->off, rseq->striping_info, rseq->group_step_size) ==
               rseq->iter);

  if (rseq->write_buf == NULL) {
    for (i = 0; i < rseq->nprocs; i++)
      if (rseq->others_req[i].count) {
        rseq->write_buf = ADIOI_Malloc(rseq->striping_info[0]);
        break;
      }
  }

  for (i = 0; i < rseq->nprocs; i++) {
    if (rseq->my_req[i].count &&
        send_curr_offlen_ptr[i] < rseq->my_req[i].count) {
      rseq->this_buf_idx[i] = rseq->buf_idx[i][send_curr_offlen_ptr[i]];
      for (j = send_curr_offlen_ptr[i]; j < rseq->my_req[i].count; j++) {
        send_off = rseq->my_req[i].offsets[j];
        if (ITER_N(send_off, rseq->striping_info, rseq->group_step_size) ==
            rseq->iter) {
          rseq->send_size[i] += rseq->my_req[i].lens[j];
        } else {
          break;
        }
      }
      send_curr_offlen_ptr[i] = j;
    }
    if (rseq->others_req[i].count
        && recv_curr_offlen_ptr[i] < rseq->others_req[i].count) {
      rseq->recv_start_pos[i] = recv_curr_offlen_ptr[i];
      for (j = recv_curr_offlen_ptr[i]; j < rseq->others_req[i].count; j++) {
        req_off = rseq->others_req[i].offsets[j];
        if (ITER_N(req_off, rseq->striping_info, rseq->group_step_size) ==
            rseq->iter) {
          rseq->recv_count[i]++;
          ADIOI_Assert(req_off - rseq->off < stripe_size);
          ADIOI_Assert((((ADIO_Offset)(MPIR_Upint)rseq->write_buf)+req_off-rseq->off) ==
                       (ADIO_Offset)(MPIR_Upint)(rseq->write_buf+req_off-rseq->off));
          MPI_Address(rseq->write_buf + req_off - rseq->off,
                      &(rseq->others_req[i].mem_ptrs[j]));
          rseq->recv_size[i] += rseq->others_req[i].lens[j];
        } else {
          break;
        }
      }
      recv_curr_offlen_ptr[i] = j;
    }
  }
  rseq->hole = rseq->fd->hints->fs_hints.lustre.ds_in_coll;
  ADIOI_LUSTRE_W_Exchange_data(rseq->fd, rseq->buf, rseq->write_buf,
                               rseq->flat_buf, rseq->offset_list,
                               rseq->len_list, rseq->send_size,
                               rseq->recv_size, rseq->off, rseq->real_size,
                               rseq->recv_count, rseq->recv_start_pos,
                               sent_to_proc, rseq->nprocs, rseq->myrank,
                               rseq->buftype_is_contig,
                               rseq->contig_access_count,
                               rseq->striping_info, rseq->others_req,
                               rseq->send_buf_idx, rseq->curr_to_proc,
                               rseq->done_to_proc, &(rseq->hole), rseq->iter,
                               rseq->buftype_extent, rseq->this_buf_idx,
                               &(rseq->srt_off), &(rseq->srt_len),
                               &(rseq->srt_num),
                               rseq->sub_requests, &(rseq->num_sub_requests),
                               rseq->send_buf, &(rseq->error_code));
  if (rseq->error_code || rseq->num_sub_requests == 0) {
    rseq->phase = 3;
    rseq->num_sub_requests = 0;
  }
  return rseq->error_code;
}

static int poll_exch(struct ExchAndWrite_RequestSeq *rseq)
{
  int exch_complete, flag, i;
  ADIO_Offset block_offset;
  int block_len;

  ADIOI_Assert(rseq->phase == 1);
  ADIOI_Assert(rseq->error_code == MPI_SUCCESS);

  if (rseq->num_sub_requests == 0)
    exch_complete = TRUE;
  else
    rseq->error_code = MPI_Testall(rseq->num_sub_requests, rseq->sub_requests,
                                   &exch_complete, MPI_STATUSES_IGNORE);
  if (rseq->error_code || exch_complete)
    rseq->num_sub_requests = 0;
  if (!rseq->error_code && exch_complete) {
    for (i = 0; i < rseq->nprocs; i++) {
      if (rseq->recv_count[i])
        rseq->phase = 2;
      if (rseq->send_buf[i]) {
        ADIOI_Free(rseq->send_buf[i]);
        rseq->send_buf[i] = NULL;
      }
    }
    if (rseq->phase == 2) {
      /* Although we have recognized the data according to OST index,
       * a read-modify-write will be done if there is a hole between
       * the data. For example: if blocksize=60, xfersize=30 and
       * stripe_size=100, then rank0 will collect data [0, 30] and
       * [60, 90] then write. There is a hole in [30, 60], which will
       * cause a read-modify-write in [0, 90].
       *
       * To reduce its impact on the performance, we can disable data
       * sieving by hint "ds_in_coll".
       */
      /* check the hint for data sieving */
      if (rseq->fd->hints->fs_hints.lustre.ds_in_coll == ADIOI_HINT_ENABLE) {
#ifdef USE_IWRITE_FUNCTIONS
        rseq->num_sub_requests = 1;
        ADIO_IwriteContig(rseq->fd, rseq->write_buf, rseq->real_size,
                          MPI_BYTE, ADIO_EXPLICIT_OFFSET, rseq->off,
                          &(rseq->sub_requests[0]), &(rseq->error_code));
#else /* !USE_IWRITE_FUNCTIONS */
        rseq->phase = 3;
        ADIO_WriteContig(rseq->fd, rseq->write_buf, rseq->real_size,
                         MPI_BYTE, ADIO_EXPLICIT_OFFSET, rseq->off,
                         MPI_STATUS_IGNORE, &(rseq->error_code));
#endif /* USE_IWRITE_FUNCTIONS */
      } else {
        /* if there is no hole, write data in one time;
         * otherwise, write data in several times */
        if (!rseq->hole) {
#ifdef USE_IWRITE_FUNCTIONS
          rseq->num_sub_requests = 1;
          ADIO_IwriteContig(rseq->fd, rseq->write_buf, rseq->real_size,
                            MPI_BYTE, ADIO_EXPLICIT_OFFSET, rseq->off,
                            &(rseq->sub_requests[0]), &(rseq->error_code));
#else /* !USE_IWRITE_FUNCTIONS */
          rseq->phase = 3;
          ADIO_WriteContig(rseq->fd, rseq->write_buf, rseq->real_size,
                           MPI_BYTE, ADIO_EXPLICIT_OFFSET, rseq->off,
                           MPI_STATUS_IGNORE, &(rseq->error_code));
#endif /* USE_IWRITE_FUNCTIONS */
        } else {
          // Don't do async io in this case. TODO: maybe
          // allow a sequence of async ios, if the number of
          // blocks is small enough.
          rseq->phase = 3;
          block_offset = -1;
          block_len = 0;
          for (i = 0; i < rseq->srt_num; ++i) {
            if (rseq->srt_off[i] < rseq->off + rseq->real_size &&
                rseq->srt_off[i] >= rseq->off) {
              if (block_offset == -1) {
                block_offset = rseq->srt_off[i];
                block_len = rseq->srt_len[i];
              } else {
                if (rseq->srt_off[i] == block_offset + block_len) {
                  block_len += rseq->srt_len[i];
                } else {
                  ADIO_WriteContig(rseq->fd,
                                   rseq->write_buf + block_offset - rseq->off,
                                   block_len,
                                   MPI_BYTE, ADIO_EXPLICIT_OFFSET,
                                   block_offset, MPI_STATUS_IGNORE,
                                   &(rseq->error_code));
                  if (rseq->error_code != MPI_SUCCESS) {
                    block_offset = -1;
                    break;
                  }
                  block_offset = rseq->srt_off[i];
                  block_len = rseq->srt_len[i];
                }
              }
            }
          }
          if (block_offset != -1) {
            ADIO_WriteContig(rseq->fd,
                             rseq->write_buf + block_offset - rseq->off,
                             block_len,
                             MPI_BYTE, ADIO_EXPLICIT_OFFSET,
                             block_offset, MPI_STATUS_IGNORE,
                             &(rseq->error_code));
          }
        }
      }
    }
    else
      rseq->phase = 3;
  }
  return rseq->error_code;
}

static int poll_write(struct ExchAndWrite_RequestSeq *rseq)
{
  int write_complete;

  ADIOI_Assert(rseq->phase == 2);
  ADIOI_Assert(rseq->error_code == MPI_SUCCESS);

  if (rseq->num_sub_requests == 0)
    write_complete = TRUE;
  else
    rseq->error_code = MPI_Testall(rseq->num_sub_requests, rseq->sub_requests,
                                   &write_complete, MPI_STATUSES_IGNORE);
  if (rseq->error_code || write_complete)
    rseq->num_sub_requests = 0;
  if (!rseq->error_code && write_complete)
    rseq->phase = 3;
  return rseq->error_code;
}
