/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement non-blocking 1D communication scheme
//       along X axis.
//
// SUMMARY:
//     - 1D splitting along X
// NEW:
//     - >>> Non-blocking communications <<<
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex3(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation then ex1
	lbm_comm_init_ex1(comm, total_width, total_height);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex3(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 1D communication with non-blocking MPI functions.
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[9] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)
	
	//example to access cell
	//double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	//double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);

	// Computing neighbors rank
	int left_rank = (comm->rank_x > 0) ? comm->rank_x - 1 : MPI_PROC_NULL;
	int right_rank = (comm->rank_x < comm->nb_x - 1) ? comm->rank_x + 1 : MPI_PROC_NULL;

	// Standardize tag to factorize
	int tag_to_left = 0;
	int tag_to_right = 1;

	// Size of a column, accounting the fact that each cell contains DIRECTIONS times doubles
	int col_size = comm->height * DIRECTIONS;

	// Initializing request and status array
	int request_count = 0;
	MPI_Request requests[4];
	MPI_Status statuses[4];

	/* ---------------------------------------------------------------
	left to right : send from left neighbor | reception in the right neighbor
		Note : we send from column width-2 as column width - 1 is a ghost one
			   we receive from column width = 0 as its the left ghost cells column
	--------------------------------------------------------------- */

	// We send directly the whole column ! It's way more elegant and efficient that using a loop like I first wanted to do
	MPI_Isend(lbm_mesh_get_cell(mesh, comm->width - 2, 0),
			 col_size, MPI_DOUBLE,
			 right_rank, tag_to_right,
			 comm->communicator, &requests[request_count++]);

	MPI_Irecv(lbm_mesh_get_cell(mesh, 0, 0),
			 col_size, MPI_DOUBLE,
			 left_rank, tag_to_right,
			 comm->communicator, &requests[request_count++]);

	/* ---------------------------------------------------------------
	right to left : send from right neighbor | reception in the left neighbor
	Note : we send from column 1 as column width = 0 is a ghost one
		   we receive from column width - 1 as its the right ghost cells column
	--------------------------------------------------------------- */
	MPI_Isend(lbm_mesh_get_cell(mesh, 1, 0),
			 col_size, MPI_DOUBLE,
			 left_rank, tag_to_left,
			 comm->communicator,&requests[request_count++]);

	MPI_Irecv(lbm_mesh_get_cell(mesh, comm->width - 1, 0),
			 col_size, MPI_DOUBLE,
			 right_rank, tag_to_left,
			 comm->communicator, &requests[request_count++]);

	MPI_Waitall(request_count, requests, statuses);
}
