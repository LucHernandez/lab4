/*****************************************************
	AUTHOR  : Sébastien Valat
	MAIL    : sebastien.valat@univ-grenoble-alpes.fr
	LICENSE : BSD
	YEAR    : 2021
	COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
//
// GOAL: Implement a 1D communication scheme along
//       X axis with blocking communications.
//
// SUMMARY:
//     - 1D splitting along X
//     - Blocking communications
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex1(lbm_comm_t *comm, int total_width, int total_height)
{
	int rank;
	int comm_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	//
	// TODO: calculate the splitting parameters for the current task.
	//
	// HINT: You can look in exercise_0.c to get an example for the sequential case.
	//
	if (total_width % comm_size != 0)
	{
		fatal("Total width not dividable by comm_size !");
	}

	// TODO: calculate the number of tasks along X axis and Y axis.
	comm->nb_x = comm_size;
	comm->nb_y = 1;

	// TODO: calculate the current task position in the splitting
	comm->rank_x = rank;
	comm->rank_y = 0;

	// TODO : calculate the local sub-domain size (do not forget the
	//        ghost cells). Use total_width & total_height as starting
	//        point.
	comm->width = (total_width / comm_size) + 2;
	comm->height = total_height + 2;

	// TODO : calculate the absolute position in the global mesh.
	//        without accounting the ghost cells
	//        (used to setup the obstable & initial conditions).
	comm->x = comm->rank_x * (total_width / comm_size);
	comm->y = 0;

	comm->communicator = MPI_COMM_WORLD;

	MPI_Comm newcomm;
	MPI_Comm_split(comm->communicator, comm->rank_y, comm->rank_x, &newcomm);

	int newrank;
	MPI_Comm_rank(newcomm, &newrank);

	// if debug print comm
	lbm_comm_print(comm);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex1(lbm_comm_t *comm, lbm_mesh_t *mesh)
{
	//
	// TODO: Implement the 1D communication with blocking MPI functions (MPI_Send & MPI_Recv)
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[DIRECTIONS] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)

	// example to access cell
	// double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	// double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);

	// Computing neighbors rank
	int left_rank = (comm->rank_x > 0) ? comm->rank_x - 1 : MPI_PROC_NULL;
	int right_rank = (comm->rank_x < comm->nb_x - 1) ? comm->rank_x + 1 : MPI_PROC_NULL;

	// Standardize tag to factorize
	int tag_to_left = 0;
	int tag_to_right = 1;

	// Size of a column, accounting the fact that each cell contains DIRECTIONS times doubles
	int col_size = comm->height * DIRECTIONS;

	// Necessary to allow recv to return informations
	MPI_Status status;

	/* ---------------------------------------------------------------
	left to right : send from left neighbor | reception in the right neighbor
		Note : we send from column width-2 as column width - 1 is a ghost one
			   we receive from column width = 0 as its the left ghost cells column
	--------------------------------------------------------------- */

	// We send directly the whole column ! It's way more elegant and efficient that using a loop like I first wanted to do
	MPI_Send(lbm_mesh_get_cell(mesh, comm->width - 2, 0),
			 col_size, MPI_DOUBLE,
			 right_rank, tag_to_right,
			 comm->communicator);

	MPI_Recv(lbm_mesh_get_cell(mesh, 0, 0),
			 col_size, MPI_DOUBLE,
			 left_rank, tag_to_right,
			 comm->communicator, &status);

	/* ---------------------------------------------------------------
	right to left : send from right neighbor | reception in the left neighbor
	Note : we send from column 1 as column width = 0 is a ghost one
		   we receive from column width - 1 as its the right ghost cells column
	--------------------------------------------------------------- */
	MPI_Send(lbm_mesh_get_cell(mesh, 1, 0),
			 col_size, MPI_DOUBLE,
			 left_rank, tag_to_left,
			 comm->communicator);

	MPI_Recv(lbm_mesh_get_cell(mesh, comm->width - 1, 0),
			 col_size, MPI_DOUBLE,
			 right_rank, tag_to_left,
			 comm->communicator, &status);
}