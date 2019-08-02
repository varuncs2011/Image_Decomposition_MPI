#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "pgmio.h"
#include "arralloc.h"


#define THRESHOLD 0.1
#define MAXITER   1500
#define PRINTFREQ  200
#define DEFAULT 1
#define NDIMS 2
#define DIM_X 0
#define DIM_Y 1

#define LEFT_NEIGHBOR 0
#define RIGHT_NEIGHBOR 1
#define UP_NEIGHBOR 2
#define DOWN_NEIGHBOR 3
#define OFFSET 2



MPI_Comm comm, cart_comm;
int size, rank;
int neighbor[4];
int coords[NDIMS];
int dims[NDIMS];
int length_buf, width_buf;
int x_size, y_size;
double start, end;
char *input_image, *output_image;

double **edge, **old, **new;



void write(int, int );
double boundaryval(int , int );
void readedges(int , int );
void processing(int , int);
void setboundary(int, int);
void finalize(double **, double *, double **, double *, double **, double *);




int main(int argc, char **argv)
{
	int length, width;       // original image size
	int curlength, curwidth; // current process image chunk size
	double *edge_content, *old_content, *new_content;


	//input_image = argv[1];
	//output_image = argv[2];

	input_image = argv[1];
	output_image= argv[2];
	//printf("Image file name %s",input_image);

	// calculate the source image size using pgmsize method
	pgmsize(input_image, &length, &width);

	length_buf = length;
	width_buf = width;
	if (rank == 0)
	{
	//	printf("Image size: %d x %d.\n", length, width);
	//	printf("Max iteration :%d\n", MAXITER);
	}

	// Call MPI INIT
	MPI_Init(NULL, NULL);
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	// Periodicity is only in vertical direction. So horizontal is 0 and vertical is 1
	int periods[NDIMS] = {0, 1};

	// create dimension (i*j) based on the number of process and dimensions
	MPI_Dims_create(size, NDIMS, dims);

	//Create cartision topology based on the dimenstion and periodicity
	MPI_Cart_create(comm, NDIMS, dims, periods, 0, &cart_comm);
	MPI_Cart_coords(cart_comm, rank, NDIMS, coords);

	if (rank == 0)
	{
//		printf("Mapping processes into %d x %d grids.\n", dims[DIM_X], dims[DIM_Y]);
	}
//	printf("Rank %d in coords (%d , %d)\n", rank, coords[DIM_X], coords[DIM_Y]);

	// find the neighbors

	// In 1st direction left & right
	MPI_Cart_shift(cart_comm, 0, 1, &neighbor[LEFT_NEIGHBOR], &neighbor[RIGHT_NEIGHBOR]);

	// in 2nd direction down and up
	MPI_Cart_shift(cart_comm, 1, 1, &neighbor[DOWN_NEIGHBOR], &neighbor[UP_NEIGHBOR]);


	//divide the image based on the dimensions calculated earlier
	x_size = length / dims[DIM_X];
	y_size = width / dims[DIM_Y];
	curlength = x_size;
	curwidth = y_size;

//	printf("Rank %d curlength %d, curwidth %d\n", rank, curlength, curwidth);


	allocatememory(&new, &new_content, curlength + 2, curwidth + 2);
	allocatememory(&old, &old_content, curlength + 2, curwidth + 2);
	allocatememory(&edge, &edge_content, curlength + 2, curwidth + 2);

	// readEdge file
	readedges(curlength, curwidth);

	//set boundary conditions
	setboundary(curlength, curwidth);

	// Construct image
	processing(curlength, curwidth);

	// collect data and write to disc
	write(curlength, curwidth);

	// Finalize
         finalize(edge, edge_content, old, old_content, new, new_content);

	return 0;
}


// read the edge file.
void readedges(int curlength, int curwidth)
{

	pgmreadchunk(input_image, &(edge[1][1]), curlength, curwidth, coords[DIM_X] * x_size, coords[DIM_Y] * y_size, OFFSET);

	// Initialize old to 255 (all white)
	for (int i = 0; i <= curlength + 1; ++i)
	{
		for (int j = 0; j <= curwidth + 1; ++j)
		{
			old[i][j] = 255.0;
		}
	}
}

// set boundary conditions.
void setboundary(int curlength, int curwidth)
{
	// Boundary conditions for bottom and top
	double value;
	if (MPI_PROC_NULL == neighbor[UP_NEIGHBOR])
	{
		for (int i = 1; i <= curlength; ++i)
		{
			value = boundaryval(coords[DIM_X] * x_size + i, length_buf);
			old[i][curwidth + 1] = (int)(255.0 * (1.0 - value));
		}
	}
	if (MPI_PROC_NULL == neighbor[DOWN_NEIGHBOR])
	{
		for (int i = 1; i <= curlength; ++i)
		{
			value = boundaryval(coords[DIM_X] * x_size + i, length_buf);
			old[i][0] = (int)(255.0 * value);
		}
	}
}


// This method does the actual processing for the image construction by swapping halos with the neighbour process.
void processing(int length, int width)
{

	// ROW datatype declaration. ROW is vector which is used to send halos.
	MPI_Datatype  row;
	MPI_Type_vector(length + 2, 1, width + 2, MPI_DOUBLE, &row);
	MPI_Type_commit(&row);
	MPI_Request tracking_id[4];
	MPI_Status status[4];

	// All the processes stops here. This is used to calculate the start time.
	MPI_Barrier(comm);
	start = MPI_Wtime();

	// Send and receive tags should match.
	int  sendLeft= 2;
	int  recvRight= 2;
	int  sendRight = 3;
	int  recvLeft= 3;

//	printf("Rank  %d iteration \n", rank);

	for (int iter = 1; iter <= MAXITER; ++iter)
	{
		// call Non blocking send(Isend) for all the neighbours of this process.


		// send all the elements equal to the width +2 at once to left . Since the halos are located in adjacent memeory location.we dont neet vector for left and right send.
		MPI_Issend(&(old[1][0]), width + 2,MPI_DOUBLE , neighbor[LEFT_NEIGHBOR], sendLeft, comm, &tracking_id[LEFT_NEIGHBOR]);
		MPI_Issend(&(old[length][0]), width+2, MPI_DOUBLE, neighbor[RIGHT_NEIGHBOR], sendRight, comm, &tracking_id[RIGHT_NEIGHBOR]);

		//  send one row at a time. the halo cells are not present in continuous memory location. We need to use Vector to send data up and down .
		MPI_Issend(&(old[0][1]), 1, row, neighbor[DOWN_NEIGHBOR], 1, comm, &tracking_id[DOWN_NEIGHBOR]);
		MPI_Issend(&(old[0][width]), 1, row, neighbor[UP_NEIGHBOR], 1, comm, &tracking_id[UP_NEIGHBOR]);

		// Calculation can be done in parallel as the send routine was Non-Blocking
		for (int i = 2; i <= length - 1; ++i)
		{
			for (int j = 2; j <= width - 1; ++j)
			{
				new[i][j] = 0.25 * (old[i - 1][j] + old[i + 1][j] + old[i][j - 1] + old[i][j + 1] - edge[i][j]);
			}
		}

		// Blokcing Receive
		MPI_Recv(&(old[0][0]), width + 2, MPI_DOUBLE, neighbor[LEFT_NEIGHBOR], recvLeft, comm, &status[LEFT_NEIGHBOR]);
		MPI_Recv(&(old[length + 1][0]), width + 2, MPI_DOUBLE, neighbor[RIGHT_NEIGHBOR], recvRight, comm, &status[RIGHT_NEIGHBOR]);

		// Similar to send, we can use row vector to receive the data.
		MPI_Recv(&(old[0][0]), 1, row, neighbor[DOWN_NEIGHBOR], 1, comm, &status[DOWN_NEIGHBOR]);
		MPI_Recv(&(old[0][width + 1]), 1, row, neighbor[UP_NEIGHBOR], 1, comm, &status[UP_NEIGHBOR]);


		// Halo calculation
		for (int j = 1; j <= width; ++j)
		{
			new[1][j] = 0.25 * (old[0][j] + old[2][j] + old[1][j - 1] + old[1][j + 1] - edge[1][j]);
			new[length][j] = 0.25 * (old[length - 1][j] + old[length + 1][j] + old[length][j - 1] + old[length][j + 1] - edge[length][j]);
		}

		for (int i = 1; i <= length; ++i)
		{
			new[i][1] = 0.25 * (old[i - 1][1] + old[i + 1][1] + old[i][0] + old[i][2] - edge[i][1]);
			new[i][width] = 0.25 * (old[i - 1][width] + old[i + 1][width] + old[i][width - 1] + old[i][width + 1] - edge[i][width]);
		}

		// Call MPI wait ALL and pass tracking ids and status array
		MPI_Waitall(2,tracking_id, status);

		// This code is to check if we can stop the iteration if we are not getting a change in any of the pixel . maxium absolute change of any pixel
		double delta = 0, total_delta = 0;
		double sum = 0, total_sum = 0;
		for (int i = 1; i <= length; ++i)
		{
			for (int j = 1; j <= width; ++j)
			{
				float d = fabs(new[i][j] - old[i][j]);
				if(d > delta)
				{
					delta = d;
				}
				sum += new[i][j];

				// copy new to old
				old[i][j] = new[i][j];
			}
		}

		// Call MPI All reduce
		MPI_Allreduce(&delta, &total_delta, 1, MPI_DOUBLE, MPI_MAX, comm);

		if (total_delta <= THRESHOLD)
		{
			if (rank == 0)
			{
				printf("Finished at %d  for delta : %.5f\n", iter, total_delta);
			}
			break;
		}

		if (0 == iter % PRINTFREQ)
		{
			MPI_Allreduce(&sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
			if (rank == 0)
			{
				printf(" Iteration %d : Average pixels %.3f\n", iter, total_sum / (length_buf * width_buf));
			}
		}
	}
	end = MPI_Wtime();
	MPI_Barrier(comm);

	if (rank == 0)
	{
		printf("Iteration time %.5f\n", end - start);
	}
}


// This method accumulate the partial images from all the process into master process at rank 0

void write(int curlength, int curwidth)
{
	//
	if (rank != 0)
	{
		// all the process except rank one sends their partial image to rank 0

		MPI_Datatype send_to_master;

		// define a vector datatype to send subsection of the 2d array.
		MPI_Type_vector(curlength, curwidth, curwidth + 2, MPI_DOUBLE, &send_to_master);
		MPI_Type_commit(&send_to_master);

		MPI_Ssend(&(old[1][1]), 1, send_to_master, 0, 1, comm);

	}
	if (rank == 0)
	{
		double **masterbuf;
		double *master_chunk;
		allocatememory(&masterbuf, &master_chunk, length_buf, width_buf);

		// copy the  partial image of master process at rank0 into masterbuf
		for (int i = 0; i <= curlength - 1; ++i)
		{
			for (int j = 0; j <= curwidth - 1; ++j)
			{
				masterbuf[i][j] = old[i + 1][j + 1];
			}
		}

		// call MPI receive for all the process
		for (int r = 1; r <= size - 1; ++r)
		{
			int recv_coords[NDIMS]; // coords of rank 0.
			MPI_Cart_coords(cart_comm, r, NDIMS, recv_coords);

			int recv_length = x_size;
			int recv_width = y_size;

			MPI_Datatype recv_from_slaves;
			MPI_Type_vector(recv_length, recv_width, width_buf, MPI_DOUBLE, &recv_from_slaves);
			MPI_Type_commit(&recv_from_slaves);

			// copy the partial images in masterbuf copy into 0,1,Width for process 0
			MPI_Recv(&(masterbuf[recv_coords[DIM_X] * x_size][recv_coords[DIM_Y] * y_size]), 1, recv_from_slaves, r, 1, comm, MPI_STATUS_IGNORE);
		}

		//output_image=inpugsize
		// write the output to output file
		pgmwrite(output_image, master_chunk, length_buf, width_buf);

		// release the memory for masterbuf
		releasememory(&masterbuf, &master_chunk);
	}
}


double boundaryval(int i, int m)
{
	double val;

	val = 2.0 * ((double)(i - 1)) / ((double)(m - 1));
	if (i >= m / 2 + 1)
		val = 2.0 - val;

	return val;
}



// Release the memory allocated to 2D array and the the array of pointers which is storing the base address of 1D array which is subsection of 2d Array.
void finalize(double **edge, double *edge_content, double **old, double *old_content, double **new, double *new_content)
{
	// release the memory
	releasememory(&new, &new_content);
	releasememory(&old, &old_content);
	releasememory(&edge,&edge_content);

	// MPI Finalize
	MPI_Finalize();
}
