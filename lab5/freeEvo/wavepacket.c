/********************************************************************

 Time evolution of a quantum wave packet, in 1D, 2D or 3D

 MPI parallel version compiled using flag -DMPI

 Version 1.0 (April 2013)

 Copyright 2013
 C.M. Dion, A. Hashemloo, and G. Rahali

********************************************************************/

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>
#include "wavepacket.h"
#include "wavepacket_private.h"

#ifdef MPI
#include "mpi.h"
#include <fftw3-mpi.h>
#endif


/*** Main program ***/

int
main (int argc, char *argv[])
{
  size_t alloc_local;
  unsigned int nt_inner, nt_outer;
  double t;
  parameters params;
  output_flags output;
  char results_file[80], wf_text[80], wf_bin[80];

  FILE *fp_results;

#ifdef MPI
  // Initialize MPI
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &params.size);
  MPI_Comm_rank (MPI_COMM_WORLD, &params.rank);
  fftw_mpi_init ();
#else
  // Default values for serial version
  params.size = 1;
  params.rank = 0;
#endif
  
  if (params.rank == 0)  // Input and output is managed by processor 0
    {
      // Get name of parameter file
      if (argc < 2)
	{
	  fprintf (stderr, "usage: %s parameter_file [user-defined arguments]\n", argv[0]);
	  abort ();
	}
      
      // Get parameters for the simulation
      read_parameters (&params, &output, &nt_inner, &nt_outer, argv[1], 
		       results_file, wf_text, wf_bin);
    }

#ifdef MPI
  // Transfer parameters to all processors
  distribute_parameters (&params, &output, &nt_inner, &nt_outer);
#endif

  // Determine if the simulation is 1D
  if (params.ny > 1)
    {
      params.is1D = false;
    }
  else
    {
      params.is1D = true;
    }

  // Instantiate carriers of information for kinetic and energy operators
  kin_data kd;
  ener_data ed;

#ifdef MPI
  ptrdiff_t local_n0, local_0_start, local_n1, local_1_start;  
  
  // Determine the part of each processor
  if (params.is1D)
    {
      // 1D system 
      ptrdiff_t alloc_in = fftw_mpi_local_size_1d
	(params.nx, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE,
	 &local_n0, &local_0_start, &local_n1, &local_1_start);
      
      ptrdiff_t alloc_out = fftw_mpi_local_size_1d
	(params.nx, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE,
	 &local_n1, &local_1_start, &local_n0, &local_0_start);

      alloc_local = MAX(alloc_in, alloc_out);

      ed.n1 = kd.n1 = params.nx;
      ed.n1_local = kd.n1_local = local_n1;
      ed.n1_0 = kd.n1_0 = local_1_start;
      ed.n2 = kd.n2 = 1;
      ed.n3 = kd.n3 = 1;
    }
  else
    {
      // 2D-3D system
      alloc_local = fftw_mpi_local_size_3d_transposed 
	(params.nx, params.ny, params.nz, MPI_COMM_WORLD, 
	 &local_n0, &local_0_start, &local_n1, &local_1_start);      

      ed.n1 = kd.n1 = params.ny;
      ed.n1_local = kd.n1_local = local_n1;
      ed.n1_0 = kd.n1_0 = local_1_start;
      ed.n2 = kd.n2 = params.nx;
      ed.n3 = kd.n3 = params.nz;
    }

  params.nx_local = local_n0;
  params.nx0 = local_0_start;
#else
  // Parameters for single processor are the same as global parameters
  alloc_local = params.n;
  params.nx_local = params.nx;
  params.nx0 = 0;
  
  ed.n1 = kd.n1 = params.nx;
  ed.n1_local = kd.n1_local = params.nx_local;
  ed.n1_0 = kd.n1_0 = 0;
  ed.n2 = kd.n2 = params.ny;
  ed.n3 = kd.n3 = params.nz;
#endif  

  // Total local number of grid points 
  params.n_local = params.nx_local * params.ny * params.nz;

  // Make spatial grid
  make_grid (&params);

  // Create arrays for the wave function 
  fftw_complex *psi = fftw_alloc_complex (alloc_local);
  fftw_complex *psitmp = fftw_alloc_complex (alloc_local);
  
  if (psi == NULL || psitmp == NULL)
    {
      fprintf (stderr, "Out of memory while allocating wave functions.\n");
      abort ();
    }

  // Array to store the potential
  double * const pot = (double *) malloc (params.n_local * sizeof (double));
  
  if (pot == NULL)
    {
      fprintf (stderr, "Out of memory while allocating potential array.\n");
      abort ();
    }

  // Recurring factor in the exponential of the evolution operator
  double complex idt_hbar = -I * params.dt / params.hbar;

  // Initialize kinetic energy operator
  initialize_kinetic (params, psi, &kd);

  // Initialize the potential function (user-supplied routine)
  initialize_potential (params, argc, argv);

  // Initialize calculation of the energy
  initialize_energy (params, psi, psitmp, &ed);

  // Initialize wave function (user-supplied routine)
  initialize_wf (params, argc, argv, psi);

  // Initialize user-defined observables
  if (output.user)
    initialize_user_observe (params, argc, argv);

  // Initial wave function must be kept if calculating autocorrelation function
  double complex *psii;
  if (output.autoc)
    {
      psii = (double complex *) malloc 
	(params.n_local * sizeof (double complex));
      
      if (psii == NULL)
	{
	  fprintf (stderr, "Out of memory while allocating initial wave functions.\n");
	  abort ();
	}
      
      for (size_t l = 0; l < params.n_local; ++l)
	psii[l] = psi[l];
    }

  // Header of results file for values that will be printed out
  if (params.rank == 0)
    print_header (params, output, results_file, &fp_results);
  
  // Calculate observables for initial wave function
  if (output.yes)
    observe (params, output, 0., ed, pot, psii, psi, psitmp, fp_results);

  // Flush buffer before time evolution (in case of problems / for testing)
  if (params.rank == 0)
    {
      fflush (fp_results);
    }

  /* Time evolution */
  t = 0.;

  for (size_t outer = 0; outer < nt_outer; ++outer)  // Main loop
    {
      // half kinetic operator
      kinetic_operator (kd, psi, HALF_STEP);
      
      t -= 0.5 * params.dt;  // set time for inner loop
      
      for (size_t inner = 0; inner < nt_inner; ++inner)  // inner loop
	{
	  t += params.dt;  // accumulate time

	  // potential operator
	  potential_operator (params, idt_hbar, t, pot, psi);
	  
	  // double kinetic operator
	  kinetic_operator (kd, psi, FULL_STEP);
	}
      
      t += params.dt; // accumulate time 
      
      // potential operator
      potential_operator (params, idt_hbar, t, pot, psi);
      
      // half Kinetic operator
      kinetic_operator (kd, psi, HALF_STEP);

      t += 0.5 * params.dt; // Current time
      
      if (output.yes) 
	{
	  // Print out intermediate results
	  observe (params, output, t, ed, pot, psii, psi, psitmp, fp_results);
	}
    }
    
  // Clean up

  if (params.rank == 0)
    {
      fclose (fp_results);
    }
  
  fftw_destroy_plan (kd.forward);
  fftw_destroy_plan (kd.backward);
  fftw_destroy_plan (ed.forward);
  fftw_destroy_plan (ed.backward);
  fftw_free (psitmp);
  if (output.autoc)
    {
      free (psii);
    }
  free (kd.ek1);
  free (kd.ek2);
  free (kd.ek3);
  free (kd.ek1_2);
  free (kd.ek2_2);
  free (kd.ek3_2);
  free (ed.k1);
  free (ed.k2);
  free (ed.k3);
  free (pot);
  free (params.x2);
  free (params.y2);
  free (params.z2);

  // Print out final wave function
  if (output.wf == 1 || output.wf == 3)
    print_wf (params, psi, wf_text);

  if (output.wf == 2 || output.wf == 3)
    print_wf_bin (params, psi, wf_bin);

  // Finish clean up
  fftw_free (psi);

  free (params.x);
  free (params.y);
  free (params.z);

#ifdef MPI
  MPI_Finalize ();
  fftw_mpi_cleanup ();
#else
  fftw_cleanup ();
#endif

  return 0;
}



/*** Routines for input and output ***/

void 
read_parameters (parameters * const params, output_flags * const output, 
		 unsigned int * const nt_inner, unsigned int * const nt_outer, 
		 const char * const parameter_file, 
		 char * const results_file, char * const wf_text, 
		 char * const wf_bin)
{
  /** Read parameters from the file specified on the command line **/

  if (params->rank != 0)
    {
      /* (This function should only be called by the rank=0 processor in a 
	 multi-processor environment) */
      fprintf (stderr, "Function read_parameters called by a processor other than rank 0 (rank = %d).\n", params->rank);
      abort ();
    }

  bool error = false, have_results_file = false, have_wf_output_text = false, 
    have_wf_output_bin = false, have_mass = false, 
    have_nx = false, have_ny = false, have_nz = false, have_x_min = false, 
    have_x_max = false, have_y_min = false, have_y_max = false, 
    have_z_min = false, have_z_max = false, have_dt = false, have_nt = false, 
    have_nprint = false;
  int units = 0, nt, nprint;
  
  // Open file for reading
  FILE *fp_in = fopen (parameter_file, "r");
  if (fp_in == NULL)
    {
      fprintf (stderr, "Could not open parameter file %s\n", 
	       parameter_file);
      abort ();
    }        
  
  // Set all output flags to false before reading parameters
  output->yes = false;
  output->norm = false;
  output->energy = false;
  output->x_avg = false;
  output->y_avg = false;
  output->z_avg = false;
  output->sx = false;
  output->sy = false;
  output->sz = false;
  output->autoc = false;
  output->user = false;

  /* Read parameters */
  char key[30], value[80];
  while (fscanf (fp_in, "%s = %s", key, value) != EOF) 
    {
      if (!strcmp (key, "results_file"))
	{
	  strcpy (results_file, value);
	  have_results_file = true;
	}
      else if (!strcmp (key, "wf_output_text"))
	{
	  strcpy (wf_text, value);
	  have_wf_output_text = true;
	}
      else if (!strcmp (key, "wf_output_binary"))
	{
	  strcpy (wf_bin, value);
	  have_wf_output_bin = true;
	}
      else if (!strcmp (key, "units"))
	{
	  if (!strcmp (value, "AU"))
	    {
	      units = 0;
	    }
	  else if (!strcmp (value, "SI"))
	    {
	      units = 1;
	    }
	  else 
	    {
	      fprintf (stderr, "In parameter file %s, value %s for key units unknown, ignored\n", parameter_file, value);
	    }
	}
      else if (!strcmp (key, "mass"))
	{ 
	  params->mass = atof (value);
	  have_mass = true;
	}
      else if (!strcmp (key, "nx"))
	{
	  params->nx = atoi (value);
	  have_nx = true;
	}
      else if (!strcmp (key, "ny"))
	{
	  params->ny = atoi (value);
	  have_ny = true;
	}
      else if (!strcmp (key, "nz"))
	{
	  params->nz = atoi (value);
	  have_nz = true;
	}
      else if (!strcmp (key, "x_min"))
	{
	  params->x_min = atof (value);
	  have_x_min = true;
	}
      else if (!strcmp (key, "x_max"))
	{
	  params->x_max = atof (value);
	  have_x_max = true;
	}
      else if (!strcmp (key, "y_min"))
	{
	  params->y_min = atof (value);
	  have_y_min = true;
	}
      else if (!strcmp (key, "y_max"))
	{
	  params->y_max = atof (value);
	  have_y_max = true;
	}
      else if (!strcmp (key, "z_min"))
	{
	  params->z_min = atof (value);
	  have_z_min = true;
	}
      else if (!strcmp (key, "z_max"))
	{
	  params->z_max = atof (value);
	  have_z_max = true;
	}
      else if (!strcmp (key, "dt"))
	{
	  params->dt = atof (value);
	  have_dt = true;
	}
      else if (!strcmp (key, "nt"))
	{
	  nt = atoi (value);
	  have_nt = true;
	}
      else if (!strcmp (key, "nprint"))
	{
	  nprint = atoi (value);
	  have_nprint = true;
	}
      else if (!strcmp (key, "output"))      
	{
	  if (!strcmp (value, "norm"))
	    {
	      output->norm = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "energy"))
	    {
	      output->energy = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "x_avg"))
	    {
	      output->x_avg = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "y_avg"))
	    {
	      output->y_avg = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "z_avg"))
	    {
	      output->z_avg = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "sx"))
	    {
	      output->sx = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "sy"))
	    {
	      output->sy = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "sz"))
	    {
	      output->sz = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "autocorrelation"))
	    {
	      output->autoc = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "user_defined"))
	    {
	      output->user = true;
	      output->yes = true;
	    }
	  else
	    {
	      fprintf (stderr, "In parameter file %s, output value '%s' unknown, ignored\n", parameter_file, value);
	    }
	}
      else
	{
	  fprintf (stderr, "In parameter file %s, parameter %s unknown, ignored\n", parameter_file, key);
	}
    } 

  fclose (fp_in);
  
  /* Check for missing parameters */
  if (!have_mass)
    {
      fprintf (stderr, "Parameter mass missing from file %s\n",
	       parameter_file);
      error = true;
    }

  // 1D means x
  if (!have_nx)
    {
      fprintf (stderr, "Parameter nx missing from file %s\n",
	       parameter_file);
      fprintf (stderr, "(1D system must be defined along x, 2D system in xy plane)\n");
      error = true;
    }
  else if (!have_ny && !have_nz)
    {
      fprintf (stdout, "Parameter ny and nz missing from file %s, values set to 1\n",
	       parameter_file);
      params->ny = 1;
      params->nz = 1;
    }
  else if (!have_ny)
    {
      fprintf (stdout, "Parameter ny missing from file %s, value set to 1\n",
	       parameter_file);
      params->ny = 1;
    }
  else if (!have_nz)
    {
      fprintf (stdout, "Parameter nz missing from file %s, value set to 1\n",
	       parameter_file);
      params->nz = 1;
    }
  
  // 2D means xy
  if (have_nx && params->ny <= 1 && params->nz > 1)
    {
      fprintf (stderr, "2D system must be defined along the x and y axes\n");
      error = true;
    }

  if (!have_x_min)
    {
      fprintf (stderr, "Parameter x_min missing from file %s\n", 
	       parameter_file);
      error = true;
    } 
  
  if (!have_x_max)
    {
      fprintf (stderr, "Parameter x_max missing from file %s\n", 
	       parameter_file);
      error = true;
    }
  
  if (!have_y_min)
    {
      if(params->ny > 1)
	{
	  fprintf (stderr, "Parameter y_min missing from file %s\n", parameter_file);
	  error = true;
	} 
      else
	{
	  fprintf (stdout, "Parameter y_min missing from file %s, value set to 0.\n", parameter_file);
	  params->y_min = 0.;
	}
    }
  
  if (!have_y_max)
    {
      if (params->ny > 1)
	{
	  fprintf (stderr, "Parameter y_max missing from file %s\n",
		   parameter_file);
	  error = true;
	}
      else
	{
	  params->y_max = params->y_min;
	}
    }


  if (!have_z_min)
    {
      if(params->nz > 1)
	{
	  fprintf (stderr, "Parameter z_min missing from file %s\n",
		   parameter_file);
	  error = true;
	} 
      else
	{
	  fprintf (stdout, "Parameter z_min missing from file %s, value set to 0.\n", parameter_file);
	  params->z_min = 0.;
	}
    }

  if (!have_z_max)
    {
      if (params->nz > 1)
	{
	  fprintf (stderr, "Parameter z_max missing from file %s\n",
		   parameter_file);
	  error = true;
	}
      else
	{
	  params->z_max = params->z_min;
	}
    }
  
  if (!have_dt)
    {
      fprintf (stderr, "Parameter dt missing from file %s\n",
	       parameter_file);
      error = true;
    }

  if (!have_nt)
    {
      fprintf (stderr, "Parameter nt missing from file %s\n",
	       parameter_file);
      error = true;
    }

  if (!have_nprint && output->yes)
    {
      fprintf (stderr, "Parameter nprint missing from file %s\n",
	       parameter_file);
      error = true;
    }
  
  if (error) abort ();

  if (!have_results_file)
    {
      fprintf (stdout, "No results_file specified, using default value 'results'\n");
      strcpy (results_file, "results");
    }

  /* Check validity of the parameters */
  if (params->mass <= 0.)
    {
      fprintf (stderr, "mass = %g, must be greater than zero.\n", params->mass);
      error = true;
    }

  if (params->nx == 0)
    {
      fprintf (stderr, "nx has to be greater than zero\n");
      fprintf (stderr, "(1D system must be defined along x, 2D system in xy plane))\n");
      error = true;
    }
#ifdef MPI
  else if (params->nx < params->size)
    {
      fprintf (stderr, "fewer grid points in x than number of processors\n");
      error = true;
    }
#endif

  if (params->ny == 0)
    {
      if (params->nz > 1)
	{
	  fprintf (stderr, "ny has to be greater than zero\n");
	  fprintf (stderr, "(2D system must be defined in xy plane))\n");
	  error = true;
	}
      else
	{
	  fprintf (stdout, "ny should be greater than zero\n");
	  fprintf (stdout, "Assuming 1D system and setting ny = 1\n");
	  params->ny = 1;
	}
    }
  
  if (params->nz == 0)
    {
      fprintf (stdout, "nz should be greater than zero\n");
      fprintf (stdout, "Assuming reduced dimensionality and setting nz = 1\n");
      params->nz = 1;
    }
  
  if (params->x_max <= params->x_min)
    {
      fprintf (stderr, "x_max = %lg must be greater than x_min = %lg\n", 
	       params->x_max, params->x_min);
      error = true;
    }

  if (params->ny > 1 && params->y_max <= params->y_min)
    {
      fprintf (stderr, "y_max = %lg must be greater than y_min = %lg\n", 
	       params->y_max, params->y_min);
      error = true;
    }

  if (params->nz > 1 && params->z_max <= params->z_min)
    {
      fprintf (stderr, "z_max = %lg must be greater than z_min = %lg\n", 
	       params->z_max, params->z_min);
      error = true;
    }
  
  if (error) abort ();

  if (output->yes && nprint > nt)
    {
      fprintf (stdout, "* Warning *\n");
      fprintf (stdout, "nprint = %d is greater than nt = %d.\n",
	       nprint, nt);
      output->yes = 0;
    }
  else if (output->yes && nt % nprint != 0)
    {
      fprintf (stdout, "* Warning *\n");
      fprintf (stdout, "nprint = %d is not an integer multiple of nt = %d.\n",
	       nprint, nt);
      
      nt = (nt / nprint + 1) * nprint;

      fprintf (stdout, "nt is modified to %d.\n", nt);
    }

  // Figure out the final output of the wave function
  output->wf = 0;
  
  if (have_wf_output_text)
    output->wf = 1;
  
  if (have_wf_output_bin)
    output->wf += 2;
  
  if (!output->yes)
    {
      if (output->wf == 0)
	{
	  fprintf (stderr, "Program will produce no output.\nAborting...\n");
	  abort ();
	}
      else
	{
	  fprintf (stdout, "No output during time evolution,\nonly final wave function will be printed out\n");
	  
	  *nt_inner = nt - 1;
	  *nt_outer = 1;
	}
    }
  else
    {
      *nt_inner = nprint - 1;
      *nt_outer = nt/nprint;
    }
  
  if (output->y_avg && params->ny == 1)
    {
      fprintf (stdout, "No grid along y, so y_avg will not be calculated\n");
      output->y_avg = 0;
    }

  if (output->z_avg && params->nz == 1)
    {
      fprintf (stdout, "No grid along z, so z_avg will not be calculated\n");
      output->z_avg = 0;
    }

  if (output->sy && params->ny == 1)
    {
      fprintf (stdout, "No grid along y, so sy will not be calculated\n");
      output->sy = 0;
    }

  if (output->sz && params->nz == 1)
    {
      fprintf (stdout, "No grid along z, so sz will not be calculated\n");
      output->sz = 0;
    }

  /* Derived parameters */
  params->n = params->nx * params->ny * params->nz;

  if (units == 1)
    {
      params->hbar = 1.054571726e-34; // Planck constant over 2 pi (2010 CODATA)
    }
  else
    {
      params->hbar = 1.;
    }
  
  /* Print out all parameters */
  fprintf (stdout, "*** Parameters ***\n");
  if (units == 0)
    {
      fprintf (stdout, "units = AU\n");
    }
  else
    {
      fprintf (stdout, "units = SI\n");
    }
  fprintf (stdout, "mass = %.16g\n", params->mass);
  fprintf (stdout, "nx = %zu\n", params->nx);
  fprintf (stdout, "ny = %zu\n", params->ny);
  fprintf (stdout, "nz = %zu\n", params->nz);
  fprintf (stdout, "x_min = %.16g\n", params->x_min);
  fprintf (stdout, "x_max = %.16g\n", params->x_max);
  fprintf (stdout, "y_min = %.16g\n", params->y_min);
  fprintf (stdout, "y_max = %.16g\n", params->y_max);
  fprintf (stdout, "z_min = %.16g\n", params->z_min);
  fprintf (stdout, "z_max = %.16g\n", params->z_max);
  fprintf (stdout, "dt = %.16g\n", params->dt);
  fprintf (stdout, "nt = %d\n", nt);
  fprintf (stdout, "nprint = %d\n\n", nprint);
  
  if (output->yes)
    fprintf (stdout, "results_file = %s\n", results_file);

  if (output->wf == 1 || output->wf == 3)
    fprintf (stdout, "wf_output_text = %s\n", wf_text);

  if (output->wf == 2 || output->wf == 4)
    fprintf (stdout, "wf_output_bin = %s\n", wf_bin);

  if (output->yes)
    {
      fprintf (stdout, "\n");
      if (output->norm)
	fprintf (stdout, "output = norm\n");
      if (output->energy)
	fprintf (stdout, "output = energy\n");
      if (output->x_avg)
	fprintf (stdout, "output = x_avg\n");
      if (output->y_avg)
	fprintf (stdout, "output = y_avg\n");
      if (output->z_avg)
	fprintf (stdout, "output = z_avg\n");
      if (output->sx)
	fprintf (stdout, "output = sx\n");
      if (output->sy)
	fprintf (stdout, "output = sy\n");
      if (output->sz)
	fprintf (stdout, "output = sz\n");
      if (output->autoc)
	fprintf (stdout, "output = autocorrelation\n");
      if (output->user)
	fprintf (stdout, "output = user_defined\n");
    }

  return;
}



void
print_header (const parameters params, const output_flags output, 
	      const char * const results_file, FILE **fp_results)
{
  /** Header of the results file for values that will be printed out **/
  
  if (params.rank != 0)
    {
      /* (This function should only be called by the rank=0 processor in a 
	 multi-processor environment) */
      fprintf (stderr, "Function print_header called by a processor other than rank 0 (rank = %d).\n", params.rank);
      abort ();
    }

  // Open file for results
  *fp_results = fopen (results_file, "w");
  if (*fp_results == NULL)
    {
      fprintf (stderr, "Could not open file %s for output.\n", results_file);
      abort ();
    }   

  fprintf (*fp_results, "# time      ");

  if (output.norm)
    fprintf (*fp_results, " norm        ");

  if (output.energy)
    fprintf (*fp_results, "  energy      ");

  if (output.x_avg)
    fprintf (*fp_results, "  x_avg       ");

  if (output.y_avg)
    fprintf (*fp_results, "  y_avg       ");

  if (output.z_avg)
    fprintf (*fp_results, "  z_avg       ");

  if (output.sx)
    fprintf (*fp_results, " sx          ");

  if (output.sy)
    fprintf (*fp_results, " sy          ");

  if (output.sz)
    fprintf (*fp_results, " sz          ");

  if (output.autoc)
    fprintf (*fp_results, " autoc       ");

  fprintf (*fp_results, "\n");

  return;
}



void 
print_wf (const parameters params, const fftw_complex * const psi,
	  const char * const wf_text)
{
  /** Print out the wave function psi in a text file **/
  
  size_t index_x, index_xy, index;
  FILE *fp;

  if (params.rank == 0)
    {
      fp = fopen (wf_text, "w");
      if (fp == NULL)
	{
	  fprintf (stderr, "Could not open file %s for output.\n", wf_text);
	  abort ();
	}   
    }      
  
  bool print_y = (params.ny > 1);
  bool print_z = (params.nz > 1);
  
  if (params.rank == 0)
    {
      // Print header
      fprintf (fp, "# x");
      
      if (print_y)
	fprintf (fp, " y");
      
      if (print_z)
	fprintf (fp, " z");
      
      fprintf (fp, " psi_r psi_i psi2 \n");

      // Processor zero prints its wave function
      
      for (size_t i = 0; i < params.nx_local; ++i)
	{
	  index_x = i * params.nz * params.ny;
	  for (size_t j = 0; j < params.ny; ++j)
	    {
	      index_xy = j * params.nz + index_x;
	      for (size_t k = 0; k < params.nz; ++k)
		{
		  index = k + index_xy;
		  
		  fprintf (fp, "%13.6e ", params.x[i]);
		  
		  if (print_y)
		    fprintf (fp, "%13.6e ", params.y[j]);
		  
		  if (print_z)
		    fprintf (fp, "%13.6e ", params.z[k]);
		  
		  double psi2 = cabs (psi[index]);
		  psi2 *= psi2;
		  fprintf (fp, "%13.6e %13.6e %.6e\n",
			   creal (psi[index]), cimag (psi[index]), psi2);
		}
	    }
	}
    }

#ifdef MPI
  // Processor 0 gathers the wave function from the other processors to print
  int nx, npsi;
  MPI_Status status;

  if (params.rank == 0)
    {
      for (int core = 1; core < params.size; ++core)
	{ 
	  // Get size of the data
	  MPI_Recv (&nx, 1, MPI_INT, core, 1, MPI_COMM_WORLD, &status);

	  // Allocate storage for x and wave function
	  double *x = (double *) malloc (nx * sizeof (double));
	  if (x == NULL)
	    {
	      fprintf (stderr, "Out of memory while allocating x in print_wf.\n");
	      abort ();
	    }

	  npsi = 2 * nx * params.ny * params.nz;
	  double *psii = (double *) malloc (npsi * sizeof (double));
	  if (psii == NULL)
	    {
	      fprintf (stderr, "Out of memory while allocating psii in print_wf.\n");
	      abort ();
	    }

	  // Get x
	  MPI_Recv (x, nx, MPI_DOUBLE, core, 2, MPI_COMM_WORLD, &status);

	  // Get wave function
	  MPI_Recv (psii, npsi, MPI_DOUBLE, core, 3, MPI_COMM_WORLD, &status);

	  // Print data
	  for (size_t i = 0; i < nx; ++i)
	    {
	      index_x = i * params.nz * params.ny;
	      for (size_t j = 0; j < params.ny; ++j)
		{
		  index_xy = j * params.nz + index_x;
		  for (size_t k = 0; k < params.nz; ++k)
		    {
		      index = k + index_xy;
		      
		      fprintf (fp, "%13.6e ", x[i]);
		      
		      if (print_y)
			fprintf (fp, "%13.6e ", params.y[j]);
		      
		      if (print_z)
			fprintf (fp, "%13.6e ", params.z[k]);
		      
		      double psii2 = (SQUARE(psii[2*index]) 
				      + SQUARE(psii[2*index+1]));
		      fprintf (fp, "%13.6e %13.6e %13.6e\n",
			       psii[2*index], psii[2*index+1], psii2);
		    }
		}
	    }
	  // Free memory
	  free (x);
	  free (psii);
	}
    }
  else
    {
      // Send size of local data
      nx = params.nx_local;
      MPI_Send (&nx, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD); 

      // Send x
      MPI_Send (params.x, nx, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);

      // Send wave function
      npsi = 2 * nx * params.ny * params.nz;
      MPI_Send ((double *)psi, npsi, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
    }
#endif

  // Close file
  if (params.rank == 0)
    fclose (fp);
  
  return;
}



void 
print_wf_bin (const parameters params, const double complex * const psi,
	  const char * const wf_bin)
{
  /** Print out the wave function psi in a binary file **/
  
  FILE *fp;

  if (params.rank == 0)
    {
      fp = fopen (wf_bin, "wb");
      if (fp == NULL)
	{
	  fprintf (stderr, "Could not open file %s for output.\n", wf_bin);
	  abort ();
	}   
      
      // Write size of the wave function
      fwrite (&params.nx, sizeof (size_t), 1, fp);
      fwrite (&params.ny, sizeof (size_t), 1, fp);
      fwrite (&params.nz, sizeof (size_t), 1, fp);

      // Processor zero prints its wave function
      fwrite (psi, sizeof (double complex), params.n_local, fp); 
    }

#ifdef MPI
  // Processor 0 gathers the wave function from the other processors to print
  unsigned int nx;
  MPI_Status status;

  if (params.rank == 0)
    {
      for (int core = 1; core < params.size; ++core)
	{ 
	  // Get size of the data
	  MPI_Recv (&nx, 1, MPI_UNSIGNED, core, 1, MPI_COMM_WORLD, &status);
      
	  // Allocate storage for x and wave function
	  double *psii = (double*) malloc (2 * nx * sizeof (double));

	  // Get wave function
	  MPI_Recv (psii, 2*nx, MPI_DOUBLE, core, 3, MPI_COMM_WORLD, &status);

	  // Print data
	  fwrite (psii, sizeof (double), 2*nx, fp); 

	  // Free memory
	  free (psii);
	}
    }
  else
    {
      // Send size of local data
      nx = params.nx_local;
      MPI_Send (&nx, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD); 

      // Send wave function
      MPI_Send ((double *)psi, 2*nx, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
    }
#endif

  if (params.rank == 0)
    fclose (fp);
  
  return;
}



void 
read_wf_bin (const parameters params, const char * const wf_bin,
	     double complex * const psi)
{
  /** Read the wave function psi from a binary file **/

  double complex *psitmp;
  FILE *fp;
  
  if (params.rank == 0)
    {
      fp = fopen (wf_bin, "rb");
      if (fp == NULL)
	{
	  fprintf (stderr, "Could not open file %s for input.\n", wf_bin);
	  abort ();
	}   
      
      // Read size of the stored wave function
      int check;
      size_t n[3];
      check = fread (n, sizeof (size_t), 3, fp);

      if (check < 3)
	{
	  fprintf (stderr, "Error reading grid size in file %s.\n", wf_bin);
	  abort ();
	}
      else if (n[0] != params.nx || n[1] != params.ny || n[2] != params.nz)
	{
	  fprintf (stderr, "Grid size in file %s, nx = %ld, ny = %ld, nz = %ld,\n",
		   wf_bin, n[0], n[1], n[2]);
	  fprintf (stderr, "not compatible with current parameters, nx = %ld, ny = %ld, nz = %ld.\n", 
		   params.nx, params.ny, params.nz);
	  abort ();
	}

      // Read the wave function      
#ifdef MPI
      psitmp = (double complex *) malloc (params.n * sizeof (double complex));
#else
      psitmp = psi;
#endif
      check = fread (psitmp, sizeof (double complex), params.n, fp); 

      if (check < params.n)
	{
	  fprintf (stderr, "Error reading wave function in file %s.\n", wf_bin);
	  abort ();
	}

      fclose (fp);
    }
  
#ifdef MPI
  // Share the wave function among processors
  distribute_wf (params, psitmp, psi);

  if (params.rank == 0)
    free (psitmp);
#endif
  
  return;
}



/*** Routines for initialization ***/

void
make_grid (parameters *params)
{
  /** Construct the three-dimensional spatial grid **/

  params->x = (double *) malloc (params->nx_local * sizeof (double));
  params->y = (double *) malloc (params->ny * sizeof (double));
  params->z = (double *) malloc (params->nz * sizeof (double));
  params->x2 = (double *) malloc (params->nx_local * sizeof (double));
  params->y2 = (double *) malloc (params->ny * sizeof (double));
  params->z2 = (double *) malloc (params->nz * sizeof (double));

  if (params->x == NULL || params->y == NULL || params->z == NULL
      || params->x2 == NULL || params->y2 == NULL || params->z2 == NULL)
    {
      fprintf (stderr, "Out of memory while creating spatial grid.\n");
      abort ();
    }

  params->dx = (params->x_max - params->x_min) / (double)(params->nx - 1);

  if (params->ny == 1) 
    params->dy = 1.;
  else
    params->dy = (params->y_max - params->y_min) / (double)(params->ny - 1);

  if (params->nz == 1)   
    params->dz = 1.;
  else
    params->dz = (params->z_max - params->z_min) / (double)(params->nz - 1);

  for (size_t i = 0; i < params->nx_local; ++i)
    {
      params->x[i] = params->x_min + (double)(i + params->nx0) * params->dx;
      params->x2[i] = SQUARE(params->x[i]);
    }

  for (size_t j = 0; j < params->ny; ++j)
    {
      params->y[j] = params->y_min + (double)j * params->dy;
      params->y2[j] = SQUARE(params->y[j]);
    }

  for (size_t k = 0; k < params->nz; ++k)
    {
      params->z[k] = params->z_min + (double)k * params->dz;
      params->z2[k] = SQUARE(params->z[k]);
    }

  return;
}



void
initialize_kinetic (const parameters params, fftw_complex *psi, kin_data *kd)
{
  /** Initialize kinetic energy operator **/
  size_t last, mid;

  kd->ek1 = (double complex*) malloc (kd->n1_local * sizeof (double complex)); 
  kd->ek2 = (double complex*) malloc (kd->n2 * sizeof (double complex));
  kd->ek3 = (double complex*) malloc (kd->n3 * sizeof (double complex));
  kd->ek1_2 = (double complex*) malloc (kd->n1_local * sizeof (double complex));
  kd->ek2_2 = (double complex*) malloc (kd->n2 * sizeof (double complex));
  kd->ek3_2 = (double complex*) malloc (kd->n3 * sizeof (double complex));
  
  if (kd->ek1 == NULL || kd->ek2 == NULL || kd->ek3 == NULL
      || kd->ek1_2 == NULL || kd->ek2_2 == NULL || kd->ek3_2 == NULL)

    {
      fprintf (stderr, "Out of memory while creating spatial grid.\n");
      abort ();
    }
  
  // Kinetic operator along the first dimension
  double complex factor1 = -I * PISQR * params.hbar * params.dt / params.mass;
  double complex factor2;

#ifdef MPI
  if (!params.is1D)
    {
      // x and y are transposed, so y is first dimension
      factor2 = factor1 / SQUARE((double)params.ny * params.dy);
    }
  else
#endif
    {
      // x is the first dimension
      factor2 = factor1 / SQUARE((double)params.nx * params.dx);
    }
  
  last = kd->n1_0 + kd->n1_local;
  mid = MIN(kd->n1/2 + 1, last);
  
  for (size_t i = kd->n1_0; i < mid; ++i)
    {
      kd->ek1[i - kd->n1_0] = cexp (factor2 * (double)SQUARE(i));
      kd->ek1_2[i - kd->n1_0] = cexp (2. * factor2 * (double)SQUARE(i));
    }
  
  mid = MAX(kd->n1/2 + 1, kd->n1_0);
  
  for (size_t i = mid; i < last; ++i)
    {
      kd->ek1[i - kd->n1_0] = cexp (factor2 * (double)SQUARE(kd->n1 - i));
      kd->ek1_2[i - kd->n1_0] = cexp (2. * factor2 
				      * (double)SQUARE(kd->n1 - i));
    }


// Kinetic operator along the second dimension  
#ifdef MPI
  if (!params.is1D)
    {
      // x and y are transposed, so x is 2nd dimension
      factor2 = factor1 / SQUARE((double)params.nx * params.dx);
    }
  else      
#endif
    {
      // y is 2nd dimension
      factor2 = factor1 / SQUARE((double)params.ny * params.dy);
    }

  kd->ek2[0] = 1.;
  kd->ek2_2[0] = 1.;
  for (size_t j = 1; j <= kd->n2/2; ++j)
    {
      double complex exponential = cexp (factor2 * (double)SQUARE(j));
      kd->ek2[j] = exponential;
      kd->ek2[kd->n2 - j] = exponential;
      
      exponential = cexp (2. * factor2 * (double)SQUARE(j));
      kd->ek2_2[j] = exponential;
      kd->ek2_2[kd->n2 - j] = exponential;
    }
  
  // Kinetic operator along the third dimension (z)
  factor2 = factor1 / SQUARE((double)params.nz * params.dz);
  
  // FFT scaling factor is absorbed in ek3 and ek3_2
  double fft_factor = 1. / (double)params.n; 
  
  kd->ek3[0] = fft_factor;
  kd->ek3_2[0] = fft_factor;
  for (size_t k = 1; k <= kd->n3/2; ++k)
    {
      double complex exponential = cexp (factor2 * (double)SQUARE(k)) 
	* fft_factor;
      kd->ek3[k] = exponential;
      kd->ek3[kd->n3 - k] = exponential;

      exponential = cexp (2. * factor2 * (double)SQUARE(k)) * fft_factor;
      kd->ek3_2[k] = exponential;
      kd->ek3_2[kd->n3 - k] = exponential;
    }

  // Prepare FFT
#ifdef MPI
  if (params.is1D)
    {
      // 1D, no transpose
      kd->forward = fftw_mpi_plan_dft_1d (params.nx, psi, psi, MPI_COMM_WORLD,
					  FFTW_FORWARD, FFTW_ESTIMATE);
      kd->backward = fftw_mpi_plan_dft_1d (params.nx, psi, psi, MPI_COMM_WORLD,
					   FFTW_BACKWARD, FFTW_ESTIMATE);
    }
  else
    {
      // x and y are transposed
      kd->forward = fftw_mpi_plan_dft_3d (params.nx, params.ny, params.nz, psi, 
					  psi, MPI_COMM_WORLD, FFTW_FORWARD, 
					  FFTW_ESTIMATE | 
					  FFTW_MPI_TRANSPOSED_OUT);
      kd->backward = fftw_mpi_plan_dft_3d (params.nx, params.ny, params.nz, psi,
					   psi, MPI_COMM_WORLD, FFTW_BACKWARD, 
					   FFTW_ESTIMATE | 
					   FFTW_MPI_TRANSPOSED_IN);
    }
#else
  kd->forward = fftw_plan_dft_3d (params.nx, params.ny, params.nz, 
				  psi, psi, FFTW_FORWARD, FFTW_ESTIMATE);
  kd->backward = fftw_plan_dft_3d (params.nx, params.ny, params.nz, 
				   psi, psi, FFTW_BACKWARD, FFTW_ESTIMATE);
#endif

  return;
}



void
initialize_energy (const parameters params, fftw_complex *psi, 
		   fftw_complex *psitmp, ener_data *ed)
{
  /** Initialize calculation of the energy **/
  size_t last, mid;
  double factor1, factor2, tmp;

  ed->k1 = (double *) malloc (ed->n1_local * sizeof (double));
  ed->k2 = (double *) malloc (ed->n2 * sizeof (double));
  ed->k3 = (double *) malloc (ed->n3 * sizeof (double));

  if (ed->k1 == NULL || ed->k2 == NULL || ed->k3 == NULL)
    {
      fprintf (stderr, "Out of memory in initialize_energy.\n");
      abort ();
    }  
  
  factor1 = 2. * PISQR * SQUARE(params.hbar) / params.mass;

  factor1 /= (double)params.n; // FFT scaling factor

#ifdef MPI
  if (!params.is1D)
    {
      // x and y are transposed, so y is the first dimension
      factor2 = factor1 / SQUARE((double)params.ny * params.dy);
    }
  else
#endif      
    {
      // x is the first dimension
      factor2 = factor1 / SQUARE((double)params.nx * params.dx);
    }

  last = ed->n1_0 + ed->n1_local;
  mid = MIN(ed->n1/2 + 1, last);
  
  for (size_t i = ed->n1_0; i < mid; ++i)
    {
      ed->k1[i - ed->n1_0] = factor2 * (double)SQUARE(i);
    }
  
  mid = MAX(ed->n1/2 + 1, ed->n1_0);
  
  for (size_t i = mid; i < last; ++i)
    {
      ed->k1[i - ed->n1_0] = factor2 * (double)SQUARE(ed->n1 - i);
    }      

#ifdef MPI
  if (!params.is1D)
    {
      // x and y are transposed, so x is 2nd dimension
      factor2 = factor1 / SQUARE((double)params.nx * params.dx);
    }      
  else
#endif
    {
      // y is the 2nd dimension      
      factor2 = factor1 / SQUARE((double)params.ny * params.dy);
    }
      
  ed->k2[0] = 0.;
  for (size_t j = 1; j <= ed->n2/2; ++j)
    {
      tmp = factor2 * (double)SQUARE(j);
      ed->k2[j] = tmp;
      ed->k2[ed->n2 - j] = tmp;
    }
  
  // z is the 3rd dimension
  factor2 = factor1 / SQUARE((double)params.nz * params.dz);
  
  ed->k3[0] = 0.;
  for (size_t k = 1; k <= ed->n3/2; ++k)
    {
      tmp = factor2 * (double)SQUARE(k);
      ed->k3[k] = tmp;
      ed->k3[ed->n3 - k] = tmp;
    }
  
  // Prepare FFT
#ifdef MPI
  if (params.is1D)
    {
      // 1D, no transpose
      ed->forward = fftw_mpi_plan_dft_1d (params.nx, psi, psitmp, 
					  MPI_COMM_WORLD,
					  FFTW_FORWARD, FFTW_ESTIMATE);
      ed->backward = fftw_mpi_plan_dft_1d (params.nx, psitmp, psitmp,
					   MPI_COMM_WORLD,
					   FFTW_BACKWARD, FFTW_ESTIMATE);
    }
  else
    {
      // 2D or 3D, x and y are transposed
      ed->forward = fftw_mpi_plan_dft_3d (params.nx, params.ny, params.nz, psi, 
					  psitmp, MPI_COMM_WORLD, FFTW_FORWARD, 
					  FFTW_ESTIMATE | 
					  FFTW_MPI_TRANSPOSED_OUT);
      ed->backward = fftw_mpi_plan_dft_3d (params.nx, params.ny, params.nz, 
					   psitmp, psitmp,
					   MPI_COMM_WORLD, FFTW_BACKWARD, 
					   FFTW_ESTIMATE | 
					   FFTW_MPI_TRANSPOSED_IN);
    }
#else
  ed->forward = fftw_plan_dft_3d (params.nx, params.ny, params.nz, psi, 
				  psitmp, FFTW_FORWARD, FFTW_ESTIMATE);
  ed->backward = fftw_plan_dft_3d (params.nx, params.ny, params.nz, psitmp, 
				   psitmp, FFTW_BACKWARD, FFTW_ESTIMATE);
#endif

  return;
}



/*** Routines implementing operators ***/

void
kinetic_operator (const kin_data kd, fftw_complex *psi, const int step)
{
  /** Apply the kinetic operator to wave function psi **/

  /* If step == HALF_STEP, corresponds to the kinetic operator with 
     a factor 1/2 */

  // Forward FFT to momentum space
  fftw_execute (kd.forward);   

  // Apply kinetic operator in momentum space
  double complex *ek1, *ek2, *ek3;
  if (step == HALF_STEP)
    {
      ek1 = kd.ek1;
      ek2 = kd.ek2;
      ek3 = kd.ek3;
    }
  else
    {
      ek1 = kd.ek1_2;
      ek2 = kd.ek2_2;
      ek3 = kd.ek3_2;
    }

  size_t index_1, index_2, index;
  double complex ek1ek2;
  for (size_t i = 0; i < kd.n1_local; ++i)
    {
      index_1 = i * kd.n3 * kd.n2;
      for (size_t j = 0; j < kd.n2; ++j)
	{
	  index_2 = j * kd.n3 + index_1;
	  ek1ek2 = ek1[i] * ek2[j];
	  for (size_t k = 0; k < kd.n3; ++k)
	    {
	      index = k + index_2;
	      psi[index] *= ek3[k] * ek1ek2;
	    }
	}
    }
    
  // Backward FFT to position space
  fftw_execute (kd.backward);   

  return;
}



void
potential_operator (const parameters params, const double complex idt_hbar, 
		    const double t, double * const pot, 
		    fftw_complex * const psi) 
{
  /** Apply the potential operator to wave function psi **/

  // Get potential
  potential (params, t, pot);

  for (size_t l = 0; l < params.n_local; ++l)
    {
      psi[l] *= cexp (idt_hbar * pot[l]);
    }

  return;
}



/*** Routines for calculating observables ***/

void
observe (const parameters params, const output_flags output, const double t, 
	 const ener_data ed, double * const pot, 
	 const double complex * const psii, fftw_complex *psi, 
	 fftw_complex *psitmp, FILE *fp_results)
{
  /** Calculate observables and print results to file fp_results **/

  /* Norm */
  double norm_ = norm (params, psi);
  double inv_norm2 = SQUARE(norm_);

  if (params.rank == 0)
    {
      fprintf (fp_results, "%.6e", t);
      
      if (output.norm)
	fprintf (fp_results, " %.6e", norm_);
      
    }

  /* Energy */
  if (output.energy)
    {
      double ener = energy (params, t, ed, pot, psi, psitmp);
      
      if (params.rank == 0)
	fprintf (fp_results, " %13.6e", ener);
    }
  
  /* Expectation values of position  */
  double x_avg, y_avg, z_avg;
  
  if (output.x_avg)
    {
      // <x>
      x_avg = expectation1D (params, 1, params.x, psi) * inv_norm2;
    
      if (params.rank == 0)
	fprintf (fp_results, " %13.6e", x_avg);
    }
  
  if (output.y_avg)
    {  
      // <y>
      if (params.ny != 1)
	y_avg = expectation1D (params, 2, params.y, psi);
      else
	y_avg = params.y_min;
      
      y_avg *= inv_norm2;
    
      if (params.rank == 0)
	fprintf (fp_results, " %13.6e", y_avg);
    }
  
  if (output.z_avg)
    {
      // <z>
      if (params.nz != 1)
	z_avg = expectation1D (params, 3, params.z, psi);
      else
	z_avg = params.z_min;
      
      z_avg *= inv_norm2;

      if (params.rank == 0)
	fprintf (fp_results, " %13.6e", z_avg);
    }

  /* Standard deviations (width of the wave packet) */
  if (output.sx)
    {
      // sx = sqrt (<x^2> - <x>^2)

      if (!output.x_avg)
	{
	  // <x>
	  x_avg = expectation1D (params, 1, params.x, psi) * inv_norm2;
	}

      double sx = sqrt (expectation1D (params, 1, params.x2, psi) * inv_norm2
		 - SQUARE(x_avg));

      if (params.rank == 0)
	fprintf (fp_results, " %.6e", sx);
    }
  
  if (output.sy)
    {
      // sy = sqrt (<y^2> - <y>^2)

      if (!output.y_avg)
	{
	  // <y>
	  y_avg = expectation1D (params, 2, params.y, psi) * inv_norm2;
	}

      double sy = sqrt (expectation1D (params, 2, params.y2, psi) * inv_norm2
			- SQUARE(y_avg));

      if (params.rank == 0)
	fprintf (fp_results, " %.6e", sy);
    }

  if (output.sz)
    {
      // sz = sqrt (<z^2> - <z>^2)

      if (!output.z_avg)
	{
	  // <z>
	  z_avg = expectation1D (params, 3, params.z, psi) * inv_norm2;
	}

      double sz = sqrt (expectation1D (params, 3, params.z2, psi) * inv_norm2
			- SQUARE(z_avg));
      
      if (params.rank == 0)
	fprintf (fp_results, " %.6e", sz);
    }
  
  /* Autocorrelation function */
  if (output.autoc)
    {
      double ac = cabs (integrate3D (params, psii, psi));
      ac *= ac;

      if (params.rank == 0)
	fprintf (fp_results, " %.6e", ac);
    }
      
  if (params.rank == 0)
    {
      fprintf (fp_results, "\n");
      fflush (fp_results);
    }

  /* User-defined observables */
  if (output.user)
    user_observe (params, t, psi);
  
  return;
}



double
energy (const parameters params, const double t, const ener_data ed,
	double * const pot, fftw_complex *psi, fftw_complex *psitmp)
{
  /** Calculates the energy of the wave function psi **/

  /* Kinetic energy */
  
  // Forward FFT to momentum space
  fftw_execute (ed.forward);

  // Apply kinetic operator in momentum space
  size_t index_1, index_2, index;
  double k1k2;

  for (size_t i = 0; i < ed.n1_local; ++i)
    {
      index_1 = i * ed.n3 * ed.n2;
      for (size_t j = 0; j < ed.n2; ++j)
	{
	  index_2 = j * ed.n3 + index_1;
	  k1k2 = ed.k1[i] + ed.k2[j];
	  for (size_t k = 0; k < ed.n3; ++k)
	    {
	      index = k + index_2;
	      psitmp[index] *= ed.k3[k] + k1k2;
	    }
	}
    }

  // Backward FFT to position space
  fftw_execute (ed.backward);   

  /* Potential energy */

  // Get potential
  potential (params, t, pot);

  for (size_t l = 0; l < params.n_local; ++l)
    {
      psitmp[l] += psi[l] * pot[l];
    }

  /*** Integrate psi* psitmp ***/
  
  // (Note: imaginary part should be zero)
  return creal (integrate3D (params, psi, psitmp));
}



double
expectation1D (const parameters params, const int dir, const double * const a, 
	       const fftw_complex * const psi)
{
  /** Calculate the expectation value of a, considering that it varies 
      only along one direction (1=x, 2=y, 3=z) **/

  size_t index_x, index_xy, index;
 
  double sum = 0., sum_total;

  if (dir == 1)
    {
      for (size_t i = 0; i < params.nx_local; ++i)
	{
	  index_x = i * params.nz * params.ny;
	  for (size_t j = 0; j < params.ny; ++j)
	    {
	      index_xy = j * params.nz + index_x;
	      for (size_t k = 0; k < params.nz; ++k)
		{
		  index = k + index_xy;
		  double psi_abs = cabs (psi[index]);
		  sum += SQUARE(psi_abs) * a[i];
		}
	    }
	}
    }
  else if (dir == 2)
    {
      for (size_t i = 0; i < params.nx_local; ++i)
	{
	  index_x = i * params.nz * params.ny;
	  for (size_t j = 0; j < params.ny; ++j)
	    {
	      index_xy = j * params.nz + index_x;
	      for (size_t k = 0; k < params.nz; ++k)
		{
		  index = k + index_xy;
		  double psi_abs = cabs (psi[index]);
		  sum += SQUARE(psi_abs) * a[j];
		}
	    }
	}
    }
  else if (dir == 3)
    {
      for (size_t i = 0; i < params.nx_local; ++i)
	{
	  index_x = i * params.nz * params.ny;
	  for (size_t j = 0; j < params.ny; ++j)
	    {
	      index_xy = j * params.nz + index_x;
	      for (size_t k = 0; k < params.nz; ++k)
		{
		  index = k + index_xy;
		  double psi_abs = cabs (psi[index]);
		  sum += SQUARE(psi_abs) * a[k];
		}
	    }
	}
    }
  else
    {
      fprintf (stderr, "Invalid argument dir = %ud in function expectation1D.\n", dir);
      abort ();
    }
  
  sum *= params.dx * params.dy * params.dz;

#ifdef MPI
  MPI_Allreduce (&sum, &sum_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  sum_total = sum;
#endif

  return sum_total;
}



double
norm (const parameters params, const fftw_complex * const psi)
{
  /** Calculates the norm of the wave function **/

  return sqrt (creal (integrate3D (params, psi, psi)));
}



void
renormalize (const parameters params, fftw_complex *psi)
{
  /** Renormalize the wave function **/

  double fact = 1. / norm (params, psi);
  
  for (size_t l = 0; l < params.n_local; ++l)
    psi[l] *= fact;

  return;
}
  


double complex
integrate3D (const parameters params, const double complex * const f1,
	     const double complex * const f2)
{
  /** Integrates f1^* f2 dtau  **/
  
  double complex sum = 0., sum_total;
  
  for (size_t l = 0; l < params.n_local; ++l)
    {
      sum += conj (f1[l]) * f2[l];
    }
  
  sum *= params.dx * params.dy * params.dz;
  
#ifdef MPI
  MPI_Allreduce (&sum, &sum_total, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  sum_total = sum;
#endif

  return sum_total;
}



/*** Utility routines ***/

void
abort ()
{
  /** Abort execution **/

  fprintf (stderr, "Execution aborted.\n");

#ifdef MPI
  MPI_Abort (MPI_COMM_WORLD, -99);
#endif

  exit(-99);
}



#ifdef MPI
void 
distribute_parameters (parameters *params, output_flags *output, 
		  unsigned int *nt_inner, unsigned int *nt_outer)
{
  /** Transfer parameters from processor 0 to others **/

  /* Since size_t may differ from one architecture to another,
     we temporarily convert size_t to unsigned long int */
  unsigned long int nx, ny, nz, n;

  if (params->rank == 0)
    {
      nx = (unsigned long int)params->nx;
      ny = (unsigned long int)params->ny;
      nz = (unsigned long int)params->nz;
      n = (unsigned long int)params->n;
    }

    MPI_Bcast (&nx, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast (&ny, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast (&nz, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast (&n, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast (&params->x_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&params->y_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&params->z_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&params->x_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&params->y_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&params->z_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&params->dx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&params->dy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&params->dz, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&params->mass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&params->dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&params->hbar, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (nt_inner, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast (nt_outer, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast (&output->yes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&output->norm, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&output->energy, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&output->x_avg, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&output->y_avg, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&output->z_avg, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&output->sx, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&output->sy, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&output->sz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&output->autoc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&output->wf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (params->rank != 0)
      {
	params->nx = (size_t)nx;
	params->ny = (size_t)ny;
	params->nz = (size_t)nz;
	params->n = (size_t)n;
      }

    return;
}
#endif



void
distribute_wf (const parameters params, double complex * const psi_in,
	       double complex * const psi_out)
{
#ifdef MPI
  /** Processor 0 sends the wave function to the other processors **/
  
  // rank 0 needs the size for each processor
  int *nx = (int *) malloc (sizeof (int) * params.size);
  if (nx == NULL)
    {
      fprintf (stderr, "Out of memory in function read_wf_bin\n");
      abort ();
    }
  
  int nx_local = 2 * params.nx_local; // (factor of 2 because psi is complex)
  
  MPI_Gather (&nx_local, 1, MPI_INT, nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  // Calculate displacement vector
  int *displs = (int *) malloc (sizeof (int) * params.size);
  if (displs == NULL)
    {
      fprintf (stderr, "Out of memory in function read_wf_bin\n");
      abort ();
    }
  
  if (params.rank == 0)
    {
      displs[0] = 0;
      for (int i = 1; i < params.size; ++i)
	{
	  displs[i] = displs[i-1] + nx[i-1];
	}
    }
  
  // Send information to all
  MPI_Scatterv (psi_in, nx, displs, MPI_DOUBLE, psi_out, nx_local, MPI_DOUBLE, 
		0, MPI_COMM_WORLD);
  
  free (nx);
  free (displs);
#else
  /* This function should be called only in the parallel version. 
     In case it is called from the serial version, it only copies psi_in 
     to psi_out.
  */
  
  for (size_t l = 0; l < params.n; ++l)
    {
      psi_out[l] = psi_in[l];
    }
#endif
  
  return;
}
