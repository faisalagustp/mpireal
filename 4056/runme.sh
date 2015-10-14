#!/bin/bash

mpicc -o vv1 vv_mult_blkstp_unf.c;
mpicc -o vv2 vv_mult_blk_cyclic.c;
mpicc -o mv1 mv_mult_checkerboard.c;
mpicc -o mv2 mv_mult_blkstp.c;
mpicc -o mm1 mpi_mm.c;
mpicc -o mm2 mm_mult_fox.c;

sbatch vv1.8;
sbatch vv1.16;
sbatch vv1.32;
sbatch vv1.64;
sbatch vv1.128;
sbatch vv1.256;

sbatch vv2.8;
sbatch vv2.16;
sbatch vv2.32;
sbatch vv2.64;
sbatch vv2.128;
sbatch vv2.256;

sbatch mv1.8;
sbatch mv1.16;
sbatch mv1.32;
sbatch mv1.64;
sbatch mv1.128;
sbatch mv1.256;

sbatch mv2.8;
sbatch mv2.16;
sbatch mv2.32;
sbatch mv2.64;
sbatch mv2.128;
sbatch mv2.256;

sbatch mm1.8;
sbatch mm1.16;
sbatch mm1.32;
sbatch mm1.64;
sbatch mm1.128;
sbatch mm1.256;

sbatch mm2.8;
sbatch mm2.16;
sbatch mm2.32;
sbatch mm2.64;
sbatch mm2.128;
sbatch mm2.256;
