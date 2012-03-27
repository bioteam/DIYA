#!/bin/bash
# Run MPI-enabled cmsearch using LSF
#BSUB -L /bin/sh

ARGS=3
NP=100
MPIRUN=/usr/mpi/gcc/openmpi-1.4.3/bin/mpirun
CMSEARCH=/Jake/apps/infernal-mpi/bin/cmsearch
WORKDIR=/home/briano/tmp

if [ $# -ne $ARGS ]  # Correct number of arguments passed to script?
then
  echo "Usage: $0 <*.cm file> <sequence file> <OUTPUTFILE>"
  exit
fi

# module load infernal-mpi

CMFILE=$1      # For example, /Jake/data/mirror/Rfam-10.0/sRNAs.cm
SEQFILE=$2     # For example, test0001.fa
OUTPUTFILE=$3

$MPIRUN -np $NP $CMSEARCH --mpi --tabfile $OUTPUTFILE -E 1.0e-20 $CMFILE $WORKDIR/$SEQFILE

# Example mpirun command:
# /usr/mpi/gcc/openmpi-1.4.3/bin/mpirun -np 4 /Jake/apps/infernal-mpi/bin/cmsearch --mpi --tabfile OUTPUTFILE -E 1.0e-20 -o test.out /Jake/data/mirror/Rfam-10.0/sRNAs.cm test0001.fa
# Using bsub:
# bsub -e err -o out -J cmsearch "/Jake/apps/DIYA/scripts/run-cmsearch-mpi.sh /Jake/data/mirror/Rfam-10.0/sRNAs.cm test0001.fa OUTPUTFILE"