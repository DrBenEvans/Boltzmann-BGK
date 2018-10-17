# MPI IMPLEMENTATION NOTES

In the original implementation, no MPI was used. A domain decomposition with
MPI was implemented later to tackle memory issues: the goal was *not*
performance, but only be able to solve a larger problem. 

A master-slave configuration was chosen, where the master would do the task
relative to IO and to the domain decomposition logic as, e.g., splitting the
mesh and set up the lookup tables for all the processes - see [main
README](README.md).

### SOME COMMENTS ON THE OLD MPI IMPLEMENTATION

In the original MPI implementation, in most cases a very inefficient
communication pattern was choosen, like the one descripted in the pseudocode
snippets below.

Here is one pattern used in initialization:
```c
for all points in the mesh 
    if on the master rank:
        find the process IG the point belongs to
    
    MPI_BCAST of IG to all processes

    if on the master process:      MPI_SEND to the IG-th process
    if on the IG-th process:       MPI_RECV from the master process

end for on all points in the mesh 
```

Another variant, for inter process comunication of border data:

```c
for IG in all processes 
    for all points in the partition IG
        if on the IG-th process:
            find the process IG2 to which data must be sent

        MPI_BCAST of IG2 to all other processes

    if on the IG-th process:     MPI_SEND to IG2-th process
    if on the IG2-th process:    MPI_RECV from the IG-th process


    end for on all points in the partition
end for on all processes
```

This is computationally inefficient, since for each step of the time marching 
scheme millions of MPI_BCAST calls are performed. 
This scheme has been modified only where it would severely decrease performance
(i.e., during the time steps and in the [OUTPUT](../src/OUTPUT.f90) subroutine).

## VELOCITY SPACE PARTITIONING.

In order to do velocity space partitioning, the MPI communicator layout is the
one depicted in this figure: ![communicator_layout.](communicators_layout.svg)
(At the moment on Gitlab the svg is not shown in the text but you can
visualise it separately).

Some remarks:
* Velocity space has been partitioned replicating the master-slave layout in
  the original code. In order to remove the performance penalty associated to 
  the lazy master process, a share of the work of the slaves has been assigned
  to the master.

List of relevant variables:
* **VSPACE_FIRST**: First index of vspace that is in the responsibility of the
  current process.
* **VSPACE_LAST**: Last index of vspace that is in the responsibility of the
  current process.
* **VCORD**: the *u* and *v* coordinates associated to every point (labeled,
  e.g. IV) in the velocity space. It has size 1:VNPNT on every rank instead of
just the minimum possible which would be VSPACE_FIRS:VSPACE_LAST (it does not 
contain a big amount of data, by the way).

As stated in [README](./README.md), there are some constraints on the choice of 
the number of partitions in velocity space (see discussion around MPI_SIZE_V).

### OUTPUT 

The **OUTPUT** subroutine prints to file the content of the DISNF array, which
contains the particle density in che cartesian product of partition and
velocity space. In this subroutine, a communication pattern described in the
[previous section](#some-comments-on-the-old-mpi-implementation) was used.
In the new implementation, the subroutine consists of the following steps:
* All the "master" ranks (having MPI_RANK_P equal to 0) collect and reorder the
  data pertaining to their velocity space partition from the ranks in the
same MPI_COMM_P communicator:
    *  Each rank sends to the master in its MPI_COMM_P its whole DISNF_PP
       array (which has size 3x[VSPACE_FIRST:VSPACE_LAST]xNELEM_PP).
    *  The master ranks in every communicator MPI_COMMP_P make use of the ELGRP
       dictionary to dispatch the content of each DISNF_PP received to the
right location in DISNF (which has size  3x[VSPACE_FIRST:VSPACE_LAST]xNELEM).
* The master rank in MPI_COMM_WORLD (that is also the master rank in the
  MPI_COMM_V communicator between the ranks having MPI_RANK_P equal to zero)
collects the DISNF arrays coming from each other rank in the MPI_COMM_V
communicator and writes the content to a file.

### NOTES

* The **INICON** and **INICON2** subroutines have been modified for VSPACE
  partitioning, but the inefficient communication pattern was kept in
place, since this is not a performance-critical part of the program.
