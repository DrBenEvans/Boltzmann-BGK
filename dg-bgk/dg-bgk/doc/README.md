# Introduction

This documentation is only a collection of technical information that I 
found useful to implement the optimization I needed, and is provided to future 
generations in order to facilitate the process of understanding how the code
works.

My work has been mainly related to the optimization of the multi process aspect
of the program, through MPI. Some notes can be found 
[here](../doc/mpi_structure.md).

In what follows, the claims made on lines marked with (?) would require some 
verification.



# Random notes, but important

*  The program requires intel fortran to be compiled and work with the input 
data provided (see file [example.con](../doc/example.con)). As an example, 
gfortran would need quotes around strings in the configuration file. This does 
not exclude that the program can be compiled successfully with gnu tools and 
execute correctly - I just cannot guarantee that (also due to my ignorance). 
* The point indices for element definition in the mesh file are given in
clockwise order, otherwise the algorithm used to fill **ISIDE** would not 
work.
* The mesh was partitioned starting from a node-base partitioning file, and 
elements were assigned to partitions using an arbitrary rule, that is - the
element is assigned to the partition, among the partitions to which its nodes
are assigned, with the highest index. It is also possible to use native 
element-based partitioning and that would remove the need to do this 
potentially suboptimal step.
* In the following, notice the difference between *edge* and *boundary* 



# Info on variables and arrays (on the master rank)

* **NELEM**: Number of (triangular) elements in the mesh
* **NPOIN**: Number of nodes in the mesh
* **NBOUN**: Number of edges in the mesh *on the boundary*.
* **MXSID**: Maximum number of sides. In the code is set to **2x(NELEM+NBOUN)**, 
  but I think that **(3xNELEM+NBOUN)/2** should be enough.
* **NSIDE**: Actualy number of edges in the mesh (I think it should be 
(3xNELEM+NBOUN)/2).

In what follows, a 2D mesh is assumed.

* **INTMA(3,NELEM)** : 
   "Element array" (in the opengl sense). Contains the point indices (integers) 
   for the 3 nodes of every element in the mesh.
                 
* **BSIDO(4+NDIMN,NBOUN)**: 
   Boundary data. NDIMN is actually 2, so this is a 6xNBOUN 2D array 
   of integers. For each row, we have:
    1. Point index
    2. Second point index
    3. Element index
    4. Type of edge (1 inflow, 4 object surface, 2 (?) ).
    5. (?)
    6. (?)

* **RSIDO(1+NDIMN,NBOUN)**:
   Other boundary data. NDIMN is acually 2. Let **n** be the normal to the 
   edge and **l** the lenght of the edge, then for each row we have:
    1. n<sub>x</sub> l <sup>2</sup>
    2. n<sub>y</sub> l <sup>2</sup>
    3. 1/l

* **LHOWM(NPOIN)**: How many elements join on a given node
* **LWHERE(NPOIN)**: A scan of LHOWM (that is, LWHERE(i) is the sum of
  LHOWM(j) for j < i, and as a consequence LWHERE(1)=0)

* **ICONE(3xNELEM)**: Index of the element associated with every *duplicated* 
   node (a node *i* is duplicated LHOWM(i) times).

* **ISIDE(8,MXSID)**: Array containing informations on sides. 
    1. First node index (IP1).
    2. Second node index (IP2).
    3. Index of the element on the left of the edge (as seen from the first 
        node looking to the second node). *Cannot be zero*.
    4. Index of the element on the right of the edge (same convention as above)
    5. First point index in the left element (between 1 and 3 included)
    6. Second point index in the left element (between 1 and 3 included)
    7. First point index in the right element (between 1 and 3 included)
    8. Second point index in the right element (between 1 and 3 included)


    Some remarks:
    * The element index in ISIDE(4) can be zero if the edge is on the border of 
the whole domain.
    * The first node index is always smaller than the second node index 
*except on the boundaries*, where the two points are flipped so that only
ISIDE(4) can be zero.
    * When ISIDE(4) is zero, then also ISIDE(7) and ISIDE(8) are.
    * MXSID is set to 2*(NELEM+NBOUN), but I think it can be proven that 
      (3*NELEM+NBOUN)/2 is enough. Anyway, oversizing ISIDE is not a problem.

* **NX(NSIDE)**: Contains the *x* components of the normal to each edge in the
  mesh. The normals are taken 90 degrees clockwise with respect to the 
direction of the edge.
* **NY(NSIDE)**: Contains the *y* components of the normal to each edge in the
  mesh. The normals are taken 90 degrees clockwise with respect to the 
direction of the edge.


### Macroscopic quantities
* **RHO(NPOIN)**: Density for each point.
* **UVEL(NPOIN)**: Average x-velocity for each point.
* **VVEL(NPOIN)**: Average y-velocity for each point.
* **PS(NPOIN)**: Pressure for each point (as the average of the pressure tensor).
* **TEMP(NPOIN)**: Temperature for each point.

Notice that these quantities are somehow averaged between the various copies of
the nodes belonging to each element the node is part of (because it's 
Discontinuous Galerkin, right?)

### Not fully understood stuff

* **IELSI(2,NELEM)**: (?) Integer values.
* **COORD(2,NPOINT)**: Coordinates of the nodes.
* **GEOME(7,NELEM)**: (?) Real values. Something related to the normals for the 
3 edges of the element, multiplied by l/2A 
* **MMAT(3,NELEM)**: Either this or CMAT is used. Related with GEOME(7,...)
* **CMAT(3,NELEM)**: Either this or MMAT is used. Related with GEOME(1,...)

# Multi-rank related 
In what follows, the words *rank*, *partition*, *process* and *group* might be 
used interchangeably depending on my current mood. Their real meaning is of
course tightly connected.
## Variables computed by master and communicated to slaves (not specific to each slave)

* **NGRPS**: Number of slaves (groups of points).

* **ELGRP(NELEM,2)**: A dictionary for element indices between master and slaves.
   It is computed in the [PTMESH](../src/PTMESH.f90)
   subroutine on the master rank and broadcast to slaves.
    1. Partition (group) the element belongs to.
    2. Index IEG of the element in its partition ( 1 <= IEG <= NELEM_PP )

* **NEGRP(NGRPS)**: Number of elements in each partition (group). Each nuber in
  this array is sent to the respective rank as NELEM_PP.
* **GCOMM(NGRPS,NGRPS)**: Symmetric matrix that contains the count of the edges
  shared between two partitions (groups). Note: GCOMM(i,i) = 0. As for ELGRP,
this is computed by the master rank in the [PTMESH](../src/PTMESH.f90) routine 
and broadcasted.

About MPI-related stuff (see also [mpi_structure.md](mpi_structure.md)):

* **MPI_COMM_P**: MPI communicator which is common to all the processes sharing
  the same portion velocity space (P-space-spanning communicator).
* **MPI_RANK_P**: rank number for the current process associated to the 
MPI_COMM_P communicator.
* **MPI_SIZE_P**: Number of ranks in each P-space-spanning communicator.
 It's equal to the number of partitions in position space.

* **MPI_COMM_V**: MPI communicator which is common to all the processes sharing
  the same portion of position space (V_space-spanning communicator).
* **MPI_RANK_V**: rank number for the current process associated to the 
MPI_COMM_V communicator.
* **MPI_SIZE_V**: Number of ranks in each V-space-spanning communicator. 
**NOTE:** Because of reflections of particles on the object surface, at each 
step of the algorithm some particles bounce across velocity space. Fortunately, 
energy wrt the object surface is conserved, which allows us to claim that the 
modulus of speed does not change in the collision. This allows us not to 
modify the [GETBOU](../src/GETBOU.f90) subroutine, as long as the portion of the velocity
space with fixed **|v|** is not split across ranks. This is guaranteed only if 
MPI_SIZE_V is a divisor of **RORDER**: for this reason, this constraint is 
asserted at runtime (i.e. the program will crash if this is not true). See also
the example [example.con](./example.con)


## Variables computed by master and communicated to slaves (specific to each slave)
Usually the variables which end with "**_PP**" fall in this cathegory. 
Some notes:
* The useful part of arrays declared with a lenght of maxNPOIN_PP is only 
NPOIN_PP-big.
* The useful part of arrays declared with a lenght of maxNELEM_PP is only 
NELEM_PP-big, that is NEGRP(IG)-big.
* *Mutatis mutandis*, this holds for maxNBOUN_PP, and similar cases. 

Now to the list of variables and arrays:

* **NELEM_PP**: Number of elements in a partition (process specific).
* **NCOMM_PP**: Number of edges shared with all the other ranks-processes.

* **IPCOM_PP(maxNPOIN_PP)**: Dictionary between local (process-specific) node 
index and global node index.
 The local point index is obtained incrementing a counter in a 
cycle over the elements assigned to that process-partition. 

* **INTMA_PP(3,maxNELEM_PP)**: Like INTMA on the master rank, but using the
  IPCOM_PP dictionary to translate node indices and the ELGRP dictionary to
translate element indices.
* **COORD_PP(2,maxNPOIN_PP)**: Translation of COORD using the IPCOM_PP
  dictionary for the point indices. 

* **IBCOM_PP(maxNBOUND_PP)**: Dictionary between local *boundary* index and global
  *boundary* index, obtained by cycling on all the sides globally  and 
incrementing a counter when the side is in the partition considered.

* **BSIDO_PP(6,maxNBOUN_PP)**: 
    1. Local point 1 index, translated from global indices using IPCOM_PP
    2. Local point 2 index, translated from global indices using IPCOM_PP
    3. Local element index, translated from global indices using ELGRP

    4,5,6 are copied from BSIDO as they are. The border index in BSIDO_PP is
related to the one in BSIDO by the dictionary IBCOM_PP.

* **RSIDO_PP(3,maxNBOUN_PP)**: Content is copied as is from the global RSIDO,
  and the border index is related to the one in RSIDO by the dictionary
IBCOM_PP as in the BSIDO_PP/BSIDO case.

* **ISCOM_PP(maxNSIDE_PP)**: Dictionary between local and global *edge* index.

* **ISIDE_PP(8,maxNSIDE_PP)**: Translation of ISIDE. Point indices 1 and 2 are
  translated using IPCOM_PP, element indices 3 and 4 are translated using
ELGRP(...,2), and the 5,6,7,8 indices are copied as they are. 

* **SDCOM_PP(3,NCOM_PP)**: Informations about shared edges. This array is set up
  in the EDGCOM routine.
    1. Local edge index (in the range from 1 to NSIDE_PP)
    2. Local edge index on opposite rank
    3. Opposite rank 

* **LCOMM_PP(maxNSIDE_PP)**: This array states if an edge is shared or not. If
  edge *i* is shared, LCOMM_PP(i) is equal to the other rank that shares the
edge, otherwise (that is, in the case of an internal edge) it is equal to zero.

* **DISNF_PP(NNODE,VNPNT_PART*MPI_RANK_V+1:VNPNT_PART*(MPI_RANK_V+1),maxNELEM_PP)**:
    The density function discretized in position and velocity space. In order
to implement velocity space partitioning, the second index is chosen to go from
VNPNT_PART*MPI_RANK_V+1 to VNPNT_PART*(MPI_RANK_V+1). 

### Not fully understood stuff

The content of **GEOME_PP**, **MMAT_PP** (or **CMAT_PP**) are copied from the
global arrays, using the ELGRP(...,2) dictionary to translate between global and
local element index.

