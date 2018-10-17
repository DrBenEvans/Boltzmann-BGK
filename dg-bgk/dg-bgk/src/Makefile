all = BGK

FFLAGS= -fc=ifort -heap-arrays -g -traceback -check all 
#FFLAGS= -fc=ifort -heap-arrays -check all 
#FFLAGS= -fc=ifort -heap-arrays -O2

OBJ = ADDVEC.o ADTHEM.o ADVNCE.o ALOTIM2.o ALOTIM.o BCASTP.o BCASTV.o CMMINC.o COMPMU.o CONMAT.o DG_lin_convec.o EDGCOM.o EDGFLX.o EDGFLXOPT.o ELNODE.o FLUCON2.o GETAVR.o GETBOU.o GETELC.o GETFLA.o GETFLU.o GETFOR.o GETGEO.o GETINC.o GETINF.o GETLEN.o GETLOC.o GETMAC.o GETMAT.o GETNOR.o GETRES.o GETRHO.o GETROW.o GETSID.o GTEQNF.o GTINPT.o IDENTA.o IDENTM.o IFILLA.o IFILLM.o IFILLV.o INFLOW.o INICON2.o INICON.o MACROS.o normal.o OUTPUT.o PTMESH.o RFILLA.o RFILLM.o RFILLV.o spectralBGK.o VSPACE.o

%.o : %.f90 
	mpif90 $(FFLAGS) -c $<

BGK : $(OBJ)
	mpif90 $(FFLAGS) -o $@ $^

clean:
	rm BGK $(OBJ)
