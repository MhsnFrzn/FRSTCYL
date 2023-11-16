c  FLTIMER = whether to do timing
c  FLAVER  = whether to do axial average
c  FLTHRM  = whether to include temperature
c  FLTHRM2 = whether to include concentration 
c  FLSTLE  = whether tt has settling velocity
c  FLGVARY = whether gravity vector changes with y 
c  FLPGVAR = whether pressure gradient changes with time 
c  FLIBM   = whether to account for immersed boundary method 
c  FLMOV   = whether immersed boundary can move 
c  FLSURF  = whether to read IBM surface 

      LOGICAL FLTIMER,FLAVER,FLTHRM,FLSTLE,FLTHRM2,
     *        FLSV_YES,FLSV_NO,FLGVARY,FLPGVAR,
     *        FLSV2_YES,FLSV2_NO,FLIBM,FLMOV,FLSURF
      COMMON /flags/ FLTIMER,FLAVER,FLTHRM,FLTHRM2,
     *               FLSTLE,FLIBM,FLMOV,FLSURF,
     *               FLSV_YES,FLSV_NO,
     *               FLSV2_YES,FLSV2_NO,
     *               FLGVARY,FLPGVAR

