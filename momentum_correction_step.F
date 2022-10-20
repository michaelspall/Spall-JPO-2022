C $Header: /u/gcmpack/MITgcm/model/src/momentum_correction_step.F,v 1.10 2011/11/13 01:44:14 jmc Exp $
C $Name:  $

#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: MOMENTUM_CORRECTION_STEP
C     !INTERFACE:
      SUBROUTINE MOMENTUM_CORRECTION_STEP( myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE MOMENTUM_CORRECTION_STEP
C     *==========================================================*
C     |1rst Part : Update U,V.
C     |
C     | The arrays used for time stepping are cycled.
C     | Momentum:
C     |           V(n) = Gv(n) - dt * grad Eta
C     |
C     |part1: update U,V
C     |  U*,V* (contained in gU,gV) have the surface
C     |     pressure gradient term added and the result stored
C     |     in U,V (contained in uVel, vVel)
C     |
C     |part2: Adjustments
C     |   o Filter  U,V (Shapiro Filter, Zonal_Filter)
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"

#ifdef ALLOW_SHAP_FILT
#include "SHAP_FILT.h"
#endif
#ifdef ALLOW_ZONAL_FILT
#include "ZONAL_FILT.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myTime :: Current time in simulation
C     myIter :: Current iteration number in simulation
C     myThid :: my Thread Id. number
      _RL myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables
      _RL phiSurfX(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL phiSurfY(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gUtmp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gVtmp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin,iMax
      INTEGER jMin,jMax
      INTEGER bi,bj
      INTEGER k,i,j
CEOP

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

C--     Set up work arrays that need valid initial values
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          phiSurfX(i,j)=0.
          phiSurfY(i,j)=0.
         ENDDO
        ENDDO

C       Loop range: Gradients of Eta are evaluated so valid
C       range is all but first row and column in overlaps.
        iMin = 1-OLx+1
        iMax = sNx+OLx
        jMin = 1-OLy+1
        jMax = sNy+OLy

C-      Calculate gradient of surface Potentiel
        CALL CALC_GRAD_PHI_SURF(
     I       bi,bj,iMin,iMax,jMin,jMax,
     I       etaN,
     O       phiSurfX,phiSurfY,
     I       myThid )

        IF (momStepping) THEN
C--       Loop over all layers, top to bottom
          DO k=1,Nr
         DO j=jMin,jMax
         DO i=iMin,iMax
          gUtmp(i,j) = gU(i,j,k,bi,bj)
          gVtmp(i,j) = gV(i,j,k,bi,bj)
         ENDDO
        ENDDO

C-        Update velocity fields:  V(n) = V** - dt * grad Eta
            CALL CORRECTION_STEP(
     I         bi,bj,iMin,iMax,jMin,jMax,k,
     I         phiSurfX,phiSurfY,myTime,myThid )

c MAS
         DO j=jMin,jMax
         DO i=iMin,iMax
        gUtmp(i,j) = (uVel(i,j,k,bi,bj)-gUtmp(i,j))/deltaTmom
        gVtmp(i,j) = (vVel(i,j,k,bi,bj)-gVtmp(i,j))/deltaTmom
         ENDDO
        ENDDO
#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics  )
     &    CALL DIAGNOSTICS_FILL(gUtmp,'UDIAG1  ',k,1,2,bi,bj,myThid)
#endif /* ALLOW_DIAGNOSTICS */
#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics  )
     &    CALL DIAGNOSTICS_FILL(gVtmp,'UDIAG2  ',k,1,2,bi,bj,myThid)
#endif /* ALLOW_DIAGNOSTICS */

C--       End DO k=1,Nr
          ENDDO

        ENDIF

#ifdef ALLOW_OBCS
        IF (useOBCS) THEN
          CALL OBCS_APPLY_UV( bi, bj, 0, uVel, vVel, myThid )
        ENDIF
#endif /* ALLOW_OBCS */

C--    End of 1rst bi,bj loop
       ENDDO
      ENDDO

C--- 2nd Part : Adjustment.

C--   Filter (and exchange)
#ifdef ALLOW_SHAP_FILT
      IF (useSHAP_FILT) THEN
       IF ( .NOT.shap_filt_uvStar ) THEN
        CALL TIMER_START('SHAP_FILT_UV       [MOM_CORR_STEP]',myThid)
        CALL SHAP_FILT_APPLY_UV( uVel, vVel, myTime, myIter, myThid )
        CALL TIMER_STOP ('SHAP_FILT_UV       [MOM_CORR_STEP]',myThid)
       ENDIF
      ENDIF
#endif
#ifdef ALLOW_ZONAL_FILT
      IF (useZONAL_FILT) THEN
       IF ( .NOT.zonal_filt_uvStar ) THEN
        CALL TIMER_START('ZONAL_FILT_UV      [MOM_CORR_STEP]',myThid)
        CALL ZONAL_FILT_APPLY_UV( uVel, vVel, myThid )
        CALL TIMER_STOP ('ZONAL_FILT_UV      [MOM_CORR_STEP]',myThid)
       ENDIF
      ENDIF
#endif

      RETURN
      END
