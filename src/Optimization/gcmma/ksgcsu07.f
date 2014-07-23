C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE MMASUG(ITER,M,N,GEPS,IYFREE,XVAL,XMMA,
     1                  XMIN,XMAX,XLOW,XUPP,ALFA,BETA,
     2                  A,B,C,Y,Z,RAA0,RAA,ULAM,
     3                  F0VAL,FVAL,F0APP,FAPP,FMAX,DF0DX,DFDX,
     4                  P,Q,P0,Q0,UU,GRADF,DSRCH,HESSF)
C
C       Version "August 2007".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     MMASUG generates and solves the MMA subproblem within
C     the globally convergent version of MMA.
C     XMMA,Y,Z = optimal solution of the generated subproblem.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION XVAL(1),XMMA(1),XMIN(1),XMAX(1),XLOW(1),XUPP(1),
     1          ALFA(1),BETA(1),A(1),B(1),C(1),Y(1),RAA(1),ULAM(1),
     2          FVAL(1),FAPP(1),FMAX(1),DF0DX(1),DFDX(1),
     3          P(1),Q(1),P0(1),Q0(1),UU(1),GRADF(1),DSRCH(1),
     4          HESSF(1)
      INTEGER IYFREE(1)
C
C********+*********+*********+*********+*********+*********+*********+
C
      CALL GENSUG(M,N,XVAL,XMIN,XMAX,F0VAL,DF0DX,FMAX,FVAL,
     1            DFDX,P,Q,B,P0,Q0,R0,XLOW,XUPP,RAA0,RAA)
C
C***** GENSUG generates the MMA subproblem by calculating the
C***** coefficients P(i,j),Q(i,j),B(i),P0(j),Q0(j) and R0.
C     
      CALL MAXIM(M,N,GEPS,IYFREE,GRADF,DSRCH,HESSF,XMMA,Y,Z,
     1           ULAM,UU,XLOW,XUPP,ALFA,BETA,A,B,C,P,Q,P0,Q0)
C
C***** MAXIM solves the dual problem of the MMA subproblem.
C***** ULAM = optimal solution of this dual problem.
C***** XMMA,Y,Z = optimal solution of the MMA subproblem.
C
      CALL APPRF(M,N,XMMA,P,Q,B,P0,Q0,R0,XLOW,XUPP,FMAX,F0APP,FAPP)
C
C***** F0APP, FAPP = values of the approximating functions at XMMA.
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE ASYMPG(ITER,M,N,XVAL,XMIN,XMAX,XOLD1,XOLD2,
     1                  XLOW,XUPP,ALFA,BETA)
C
C       Version "August 2007".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     ASYMPG calculates the asymptotes XLOW and XUPP,
C     and the bounds ALFA and BETA, for the current subproblem.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XVAL(1),XMIN(1),XMAX(1),XOLD1(1),XOLD2(1),
     1          XLOW(1),XUPP(1),ALFA(1),BETA(1)
C
      ALBEFA=0.1
      GHINIT=0.5
      GHDECR=0.7
      GHINCR=1.2
      IF(ITER.GE.3) GOTO 350
C
C***  Here ITER = 1 or 2 .
      DO 200 J=1,N
      XMMJ=XMAX(J)-XMIN(J)
      IF(XMMJ.LT.0.00001) XMMJ=0.00001
      XLOW(J)=XVAL(J)-GHINIT*XMMJ
      XUPP(J)=XVAL(J)+GHINIT*XMMJ
  200 CONTINUE
      GOTO 500
C
C***  Here ITER is greater than 2.
  350 CONTINUE
C
      DO 400 J=1,N
      XTEST=(XVAL(J)-XOLD1(J))*(XOLD1(J)-XOLD2(J))
      FAK=1.0
      IF(XTEST.LT.0.) FAK=GHDECR
      IF(XTEST.GT.0.) FAK=GHINCR
      XLOW(J)=XVAL(J)-FAK*(XOLD1(J)-XLOW(J))
      XUPP(J)=XVAL(J)+FAK*(XUPP(J)-XOLD1(J))
      XMMJ=XMAX(J)-XMIN(J)
      IF(XMMJ.LT.0.00001) XMMJ=0.00001
      GMINJ = XVAL(J)-10.0*XMMJ
      GMAXJ = XVAL(J)-0.01*XMMJ
      HMINJ = XVAL(J)+0.01*XMMJ
      HMAXJ = XVAL(J)+10.0*XMMJ
      IF(XLOW(J).LT.GMINJ) XLOW(J)=GMINJ
      IF(XLOW(J).GT.GMAXJ) XLOW(J)=GMAXJ
      IF(XUPP(J).LT.HMINJ) XUPP(J)=HMINJ
      IF(XUPP(J).GT.HMAXJ) XUPP(J)=HMAXJ
  400 CONTINUE
C
  500 CONTINUE
C
      DO 600 J=1,N
      XMIJ=XMIN(J)-0.000001
      XMAJ=XMAX(J)+0.000001
      IF(XVAL(J).GE.XMIJ) GOTO 550
      XLOW(J)=XVAL(J)-(XMAJ-XVAL(J))/0.9
      XUPP(J)=XVAL(J)+(XMAJ-XVAL(J))/0.9
      GOTO 600
  550 CONTINUE
      IF(XVAL(J).LE.XMAJ) GOTO 600
      XLOW(J)=XVAL(J)-(XVAL(J)-XMIJ)/0.9
      XUPP(J)=XVAL(J)+(XVAL(J)-XMIJ)/0.9
  600 CONTINUE
C
      DO 700 J=1,N
      ALFA(J)=XLOW(J)+ALBEFA*(XVAL(J)-XLOW(J))
      BETA(J)=XUPP(J)-ALBEFA*(XUPP(J)-XVAL(J))
      IF(ALFA(J).LT.XMIN(J)) ALFA(J)=XMIN(J)
      IF(BETA(J).GT.XMAX(J)) BETA(J)=XMAX(J)
  700 CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE GENSUG(M,N,XVAL,XMIN,XMAX,F0VAL,DF0DX,FMAX,FVAL,
     1                  DFDX,P,Q,B,P0,Q0,R0,XLOW,XUPP,RAA0,RAA)
C
C       Version "August 2007".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     GENSUG calculates P( ),Q( ),B( ),P0( ),Q0( ) and R0
C     for the current subproblem.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XVAL(1),XMIN(1),XMAX(1),DF0DX(1),FMAX(1),FVAL(1),
     1          DFDX(1),P(1),Q(1),B(1),P0(1),Q0(1),XLOW(1),XUPP(1),
     2          RAA(1)
C
      R0=F0VAL
      DO 20 I=1,M
      B(I)=FMAX(I)-FVAL(I)
   20 CONTINUE
C
      DO 50 J=1,N
      MJ=M*(J-1)
      UJLJ=XUPP(J)-XLOW(J)
      UJXJ=XUPP(J)-XVAL(J)
      XJLJ=XVAL(J)-XLOW(J)
      UJXJ2=UJXJ*UJXJ
      XJLJ2=XJLJ*XJLJ
      XMMJ=XMAX(J)-XMIN(J)
      IF(XMMJ.LT.0.00001) XMMJ=0.00001
      P0J=RAA0/XMMJ
      Q0J=RAA0/XMMJ
      IF(DF0DX(J).GT.0.) P0J=P0J+1.001*DF0DX(J)
      IF(DF0DX(J).GT.0.) Q0J=Q0J+0.001*DF0DX(J)
      IF(DF0DX(J).LT.0.) Q0J=Q0J-1.001*DF0DX(J)
      IF(DF0DX(J).LT.0.) P0J=P0J-0.001*DF0DX(J)
      P0J=P0J*UJXJ2
      Q0J=Q0J*XJLJ2
      P0(J)=P0J
      Q0(J)=Q0J
      R0=R0-P0J/UJXJ-Q0J/XJLJ
C
      DO 40 I=1,M
      IJ=MJ+I
      PIJ=RAA(I)/XMMJ
      QIJ=RAA(I)/XMMJ
      DFIJ=DFDX(IJ)
      IF(DFIJ.GT.0.) PIJ=PIJ+1.001*DFIJ
      IF(DFIJ.GT.0.) QIJ=QIJ+0.001*DFIJ
      IF(DFIJ.LT.0.) QIJ=QIJ-1.001*DFIJ
      IF(DFIJ.LT.0.) PIJ=PIJ-0.001*DFIJ
      PIJ=PIJ*UJXJ2
      QIJ=QIJ*XJLJ2
      P(IJ)=PIJ
      Q(IJ)=QIJ
      B(I)=B(I)+PIJ/UJXJ+QIJ/XJLJ
C
   40 CONTINUE
   50 CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE RAASTA(M,N,RAA0,RAA,XMIN,XMAX,DF0DX,DFDX)
C
C       Version "July 2007".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     RAASTA calculates starting values on RAA0 and RAA for the
C     current outer iteration in the globally convergent version
C     of MMA.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RAA(1),XMIN(1),XMAX(1),DF0DX(1),DFDX(1)
C
      ZN=N
      RAAMIN=0.000001
      RAA0=0.
      DO 10 I=1,M
      RAA(I)=0.
 10   CONTINUE
C
      DO 30 J=1,N
      MJ=M*(J-1)
      XMMJ=XMAX(J)-XMIN(J)
      IF(XMMJ.LT.0.00001) XMMJ=0.00001
      XMMJN = 0.1*XMMJ/ZN
      IF(DF0DX(J).GT.0.) RAA0=RAA0+XMMJN*DF0DX(J)
      IF(DF0DX(J).LT.0.) RAA0=RAA0-XMMJN*DF0DX(J)
      DO 20 I=1,M
      IJ=MJ+I
      DFIJ=DFDX(IJ)
      IF(DFIJ.GT.0.) RAA(I)=RAA(I)+XMMJN*DFIJ
      IF(DFIJ.LT.0.) RAA(I)=RAA(I)-XMMJN*DFIJ
 20   CONTINUE
 30   CONTINUE
C
      IF(RAA0.LT.RAAMIN) RAA0=RAAMIN
      DO 40 I=1,M
      IF(RAA(I).LT.RAAMIN) RAA(I)=RAAMIN
 40   CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE RAAUPD(M,N,GEPS,XMMA,XVAL,XMIN,XMAX,XLOW,XUPP,
     1                  F0NEW,FNEW,F0APP,FAPP,RAA0,RAA)
C
C       Version "July 2007".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     RAAUPD updates the values on RAA0 and RAA.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XMMA(1),XVAL(1),XMIN(1),XMAX(1),XLOW(1),XUPP(1),
     1          FNEW(1),FAPP(1),RAA(1)
C
      RCOMIN = 0.01*GEPS
      RAACOF = 0.
      DO 15 J=1,N
      UXXL = (XUPP(J)-XMMA(J))*(XMMA(J)-XLOW(J))
      XXXX = (XMMA(J)-XVAL(J))*(XMMA(J)-XVAL(J))
      UJLJ=XUPP(J)-XLOW(J)
      XMMJ=XMAX(J)-XMIN(J)
      IF(XMMJ.LT.0.00001) XMMJ=0.00001
      ULXX=UJLJ/XMMJ
      RAACOF = RAACOF + ULXX*XXXX/UXXL
 15   CONTINUE
C
      IF(RAACOF.LT.RCOMIN) RAACOF=RCOMIN
C
      F0APPE = F0APP + 0.5*GEPS
      IF(F0NEW.LE.F0APPE) GOTO 25
      DELTA0 = (F0NEW-F0APP)/RAACOF
      RAANE0 = 1.1*(RAA0 + DELTA0)
      RAA10 = 10.*RAA0
      IF(RAANE0.GT.RAA10) RAANE0 = RAA10
      RAA0 = RAANE0
 25   CONTINUE
C
      DO 80 I=1,M
      FAPPE=FAPP(I) + 0.5*GEPS
      IF(FNEW(I).LE.FAPPE) GOTO 75
      DELTAI = (FNEW(I)-FAPP(I))/RAACOF
      RAANEI = 1.1*(RAA(I) + DELTAI)
      RAA10 = 10.*RAA(I)
      IF(RAANEI.GT.RAA10) RAANEI = RAA10
      RAA(I) = RAANEI
 75   CONTINUE
 80   CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE FUPDAT(M,F0NEW,FNEW,F0VAL,FVAL)
C
C       Version "July 2007".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FNEW(1),FVAL(1)
C
      F0VAL=F0NEW
C
      DO 10 I=1,M
      FVAL(I)=FNEW(I)
 10   CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE APPRF(M,N,XMMA,P,Q,B,P0,Q0,R0,XLOW,XUPP,
     1                 FMAX,F0APP,FAPP)
C
C       Version "July 2007".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     APPRF calculates F0APP and FAPP.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XMMA(1),P(1),Q(1),B(1),P0(1),Q0(1),
     1          XLOW(1),XUPP(1),FMAX(1),FAPP(1)
C
      F0APP=R0
C
      DO 10 I=1,M
      FAPP(I)=FMAX(I)-B(I)
 10   CONTINUE
C
      DO 50 J=1,N
      UJXJ=XUPP(J)-XMMA(J)
      XJLJ=XMMA(J)-XLOW(J)
      F0APP=F0APP+P0(J)/UJXJ+Q0(J)/XJLJ
      MJ=M*(J-1)
      DO 40 I=1,M
      IJ=MJ+I
      FAPP(I)=FAPP(I)+P(IJ)/UJXJ+Q(IJ)/XJLJ
 40   CONTINUE
 50   CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE CONSER(M,ICONSE,GEPS,F0NEW,F0APP,FNEW,FAPP)
C
C       Version "July 2007".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     CONCER calculates ICONSE which is =1 if the current
C     approximations are conservative and =0 otherwise.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FNEW(1),FAPP(1)
C
      ICONSE=1
      F0APPE=F0APP+GEPS
      IF(F0NEW.GT.F0APPE) ICONSE=0
      DO 50 I=1,M
      FAPPE=FAPP(I)+GEPS
      IF(FNEW(I).GT.FAPPE) ICONSE=0
 50   CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE KKT(M,N,X,Y,Z,ULAM,XMIN,XMAX,DF0DX,FVAL,DFDX,FMAX,
     1               A,C,RESNOR,RESMAX)
C
C       Version "July 2007".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     KKT calculates the 2-norm RESNOR and the max-norm RESMAX
C     of the residual vector of the KKT conditions.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(1),Y(1),ULAM(1),XMIN(1),XMAX(1),DF0DX(1),
     1          FVAL(1),DFDX(1),FMAX(1),A(1),C(1)
C
      RESNOR=0.
      RESMAX=0.
C
      DO 40 J=1,N
      DLDXJ=DF0DX(J)
      MJ1=M*(J-1)
      DO 30 I=1,M
      IJ=MJ1+I
      DLDXJ=DLDXJ+ULAM(I)*DFDX(IJ)
 30   CONTINUE
      XSIJ=0.
      IF(DLDXJ.GT.0.) XSIJ=DLDXJ
      RESI=XSIJ*(X(J)-XMIN(J))
      IF(RESI.GT.RESMAX) RESMAX=RESI
      RESNOR=RESNOR+RESI**2
      ETAJ=0.
      IF(DLDXJ.LT.0.) ETAJ=-DLDXJ
      RESI=ETAJ*(XMAX(J)-X(J))
      IF(RESI.GT.RESMAX) RESMAX=RESI
      RESNOR=RESNOR+RESI**2
 40   CONTINUE
C
      DO 50 I=1,M
      FFF=FMAX(I)+A(I)*Z+Y(I)-FVAL(I)
      RESI=FFF*ULAM(I)
      IF(RESI.GT.RESMAX) RESMAX=RESI
      RESNOR=RESNOR+RESI**2
 50   CONTINUE
      RESNOR=DSQRT(RESNOR)
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
