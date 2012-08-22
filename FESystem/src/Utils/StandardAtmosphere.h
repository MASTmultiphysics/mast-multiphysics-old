//
//  StandardAtmosphere.h
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

#ifndef FESystem_StandardAtmosphere_h
#define FESystem_StandardAtmosphere_h


subroutine stdatmFORTRAN(z,t,p,r,a,mu,ts,rr,pp,rm,qm,kd,kk)
c
c   *********** 1976 STANDARD ATMOSPHERE SUBROUTINE **********
c
c     Mason's BASIC program, converted to FORTRAN - Sept. 1, 1989
c
c     kd -   = 0 - metric units
c           <> 0 - English units
c
c     kk - 0 - good return
c          1 - error: altitude out of table,
c              do not use output
c
c     z  - input altitude, in feet or meters (depending on kd)
c
c     output:
c
c     t  - temp.
c     p  - pressure
c     r  - density
c     a  - speed of sound
c     mu - viscosity
c
c     ts - t/t at sea level
c     rr - rho/rho at sea level
c     pp - p/p at sea level
c
c     rm - Reynolds number per Mach per unit of length
c     qm - dynamic pressure/Mach^2
c
implicit none
integer KK, KD
double precision k, h, t, p, r, a, mu, ml, c1, tl, pl, rl, al, bt,
$     ts , z,pp, rm, qm, rr

KK = 0
K  = 34.163195
C1 = 3.048D-04
IF (KD .eq. 0) goto 1240
TL = 518.67
PL = 2116.22
RL = .0023769
AL = 1116.45
ML = 3.7373E-07
BT = 3.0450963E-08
GOTO 1260
1240 TL = 288.15
PL = 101325
RL = 1.225
C1 = .001
AL = 340.294
ML = 1.7894E-05
BT = 1.458E-06
1260 H = C1 * Z / (1 + C1 * Z / 6356.766)
IF (H .gt. 11.0) goto 1290
T = 288.15 - 6.5 * H
PP = (288.15 / T) ** ( - K / 6.5)
GOTO 1420
1290 IF (H .gt. 20.0) goto 1310
T = 216.65
PP = .22336 *  EXP ( - K * (H - 11) / 216.65)
GOTO 1420
1310  IF (H .gt. 32.0) goto 1330
T = 216.65 + (H - 20)
PP = .054032 * (216.65 / T) ** K
GOTO 1420
1330  IF (H .gt. 47.0) goto 1350
T = 228.65 + 2.8 * (H - 32)
PP = .0085666 * (228.65 / T) ** (K / 2.8)
GOTO 1420
1350  IF( H .gt. 51.0) goto 1370
T = 270.65
PP = .0010945 *  EXP ( - K * (H - 47) / 270.65)
GOTO 1420
1370  IF (H .gt. 71.) goto 1390
T = 270.65 - 2.8 * (H - 51)
PP = .00066063 * (270.65 / T) ** ( - K / 2.8)
GOTO 1420
1390  IF (H .gt. 84.852) THEN
kk = 1
write(6,200) H
return
END IF
T = 214.65 - 2 * (H - 71)
PP = 3.9046E-05 * (214.65 / T) ** ( - K / 2)
1420  RR = PP / (T / 288.15)
MU = BT * T**1.5 / (T + 110.4)
TS = T / 288.15
A  = AL *  SQRT (TS)
T  = TL * TS
R  = RL * RR
P  = PL * PP
RM = R * A / MU
QM = .7 * P
c
200 format('   Out of Table in StdAtm- too high !'//
           1        4x,'H =',f12.3,'  > 84.852 km'/)
c
return
end

#endif
