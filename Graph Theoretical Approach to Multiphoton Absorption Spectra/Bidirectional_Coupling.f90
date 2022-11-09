Program Bidirectional_Coupling

IMPLICIT NONE
integer,parameter                   :: totstep=20000
integer,parameter                   :: MM=31 !Must be odd
integer,parameter                   :: DD=2
integer,parameter                   :: NN=2*MM
integer,parameter                   :: LDA=NN
integer,parameter                   :: LWORK=64*NN
double precision,dimension(LWORK)   :: WORK
integer,parameter                   :: LWORK2=64*NN*MM
double precision,dimension(LWORK2)  :: WORK2
double precision,parameter          :: pi=3.1415926D0 
!double precision,parameter         :: eof=0.0856454d0  !563.5nm   
double precision,dimension(NN)      :: EIG
double precision,dimension(NN,NN)   :: A
double precision,dimension(NN)      :: Arow
CHARACTER(1)                        :: UPLO
double precision,dimension(2,2)     :: TWO_LEVEL1
double precision,dimension(MM,MM)   :: IDENT_TWO
double precision,dimension(NN,NN)   :: TWO_IDENT1
double precision,dimension(MM,MM)   :: BOSON
double precision,dimension(2,2)     :: IDENT_BOS
double precision,dimension(NN,NN)   :: BOSON_IDENT
double precision,dimension(NN,NN)   :: BOSON_IDENTD
double precision,dimension(NN,NN)   :: BOSON_IDENTO
double precision,dimension(NN,NN)   :: GRAPH
double complex,dimension(3*totstep+100)   :: COMM
integer ifail,INFO
integer i,j,k,nt,iter,jjj,kkk,iii
integer d,n1,m1,P1,Q1,ii,jj,ip,fff
integer n2,m2,P2,Q2,n,ppp,i_QE
double precision w12,eof,prob,temp,norm,prob2,prob1,dt,Elip
double precision nn1,mm1,ll1,kk1,ss1,coupling,RF,delta,intensity
double precision  noc1,Time,prob23,prob24,norm_max,prob0,normal,RFD,RFO,Intensity1,intensity2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

eof=1.0d0 !0.0856454d0

noc1=16
dt = (noc1*2.0D0*pi/eof)/dble(totstep) 
Time= 2.d0*pi/eof

intensity1=5.0d0*eof
!RFD= 0.25d0*Intensity1

intensity2=1.0d0*eof
!RFO= 0.25d0*Intensity2

delta=0.5d0*eof
!coupling=0.5d0*delta


DO jjj=0,200


IF(jjj.LE.100) THEN
  Elip = dble(jjj)*0.01d0
  RFO=0.25d0*eof
  coupling = Elip*(0.5d0*delta)
  RFD=Elip*(0.25d0*Intensity1)
ELSE IF(jjj.GT.100) THEN
  fff=jjj-100 
  Elip = dble(fff)*0.01d0
  RFO=(1.0d0-Elip)*0.25d0*eof
  coupling = 0.5d0*delta
  RFD=0.25d0*Intensity1
END IF


DO iii=1,totstep

w12 =dble(iii)*0.0005d0  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! FIRST TWO-LEVEL SYSTEM Laplacian

DO i=1,2
DO j=1,2
 IF(i.NE.j) THEN
  TWO_LEVEL1(i,j)=-coupling
 ELSE IF((i.EQ.1).AND.(j.EQ.1)) THEN
  TWO_LEVEL1(i,j)=0.5d0*W12
 ELSE IF((i.EQ.2).AND.(j.EQ.2)) THEN
  TWO_LEVEL1(i,j)=-0.5d0*W12
 END IF
END DO
END DO


!! IDENTITY MATRIX

DO i=1,MM
DO j=1,MM
 IF(i.EQ.j) THEN
  IDENT_TWO(i,j)=1.0d0
 ELSE
  IDENT_TWO(i,j)=0.0d0
 END IF
END DO
END DO


!! KRONECKER PRODUCT

N1=2
M1=2
P1=MM
Q1=MM


DO i=1,N1
DO j=1,M1
DO k=1,P1
DO d=1,Q1
  TWO_IDENT1((i-1)*P1+k,(j-1)*Q1+d)=TWO_LEVEL1(i,j)*IDENT_TWO(k,d)
END DO
END DO
END DO
END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! BOSON LAPLACIAN

DO i=1, MM
DO j=1, MM
      BOSON(i,j)=0.0d0
END DO
END DO


DO i=1, MM
DO j=1, MM
 IF(i.EQ.j) THEN
   BOSON(i,j)=(   dble(i)-dble(MM)+((dble(MM)-1.0d0)/2.0d0) ) *eof
        BOSON(i,i+1)=RFD
        BOSON(i,i-1)=RFD
 END IF
END DO 
END DO


!! IDENTITY MATRIX

DO i=1,2
DO j=1,2
 IF(i.EQ.j) THEN
  IDENT_BOS(i,j)=1.0d0
 ELSE
  IDENT_BOS(i,j)=0.0d0
 END IF
END DO
END DO



!! KRONECKER PRODUCT

N2=2
M2=2
P2=MM
Q2=MM


DO i=1,N2
DO j=1,M2
DO k=1,P2
DO d=1,Q2
  BOSON_IDENTD((i-1)*P2+k,(j-1)*Q2+d)=IDENT_BOS(i,j)*BOSON(k,d)
END DO
END DO
END DO
END DO


DO i=1, NN
DO j=1, NN
  IF(i.EQ.j) THEN
   IF(i.GT.MM) THEN
     BOSON_IDENTD(i,i+1)=-BOSON_IDENTD(i,i+1)
     BOSON_IDENTD(i,i-1)=-BOSON_IDENTD(i,i-1)
   END IF
  END IF
END DO
END DO


!!!!!!!!!!!!!!!! TRANSVERSE COUPLING !!!!!!!!!!!!!!!!


DO i=1, NN
DO j=1, NN
     BOSON_IDENTO(i,j)=BOSON_IDENTD(i,j)
     BOSON_IDENTO(i,j)=BOSON_IDENTD(i,j)
END DO
END DO

DO i=1, NN
DO j=1, NN
     BOSON_IDENTO(i,j)=0.0d0
     BOSON_IDENTO(i,j)=0.0d0
END DO
END DO


DO i=1, NN
DO j=1, NN
  IF(i.LE.MM) THEN
        BOSON_IDENTO(i,i+MM+1)=RFO
        BOSON_IDENTO(i,i+MM-1)=RFO
  END IF
END DO
END DO

BOSON_IDENTO(1,MM)=0.0d0


DO i=1, NN
DO j=1, NN
  IF(i.GT.MM) THEN
        BOSON_IDENTO(i,i-MM+1)=RF
        BOSON_IDENTO(i,i-MM-1)=RF 
  END IF
END DO
END DO

BOSON_IDENTO(NN,MM+1)=0.0d0



DO i=1, NN
DO j=1, NN
        BOSON_IDENT(i,j)=BOSON_IDENTD(i,j)+BOSON_IDENTO(i,j)
END DO
END DO


!!!!! GRAPH ADDITION

DO i=1, NN
DO j=1, NN
   GRAPH(i,j)=BOSON_IDENT(i,j)+TWO_IDENT1(i,j)
END DO
END DO



!!!! SOLVING THE SYSTEM

UPLO='L'
   ifail=0

     DO i=1,NN
           DO j=1,NN
                  a(i,j)=GRAPH(i,j)
           END DO
     END DO

CALL DSYEV('V','U',NN,A,NN,EIG,WORK,LWORK,INFO)
!WRITE(6,*) 'INFO',INFO


IF(jjj.EQ.0) THEN
  DO i=1,NN
    WRITE(120,*) w12/eof, EIG(i)/eof
  END DO
ELSE IF(jjj.EQ.50) THEN
  DO i=1,NN
    WRITE(121,*) w12/eof, EIG(i)/eof
  END DO
ELSE IF(jjj.EQ.100) THEN
  DO i=1,NN
    WRITE(122,*) w12/eof, EIG(i)/eof
  END DO
ELSE IF(jjj.EQ.150) THEN
  DO i=1,NN
    WRITE(123,*) w12/eof, EIG(i)/eof
  END DO
ELSE IF(jjj.EQ.200) THEN
  DO i=1,NN
    WRITE(124,*) w12/eof, EIG(i)/eof
  END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!! TRANSITION PROBABILITY !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


DO i=1, NN
DO j=1,NN
  IF(j.EQ.MM-2) THEN
    Arow(i)=A(i,j)
  END IF
END DO
END DO



prob1=0.0d0
DO i=1,MM
  DO j=MM+1,NN
   IF(abs(i+MM-j).LE.(MM-1)) THEN
    prob1= prob1+(Arow(i)*Arow(j))**2
   END IF
  END DO 
END DO


IF(jjj.LE.100) THEN
      write(780,*) Elip, w12/eof , 2.0d0*prob1
ELSE IF(jjj.GT.100) THEN
      write(790,*) Elip, w12/eof , 2.0d0*prob1
END IF


IF(jjj.EQ.0) THEN
      write(78,*) w12/eof , 2.0d0*prob1
ELSE IF(jjj.EQ.50) THEN
      write(79,*) w12/eof , 2.0d0*prob1
ELSE IF(jjj.EQ.100) THEN
      write(80,*) w12/eof , 2.0d0*prob1
ELSE IF(jjj.EQ.150) THEN
      write(81,*) w12/eof , 2.0d0*prob1
ELSE IF(jjj.EQ.200) THEN
      write(82,*) w12/eof , 2.0d0*prob1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



END DO   ! for iii loop


write(780,*)
write(790,*)



END DO   ! for jjj loop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END