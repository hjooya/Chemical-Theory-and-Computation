Program TD_Equal_Spacing

IMPLICIT NONE
integer,parameter                   :: totstep=400000
integer,parameter                   :: MM=21  ! Must be odd
integer,parameter                   :: DD=2
integer,parameter                   :: NN=2*MM
integer,parameter                   :: PP=MM*MM
integer,parameter                   :: QQ=NN*MM
integer,parameter                   :: LDA=NN
double precision,dimension(NN,NN)     :: H
integer,parameter                   :: LWORK=64*NN
double precision,dimension(LWORK)   :: WORK
integer,parameter                   :: LWORK2=64*QQ
double precision,dimension(LWORK2)   :: WORK2
double precision,parameter          :: pi=3.1415926D0 
!double precision,parameter          :: eof=0.0856454d0  !563.5nm   
integer i,j,k,nt,iter,jjj,kkk,iii
double precision,dimension(NN)       :: EIG
double precision,dimension(NN,NN)   :: A
double precision,dimension(NN)   :: Arow

double precision,dimension(QQ)       :: EIG2
double precision,dimension(QQ,QQ)   :: A2
double precision,dimension(QQ)   :: Arow2

double precision,dimension(NN)   :: wf0
double precision,dimension(NN)   :: qwf
double precision,dimension(NN)   :: QE
double precision,dimension(NN,NN)   :: S1
double precision,dimension(NN)   :: p_d
double precision,dimension(NN,NN)   :: S2
double precision,dimension(NN,NN)   :: INITIAL
double precision,dimension(NN,NN)   :: FINAL
double precision,dimension(NN*MM,NN*MM)   :: AA
integer ifail,INFO,ipp
CHARACTER(1)                        :: UPLO
double precision,dimension(2,2)    :: TWO_LEVEL1
double precision,dimension(2,2)    :: TWO_LEVEL2
double precision,dimension(MM,MM)    :: IDENT_TWO
double precision,dimension(NN,NN)  :: TWO_IDENT1
double precision,dimension(NN,NN)  :: TWO_QUBITS
double precision,dimension(MM,MM)          :: IDENT_COUP
double precision,dimension(NN*MM,NN*MM)    :: COUP_IDENT

double precision,dimension(MM,MM)    :: BOSON1
double precision,dimension(2,2)    :: IDENT_BOS
double precision,dimension(NN,NN)  :: BOSON1_IDENT
double precision,dimension(NN,NN)  :: GRAPH1

double precision,dimension(MM,MM)    :: BOSON2
double precision,dimension(NN,NN)    :: IDENT_BOS2
double precision,dimension(QQ,QQ)  :: BOSON2_IDENT
double precision,dimension(QQ,QQ)  :: BOSON1T_IDENT
double precision,dimension(QQ,QQ)  :: BOSONTotalT
double precision,dimension(QQ,QQ)  :: GRAPH2
double precision,dimension(QQ,QQ)  :: TWO_IDENT2
integer d,n1,m1,P1,Q1,ii,jj,ip,N_S,half1,half2,j1,j2
integer n2,m2,P2,Q2,n,ppp,i_QE
double complex,dimension(totstep)         :: sprob1
double complex,dimension(3*totstep+100)            :: COMM
double precision w12,eof,prob,temp,norm,prob2,prob0,dt,prob1,prob4
double precision nn1,mm1,ll1,kk1,ss1,coupling,RF1,RF2,delta,intensity
double precision  noc1,Time,prob23,prob24,norm_max,normal,Sum,W1,W2,eof1,eof2,RF2set,ELIP,time1
integer,parameter                   :: ipptotal=10000
double precision,dimension(ipptotal)   :: Eext
double precision tt,p_time,Phase
complex*16                          :: a_t
complex*16                          :: b_t
double precision EIGChos1,EIGChos2


noc1=16
dt = (noc1*2.0D0*pi/eof)/dble(totstep) 
Time= 2.d0*pi/eof

intensity=1.0d0*eof

RF1= (0.25d0/4.0d0)/80.d0
RF2 = RF1



W12=2.0d0


DO iii=0,200

eof1=1.97d0+dble(iii)*0.0003d0


DO jjj=0,200

eof2=1.97d0+dble(jjj)*0.0003d0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! FIRST TWO-LEVEL SYSTEM Laplacian

DO i=1,2
DO j=1,2
 IF(i.NE.j) THEN
  TWO_LEVEL1(i,j)=0.0d0
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


!! FIRST BOSON LAPLACIAN

DO i=1, MM
DO j=1, MM
      BOSON1(i,j)=0.0d0
END DO
END DO



DO i=1, MM
DO j=1, MM
 IF(i.EQ.j) THEN
   BOSON1(i,j)=(   dble(i)-dble(MM)+((dble(MM)-1.0d0)/2.0d0) ) *eof1
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
  BOSON1_IDENT((i-1)*P2+k,(j-1)*Q2+d)=IDENT_BOS(i,j)*BOSON1(k,d)
END DO
END DO
END DO
END DO


DO i=1, NN
DO j=1, NN
  IF(i.LE.MM) THEN
        BOSON1_IDENT(i,i+MM+1)=RF1
        BOSON1_IDENT(i,i+MM-1)=RF1
  END IF
END DO
END DO

BOSON1_IDENT(1,MM)=0.0d0


DO i=1, NN
DO j=1, NN
  IF(i.GT.MM) THEN
        BOSON1_IDENT(i,i-MM+1)=RF1
        BOSON1_IDENT(i,i-MM-1)=RF1
  END IF
END DO
END DO

BOSON1_IDENT(NN,MM+1)=0.0d0



!!!!! GRAPH ADDITION

DO i=1, NN
DO j=1, NN
   GRAPH1(i,j)=BOSON1_IDENT(i,j)+TWO_IDENT1(i,j)
END DO
END DO


!!! SOLVING THE SYSTEM

UPLO='L'
   ifail=0

     DO i=1,NN
           DO j=1,NN
                  a(i,j)=GRAPH1(i,j)
           END DO
     END DO


CALL DSYEV('V','U',NN,A,NN,EIG,WORK,LWORK,INFO)
!WRITE(6,*) 'INFO',INFO


!DO i=1,NN
!    WRITE(120,*) w12/eof, EIG(i)/eof
!END DO

!DO i=2,NN,2
!  IF( (EIG(i).GE.-1).AND.(EIG(I).LE.5)  ) THEN
!    WRITE(121,*) w12/eof, EIG(i)/eof
!  END IF
!END DO

!DO i=1,NN,2
!  IF( (EIG(i).GE.-1).AND.(EIG(I).LE.5)  ) THEN
!    WRITE(122,*) w12/eof, EIG(i)/eof
!  END IF
!END DO

!!!!!!!!!!!!!!!!!!!!!!!! TRANSITION PROBABILITY !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO i=1, NN
DO j=1,NN
  IF(j.EQ.MM-2) THEN
    Arow(i)=A(i,j)
    EIGChos1=EIG(j)
!    write(101,*) i, Arow(i)
  END IF
END DO
END DO


IF(iii.EQ.1) THEN
DO i=1,MM
  DO j=MM+1,NN
   IF(abs(i+MM-j).LE.(MM-1)) THEN
!   write(157,*) i,j, abs(i+MM-j)
!   write(158,*) ABS(BOSON_IDENT(i,i)-BOSON_IDENT(j,j))

   END IF
  END DO 
END DO
END IF



prob1=0.0d0
DO i=1,MM
  DO j=MM+1,NN
   IF(abs(i+MM-j).LE.(MM-1)) THEN
    prob1= prob1+(Arow(i)*Arow(j))**2
   END IF
  END DO 
END DO


!      write(77,*) w12/eof1 , 2.0d0*prob1



!!!!!!!!!!!!!!!!!!!!!!!SECOND COLOR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! KRONECKER PRODUCT

N1=NN
M1=NN
P1=MM
Q1=MM


DO i=1,N1
DO j=1,M1
DO k=1,P1
DO d=1,Q1
  TWO_IDENT2((i-1)*P1+k,(j-1)*Q1+d)=GRAPH1(i,j)*IDENT_TWO(k,d)
END DO
END DO
END DO
END DO



!! SECOND BOSON LAPLACIAN

DO i=1, MM
DO j=1, MM
      BOSON2(i,j)=0.0d0
END DO
END DO



DO i=1, MM
DO j=1, MM
 IF(i.EQ.j) THEN
   BOSON2(i,j)=(   dble(i)-dble(MM)+((dble(MM)-1.0d0)/2.0d0) ) *eof2
 END IF
END DO 
END DO


!! IDENTITY MATRIX

DO i=1,NN
DO j=1,NN
 IF(i.EQ.j) THEN
  IDENT_BOS2(i,j)=1.0d0
 ELSE
  IDENT_BOS2(i,j)=0.0d0
 END IF
END DO
END DO



!! KRONECKER PRODUCT

N2=NN
M2=NN
P2=MM
Q2=MM


DO i=1,N2
DO j=1,M2
DO k=1,P2
DO d=1,Q2
  BOSON2_IDENT((i-1)*P2+k,(j-1)*Q2+d)=IDENT_BOS2(i,j)*BOSON2(k,d)
END DO
END DO
END DO
END DO


DO i=1, QQ
DO j=1, QQ
  IF(i.LE.PP) THEN
        BOSON2_IDENT(i,i+PP+1)=RF2
        BOSON2_IDENT(i,i+PP-1)=RF2
  END IF
END DO
END DO

BOSON2_IDENT(1,PP)=0.0d0


DO i=1, QQ
DO j=1, QQ
  IF(i.GT.PP) THEN
        BOSON2_IDENT(i,i-PP+1)=RF2
        BOSON2_IDENT(i,i-PP-1)=RF2 
  END IF
END DO
END DO

BOSON2_IDENT(QQ,PP+1)=0.0d0





!!!!!!!!!!!!!!!!! WE NEED BOSON1 for Time-Dependnet !!!!!!!!!!!!!!!!!!!!!!!

!! KRONECKER PRODUCT

N2=NN
M2=NN
P2=MM
Q2=MM


DO i=1,N2
DO j=1,M2
DO k=1,P2
DO d=1,Q2
  BOSON1T_IDENT((i-1)*P2+k,(j-1)*Q2+d)=IDENT_BOS2(i,j)*BOSON1(k,d)
END DO
END DO
END DO
END DO


DO i=1, QQ
DO j=1, QQ
  IF(i.LE.PP) THEN
        BOSON1T_IDENT(i,i+PP+1)=0.0d0
        BOSON1T_IDENT(i,i+PP-1)=0.0d0
  END IF
END DO
END DO

BOSON1T_IDENT(1,PP)=0.0d0


DO i=1, QQ
DO j=1, QQ
  IF(i.GT.PP) THEN
        BOSON1T_IDENT(i,i-PP+1)=0.0d0
        BOSON1T_IDENT(i,i-PP-1)=0.0d0
  END IF
END DO
END DO

BOSON1T_IDENT(QQ,PP+1)=0.0d0



DO i=1, QQ
DO j=1, QQ
  IF(i.EQ.j) THEN
   BOSONTotalT(i,j)=BOSON1T_IDENT(i,j) + BOSON2_IDENT(i,j)
  END IF
END DO
END DO



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! GRAPH ADDITION

DO i=1, QQ
DO j=1, QQ
   GRAPH2(i,j)=BOSON2_IDENT(i,j)+TWO_IDENT2(i,j)
END DO
END DO


!!!! SOLVING THE SYSTEM

UPLO='L'
   ifail=0

     DO i=1,QQ
           DO j=1,QQ
                  A2(i,j)=GRAPH2(i,j)
           END DO
     END DO


CALL DSYEV('V','U',QQ,A2,QQ,EIG2,WORK2,LWORK2,INFO)
!WRITE(6,*) 'INFO',INFO


!DO i=1,NN
!    WRITE(120,*) w12/eof, EIG(i)/eof
!END DO

!DO i=2,NN,2
!  IF( (EIG(i).GE.-1).AND.(EIG(I).LE.5)  ) THEN
!    WRITE(121,*) w12/eof, EIG(i)/eof
!  END IF
!END DO

!DO i=1,NN,2
!  IF( (EIG(i).GE.-1).AND.(EIG(I).LE.5)  ) THEN
!    WRITE(122,*) w12/eof, EIG(i)/eof
!  END IF
!END DO

!!!!!!!!!!!!!!!!!!!!!!!! TRANSITION PROBABILITY !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO i=1, QQ
DO j=1,QQ
  IF(j.EQ.PP-2) THEN
    Arow2(i)=A2(i,j)
    EIGChos2=EIG2(j)
!    write(101,*) i, Arow2(i)
  END IF
END DO
END DO



prob1=0.0d0
DO i=1,PP
  DO j=PP+1,QQ
!   IF(abs(i+PP-j).LE.(PP-1)) THEN
    prob1= prob1+(Arow2(i)*Arow2(j))**2
!   END IF
  END DO 
END DO


      write(10,*) eof1, eof2, 2.0d0*prob1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END DO   ! iii


write(10,*)


END DO   ! jjj


END   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

