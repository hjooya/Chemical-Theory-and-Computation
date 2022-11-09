Program Transverse_Coupling

IMPLICIT NONE
integer,parameter                  :: totstep=40000
integer,parameter                  :: MM=31  ! Must be odd
integer,parameter                  :: DD=2
integer,parameter                  :: NN=2*MM
integer,parameter                  :: LDA=NN
double precision,dimension(NN,NN)  :: H
integer,parameter                  :: LWORK=64*NN
double precision,dimension(LWORK)  :: WORK
integer,parameter                  :: LWORK2=64*NN*MM
double precision,dimension(LWORK2) :: WORK2
double precision,dimension(NN)     :: EIG
double precision,dimension(NN,NN)  :: A
double precision,dimension(NN)     :: Arow
CHARACTER(1)                       :: UPLO
double precision,dimension(2,2)    :: TWO_LEVEL1
double precision,dimension(MM,MM)  :: IDENT_TWO
double precision,dimension(NN,NN)  :: TWO_IDENT1
double precision,dimension(MM,MM)  :: BOSON
double precision,dimension(2,2)    :: IDENT_BOS
double precision,dimension(NN,NN)  :: BOSON_IDENT
double precision,dimension(NN,NN)  :: GRAPH
double complex,dimension(3*totstep+100) :: COMM
integer ifail,INFO
integer d,n1,m1,P1,Q1,ii,jj
integer i,j,k,jjj,kkk,iii
integer n2,m2,P2,Q2,n
double precision w12,eof,prob1,RF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


eof=1.0d0 !0.0856454d0
RF= 0.25d0*eof


DO iii=1,totstep

w12 =dble(iii)*0.00025d0

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
  BOSON_IDENT((i-1)*P2+k,(j-1)*Q2+d)=IDENT_BOS(i,j)*BOSON(k,d)
END DO
END DO
END DO
END DO


DO i=1, NN
DO j=1, NN
  IF(i.LE.MM) THEN
        BOSON_IDENT(i,i+MM+1)=RF
        BOSON_IDENT(i,i+MM-1)=RF
  END IF
END DO
END DO

BOSON_IDENT(1,MM)=0.0d0


DO i=1, NN
DO j=1, NN
  IF(i.GT.MM) THEN
        BOSON_IDENT(i,i-MM+1)=RF
        BOSON_IDENT(i,i-MM-1)=RF 
  END IF
END DO
END DO

BOSON_IDENT(NN,MM+1)=0.0d0



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


!!!!! QuasiEnergies

DO i=1,NN
    WRITE(120,*) w12/eof, EIG(i)/eof
END DO

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


  write(77,*) w12/eof , 2.0d0*prob1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END DO


END   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
