Program Longitudinal_Coupling

IMPLICIT NONE
integer,parameter                   :: totstep=40000
integer,parameter                   :: MM=51  ! Must be odd
integer,parameter                   :: DD=2
integer,parameter                   :: NN=2*MM
integer,parameter                   :: LDA=NN
double precision,dimension(NN,NN)   :: H
integer,parameter                   :: LWORK=64*NN
double precision,dimension(LWORK)   :: WORK
integer,parameter                   :: LWORK2=64*NN*MM
double precision,dimension(LWORK2)  :: WORK2
double precision,dimension(NN)      :: EIG
double precision,dimension(NN,NN)   :: A
double precision,dimension(NN)      :: Arow
double precision,dimension(2,2)     :: TWO_LEVEL1
double precision,dimension(MM,MM)   :: IDENT_TWO
double precision,dimension(NN,NN)   :: TWO_IDENT1
double precision,dimension(MM,MM)   :: BOSON
double precision,dimension(2,2)     :: IDENT_BOS
double precision,dimension(NN,NN)   :: BOSON_IDENT
double precision,dimension(NN,NN)   :: GRAPH
CHARACTER(1)                        :: UPLO
integer ifail,INFO
integer i,j,k,jjj,kkk,iii
integer d,n1,m1,P1,Q1,ii,jj,ip
integer n2,m2,P2,Q2,n,ppp
double complex,dimension(3*totstep+100)            :: COMM
double precision w12,eof,prob1,coupling,RF,delta,intensity


eof=1.0d0

intensity=5.0d0*eof
RF= 0.5d0*Intensity

coupling=0.1d0*eof



DO iii=1,totstep

w12 =dble(iii)*0.00025d0 

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
        BOSON(i,i+1)=RF
        BOSON(i,i-1)=RF
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


!!!!! GRAPH ADDITION

DO i=1, NN
DO j=1, NN
   GRAPH(i,j)=BOSON_IDENT(i,j)+TWO_IDENT1(i,j)
END DO
END DO

DO i=1, NN
DO j=1, NN
  IF(i.EQ.j) THEN
   IF(i.GT.MM) THEN
     GRAPH(i,i+1)=-GRAPH(i,i+1)
     GRAPH(i,i-1)=-GRAPH(i,i-1)
   END IF
  END IF
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


DO i=1,NN
IF ( (EIG(i).GE.-2.0d0).AND.(EIG(i).LE.2.0d0) ) THEN
    WRITE(120,*) w12/eof, EIG(i)/eof
END IF
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
IF (  (j.NE.(i+MM-3))  .AND. (j.NE.(i+MM+3))    ) THEN  ! This removes three photon transfers

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