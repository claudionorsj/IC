C
C	This is a pack of routines to be used with TRALAM and variations.
C	It includes: QTOW, WTOQ, QTOX, XTOQ, DERIVS, FITORO, ROTOFI, RK4,
C	RKQC, ODEINT, JACOBI. It is required for swing-by programs.
C
C=======================================================================
      SUBROUTINE DERIVS(U,NN,TAU,W,DWDTAU)
C=======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION W(5),DWDTAU(5)
C
      common/actejac/C0,C
C
      W1=W(1)
      W2=W(2)
      W3=W(3)
      W4=W(4)
C
      W1s=w1*w1
      W13=W1s*w1
      W14=W1s*W1s
      W2s=w2*w2
      W23=W2s*w2
      W24=W2s*W2s
      W1W2s=W1s*W2s
C
       fl2m = ((4.0d0*W13*w2 - 4.0d0*w1*W23)**2 +
     *(-1.0d0 + W14 - 6.0d0*W1W2s + W24)**2)/(4.0d0*(W1s + W2s)**3)
C
      dodw1f = ((1.0d0-u)*(8.0d0*w1*W2s + 4.0d0*w1*(1.0d0 + W1s - W2s))+
     *u*(8.0d0*w1*W2s - 4.0d0*w1*(1.0d0 - W1s + W2s)))/(W1s + W2s)**2 - 
     *4.0d0*w1*((1.0d0-u)*(4.0d0*W1W2s + (1.0d0 + W1s - W2s)**2) + 
     *u*(4.0d0*W1W2s + (1.0d0- W1s + W2s)**2))/(W1s + W2s)**3 + 
     *((4.0d0*W13*w2-4.0d0*w1*W23)**2+(-1.0d0+W14-6.0d0*W1W2s+W24)**2)*
     *(-64.0d0*c*w1*(W1s + W2s) + 
     *2.0d0*(8.0d0*w1*W2s - 4.0d0*w1*(1.0d0 - W1s + W2s))*
     *(4.0d0*W1W2s + (1.0d0 - W1s + W2s)**2) - 
     *u*(-2.0d0*(8.0d0*w1*W2s + 4.0d0*w1*(1.0d0 + W1s - W2s))*
     *(4.0d0*W1W2s + (1.0d0 + W1s - W2s)**2) + 
     *2.0d0*(8.0d0*w1*W2s - 4.0d0*w1*(1.0d0 - W1s + W2s))*
     *(4.0d0*W1W2s+(1.0d0-W1s+W2s)**2)))/(128.0d0*(W1s+W2s)**5)
     *+(2.0d0*(12.0d0*W1s*w2 - 4.0d0*W23)*(4.0d0*W13*w2 - 4.0d0*w1*W23)+
     *2.0d0*(4.0d0*W13 - 12.0d0*w1*W2s)*(-1.0d0+W14-6.0d0*W1W2s + W24))*
     *(-16.0d0*c*(W1s + W2s)**2 + 
     *(4.0d0*W1W2s + (1.0d0 - W1s + W2s)**2)**2 - 
     *u*(-((4.0d0*W1W2s + (1.0d0 + W1s - W2s)**2)**2) + 
     *(4.0d0*W1W2s + (1.0d0 - W1s + W2s)**2)**2))/
     *(128.0d0*(W1s + W2s)**5) - 5.0d0*w1*
     *((4.0d0*W13*w2-4.0d0*w1*W23)**2+(-1.0d0+W14-6.0d0*W1W2s+W24)**2)*
     *(-16.0d0*c*(W1s + W2s)**2 + 
     *(4.0d0*W1W2s + (1.0d0 - W1s + W2s)**2)**2 - 
     *u*(-((4.0d0*W1W2s + (1.0d0 + W1s - W2s)**2)**2) + 
     *(4.0d0*W1W2s + (1.0d0 - W1s + W2s)**2)**2))/
     *(64.0d0*(W1s + W2s)**6)
C
       dodw2f = ((1.0d0-u)*(8.0d0*W1s*w2-4.0d0*w2*(1.0d0+ W1s - W2s)) + 
     *u*(8.0d0*W1s*w2+4.0d0*w2*(1.0d0-W1s + W2s)))/(W1s + W2s)**2 - 
     *4.0d0*w2*((1 - u)*(4.0d0*W1W2s + (1.0d0 + W1s - W2s)**2) + 
     *u*(4.0d0*W1W2s + (1.0d0-W1s + W2s)**2))/(W1s + W2s)**3 + 
     *((4.0d0*W13*w2-4.0d0*w1*W23)**2+(-1.0d0+W14-6.0d0*W1W2s+W24)**2)*
     *(-64.0d0*c*w2*(W1s + W2s) + 
     *2.0d0*(8.0d0*W1s*w2 + 4.0d0*w2*(1.0d0 - W1s + W2s))*
     *(4.0d0*W1W2s + (1.0d0 - W1s + W2s)**2) - 
     *u*(-2.0d0*(8.0d0*W1s*w2 - 4.0d0*w2*(1.0d0 + W1s - W2s))*
     *(4.0d0*W1W2s + (1.0d0 + W1s - W2s)**2) + 
     *2.0d0*(8.0d0*W1s*w2 + 4.0d0*w2*(1.0d0 - W1s + W2s))*
     *(4.0d0*W1W2s+(1.0d0-W1s+W2s)**2)))/(128.0d0*(W1s+W2s)**5)
     *+ (2.0d0*(4.0d0*W13 - 12.0d0*w1*W2s)*(4.0d0*W13*w2-4.0d0*w1*W23)+
     *2.0d0*(-12.0d0*W1s*w2+4.0d0*W23)*(-1.0d0+W14-6.0d0*W1W2s+W24))*
     *(-16.0d0*c*(W1s + W2s)**2 + 
     *(4.0d0*W1W2s + (1.0d0-W1s + W2s)**2)**2 - 
     *u*(-((4.0d0*W1W2s + (1.0d0 + W1s - W2s)**2)**2) + 
     *(4.0d0*W1W2s + (1.0d0 - W1s + W2s)**2)**2))/
     *(128.0d0*(W1s + W2s)**5) - 5.0d0*w2*
     *((4.0d0*W13*w2-4.0d0*w1*W23)**2+(-1.0d0+W14-6.0d0*W1W2s+W24)**2)*
     *(-16.0d0*c*(W1s + W2s)**2 + 
     *(4.0d0*W1W2s + (1.0d0 - W1s + W2s)**2)**2 - 
     *u*(-((4.0d0*W1W2s + (1.0d0 + W1s - W2s)**2)**2) + 
     *(4.0d0*W1W2s + (1.0d0 - W1s + W2s)**2)**2))/
     *(64.0d0*(W1s + W2s)**6)
C
      DWDTAU(1)=W3
      DWDTAU(2)=W4
      DWDTAU(3)=2.0D0*W4*FL2M+DODW1F
      DWDTAU(4)=-2.0D0*W3*FL2M+DODW2F
      DWDTAU(5)=FL2M
C
      RETURN
      END
C
C=======================================================================
      SUBROUTINE QTOW(Q,W)
C=======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 WA,QA,FL,WAV,QAV
      DIMENSION W(5),Q(5)
      QA=DCMPLX(Q(1),Q(2))
      QAV=DCMPLX(Q(3),Q(4))
      WA=CDSQRT(2.0D0*QA+CDSQRT(4.0D0*QA*QA-1.0D0))
      W(1)=DREAL(WA)
      W(2)=DIMAG(WA)
      FL=(WA-1.0D0/(WA*WA*WA))/2.0D0
      WAV=(CDABS(FL)*CDABS(FL))*QAV/FL
      W(3)=DREAL(WAV)
      W(4)=DIMAG(WAV)
      RETURN
      END
C
C=======================================================================
      SUBROUTINE QTOX(Q,U,X)
C=======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Q(5),X(5)
      X(1)=Q(1)-0.5D0+U
      X(2)=Q(2)
      X(3)=Q(3)
      X(4)=Q(4)
      RETURN
      END
C
C=======================================================================
      SUBROUTINE ROTOFI(XR,T,XF)
C=======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XF(5),XR(5)
      CT=DCOS(T)
      ST=DSIN(T)
      XF(1)=XR(1)*CT-XR(2)*ST
      XF(2)=XR(2)*CT+XR(1)*ST
      XF(3)=-XR(1)*ST-XR(4)*ST+XR(3)*CT-XR(2)*CT
      XF(4)=XR(3)*ST-XR(2)*ST+XR(1)*CT+XR(4)*CT
      RETURN
      END
C
C=======================================================================
      SUBROUTINE WTOQ(W,Q)
C=======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 WA,QA,FL,WAV,QAV
      DIMENSION W(5),Q(5)
      WA=DCMPLX(W(1),W(2))
      WAV=DCMPLX(W(3),W(4))
      QA=(WA*WA+1.0D0/(WA*WA))/4.0D0
      Q(1)=DREAL(QA)
      Q(2)=DIMAG(QA)
      FL=(WA-1.0D0/(WA*WA*WA))/2.0D0
      QAV=FL*WAV/(CDABS(FL)*CDABS(FL))
      Q(3)=DREAL(QAV)
      Q(4)=DIMAG(QAV)
      RETURN
      END
C
C=======================================================================
      SUBROUTINE XTOQ(X,U,Q)
C=======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Q(5),X(5)
      Q(1)=X(1)+0.5D0-U
      Q(2)=X(2)
      Q(3)=X(3)
      Q(4)=X(4)
      RETURN
      END
C
C=======================================================================
      SUBROUTINE JACOBI(W,U,C1)
C=======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION W(5)
      COMPLEX*16 WA
      WA=DCMPLX(W(1),W(2))
      UL=W(3)
      VL=W(4)
      OME=(UL*UL+VL*VL)/2.0D0
      WM=CDABS(WA)
      WP1=CDABS(WA*WA+1.0D0)
      WM1=CDABS(1.0D0-WA*WA)
      AU=(U*WM1*WM1+WP1*WP1*(1.0D0-U))/WM**4
      AM=(CDABS(1.0D0-WA**4))**2/128.0D0/WM**10
      AU2=WM1**4-U*(WM1**4-WP1**4)
      C1=((OME-AU)/AM-AU2)/(-16.0D0*WM**4)
      RETURN
      END
C
C=======================================================================
	   SUBROUTINE JACOBIW(W,U,C)
C=======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XR(5),W(6),Q(5)
      CALL WTOQ(W,Q)
      CALL QTOX(Q,U,XR)
      V=XR(3)*XR(3)+XR(4)*XR(4)
      R12=(XR(1)-U)**2+XR(2)*XR(2)
      R22=(XR(1)+1.0D0-U)**2+XR(2)*XR(2)
      R1=DSQRT(R12)
      R2=DSQRT(R22)
      OME=(1.0D0-U)/R1+U/R2+((1.0D0-U)*R1*R1+U*R2*R2)/2.0D0
      C=2.0D0*OME-V
      RETURN
      END
C
C=======================================================================
      SUBROUTINE FITORO(XF,T,XR)
C=======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XF(5),XR(5)
      XR(1)=XF(1)*DCOS(T)+XF(2)*DSIN(T)
      XR(2)=XF(2)*DCOS(T)-XF(1)*DSIN(T)
      XR(4)=-XF(1)*DCOS(T)+XF(4)*DCOS(T)-XF(3)*DSIN(T)-XF(2)*DSIN(T)
      XR(3)=XF(3)*DCOS(T)+XF(2)*DCOS(T)-XF(1)*DSIN(T)+XF(4)*DSIN(T)
      RETURN
      END
