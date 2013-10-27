      program p_zvc
c-------------------------
c     Calculo das curvas de velocidade zero dado valor de C ou um ponto 
c     do espaco de configuracao (CPR3BP), usando as posicoes de primarios 
c     da Caltech
c-------------------------
      use m_sysconst
      use m_hill
   
      implicit real*8 (A-H,O-Z)
      dimension y(4)
      
      write (*,*)'Entre com C:'
      read (*,*) C
c-----------------------------------------------------------------------
      call zvc(C)
      STOP
      END
c      include 's_zvc_bacia.f'
c=======================================================================
      function CteJac (x1,x2,x3,x4)
c=======================================================================
      use m_sysconst
      implicit real*8 (A-H,O-Z)
      
      R1=dsqrt((x1-xP1)**2.0d0+x2**2.0d0)
      R2=dsqrt((x1-xP2)**2.0d0+x2**2.0d0)
      
      Omega=0.5d0*(x1*x1+x2*x2)+((1.0d0-u)/R1)+
     &     (u/R2)+0.5d0*u*(1.0d0-u)
      CteJac=2.0d0*Omega-(x3*x3+x4*x4)
      
      return
      end
c=======================================================================
      subroutine zvc (CCC)
c-----------------------------------------------------------------------
c     C = Const de Jac fixa e xponto = yponto = 0 e dado x, quem é o y
c     correspondente. Queremos resolver a equacao C-2*OMEGA=0
c-----------------------------------------------------------------------
      use m_sysconst
      use m_hill
      implicit real*8 (A-H,O-Z)
      external rtsafe,fundc
      dimension p(100)

      open(300, status='unknown', file ='Dados/c_zvc_neg1.dat')
      open(302, status='unknown', file ='Dados/c_zvc_neg2.dat')
      open(303, status='unknown', file ='Dados/c_zvc_pos1.dat')
      open(304, status='unknown', file ='Dados/c_zvc_pos2.dat')

      open(200, status='unknown', file ='DadosBarcelona/c_zvc_neg1.dat')
      open(202, status='unknown', file ='DadosBarcelona/c_zvc_neg2.dat')
      open(203, status='unknown', file ='DadosBarcelona/c_zvc_pos1.dat')
      open(204, status='unknown', file ='DadosBarcelona/c_zvc_pos2.dat')

      open(3, status='unknown',file = 'Dados/c_lagrangianos.dat')
      write(3,*)xL1,yL1
      write(3,*)xL2,yL2
      write(3,*)xL3,yL3
      write(3,*)xL4,yL4
      write(3,*)xL5,yL5
      close(3)

      open(3,status='unknown',file ='DadosBarcelona/c_lagrangianos.dat')
      write(3,*)-xL1,yL1
      write(3,*)-xL2,yL2
      write(3,*)-xL3,yL3
      write(3,*)-xL4,yL4
      write(3,*)-xL5,yL5
      close(3)
      tol=1.0d-8
c-----------------------------------------------------------------------
      v1=0.0d0
      v2=0.0d0
      b1=-2.0d0
      b2=2.0d0
      dbx=0.0001d0
      dby=0.01d0
      zero=0.0d0
      
      x=b1
      write(*,*)'Construindo zvc'

      do while (x.le.b2) 
c         write(*,*)'x=',x
         x1=x  
         y=b1
         kcont=0
         Ncont=0
         C2=CteJac(x1,y,v1,v2)
         do while (y.le.zero)
c            write(*,*)y
            C1=C2
            ynovo=y+dby
            C2=CteJac(x1,ynovo,v1,v2)
            prod=(CCC-C1)*(CCC-C2)
            if (prod.le.0.0d0)then
               x2=rtsafe(fundc,y,ynovo,tol)
c               write(*,*)'zero em :',x2
               kcont=kcont+1
               if (kcont.eq.1) then
                  write(300,*)x1,x2
                  write(303,*)x1,-x2
c
                  write(200,*)-x1,x2
                  write(203,*)-x1,-x2
               else if (kcont.eq.2) then
                  write(302,*)x1,x2
                  write(304,*)x1,-x2
c
                  write(202,*)-x1,x2
                  write(204,*)-x1,-x2
               endif
            endif
            y=ynovo
         enddo
         if (Ncont.gt.kcont) then
            write(302,*)
            write(304,*)
            write(202,*)
            write(204,*)
         endif
         Ncont=kcont
         x=x+dbx
      enddo
      close(300)
      close(302)
      close(303)
      close(304)
      close(200)
      close(202)
      close(203)
      close(204)
      write(*,*)'fim da construcao da zvc'
      return
      end subroutine zvc   
c-----------------------------------------------------------------------
      subroutine fundc(y,efe,de)
      use m_sysconst
      use m_hill
      implicit real*8 (A-H,O-Z)
      
      r13=dsqrt((x1-xP1)**2.0d0+y**2.0d0)
      r23=dsqrt((x1-xP2)**2+y**2.0d0)

      exp1=2.0d0*y
      exp2=-0.5d0*(2.0d0*y)/(r13**3.0d0)
      exp3=-0.5d0*(2.0d0*y)/(r23**3.0d0)

      OMEGA=0.5d0*(x1**2.0d0+y**2.0d0)+((1.0d0-u)/r13)+
     & (u/r23)+0.5d0*u*(1.0d0-u)
      OMEGAPT=0.5d0*exp1+(1.0d0-u)*exp2+u*exp3

      efe=C-2.0*OMEGA
      de=-2*OMEGAPT

      return
      end subroutine fundc
      
c-----------------------------------------------------------------------
      function rtsafe(funcd,x1,x2,xacc)
      implicit real*8 (A-H,O-Z)
      INTEGER MAXIT
      DOUBLE PRECISION rtsafe,x1,x2,xacc
      external funcd
      parameter (MAXIT=100)
      INTEGER j
      DOUBLE PRECISION df,dx,dxold,f,fh,fl,temp,xh,xl
      call funcd(x1,fl,df)
      call funcd(x2,fh,df)
      if((fl.GT.0.d0.and.fh.GT.0.0d0).or.(fl.LT.0.d0.and.fh.LT.0.0d0))
     &     then
         write(*,*) 'root must be bracketed in rtsafe'
         read(*,*)bobo
      endif
      if(fl.EQ.0.0d0)then
        rtsafe=x1
        return
      else if(fh.EQ.0.d0)then
        rtsafe=x2
        return
      else if(fl.LT.0.d0)then
        xl=x1
        xh=x2
      else
        xh=x1
        xl=x2
      endif
      rtsafe=.5d0*(x1+x2)
      dxold=abs(x2-x1)
      dx=dxold
      call funcd(rtsafe,f,df)
      do 11 j=1,MAXIT
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).GE.0.0d0.or.
     & abs(2.d0*f).GT.abs(dxold*df) ) then
          dxold=dx
          dx=0.5d0*(xh-xl)
          rtsafe=xl+dx
          if(xl.EQ.rtsafe)return
        else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if(temp.EQ.rtsafe)return
        endif
        if(abs(dx).LT.xacc) return
        call funcd(rtsafe,f,df)
        if(f.LT.0.0d0) then
          xl=rtsafe
        else
          xh=rtsafe
        endif
11    continue
      write(*,*) 'rtsafe exceeding maximum iterations'
      read(*,*)bobo
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software 7L`2.d0 
c-----------------------------------------------------------------------
      
