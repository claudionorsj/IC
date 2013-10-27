c=======================================================================
      program integrar
c=======================================================================
c     Integra condições iniciais no PR3CP e constroi secao de Poincare. 
c=======================================================================
      implicit real*8 (a-h,o-z)      ! Definindo variáveis reais
      parameter (ndim=4)             ! Definindo parametro (dimensao do espaco de fases)
c     Determinando as dimensoes dos vetores
      dimension x(ndim),xold(ndim),xsec(ndim),ci(ndim)
      dimension relerr(ndim),abserr(ndim),work(56) ! no minimo 14*ndim

      common/actejac/C0,C

      external deriv,derivsec ! Declarando subrotinas que serao argumentos do integrador

c     Abrindo arquivos de saida
      open(unit=1,file='traj.dat',status='unknown')
      open(unit=2,file='sec.dat',status='unknown')

c     Parametro de massa
      rmu=0.0121506683d0

c     Parametros do integrador
      do i=1,ndim; relerr(i)=1.0d-14; abserr(i)=1.0d-15; enddo
      dtmax=0.1d0; epsJac=1.d-12

c     Atribuindo o valor da Constante de Jacobi
      C0=3.1d0

c     Definindo a linha de cis e setando tempo final de integracao
      x1_ini=-1.0d0+rmu-0.1; x1_fim=rmu-0.2d0
      x1_step=abs(x1_fim-x1_ini)/1.0d1
      ci(1)=x1_ini; tf=50.0d0; sentido=1.0d0 ! Integracao direta (+1.0d0) ou retrograda (-1.0d0)

      do 1 while (ci(1).le.x1_fim)
         write(*,*)'ci',ci(1)
         ci(2)=0.0d0; ci(3)=0.0d0
         call Vely (rmu,ndim,ci,C0,kvel)
         if(kvel.ne.0)then; write(*,*)'v2<0!'; goto 2; endif
         do i=1,ndim; x(i)=ci(i); enddo ! Atribuindo valor inicial ao vetor de estado
         ti=0.0d0; call evolucao (rmu,ndim,x,ti,tf,sentido)
 2       ci(1)=ci(1)+x1_step
         write(1,*); write(2,*); write(3,*)
 1    continue

      stop; end program integrar ! Finalizando a rotina principal

c=======================================================================
      subroutine evolucao (rmu,neq,y0,ti,tf,sentido)
c=======================================================================
      implicit real*8 (A-H,O-Z)
      dimension y0(neq),y(neq),yold(neq)

      do i=1,neq; y(i)=y0(i); yold(i)=y(i); enddo
      t=ti; told=t+1.0d0 ! dummy diferente de t
      dt=sentido*0.001d0; nfuros=0

      call testar (rmu,neq,told,t,yold,y,nfuros,kopcao); tc=t

      do 10 while (kopcao.le.1)
         if(kopcao.eq.0)then
            call evoluinorm (rmu,neq,tf,t,told,y,yold,dt,
     &           nfuros,kopcao)
            tc=t
         else if (kopcao.eq.1) then ! integra com regularizacao
            call evoluiregu (rmu,neq,tf,tc,told,y,yold,dt,
     &           nfuros,kopcao)
         endif
 10   continue
      return; end subroutine evolucao

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
c     kopcao = 0 integracao sem regularizacao
c              1 integracao com regularizacao
c              5 colisao com o primario menor
c              6 colisao com o primario maior
c              4 cte de Jacobi nao conservada
c              3 tempo final atingido
c              7 alguma secao é criterio de parada
c              2 erro de integracao
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c=======================================================================
      subroutine testar (rmu,neq,told,time,yold,y,nfuros,kopcao)
c=======================================================================
      implicit real*8 (A-H,O-Z)
      dimension y(neq),yold(neq),ysec(neq)
      logical ok
      parameter (epsJac=1.0d-12,epsReg1=4.0d-2,epsReg2=1.0d-2)
      parameter (dpri=384.4d3,rP1=6378.136d0/dpri,rP2=1738.0d0/dpri)
      common/actejac/C0,C
      d1=dsqrt((y(1)-rmu)**2.0d0+y(2)**2.0d0)
      d2=dsqrt((y(1)-rmu+1.0d0)**2.0d0+y(2)**2.0d0)
!!!!!!Confere a cte de Jacobi!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      C=CteJac(rmu,neq,y); write(3,*)time,dabs(C-C0)
      if(dabs(C-C0).gt.epsJac)then
         if(kopcao.ne.1)then !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            write(*,*)'dC>eps. kopcao =',kopcao,'t =',time
            write(*,*)'C,dC =',C,dabs(C-C0)
            write(*,*)'d1= ',d1,' d2= ',d2 !; read(*,*)bobo
            kopcao=4; return
         endif !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      kopcao=0
      if (d2.gt.rP2.and.d2.le.epsReg2) then
         kopcao=1
      else if (d2.le.rP2) then
         write(*,*)'Colisao com o secundario. d2,rP2,t =',d2,rP2,time
         kopcao=5; return
      else if (d1.gt.rP1.and.d1.le.epsReg1) then
         kopcao=1
      else if (d1.le.rP1) then
         write(*,*)'Colisao com o primario. d1,rP1,t =',d1,rP1,time
         kopcao=6; return
      endif
      ok=.false.; call secao (rmu,neq,told,time,yold,y,tsec,ysec,ok)
      if(ok)then
         nfuros=nfuros+1
         write(2,12)(ysec(i),i=1,4),tsec ! Escrevendo os pontos na secao
         if(nfuros.eq.10)then; kopcao=7; return; endif
      endif
      write(1,12)(y(i),i=1,4),time ! Escrevendo a trajetoria
 12   format(5f24.16)
      return; end subroutine testar

c=======================================================================
      subroutine secao (rmu,neq,told,time,yold,y,tsec,ysec,ok)
c=======================================================================
      implicit real*8 (A-H,O-Z)
      dimension y(neq),yold(neq),ysec(neq)
      external derivsec
         if(yold(2).lt.0.0d0.and.y(2).ge.0.0d0)then
            do i=1,ndim; ysec(i)=yold(i); enddo
            passo=0.0d0-yold(2); tdum=0.0d0
            call rk4 (rmu,derivsec,neq,ysec,tdum,passo)
            tsec=told+tdum*passo; ok=.true.
         endif
      return; end subroutine secao

c=======================================================================
      subroutine evoluinorm (rmu,neq,tfim,time,told,y,yold,dt,
     &           nfuros,kopcao)
c=======================================================================
      implicit real*8 (A-H,O-Z)
      dimension y(neq),yold(neq)
      dimension relerr(neq),abserr(neq),work(56) ! no minimo 14*neq
      external deriv
      do i=1,neq; relerr(i)=1.0D-14; abserr(i)=1.0D-15; enddo
      iflag=1; tout=time+dt; dtmax=0.1d0
c-----
      do 100 while (dabs(time).le.dabs(tfim).and.kopcao.eq.0)
         do i=1,neq; yold(i)=y(i); enddo; told=time
         call rkf78(rmu,deriv,neq,y,time,tout,relerr,abserr,
     &        iflag,work,dt)
         if(iflag.ne.2)then 
            kopcao=2; write(*,11)'iflag=',iflag,' em t=',t; return
         else ! t=tout
            if(dt.lt.dtmax)then; tout=tout+dt; else
            tout=tout+dtmax; endif
            if(dabs(tout).gt.dabs(tfim))then
               if(dabs(time).lt.dabs(tfim))then
               tout=tfim; else; kopcao=3
            return; endif; endif
         endif
         call testar (rmu,neq,told,time,yold,y,nfuros,kopcao)
 100  continue
 11   format(a6,i4,a6,f6.3)
      return; end subroutine evoluinorm

c=======================================================================
      subroutine evoluiregu (rmu,neq,tfim,time,told,x,xold,dt,
     &           nfuros,kopcao)
c=======================================================================
      implicit real*8 (A-H,O-Z)
      dimension x(neq),xold(neq),y(5),yold(5)
      dimension relerr(5),abserr(5),work(70) ! no minimo 14*5=70
      dimension W(5),Q(5),Wold(5),DWDTAU(5)
      external derivs
c-----
      do i=1,neq; y(i)=x(i); yold(i)=xold(i); enddo
      y(5)=time; yold(5)=told
c-----
      do i=1,5; relerr(i)=1.0D-14; abserr(i)=1.0D-15; enddo
      iflag=1; tau=0.0d0; dtantes=dt; dt=dt*0.1d0; tout=tau+dt
      dtmax=dt*10.d0
c-----Levando para variaveis regularizadas------------------------------
      call XTOQ(Y,RMU,Q) ! Lembrar que C já deve ter sido
      call QTOW(Q,W)   ! calculado pelo menos uma vez
      call DERIVS(RMU,5,TAU,W,DWDTAU); FLW2=DWDTAU(5)
c-----------------------------------------------------------------------
      do 100 while (dabs(time).le.dabs(tfim).and.kopcao.eq.1)
         do i=1,5; yold(i)=y(i); Wold(i)=W(i); enddo
         told=time; tauold=tau; FLW2old=FLW2
         call rkf78(rmu,derivs,5,W,tau,tout,relerr,abserr,
     &   iflag,work,dt)
         if (iflag.ne.2) then
            kopcao=2; write(*,11)'iflag=',iflag,' em t=',t; return
         else
            if(dt.lt.dtmax)then; tout=tout+dt; else
            tout=tout+dtmax; endif
         endif    
c--------Trazendo para variaveis fisicas--------------------------------
         call DERIVS(RMU,5,TAU,W,DWDTAU); FLW2=DWDTAU(5)
!        em lugar de time=time+(tau-tauold)*flw2
         time=time+0.5d0*(tau-tauold)*(FLW2+FLW2old)
         call WTOQ(W,Q)
         call QTOX(Q,RMU,Y)
c-----------------------------------------------------------------------
         do i=1,neq; x(i)=y(i); xold(i)=yold(i); enddo
         if(dabs(time).gt.dabs(tfim))then; kopcao=3; return; endif
         call testar (rmu,neq,told,time,xold,x,nfuros,kopcao)
 100  continue
      dt=dtantes
 11   format(a6,i4,a6,f6.3)
      return; end subroutine evoluiregu

cccccccccccccccccccccccccccccccccccccccccccccccccc
c     E se usarmos uma subrotina do totfast?     c
c     call JACOBIW(W,U,Cw)                       c    
c     write(*,*)'C,Cw',C,Cw,dabs(C-Cw)           c
c     read(*,*)                                  c
c     OK! Cw confere com nosso C                 c
c     Se usarmos outra subrotina do totfast?     c
c     call JACOBI(W,U,C1)                        c
c     write(*,*)'C,C1',C,C1,dabs(C-C1)           c
c     read(*,*)                                  c
c     RUIM! difere um pouco do nosso C           c
cccccccccccccccccccccccccccccccccccccccccccccccccc

c=======================================================================
      real*8 function CteJac (u,neq,y)
c=======================================================================
c     Esta funcao fornece (na variavel real*8 CteJac) o valor da
c     Constante de Jacobi de uma particula no ponto y(1),y(2) com 
c     velocidade y(3),y(4)
c=======================================================================
      implicit real*8 (a-h,o-z)
      dimension y(neq)
      x1=y(1); x2=y(2); x3=y(3); x4=y(4)
      R1=dsqrt((x1-u)*(x1-u)+x2*x2)
      R2=dsqrt((x1+1.0d0-u)*(x1+1.0d0-u)+x2*x2)
      Omega=0.5d0*(x1*x1+x2*x2)+(1.0d0-u)/R1+u/R2+0.5d0*u*(1.0d0-u)
      CteJac=2.0d0*Omega-(x3*x3+x4*x4)
      return; end function CteJac

c=======================================================================
      subroutine Vely (u,neq,y,C,kvel)
c=======================================================================
c     Esta funçao fornece (na quarta componente do vetor de 
c     variaveis y, y(4)) o valor da segunda componente da velocidade de 
c     uma particula no ponto y(1),y(2) com a primeira componente da 
c     velocidade dada por y(3) no nivel de energia correspondente a C
c     Se o calculo de y(4) é bem sucedido iflag retorna 0. Senão, 
c     iflag retorna 1
c=======================================================================
      implicit real*8 (a-h,o-z)
      dimension y(neq)
      x1=y(1); x2=y(2); x3=y(3)
      R1=dsqrt((x1-u)*(x1-u)+x2*x2)
      R2=dsqrt((x1+1.0d0-u)*(x1+1.0d0-u)+x2*x2)
      Omega=0.5d0*(x1*x1+x2*x2)+(1.0d0-u)/R1+u/R2+0.5d0*u*(1.0d0-u)
      Vely2=2.0d0*Omega-x3*x3-C; if(Vely2.ge.0.0d0)then
      kVel=0; y(4)=sqrt(Vely2); else; kVel=1; endif
      return; end subroutine Vely

c=======================================================================
      subroutine deriv (u,neq,t,y,yprime)
c=======================================================================
      implicit real*8 (a-h,o-z)
      dimension y(neq),yprime(neq)
      x1=y(1); x2=y(2); x3=y(3); x4=y(4)
      R1=dsqrt((x1-u)*(x1-u)+x2*x2)
      R2=dsqrt((x1+1.0d0-u)*(x1+1.0d0-u)+x2*x2)
      Omega_x1=x1-(1.0d0-u)*(x1-u)/R1**3.0d0-
     &         u*(x1+1.0d0-u)/R2**3.0d0
      Omega_x2=x2-(1.0d0-u)*x2/R1**3.0d0-u*x2/R2**3.0d0
c-----sistema dinamico--------------------------------------------------
      yprime(1)=x3
      yprime(2)=x4
      yprime(3)=Omega_x1+2.0d0*x4
      yprime(4)=Omega_x2-2.0d0*x3
c-----------------------------------------------------------------------
      return; end subroutine deriv

c=======================================================================
      subroutine derivsec (u,neq,t,y,yprime)
c=======================================================================
      implicit real*8 (a-h,o-z)
      dimension y(neq),yprime(neq)

      x1=y(1); x2=y(2); x3=y(3); x4=y(4)

      R1=dsqrt((x1-u)*(x1-u)+x2*x2)
      R2=dsqrt((x1+1.0d0-u)*(x1+1.0d0-u)+x2*x2)

      Omega_x1=x1-(1.0d0-u)*(x1-u)/R1**3.0d0-
     &         u*(x1+1.0d0-u)/R2**3.0d0
      Omega_x2=x2-(1.0d0-u)*x2/R1**3.0d0-u*x2/R2**3.0d0

c-----sistema dinamico--------------------------------------------------
      yprime(1)=x3
      yprime(2)=x4
      yprime(3)=Omega_x1+2.0d0*x4
      yprime(4)=Omega_x2-2.0d0*x3
c-----------------------------------------------------------------------

      aux=1.0d0/yprime(2)
      do i=1,neq; yprime(i)=yprime(i)*aux; enddo; t=aux

      return; end subroutine derivsec

c=======================================================================
      subroutine RK4 (rmu,f,neq,y,t,h)
c=======================================================================
      implicit real*8 (a-h,o-z)
      dimension dydx(neq),y(neq),dym(neq),dyt(neq),yt(neq)
      hh=h*0.5d0; h6=h/6.0d0; th=t+hh
      call f(rmu,neq,t,y,dydx)
      do 11 i=1,neq
         yt(i)=y(i)+hh*dydx(i)
 11   continue
      call f(rmu,neq,th,yt,dyt)
      do 12 i=1,neq
         yt(i)=y(i)+hh*dyt(i)
 12   continue
      call f(rmu,neq,th,yt,dym)
      do 13 i=1,neq
         yt(i)=y(i)+h*dym(i)
 13   continue
      do 14 i=1,neq
         dym(i)=dyt(i)+dym(i)
 14   continue
      call f(rmu,neq,t+h,yt,dyt)
      do 15 i=1,neq
         y(i)=y(i)+h6*(dydx(i)+dyt(i)+(2.0d0*dym(i)))
 15	continue
      return; end subroutine RK4
