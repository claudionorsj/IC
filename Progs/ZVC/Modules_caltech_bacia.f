c=======================================================================
      module m_sysconst
c=======================================================================
      implicit none
c     Modelo: CPR3BP
c     Configuracao escolhida: L2 P2 L1 P1 L3
c     P1: Terra (m1=1-mu)
c     P2: Lua (m2=mu)
c     Parametro de massa: mu=u
      real*8,parameter :: u=0.0121506683d0
c     Coordenadas dos primarios:
      real*8,parameter :: xP1=-u, xP2=1.0d0-u, yP1=0.0d0, yp2=0.0d0
c     Coordenadas dos pontos lagrangeanos:
      real*8,parameter :: xL1=0.83691471889320d0, yL1=0.0d0
      real*8,parameter :: xL2=1.15568248347861d0, yL2=0.0d0
      real*8,parameter :: xL3=-1.00506268026259d0, yL3=0.0d0
      real*8,parameter :: xL4=0.4878493317d0, yL4=0.866025403784439d0
      real*8,parameter :: xL5=0.4878493317d0, yL5=-0.86602540378443d0
c     Distancia entre os primarios em km:
      real*8,parameter :: dpri=384.4d3
c     Raio equatorial da Terra normalizado:
      real*8,parameter :: rP1=6378.136d0/dpri
c     Raio medio da Lua normalizado:
!      real*8,parameter :: rP2=1.0d-10
      real*8,parameter :: rP2=1738.0d0/dpri
c     Definicao da caixa de CIs
      real*8,parameter :: xi=xL1, xf=xL2, dx=0.001d0
      real*8,parameter :: yi=-0.3d0, yf=0.3d0, dy=0.0025d0
c     Definicao da caixa de gracao d curvaplano
      real*8,parameter :: yyi=-2.0d0, yyf=2.0d0, dyy=0.001d0  
      end module m_sysconst

c=======================================================================
      module m_hill
c=======================================================================
      implicit real*8 (A-H,O-Z)
      real*8 C,x1,v1,v2

      end module m_hill
c=======================================================================
