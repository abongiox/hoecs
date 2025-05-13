!
! This file is distributed under the terms of the GNU General Public 
! License: http://www.gnu.org/copyleft/gpl.txt .
!
      module constants
!
        real*8  , parameter :: zero=0.d0, one=1.d0, two=2.d0
!
        real*8 , parameter :: bohr  = 0.529177210544d0
        real*8 , parameter :: au2ev = 27.211386245981d0
        real*8 , parameter :: au2j  = 4.3597482d-18
        real*8 , parameter :: au2pas= au2j/bohr**3*1.d30
        real*8 , parameter :: au2gpa= au2pas/1.d9
        real*8 , parameter :: ry2gpa= au2gpa/two
!
      end module constants
!
!
      module controls 
!
        character fname*255,dname*255,oname*255
        logical tdeform,tpk2ecs,tout,tonlyi,tperm,thomo,tone
        integer iecs,ndef,nvec,ione
        real*8 pmst, numthr
!
        integer , allocatable :: idef(:,:), grid(:,:), isym(:,:)
!
        integer , parameter :: oinp=1, oist=2, oxyz=3, oout=100
        integer , parameter :: opio=11, oecs=10, oall=12, ocom=13, ohom=23
!
      end module controls 
!
      module crystal
!
        integer natom
        real*8 box(3,3),boxit(3,3),vol
        character*2, allocatable :: id(:)
        real*8, allocatable :: x(:,:),xf(:,:)
!
        logical tcry(48)
        real*8 sym(3,3,48) 
        integer nsym,ncry,sabc(3,3,48),imat(48)
        character sname(48)*45
!
      end module crystal
!
      module permutations
!
        implicit none
        integer nxp,np
        character , allocatable :: perm(:)*8
        integer , allocatable :: ind(:,:)
!
      end module permutations
!
      module tensors
        implicit none
!
        real*8 , parameter :: ectoll=1.d0
        integer iref
        integer nc2,nc3,nc4,nc5,nc6,nc7,nc8
        logical tc2,tc3,tc4,tc5,tc6,tc7,tc8
!
        integer, allocatable :: imap(:,:)
        real*8 , allocatable :: cau(:,:,:),pio(:,:,:),cstr(:,:)
!
        real*8 , allocatable :: c2(:,:,:),c3(:,:,:,:),c4(:,:,:,:,:)
        real*8 , allocatable :: c5(:,:,:,:,:,:), c6(:,:,:,:,:,:,:)
        real*8 , allocatable :: c7(:,:,:,:,:,:,:,:)
        real*8 , allocatable :: c8(:,:,:,:,:,:,:,:,:)
!
        logical, allocatable :: t3(:,:,:),t4(:,:,:,:),t5(:,:,:,:,:)
        logical, allocatable :: t6(:,:,:,:,:,:), t7(:,:,:,:,:,:,:)
        logical, allocatable :: t8(:,:,:,:,:,:,:,:)
!
        real*8 , allocatable :: c2r(:,:),c3r(:,:,:),c4r(:,:,:,:)
        real*8 , allocatable :: c5r(:,:,:,:,:),c6r(:,:,:,:,:,:)
        real*8 , allocatable :: c7r(:,:,:,:,:,:,:),c8r(:,:,:,:,:,:,:,:)
!
!
        real*8 c2p(6,6), c2a(6,6)
        real*8 c3p(6,6,6), c3a(6,6,6)
        real*8 c4p(6,6,6,6), c4a(6,6,6,6)
        real*8 c5p(6,6,6,6,6), c5a(6,6,6,6,6)
        real*8 c6p(6,6,6,6,6,6), c6a(6,6,6,6,6,6)
        real*8 c7p(6,6,6,6,6,6,6), c7a(6,6,6,6,6,6,6)
        real*8 c8p(6,6,6,6,6,6,6,6), c8a(6,6,6,6,6,6,6,6)
!
      end module tensors
!
!
!-----------------------------------------------------------------------
      program hoecs
!-----------------------------------------------------------------------
      implicit none
!
!     definition of the finite-difference operator
      call finite_difference_operator
!
!     read control/input parameters as arguments
      call read_arguments
!
!     the reference state
      call read_reference_state
!
!     find the symmetry operations 
      call symmetry_operations
!
!     apply the strain and read the stress, if present
      call strain_stress_vectors
!
!     caculate the HOECs
      call divided_differences
!
!     symmetrization
      call impose_material_symmetries
      stop
      end
!-----------------------------------------------------------------------
      subroutine finite_difference_operator
!-----------------------------------------------------------------------
      use controls
      implicit none
!
!     first-order central finite difference
      ndef=12
      allocate(idef(6,ndef))
      idef(:,:)=0
!
      idef(1,1) = 1 ; idef(1,2) =-1
      idef(2,3) = 1 ; idef(2,4) =-1
      idef(3,5) = 1 ; idef(3,6) =-1
      idef(4,7) = 1 ; idef(4,8) =-1
      idef(5,9) = 1 ; idef(5,10)=-1
      idef(6,11)= 1 ; idef(6,12)=-1
!
      return
      end
!-----------------------------------------------------------------------
      subroutine single_deformation_mode
!-----------------------------------------------------------------------
      use constants
      use controls
      implicit none
      integer i,k,it(6),itdef(6)
      real*8 a,ide(3,3),tmp(3,3)
!
!
      if (tpk2ecs) goto 56
!
      itdef(:)=0
      if (ione .eq.-1) then
        itdef(1)=1 ; itdef(2)=1 ; itdef(3)=1
      else if (ione.ge.1.and.ione.le.6) then
        itdef(ione)=1
      else if (ione.eq.7) then
        itdef(1)=1 ; itdef(2)=1
      else if (ione.eq.8) then
        itdef(2)=1 ; itdef(3)=1
      else
        write(6,*) 'Wrong strain identifier!'
        stop
      endif
!
      ide(:,:)=zero
      ide(1,1)=one
      ide(2,2)=one
      ide(3,3)=one
!
      k=0
      it(:)=0
      do i=-iecs,iecs
        k=k+1
        it(:)=i*itdef(:)
        call deformation(k,it,i,1,ide,tmp,tmp)
      enddo
      return
!
56    continue
      stop
!
      return
      end
!-----------------------------------------------------------------------
      subroutine read_reference_state
!-----------------------------------------------------------------------
      use constants
      use controls
      use crystal
      implicit none
      integer i,j,k,n,iux
      real*8 a
      character ida*2,line*255
      logical matches
!
!
      open(1,file=trim(fname),iostat=i,status='old')
      if (i.ne.0) stop 328
      read(1,'(A)',end=100,err=100)line
      if (.not.matches('CELL_PARAMETERS',line)) then
        write(6,*)'something is wrong ...'
        write(6,*)line
        stop
      endif
      a=bohr
      if (matches('ANGSTROM',line).or.matches('angstrom',line)) a=one
      if (matches('Angstrom',line)) a=one
      do i=1,3
        read(1,*,end=100,err=100) (box(k,i),k=1,3)
        do k=1,3
          box(k,i)=box(k,i)*a
        enddo
      enddo
!
      read(1,'(A)',end=100,err=100) line
      if (.not.matches('ATOMIC_POSITIONS',line)) then
        write(6,*)'something is wrong ...'
        write(6,*)line
        stop
      endif
      n=0
10    read(1,*,end=99,err=99) ida
      n=n+1
      goto 10
99    natom=n
      rewind(1)
      do i=1,5
        read(1,*,end=100,err=100)
      enddo
!
      a=bohr
      if (matches('ANGSTROM',line)) a=one
      if (matches('angstrom',line)) a=one
      if (matches('Angstrom',line)) a=one
      if (matches('crystal',line)) a=one
      if (matches('Crystal',line)) a=one
      if (matches('CRYSTAL',line)) a=one
!
      allocate(id(natom),x(3,natom),xf(3,natom))
      do i=1,natom
        read(1,*,end=100,err=100) id(i),(x(k,i),k=1,3)
        do k=1,3
          x(k,i)=x(k,i)*a ! to A if iux=0
        enddo
      enddo
      close(1)
!
      call lattices(box,vol,boxit)
!
      iux=1
      if (matches('crystal',line)) iux=2
      if (matches('Crystal',line)) iux=2
      if (matches('CRYSTAL',line)) iux=2
      if (iux.eq.2) then
        do i=1,3*natom
          xf(i,1)=x(i,1)
          x(i,1)=zero
        enddo
        call fract2cart(box,xf,x,natom)
      else
        call cart2fract(boxit,x,xf,natom)
      endif
!
!
      return
100   stop 100
      end
!-----------------------------------------------------------------------
      subroutine symmetry_operations
!-----------------------------------------------------------------------
      use controls
      use crystal
      implicit none
      integer i,j,k,ou
!
      call set_sym_bl
      call subgroup(nsym,sabc,natom,xf,id,box,tcry)
      ncry=0
      do i=1,nsym
        if (tcry(i)) ncry=ncry+1
      enddo
!
      if (tout) then
        ou=16
        open(ou,file=trim(dname)//'symmetries')
        write(ou,'(2(1x,i3))') nsym,ncry
        do i=1,nsym
          write(ou,'(1x,i3,1x,l1,1x,a45)') i,tcry(i),sname(i)
          do j=1,3
            write(ou,*) (sym(j,k,i),k=1,3)!,(sabc(j,k,i),k=1,3)
          enddo
        enddo
        close(ou)
      endif
!
      return
      end
!-----------------------------------------------------------------------
      subroutine strain_stress_vectors
!-----------------------------------------------------------------------
      use constants
      use controls
      use crystal
      use tensors
      implicit none
      integer ii(6),idist,necs
      integer i,j,k,l,m,n,ou,nt,ji,is
      real*8 a,ca(3,3),pk(3,3),c(3,3),tmp(3,3),s(3,3)
      character com*255
!
!
      n=1
      do i=1,iecs-1
        n=n*ndef
      enddo
      n=n+1
      allocate(grid(6,n))
      do i=1,6*n
        grid(i,1)=0
      enddo
!
      n=1                  ! reference state
      nt=1
      necs=1
1     necs=necs+1
      do i=1,n
        do l=1,ndef
          do j=1,6
            ii(j)=grid(j,i)+idef(j,l)
          enddo
          do k=1,nt
            j=idist(ii,grid(1,k),6)
            if (j.eq.0) goto 2
          enddo
          nt=nt+1
          do k=1,6
            grid(k,nt)=ii(k)
          enddo
2         continue
        enddo
      enddo
      n=nt
      if (necs.lt.iecs) goto 1
!
      nvec=n
      allocate(isym(2,n))
      call symmetry_redundancy
!
!
      if (tpk2ecs) then
        if (tout) open(opio,file=trim(dname)//'stress.dat')
        allocate(cau(3,3,nvec),pio(3,3,nvec))
      endif
!
!
      do i=1,n
        ji=isym(1,i)
        is=isym(2,i)
!
        if (tonlyi) then
          if (ji.eq.i.and.is.eq.1) call deformation(i,grid(1,i),ji,is,sym(1,1,is),ca,pk)
        else
          call deformation(i,grid(1,i),ji,is,sym(1,1,is),ca,pk)
        endif
!
        if (tpk2ecs) then
!
          if (ji.eq.i.and.is.eq.1) then
            cau(:,:,i)=ca(:,:)      
            pio(:,:,i)=pk(:,:)      
            com=' :irreducible member'
          else
            tmp(:,:)=zero
            c(:,:)=zero
            s(:,:)=sym(:,:,is)
            call dgemm('t','n',3,3,3,one,s,3,cau(1,1,ji),3,zero,tmp,3)
            call dgemm('n','n',3,3,3,one,tmp,3,s,3,zero,c,3)
!
            if (.not.tonlyi) then
              a=zero
              do m=1,9
                a=a+abs(ca(m,1)-c(m,1))
              enddo
              if (a.gt.numthr) then
                write(6,*) 'Warning ...'
                write(6,*) i,ji,is,a
                do m=1,3
                  write(6,77) (ca(m,l),l=1,3),(c(m,l),l=1,3)
77                format(3(1x,f8.5),4x,3(1x,f8.5))
                enddo
              endif
              cau(:,:,i)=ca(:,:)      
              pio(:,:,i)=pk(:,:)      
              com=' :redundant, stress read from file!'
            else
              cau(:,:,i)=c(:,:)      
              tmp(:,:)=zero
              c(:,:)=zero
              call dgemm('t','n',3,3,3,one,s,3,pio(1,1,ji),3,zero,tmp,3)
              call dgemm('n','n',3,3,3,one,tmp,3,s,3,zero,c,3)
              pio(:,:,i)=c(:,:)
              com=' :redundant, stress obtained from sym. op.!'
            endif
          endif

          if (tout) then
            write(opio,'(a1)')'-'
            write(opio,*)i,ji,is,trim(com)
            write(opio,96)(grid(l,i),l=1,6)
            do l=1,3
              write(opio,99)(cau(l,k,i),k=1,3),(pio(l,k,i),k=1,3)
            enddo
96          format(6(1x,i3))
99          format(3(1x,f10.4),3x,3(1x,f10.4))
          endif
!        
        endif
      enddo
      if (tpk2ecs.and.tout) close(opio)
!
!
      return
      end
!-----------------------------------------------------------------------
      subroutine symmetry_redundancy
!-----------------------------------------------------------------------
      use constants
      use controls
      use crystal
      implicit none
      integer ou,nx
      integer i,j,k,l,m,n
      real*8 tmp(3,3),ost(3,3),vo(6)
      real*8 , allocatable :: str(:,:,:)
      logical , allocatable :: tsym(:)
!
!
      nx=nvec
      allocate(str(3,3,nx),tsym(nx))
      do i=1,nx
        vo(:)=pmst*dble(grid(:,i))
        call voigt2str(vo,str(1,1,i))
        tsym(i)=.false.
        isym(1,i)=0
        isym(2,i)=0
      enddo
!
      ou=12
      if (tout) open(ou,file=trim(dname)//'redundacies')
      do i=1,nx
!
        if (tsym(i)) goto 30
        tsym(i)=.true.
        isym(1,i)=i
        isym(2,i)=1        
!
        do j=i+1,nx
!
          if (tsym(j)) goto 35
          do k=1,nsym
            if (.not.tcry(k)) goto 10  ! only the relevant symmetry operations
            tmp(:,:)=zero 
            ost(:,:)=zero 
           ! AT * Strain
            call dgemm('t','n',3,3,3,one,sym(1,1,k),3,str(1,1,i),3,zero,tmp,3)
           ! Strain * A
            call dgemm('n','n',3,3,3,one,tmp,3,sym(1,1,k),3,zero,ost,3)
           ! At * S * A = S ?
!
            do l=1,9
              if (abs(str(l,1,j)-ost(l,1)).gt.numthr) goto 10
            enddo
            tsym(j)=.true.
            isym(1,j)=i
            isym(2,j)=k
!
            if (tout) then
              write(ou,*) i,'=ref',j,'=to ref throu sym op.',k
              do m=1,3
                write(ou,'(3(1x,f7.4),3x,3(1x,f7.4),3x,3(1x,f7.4))') &
     &          (str(m,l,i),l=1,3),(str(m,l,j),l=1,3),(ost(m,l),l=1,3)
              enddo
            endif
!
            goto 35
10          continue ! not this symmetry operation
          enddo
35        continue ! already redundant
!
        enddo
30      continue   ! already redundant
      enddo
!
!     number of irriducible strained configurations
      n=0
      do i=1,nx
        k=isym(2,i)
        if (k.eq.1) n=n+1
      enddo
!
      if (tout) then
      write(ou,*) '-------------------------------------------'
      write(ou,*) 'Number of irriducible deformations: ',n
      write(ou,*) '-------------------------------------------'
      do i=1,nx
        k=isym(2,i)
        write(ou,100) i,isym(1,i),k,sname(k)
100     format(1x,i7,2x,i5,2x,i3,2x,a45)
      enddo
      close(ou)
      endif
!
      deallocate(tsym,str)
      return
      end
!-----------------------------------------------------------------------
      subroutine deformation(itag,ivoigt,isym1,isym2,symop,ca,pk)
!-----------------------------------------------------------------------
      use constants
      use controls
      use crystal
      use tensors , only : iref
      implicit none
      integer i,j,k,l,i1,i2,io,ivoigt(6),itag,isym1,isym2
      integer idist,ii(6)
      real*8 vo(6),dbox(3,3),symop(3,3),a,ca(3,3),pk(3,3),ins(3,3)
      character tag*6
!
!
!
      ca(:,:)=zero
      pk(:,:)=zero
      ins(:,:)=zero
!
      write(tag,'(i6.6)') itag
!
      ii(:)=0
      k=idist(ii,ivoigt,6)
      if (k.eq.0) iref=itag
!
      do k=1,6
        vo(k)=dble(ivoigt(k))*pmst
      enddo
      i=1
      if (tout) i=oout
      call strain2dbox(i,vo,box,dbox)
!
!
      if (tdeform) then
        open(oist,file=trim(dname)//'istr.'//tag,status='new',iostat=io)
        if (io.ne.0) then
          write(6,*) 'Cannot open file: istr.'//tag//'. Aborting.'
          stop
        endif
        write(oist,25) iecs,itag,(ivoigt(i),i=1,6),pmst,isym1,isym2
25      format(1x,i3,1x,i7,2x,6(1x,i3),2x,f13.7,2x,i7,1x,i2)
        do i=1,3
          write(oist,*) (symop(i,k),k=1,3)
        enddo
        close(oist)
      else
        open(oist,file=trim(dname)//'istr.'//tag,status='old',iostat=io)
        if (io.ne.0) then
          write(6,*) 'Cannot open file: istr.'//tag//'. Aborting.'
          stop
        endif
        read(oist,*,err=902) i,j,(ii(i),i=1,6),a,i1,i2
        k=idist(ii,ivoigt,6)
        if (a.ne.pmst.or.i.ne.iecs.or.j.ne.itag.or.k.ne.0.or.i1.ne.isym1.or.i2.ne.isym2) then
          write(6,*) 'Inconsistency in istr.'//tag//'. Aborting.'
          stop
        endif
        do i=1,3
          read(oist,*,err=902) (ins(i,k),k=1,3)
        enddo
        a=zero
        do l=1,3*3
          a=a+abs(symop(l,1)-ins(l,1))
        enddo
        if (a.gt.numthr) then
          write(6,*) 'Inconsistent symmetry op. in istr.'//tag//'. Aborting.'
          write(6,*) 'Expected sym. op.:'
          do i=1,3
            write(6,*)(symop(i,k),k=1,3)
          enddo
          write(6,*) 'Read from file:'
          do i=1,3
            write(6,*)(ins(i,k),k=1,3)
          enddo
          stop
        endif
        do l=1,3
          read(oist,*,err=902,end=902) (ca(l,k),k=1,3)
        enddo
        ca(:,:)=ca(:,:)*ry2gpa
        close(oist)
!
        call cauchy2piola(vo,ca,pk)
!
      endif
!
!
!
!
      if (tpk2ecs) return
      open(oinp,file=trim(dname)//'ipos.'//tag)
      write(oinp,10)'CELL_PARAMETERS (bohr)'
10    format(a22)
      do i=1,3
        write(oinp,*) (dbox(k,i)/bohr,k=1,3)
      enddo
      write(oinp,11)'ATOMIC_POSITIONS (crystal)'
11    format(a26)
      do i=1,natom
        write(oinp,*) id(i),(xf(k,i),k=1,3)
      enddo
      close(oinp)
!
!
!
      if (.not.tout) return
      open(oxyz,file=trim(dname)//'xyz.'//tag)
      call fract2cart(dbox,xf,x,natom)
      write(oxyz,*) natom
      write(oxyz,25) (ivoigt(k),k=1,6)
      do i=1,natom
        write(oxyz,16)id(i),(x(k,i),k=1,3)
16      format(a2,3(1x,f12.5))
      enddo
      close(oxyz)
      return
902   write(6,*) 'Reading error in istr.'//tag//'. Aborting.'
      stop
      end
!-----------------------------------------------------------------------
      subroutine defgrad2voigt(br,b,vo)
!-----------------------------------------------------------------------
      use constants
      implicit none
      real*8 vo(6),br(3,3),b(3,3),brit(3,3)
      real*8 tmp(3,3),cho(3,3),vol
      integer i,j,k
!
!
      call lattices(br,vol,brit)
      call dgemm('n','t',3,3,3,one,b,3,brit,3,zero,cho,3)
!
      do i=1,3*3
        tmp(i,1)=zero
      enddo
      call dgemm('t','n',3,3,3,one,cho,3,cho,3,zero,tmp,3)
      do i=1,3
        tmp(i,i)=tmp(i,i)-one
      enddo
      do i=1,9
        tmp(i,1)=tmp(i,1)/two
      enddo
      do i=1,3
        vo(i)=tmp(i,i)
      enddo
      vo(4)=tmp(2,3)*two
      vo(5)=tmp(1,3)*two
      vo(6)=tmp(1,2)*two
!
      return
      end
!-----------------------------------------------------------------------
      subroutine voigt2str(v,str)
!-----------------------------------------------------------------------
      use constants
      implicit none
      real*8 str(3,3),v(6)
!
      str(1,1)=v(1)
      str(2,2)=v(2)
      str(3,3)=v(3)
      str(1,2)=v(6)/two
      str(2,1)=str(1,2)
      str(1,3)=v(5)/two
      str(3,1)=str(1,3)
      str(2,3)=v(4)/two
      str(3,2)=str(2,3)
!
      return
      end
!-----------------------------------------------------------------------
      subroutine strain2dbox(ou,vo,box,dbox)
!-----------------------------------------------------------------------
      use constants
      implicit none
      integer , parameter :: lw=20
      real*8 vo(6),cho(3,3)
      integer i,j,k,ou
      real*8 str(3,3),mat(3,3),chot(3,3),tmp(3,3),w(3),wrk(lw)
      real*8 dbox(3,3),box(3,3)
!
!
      if (ou.gt.6) then
        write(ou,*)'#'
        write(ou,*)'voigt'
        write(ou,12)(vo(i),i=1,6)
12      format(6(1x,f12.6))
      endif

      str(1,1)=vo(1)
      str(2,2)=vo(2)
      str(3,3)=vo(3)
      str(1,2)=vo(6)/two
      str(2,1)=str(1,2)
      str(1,3)=vo(5)/two
      str(3,1)=str(1,3)
      str(2,3)=vo(4)/two
      str(3,2)=str(2,3)
!
      if (ou.gt.6) then
        write(ou,*)'strain'
        do i=1,3
          write(ou,12) (str(k,i),k=1,3)
        enddo
      endif
!
      do i=1,9
        cho(i,1)=two*str(i,1)
      enddo
      do i=1,3
        cho(i,i)=cho(i,i)+one
      enddo
!
      call chole(3, cho, mat)
!
!     traspose
!
      do i=1,3
        w(i)=zero
        do j=1,3
          cho(i,j)=mat(j,i)
          tmp(i,j)=zero
          chot(i,j)=zero
        enddo
      enddo
!
!     single value decomposition
!
      call dgesvd('A','A',3,3,cho,3,w,tmp,3,chot,3,wrk,lw,i)
      if (i.ne.0) stop 6780
      do i=1,9 
        tmp(i,1)=zero
        mat(i,1)=zero
      enddo
!
!     polar decomposition
!
      do i=1,3
        mat(i,i)=w(i)
      enddo
      call dgemm('t','n',3,3,3,one,chot,3,mat,3,zero,tmp,3)
      do i=1,9 
        mat(i,1)=zero
      enddo
      call dgemm('n','n',3,3,3,one,tmp,3,chot,3,zero,mat,3)
      do i=1,3*3
        cho(i,1)=mat(i,1)
      enddo
!
      if (ou.gt.6) then
        write(ou,*)'deformation gradient'
        do i=1,3
          write(ou,12) (cho(k,i),k=1,3)
        enddo
      endif
!
!     make sure everything is ok
!
      do i=1,3*3
        tmp(i,1)=zero
      enddo
      call dgemm('t','n',3,3,3,one,cho,3,cho,3,zero,tmp,3)
      do i=1,3
        tmp(i,i)=tmp(i,i)-one
      enddo
      do i=1,9
        tmp(i,1)=tmp(i,1)/two
      enddo
!
      if (ou.gt.6) then
        write(ou,*)'strain, from def grad'
        do i=1,3
          write(ou,12) (tmp(k,i),k=1,3)
        enddo
      endif
!
      dbox(:,:)=zero
      call dgemm('n','n',3,3,3,one,cho,3,box,3,zero,dbox,3)
!
      return
      end
!-----------------------------------------------------------------------
      subroutine lattices(box,omega,boxi)
!-----------------------------------------------------------------------
      use constants
      implicit none
      real*8 , parameter :: toll=1.d-4
      real*8 box(3,3),alat
      real*8 a1(3),a2(3),a3(3),omega,boxi(3,3)
      real*8 b1(3),b2(3),b3(3),den,s,ide(3,3)
      integer i,j,k,l,iperm,ir
!
!
      alat=zero
      do k=1,3
        a1(k)=box(k,1)
        alat=alat+a1(k)**2
        a2(k)=box(k,2)
        a3(k)=box(k,3)
      enddo
      alat=sqrt(alat)
!
      den=zero
      i=1
      j=2
      k=3
      s=one
1     continue
      do iperm=1,3
         den=den+s*a1(i)*a2(j)*a3(k)
         l=i
         i=j
         j=k
         k=l
      enddo
      i=2
      j=1
      k=3
      s=-s
      if (s.lt.zero) goto 1
      i=1
      j=2
      k=3
      omega=abs(den)
      den=one/omega
!
      do ir=1,3
        b1(ir)=den*(a2(j)*a3(k)-a2(k)*a3(j))
        b2(ir)=den*(a3(j)*a1(k)-a3(k)*a1(j))
        b3(ir)=den*(a1(j)*a2(k)-a1(k)*a2(j))
        l=i
        i=j
        j=k
        k=l
      enddo
!
      do i=1,3
        boxi(1,i)=b1(i)
        boxi(2,i)=b2(i)
        boxi(3,i)=b3(i)
      enddo
!
      do i=1,3*3
        ide(i,1)=zero
      enddo
      call dgemm('n','n',3,3,3,one,box,3,boxi,3,zero,ide,3)
      k=0
      do i=1,3
        if (abs(ide(i,i)-one).gt.toll) k=1
        do j=i+1,3
          if (abs(ide(i,j)).gt.toll) k=1
          if (abs(ide(j,i)).gt.toll) k=1
        enddo
      enddo
      if (k.eq.0) return
!
!
      write(6,*)'BOX'
      do l=1,3
        write(6,12) (box(j,l),j=1,3)
      enddo
      write(6,*)'BOXI'
      do l=1,3
        write(6,12) (boxi(j,l),j=1,3)
      enddo
      write(6,*)'Identity'
      do l=1,3
        write(6,12) (ide(j,l),j=1,3)
      enddo
      stop
12    format(3(1x,f12.5))
      return
      end
!-----------------------------------------------------------------------
      function pbcr(r1,r2,box,boxi,rout)
!-----------------------------------------------------------------------
      use constants
      implicit none
      real*8 pbcr,box(3,3),boxi(3,3),rout(3),r
      real*8 r1(3),x1,y1,z1,r2(3),x2,y2,z2,x,y,z
!
      x1 = boxi(1,1)*r1(1)+boxi(2,1)*r1(2)+boxi(3,1)*r1(3)
      y1 = boxi(1,2)*r1(1)+boxi(2,2)*r1(2)+boxi(3,2)*r1(3)
      z1 = boxi(1,3)*r1(1)+boxi(2,3)*r1(2)+boxi(3,3)*r1(3)
!
      x2 = boxi(1,1)*r2(1)+boxi(2,1)*r2(2)+boxi(3,1)*r2(3)
      y2 = boxi(1,2)*r2(1)+boxi(2,2)*r2(2)+boxi(3,2)*r2(3)
      z2 = boxi(1,3)*r2(1)+boxi(2,3)*r2(2)+boxi(3,3)*r2(3)
!
      x=mod(x1-x2,one)
      x=x-dint(x*two)
      y=mod(y1-y2,one)
      y=y-dint(y*two)
      z=mod(z1-z2,one)
      z=z-dint(z*two)
!
      rout(1) = x*box(1,1)+y*box(1,2)+z*box(1,3)
      rout(2) = x*box(2,1)+y*box(2,2)+z*box(2,3)
      rout(3) = x*box(3,1)+y*box(3,2)+z*box(3,3)
      r=sqrt(rout(1)**2+rout(2)**2+rout(3)**2)
      pbcr=r
      return
      end
!-----------------------------------------------------------------------
      subroutine born_von_karman(box,r1,f1,boxi)
!-----------------------------------------------------------------------
      use constants
      implicit none
      real*8 box(3,3),boxi(3,3),r1(3),f1(3),x1,y1,z1,x,y,z
!
      x1 = boxi(1,1)*r1(1)+boxi(2,1)*r1(2)+boxi(3,1)*r1(3)
      y1 = boxi(1,2)*r1(1)+boxi(2,2)*r1(2)+boxi(3,2)*r1(3)
      z1 = boxi(1,3)*r1(1)+boxi(2,3)*r1(2)+boxi(3,3)*r1(3)
!
      x=mod(x1,one)
      x=x-(sign(one,x)-one)/two
      y=mod(y1,one)
      y=y-(sign(one,y)-one)/two
      z=mod(z1,one)
      z=z-(sign(one,z)-one)/two
!
      r1(1) = x*box(1,1)+y*box(1,2)+z*box(1,3)
      r1(2) = x*box(2,1)+y*box(2,2)+z*box(2,3)
      r1(3) = x*box(3,1)+y*box(3,2)+z*box(3,3)
      f1(1) = x
      f1(2) = y
      f1(3) = z
!
      return
      end
!-----------------------------------------------------------------------
      subroutine incell(box,x,xf,n,boxi)
!-----------------------------------------------------------------------
      implicit none
      integer i,n
      real*8 x(3,n),xf(3,n),box(3,3),boxi(3,3)
      do i=1,n
        call born_von_karman(box,x(1,i),xf(1,i),boxi)
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine fract2cart(box,xf,x,n)
!-----------------------------------------------------------------------
      implicit none
      integer i,n
      real*8 x(3,n),xf(3,n),box(3,3),a,b,c
      do i=1,n
        a=xf(1,i)
        b=xf(2,i)
        c=xf(3,i)
        x(1,i) = a*box(1,1)+b*box(1,2)+c*box(1,3)
        x(2,i) = a*box(2,1)+b*box(2,2)+c*box(2,3)
        x(3,i) = a*box(3,1)+b*box(3,2)+c*box(3,3)
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine cart2fract(boxi,x,xf,n)
!-----------------------------------------------------------------------
      use constants
      implicit none
      integer i,n
      real*8 x(3,n),xf(3,n),boxi(3,3),x1,y1,z1,a,b,c
      do i=1,n
        a=x(1,i)
        b=x(2,i)
        c=x(3,i)
        x1 = boxi(1,1)*a+boxi(2,1)*b+boxi(3,1)*c
        y1 = boxi(1,2)*a+boxi(2,2)*b+boxi(3,2)*c
        z1 = boxi(1,3)*a+boxi(2,3)*b+boxi(3,3)*c
        xf(1,i)=mod(x1,one)
        xf(1,i)=xf(1,i)-(sign(one,xf(1,i))-one)/two
        xf(2,i)=mod(y1,one)
        xf(2,i)=xf(2,i)-(sign(one,xf(2,i))-one)/two
        xf(3,i)=mod(z1,one)
        xf(3,i)=xf(3,i)-(sign(one,xf(3,i))-one)/two
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine invmat (n, a, a_inv, da)
!-----------------------------------------------------------------------
! computes the inverse "a_inv" of matrix "a", both dimensioned (n,n)
! if the matrix is dimensioned 3x3, it also computes determinant "da"
! matrix "a" is unchanged on output - LAPACK
!
!USE kinds
      implicit none
      integer :: n
      real*8, DIMENSION (n,n) :: a, a_inv
      real*8 :: da
!
      integer :: info, lda, lwork, ipiv (n)
! info=0: inversion was successful
! lda   : leading dimension (the same as n)
! ipiv  : work space for pivoting (assumed of length lwork=n)
      real*8 :: work (n)
! more work space
!
      lda = n
      lwork=n
!
      a_inv(:,:) = a(:,:)
!
      call DGETRF (n, n, a_inv, lda, ipiv, info)
      if (info.ne.0) stop 21 !call errore ('invmat', 'error in DGETRF', abs (info) )
      call DGETRI (n, a_inv, lda, ipiv, work, lwork, info)
      if (info.ne.0) stop 22 !call errore ('invmat', 'error in DGETRI', abs (info) )
!
      if (n == 3) then
        da = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) + &
             a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) + &
             a(1,3)*(a(2,1)*a(3,2)-a(3,1)*a(2,2))
        IF (ABS(da) < 1.d-10) stop 23 !CALL errore(' invmat ',' singular matrix ', 1)
      else
        da = 0.d0
      end if

      return
      end subroutine invmat
!
!-----------------------------------------------------------------------
      function matches( string1, string2 )
!-----------------------------------------------------------------------
!
! ... .TRUE. if string1 is contained in string2, .FALSE. otherwise
!
      IMPLICIT NONE
!
      CHARACTER (LEN=*), INTENT(IN) :: string1, string2
      LOGICAL                       :: matches
      INTEGER                       :: len1, len2, l
!
!
      len1 = LEN_TRIM( string1 )
      len2 = LEN_TRIM( string2 )
      DO l = 1, ( len2 - len1 + 1 )
        IF ( string1(1:len1) == string2(l:(l+len1-1)) ) THEN
          matches = .TRUE.
          RETURN
        END IF
      END DO
      matches = .FALSE.
      RETURN
      end
!-----------------------------------------------------------------------
      subroutine subgroup(nrot,s,nat,xau,ityp,at,sym)
!-----------------------------------------------------------------------
!     Given the point group of the Bravais lattice, this routine finds 
!     the subgroup which is the point group of the considered crystal.

!     Non symmorphic groups are allowed, provided that fractional
!     translations are allowed (nofrac=.false), that the unit cell is
!     not a supercell, and that they are commensurate with the FFT grid
!
!     On output, the array sym is set to .true.. for each operation
!     of the original point group that is also a symmetry operation
!     of the crystal symmetry point group
!
      implicit none
      integer nrot, nat
      real*8 xau(3,nat),at(3,3)
      character*2 ityp(nat)
      integer s(3,3,48)
!
      integer irt (48, nat)
      logical sym (48)
      integer na, nb, irot, i, j, k
      real*8 , allocatable :: rau(:,:)
      logical fractional_translations,ttag
      real*8 ft (3)
!
      allocate(rau(3,nat))
      do i=1,48
        sym(i)=.false.
      enddo
!
      nb = 1
      irot = 1
!
      ttag=.true.
      fractional_translations = .true.
      do na = 2, nat
        if ( fractional_translations ) then
          if (ityp(nb)==ityp(na)) then
            ft (:) = xau(:,na) - xau(:,nb) - nint( xau(:,na) - xau(:,nb) )
            call checksym (irot, nat, ityp, xau, xau, ft, sym, irt)
!
            if (sym(irot).and.(abs (ft(1)**2+ft(2)**2+ft(3)**2) < 1.d-8)) then
              write(6,*)'Overlapping atoms!'
              stop
            endif
!
            if (sym(irot).and.ttag) then
              ttag=.false.
        !     fractional_translations = .false.
            write(6,'(5x,"Found symmetry operation: I + (",&
     &             3f8.4, ")",/,5x,"This is a supercell!")') ft
            write(6,'(5x,"Nothing to worry about. Just noticing.")')
            endif
          endif
        endif
      enddo
!
      do irot = 1, nrot
        do na = 1, nat
        ! rau = rotated atom coordinates
          rau (:, na) = s (1,:, irot) * xau (1, na) + &
                        s (2,:, irot) * xau (2, na) + &
                        s (3,:, irot) * xau (3, na)
        enddo
        !
        ! first attempt: no fractional translation
        !
        ft (:) = 0.d0
        call checksym (irot,nat,ityp,xau,rau,ft,sym,irt)
        !
        if (.not.sym (irot).and.fractional_translations) then
          nb = 1
          do na = 1, nat
            if (ityp(nb)==ityp(na)) then
            !
            ! second attempt: check all possible fractional translations
            !
              ft (:) = rau(:,na)-xau(:,nb)-nint(rau(:,na)-xau(:,nb))
              call checksym (irot,nat,ityp,xau,rau,ft,sym,irt)
              if (sym(irot)) goto 100
            endif
          enddo
        endif
100     continue
      enddo
      deallocate (rau)
      return
      end
!-----------------------------------------------------------------------
      subroutine checksym (ir,nat,ityp,xau,rau,ft,sym,irt)
!-----------------------------------------------------------------------
!
!   This routine receives as input all the atomic positions xau,
!   and the rotated rau by the symmetry operation ir. It sets to true
!   sym(ir) if for each atom na, it is possible to find an atom nb
!   which is of the same type of na, and coincide with it after the
!   symmetry operation. Fractional translations are allowed.
!
!   Revised layout 1 may 1995 by A. Dal Corso
!
      implicit none
      integer nat, irt (48, nat), ir
      character*2 ityp (nat)
      real*8 xau (3, nat), rau (3, nat), ft (3)
      logical sym (48)
      integer na, nb
      logical eqvect
!
      do na = 1, nat
        do nb = 1, nat
          sym(ir)=ityp(na).eq.ityp(nb).and.eqvect(rau(1,na),xau(1,nb),ft)
          if (sym (ir)) then
            ! the rotated atom does coincide with one of the like atoms
            ! keep track of which atom the rotated atom coincides with
            irt (ir, na) = nb
            goto 10
          endif
        enddo
        ! the rotated atom does not coincide with any of the like atoms
        ! s(ir) + ft is not a symmetry operation
        return
10      continue
      enddo
!
! s(ir) + ft is a symmetry operation
!
      return
      end
!-----------------------------------------------------------------------
      logical function eqvect (x, y, f)
!-----------------------------------------------------------------------
!
!   This function test if the difference x-y-f is an integer.
!   x, y = 3d vectors in crystal axis, f = fractionary translation
!
      use constants
      use controls , only : numthr
      implicit none
      real*8 x(3),y(3),f(3)
      real*8 accep

      accep=numthr
      eqvect = abs(x(1)-y(1)-f(1)-nint(x(1)-y(1)-f(1))) < accep .and. &
     &     abs( x(2)-y(2)-f(2) - nint(x(2)-y(2)-f(2)) ) < accep .and. &
     &     abs( x(3)-y(3)-f(3) - nint(x(3)-y(3)-f(3)) ) < accep
!
      return
      end
!-----------------------------------------------------------------------
      subroutine divided_differences
!-----------------------------------------------------------------------
      use constants
      use controls
      use crystal
      use tensors
      implicit none
      integer i,j,k,l,m,n,nx,idist,ou,necs
      integer nm,np,ii(6),jdef(ndef)
      integer, allocatable :: im(:),ip(:)
      character tag*4
!
!
!
      tc2=.false. ; tc3=.false. ; tc4=.false. ; tc5=.false.
      tc6=.false. ; tc7=.false. ; tc8=.false.
!
      nx=nvec
      if (tout) then
        ou=11
        open(ou,file=trim(dname)//'indexing')
        write(ou,87) nx,pmst,iecs
        do i=1,ndef
          write(ou,88) i,(idef(k,i),k=1,6)
        enddo
87      format(1x,i7,1x,f11.6,1x,i3)
88      format(1x,i7,2x,6(1x,i3),3x,i7,1x,i3,1x,l1)
        do i=1,nx
          k=isym(2,i)
          if (tout) write(ou,88) i,(grid(j,i),j=1,6),isym(1,i),k,tcry(k)
        enddo
      endif
!
!
      allocate(im(nx),ip(nx),imap(ndef,nx))
      nm=nx
      do i=1,nx
        im(i)=i
        ip(i)=0
      enddo
      do i=1,nx*ndef
        imap(i,1)=0
      enddo
!
      necs=1
1     necs=necs+1
      np=0
      do i=1,nm 
        m=0
        n=im(i)
        do l=1,ndef
          do j=1,6
            ii(j)=grid(j,n)+idef(j,l)
          enddo
          do j=1,nm
            k=idist(ii,grid(1,im(j)),6)
            if (k.eq.0) then
              m=m+1
              jdef(l)=j
              goto 15
            endif
          enddo
15        continue
        enddo
        if (m.eq.ndef) then
          np=np+1
          ip(np)=n
          do k=1,ndef
            imap(k,np)=jdef(k)
          enddo
        endif
      enddo
!
      if (tout) then
        write(ou,89) np,necs
        do i=1,np
          write(ou,89) i,(imap(k,i),k=1,ndef)
        enddo
89      format(13(1x,i7))
      endif
!
      if (tpk2ecs) call hoecsdd(necs,np)
!
      do i=1,np
        im(i)=ip(i)
        ip(i)=0
        do j=1,ndef
          imap(j,i)=0
        enddo
      enddo
      nm=np
      if (necs.lt.iecs) goto 1
!
      if (tout) close(ou)
      deallocate(im,ip,imap)
!
!
      if (.not.tpk2ecs) return
!
      tag='ecsn'
      call print_hoecs(tag)
!
      open(oall,file=trim(dname)//'hoecs.'//tag)
      write(oall,*)'CELL_PARAMETER (bohr)'
      do i=1,3
        write(oall,*)(box(k,i)/bohr,k=1,3)
      enddo
      write(oall,*)'CAUCHY (GPa)'
      do i=1,3
        write(oall,*)(cau(k,i,iref),k=1,3)
      enddo
      if (.not.tc2) return
      write(oall,*)'2ECS (GPa)'
      write(oall,*)(c2r(i,1),i=1,6*6)
!
    ! necs=2
    ! call figure_of_merit(necs)
!
      if (.not.tc3) return
      write(oall,*)'3ECS (GPa)'
      write(oall,*)(c3r(i,1,1),i=1,6*6*6)
!
    ! necs=3
    ! call figure_of_merit(necs)
!
      if (.not.tc4) return
      write(oall,*)'4ECS (GPa)'
      write(oall,*)(c4r(i,1,1,1),i=1,6*6*6*6)
!
    ! necs=4
    ! call figure_of_merit(necs)
!
      if (.not.tc5) return
      write(oall,*)'5ECS (GPa)'
      write(oall,*)(c5r(i,1,1,1,1),i=1,6*6*6*6*6)
!
    ! necs=5
    ! call figure_of_merit(necs)
!
      if (.not.tc6) return
      write(oall,*)'6ECS (GPa)'
      write(oall,*)(c6r(i,1,1,1,1,1),i=1,6*6*6*6*6*6)
!
    ! necs=6
    ! call figure_of_merit(necs)
!
      if (.not.tc7) return
      write(oall,*)'7ECS (GPa)'
      write(oall,*)(c7r(i,1,1,1,1,1,1),i=1,6*6*6*6*6*6*6)
!
    ! necs=7
    ! call figure_of_merit(necs)
!
      if (.not.tc8) return
      write(oall,*)'8ECS (GPa)'
      write(oall,*)(c8r(i,1,1,1,1,1,1,1),i=1,6*6*6*6*6*6*6*6)
!
    ! necs=8
    ! call figure_of_merit(necs)
!
      return
      end
!-----------------------------------------------------------------------
      subroutine figure_of_merit(necs)
!-----------------------------------------------------------------------
      use constants
      use controls
      use tensors
      implicit none
      integer i,j,k,l,necs
      real*8 vo(6),pc(6),pr(6),a,b
      character tag*255,num*1
!
      if (.not.tout) return
!
      write(num,'(i1.1)') necs
      tag=trim(dname)//'fom.'//num
      
      open(ocom,file=trim(tag))
      do i=1,nvec
!
        b=zero
        do l=1,6
          vo(l)=dble(grid(l,i))*pmst
          b=b+vo(l)**2
        enddo
        b=sqrt(b)
!
        call e2pk(vo,pc)
        call p2v(pio(1,1,i),pr)
!
        a=zero
        do j=1,6
          a=a+(pc(j)-pr(j))**2
        enddo
        a=sqrt(a)
!
        write(ocom,291) a,b,(grid(k,i),k=1,6)
291     format(2(1x,f12.4),1x,6(1x,i2))
      enddo
      close(ocom)
      return
      end
!-----------------------------------------------------------------------
      subroutine e2pk(e,p)
!-----------------------------------------------------------------------
      use constants
      use tensors
      implicit none
      real*8 , parameter :: three=3.d0, four=4.d0, five=5.d0, six=6.d0
      real*8 , parameter :: seven=7.d0
      integer i,j,k,l,m,n,q,r,s,t
      real*8 e(6),p(6),f5,f6,f7,f8
!
      f5=one/two/three/four
      f6=f5/five
      f7=f6/six
      f8=f7/seven
!
      n=6
      call p2v(pio(1,1,iref),p)
      do i=1,n
        p(i)=-p(i)
      enddo
!
      if (.not.tc2) goto 222
      do i=1,n ; do j=1,n
        p(i)=p(i)+c2r(i,j)*e(j)
      enddo ; enddo
!
      if (.not.tc3) goto 222
      do i=1,n ; do j=1,n ; do k=1,n
        p(i)=p(i)+c3r(i,j,k)*e(j)*e(k)/two
      enddo ; enddo ; enddo
!
      if (.not.tc4) goto 222
      do i=1,n ; do j=1,n ; do k=1,n ; do l=1,n
        p(i)=p(i)+c4r(i,j,k,l)*e(j)*e(k)*e(l)/two/three
      enddo ; enddo ; enddo ; enddo
!
      if (.not.tc5) goto 222
      do i=1,n ; do j=1,n ; do k=1,n ; do l=1,n ; do m=1,n
        p(i)=p(i)+c5r(i,j,k,l,m)*e(j)*e(k)*e(l)*e(m)*f5
      enddo ; enddo ; enddo ; enddo ; enddo
!
      if (.not.tc6) goto 222
      do i=1,n ; do j=1,n ; do k=1,n ; do l=1,n ; do m=1,n ; do q=1,n
        p(i)=p(i)+c6r(i,j,k,l,m,q)*e(j)*e(k)*e(l)*e(m)*e(q)*f6
      enddo ; enddo ; enddo ; enddo ; enddo ; enddo
!
      if (.not.tc7) goto 222
      do i=1,n ; do j=1,n ; do k=1,n ; do l=1,n ; do m=1,n
      do q=1,n ; do r=1,n
        p(i)=p(i)+c7r(i,j,k,l,m,q,r)*e(j)*e(k)*e(l)*e(m)*e(q)*e(r)*f7
      enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo
!
      if (.not.tc8) goto 222
      do i=1,n ; do j=1,n ; do k=1,n ; do l=1,n ; do m=1,n
      do q=1,n ; do r=1,n ; do t=1,n
        p(i)=p(i)+c8r(i,j,k,l,m,q,r,t)*e(j)*e(k)*e(l)*e(m)*e(q)*e(r)*e(t)*f8
      enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo
!
222   continue
      do i=1,6
        p(i)=-p(i)
      enddo
      return
      end
!-----------------------------------------------------------------------
      SUBROUTINE set_sym_bl
!---------------------------------------------------------------------
     !! Provides symmetry operations for all bravais lattices.  
     !! Tests the 24 proper rotations for the cubic lattice first, then
     !! the 8 rotations specific for the hexagonal axis (special axis c), 
     !! then inversion is added.
     !
      use constants
      use crystal
      use controls , only : numthr
      IMPLICIT NONE
     !
     ! sin3 = sin(pi/3), cos3 = cos(pi/3), msin3 = -sin(pi/3), mcos3 = -cos(pi/3)
     !
      REAL*8, PARAMETER :: sin3 = 0.866025403784438597d0,  cos3 =  0.5d0, &
                           msin3 =-0.866025403784438597d0, mcos3 = -0.5d0
     !
      integer i,k,l,nrot
      real*8 at(3,3)
      REAL*8 :: s0(3,3,32), overlap(3,3), rat(3), rot(3,3), value
     ! s0: the s matrices in cartesian axis
     ! overlap: inverse overlap matrix between direct lattice
     ! rat: the rotated of a direct vector ( cartesian )
     ! rot: the rotated of a direct vector ( crystal axis )
     ! value: component of the s matrix in axis basis
     INTEGER :: jpol, kpol, mpol, irot
     ! counters over the polarizations and the rotations
     !
     CHARACTER(LEN=45) :: s0name(64)
     ! full name of the rotational part of each symmetry operation
     !
     DATA s0/ 1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
             -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
             -1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
              1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
              0.d0,  1.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
              0.d0, -1.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
              0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
              0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
              0.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  0.d0, &
              0.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0,  0.d0, &
              0.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0,  0.d0, &
              0.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  0.d0, &
             -1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0, &
             -1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0, &
              1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0, &
              1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0, &
              0.d0,  0.d0,  1.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, &
              0.d0,  0.d0, -1.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, &
              0.d0,  0.d0, -1.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, &
              0.d0,  0.d0,  1.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, &
              0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  1.d0,  0.d0,  0.d0, &
              0.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  1.d0,  0.d0,  0.d0, &
              0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0, -1.d0,  0.d0,  0.d0, &
              0.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0, -1.d0,  0.d0,  0.d0, &
              cos3,  sin3, 0.d0, msin3,  cos3, 0.d0, 0.d0, 0.d0,  1.d0, &
              cos3, msin3, 0.d0,  sin3,  cos3, 0.d0, 0.d0, 0.d0,  1.d0, &
             mcos3,  sin3, 0.d0, msin3, mcos3, 0.d0, 0.d0, 0.d0,  1.d0, &
             mcos3, msin3, 0.d0,  sin3, mcos3, 0.d0, 0.d0, 0.d0,  1.d0, &
              cos3, msin3, 0.d0, msin3, mcos3, 0.d0, 0.d0, 0.d0, -1.d0, &
              cos3,  sin3, 0.d0,  sin3, mcos3, 0.d0, 0.d0, 0.d0, -1.d0, &
             mcos3, msin3, 0.d0, msin3,  cos3, 0.d0, 0.d0, 0.d0, -1.d0, &
             mcos3,  sin3, 0.d0,  sin3,  cos3, 0.d0, 0.d0, 0.d0, -1.d0 /
     !
     DATA s0name/  'identity                                     ',&
                   '180 deg rotation - cart. axis [0,0,1]        ',&
                   '180 deg rotation - cart. axis [0,1,0]        ',&
                   '180 deg rotation - cart. axis [1,0,0]        ',&
                   '180 deg rotation - cart. axis [1,1,0]        ',&
                   '180 deg rotation - cart. axis [1,-1,0]       ',&
                   ' 90 deg rotation - cart. axis [0,0,-1]       ',&
                   ' 90 deg rotation - cart. axis [0,0,1]        ',&
                   '180 deg rotation - cart. axis [1,0,1]        ',&
                   '180 deg rotation - cart. axis [-1,0,1]       ',&
                   ' 90 deg rotation - cart. axis [0,1,0]        ',&
                   ' 90 deg rotation - cart. axis [0,-1,0]       ',&
                   '180 deg rotation - cart. axis [0,1,1]        ',&
                   '180 deg rotation - cart. axis [0,1,-1]       ',&
                   ' 90 deg rotation - cart. axis [-1,0,0]       ',&
                   ' 90 deg rotation - cart. axis [1,0,0]        ',&
                   '120 deg rotation - cart. axis [-1,-1,-1]     ',&
                   '120 deg rotation - cart. axis [-1,1,1]       ',&
                   '120 deg rotation - cart. axis [1,1,-1]       ',&
                   '120 deg rotation - cart. axis [1,-1,1]       ',&
                   '120 deg rotation - cart. axis [1,1,1]        ',&
                   '120 deg rotation - cart. axis [-1,1,-1]      ',&
                   '120 deg rotation - cart. axis [1,-1,-1]      ',&
                   '120 deg rotation - cart. axis [-1,-1,1]      ',&
                   ' 60 deg rotation - cryst. axis [0,0,1]       ',&
                   ' 60 deg rotation - cryst. axis [0,0,-1]      ',&
                   '120 deg rotation - cryst. axis [0,0,1]       ',&
                   '120 deg rotation - cryst. axis [0,0,-1]      ',&
                   '180 deg rotation - cryst. axis [1,-1,0]      ',&
                   '180 deg rotation - cryst. axis [2,1,0]       ',&
                   '180 deg rotation - cryst. axis [0,1,0]       ',&
                   '180 deg rotation - cryst. axis [1,1,0]       ',&
                   'inversion                                    ',&
                   'inv. 180 deg rotation - cart. axis [0,0,1]   ',&
                   'inv. 180 deg rotation - cart. axis [0,1,0]   ',&
                   'inv. 180 deg rotation - cart. axis [1,0,0]   ',&
                   'inv. 180 deg rotation - cart. axis [1,1,0]   ',&
                   'inv. 180 deg rotation - cart. axis [1,-1,0]  ',&
                   'inv.  90 deg rotation - cart. axis [0,0,-1]  ',&
                   'inv.  90 deg rotation - cart. axis [0,0,1]   ',&
                   'inv. 180 deg rotation - cart. axis [1,0,1]   ',&
                   'inv. 180 deg rotation - cart. axis [-1,0,1]  ',&
                   'inv.  90 deg rotation - cart. axis [0,1,0]   ',&
                   'inv.  90 deg rotation - cart. axis [0,-1,0]  ',&
                   'inv. 180 deg rotation - cart. axis [0,1,1]   ',&
                   'inv. 180 deg rotation - cart. axis [0,1,-1]  ',&
                   'inv.  90 deg rotation - cart. axis [-1,0,0]  ',&
                   'inv.  90 deg rotation - cart. axis [1,0,0]   ',&
                   'inv. 120 deg rotation - cart. axis [-1,-1,-1]',&
                   'inv. 120 deg rotation - cart. axis [-1,1,1]  ',&
                   'inv. 120 deg rotation - cart. axis [1,1,-1]  ',&
                   'inv. 120 deg rotation - cart. axis [1,-1,1]  ',&
                   'inv. 120 deg rotation - cart. axis [1,1,1]   ',&
                   'inv. 120 deg rotation - cart. axis [-1,1,-1] ',&
                   'inv. 120 deg rotation - cart. axis [1,-1,-1] ',&
                   'inv. 120 deg rotation - cart. axis [-1,-1,1] ',&
                   'inv.  60 deg rotation - cryst. axis [0,0,1]  ',&
                   'inv.  60 deg rotation - cryst. axis [0,0,-1] ',&
                   'inv. 120 deg rotation - cryst. axis [0,0,1]  ',&
                   'inv. 120 deg rotation - cryst. axis [0,0,-1] ',&
                   'inv. 180 deg rotation - cryst. axis [1,-1,0] ',&
                   'inv. 180 deg rotation - cryst. axis [2,1,0]  ',&
                   'inv. 180 deg rotation - cryst. axis [0,1,0]  ',&
                   'inv. 180 deg rotation - cryst. axis [1,1,0]  ' /
     !
     ! ... compute the overlap matrix for crystal axis
     !
     at(:,:)=box(:,:)
     !
     DO jpol = 1, 3
        DO kpol = 1, 3
           rot(kpol,jpol) = at(1,kpol)*at(1,jpol) + &
                            at(2,kpol)*at(2,jpol) + &
                            at(3,kpol)*at(3,jpol)
        ENDDO
     ENDDO
     !
     ! ... then its inverse (rot is used as work space)
     CALL invmat( 3, rot, overlap, value )
     !
     nrot = 1
     !
     DO irot = 1, 32
        !
        ! ... for each possible symmetry
        DO jpol = 1, 3
           DO mpol = 1, 3
              !
              ! ... compute, in cartesian coordinates the rotated vector
              rat(mpol) = s0(mpol,1,irot)*at(1,jpol) + &
                          s0(mpol,2,irot)*at(2,jpol) + &
                          s0(mpol,3,irot)*at(3,jpol)
           ENDDO

           DO kpol = 1, 3
              !
              ! ... the rotated vector is projected on the direct lattice
              rot(kpol,jpol) = at(1,kpol)*rat(1) + &
                               at(2,kpol)*rat(2) + &
                               at(3,kpol)*rat(3)
           ENDDO
        ENDDO
        !
        ! ... and the inverse of the overlap matrix is applied
        DO jpol = 1,3
           DO kpol = 1,3
              value = overlap(jpol,1)*rot(1,kpol) + &
                      overlap(jpol,2)*rot(2,kpol) + &
                      overlap(jpol,3)*rot(3,kpol)
              !
              IF ( ABS(DBLE(NINT(value))-value) > numthr ) THEN
                 !
                 ! ... if a noninteger is obtained, this implies that this operation
                 ! is not a symmetry operation for the given lattice
                 !
                 GOTO 10
              ENDIF
              !
              sabc(kpol,jpol,nrot) = nint(value)
           ENDDO
        ENDDO
        !
        sname(nrot) = s0name(irot)
        imat(nrot) = irot
        nrot = nrot+1
        !
   10   CONTINUE
        !
     ENDDO
     !
     nrot = nrot-1
     !
     IF ( nrot /= 1 .AND. nrot /= 2 .AND. nrot /= 4 .AND. nrot /= 6 .AND. &
          nrot /= 8 .AND. nrot /=12 .AND. nrot /=24 ) THEN
          WRITE (6,*)'Bravais lattice has wrong number of symmetries', nrot
         nrot = 1
     ENDIF
     !
     ! ... set the inversion symmetry (Bravais lattices have always inversion symmetry)
     DO irot = 1, nrot
        sname(irot+nrot) = s0name(imat(irot)+32)
        imat(irot+nrot) = -imat(irot)
        DO kpol = 1, 3
           DO jpol = 1, 3
              sabc(kpol,jpol,irot+nrot) = -sabc(kpol,jpol,irot)
           ENDDO
        ENDDO
        i=imat(irot)
        sym(:,:,irot)     = s0(:,:,i)
        sym(:,:,irot+nrot)=-s0(:,:,i)
     ENDDO
     !
     nrot = 2*nrot
     nsym=nrot
     !
     RETURN
     END
!-----------------------------------------------------------------------
      subroutine read_arguments
!-----------------------------------------------------------------------
      use constants
      use controls
      use tensors , only : iref
      implicit none
      real*8 d(3)
      character arg*255,argv*255,tag*6
      integer i,n,istat,ii(6),k
      logical ttmp
!
!
      pmst=0.01d0
      iecs=2
      fname='refstate'
      dname='./'
      tonlyi=.true.
      numthr = 1.d-4
      tout=.false.
      tdeform=.false.
      tpk2ecs=.false.
      tperm=.true.
      thomo=.false.
      tone=.false.
      ione=0
!
      n=command_argument_count()
      i=0
10    i=i+1
      if (i.gt.n) goto 1000
      call get_command_argument(i, arg)
!
      select case (trim(arg))
!
            case ('str2pk')
                if (i.ne.1) goto 900
                tdeform=.true.

            case ('pk2ecs')
                if (i.ne.1) goto 900
                tpk2ecs=.true.

            case ('-h','--help')
                if (i.ne.1) then
                  write(6,*)' -h : ignored ..'
                else
                  call help_routine
                  stop
                endif

! strain parameter
            case ('-s','--strpar')
                i=i+1
                call get_command_argument(i,argv)
                read(argv,*,err=790) pmst
                if (tpk2ecs) write(6,*)' -s: ignored ..'

! highest order of the elastic constants
            case ('-n','--norder')
                i=i+1
                call get_command_argument(i,argv)
                read(argv,*) iecs
                if (tpk2ecs) write(6,*)' -n: ignored ..'

            case ('-r','--refstate')
                i=i+1
                call get_command_argument(i,argv)
                fname=trim(argv)
                if (tpk2ecs) write(6,*)' -r: ignored ..'
      
            case ('-d','--dir')
                i=i+1
                call get_command_argument(i,argv)
                dname=trim(argv)
      
            case ('-nosym')
                tonlyi=.false.
                if (tpk2ecs) write(6,*)' -nosym: ignored ...'

            case ('-noper')
                tperm=.false.
                if (tdeform) write(6,*)' -noper: ignored ...'

            case ('-matsym')
                thomo=.true.
                if (tdeform) then
                  write(6,*)' -matsym: ignored ...'
                  thomo=.false.
                endif

            case ('-t','--thr')
                i=i+1
                call get_command_argument(i,argv)
                read(argv,*) numthr

            case ('-v','--verbose')
                tout=.true.
                oname='debug.out'

            case ('-p','--pencil')
                i=i+1
                call get_command_argument(i,argv)
                read(argv,*,err=800) ione 
                tone=.true.

            case default
                write(6,*) ' Argument not recognized'

      end select
      goto 10
1000  continue
!
!
!
      if (.not.(tdeform.or.tpk2ecs)) goto 900
!
      if (trim(dname).ne.'./') then
        k=len_trim(dname)
        if (dname(k:k).ne.'/') dname=trim(dname)//'/'
      endif
      inquire(file=trim(dname), exist=ttmp)
!
!
      if (tpk2ecs) then
!
        if (.not.ttmp) then
           write(6,*)'Directory not found. Aborting.'
           stop
        endif
!
        iref=1                     ! Default reference state index
        write(tag,'(i6.6)') iref
!
        fname=trim(dname)//'ipos.'//tag
        inquire(file=trim(fname), exist=ttmp)
        if (.not.ttmp) then
          write(6,*)'Directory does not contain "ipos.000001". Aborting'
          stop
        endif
!
        arg=trim(dname)//'istr.'//tag
        inquire(file=trim(arg), exist=ttmp)
        if (.not.ttmp) then
          write(6,*)'Directory does not contain "istr.000001". Aborting'
          stop
        endif
!
        open(oist,file=trim(arg),status='old',iostat=i)
        if (i.ne.0) then
          write(6,*)'Error opening file "istr.000001". Aborting.'
          stop
        endif
        read(oist,*,err=950,end=950) iecs,i,(ii(i),i=1,6),pmst,i,i
        do i=1,6
          read(oist,*,err=950,end=950) (d(k),k=1,3)
        enddo
        close(oist)
!
      endif
!
!
      if (tdeform) then
!
        if (.not.ttmp) then
          call system('mkdir '//trim(dname),istat)
          if (istat.ne.0) then
            write(6,*)'Cannot create directory. Aborting.'
            stop
          endif
        endif
!
        if (ttmp.and..not.(trim(dname).eq.'./')) then
          write(6,*)'Directory exits. Aborting.'
          stop
        endif
!
        inquire(file=trim(fname), exist=ttmp)
        if (.not.ttmp) then
          write(6,*)'Filename of the reference state "',trim(fname),'" not existing. Aborting'
          write(6,*)
          write(6,*)'Sample input file (pwscf format):'
          write(6,*)
          write(6,*)'CELL_PARAMETERS (bohr)'
          write(6,*)'0.000000000  5.131553337  5.131553337'
          write(6,*)'5.131553337  0.000000000  5.131553337'
          write(6,*)'5.131553337  5.131553337  0.000000000'
          write(6,*)'ATOMIC_POSITIONS (crystal)'
          write(6,*)'Si  0.000  0.000  0.000'
          write(6,*)'Si  0.250  0.250  0.250'
          stop
        endif
!
      endif
!
      iecs=max(iecs,2)
      if (iecs.ge.8) then
        write(6,*) 'Unless the stress is calculated with (super) high precision,'
        write(6,*) 'this method yields elastic constants of order larger than 7 '
        write(6,*) 'that are, most likely, inaccurate.'
      endif
!
      if (tout) then
        open(oout,file=trim(dname)//trim(oname))
        write(oout,*)'----------------------------------------------'
        write(oout,'(a9,1x,l1)') 'tdeform: ', tdeform
        write(oout,'(a9,1x,l1)') 'tpk2ecs: ', tpk2ecs
        write(oout,'(a9,1x,a)') 'dname  : ', trim(dname)
        write(oout,'(a9,1x,a)') 'fname  : ', trim(fname)
        write(oout,'(a9,1x,i1)') 'iecs   : ', iecs
        write(oout,'(a9,1x,e10.4)') 'pmst   : ', pmst
        write(oout,'(a9,1x,l1)') 'tonlyi : ', tonlyi
        write(oout,'(a9,1x,l1)') 'tperm  : ', tperm
        write(oout,'(a9,1x,l1)') 'thomo  : ', thomo
        write(oout,'(a9,1x,e10.4)') 'numthr : ', numthr
        write(oout,*)'----------------------------------------------'
      endif
!
      return
790   write(6,*) 'Error reading the strain parameter. Aborting!'
      stop
800   write(6,*) 'Error reading the Voigt strain vector. Aborting!'
      stop
900   write(6,*)' First argument must be either "str2pk" or "pk2ecs".'
      write(6,*)
      write(6,*)' Examples:'                       
      write(6,*)' ./hoecs.x str2pk -n 3 -r rstate.dat'                       
      write(6,*)' ./hoecs.x pk2ecs -d ./files'                       
      stop
950   write(6,*) 'Error reading "istr.000001". Aborting!'
      stop
      end
!-----------------------------------------------------------------------
      subroutine help_routine
!-----------------------------------------------------------------------
      implicit none

      write(6,*)
      write(6,*)'This module applies finite strains to a reference state'
      write(6,*)'state (mode "str2pk"), and  after  computing the stress'
      write(6,*)'tensors resulting from the deformations, it can be used'
      write(6,*)'(mode "pk2ecs") to  calculate  the elastic constants up'
      write(6,*)'to a  given  order.  The  module accepts arguments. The' 
      write(6,*)'first  argument  specifying the  "mode" of operation is'
      write(6,*)'mandatory. The rest of the arguments are  optional (see'
      write(6,*)'below their default value in parantesis).'
      write(6,*)
      write(6,*)'Mode "str2pk" (1st argument):'
      write(6,*)
      write(6,*)' -r, --refstate [<str>] : Reference state ("refstate")'
      write(6,*)' -s, --strpar   [<real>]: Strain parameter (0.01)'
      write(6,*)' -n, --norder   [<int>] : Maximum EC order (2)'
      write(6,*)' -nosym                 : Turn off symmetry (ON)'
      write(6,*)
      write(6,*)'Mode "pk2ecs" (1st argument):' 
      write(6,*)
      write(6,*)' -noper                 : Turn off Voigt symmetr. (ON)'
      write(6,*)' -matsym                : Turn on tensor symmetr. (ON)'
      write(6,*)
      write(6,*)'Common arguments:'
      write(6,*)
      write(6,*)' -h, --help             : This message'
      write(6,*)' -d, --dir      [<str>] : Working directory ("./")'
      write(6,*)' -v, --verbose          : Extensive output (OFF)'
      write(6,*)' -t, --thr      [<real>]: Tolerance (1.e-4)'
      write(6,*)
      write(6,*)
      write(6,*)'Example 1: ./hoecs.x str2pk -r rfile -s 0.005 -n 5'
      write(6,*)'Example 2: ./hoecs.x pk2ecs -d ./files -matsym'
      write(6,*)
!
      return
      end
!-----------------------------------------------------------------------
      function idist(i,j,nx)
!-----------------------------------------------------------------------
      implicit none
      integer i(nx),j(nx),jj,n,nx,idist
      jj=0
      do n=1,nx
        jj=jj+(i(n)-j(n))**2
      enddo
      idist=jj
      return
      end
!-----------------------------------------------------------------------
      subroutine cauchy2piola(vo,pc,pp)
!-----------------------------------------------------------------------
      implicit none
      real*8 , parameter :: zero=0.d0,one=1.d0,two=2.d0
      integer i,j,k,ou,ipov(3)
!
      logical tout
      real*8 str(3,3),vo(6),pc(3,3),pp(3,3),cho(3,3)
      real*8 mat(3,3),tmp(3,3),detf
      real*8 chot(3,3),chom(3,3),chotm(3,3)
      integer , parameter :: lw=20
      real*8 w(3),wrk(lw)
!
!
      tout=.false.
      ou=100
!
12    format(6(1x,f12.6))
      if (tout) then
        write(ou,*)'#'
        write(ou,*)'Voigt Lagrangian strain'
        write(ou,12)(vo(i),i=1,6)
        write(ou,*)'Cauchy stress'
        do i=1,3
          write(ou,12) (pc(k,i),k=1,3)
        enddo
      endif
!
      str(1,1)=vo(1)
      str(2,2)=vo(2)
      str(3,3)=vo(3)
      str(1,2)=vo(6)/two
      str(2,1)=str(1,2)
      str(1,3)=vo(5)/two
      str(3,1)=str(1,3)
      str(2,3)=vo(4)/two
      str(3,2)=str(2,3)
!
      if (tout) then
        write(ou,*)'Strain matrix - input'
        do i=1,3
          write(ou,12) (str(k,i),k=1,3)
        enddo
      endif
!
      do i=1,9
        cho(i,1)=two*str(i,1)
      enddo
      do i=1,3
        cho(i,i)=cho(i,i)+one
      enddo
!
      call chole(3, cho, mat)
      do i=1,3
      do j=1,3
        cho(i,j)=mat(j,i)
      enddo
      enddo

      call dgesvd('A','A',3,3,cho,3,w,tmp,3,chot,3,wrk,lw,i)
      if (i.ne.0) stop 6780
      do i=1,9 
        tmp(i,1)=zero
        mat(i,1)=zero
        chotm(i,1)=zero
      enddo
      do i=1,3
        mat(i,i)=w(i)
      enddo
      call dgemm('t','n',3,3,3,one,chot,3,mat,3,zero,tmp,3)
      call dgemm('n','n',3,3,3,one,tmp,3,chot,3,zero,chotm,3)
!
      do i=1,3*3
        cho(i,1)=chotm(i,1)
      enddo
      call dgemm('t','n',3,3,3,one,cho,3,cho,3,zero,tmp,3)
      do i=1,3
        tmp(i,i)=tmp(i,i)-one
      enddo
      do i=1,9
        tmp(i,1)=tmp(i,1)/two
      enddo
!
      if (tout) then
        write(ou,*)'Strain matrix - recalculated'
        do i=1,3
          write(ou,12) (tmp(k,i),k=1,3)
        enddo
      endif
!
      do i=1,3*3
        mat(i,1)=cho(i,1)
      enddo
      call dgetrf( 3, 3, mat, 3, ipov ,i)
      detf=one
      do i=1,3
        detf=detf*mat(i,i)
      enddo
      detf=abs(detf)
!
      do i=1,3*3
        chom(i,1)=zero
        chotm(i,1)=zero
      enddo
      call invert3x3(cho,chom)
      do i=1,3
        do j=1,3
          tmp(i,j)=cho(j,i)
        enddo
      enddo
      call invert3x3(tmp,chotm)
!
      call dgemm('n','n',3,3,3,one,cho,3,chom,3,zero,mat,3)
!
      if (tout) then
        write(ou,*)'Check unity'
        do i=1,3
          write(ou,12) (mat(k,i),k=1,3)
        enddo
      endif
!
      call dgemm('n','n',3,3,3,one,tmp,3,chotm,3,zero,mat,3)
!
      if (tout) then
        write(ou,*)'Check unity - again'
        do i=1,3
          write(ou,12) (mat(k,i),k=1,3)
        enddo
      endif
!
      do i=1,3*3
        tmp(i,1)=zero
      enddo
      call dgemm('n','n',3,3,3,one,pc,3,chotm,3,zero,tmp,3)
      do i=1,3*3
        pp(i,1)=zero
      enddo
      call dgemm('n','n',3,3,3,detf,chom,3,tmp,3,zero,pp,3)
!
      if (tout) then
        write(ou,*)'Finally .. the 2PK stress'
        do i=1,3
          write(ou,12) (pp(k,i),k=1,3)
        enddo
        write(ou,*)
      endif
!
      return
      end
!-----------------------------------------------------------------------
      subroutine invert3x3(box,boxi)
!-----------------------------------------------------------------------
      implicit none
      real*8 , parameter :: zero = 0.d0 , one = 1.d0
      real*8 box(3,3)
      real*8 a1(3),a2(3),a3(3),omega,boxi(3,3)
      real*8 b1(3),b2(3),b3(3)
      real*8 den,s
      integer i,j,k,l,iperm,ir
!
!
      do k=1,3
        a1(k)=box(k,1)
        a2(k)=box(k,2)
        a3(k)=box(k,3)
      enddo
!
      den=zero
      i=1
      j=2
      k=3
      s=one
1     continue
      do iperm=1,3
         den=den+s*a1(i)*a2(j)*a3(k)
         l=i
         i=j
         j=k
         k=l
      enddo
      i=2
      j=1
      k=3
      s=-s
      if (s.lt.zero) goto 1
      i=1
      j=2
      k=3
      omega=abs(den)
      den=one/omega
!
      do ir=1,3
        b1(ir)=den*(a2(j)*a3(k)-a2(k)*a3(j))
        b2(ir)=den*(a3(j)*a1(k)-a3(k)*a1(j))
        b3(ir)=den*(a1(j)*a2(k)-a1(k)*a2(j))
        l=i
        i=j
        j=k
        k=l
      enddo
!
     !do i=1,3
     !  boxi(i,1)=b1(i)
     !  boxi(i,2)=b2(i)
     !  boxi(i,3)=b3(i)
     !enddo
      do i=1,3
        boxi(1,i)=b1(i)
        boxi(2,i)=b2(i)
        boxi(3,i)=b3(i)
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine chole(n,a,g)
!-----------------------------------------------------------------------
      use constants
      implicit none
      integer n,i,j
      real*8 a(n,n),g(n,n)
!
      g(:,:)=zero
      do j = 1,n
        g(j,j) = sqrt( a(j,j) - dot_product(g(j,1:j-1),g(j,1:j-1)) )
         do i = j+1, n
            g(i,j)=(a(i,j)-dot_product(g(i,1:j-1),g(j,1:j-1)))/g(j,j)
         enddo
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine hoecsdd(necs,np)
!-----------------------------------------------------------------------
      use constants
      use controls
      use crystal
      use tensors
      implicit none
      integer i,j,k,l,m,n,p,q,r,s,ii(6)
      integer necs,np,ip,im
      real*8 c(3,3),a,zeros
!
!
      go to (1002,1003,1004,1005,1006,1007,1008), necs-1
!
!     SOECS
1002  tc2=.true.
      nc2=np
      allocate(c2(6,6,nc2))
      do n=1,nc2
        do i=1,6
          ip=imap(2*i-1,n)
          im=imap(2*i  ,n)
          do l=1,9
            c(l,1)=-(pio(l,1,ip)-pio(l,1,im))/two/pmst
          enddo
          c2(1,i,n)=c(1,1)
          c2(2,i,n)=c(2,2)
          c2(3,i,n)=c(3,3)
          c2(4,i,n)=(c(3,2)+c(2,3))/two
          c2(5,i,n)=(c(3,1)+c(1,3))/two
          c2(6,i,n)=(c(1,2)+c(2,1))/two
        enddo
      enddo
!
      if (tperm) then
        do n=1,nc2
          do i=1,6
            do j=i+1,6
              a=(c2(i,j,n)+c2(j,i,n))/two
              c2(i,j,n)=a
              c2(j,i,n)=a
            enddo
          enddo
        enddo
      else
        n=iref
        do i=1,6
          do j=i+1,6
            a=(c2(i,j,n)+c2(j,i,n))/two
            c2(i,j,n)=a
            c2(j,i,n)=a
          enddo
        enddo
      endif
!
      allocate(c2r(6,6))
      c2r(:,:) = c2(:,:,iref)
      return
!
!     TOECS
!
1003  tc3=.true.
      if (.not.tc2) then
        write(6,*)'Funny business .. ', necs
        stop
      endif
      nc3=np
      allocate(c3(6,6,6,nc3))
      do n=1,nc3
        do i=1,6
          ip=imap(2*i-1,n)
          im=imap(2*i  ,n)
          do l=1,6*6
            c3(l,1,i,n)=(c2(l,1,ip)-c2(l,1,im))/two/pmst
          enddo
        enddo
      enddo
!
      if (tperm) then
        allocate(t3(6,6,6))
        do n=1,nc3
          t3(:,:,:)=.true.
          do i=1,6 ; ii(1)=i ; do j=1,6 ; ii(2)=j ; do k=1,6 ; ii(3)=k
            if (t3(i,j,k)) call permuta(ii(1),3,c3(1,1,1,n),t3(1,1,1))
          enddo ; enddo ; enddo
        enddo
        deallocate(t3)
      else
        allocate(t3(6,6,6))
        n=iref
        t3(:,:,:)=.true.
        do i=1,6 ; ii(1)=i ; do j=1,6 ; ii(2)=j ; do k=1,6 ; ii(3)=k
          if (t3(i,j,k)) call permuta(ii(1),3,c3(1,1,1,n),t3(1,1,1))
        enddo ; enddo ; enddo
        deallocate(t3)
      endif
!
      deallocate(c2)
      allocate(c3r(6,6,6))
      c3r(:,:,:) = c3(:,:,:,iref)
      return
!
!     4OECS
!
1004  tc4=.true.
      nc4=np
      if (.not.tc3) then
        write(6,*)'Funny business .. ', necs
        stop
      endif
!
      allocate(c4(6,6,6,6,nc4))
      do n=1,nc4
        do i=1,6
          ip=imap(2*i-1,n)
          im=imap(2*i  ,n)
          do l=1,6*6*6
            c4(l,1,1,i,n)=(c3(l,1,1,ip)-c3(l,1,1,im))/two/pmst
          enddo
        enddo
      enddo
!
      if (tperm) then
        allocate(t4(6,6,6,6))
        do n=1,nc4
          t4(:,:,:,:)=.true.
          do i=1,6 ; ii(1)=i ; do j=1,6 ; ii(2)=j ; do k=1,6 ; ii(3)=k
          do l=1,6 ; ii(4)=l
            if (t4(i,j,k,l)) call permuta(ii,4,c4(1,1,1,1,n),t4(1,1,1,1))
          enddo ; enddo ; enddo ; enddo
        enddo
        deallocate(t4)
      else
        allocate(t4(6,6,6,6))
        n=iref
        t4(:,:,:,:)=.true.
        do i=1,6 ; ii(1)=i ; do j=1,6 ; ii(2)=j ; do k=1,6 ; ii(3)=k
        do l=1,6 ; ii(4)=l
          if (t4(i,j,k,l)) call permuta(ii,4,c4(1,1,1,1,n),t4(1,1,1,1))
        enddo ; enddo ; enddo ; enddo
        deallocate(t4)
      endif
!
      deallocate(c3)
      allocate(c4r(6,6,6,6))
      c4r(:,:,:,:) = c4(:,:,:,:,iref)
      return
!
!
!     5OECS
!
1005  tc5=.true.
      nc5=np
      if (.not.tc4) then
        write(6,*)'Funny business .. ', necs
        stop
      endif
!
      allocate(c5(6,6,6,6,6,nc5))
      do n=1,nc5
        do i=1,6
          ip=imap(2*i-1,n)
          im=imap(2*i  ,n)
          do l=1,6*6*6*6
            c5(l,1,1,1,i,n)=(c4(l,1,1,1,ip)-c4(l,1,1,1,im))/two/pmst
          enddo
        enddo
      enddo
!
      if (tperm) then
        allocate(t5(6,6,6,6,6))
        do n=1,nc5
          t5(:,:,:,:,:)=.true.
          do i=1,6 ; ii(1)=i ; do j=1,6 ; ii(2)=j ; do k=1,6 ; ii(3)=k 
          do l=1,6 ; ii(4)=l ; do m=1,6 ; ii(5)=m
            if (t5(i,j,k,l,m)) call permuta(ii,5,c5(1,1,1,1,1,n),t5(1,1,1,1,1))
          enddo ; enddo ; enddo ; enddo ; enddo
        enddo
        deallocate(t5)
      else
        allocate(t5(6,6,6,6,6))
        n=iref
        t5(:,:,:,:,:)=.true.
        do i=1,6 ; ii(1)=i ; do j=1,6 ; ii(2)=j ; do k=1,6 ; ii(3)=k
        do l=1,6 ; ii(4)=l ; do m=1,6 ; ii(5)=m
          if (t5(i,j,k,l,m)) call permuta(ii,5,c5(1,1,1,1,1,n),t5(1,1,1,1,1))
        enddo ; enddo ; enddo ; enddo ; enddo
        deallocate(t5)
      endif
!
      deallocate(c4)
      allocate(c5r(6,6,6,6,6))
      c5r(:,:,:,:,:) = c5(:,:,:,:,:,iref)
      return
!
!
!     6OECS
!
1006  tc6=.true.
      nc6=np
      if (.not.tc5) then
        write(6,*)'Funny business .. ', necs
        stop
      endif
!
      allocate(c6(6,6,6,6,6,6,nc6))
      do n=1,nc6
        do i=1,6
          ip=imap(2*i-1,n)
          im=imap(2*i  ,n)
          do l=1,6*6*6*6*6
            c6(l,1,1,1,1,i,n)=(c5(l,1,1,1,1,ip)-c5(l,1,1,1,1,im))/two/pmst
          enddo
        enddo
      enddo
!
      if (tperm) then
        allocate(t6(6,6,6,6,6,6))
        do n=1,nc6
          t6(:,:,:,:,:,:)=.true.
          do i=1,6 ; ii(1)=i ; do j=1,6 ; ii(2)=j ; do k=1,6 ; ii(3)=k
          do l=1,6 ; ii(4)=l ; do m=1,6 ; ii(5)=m ; do p=1,6 ; ii(6)=p
            if (t6(i,j,k,l,m,p)) then 
              call permuta(ii,6,c6(1,1,1,1,1,1,n),t6(1,1,1,1,1,1))
            endif
          enddo ; enddo ; enddo ; enddo ; enddo ; enddo
        enddo
        deallocate(t6)
      else
        allocate(t6(6,6,6,6,6,6))
        n=iref
        t6(:,:,:,:,:,:)=.true.
        do i=1,6 ; ii(1)=i ; do j=1,6 ; ii(2)=j ; do k=1,6 ; ii(3)=k
        do l=1,6 ; ii(4)=l ; do m=1,6 ; ii(5)=m ; do p=1,6 ; ii(6)=p
          if (t6(i,j,k,l,m,p)) then
            call permuta(ii,6,c6(1,1,1,1,1,1,n),t6(1,1,1,1,1,1))
          endif
        enddo ; enddo ; enddo ; enddo ; enddo ; enddo
        deallocate(t6)
      endif
!
      deallocate(c5)
      allocate(c6r(6,6,6,6,6,6))
      c6r(:,:,:,:,:,:) = c6(:,:,:,:,:,:,iref)
      return
!
!
!     7OECS
!
1007  tc7=.true.
      nc7=np
      if (.not.tc6) then
        write(6,*)'Funny business .. ', necs
        stop
      endif
!
      allocate(c7(6,6,6,6,6,6,6,nc7))
      do n=1,nc7
        do i=1,6
          ip=imap(2*i-1,n)
          im=imap(2*i  ,n)
          do l=1,6*6*6*6*6*6
            c7(l,1,1,1,1,1,i,n)=(c6(l,1,1,1,1,1,ip)-c6(l,1,1,1,1,1,im))/two/pmst
          enddo
        enddo
      enddo
!
      if (tperm) then
        allocate(t7(6,6,6,6,6,6,6))
        do n=1,nc7
          t7(:,:,:,:,:,:,:)=.true.
          do i=1,6 ; ii(1)=i ; do j=1,6 ; ii(2)=j ; do k=1,6 ; ii(3)=k
          do l=1,6 ; ii(4)=l ; do m=1,6 ; ii(5)=m ; do p=1,6 ; ii(6)=p
          do q=1,6 ; ii(7)=q
!
          if (t7(i,j,k,l,m,p,q)) then 
            call permuta(ii,7,c7(1,1,1,1,1,1,1,n),t7(1,1,1,1,1,1,1))
          endif
!
        enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo
        enddo
        deallocate(t7)
      else
        allocate(t7(6,6,6,6,6,6,6))
        n=iref
        t7(:,:,:,:,:,:,:)=.true.
        do i=1,6 ; ii(1)=i ; do j=1,6 ; ii(2)=j ; do k=1,6 ; ii(3)=k
        do l=1,6 ; ii(4)=l ; do m=1,6 ; ii(5)=m ; do p=1,6 ; ii(6)=p
        do q=1,6 ; ii(7)=q
!
        if (t7(i,j,k,l,m,p,q)) then
          call permuta(ii,7,c7(1,1,1,1,1,1,1,n),t7(1,1,1,1,1,1,1))
        endif
!
        enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo
        deallocate(t7)
      endif
!
      deallocate(c6)
      allocate(c7r(6,6,6,6,6,6,6))
      c7r(:,:,:,:,:,:,:) = c7(:,:,:,:,:,:,:,iref)
      return
!
!
!     8OECS
!
1008  tc8=.true.
      nc8=np
      if (.not.tc7) then
        write(6,*)'Funny business .. ', necs
        stop
      endif
!
      allocate(c8(6,6,6,6,6,6,6,6,nc8))
      do n=1,nc8
        do i=1,6
          ip=imap(2*i-1,n)
          im=imap(2*i  ,n)
          do l=1,6*6*6*6*6*6*6
      c8(l,1,1,1,1,1,1,i,n)=(c7(l,1,1,1,1,1,1,ip)-c7(l,1,1,1,1,1,1,im))/two/pmst
          enddo
        enddo
      enddo
!
      if (tperm) then
        allocate(t8(6,6,6,6,6,6,6,6))
        do n=1,nc8
          t8(:,:,:,:,:,:,:,:)=.true.
          do i=1,6 ; ii(1)=i ; do j=1,6 ; ii(2)=j ; do k=1,6 ; ii(3)=k
          do l=1,6 ; ii(4)=l ; do m=1,6 ; ii(5)=m ; do p=1,6 ; ii(6)=p
          do q=1,6 ; ii(7)=q ; do r=1,6 ; ii(8)=r
!
            if (t8(i,j,k,l,m,p,q,r)) then 
              call permuta(ii,8,c8(1,1,1,1,1,1,1,1,n),t8(1,1,1,1,1,1,1,1))
            endif
          enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo
        enddo
        deallocate(t8)
      else
        allocate(t8(6,6,6,6,6,6,6,6))
        n=iref
        t8(:,:,:,:,:,:,:,:)=.true.
        do i=1,6 ; ii(1)=i ; do j=1,6 ; ii(2)=j ; do k=1,6 ; ii(3)=k
        do l=1,6 ; ii(4)=l ; do m=1,6 ; ii(5)=m ; do p=1,6 ; ii(6)=p
        do q=1,6 ; ii(7)=q ; do r=1,6 ; ii(8)=r
!
          if (t8(i,j,k,l,m,p,q,r)) then
            call permuta(ii,8,c8(1,1,1,1,1,1,1,1,n),t8(1,1,1,1,1,1,1,1))
          endif
        enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo
        deallocate(t8)
      endif
!
      deallocate(c7)
      allocate(c8r(6,6,6,6,6,6,6,6))
      c8r(:,:,:,:,:,:,:,:) = c8(:,:,:,:,:,:,:,:,iref)
      deallocate(c8)
      return
      end
!-----------------------------------------------------------------------
      subroutine print_hoecs(tag)
!-----------------------------------------------------------------------
      use constants
      use controls
      use crystal
      use tensors
      implicit none
      integer i,j,k,l,m,p,q,r
      character tag*4
!
!
      if (.not.tc2) return
      call bulkmod(tag,pio(1,1,iref),c2r)
!
      open(oecs,file=trim(dname)//tag//'.2')
      do i=1,6 ; do j=1,6
        if (abs(c2r(i,j)).gt.ectoll) write(oecs,112) c2r(i,j),i,j
112     format(1x,f12.2,2x,2i1)
      enddo ; enddo
      close(oecs)
!
      if (.not.tc3) return
      open(oecs,file=trim(dname)//tag//'.3')
      do i=1,6 ; do j=1,6 ; do k=1,6
        if (abs(c3r(i,j,k)).gt.ectoll) write(oecs,113) c3r(i,j,k),i,j,k
113     format(1x,f15.2,2x,3i1)
      enddo ; enddo ; enddo
      close(oecs)
!
      if (.not.tc4) return
      open(oecs,file=trim(dname)//tag//'.4')
      do i=1,6 ; do j=1,6 ; do k=1,6 ; do l=1,6
        if (abs(c4r(i,j,k,l)).gt.ectoll) write(oecs,114) c4r(i,j,k,l),i,j,k,l
114     format(1x,f15.2,2x,4i1)
      enddo ; enddo ; enddo ; enddo
      close(oecs)
!
      if (.not.tc5) return
      open(oecs,file=trim(dname)//tag//'.5')
      do i=1,6 ; do j=1,6 ; do k=1,6 ; do l=1,6 ; do m=1,6
        if (abs(c5r(i,j,k,l,m)).gt.ectoll) write(oecs,115) c5r(i,j,k,l,m),i,j,k,l,m
115     format(1x,f15.2,2x,5i1)
      enddo ; enddo ; enddo ; enddo ; enddo
      close(oecs)
!
      if (.not.tc6) return
      open(oecs,file=trim(dname)//tag//'.6')
      do i=1,6 ; do j=1,6 ; do k=1,6 ; do l=1,6 ; do m=1,6 ; do p=1,6
        if (abs(c6r(i,j,k,l,m,p)).gt.ectoll) write(oecs,116) c6r(i,j,k,l,m,p),i,j,k,l,m,p
116     format(1x,f15.2,2x,6i1)
      enddo ; enddo ; enddo ; enddo ; enddo ; enddo
      close(oecs)
!
      if (.not.tc7) return
      open(oecs,file=trim(dname)//tag//'.7')
      do i=1,6 ; do j=1,6 ; do k=1,6 ; do l=1,6 ; do m=1,6 ; do p=1,6 ; do q=1,6
        if (abs(c7r(i,j,k,l,m,p,q)).gt.ectoll) write(oecs,117) c7r(i,j,k,l,m,p,q),i,j,k,l,m,p,q
117     format(1x,f15.2,2x,7i1)
      enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo
      close(oecs)
!
      if (.not.tc8) return
      open(oecs,file=trim(dname)//tag//'.8')
      do i=1,6 ; do j=1,6 ; do k=1,6 ; do l=1,6 ; do m=1,6 ; do p=1,6 ; do q=1,6 ; do r=1,6
        if (abs(c8r(i,j,k,l,m,p,q,r)).gt.ectoll) then
          write(oecs,118) c8r(i,j,k,l,m,p,q,r),i,j,k,l,m,p,q,r
118       format(1x,f15.2,2x,8i1)
        endif
      enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo
      close(oecs)
!
      return
      end
!-----------------------------------------------------------
      subroutine permuta(ii,nx,c,tc)
!-----------------------------------------------------------
      use permutations
      use controls , only : tout,oout
      implicit none
      real*8 , parameter :: zero =0.d0, one=1.d0, two=2.d0, ten=10.d0
      real*8 , parameter :: toll=1.d-2
      integer i,j,k,n,m,nj,nx
      integer ii(nx),jj(nx),ipoint,idist
      character l*8,t
      real*8 c(6**nx),a,av,b,aa
      logical tc(6**nx)
!
!     ii(1:nx) is the list of input indices
!
      if (nx.gt.8) stop 8
      if (nx.ne.nxp) then
        nxp=nx
!
        l=''
        k=1
        do i=1,nx
          write(l(i:i),'(i1)') i
          k=k*i
        enddo
        np=k
        if (allocated(perm)) deallocate(perm)
        allocate(perm(np))
!
        m=1
        perm(m)=l(1:nx)
1       n=nx
2       t = l(1:1) ! first index ...
        do i=1,n-1
          l(i:i)=l(i+1:i+1)
        enddo
        l(n:n) = t ! ... is now the last one.
!
        do i=1,n-1
          if (l(n:n).le.l(i:i)) then
            m=m+1
            perm(m)=l(1:nx)
            goto 1
          endif
        enddo
        if (n.ne.2) then
          n=n-1
          goto 2
        endif
        if (m.ne.np) stop 89
!
        if (allocated(ind)) deallocate(ind)
        allocate(ind(nx,np))
!
!       write(1000,*) nx,np,' #'
!       do n=1,np
!         write(1000,*) n,trim(perm(n))
!       enddo
      endif
!
!
      nj=0
      do i=1,np
        do j=1,nx
          read(perm(i)(j:j),'(i1)') k
          jj(j)=ii(k)
        enddo
        do j=1,nj
          k=idist(ind(1,j),jj,nx)
          if (k.eq.0) goto 4
        enddo
        nj=nj+1
        ind(:,nj)=jj(:)
!       write(6,11) nj, perm(i), (ind(j,nj),j=1,nx)
4       continue
      enddo
!
11    format(1x,i6,2x,a8,2x,8i1)
12    format(2(f17.2,2x),8i1)
      a=zero
      do i=1,nj
        k=ipoint(ind(1,i),nx,6)
        a=a+c(k)
        tc(k)=.false.
      enddo
      av=a/dble(nj)
!
      if (tout) then
        a=zero
        do i=1,nj
          k=ipoint(ind(1,i),nx,6)
          a=a+(c(k)-av)**2
        enddo
        a=sqrt(a/dble(nj))
        if (abs(av).gt.toll) then
          b=a/av*ten*two
          if (b.gt.ten*ten) then
            do i=1,nj
              k=ipoint(ind(1,i),nx,6)
              write(oout,16) c(k),av,a,(ind(j,i),j=1,nx)
            enddo
          endif
        endif
      endif
!
      do i=1,nj
        k=ipoint(ind(1,i),nx,6)
        c(k)=av
      enddo
!
16    format(3(1x,f15.2),1x,8i1)
      return
      end
!----------------------------------------------------------------------
      function ipoint(ii,n,ij)
!----------------------------------------------------------------------
      implicit none
      integer n,ii(n),ij,k,i,ipoint
!
      k=ii(1)
      do i=2,n
        k=k+(ii(i)-1)*ij**(i-1)
      enddo
      ipoint=k
!
      return
      end
!-----------------------------------------------------------------------
      subroutine bulkmod(tag,st,c)
!-----------------------------------------------------------------------
      use constants
      use controls
      implicit none
      real*8 st(3,3),c(6,6)
      integer i,j,k,m,ou,ipiv(6)
      real*8 cijd(6,6),c2(6,6),a,bm,eig(6),vst(6)
      integer, parameter :: lwork = 3*6-1
      real*8 work(lwork)
      character tag*4
!
      do k=1,6*6
        cijd(k,1)=zero
      enddo
      call p2v(st,vst)
      call birch_coefficients(vst,cijd)
!
      do i=1,6*6
        c2(i,1)=c(i,1)-cijd(i,1)
      enddo
!
!     Compliance
      do k=1,6*6
        cijd(k,1)=c2(k,1)
      enddo
      i=6
      call dgetrf(i,i,cijd,i,ipiv,k)
      if (k.ne.0) stop 6761
      call dgetri(i,cijd,i,ipiv,work,lwork,k)
      if (k.ne.0) stop 6762
!
      bm=cijd(1,1)+cijd(2,2)+cijd(3,3)+ &
     &   cijd(1,2)+cijd(1,3)+cijd(2,3)+ &
     &   cijd(2,1)+cijd(3,1)+cijd(3,2)
      bm=one/bm
!
      ou=18
      open(ou,file=trim(dname)//'moduli.'//tag)
      write(ou,*) 'Bulk modulus'
      write(ou,20) bm
      write(ou,*) 'E1 E2 E3 G1 G2 G3'
      write(ou,20) (one/cijd(k,k),k=1,6)
      write(ou,*) 'Vxy Vyx Vxz Vzx Vyz Vzy'
      write(ou,20)-cijd(1,2)/cijd(2,2),-cijd(2,1)/cijd(1,1), &
     &            -cijd(1,3)/cijd(3,3),-cijd(3,1)/cijd(1,1), &
     &            -cijd(2,3)/cijd(2,2),-cijd(3,2)/cijd(3,3)

      eig(:)=zero
      do i=1,6
        do j=i,6
          a=(c2(i,j)+c2(j,i))/two
          cijd(i,j)=a
          cijd(j,i)=a
        enddo
      enddo
      call dsyev('V','U', 6, cijd, 6, eig , work, lwork, i)
      if (i.ne.0) stop 6464
!
      do i=1,6
        do j=i,6
          a=(c2(i,j)+c2(j,i))/two
          c2(i,j)=a
          c2(j,i)=a
        enddo
      enddo
!
      write(ou,*) 'Stress tensor of the reference state'
      write(ou,20) (vst(i),i=1,6)
      write(ou,*) 'SOECs'
      do i=1,6
        write(ou,20) (c(i,k),k=1,6)
      enddo
      write(ou,*) 'Symmetrized Birch tensor'
      do i=1,6
        write(ou,20) (c2(i,k),k=1,6)
      enddo
      write(ou,*) 'Eigenvalues/vectors of symmetrized Birch tensor'
      do m=1,6
        write(ou,'(1x,f10.2,4x,6(1x,f8.3))') eig(m),(cijd(i,m),i=1,6)
      enddo
      close(ou)
20    format(6(1x,f10.2))
!
      return
      end
!-----------------------------------------------------------------------
      subroutine birch_coefficients(s,bs)
!-----------------------------------------------------------------------
      implicit none
      integer ia,ja,i,j,k,l
      real*8 s(6),c1(3,3),bcoef,a(2),bs(6,6)
!
      c1(1,1)=s(1)
      c1(2,2)=s(2)
      c1(3,3)=s(3)
      c1(2,3)=s(4)
      c1(1,3)=s(5)
      c1(1,2)=s(6)
      c1(3,2)=c1(2,3)
      c1(3,1)=c1(1,3)
      c1(2,1)=c1(1,2)
!
      do ia=1,6
        call ij(ia,i,j)
        do ja=ia,6
          call ij(ja,k,l)
!           
          a(1)=bcoef(i,j,k,l,c1)
          a(2)=bcoef(k,l,i,j,c1)
          bs(ia,ja)=a(1)
          bs(ja,ia)=a(2)
!
        enddo
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine ij(ia,i,j)
!-----------------------------------------------------------------------
      implicit none
      integer ia,i,j
!
      if (ia.lt.1.or.ia.gt.6) stop 1000
      goto (10,20,30,40,50,60),ia
10    i=1
      j=1
      goto 70
20    i=2
      j=2
      goto 70
30    i=3
      j=3
      goto 70
40    i=3
      j=2
      goto 70
50    i=3
      j=1
      goto 70
60    i=2
      j=1
      goto 70
70    continue
      return
      end
!-----------------------------------------------------------------------
      function delta(i,j)
!-----------------------------------------------------------------------
      implicit none
      integer i,j
      real*8 delta
      if (i.eq.j) then
        delta=1.d0
      else
        delta=0.d0
      endif
      return
      end
!-----------------------------------------------------------------------
      function bcoef(i,j,k,l,s)
!-----------------------------------------------------------------------
      implicit none
      integer i,j,k,l
      real*8 delta,s(3,3),bcoef
!
      bcoef=s(i,k)*delta(j,l)+ &
     &      s(i,l)*delta(j,k)+ &
     &      s(j,k)*delta(i,l)+ &
     &      s(j,l)*delta(i,k)-2.d0*s(i,j)*delta(k,l)
      bcoef=bcoef*0.5d0
!
      return
      end
!-----------------------------------------------------------------------
      subroutine p2v(p,v)
!-----------------------------------------------------------------------
      implicit none
      real*8 v(6),p(3,3)
      integer i
!
      do i=1,3
        v(i)=p(i,i)
      enddo
      v(4)=p(2,3)
      v(5)=p(1,3)
      v(6)=p(1,2)
!
      return
      end
!----------------------------------------------------------------------
      subroutine impose_material_symmetries
!----------------------------------------------------------------------
      use constants
      use controls
      use crystal
      use tensors
      implicit none
      integer i,k,n,ou
      real*8 , allocatable :: s(:,:,:)
      character tag*4
!
      if (.not.thomo) return
!
      n=0
      do k=1,nsym
        if (tcry(k)) n=n+1
      enddo
      allocate(s(3,3,n))
!
      i=0
      do k=1,nsym
        if (tcry(k)) then
          i=i+1
          s(:,:,i)=sym(:,:,k)
        endif
      enddo
      if (i.ne.n) stop 1111
!
!
      call symmetrization(n,s)
!
      tag='aecs'
      call print_hoecs(tag)
!
      ou=ohom
      open(ou,file=trim(dname)//'hoecs.'//tag)
      write(ou,*)'CELL_PARAMETER (bohr)'
      do i=1,3
        write(ou,*)(box(k,i)/bohr,k=1,3)
      enddo
      write(ou,*)'CAUCHY (GPa)'
      do i=1,3
        write(ou,*)(cau(k,i,iref),k=1,3)
      enddo
!
      if (.not.tc2) return
      write(ou,*)'2ECS (GPa)'
      write(ou,*)(c2r(i,1),i=1,6*6)
!
      if (.not.tc3) return
      write(ou,*)'3ECS (GPa)'
      write(ou,*)(c3r(i,1,1),i=1,6*6*6)
!
      if (.not.tc4) return
      write(ou,*)'4ECS (GPa)'
      write(ou,*)(c4r(i,1,1,1),i=1,6*6*6*6)
!
      if (.not.tc5) return
      write(ou,*)'5ECS (GPa)'
      write(ou,*)(c5r(i,1,1,1,1),i=1,6*6*6*6*6)
!
      if (.not.tc6) return
      write(ou,*)'6ECS (GPa)'
      write(ou,*)(c6r(i,1,1,1,1,1),i=1,6*6*6*6*6*6)
!
      if (.not.tc7) return
      write(ou,*)'7ECS (GPa)'
      write(ou,*)(c7r(i,1,1,1,1,1,1),i=1,6*6*6*6*6*6*6)
!
      if (.not.tc8) return
      write(ou,*)'8ECS (GPa)'
      write(ou,*)(c8r(i,1,1,1,1,1,1,1),i=1,6*6*6*6*6*6*6*6)
!
      return
      end
!----------------------------------------------------------------------
      subroutine symmetrization(n,s)
!----------------------------------------------------------------------
      use constants
      use tensors
      implicit none
      integer i,n
      real*8 s(3,3,n)
!
      if (tc2) c2a(:,:)=zero
      if (tc3) c3a(:,:,:)=zero
      if (tc4) c4a(:,:,:,:)=zero
      if (tc5) c5a(:,:,:,:,:)=zero
      if (tc6) c6a(:,:,:,:,:,:)=zero
      if (tc7) c7a(:,:,:,:,:,:,:)=zero
      if (tc8) c8a(:,:,:,:,:,:,:,:)=zero
      do i=1,n
        call arotation(s(1,1,i))
        if (tc2) c2a(:,:)             = c2a(:,:)            +c2p(:,:)
        if (tc3) c3a(:,:,:)           = c3a(:,:,:)          +c3p(:,:,:)
        if (tc4) c4a(:,:,:,:)         = c4a(:,:,:,:)        +c4p(:,:,:,:)
        if (tc5) c5a(:,:,:,:,:)       = c5a(:,:,:,:,:)      +c5p(:,:,:,:,:)
        if (tc6) c6a(:,:,:,:,:,:)     = c6a(:,:,:,:,:,:)    +c6p(:,:,:,:,:,:)
        if (tc7) c7a(:,:,:,:,:,:,:)   = c7a(:,:,:,:,:,:,:)  +c7p(:,:,:,:,:,:,:)
        if (tc8) c8a(:,:,:,:,:,:,:,:) = c8a(:,:,:,:,:,:,:,:)+c8p(:,:,:,:,:,:,:,:)
      enddo
      if (tc2) c2r(:,:)               = c2a(:,:)/dble(n)
      if (tc3) c3r(:,:,:)             = c3a(:,:,:)/dble(n)
      if (tc4) c4r(:,:,:,:)           = c4a(:,:,:,:)/dble(n)
      if (tc5) c5r(:,:,:,:,:)         = c5a(:,:,:,:,:)/dble(n)
      if (tc6) c6r(:,:,:,:,:,:)       = c6a(:,:,:,:,:,:)/dble(n)
      if (tc7) c7r(:,:,:,:,:,:,:)     = c7a(:,:,:,:,:,:,:)/dble(n)
      if (tc8) c8r(:,:,:,:,:,:,:,:)   = c8a(:,:,:,:,:,:,:,:)/dble(n)
!
      return
      end
!----------------------------------------------------------------------
      subroutine arotation(mat)
!----------------------------------------------------------------------
      use constants
      use tensors
      implicit none
      integer i,j,k,l,m,n,ii,jj
      integer p,q,r,s,t,u,pp,qq
      real*8 a,mat(3,3),vq(6,6)
!
!
      vq(:,:)=zero
      do i=1,3
        do j=1,3
          vq(i,j)=mat(i,j)**2
        enddo
      enddo
      vq(1,4)=two*mat(1,2)*mat(1,3)
      vq(1,5)=two*mat(1,3)*mat(1,1)
      vq(1,6)=two*mat(1,1)*mat(1,2)
      vq(2,4)=two*mat(2,2)*mat(2,3)
      vq(2,5)=two*mat(2,3)*mat(2,1)
      vq(2,6)=two*mat(2,1)*mat(2,2)
      vq(3,4)=two*mat(3,2)*mat(3,3)
      vq(3,5)=two*mat(3,3)*mat(3,1)
      vq(3,6)=two*mat(3,1)*mat(3,2)
!
      vq(4,1)=    mat(2,1)*mat(3,1)
      vq(5,1)=    mat(3,1)*mat(1,1)
      vq(6,1)=    mat(1,1)*mat(2,1)
      vq(4,2)=    mat(2,2)*mat(3,2)
      vq(5,2)=    mat(3,2)*mat(1,2)
      vq(6,2)=    mat(1,2)*mat(2,2)
      vq(4,3)=    mat(2,3)*mat(3,3)
      vq(5,3)=    mat(3,3)*mat(1,3)
      vq(6,3)=    mat(1,3)*mat(2,3)
!
      vq(4,4)=mat(2,2)*mat(3,3)+mat(2,3)*mat(3,2)
      vq(4,5)=mat(2,3)*mat(3,1)+mat(2,1)*mat(3,3)
      vq(4,6)=mat(2,1)*mat(3,2)+mat(2,2)*mat(3,1)
      vq(5,4)=mat(3,2)*mat(1,3)+mat(3,3)*mat(1,2)
      vq(5,5)=mat(3,3)*mat(1,1)+mat(3,1)*mat(1,3)
      vq(5,6)=mat(3,1)*mat(1,2)+mat(1,1)*mat(3,2)
      vq(6,4)=mat(1,2)*mat(2,3)+mat(1,3)*mat(2,2)
      vq(6,5)=mat(1,3)*mat(2,1)+mat(1,1)*mat(2,3)
      vq(6,6)=mat(1,1)*mat(2,2)+mat(1,2)*mat(2,1)
!
      if (.not.tc2) return
      do i=1,6 ; do j=1,6
        a=zero
        do p=1,6 ; do q=1,6
          a=a+vq(i,p)*vq(j,q)*c2r(p,q)
        enddo ; enddo
        c2p(i,j)=a
      enddo ; enddo
!
      if (.not.tc3) return
      do i=1,6 ; do j=1,6 ; do k=1,6
        a=zero
        do p=1,6 ; do q=1,6 ; do r=1,6
          a=a+vq(i,p)*vq(j,q)*vq(k,r)*c3r(p,q,r)
        enddo ; enddo ; enddo
        c3p(i,j,k)=a
      enddo ; enddo ; enddo
!
      if (.not.tc4) return
      do i=1,6 ; do j=1,6 ; do k=1,6 ; do l=1,6
        a=zero
        do p=1,6 ; do q=1,6 ; do r=1,6 ; do s=1,6
          a=a+vq(i,p)*vq(j,q)*vq(k,r)*vq(l,s)*c4r(p,q,r,s)
        enddo ; enddo ; enddo ; enddo
        c4p(i,j,k,l)=a
      enddo ; enddo ; enddo ; enddo
!
      if (.not.tc5) return
      do i=1,6 ; do j=1,6 ; do k=1,6 ; do l=1,6 ; do m=1,6
        a=zero
        do p=1,6 ; do q=1,6 ; do r=1,6 ; do s=1,6 ; do t=1,6
          a=a+vq(i,p)*vq(j,q)*vq(k,r)*vq(l,s)*vq(m,t)*c5r(p,q,r,s,t)
        enddo ; enddo ; enddo ; enddo ; enddo
        c5p(i,j,k,l,m)=a
      enddo ; enddo ; enddo ; enddo ; enddo
!
      if (.not.tc6) return
      do i=1,6 ; do j=1,6 ; do k=1,6 ; do l=1,6 ; do m=1,6 ; do n=1,6
        a=zero
        do p=1,6 ; do q=1,6 ; do r=1,6 ; do s=1,6 ; do t=1,6 ; do u=1,6
          a=a+vq(i,p)*vq(j,q)*vq(k,r)*vq(l,s)*vq(m,t)*vq(n,u)*c6r(p,q,r,s,t,u)
        enddo ; enddo ; enddo ; enddo ; enddo ; enddo
        c6p(i,j,k,l,m,n)=a
      enddo ; enddo ; enddo ; enddo ; enddo ; enddo
!
      if (.not.tc7) return
      do i=1,6 ; do j=1,6 ; do k=1,6 ; do l=1,6 ; do m=1,6 ; do n=1,6 ; do ii=1,6
        a=zero
        do p=1,6 ; do q=1,6 ; do r=1,6 ; do s=1,6 ; do t=1,6 ; do u=1,6 ; do pp=1,6
          a=a+vq(i,p)*vq(j,q)*vq(k,r)*vq(l,s)*vq(m,t)*vq(n,u)*vq(ii,pp)*c7r(p,q,r,s,t,u,pp)
        enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo
        c7p(i,j,k,l,m,n,ii)=a
      enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo
!
      if (.not.tc7) return
      do i=1,6 ; do j=1,6 ; do k=1,6 ; do l=1,6 ; do m=1,6 ; do n=1,6 ; do ii=1,6 ; do jj=1,6
        a=zero
        do p=1,6 ; do q=1,6 ; do r=1,6 ; do s=1,6 ; do t=1,6 ; do u=1,6 ; do pp=1,6 ; do qq=1,6
          a=a+vq(i,p)*vq(j,q)*vq(k,r)*vq(l,s)*vq(m,t)*vq(n,u)*vq(ii,pp)*vq(jj,qq)*c8r(p,q,r,s,t,u,pp,qq)
        enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo
        c8p(i,j,k,l,m,n,ii,jj)=a
      enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo
      return
      end
