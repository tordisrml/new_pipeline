c-----------------------------------------------------------------------
      program accurac1
c.......................................................................
c     Description: Prepares three datafiles for the calculation of
c                  approximate accuracy of genetic evaluation under an
c                  animal model.
c.......................................................................
c           Input: radnrkodi
c                  dmuped.txt
c                  uppl.txt
c          Output: accuraci.dat
c                  accuracs.dat
c                  accuracd.dat
c.......................................................................
c      References: Meyer, K. 1989. Approximate Accuracy of Genetic
c                  Evaluations under an Animal Model. Livestock
c                  Production Science, 21, 87-100.
c.......................................................................
c      Written by Agust Sigurdsson, Department of animal breeding and
c                 genetics. Swedish University of Agric. Sci., Uppsala.
c-----------------------------------------------------------------------
c     Modified by Jon Hjalti Eiriksson, Radgjafamidstod landbunadarins
c     May 2018 for test day model.
!      parameter(noani = 453604)
      parameter(nofix =    3)
!      parameter(nogrp =    36)

      integer*4, allocatable :: ix(:),s(:),d(:),off1(:)
     +           ,off2(:),ps(:),fix(:,:)
      integer*2, allocatable :: sex(:)
	  integer*2 hjord
      character*7 yearmonth
!      character pth*29
!      pth='/home/agust/skuggi_kyr/mjolk/'


      open(100,file='control.txt')  !control skráin er nauðsynleg. Er með upplýsingum og til dæmis hversu mörgum mælingum má búast við til þess að skilgreina vigra og fylki rétt
      read(100,*)yearmonth
      close(100)

	  open(13,file='uppl.txt')
	    read(13,*)noani
	  close(13)

	  allocate(ix(noani),s(noani),d(noani),off1(noani)
     +           ,off2(noani),ps(noani))
	  allocate(fix(noani,nofix),sex(noani))

c.....NULLUN
!      call nolli(ix,noani)
!      call nolli(s,noani)
!      call nolli(d,noani)
!      call nolli(off1,noani)
!      call nolli(off2,noani)
!      call nolli(ps,noani)
!      call nolli2(sex,noani)
      do i=1,noani
	    ix(i)=0
		s(i)=0
		d(i)=0
		off1(i)=0
		off2(i)=0
		ps(i)=0
       do j=1,nofix
        fix(i,j)=0
       enddo
      enddo
      hjord=0
c................................
      open(10,file='../'//yearmonth//'/dmu_data/radnrkodi')
      open(11,file='../'//yearmonth//'/dmu_data/dmu_ped.txt')
!      open(12,file=pth//'impl/data/sex.txt',status='old')

      open(20,file='../'//yearmonth//'/dmu_data/accuraci.dat')
      open(21,file='../'//yearmonth//'/dmu_data/accuracs.dat')
      open(22,file='../'//yearmonth//'/dmu_data/accuracd.dat')

      write(*,*)'Reading data...'
c...vegna gen.gr.
!      do i=1,nogrp
!        read(10,*)hjord,(fix(i,j),j=1,nofix)
!        read(11,*)ixxj,s(i),d(i)
!      enddo
  111 format(30x,3(1x,i5),i2)
      do 199 i=1,noani
        read(10,111)(fix(i,j),j=1,nofix),sex(i)
        read(11,*)ixxj,s(i),d(i)
        if(s(i).le.0)s(i)=0
		if(mod(i,10000).eq.0)write(*,*)d(i),s(i)
!        if(s(i).gt.0)s(i)=s(i)-nogrp
        if(d(i).le.0)d(i)=0
!        if(d(i).gt.0)d(i)=d(i)-nogrp
!        read(12,497)
        ix(i)=i
!        off1(i)=0
!        off2(i)=0
  199 continue
!  497 format(i15,1x,i8,1x,i2,i3,3(1x,5i),i2)
      write(*,*)'        ...done'
      write(*,*)'Counting offsprings...'
      do 188 i=1,noani
        if(mod(i,10000).eq.0)then
           print *,'Working with animal:',i
        endif
        if(sex(i).eq.1)then
          do 187 is=1,noani
            if(s(is).eq.ix(i))then
              if(d(is).gt.0)then
                off1(i)=off1(i)+1
              else
                off2(i)=off2(i)+1
              endif
            endif
  187     continue
        else
          do 186 id=1,noani
            if(d(id).eq.ix(i))then
              if(s(id).gt.0)then
                off1(i)=off1(i)+1
              else
                off2(i)=off2(i)+1
              endif
            endif
  186     continue
        endif
  188 continue
      write(*,*)'               ...done'
      write(*,*)'Creating - accuraci.dat'
      do 179 i=1,noani
        if(fix(i,1).gt.0)then
          irec=1
        else
          irec=0
        endif
      write(20,399)irec,fix(i,1),fix(i,2),fix(i,3),
     +              ix(i),off1(i),
     +   off2(i),s(i),d(i),sex(i)
  179 continue
      write(*,*)'                ...done'
      write(*,*)'Creating - accuracs.dat'
      do 169 i=1,noani
  169   ps(i)=s(i)

      call sort2b(noani,ps,ix)

      do 168 i=1,noani
        if(s(ix(i)).gt.0)then
          if(fix(ix(i),1).gt.0)then
            irec=1
          else
            irec=0
          endif
        write(21,399)irec,fix(ix(i),1),fix(ix(i),2),fix(ix(i),3),
     +    ix(i)
     +    ,off1(ix(i)),off2(ix(i)),s(ix(i)),d(ix(i)),sex(ix(i))
        endif
  168 continue
      write(*,*)'                ...done'
      write(*,*)'Creating - accuracd.dat'
      do 159 i=1,noani
        ix(i)=i
        ps(i)=d(i)
  159 continue
      call sort2b(noani,ps,ix)

      do 158 i=1,noani
        if(d(ix(i)).gt.0)then
          if(fix(ix(i),1).gt.0)then
            irec=1
          else
            irec=0
          endif
        write(22,399)irec,fix(ix(i),1),fix(ix(i),2),fix(ix(i),3),
     +  ix(i),off1(ix(i)),
     +   off2(ix(i)),s(ix(i)),d(ix(i)),sex(ix(i))
        endif
  158 continue
      write(*,*)'                ...done'
  399 format(i1,i5,2i2,i6,2i4,2i7,i1)
      close(10)
      close(11)
      close(12)
      close(20)
      close(21)
      close(22)
      stop
      end

	  subroutine sort2b(n,ra,rb)
*      dimension ra(n),rb(n)
      integer*4 ra(n),rb(n)
      l=n/2+1
      ir=n
   10 continue
         if(l.gt.1)then
            l=l-1
            rra=ra(l)
            rrb=rb(l)
         else
            rra=ra(ir)
            rrb=rb(ir)
            ra(ir)=ra(1)
            rb(ir)=rb(1)
            ir=ir-1
            if(ir.eq.1)then
               ra(1)=rra
               rb(1)=rrb
               return
            endif
         endif
         i=l
         j=l+l
   20    if(j.le.ir)then
            if(j.lt.ir)then
               if(ra(j).lt.ra(j+1))then
                 j=j+1
               endif
            endif
            if(rra.lt.ra(j))then
               ra(i)=ra(j)
               rb(i)=rb(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
            go to 20
         endif
         ra(i)=rra
         rb(i)=rrb
         go to 10
      end

!      include '/home/agust/agusts/assub/sorting.f'
!      include '/home/agust/agusts/assub/nolli.f'
!      include '/home/agust/agusts/assub/nolli2.f'
