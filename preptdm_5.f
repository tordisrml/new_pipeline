! Program for handling TDM data from Huppa and prepare for breeding value estimation.
! JHE March 2018
!
! Files needed:
! *.csv  Records. Does not need to be sorted.
!  pedegree file form Huppa sorted.
!
! Files written:
! tdm1.dat      For MY evaluation
! tdm2.dat      For the FY, PY and SCS
! tdm3.dat      For FP and PP
! dmuped.txt    The Pedigree file for DMU
! radnrkodi     Long id to new id and more
! uppl.txt      Information for solsaman.f

      implicit none !Svona byrja flest fortran forrit

    !Skilgreina breytur
	  integer*8 inid,idnow,keep1,keep2
	  integer*8, allocatable :: id(:),cow(:)     !Þarf lengri integer fyrir einstaklingsnumerin
	  integer recmax,cowmax,herdmax,yearmax,firstyear,pedmax,stdyear  !Fyrir stuðla úr control
	  integer tdmax,ruglpedmax, mlr
	  integer inherd, inlact, bdy, bdm, bdd, inn  !Fyrir breytur sem eru lesnar inn
	  integer, allocatable :: herdv(:),indexa(:)
      integer inbdy,inbdm,inbdd,intdy,intdm,intdd,incdy,incdm,incdd  !Dagsetningar teknar inn i butum
      integer idim, icd, ibd, ica, itd, gr, grs,nu
	  integer nuh,pl,i,j,k,l,herdnow,ql,last,lasth,t1,t2,t3
	  integer mdays(12),dphg(50),sphg(50),ndphg,nsphg,nphg(2,50)
	  integer, allocatable :: dim(:),bd(:),cy(:),ca(:),cm(:)
	  integer, allocatable :: ty(:),tm(:)
	  integer, allocatable :: herd(:),lact(:),cd(:),td(:),htds(:)
	  integer, allocatable :: hy(:),hys(:),htd(:)
	  integer, allocatable :: cc1(:,:),cc2(:,:),cc3(:,:)
	  integer, allocatable :: scount(:,:),ccount(:,:),tdscnow(:)
	  integer, allocatable :: yeco(:,:,:),yesco(:,:,:),tdinhco(:,:)
	  integer ncoun,sncount,ci1,ci2,tempc(3),tdnr(3),tdsnr(3)
	  integer year,yc,res,r,ysc,ress,rs,hynr(3),hysnr(3)
	  integer, allocatable :: vis_td(:), td_td(:), mjsk_td(:)
	  integer recbumax,dag,tdnu,taltdhop(4),tdnus,taltdhops(4)

      real inmy, infy, inpy, infper, inpper, inscc   !Fyrir innlestur real (float) breyta
	  real, allocatable :: my(:),fy(:),py(:),fp(:)
	  real, allocatable :: scs(:),pp(:)
	  character*10 inbd,intd,incd
	  character*30 filename1, filename2

	  logical, allocatable :: milkok(:), samok(:),samok_td(:)
      logical binsok, langsok, samplok
	  logical lac1ok,lac2ok,lac3ok
	  logical ok

	  integer*8 extra(100)
	  integer*8, allocatable :: pedid(:),dam(:),sire(:)
	  integer, allocatable :: ar(:),kyn(:),recno(:),fix(:,:)
	  integer land,extra_year(100)
	  real stdim,lp1,lp2,lp3,lp4,myout(3),fyout(3),pyout(3),wil
	  real scsout(3),ppout(3),fpout(3)

      integer*8 einst(3), arsums, arsumd
      integer*8 gognid,biljon
	  integer*8, allocatable :: draugur(:),auka(:)
	  integer, allocatable :: draug_ar(:),raddraug(:)
      integer taldraug,nrped,stdan
      integer einstut(3),allsped
      integer damtal, siretal,bil,utdraugur,aukatal,stdanno
	  integer, allocatable :: by(:),hygc(:),htdgc(:),hysgc(:)
	  integer aukaut,faertal,fjoldi,hyp,stadar
      real meddam, medsire
      logical logdraug

      !Hérna er búið að skilgreina breyturnar

      open(10,file='control.txt')  !control skráin er nauðsynleg. Er með upplýsingum og til dæmis hversu mörgum mælingum má búast við til þess að skilgreina vigra og fylki rétt
  	  read(10,'(i8)')recmax       ! hámarksfjöldi mælinga. Þarf sumsé bara öruggega að vera hærri tala en fjöldi mælinga
  	  read(10,'(i8)')cowmax       !Hámarksfjöldi kúa
  	  read(10,'(i8)')herdmax      !hámarksfjöldi búa
  	  read(10,'(i8)')yearmax      !Hámarksfjöldi ára
  	  read(10,'(i8)')firstyear    !Sleppir gögnum frá því fyrir þetta ártal.
  	  read(10,'(i8)')pedmax       !Hve stór ættskráin má vera
  	  read(10,'(i8)')stdyear      !Ár til að miða staðalhópinn við. Er með í þessu forriti vegna þess að það merkir inni í radnrkodi hvaða kýr eiga að vera með í staðalárgangi byggt á fæðingarári og hvort þær eru með mælingar
  	  read(10,'(i8)')ruglpedmax   !Hversu marga aukaeinstaklinga má búa til ef það er eitthvert rugl, til dæmis foreldrar sem eru ekki sjálfir til
  	  read(10,'(i8)')mlr	        !Hversu margar mælingar þurfa að vera að lágmarki á mjaltaskeiði til að það sé tekið með. Ég hef verið með 2. Var með einhver rök fyrir því, er ekki viss um að það sé eitthvert vit í þeim
  	  read(10,*)filename1         !Nafnið á gagnaskránni
  	  read(10,*)filename2	        !Nafnið á ættskránni
  	  close(10)

	  open(11,file='phg/damphg.txt') !les hvaða ár á að miða við í draugahópunum fyrir mæður úr sér skrá.
    !Lesa draugahópaárin
	  i=0
	  do
	    i=i+1
	    read(11,*,end=97)dphg(i)
	  enddo
   97 close(11)
      ndphg=i-1
      i=0
	  open(12,file='phg/sirephg.txt') !Draugahópar feður
	  do
	    i=i+1
	    read(12,*,end=98)sphg(i)
      enddo
   98 close(12)
      nsphg=i-1

	  write(*,*)'Data from ',filename1
	  tdmax=yearmax*12

      write(*,*)recmax,cowmax,herdmax,yearmax,firstyear,pedmax,stdyear

	  allocate(id(recmax))
	  allocate(dim(recmax),bd(recmax),cy(recmax),ca(recmax))   !Tekur frá minni. Max tölurnar úr control.txt notaðar til að ákvarða hversu mikið á að taka frá
	  allocate(ty(recmax),tm(recmax))
	  allocate(herd(recmax),lact(recmax),cd(recmax),td(recmax))
	  allocate(hy(recmax),hys(recmax),cm(recmax))
	  allocate(indexa(recmax),my(recmax),fy(recmax),py(recmax))
	  allocate(fp(recmax),pp(recmax),scs(recmax))
	  allocate(milkok(recmax),samok(recmax),by(recmax))


! A few standards
! Days of the year at the beginning of each month.
	  mdays(1)=0
	  mdays(2)=31
	  mdays(3)=59
	  mdays(4)=90
	  mdays(5)=120
	  mdays(6)=151
	  mdays(7)=181
	  mdays(8)=212
	  mdays(9)=243
	  mdays(10)=273
	  mdays(11)=304
	  mdays(12)=334


	  open(10,file=filename1) !The main data file

      inn=0      !number of lines read
	  gr=0       !number of useful milk records
	  grs=0      !number of useful milk sample records


! Loop for reading data file
! Only saves records if they are within limits
	  do
	    inmy=0
		infy=0
		inpy=0
		infper=0
		inpper=0
		inscc=0

		if(gr.eq.recmax)write(*,*)'Error: too many records',gr
		if(gr.eq.recmax)goto 999
        inn=inn+1
    !Lesa inn úr csv skránni. ID, fæðingardagur, bú, mjaltaskeið, mælidagur, mjólk, fita, prótein,fitu%, prótein%, frumutala, burðardagur
		read(10,*, end=20)inid,inbd,inherd,inlact    !ATH gagnaskrárnar eru með auka dálki aftast sem er "næsti burður" til að reikna bil á milli burða
     +  ,intd,inmy,infy,inpy,infper,inpper,     !Þarf að vita fyrir python innlestur!
     +  inscc,incd
        samplok=.true.                  ! boolean breyta um það hvort það á að nota sýnaniðurstöðurnar. Er breytt í .false. ef eitthvað er athugavert.

	    if(mod(inn,500000).eq.0)then     ! Þessi lykkja er bara til að skrifa eitthvað á skjáinn fyrir hverja fimmhundruðþúsundustu mælingu
          write(*,*)inn, 'lines read '
	      write(*,*)gr,'milkyield records kept '
          write(*,*)grs,'records from milksamples kept'
		  write(*,*)inid,inbd,inherd,inlact
     +      ,intd,inmy,infy,inpy,infper,inpper,
     +      inscc,incd
        endif

! Working with dates, bd=birthday, cd=calvingday, td=testday
		read(inbd,'(i4,x,i2,x,i2)')inbdy,inbdm,inbdd
		read(intd,'(i4,x,i2,x,i2)')intdy,intdm,intdd
		read(incd,'(i4,x,i2,x,i2)')incdy,incdm,incdd
		ibd=(inbdy-1993)*365+mdays(inbdm)+inbdd !days from 1.1.1993    Býr til dagafjölda frá 1. jan 1993 í stað dagsetninganna
		icd=(incdy-1993)*365+mdays(incdm)+incdd
		itd=(intdy-1993)*365+mdays(intdm)+intdd
! Hlauparsleidretting
        if(ibd.gt.1154) ibd=ibd+1 !1996
		if(ibd.gt.2615) ibd=ibd+1 !2000
		if(ibd.gt.4076) ibd=ibd+1 !2004
		if(ibd.gt.5537) ibd=ibd+1 !2008
		if(ibd.gt.6998) ibd=ibd+1 !2012
		if(ibd.gt.8459) ibd=ibd+1 !2016
		if(ibd.gt.9920) ibd=ibd+1 !2020
		if(ibd.gt.11381) ibd=ibd+1 !2024

        if(icd.gt.1154) icd=icd+1 !1996
		if(icd.gt.2615) icd=icd+1 !2000
		if(icd.gt.4076) icd=icd+1 !2004
		if(icd.gt.5537) icd=icd+1 !2008
		if(icd.gt.6998) icd=icd+1 !2012
		if(icd.gt.8459) icd=icd+1 !2016
		if(icd.gt.9920) icd=icd+1 !2020
		if(icd.gt.11381) icd=icd+1 !2024

        if(itd.gt.1154) itd=itd+1 !1996
		if(itd.gt.2615) itd=itd+1 !2000
		if(itd.gt.4076) itd=itd+1 !2004
		if(itd.gt.5537) itd=itd+1 !2008
		if(itd.gt.6998) itd=itd+1 !2012
		if(itd.gt.8459) itd=itd+1 !2016
		if(itd.gt.9920) itd=itd+1 !2020
		if(itd.gt.11381) itd=itd+1 !2024


!Þessar línur eru síur. Ef það er cycle er hoppað út úr lykkjunni og þessi mæling fær ekki að vera með. Ef mjólkurmagnið virðist í lagi en ekki sýnið er samplok sett á .false.
		idim=itd-icd               !days in milk    !dagar frá burði til mælingar reiknaðir
        ica=icd-ibd                !calving age    !burðaraldur í d-gum reiknaður
		if(inmy.lt.1.or.inmy.gt.55)cycle    !mjólkurmagn
		if(infy.lt.0.05.or.infy.gt.3)samplok=.false.     !fita kg
		if(inpy.lt.0.05.or.inpy.gt.2.4)samplok=.false.     !prótein kg
		if(infper.lt.2.or.infper.gt.8)samplok=.false.     !fitu%
		if(inpper.lt.2.or.inpper.gt.6)samplok=.false.    !protein %
		if(inscc.lt.1.or.inscc.gt.10000)samplok=.false.    !frumutala
        if(inherd.gt.7000000)cycle    !bæjarnúmer, mig minnir að allir bæir með hærra númer væri í færeyjum
        if(idim.lt.5.or.idim.gt.305)cycle     !dagar frá burði 5-305
		if(inlact.eq.1.and.ica.lt.540)cycle    !burðaraldur, sérstakt fyrir hvert mjaltaskeið
		if(inlact.eq.1.and.ica.gt.1350)cycle
		if(inlact.eq.2.and.ica.lt.840)cycle
		if(inlact.eq.2.and.ica.gt.1800)cycle
		if(inlact.eq.3.and.ica.lt.1140)cycle
		if(inlact.eq.3.and.ica.gt.2250)cycle
		if(incdy.lt.firstyear)cycle         !Sleppir mælingum frá ári fyrir það sem var skilgreint í control.txt

        gr=gr+1                !Telur fjölda mælinga sem fá að vera með
        if(samplok)grs=grs+1            !Telur fjölda efnamælinga sem fá að vera með
!       Put records into vectors
        samok(gr)=samplok       !Vigur tengdur hverri mælingu sem segir til um hvort efnamælingin má vera með
		milkok(gr)=.true.    !Vigur tengdur hverri mælingu sem segir til um hvort mælingin má vera með. Er sett á .true. fyrir allar núna, en er sett á .false. seinna ef eitthvað er að
		id(gr)=inid
		my(gr)=inmy
        if(samplok)then
  		  fy(gr)=infy
		  py(gr)=inpy
		  fp(gr)=infper
		  pp(gr)=inpper
          scs(gr)=log(inscc/100)/log(2.00)+3
        endif
		dim(gr)=idim
		cy(gr)=incdy
		cm(gr)=incdm
		ty(gr)=intdy
		tm(gr)=intdm
		cd(gr)=icd
		herd(gr)=inherd
		lact(gr)=inlact
		td(gr)=itd
		by(gr)=inbdy

! changing calving age to 30 day groups:
		ica=ica/30
		  if(inlact.eq.1)then  !ef mjaltaskeið 1
		    if(ica.lt.21)ica=20   !sameinar yngstu hópana
			if(ica.gt.40)ica=41  !sameinar elstu hópana
			ica=ica-19
		  endif
		  if(inlact.eq.2)then
		    if(ica.lt.33)ica=32
			if(ica.gt.54)ica=55
			ica=ica-31
		  endif
		  if(inlact.eq.3)then
		    if(ica.lt.44)ica=43
			if(ica.gt.65)ica=66
			ica=ica-42
		  endif
		ca(gr)=ica

      enddo

   20  close(10)
	  write(*,*)'Reading of data file finished'
	  write(*,*)inn, 'lines read in total'
	  write(*,*)gr,'milkyield records kept, the maximum is ',recmax

	  allocate(herdv(herdmax),cow(cowmax))

!    Get lists of cows and herds
!  For checking number of records and more
!Næsta lykkja rennir í gegn um mælingarnar og skilar lista yfir kýr (cow) og bú (herdv), auk fjölda þeirra (nu og nuh).
!Ég veit ekki alveg afhverju þetta er svona flókinn kóði, en ég held að þetta sé allt sem hann gerir
      nu=1
      nuh=1
	  last=1
	  lasth=1
      do i=1, gr
		if(nu.eq.cowmax)write(*,*)'Error: too many cows',nu
		if(nuh.eq.herdmax)write(*,*)'Error: too many herds',nuh
		if(nu.eq.cowmax.or.nuh.eq.herdmax)goto 999
        if(mod(i,100000).eq.0)then
		write(*,*)i,' records, ',nu,' cows and ',
     +     nuh,' herds.'
	    endif
        idnow=id(i)
        if(idnow.eq.cow(last))then
          goto 21
		else
          if(idnow.gt.cow(nu))then
		    nu=nu+1
            cow(nu)=idnow
            last=nu
            goto 21
          endif
		  if(langsok(cow,nu,idnow,ql))then
		    last=ql
		    goto 21
		  else
		    keep1=idnow
		    do j=ql,nu
			  keep2=keep1
			  keep1=cow(j)
			  cow(j)=keep2
            enddo
            nu=nu+1
            cow(nu)=keep1
            last=ql
		  endif
		endif

   21   herdnow=herd(i)
        if(herdnow.eq.herdv(lasth))then
          cycle
		else
		  if(herdnow.gt.herdv(nuh))then
		    nuh=nuh+1
			herdv(nuh)=herdnow
			last=nuh
			cycle
	      endif
		  if(binsok(herdv,nuh,herdnow,pl))then
		    lasth=pl
		    cycle
		  else
		    keep1=herdnow
		    do j=pl,nuh
			  keep2=keep1
			  keep1=herdv(j)
			  herdv(j)=keep2
            enddo
            nuh=nuh+1
            herdv(nuh)=keep1
            lasth=pl
		  endif
		endif
      enddo

 	  write(*,*)nu, ' cows with records'
	  write(*,*)'from ', nuh, ' herds'


 ! Allocate memory for counting records of each cow

	  allocate (ccount(nu,3),scount(nu,3))
	  allocate (cc1(nu,11),cc2(nu,11),cc3(nu,11))

! fill with 0

      do i=1, nu
	   ccount(i,1)=0
	   ccount(i,2)=0
	   ccount(i,3)=0
	   scount(i,1)=0
	   scount(i,2)=0
	   scount(i,3)=0
      enddo
      write(*,*)'Counting records of cows'

!Loop for counting records of cows. Also finds if two records in same day, puts in an average.
	  do i=1, gr                               !lykkja yfir mælingar
	   if(mod(i,300000).eq.0)write(*,*)i
	   if(langsok(cow,nu,id(i),ql))then        !finna hver kýrin er í kúalistanum
	     if(lact(i).eq.1)then                  ! ef fyrsta mjaltaskeið
		   ncoun=ccount(ql,1)+1                   !ncount er til að telja númer hvað þessi mæling verður fyrir þessa kú á þessu mjaltaskeiði. ccount geymir hve mikið á komið áður fyrir kúna
		   if(samok(i))sncount=scount(ql,1)+1    !talning fyrir efnaeiginleika
		   if(ncoun.gt.11)then                   !Eyðir öllu ef komnar eru fleiri en 11 mælingar vegna þess að þá er eitthvað í ólagi. Var helst í elstu gögnunum, að ég held því að mælingarnar voru skráðar á vitlausar kýr, svo sama kýrin fékk mælingar frá tveimur
			 milkok(i)=.false.                   !þetta setur mælinguna sem unnið er með núna á .false., sumsé ekki með
		     do j=1,11                           !Þessi lykkja finnur svo aðrar mælingar þessarar kýr fyrir þetta mjaltaskeið og setur á .false.
      		   milkok(cc1(ql,j))=.false.
			   samok(cc1(ql,j))=.false.
			 enddo
			 cycle
		   endif
		   if(ncoun.eq.1)goto 22
		   do j=1,ncoun-1
		     if(dim(cc1(ql,j)).eq.dim(i))then     !Ef nýja mælingin er skrifuð á sama dag og fyrri mæling kýrinnar er sett inn meðaltal
			   my(cc1(ql,j))=(my(cc1(ql,j))+my(i))/2
			   if(samok(i))then
			     fy(cc1(ql,j))=(fy(cc1(ql,j))+fy(i))/2
			     py(cc1(ql,j))=(py(cc1(ql,j))+py(i))/2
			     fp(cc1(ql,j))=(fp(cc1(ql,j))+fp(i))/2
			     pp(cc1(ql,j))=(pp(cc1(ql,j))+pp(i))/2
			     scs(cc1(ql,j))=(scs(cc1(ql,j))+scs(i))/2
			     samok(i)=.false.
			     sncount=sncount-1
			   endif
			   milkok(i)=.false.
			   ncoun=ncoun-1
			   exit
			 endif
		   enddo
   22      if(milkok(i))cc1(ql,ncoun)=i              ! cc1 er "pointer", þ.e. geymir hvar í aðal gagnavigrunum mælingar hverrar kýr eru. cc1 er fyrir fyrsta mjaltaskeið, cc2 fyrir annað
		   ccount(ql,1)=ncoun
		   if(samok(i))scount(ql,1)=sncount
		 endif


	     if(lact(i).eq.2)then                        !annað mjaltaskeið
		   ncoun=ccount(ql,2)+1
		   if(samok(i))sncount=scount(ql,2)+1
		   if(ncoun.gt.11)then
			 milkok(i)=.false.
		     do j=1,11
      		   milkok(cc2(ql,j))=.false.
			   samok(cc2(ql,j))=.false.
			 enddo
			 cycle
		   endif
		   if(ncoun.eq.1)goto 23
		   do j=1,ncoun-1
		     if(dim(cc2(ql,j)).eq.dim(i))then
			   my(cc2(ql,j))=(my(cc2(ql,j))+my(i))/2
			   if(samok(i))then
			     fy(cc2(ql,j))=(fy(cc2(ql,j))+fy(i))/2
			     py(cc2(ql,j))=(py(cc2(ql,j))+py(i))/2
			     fp(cc2(ql,j))=(fp(cc2(ql,j))+fp(i))/2
			     pp(cc2(ql,j))=(pp(cc2(ql,j))+pp(i))/2
			     scs(cc2(ql,j))=(scs(cc2(ql,j))+scs(i))/2
			     samok(i)=.false.
			     sncount=sncount-1
			   endif
			   milkok(i)=.false.
			   ncoun=ncoun-1
			   exit
			 endif
		   enddo
   23      if(milkok(i))cc2(ql,ncoun)=i
		   ccount(ql,2)=ncoun
		   if(samok(i))scount(ql,2)=sncount
		 endif


	     if(lact(i).eq.3)then                  !þriðja mjaltaskeið
		   ncoun=ccount(ql,3)+1
		   if(samok(i))sncount=scount(ql,3)+1
		   if(ncoun.gt.11)then
			 milkok(i)=.false.
		     do j=1,11
      		   milkok(cc3(ql,j))=.false.
			   samok(cc3(ql,j))=.false.
			 enddo
			 cycle
		   endif
		   if(ncoun.eq.1)goto 24
		   do j=1,ncoun-1
		     if(dim(cc3(ql,j)).eq.dim(i))then
			   my(cc3(ql,j))=(my(cc3(ql,j))+my(i))/2
			   if(samok(i))then
			     fy(cc3(ql,j))=(fy(cc3(ql,j))+fy(i))/2
			     py(cc3(ql,j))=(py(cc3(ql,j))+py(i))/2
			     fp(cc3(ql,j))=(fp(cc3(ql,j))+fp(i))/2
			     pp(cc3(ql,j))=(pp(cc3(ql,j))+pp(i))/2
			     scs(cc3(ql,j))=(scs(cc3(ql,j))+scs(i))/2
			     samok(i)=.false.
			     sncount=sncount-1
			   endif
			   milkok(i)=.false.
			   ncoun=ncoun-1
			   exit
			 endif
		   enddo
   24      if(milkok(i))cc3(ql,ncoun)=i
		   ccount(ql,3)=ncoun
		   if(samok(i))scount(ql,3)=sncount
         endif
       else
	     write(*,*)'Error, cow not found ',id(i)
       endif
	  enddo

!Loop for checking number of records of cows, calving interval
!Fyrst for the milk sample count

      write(*,*)'Removing records if too few samples in a lactation'
	  write(*,*)'Number controled in control.txt'

    !Ny lykkja yfir kýr. Hún eyðir út mjaltaskeiðum ef færri en mlr mælingar eru komnar frá mjaltaskeiðinu og eyðir þá líka seinni mjaltaskeiðum
 	 !Þessi gerir þetta fyrir efnaeiginleikana
	 do i=1, nu
	   tempc(1)=scount(i,1)
	   tempc(2)=scount(i,2)
	   tempc(3)=scount(i,3)
	   lac1ok=.true.
	   lac2ok=.true.
	   lac3ok=.true.

	   if(tempc(1).lt.mlr)then
         lac1ok=.false.
		 lac2ok=.false.
		 lac3ok=.false.
	   endif

	   if(tempc(2).lt.mlr)then
		 lac2ok=.false.
		 lac3ok=.false.
	   endif

	   if(tempc(3).lt.mlr)then
         lac3ok=.false.
	   endif

	   if(.not.lac1ok.and.tempc(1).gt.0)then
   	     do j=1,ccount(i,1)
	       samok(cc1(i,j))=.false.
		 enddo
	   endif
	   if(.not.lac2ok.and.tempc(2).gt.0)then
   	     do j=1,ccount(i,2)
	       samok(cc2(i,j))=.false.
		 enddo
	   endif
	   if(.not.lac3ok.and.tempc(3).gt.0)then
   	     do j=1,ccount(i,3)
	       samok(cc3(i,j))=.false.
		 enddo
	   endif
      enddo

      write(*,*) 'Then for milk records. Also checking if calving ',
     +  'interval is out of limits'
!Here for the MY and ci
  !Þessi er fyrir mjólkina.
	 do i=1, nu
	   tempc(1)=ccount(i,1)
	   tempc(2)=ccount(i,2)
	   tempc(3)=ccount(i,3)
	   lac1ok=.true.
	   lac2ok=.true.
	   lac3ok=.true.

	   if(tempc(1).lt.mlr)then
          lac1ok=.false.
		  lac2ok=.false.
		  lac3ok=.false.
	   endif

	   if(tempc(2).lt.mlr)then
		  lac2ok=.false.
		  lac3ok=.false.
	   else
	     ci1=cd(cc2(i,1))-cd(cc1(i,1))   !smá auka hérna miðað við síðustu lykkju. Tjekkar hversu langt er á milli burða, eyðir mælingum á seinna mjaltaskeiðinu ef það er minna en 295 eða meira en 730 dagar á milli burða
	     if(ci1.lt.295.or.ci1.gt.730)then
		   lac2ok=.false.
		   lac3ok=.false.
		 endif
	   endif

	   if(tempc(3).lt.mlr)then
		  lac3ok=.false.
	   else
	     ci1=cd(cc3(i,1))-cd(cc2(i,1))
	     if(ci1.lt.295.or.ci1.gt.730)then
		   lac2ok=.false.
		   lac3ok=.false.
		 endif
	   endif

	   if(.not.lac1ok.and.tempc(1).gt.0)then
   	     do j=1,tempc(1)
	       milkok(cc1(i,j))=.false.
	       samok(cc1(i,j))=.false.
		 enddo
	   endif
	   if(.not.lac2ok.and.tempc(2).gt.0)then
   	     do j=1,tempc(2)
	       milkok(cc2(i,j))=.false.
	       samok(cc2(i,j))=.false.
		 enddo
	   endif
	   if(.not.lac3ok.and.tempc(3).gt.0)then
   	     do j=1,tempc(3)
	       milkok(cc3(i,j))=.false.
	       samok(cc3(i,j))=.false.
		 enddo
	   endif
	  enddo

	  write(*,*)'Calving interval and number of records check done'
      deallocate(cc1,cc2,cc3,ccount,scount)
	  allocate(ccount(nuh,3),scount(nuh,3))
      allocate(yeco(nuh,3,yearmax),yesco(nuh,3,yearmax))!,tdinhco(nuh,3))
!      allocate(tdco(nuh,3,tdmax),tdsco(nuh,3,tdmax),tdinh(nuh,3,tdmax))

!	 allocate (cc1(nuh,10000),cc2(nuh,5000),cc3(nuh,2000))

	  write(*,*)'Reallocation done'
    !þetta er bara til að fylla vigra af núllum.
      do i=1, nuh
	   ccount(i,1)=0
	   ccount(i,2)=0
	   ccount(i,3)=0
	   scount(i,1)=0
	   scount(i,2)=0
	   scount(i,3)=0

	   do j=1, yearmax
	     yeco(i,1,j)=0
	     yeco(i,2,j)=0
	     yeco(i,3,j)=0
	     yesco(i,1,j)=0
	     yesco(i,2,j)=0
	     yesco(i,3,j)=0
	   enddo

      enddo

      write(*,*)'00'
!Loop for counting records in herds	and years
!      k=0

	  do i=1, gr
	   l=lact(i)
	   if(binsok(herdv,nuh,herd(i),pl).and.milkok(i))then
         if(mod(i,500000).eq.0)write(*,*)i
		 ccount(pl,l)=ccount(pl,l)+1
         year=ty(i)-firstyear+1
		 if(samok(i))then
		   scount(pl,l)=scount(pl,l)+1
           yesco(pl,l,year)=yesco(pl,l,year)+1    !fylki sem geymir fjölda mælinga á búi(pl), mjaltaskeiði(l),ári(year) fyrir efnaeiginleikana
		 endif
         yeco(pl,l,year)=yeco(pl,l,year)+1         !fylki sem geymir fjölda mælinga á búi(pl), mjaltaskeiði(l),ári(year)
		 indexa(i)=pl

	   else
	     if(milkok(i))write(*,*)'Error, herd not found'
		 if(milkok(i))goto 999
	   endif
	  enddo

	  write(*,*)'Making herd-year groups'
      r=0
	  rs=0
    !Hérna kemur frekar flippuð lykkja
	  !Hún byr til bú-ár hópana með raðnúmerum án þess að skrá mælingarnar í hópana strax
	  !yeco og yesco fylkin eru með fjölda á bú-árum.

	  !Í þessari lykkju er fjöldatölunni skipt út fyrir númeri á hópnum.
	  !Ef það eru nógu margar mælingar (29) í hópnum er það frekar auðvelt
	  !Ef það eru færri eru hópar sameinaðir innan bús.
	  !Þetta er gert með því að renna í gegn um búin, og svo árin innan búa. Fyrst er byrjað
	  !með hóp 1, fyrsta ár á fyrsta búi fær það númer. Ef það eru nógu margar kýr í þeim hópi er númerið hækkað
	  !um 1 fyrir næsta bú-ár. Ef það eru of fáar er haldið áfram með sama númer þangað til það eru komnar 29 í hópinn.
	  !Ef það eru of fáar í síðasta ári á búinu fær það bú-ár sama númer og næsta ár á undan.
	  !Mér sýnist líka vera eitthvert öryggi hérna á heildarfjölda á búi, ef hann nær 19 er búið með, ef það er minna virðist mér sem mælingarnar séu bara með síðasta búi.
    do l=1,3
	    write(*,*)'Lactation ',l
	    do i=1,nuh
	      if(ccount(i,l).lt.1)cycle
   		  yc=0
		  ysc=0
		  if(ccount(i,l).gt.19)r=r+1
	      if(scount(i,l).gt.19)rs=rs+1
	      res=ccount(i,l)
		  ress=scount(i,l)
		  do j=1,yearmax
            yc=yc+yeco(i,l,j)
	        ysc=ysc+yesco(i,l,j)
	        res=res-yeco(i,l,j)
		    ress=ress-yesco(i,l,j)
	        yeco(i,l,j)=r
		    yesco(i,l,j)=rs
		    if(yc.gt.29.and.res.gt.29)then
			  r=r+1
			  yc=0
			endif
	        if(ysc.gt.29.and.ress.gt.29)then
			  rs=rs+1
              ysc=0
            endif
		  enddo
	    enddo
		hynr(l)=r
		hysnr(l)=rs
	  enddo

	  do l=1, 3
	    write(*,*)hynr(l),' herd-year groups in lactation ', l
		write(*,*)hysnr(l),' milk component herd-year groups',
     +      ' in lactation ', l
	  enddo

      !Þessi krúttlega lykkja setur svo bú-ár númerin við mælingarnar, hy fyrir mjólk, hys fyrir efni.
	  !hy = herd-year group for MY
	  !hys = herd-year for the other traits
	  do i=1,gr
	    if(milkok(i))then
	      year=ty(i)-firstyear+1
	      hy(i)=yeco(indexa(i),lact(i),year)
		  if(samok(i))then
		    hys(i)=yesco(indexa(i),lact(i),year)
		  else
		    hys(i)=0
		  endif
        endif
	  enddo

      !hérna byrjar svo bú-mælidagadæmið
	  !Herd-test day
	  write(*,*)'making herd-test-day'
	  recbumax=12*250*yearmax        !hámarksfjöldi mælinga á búi.
	  allocate(td_td(recbumax),mjsk_td(recbumax),vis_td(recbumax))
	  allocate(samok_td(recbumax))
	  allocate(htd(recmax),htds(recmax))

	  t1=yearmax*19
      tdnu=0
	  tdnus=0

	  do i=1,recbumax
	    samok_td(i)=.false.
		td_td(i)=0
		mjsk_td(i)=0
		vis_td(i)=0
	  enddo

      !Hérna byrjar heljarinnar lykkja yfir öll bú
	  do i=1, nuh
		t3=0
    	taltdhop(1)=0
	    taltdhop(2)=0
	    taltdhop(3)=0
	    taltdhop(4)=0

	    taltdhops(1)=0
	    taltdhops(2)=0
	    taltdhops(3)=0
	    taltdhops(4)=0
		tdnu=tdnu+1
		tdnus=tdnus+1
	    do j=1,recbumax
	      samok_td(j)=.false.
		  td_td(j)=0
		  mjsk_td(j)=0
		  vis_td(j)=0
	    enddo
	    k=0
	    if(mod(i,100).eq.0)write(*,*)i,' herds'

      !Fram að þessu bara eitthvert dund, fylla allt af núllum aftur eftir síðasta bú
!hér kemur svo lykkja yfir all mælingar sem pikkar út mælingar frá þessu búi og setur í sérsaka vigra
	    do j=1, gr
		  if(herd(j).eq.herdv(i).and.milkok(j))then   !ef sama bú og ok mæling
            k=k+1
			vis_td(k)=j      !pointer yfir í aðal gagnavigrana
			td_td(k)=td(j)
			mjsk_td(k)=lact(j)
			samok_td(k)=samok(j)
			if(samok(j))t3=t3+1
          endif
        enddo
        dag=(firstyear-1993)*365
!	    write(*,*) 'c'
		t2=k

    !hérna kemur lykkja sem er svipuð og þegar bú-ár hóparnar voru skilgreindir, nema hér fyrir bú-mælidag.
!hérna erum við innan lykkju yfir bú, en hér byrjar svo lykkja yfir tíma
!Til þess að mælidagur sé leifilegur er þurfa að vera 5 mælingar, þvert á mjaltaskeið, í hópnm.
!þar sem bú-mælidagur er slembibreyta í líkaninu er í lagi að hóparnir séu litlir. Þeir nýta líka upplýsingar á milli mjaltaskeiða því það er fylgni á milli í líkaninu
!samt sem áður er hér skilyrði að ekki megi vera aðeins 1 mæling frá mjaltaskeiði í hópnum, þá er sameinað (0 er leift). Held reyndar að það sé della að hafa það inni.

		do l=1,t1        ! tl er einhver nógu há tala til að ná yfir alla mögulega mælidaga innan bús
          do j=1, k       !lykkja yfir mælingar búsins
!		    write(*,*)'f'
            if(td_td(j).lt.dag.and.td_td(j).gt.1.and.samok_td(j))then        !dag hækkar um 20 í hverri l lykkju, svo það tekur alltaf mælingar frá næstu 20 dögum(ef einhverjar)
!			  write(*,*)'s'
              t3=t3-1
		      htds(vis_td(j))=tdnus
!			  write(*,*)'t'
			  taltdhops(mjsk_td(j))=taltdhops(mjsk_td(j))+1
!			  write(*,*)'o'
			  taltdhops(4)=taltdhops(4)+1             !taltdhops(4) er fyrir öll mjaltaskeið, 1,2,3 eru fyrir fyrsta, annað og þriðja
!			  write(*,*)'p'
	        endif
            if(td_td(j).lt.dag.and.td_td(j).gt.1)then
!			  write(*,*)'r'
              t2=t2-1
		      htd(vis_td(j))=tdnu
			  taltdhop(mjsk_td(j))=taltdhop(mjsk_td(j))+1
			  taltdhop(4)=taltdhop(4)+1
			  td_td(j)=0
            endif
          enddo
!	  if(mod(l,1).eq.0)write(*,*)'d', l
      	  ok=.false.
	      if(taltdhop(4).gt.5.and.t2.gt.5)ok=.true.
		  if(taltdhop(4).gt.2.and.taltdhop(1).ne.1.and.taltdhop(2)
     +      .ne.1.and.taltdhop(3).ne.1.and.t2.gt.5)ok=.true.
		  if(ok)then
!		    write(*,*)'nullun'
		    tdnu=tdnu+1
		    taltdhop(1)=0
		    taltdhop(2)=0
		    taltdhop(3)=0
		    taltdhop(4)=0
		  endif
          ok=.false.
	      if(taltdhops(4).gt.5.and.t3.gt.5)ok=.true.
		  if(taltdhops(4).gt.2.and.taltdhops(1).ne.1.and.taltdhops(2)
     +      .ne.1.and.taltdhops(3).ne.1.and.t3.gt.5)ok=.true.
!	      write(*,*)'g'
		  if(ok)then
		    tdnus=tdnu
		    taltdhops(1)=0
		    taltdhops(2)=0
		    taltdhops(3)=0
		    taltdhops(4)=0
          endif
		  dag=dag+20
!		 write(*,*)'e',tdnu,tdnus
        enddo

	  enddo              !bú-lykkjan búin
	  write(*,*)tdnu,' herd-test-days'



  210 format(i15,i15,i15,6x,i6,i1)
  211 format(i15,1x,i15,1x,i15,1x,i4)

      allocate(pedid(pedmax),dam(pedmax),sire(pedmax),ar(pedmax))
	  allocate(kyn(pedmax))
! Read the pedigree file
      open(10,file=filename2)

      write(*,*)'Pedigree file: ',filename2
	  write(*,*)'Maximum number of animals is set to ',pedmax
	  biljon=100000
	  biljon=biljon*1000000
	  k=0
      !lykkja sem les ætternisskrána
	  do
	    k=k+1
	    read(10,210,end=30)pedid(k),dam(k),sire(k),land,kyn(k)
        ar(k)=pedid(k)/biljon
		if(mod(k,50000).eq.0)write(*,*)k,' animals from ped file'
	    if(land.gt.900000)k=k-1          !hendir út gripum með hærra búsnúmer en 9000000
	  enddo
  30  nrped=k-1

	  allocate(recno(nrped),fix(3,nrped))

      close(10)
	  do i=1,nrped
	    recno(i)=0
		fix(1,i)=0
		fix(2,i)=0
		fix(3,i)=0
      enddo

! Change the id numbers in the data file



      write(*,*)'put running numbers to the data'
      j=0
      do i=1,gr                         !lykkja yfir mælingar
	    if(mod(i,500000).eq.0)write(*,*)i
        if(langsok(pedid,k,id(i),ql).and.milkok(i))then         !leitar að gripnum í ætternisskránni
		  id(i)=ql                                              !Skrifar númerið í ætternisskránni sem nýtt ID
		  recno(ql)=recno(ql)+1
		  if(lact(i).eq.1.and.samok(i))fix(1,ql)=hys(i)       !Þette fix dæmi er til að skrifa út fixed hópana í sér skrá sem er svo notuð í öryggisútreikningana.
		  if(lact(i).eq.1.and.samok(i))fix(2,ql)=ca(i)
		  if(lact(i).eq.1.and.samok(i))fix(3,ql)=tm(i)
        else                                         !hoppar yfir í það sem kemur hér á eftir else ef að gripurinn finnst ekki í ætternisskánni
		  if(milkok(i))then
		  do l=1,j
		    if(id(i).eq.extra(l))then        !athugar hvort þegar er búið að búa til aukaeinstakling fyrir þennan grip
              id(i)=2000000+l            !Setur þá sama raðnúmer á
 			  goto 31                               !Hoppar til 31 ef eitthvað fannst
			endif
		  enddo
!          write(*,*)j
          j=j+1
          extra(j)=id(i)        !Extra er með lista af longu ID sem fundust í gögnum en ekki í ættskránni
		  extra_year(j)=by(i)       !Geymi líka ártalið upp á draugaforeldrahópana seinna
          id(i)=2000000+j             !Nýtt ID fyrir týnda gripinn
          endif
	    endif
  31  continue
      enddo
	  aukatal=j         !halda til haga hve margir eru aukalega

	  open(20,file='../dmu_data/tdm1.dat')
	  open(22,file='../dmu_data/tdm2.dat')
	  open(23,file='../dmu_data/tdm3.dat')
!	  open(24,file='tdm_fpf.dat')

  220 format(i8,1x,i7,1x,i7,1x,i4,1x,i4,5(f8.4,1x),3(f7.2,1x),3(1x,i7))
  222 format(i8,1x,i7,1x,i7,1x,i4,1x,i4,5(f8.4,1x),15(f7.2,1x),3(1x,i7)
     +  )
  223 format(i8,1x,i7,1x,i7,1x,i4,1x,i4,5(f8.4,1x),6(f7.2,1x),3(1x,i7))
!  229 format(i8,1x,i7,1x,i7,1x,i4,1x,i4,5(f8.4,1x),18(f7.2,1x),i8)

!writing the data file with new id numbers.
	  write(*,*)j,' animals in data file not in ped file'

      do i=1, gr              !lykkja yfir mælingar til að skrifa út gagnaskrárnar, tdm1.dat os.frv.
	    if(milkok(i))then
		  stdim=real(dim(i))
		  stdim=2*((stdim-5)/300)-1
		  lp1=1.2247*stdim            !Býr til legendre fjölliðurnar fyrir aðhvorfin í DMU.(lp1,lp2,lp3,lp4)
		  lp2=-0.7906+2.3717*stdim*stdim
		  lp3=-2.8062*stdim+4.6771*stdim*stdim*stdim
		  lp4=0.7955-7.9550*stdim*stdim+9.2808*stdim*stdim*stdim*stdim
		  wil=exp(-0.05*dim(i))                   !wilmink
          if(lact(i).eq.1)then            !hver mæling fer í eina línu. Þess vegna kemur -999 (missing) fyrir öll önnur mjaltaskeið en þau sem mælingin er frá
		    myout(1)=my(i)
		    if(samok(i))then
			  fyout(1)=fy(i)*10
		      pyout(1)=py(i)*10
		      scsout(1)=scs(i)
		      ppout(1)=pp(i)
		      fpout(1)=fp(i)
			else
			  fyout(1)=-999
		      pyout(1)=-999
		      scsout(1)=-999
		      ppout(1)=-999
		      fpout(1)=-999
			endif
		    myout(2)=-999
		    fyout(2)=-999
		    pyout(2)=-999
		    scsout(2)=-999
		    ppout(2)=-999
		    fpout(2)=-999
		    myout(3)=-999
		    fyout(3)=-999
		    pyout(3)=-999
		    scsout(3)=-999
		    ppout(3)=-999
		    fpout(3)=-999
          elseif(lact(i).eq.2)then
		    myout(2)=my(i)
			if(samok(i))then
		      fyout(2)=fy(i)*10
		      pyout(2)=py(i)*10
		      scsout(2)=scs(i)
		      ppout(2)=pp(i)
		      fpout(2)=fp(i)
			else
			  fyout(2)=-999
		      pyout(2)=-999
		      scsout(2)=-999
		      ppout(2)=-999
		      fpout(2)=-999
			endif
		    myout(1)=-999
		    fyout(1)=-999
		    pyout(1)=-999
		    scsout(1)=-999
		    ppout(1)=-999
		    fpout(1)=-999
		    myout(3)=-999
		    fyout(3)=-999
		    pyout(3)=-999
		    scsout(3)=-999
		    ppout(3)=-999
		    fpout(3)=-999
          elseif(lact(i).eq.3)then
		    myout(3)=my(i)
			if(samok(i))then
		      fyout(3)=fy(i)*10
		      pyout(3)=py(i)*10
		      scsout(3)=scs(i)
		      ppout(3)=pp(i)
		      fpout(3)=fp(i)
			else
			  fyout(3)=-999
		      pyout(3)=-999
		      scsout(3)=-999
		      ppout(3)=-999
		      fpout(3)=-999
			endif
		    myout(2)=-999
		    fyout(2)=-999
		    pyout(2)=-999
		    scsout(2)=-999
		    ppout(2)=-999
		    fpout(2)=-999
		    myout(1)=-999
		    fyout(1)=-999
		    pyout(1)=-999
		    scsout(1)=-999
		    ppout(1)=-999
		    fpout(1)=-999
          endif
		  write(20,220)id(i),hy(i),htd(i),ca(i),tm(i),    !skrifa tdm1.dat ID, bú-ár, bú-mælidagur, burðaraldur, mánuður??, lp1,lp2,lp3,lp4,wilmink,mjólk1,mjólk2, mjólk3, dagar frá burði, bú, burðarmán.
     +       lp1,lp2,lp3,lp4,wil,
     +       myout(1),myout(2),myout(3)
     +       ,dim(i),herd(i),cm(i)
	     if(samok(i))write(22,222)id(i),hys(i),htds(i),ca(i),tm(i),       !Skrifa tdm2.dat. Það er með öðru en mjólk
     +       lp1,lp2,lp3,lp4,wil,
     +       fyout(1),fyout(2),fyout(3),
     +       pyout(1),pyout(2),pyout(3),fpout(1),fpout(2),fpout(3)
     +       ,ppout(1),ppout(2),ppout(3),scsout(1),scsout(2),scsout(3),
     +        dim(i),herd(i),cm(i)
         if(samok(i))write(23,223)id(i),hys(i),htds(i),ca(i),tm(i),     !tmd3.dat er bara með prósentueiginleikunum
     +       lp1,lp2,lp3,lp4,wil,
     +       fpout(1),fpout(2),fpout(3),
     +       ppout(1),ppout(2),ppout(3),
     +       dim(i),herd(i),cm(i)
!	     if(samok(i))write(24,229)id(i),hys(i),htds(i),ca(i),tm(i),
!     +       lp1,lp2,lp3,lp4,wil,
!     +       myout(1),myout(2),myout(3),
!     +       fyout(1),fyout(2),fyout(3),
!     +       pyout(1),pyout(2),pyout(3),fpout(1),fpout(2),fpout(3)
!     +       ,ppout(1),ppout(2),ppout(3),scsout(1),scsout(2),scsout(3)
!     +       ,herd(i)
		endif
      enddo

      write(*,*)'Data written to tdm.dat'

	  close(20)
	  close(22)
	  close(23)
!	  close(24)

                     !Hérna er eiginlega búið að vinna allt með mælingarnar, restin snýst um ætternisskrána.
!	  goto 999

  221 format(4(i8,1x))
  224 format(i15,1x,i8,1x,i2,i3,3(1x,i5),i2)

      open(21,file='../dmu_data/dmu_ped.txt')
	  open(24,file='../dmu_data/radnrkodi')

	  allocate(draugur(ruglpedmax),auka(ruglpedmax),
     +   draug_ar(ruglpedmax),raddraug(ruglpedmax))

c....her byrjar lykkja til ad setja radnr i stad longu
c....einstaklingsnumeranna i aetternisskra. adalpudrid
c... er i foreldra sem eru ekki i aettskranni.
      taldraug=0                !Draugur á hér við gripi sem eru með mælingu en ekki í ætternisskránni. Sumsé ekki draughópa
      do i=1, ruglpedmax
        draugur(i)=0
      enddo
	  do i=1,40
	    nphg(1,i)=0
		nphg(2,i)=0
	  enddo
      write(*,*)'000'

	  stdanno=0
      damtal=0
      siretal=0
      arsums=0
      arsumd=0
!	  faertal=0
      do i=1, nrped             !Lykkja yfir gripi í ætternisskránni
        einst(1)=pedid(i)
        einst(2)=dam(i)
        einst(3)=sire(i)
        einstut(1)=i

        if(mod(i,50000).eq.0)write(*,*)i, 'animals'

        if(langsok(pedid,nrped,einst(2),ql))then       !leita að móður gripsins í ætternissrkánni

          if(ar(ql).eq.9999)ar(ql)=ar(i)-5      !Ef ekkert ártal á foreldri er það sett 5 árum fyrir fæðingu afkvæmis
		  if(kyn(ql).eq.1)then                        !Ef móðir er naut er henni eytt
		    write(*,*)'bull ',pedid(ql),
     +       'is registered dam of ' ,pedid(i)
	        einst(2)=0
	      endif
          bil=ar(i)-ar(ql)                      !aldursmismunum móður og afkvæmis
		  if(bil.lt.1.or.bil.gt.30)einst(2)=0  !Ef meira en 30 ár er móður eytt. Einnig ef hún er jafn gömul eða yngri en afkvæmið
          arsumd=arsumd+bil
          damtal=damtal+1
          einstut(2)=ql               !ef allt var ok fær móðirin radnúmer

!		  if(bil.lt.1.or.bil.gt.30)goto 555
!		else if(dam(i).eq.0)then
        else if(einst(2).gt.0)then                 !ef að gripurinn er með móður en hún er ekki í ætternisskránni

!          write(*,*)'modir ',einst(2),'finnst ekki'
          j=1
          logdraug=.true.
          do while(logdraug)                       !athuga hvort þessi gripur hefur áður komið fram sem foreldri

             if(einst(2).eq.draugur(j))then
              einstut(2)=raddraug(j)
!             write(*,*)'tvidraugur ',einst(2)
              logdraug=.false.

            else if(draugur(j).lt.1)then                !Ef þetta er í fyrsta skipti sem þessi gripur kemur sem foreldri (og er ekki í ætternisskránni sjálfur) fær hann nýtt ID hærra en miljón
              taldraug=taldraug+1
              raddraug(j)=j+1000000
              einstut(2)=raddraug(j)
              draugur(j)=dam(i)
              logdraug=.false.
              draug_ar(j)=ar(i)-5

            endif
            j=j+1

            if(taldraug.gt.(ruglpedmax-1))then
              write(*,*)'Alvarleg villa fullur draugur'      ! ef og margir foreldrar eru í ruglinu kemur villa og forritið stoppar
              goto 999
            endif
          enddo
        endif
		if(einst(2).eq.0)then                           !Ef að móður vantar, eða henni hefur verið eytt því hún var of gömul eða eitthvað
		  if(binsok(dphg,ndphg,ar(i),pl))continue     !finnur fæðingarárið á eftir í draughópalistanum sem var lesinn úr stýriskránni í upphafi
		  pl=pl-1
!		  if(ar(i).gt.dphg(ndphg))pl=pl+1
		  einstut(2)=-pl                              !fær ID sem er mínustala
		  nphg(1,pl)=nphg(1,pl)+1                      !talið í draugforeldrahópnum
        endif

        if(langsok(pedid,nrped,einst(3),ql))then  !Hérna byrjar sama ferlið með föður
		  if(kyn(ql).eq.2)then
		    write(*,*)'cow ',pedid(ql),
     +       'is registered sire of ',pedid(i)
          einst(3)=0
	      endif
          if(ar(ql).eq.9999)ar(ql)=ar(i)-5
		  bil=ar(i)-ar(ql)
          if(bil.lt.1.or.bil.gt.30)einst(3)=0
          arsums=arsums+bil
          siretal=siretal+1
          einstut(3)=ql
!		  if(bil.lt.1.or.bil.gt.30) goto 556
!		else if(sire(i).eq.0)then

        else if(einst(3).gt.0)then
!         write(*,*)'fadir ', einst(3), 'finnst ekki'
          j=1
          logdraug=.true.
          do while(logdraug)
            if(einst(3).eq.draugur(j))then
              einstut(3)=raddraug(j)
!              write(*,*)'tvidraugur ',einst(3)
              logdraug=.false.
            else if(draugur(j).lt.1)then
              taldraug=taldraug+1
              raddraug(j)=j+1000000
              einstut(3)=raddraug(j)
              draugur(j)=sire(i)
              logdraug=.false.
              draug_ar(j)=ar(i)-5
            endif
            j=j+1
            if(taldraug.gt.(ruglpedmax-1))then
              write(*,*)'Alvarleg villa fullur draugr'
              goto 999
            endif
          enddo
        endif
		if(einst(3).eq.0)then
		  if(binsok(sphg,nsphg,ar(i),pl))continue
		  pl=pl-1
!		  if(ar(i).gt.sphg(nsphg))pl=pl+1
		  einstut(3)=-pl-1000
		  nphg(2,pl)=nphg(2,pl)+1
        endif

!Þa er hægt að skrifa út nýja ætternisskrá með radnúmerum
        write(21,221)einstut(1),einstut(3),einstut(2),ar(i)
		if(ar(i).eq.stdyear.and.kyn(i).eq.2.and.recno(i).gt.1)then
		  stdan=1
		  stdanno=stdanno+1
		else
		  stdan=0
		endif
		write(24,224)einst(1),einstut(1),stdan,recno(i),      !skrifar radnrkodi með upplýsingum
     +    (fix(j,i),j=1,3),kyn(i)

      enddo
	  write(*,*) nrped,' gripir skrifadir med radnr'
!	  close(27)
	  stdan=0

     !Þá er bara að bæta vesenisforeldrum við ættskrána
c...  Næst er að bæta foreldrum sem eru ekki sjalfir i aettskránni
c     við og setja þeirra foreldra í rétta hópa.
      do i=1, taldraug
	    if(binsok(dphg,ndphg,draug_ar(i),pl))continue
		pl=pl-1
!		if(draug_ar(i).gt.dphg(ndphg))pl=pl+1
		einstut(2)=-pl
		nphg(1,pl)=nphg(1,pl)+1
	    if(binsok(sphg,nsphg,draug_ar(i),pl))continue
		pl=pl-1
!		if(draug_ar(i).gt.sphg(nsphg))pl=pl+1
		einstut(3)=-pl-1000
		nphg(2,pl)=nphg(2,pl)+1

		utdraugur=i+1000000
        write(21,221)utdraugur,einstut(3),einstut(2),draug_ar(i)       !draug_ar eru ætluð fæðingarár þessara tilbúnu gripa
		write(24,224)draugur(i),utdraugur,stdan
      enddo
c...  Skrifar smá upplýsingar á skjáinn
      write(*,*)taldraug,' sem vantadi baett vid'

c...  Næst er að bæta i ut aettskrana þeim sem eru með mælingar
c...  en vantaði í ættskrána.
      do i=1, aukatal
	    if(binsok(dphg,ndphg,extra_year(i),pl))continue
		pl=pl-1
!		if(extra_year(i).gt.dphg(ndphg))pl=pl+1
		einstut(2)=-pl
		nphg(1,pl)=nphg(1,pl)+1
	    if(binsok(sphg,nsphg,extra_year(i),pl))continue
		pl=pl-1
!		if(extra_year(i).gt.sphg(nsphg))pl=pl+1
		einstut(3)=-pl-1000
		nphg(2,pl)=nphg(2,pl)+1

		aukaut=i+2000000
        write(21,221)aukaut,einstut(3),einstut(2),extra_year(i)
		write(24,224)extra(i),aukaut,stdan
      enddo
	  write(*,*)aukatal,' ur gagnaska sem vantadi i aettskra skrifad'

      allsped=nrped+taldraug+aukatal
      write(*,*)'samtals ', allsped
      close(21)
	  close(24)
      meddam=real(arsumd)/real(damtal)                  !Reikna kynslóðabil til gamans
      medsire=real(arsums)/real(siretal)
      write(*,*)'kynslodabil maedur: ', meddam, ' fjoldi ',damtal
      write(*,*)'kynslodabil fedur: ', medsire, ' fjoldi ',siretal

	  open(25,file='uppl.txt')                      !skrá með helstu tölum
	  write(25,*)nrped,'  animals in radnrkodi form the pedigree file'
	  write(25,*)stdanno,' cows born in ',stdyear,' for standardizing'
	  write(25,*)allsped,'  lines in radnrkodi'
	  write(25,*)hysnr(1)
	  write(25,*)22
	  write(25,*)12

	  close(25)
	  open(31,file='../dmu_data/damphgout.txt')          !Skrá með fjölda í draugahópum til upplýsinga
      open(32,file='../dmu_data/sirephgout.txt')
	  do i=1,ndphg
	    write(31,*)dphg(i),nphg(1,i)
	  enddo
	  do i=1,nsphg
	    write(32,*)sphg(i),nphg(2,i)
	  enddo
	  close(31)
      close(32)
	  write(*,*)'Check numbers in PHG in damphgout.txt',
     +     ' and sirephgout.txt'

  999 stop          !sjálft forritið búið
      end

	 !Hraðvirk leitarfunktíon frá Ágústi. Tvær útgáfur, önnur fyrir löng ID, hin fyrir tölur.
	  logical function binsok(a,n,t,pl)
      integer n,pl
      integer min,max,mitt
      integer a(n),t
      min=1
      max=n
   10 if(min.lt.max)then
         mitt=(min+max)/2
      if(t.le.a(mitt))then
         max=mitt
      else
        min=mitt+1
      endif
        goto 10
      endif

      if(t.eq.a(min))then
        binsok=.true.
        pl=min
      else
        binsok=.false.
		pl=min
      endif

      return
      end

      logical function langsok(a,n,t,ql)
      integer n,ql
      integer min,max,mitt
      integer*8 a(n),t
      min=1
      max=n
   11 if(min.lt.max)then
         mitt=(min+max)/2
      if(t.le.a(mitt))then
         max=mitt
      else
        min=mitt+1
      endif
        goto 11
      endif

      if(t.eq.a(min))then
        langsok=.true.
        ql=min
      else
        langsok=.false.
		ql=min
      endif

      return
      end
