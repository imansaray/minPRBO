        subroutine header(filenum, fermie, maxnpw, maxband, bantot, nsppol,  &
                     cplex, nspden, nspinor, nkpt, ngfft,rprimd, kptrlatt)
	implicit none
        integer :: filenum, maxnpw, maxband
        
        character*8 :: codvsn
	integer :: headform,fform
	integer :: bantot,cplex,date,icoulomb,intxc,ixc,kptopt,mband,natom,ngfft(3),&
  	           nkpt,npsp,nshiftk_orig,nshiftk,nspden,nspinor,nsppol,nsym,ntypat,occopt,&
                   pawcpxocc,pertcase,usepaw,usewvl,ipsp
	double precision :: acell(3),cellcharge,ecut,ecutdg,ecutsm,ecut_eff,etotal,&
                   fermie,nelect,qptn(3),residm,rprimd(3,3),shiftk(3),shiftk_orig(3),&
                   stmbias,tphysel,tsmear
	integer, allocatable :: istwfk(:),nband(:),&
  	           npwarr(:),symafm(:),symrel(:,:,:),so_psp(:),typat(:),nrhoijsel(:),&
  	           rhoijselect(:,:)
  	integer :: kptrlatt(3,3),kptrlatt_orig(3,3)
	double precision, allocatable :: kpt(:,:),occ(:),tnons(:,:),wtk(:),&
	           rhoij(:,:),xred(:,:),znucltypat(:),amu(:)
	character*132 :: title
	character*32 :: md5_pseudos
	double precision :: znuclpsp,zionpsp
	integer :: pspso,pspdat,pspcod,pspxc,lmn_size,temp
        
        write(6,*) "======= Reading the WFK file ======="
        read(unit=filenum) codvsn,headform,fform
        write(6,*) "codvsn=",codvsn,"headform=",headform,"fform=",fform

        read(unit=filenum) bantot,date,intxc,ixc,natom,ngfft(1:3),&
              nkpt,nspden,nspinor,nsppol,nsym,npsp,ntypat,occopt,pertcase,usepaw,&
              ecut,ecutdg,ecutsm,ecut_eff,qptn(1:3),rprimd(1:3,1:3),stmbias,&
              tphysel,tsmear,usewvl,nshiftk_orig,nshiftk,mband

         write(6,*) bantot,date,intxc,ixc,natom,ngfft(1:3)
         write(6,*) " "
         write(6,*) nkpt,nspden,nspinor,nsppol,nsym,npsp,ntypat,occopt,pertcase,usepaw
         write(6,*) " "
         write(6,*) ecut,ecutdg,ecutsm,ecut_eff
         write(6,*) ""
         write(6,*) qptn(1:3),rprimd(1:3,1:3)
         write(6,*) " "
         write(6,*) stmbias
         write(6,*) " "
         write(6,*) tphysel,tsmear,usewvl,nshiftk_orig,nshiftk,mband
         write(6,*) " usepaw = ",usepaw
         write(6,*) " Done printing the read variables"


!         allocate( rhoijselect(50,nspden))
!         write(6,*) " now deallocating nband"
!         deallocate( rhoijselect)
         allocate( istwfk(nkpt),nband(nkpt*nsppol),&
              npwarr(nkpt),symafm(nsym),symrel(3,3,nsym),&
              kpt(3,nkpt),occ(bantot),tnons(3,nsym),wtk(nkpt))
         allocate(xred(3,natom),znucltypat(ntypat),amu(ntypat),typat(natom),&
                so_psp(npsp),nrhoijsel(nspden)) !,rhoijselect(50,nspden),rhoij(100,nspden))

!         read(unit=filenum) istwfk(1:nkpt),nband(1:nkpt*nsppol)
         read(unit=filenum) istwfk(1:nkpt),nband(1:nkpt*nsppol),&
              npwarr(1:nkpt),so_psp(1:npsp),symafm(1:nsym),symrel(1:3,1:3,1:nsym),&
              typat(1:natom),kpt(1:3,1:nkpt),occ(1:bantot),tnons(1:3,1:nsym),&
              znucltypat(1:ntypat),wtk(1:nkpt)

         read(unit=filenum) residm,xred(1:3,1:natom),etotal,fermie,amu(1:ntypat)

         read(unit=filenum) kptopt,pawcpxocc,nelect,cellcharge,icoulomb,&
              kptrlatt(3,3),kptrlatt_orig(3,3),shiftk_orig(3),shiftk(3)

         do ipsp=1,npsp
    ! (npsp lines, 1 for each pseudo; npsp=ntypat, except if alchemical pseudo-atoms)
            read(unit=filenum) title,znuclpsp,zionpsp,pspso,&
                 pspdat,pspcod,pspxc,lmn_size,md5_pseudos
         enddo

!    !(in case of usepaw==1, there are some additional records)
!        if (usepaw==1)then
!           read(unit=filenum) (pawrhoij(iatom)%nrhoijsel(1:nspden),iatom=1,natom), cplex, nspden
!           read(unit=filenum) ((pawrhoij(iatom)%rhoijselect(1:nrhoijsel(ispden),ispden),ispden=1,nspden),iatom=1,natom),&
!                     ((pawrhoij(iatom)%rhoijp(1:cplex*nrhoijsel(ispden),ispden),ispden=1,nspden),iatom=1,natom)
!        endif
         
         maxnpw = maxval(npwarr)
         maxband = maxval(nband)

         deallocate( istwfk,nband,npwarr,symafm,symrel,kpt,occ,tnons,wtk)
         deallocate( xred,znucltypat,amu,typat,so_psp,nrhoijsel)
         write(6,*) " "
         write(6,*) "----- Completed the creation of all arrays ---- "         
         return 
	 end
	 
	 
	 subroutine grabwf (filenum1,filenum2,filenum3,filenum4,bantot,nspinor,nsppol,nkpt,maxnpw,maxband)
	 
	 integer ::filenum1,filenum2,filenum3,filenum4,bantot,maxnpw,maxband,&
                 iband,nband,ii,ikpt,nkpt,npw,nspinor,nsppol,isppol,index,&
                 inputstat
         integer(kind=16):: temp
	 double precision :: eigen(bantot),occ(bantot),Ha2ryd
         integer, allocatable :: kg(:,:),kgtrans(:,:)
         double precision, allocatable  :: cg(:),temp_cg(:,:)

         ! Conversion factor from Hartree to Rydberg
          Ha2ryd = 2.0d0
         ! calculating the size of cg array , a 1-D array
          temp = 0
!          temp = 2*maxnpw*nspinor*bantot*nkpt*nsppol
          temp = 2*maxnpw*nspinor
          temp = temp*bantot*nkpt
          temp = temp*nsppol
          write(6,*) "maxnpw=",maxnpw,"nspinor=",nspinor,"bantot=",bantot
          write(6,*) "nkpt=",nkpt,"nsppol=",nsppol,"temp =",temp,kind(temp)
          allocate (cg(temp))
          cg = 0.0
          eigen = 0.0
          occ = 0.0
          
          write(6,*) "------ Now reading in the WFK data ------"
          
          bantot=0                                    
          index=0                                     
          do isppol=1,nsppol
           do ikpt=1,nkpt
            write(6,*) " "
            write(6,*) "isppol =",isppol,"ikpt=",ikpt
            read(filenum1) npw,nspinor,nband 
            write(6,*) "npw=",npw,"nspinor=",nspinor,"nband=",nband
            write(filenum4) npw,nspinor,nband 
            write(6,*) " "
            allocate (kg(3,npw))
            allocate (kgtrans(npw,3))
            kg = 0  
            kgtrans = 0

            write(6,*) "==== writing kg ======"
!            read(filenum1) kg(1:3,1:npw)
            read(filenum1,IOSTAT=inputstat) kg(1:3,1:npw)
            kgtrans = transpose(kg) ! transpose, since the serial code has kg(npw,3)
            write(6,*) " shape kgtrans =",shape(kgtrans)
            write(filenum4) kgtrans 
!            write(6,*) "IOSTAT = ",inputstat
!            write(6,*) "kg = ",kg
            write(6,*) " "           
            read(filenum1,IOSTAT=inputstat) eigen(1+bantot:nband+bantot), &
                           occ(1+bantot:nband+bantot)         
!            write(6,*) "IOSTAT = ",inputstat
                        
           ! writing the eigen values and occ data for each k pt             
            write (filenum2,*) Ha2ryd*eigen(1+bantot:nband+bantot)
            write (filenum3) occ(1+bantot:nband+bantot)
            do iband=1,nband
             read(filenum1,IOSTAT=inputstat) (cg(ii+index),ii=1,2*npw*nspinor)   
             write(6,*) "IOSTAT = ",inputstat
             write (filenum4) (cg(ii+index),ii=1,2*npw*nspinor)
            enddo                         
             
            ! deallocate temp_cg and cg_cmplx
            deallocate(kg, kgtrans)
            write(6,*) " "
            write(6,*) "--- Before update ---"
            write(6,*) "bantot =",bantot,"index =",index            
            bantot=bantot+nband
            index=index+2*npw*nspinor*nband
            write(6,*) " "
            write(6,*) "--- After update ---"
            write(6,*) "bantot =",bantot,"index =",index            
           enddo
          enddo
          
          write(6,*) " Completed reading all the WFK data" 
         end
         
         
         subroutine grabden (filenum1,filenum2,cplex,nspden,ngfft)	 
	 integer ::filenum1,filenum2,cplex,ispden,nspden,ngfft(3),ir,temp,inputstat
         double precision, allocatable  :: rhor(:) 

         ! calculating the size of rhor array , a 1-D array
          temp = cplex*ngfft(1)*ngfft(2)*ngfft(3)
          write(6,*) "cplex=",cplex,"nspden=",nspden
          write(6,*) " "
          write(6,*) "ngfft(1)=",ngfft(1),"ngfft(2)=",ngfft(2),"ngfft(3)=",ngfft(3)
          allocate (rhor(temp))
          
          rhor = 0.0
          
          write(6,*) "------ Now reading in the DEN data ------"
          do ispden=1,nspden
             read(filenum1,IOSTAT=inputstat) (rhor(ir),ir=1,cplex*ngfft(1)*ngfft(2)*ngfft(3))
             write(6,*) "IOSTAT = ",inputstat
             write(filenum2,*) (rhor(ir),ir=1,cplex*ngfft(1)*ngfft(2)*ngfft(3))
          enddo
          
          ! deallocate the rhor array
          deallocate(rhor)
                    
          write(6,*) " Completed reading all the DEN data" 
         end

!  Quantum Espresso files

         subroutine read_a_wfc(ibnd, iuni,iuni2, ik, xk, nbnd, ispin, npol, gamma_only, ngw, igwx,&
                                b1,b2,b3 )
         use iso_fortran_env, ONLY: DP=> REAL64
         implicit none 
!         character (len=*), intent(in)  :: filename 
         integer, intent(in)  :: ibnd
!         complex(DP), intent(out), allocatable :: ug(:)
         complex( kind = kind(1.0d0) ),allocatable :: ug(:)
         real(dp), intent(out) :: xk(3) 
         integer, intent(out)  :: ik, nbnd, ispin, npol, ngw, igwx
         integer  :: dummy_int   
         logical  :: gamma_only
         integer, allocatable :: mill(:,:),mill_trans(:,:) 
         ! 
!         integer :: iuni = 1111, i 
         integer :: iuni,iuni2, i,j 
         real(dp),allocatable :: ug_re(:),ug_im(:),ug_long(:)
         real(dp) :: scalef
         real(dp) :: b1(3), b2(3), b3(3), dummy_real 
!         open(UNIT = iuni, FILE = trim('filename'), FORM = 'unformatted', status = 'old')

!         write(6,*) " "
!         write(6,*) "iband =",ibnd,"iuni =",iuni
!         write(6,*) " "
         read( iuni) ik, xk, ispin, gamma_only, scalef
         read (iuni) igwx,ngw, npol, nbnd
         read (iuni) b1, b2, b3 
         write(iuni2) ngw,npol,nbnd
!         write(6,*) "ik,  ispin,   gamma_only,        scalef"
!         write(6,*) ik,  ispin, gamma_only, scalef
!         write(6,*) " "
!         write(6,*) " kpt (a.u) =",xk
!         write(6,*) " "
         write(6,*) " ngw,     igwx,      npol,      nbnd"
         write(6,*) ngw, igwx, npol, nbnd
         write(6,*) " "
!         write (6,*)" b1,  (a.u)    " 
!         write (6,*) b1
!         write(6,*) " "
!         write (6,*)" b2,  (a.u)    " 
!         write (6,*) b2
!         write(6,*) " "
!         write (6,*)" b3,  (a.u)    " 
!         write (6,*) b3
!         write(6,*) " "

    ! avoid reading miller indices of G vectors below E_cut for this kpoint 
    ! if needed allocate  an integer array of dims (1:3,1:igwx) 
    !
         allocate(mill(3,ngw),mill_trans(ngw,3))
         mill = 0
         mill_trans = 0
         read (iuni) mill(1:3,1:ngw)
         mill_trans = transpose(mill)
         write(iuni2) mill_trans

!         write(6,*) " Writing out the miller indices for band ",ibnd,"kpt #",ik
!         write(6,*) " "
!         do i=1,ngw
!             write(6,*) "g = ",i
!             write(6,*) mill(:,i)
!             write(6,*)" "
!         enddo
!    !  
         allocate(ug(npol*ngw),ug_re(npol*ngw),ug_im(npol*ngw),ug_long(2*npol*ngw))
!         allocate(mill(3,igwx))
         if ( ibnd > nbnd) then 
            print '("looking for band nr. ",I7," but there are only ",I7," bands in the file")',ibnd, nbnd
            stop
         end if 
         do i = 1, nbnd
            ug = 0.0
            ug_re = 0.0
            ug_im = 0.0
            ug_long = 0.0 
!            if ( i == ibnd) then 
!            write(6,*) " Writing Fourier coeffs for band ",i,"kpt #",ik,"ngw=",ngw
            read(iuni) ug(1:npol*ngw)
            ug_re = real(ug)
            ug_im = aimag(ug) 
!            write(6,*) " "
!            do j=1,ngw
!               write(6,*) "g = ",j
!               write(6,*) "u(g) = ",ug(j),"u(g)_re =",ug_re(j),"u(g)_im=",ug_im(j)
!               write(6,*)" "
!            enddo
       
!       ! putting it in ABINIT form so I can keep all other codes the same
            ug_long(1:npol*ngw) = ug_re
            ug_long((npol*ngw)+1:2*npol*ngw) = ug_im
            write(iuni2) ug_long
         end do 
         deallocate(ug,ug_re,ug_im,ug_long,mill,mill_trans)
         end subroutine read_a_wfc 

         subroutine grabdenqe (filenum1,filenum2,gamma_only,ngm_g,nspin)	 
	 integer ::filenum1,filenum2,ngm_g,ispin,nspin,ir,temp,inputstat
         double precision, allocatable  :: rhog(:,:) 
         integer  :: dummy_int   
         logical  :: gamma_only
         double precision :: b1(3), b2(3), b3(3), dummy_real 

          write(6,*) "------ Now reading in the DEN data ------"
        ! Reading the data
          read(filenum1,IOSTAT=inputstat)gamma_only,ngm_g,nspin
          read(filenum1,IOSTAT=inputstat)b1,b2,b3
          read(filenum1,IOSTAT=inputstat)dummy_int
          write(6,*) "gamma_only=",gamma_only,"ngm_g=",ngm_g,"nspin=",nspin
          write(6,*) " "
          allocate (rhog(ngm_g,nspin))
          
          rhog = 0.0
          
          do ispin=1,nspin
             read(filenum1,IOSTAT=inputstat)rhog
             write(6,*) "IOSTAT = ",inputstat
             write(filenum2,*) (rhog(ir,ispin),ir=1,ngm_g)
          enddo
          
          ! deallocate the rhor array
          deallocate(rhog)
                    
          write(6,*) " Completed reading all the DEN data" 
         end
