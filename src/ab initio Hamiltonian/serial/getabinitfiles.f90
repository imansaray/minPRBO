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
          
