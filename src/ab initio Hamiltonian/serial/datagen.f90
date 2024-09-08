! This program takes in the Quantum Espresso wavefunction files (wfc#k.dat), density, and potential files,
! and then generates the .dat files for the wave function coeffs, uofg in 'u1.dat',
! the eigenvalues for each k point, in 'enkfile', and the occupation for all bands for a 
! given k point 'occ.dat'. It also generates the gvecs input file following Eric's original
! format written in INPUT.doc.
! It will also generate the transition matrix elements 'tmels' using Eric's original format,
! with the polarization, and without the polarization.
! hlleps are the limits for calculating the model dielectric fxn with HLL (decut file)
! OCEAN will be used for now to carry out the QE runs and obtain the wfc files for each k point,
! the charge density file, and the generaed txt file of eigen values to make the enkfile.

       program datagen
       use iso_fortran_env
       use stdlib_sorting, only: sort_index, int_size
       implicit none

       integer :: nkpt,nspinor,nsppol,isppol,nspden,cplex,ngfft(3),kptrlatt(3,3)
       integer :: iband,ikpt,maxnpw,maxband,bantot,ik,npw,nband,ig,i,j,ii,jj,itemp
       logical :: gamma_only
       integer, parameter :: u1dat = 10, enkdat = 11, occdat = 12, qewfc = 13, &
                             qeden = 14, rhoreal = 15, efermi = 16, avecbohr = 17, &
                             kmesh = 18, xmesh = 19, qunitofb = 20, brange = 21, &
                             bvecinvbohr = 22, kpoints = 23, hlleps = 24,e0 = 25, gwipt = 26,&
                             gvecs = 27, xsteps = 28,spect = 29,rhorecp = 30,valenkraw = 31,&
                             conenkraw=32,dftgap=33
       real( kind = kind(1.0d0)), parameter :: eryd = 13.6057d0, a0 = 0.529177d0,decut = 10.0 
       real( kind = kind(1.0d0)), parameter :: fermi_energy = 0.45525071572517 ! from OCEAN,in Ryd
       double precision :: fermie, rprimd(3,3), tempvec(3,3), gprimd(3,3), templen, gvec_length,&
                           ldagap
       integer, parameter :: blftl = 1, blfth = 4,brgtl = 5, brgth = 24, nkx = 16, &
                             nky = 16, nkz = 16, N_I = 2
       double precision, parameter :: diemac = 11.4, gwgap = 1.11  !, &
!                                       ldagap = 1.11 !0.577503885730103 !0.187380739*eryd! 1.11eV --> Ryd = 0.0815834302, DFT gap(HOMO-LUMO) from
                                                       ! OCEAN,0.577503885730103, in eV
                                                       !Si expt gapstatic/macroscopic dielectric
       double precision, allocatable :: zr(:,:),zi(:,:),zrlin(:),zilin(:),zlin(:),enk(:),venkdummy(:),&
                                        cenkdummy(:)
       real( kind = kind(1.0d0) ) :: xk(3),b1(3),b2(3),b3(3)
       complex( kind = kind(1.0d0) ),allocatable :: evc(:)
       character(len = 25 ) :: wfcname
       character(len = 10) :: swfc ! string for wfc index

       ! Sorting gvectors
       integer(int_size), allocatable :: indx(:) 
!       real(dp),allocatable :: gveclen(:)
       integer, allocatable :: gvec(:,:),sorted_gvec(:,:)
       double precision, allocatable :: gveclen(:)

       ! Creating a unique array from an array that has repetitions
       logical, allocatable :: mask(:)
       double precision, allocatable :: uniq_gveclen(:)
       integer, allocatable :: indx_vector(:),cnt_uniq_gveclen(:)
       integer :: gvec_cnt,dummy_int,vdum1,vdum2,cdum1,cdum2
       
       ! flags to control the different calculations in BSE, saved in niter.h
       integer, parameter :: niter = 100, bflag = 0, lflag = 0, backf = 0, aldaf = 0, &
                             qpflg = 0, bwflg = 1, bande = 1,ixmesh = 8,nonherm = 1
       ! Parameters for optical spectra calc
       double precision, parameter :: gam = 0.1,elow = 0.0, ehigh = 20.0 ! all in eV
       integer, parameter :: ne = 100
       double precision :: term

       ! default parameter, cplex, since read only with usepaw==1
       cplex = 1

        ! open the files to read the necessary data
       ngfft(1) = 48
       ngfft(2) = 48
       ngfft(3) = 48

!       open( unit=qewfc, file='wfc1.dat', form='unformatted', status='old')
!       open( unit=qeden, file='charge-density.dat', form='unformatted', status='old')
       open( unit=kpoints, file='nkpt.ipt', form='formatted', status='unknown')
!       rewind qewfc
!       rewind qeden
       rewind kpoints

        ! open the files to store the necessary data
       open( unit=u1dat, file='u1.dat', form='unformatted', status='unknown')
       open( unit=valenkraw, file='QE_EIGS.txt', form='formatted', status='unknown')
       open( unit=conenkraw, file='QE_EIGS_shift.txt', form='formatted', status='unknown')
       open( unit=enkdat, file='enkfile', form='formatted', status='unknown')
!       open( unit=occdat, file='occ.dat', form='unformatted', status='unknown')
!       open( unit=rhoreal, file='rhoofr', form='formatted', status='unknown')
!       open( unit=rhorecp, file='rhoofg', form='formatted', status='unknown')
       open( unit=efermi, file='efermiinrydberg.ipt', form='formatted', status='unknown')
       open( unit=dftgap, file='gap', form='formatted', status='unknown')
       rewind u1dat
       rewind valenkraw
       rewind conenkraw
       rewind dftgap
       rewind efermi

       ! formatted files for input to other programs
!       open( unit=efermi, file='efermiinrydberg.ipt', form='formatted', status='unknown')
       open( unit=avecbohr, file='avecsinbohr.ipt', form='formatted', status='unknown')
       open( unit=bvecinvbohr, file='bvecs', form='formatted', status='unknown')
!       open( unit=kpoints, file='nkpt.ipt', form='formatted', status='unknown'), for QE need to count this from wfc files
       open( unit=kmesh, file='kmesh.ipt', form='formatted', status='unknown')
       open( unit=xmesh, file='xmesh.ipt', form='formatted', status='unknown')
!       open( unit=qunitofb, file='qinunitsofbvectors.ipt', form='formatted', status='unknown')
       open( unit=brange, file='brange.ipt', form='formatted', status='unknown')
       open( unit=hlleps, file='decut', form='formatted', status='unknown')
       open( unit=e0, file='epsilon', form='formatted', status='unknown')
       open( unit=gwipt, file='gwipt', form='formatted', status='unknown')
       open( unit=gvecs, file='gvectors', form='formatted', status='unknown')
       open( unit=xsteps, file='stepper.h', form='formatted', status='unknown')
       open( unit=spect, file='spect.ipt', form='formatted', status='unknown')


       ! input flags for different parts of calc : niter.h
       open( unit=99, file='niter.h', form='formatted', status='unknown' )
       rewind 99
       call iputval( niter, 'niter' )
       call iputval( bflag, 'bflag' )
       call iputval( lflag, 'lflag' )
       call iputval( backf, 'backf' )
       call iputval( aldaf, 'aldaf  ' )
       call iputval( qpflg, 'qpflg' )
       call iputval( bwflg, 'bwflg' )
       call iputval( bande, 'bande' )
       call iputval( nonherm, 'nonherm' )

       close( unit = 99 )

!      Reading in the fermi energy and DFT bandgap       
       read(efermi,*) fermie
       read(dftgap,*) ldagap

       write(6,*) "---- Starting to import the Quantum Espresso density, wfc# files and generate data ----"

       ! Reading in the wavefunction coeffs, eigenenergies, occupation data, and write to file
       read(kpoints,*) nkpt

       ! Reading first line of QE_EIGS.txt and QE_EIGS_shift.txt, and get the extra bands not needed for
       ! valence and conduction band energies respectively.
       read(valenkraw,*)dummy_int,dummy_int,vdum1,vdum2,dummy_int,dummy_int
       read(conenkraw,*)cdum1,cdum2,dummy_int,dummy_int,dummy_int,dummy_int
       write(6,*) " "
       write(6,*) "val extra bands:",vdum1,vdum2
       write(6,*) " "
       write(6,*) "con extra bands:",cdum1,cdum2
       write(6,*) " "

       ! allocote these extra arrays to read them
       allocate(venkdummy(vdum1:vdum2),cenkdummy(cdum1:cdum2))

       iband = 1 
       gamma_only = .false.
       
       do ik=1,nkpt
          write(swfc,"(i10)")ik
          write(6,*) " "
          write(6,*) "wfc index is ",swfc
          wfcname = "system.save/wfc"//trim(adjustl(swfc))//".dat"
          wfcname = trim(wfcname)
          write(6,*) " "
          write(6,*) "ik = ",ik, "wfc file is ",wfcname
          write(6,*) " "
          open( unit=qewfc, file=wfcname, form='unformatted', status='old')
          call read_a_wfc(iband, qewfc,u1dat, ik, xk, nband, isppol, nsppol, gamma_only, npw, maxnpw,&
                        b1,b2,b3 )

      ! Read the eigenvalues for each kpt and write to enkfile
          allocate(enk(nband))
          enk = 0.0
          read(valenkraw,*)enk(blftl:blfth)
          read(valenkraw,*)venkdummy
          read(conenkraw,*)cenkdummy
          read(conenkraw,*)enk(brgtl:brgth)
          write(enkdat,*)enk
          deallocate(enk)
          close(unit =qewfc )
       end do
       deallocate(venkdummy,cenkdummy)

       ! Density files
!       call grabdenqe (qeden,rhorecp,gamma_only,itemp,nsppol)	 

       
       ! Printing out the key parameters
       write(6,*) " "
       write(6,*) "*******  Key Parameters *******"
       write(6,*) " "
       write(6,*) "fermie energy = ", fermie, " Rydberg "
       write(6,*) " "
       write(6,*) "DFT gap = ", ldagap, " eV "
       write(6,*) " "
       write(6,*) "maximum # of plan waves, all kpts = ", maxnpw
       write(6,*) " "
       write(6,*) "number of bands, all kpts =", nband
       write(6,*) " "
       write(6,*) "number of spin polarizations = ",nsppol
       write(6,*) " "
!       write(6,*) "nspinor =", nspinor
!       write(6,*) " "
!       write(6,*) "nspden =", nspden
!       write(6,*) " "
       write(6,*) "number of kpts =",nkpt
       write(6,*) " "
!       write(6,*) "total number of bands, all kpts", bantot
!       write(6,*) " "
!       write(6,*) "cplex =",cplex
!       write(6,*) " "
!       write(6,*) "ngfft(1) =",ngfft(1),"ngfft(2) =",ngfft(2),"ngfft(3) =",ngfft(3)
!       write(6,*) "rprimd (:,1)=",rprimd(:,1)
!       write(6,*) " "
!       write(6,*) "rprimd(:,2)=",rprimd(:,2)
!       write(6,*) " "
!       write(6,*) "rprimd(:,3)=",rprimd(:,3)
!       write(6,*) " "
!       write(6,*) " Printing the reciprocal lattice vectors in (a.u), inv. bohr"

!       write(6,*) "b1 (a.u)              "
!       write(6,*) b1
!       write(6,*) " "
!       write(6,*) "b2 (a.u)              "
!       write(6,*) b2
!       write(6,*) " "
!       write(6,*) "b3 (a.u)              " 
!       write(6,*) b3
!       write(6,*) " "


       ! Calculating the real space primitive lattice vectors
       tempvec = 0.0
       rprimd = 0.0
       gprimd = 0.0
       gprimd(:,1) = b1
       gprimd(:,2) = b2
       gprimd(:,3) = b3
       call getrecvec ( gprimd,rprimd)
       write(6,*) " Printing the a1, a2, a3 in (a.u), bohr"
       write(6,*) " "
       write(6,*) "rprimd (:,1)=",rprimd(:,1)
       write(6,*) " "
       write(6,*) "rprimd(:,2)=", rprimd(:,2)
       write(6,*) " "
       write(6,*) "rprimd(:,3)=",rprimd(:,3)
       write(6,*) " "
!       write(6,*) " "
!       write(6,*) "gprimd (:,1)=",gprimd(:,1)
!       write(6,*) b1
!       write(6,*) " "
!       write(6,*) "gprimd(:,2)=",gprimd(:,2)
!       write(6,*) b2
!       write(6,*) " "
!       write(6,*) "gprimd(:,3)=",gprimd(:,3)
!       write(6,*) b3
!       write(6,*) " "

       ! writing the other input data to file
!       write( efermi,*) fermi_energy
       write( avecbohr,*) rprimd(:,:)
!!       write( kpoints,*) nkpt
       write( kmesh,*) nkx, nky, nkz
!       write( xmesh,*) ngfft
       write( xmesh,*) ixmesh,ixmesh,ixmesh
!       write( qunitofb,*) 0.0001,0.0001,0.0001 
       write( brange, *) blftl, blfth
       write( brange, *) brgtl, brgth 
       write( bvecinvbohr,*) gprimd(1,:)
       write( bvecinvbohr,*) gprimd(2,:)
       write( bvecinvbohr,*) gprimd(3,:)
       write( hlleps,*) decut, N_I 
       write( e0,*) diemac
       write( gwipt,*) gwgap,ldagap,0.0,0.0
       write( xsteps,*) 1
       write( spect,*) gam,elow,ehigh,ne,niter


       ! Generating the gvecs input file
       write(6,*) "----- Now checking the writing of the file u1.dat -----"
       rewind u1dat
       do ik = 1,nkpt
         read(u1dat) npw,nspinor,nband
         write(6,*) " ik = ",ik,"npw =",npw,"nspinor=",nspinor,"nband=",nband
         allocate (gvec(npw,3))
         allocate (gveclen(npw))
         allocate(indx(npw))
         allocate (zr(npw,nband),zi(npw,nband))
         allocate (zlin(2*npw*nspinor))
!         allocate (zr(brgth,npw),zi(brgth,npw))
         read(u1dat) gvec    
!         write(6,*) "size gvec array =",size(gvec),"shape=",shape(gvec)
         do ig = 1, npw 
           templen = gvec_length(gvec(ig,:),gprimd)
           gveclen(ig) = templen
!           write(6,*) " gvec coeff =",gvec(ig,:),"with length =",templen
!           write(6,*) " " 
         enddo 
!         write(6,*) " "
!         write(6,*) " The array of lengths is :"
!         write(6,*) " "
!         write(6,*) gveclen
         ! Sorting the array of lengths, and saving the index of sorted array.
         call sort_index(gveclen,indx)
!         write(6,*) " The sorted array has length :",size(gveclen)
!         write(6,*) " "
!         write(6,*) gveclen

         ! Sorting gvecs (npw,3) according to the sorting of the lengths, gveclen
         do i= 1,3
           gvec(:,i) = gvec(indx(:),i)
         enddo

         ! Generating an array of lengths that holds only unique lengths
         allocate(mask(npw))
         mask = .true.
         do ig = npw,2,-1
           mask(ig) = .not.(any(gveclen(:ig-1)==gveclen(ig)))
         enddo
         allocate(indx_vector, source=pack([(ig,ig=1,npw) ],mask))
         allocate(uniq_gveclen, source=gveclen(indx_vector))
         ! deallocate the mask so it can be used for counting the freq of elements
         deallocate(mask)
!         write(6,*) " The unique sorted array has length :",size(uniq_gveclen)
!         write(6,*) " "
!         write(6,*) uniq_gveclen
!         write(6,*) " "

         ! Counting the frequency of the unique lengths in the sorted lengths array
         allocate(cnt_uniq_gveclen(size(uniq_gveclen)))
!         allocate(mask(size(uniq_gveclen)))
         do i = 1,size(uniq_gveclen)
           mask = gveclen(:) .eq. uniq_gveclen(i)
           gvec_cnt = count(mask)
           cnt_uniq_gveclen(i) = gvec_cnt
         enddo
!         write(6,*) " The count of the unique sorted array is:"
!         write(6,*) " "
!         write(6,*) cnt_uniq_gveclen
!         write(6,*) "shape zr, zi",shape(zr),shape(zi)

         ! Now saving the sorted gvecs and their stars with their multiplicity to file gvectors
         jj = 1 ! count of the global location of gvectors
         do j=1,size(uniq_gveclen)
            write(gvecs,*)j,cnt_uniq_gveclen(j),uniq_gveclen(j)
            itemp = cnt_uniq_gveclen(j)
            do ii=1,itemp
              write(gvecs,*)jj,itemp,gvec(jj,:)
              jj = jj + 1
            end do            
         end do
         do iband =1,nband
           read(u1dat) zlin
!           read(u1dat) zi(:,iband)
         enddo
         deallocate (gvec,zr,zi,zlin,gveclen,indx)
         deallocate (uniq_gveclen,indx_vector,mask,cnt_uniq_gveclen)
!         stop
        enddo  
     

       ! close the files
!       close(unit =qewfc )
!       close(unit =qeden )
       
       close(unit =u1dat )
       close(unit =valenkraw )
       close(unit =conenkraw )
       close(unit =enkdat )
!       close(unit =occdat )
!       close(unit =rhoreal )
!       close(unit =rhorecp )
       close( unit=efermi )
       close( unit=dftgap )

!       close( unit=efermi )
       close( unit=avecbohr )
       close( unit=bvecinvbohr )
       close( unit=kpoints ) 
       close( unit=kmesh )
       close( unit=xmesh )
!       close( unit=qunitofb )
       close( unit=brange )
       close( unit=hlleps )
       close( unit=e0 )
       close( unit=gwipt )
       close( unit=gvecs )
       close( unit=xsteps )       
       close( unit=spect )       

       end program datagen

! Simple function to calculate the length of a given gvec, using its integer coefficients and the bvecs
       double precision function gvec_length( gvec, bvec) result( length )
       implicit none
       integer :: gvec( 3 )
       double precision :: bvec(3,3)
       length = ( bvec(1,1) * real(gvec(1)) + bvec(2,1) * real(gvec(2)) + bvec(3,1) * real(gvec(3)) ) ** 2.0d0 &
              + ( bvec(1,2) * real(gvec(1)) + bvec(2,2) * real(gvec(2)) + bvec(3,2) * real(gvec(3)) ) ** 2.0d0 &
              + ( bvec(1,3) * real(gvec(1)) + bvec(2,3) * real(gvec(2)) + bvec(3,3) * real(gvec(3)) ) ** 2.0d0 
       end function gvec_length

!  Subroutine to take in a list, and then generate a list that only contains unique elements in the original list.

