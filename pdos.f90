program VASP_DOSCAR
implicit none

!=======================================================================
!
!      makes individual DOSCAR files for total system and each atom
!              From VASP-generated DOSCAR with LROBIT=11
!                  Soham Subhra Ghosh
!                     02/12/14
!        Modified on 11/13/2014
!        Modified on 02/24/2016
!   
!        type ./pdos --help to get the execution syntax
!
!=======================================================================

    
    integer, parameter :: outunit=44  !holds the unit number of pdos output files created 

    character(len=70) :: fn           ! variable used to build filenames 
    character(LEN=12) :: string       ! used to read variables not needed
    character(LEN=1) :: max_l         ! largest ang. mom in the system (s/p/d/f)

    integer :: NEDOS                  ! number of DOS data points
    integer :: numfiles               ! number of pdos output files == number of ions 
    integer :: filenum,i,j,ISPIN,flag1,flag2
    
    real :: energy=0.0,Efermi
    
    real,allocatable,dimension(:) :: tot_density,partial_density

    character(len=32),dimension(3) :: arg

    ! read command line argument and give out useful help if asked
    do i = 1,2
      CALL get_command_argument(i,arg(i))
      if ((len_trim(arg(i)) .eq.0) .or. (arg(i) .eq. "--help")) then
        print*," "
        print*," usage: ./pdos  ISPIN max_l < filename(=generally DOSCAR)"
        print*," "
        print*,"ISPIN is the VASP variable. max_l is the maximum angular momentum"
        print*,"max_l = 's' or 'p' or 'd' or 'f'"
        print*,"Example: ./pdos 2 d < DOSCAR"
        CALL EXIT(0)
      end if
    end do


    ! array sizes for the total_density array, and partial density array resp.  
    flag1=0;flag2=0 
    read(arg(1),*)ISPIN
    open(unit=2,file='totaldos.dat')
    ! ISPIN =1 for no spinor, ISPIN = 2 for spinors
    ! the following syntant reads string type arg(1) into integer type ISPIN
    ! max_l is the largest angular momentum in the system
    ! max_l = s,p,d or f
    max_l = arg(2)
    ! if ISPIN = 1, flag1 will hold 2 values (dos, summed_dos)
    ! if ISPIN = 2, flag1 holds 4 values (up_dos,dn_dos,smmed_up_dos,
    ! smmed_dn_dos
    flag1 = ISPIN*2 
    ! there are ISPIN*n^2 variables in the pdos array for one atom and one energy
    if (max_l.eq.'s') then
      flag2 = ISPIN*1 
    else if (max_l.eq.'p') then
      flag2 = ISPIN*4
    else if (max_l.eq.'d') then
      flag2 = ISPIN*9
    else if (max_l.eq.'f') then
      flag2 = ISPIN*16 
    else
      print*,"proper s/p/d/f value not found"
      call EXIT(0)
    end if
    ! allocate the array with known sizes
    allocate(tot_density(flag1))
    allocate(partial_density(flag2))
    ! numfiles == number of atoms == number of pdos files to be created
    read(*,*)numfiles
    ! disregard the next 4 lines
    do i = 1,4
        read(*,*) 
    end do
    ! disregard the first two variables
    read(*,*)string,string,NEDOS,Efermi
!-------------------------------------------------
! This block reads the total dos
!-------------------------------------------------
    do i = 1,NEDOS
       read(*,*)energy,tot_density
       energy = energy - Efermi
       ! make the downspin dos negative so it goes below the y axis
       if (ISPIN.eq.2) then 
         tot_density(2) = - tot_density(2)
       end if
       write(2,*)energy,tot_density
    end do
    close(2)
!-------------------------------------------------
! this block reads the site and Y_lm projected dos
!-------------------------------------------------
    do filenum=1,numfiles
        ! build filename -- i.dat
        ! this writes a file named filenum,'dat' to the string fn
        ! Where filenum is a variable
        write(fn,fmt='(i0,a)')filenum,'.dat'
    
        ! now we can open that file by calling the string 'fn'  
        ! open it with a fixed unit number
        open(unit=outunit,file=fn, form='formatted')
        read(*,*)
        do i = 1,NEDOS
          read(*,*)energy,partial_density
          energy = energy - Efermi
          ! make all the downspin dos negative so they go below the y axis
          if (ISPIN.eq.2) then
            do j = 2,flag2,2
                 partial_density(j) = -partial_density(j)
            end do
          end if
            write(outunit,*)energy,partial_density
 
        end do
        ! close it 
        close(outunit)
    end do

    deallocate(tot_density)
    deallocate(partial_density)
end program VASP_DOSCAR
