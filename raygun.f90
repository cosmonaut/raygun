program raygun

  use xml_data_config

  implicit none

  double precision, parameter               :: pi = 3.1415926
  integer                                   :: argc
  character(len=100)                        :: file
  double precision                          :: test
  !double precision, dimension(1:30, 1:3000, 1:3)  :: rays = 0
  double precision, dimension(:, :, :), allocatable :: rays
    
  print *,"O HAI WORLDS!"

  test = 1.34567834563987498374
  print *, test

  argc = IARGC()

  if(argc > 1) then
     print *, "I CAN HAZ LESS ARGUMENTS PLOX?"
     stop
  endif

  call getarg(1, file)

  if(file .ne. "") then
     print *, "USING FILE: ", file
     call read_xml_file_config(file)
  else
     print *, "USING FILE: ", "test.xml"
     call read_xml_file_config("test.xml")
  endif

  !CHARGIN MAH LAZOR
  call fill_ray_positions()

  call plot_that_action(name, lobound, hibound)

  stop

contains

  subroutine fill_ray_positions()
    integer :: i, direction

    allocate(rays(100, numrays, 3))
    
    rays(1, 1, :) = (/beamcenter(:)/)
    print *, rays(1, 1, :)

    !print *, check_bounds(rays(1,1,:), lobound, hibound)

    !print *, rayplane(1)*beamcenter(1) + rayplane(2)*beamcenter(2) &
    !     + rayplane(3)*beamcenter(3) + rayplane(4)
    !Check mah planzez
    if (rayplane(1)*beamcenter(1) + rayplane(2)*beamcenter(2) &
         + rayplane(3)*beamcenter(3) + rayplane(4) == 0) then
       print *, "HELL, YES. WE GOTS A VALID PLANE"
    else
       print *, "ERROR: CENTOR OF BEEMZ NOTS IN RAYPLANE!!!!"
       stop
    end if
    
!     if (beamcenter(1) /= 0) then
!        direction = 1
!     else if (beamcenter(2) /= 0) then
!        direction = 2
!     else
!        direction = 3
!     end if

    do i=1, numrays
       
    end do


  end subroutine fill_ray_positions

end program raygun

logical function check_bounds(point, lobound, hibound) result(answer)
  integer, dimension(3), intent(IN)          :: lobound, hibound
  double precision, dimension(3), intent(IN) :: point
  integer                                    :: i

  do i=1, 3
     if (point(i) <= hibound(i) .and. point(i) >= lobound(i)) then
        answer = .true.
     else
        answer = .false.
        print *, "ERROR, point out of bounds"
        stop
     end if
  end do
end function check_bounds

!valgrind hates you.
subroutine plot_that_action(name, lobound, hibound)
  use plplot, PI => PL_PI

  integer, dimension(3), intent (IN)       :: lobound, hibound
  character(len=40), intent(IN)            :: name
  real(plflt), dimension(1:2)              :: x, y
  real(plflt)                              :: xmin, xmax, ymin, ymax, zmin, zmax, rot
  integer                                  :: just, axis

  xmin = lobound(1)
  ymin = lobound(2)
  xmax = hibound(1)
  ymax = hibound(2)
  zmin = lobound(3)
  zmax = hibound(3)

  just = 1
  axis = 0

  x(1) = 0
  x(2) = 5
  y(1) = 0
  y(2) = 5

  rot = 1

  !We will always use Z as the optical axis.

  !I need black damnit.
  call plscol0(15, 0, 0, 0)

  !name file, set plot device, and rotate the plot.  pdfcairo seems to
  !give the best bounding box. LaTeX loves it. I would take pdfcairo
  !to prom.
  call plsfnam(trim(name) // ".pdf")
  call plsdev("pdfcairo")

  call plssub(2, 2)

  !set background, initialize, font crap.
  call plscolbg(255, 255, 255)
  call plinit()
  call plfont(2)
  !call plfontld(1)

  call plcol0(15)
  call plenv(xmin, xmax, ymin, ymax, just, axis)
  call pllab("X Axis", "Y Axis", "View from Z axis. #[0x212b]")

  call plcol0(1)
  call plline(x, y)

  call plcol0(15)
  call plenv(zmin, zmax, ymin, ymax, just, axis)
  call pllab("Z Axis", "Y Axis", "View from X axis. #[0x212b]")

  call plcol0(15)
  call plenv(zmin, zmax, xmin, xmax, just, axis)
  call pllab("Z Axis", "X Axis", "View from Y axis. #[0x212b]")

  call plcol0(15)
  call plenv(zmin, zmax, xmin, xmax, just, axis)
  call pllab("X Axis", "Y Axis", "View from detector. #[0x212b]")

  call plend()

end subroutine plot_that_action
!kthxbai

