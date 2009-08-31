program raygun

  use xml_data_config

  implicit none

  double precision, parameter               :: pi = 3.1415926535897932

  integer                                           :: argc
  character(len=100)                                :: file
  double precision, dimension(:, :, :), allocatable :: rays
    
  print *,"O HAI WORLDS!"

  print *, "PI: ", pi

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

  call plot_that_action(name, lobound, hibound, rays, numrays)

  stop

contains

  subroutine fill_ray_positions()
    logical :: check_bounds
    double precision, dimension(3,3) :: rotx, roty
    double precision :: area, gridsize, xcount = 0, ycount = 0
    logical :: done = .false.
    integer :: rayct = 1

    !rotation
    rotx(1,:) = (/1.0, 0.0, 0.0/)
    rotx(2,:) = (/0.0, cos(beamrot(1)), sin(beamrot(1))/)
    rotx(3,:) = (/0.0, -sin(beamrot(1)), cos(beamrot(1))/)

    roty(1,:) = (/cos(beamrot(2)), 0.0, -sin(beamrot(2))/)
    roty(2,:) = (/0.0, 1.0, 0.0/)
    roty(3,:) = (/sin(beamrot(2)), 0.0, cos(beamrot(2))/)

    allocate(rays(100, numrays, 3))
    rays = 0

    if (.not. check_bounds(beamcenter, lobound, hibound)) then
       print *, "ERROR: beamcenter out of bounds"
       stop
    end if

    if (beamstyle .eq. 0) then
       area = pi * beamradius**2
       gridsize = sqrt(area/numrays)
       print *, "GRID SIZE" , gridsize

       do
          if (xcount + gridsize <= sqrt(beamradius**2 - ycount**2)) then
             done = .false.
             xcount = xcount + gridsize
             rays(1, rayct, :) = (/xcount, ycount, 0.0/) 
             rayct = rayct + 1
             if (ycount /= 0) then
                rays(1, rayct, :) = (/xcount, -ycount, 0.0/) 
                rayct = rayct + 1
             end if
          else
             if (done .eqv. .true.) then
                exit
             end if
             ycount = ycount + gridsize
             xcount = 0
             done = .true.
          end if
       end do

       rayct = rayct - 1 !we added one more after final assignment
       !print *, count
       rays(1, rayct + 1:rayct*2, 1) = -rays(1, 1:rayct, 1)
       rays(1, rayct + 1:rayct*2, 2) = rays(1, 1:rayct, 2)
       
       ycount = 0
       xcount = 0
       rayct = 2*rayct + 1
       do
          if (ycount + gridsize <= sqrt(beamradius**2 - xcount**2)) then
             !print *,count
             ycount = ycount + gridsize

             rays(1, rayct, :) = (/xcount, ycount, 0.0/)
             rayct = rayct + 1

             rays(1, rayct, :) = (/xcount, -ycount, 0.0/)
             rayct = rayct + 1
          else
             exit
          end if
       end do
    end if

    rayct = rayct - 1
    print *, "ACTUAL RAYS: ", rayct !actual rays.

    print *, count(count(rays(1, :, :) .eq. 0, 2) .eq. 3)
    
!     do i = 1, 3000
!        print *, rays(1,i,2), i
!        if (rays(1,i,2) .eq. 0.0) then
!           fake = fake + 1
!        end if
!     end do
!     print *, "FAKE ", fake

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
subroutine plot_that_action(name, lobound, hibound, rays, numrays)
  use plplot, PI => PL_PI

  integer                       :: numrays
  integer, dimension(3), intent (IN)       :: lobound, hibound
  character(len=40), intent(IN)            :: name
  real(plflt), dimension(100, numrays, 3)     :: rays
  real(plflt), dimension(1:2)              :: x, y
  real(plflt)                              :: xmin, xmax, ymin, ymax, zmin, zmax
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

  print *, "PLOTTING..."
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
  !call plenv(zmin, zmax, xmin, xmax, just, axis)
  call plenv(-5.0, 5.0 , -5.0, 5.0, just, axis)
  call pllab("X Axis", "Y Axis", "View from detector. #[0x212b]")
  call plpoin(rays(1,:,1), rays(1,:,2), 95)

  call plend()

end subroutine plot_that_action
!kthxbai
