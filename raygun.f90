program raygun

  use xml_data_config

  implicit none

!   interface
!      subroutine init_optix()
       
!      end subroutine init_optix

!      subroutine fire_lazors()

!      end subroutine fire_lazors
!   end interface

  double precision, parameter     :: pi = 3.1415926535897932

  integer                                           :: argc, bounces
  character(len=100)                                :: file
  double precision, dimension(:, :, :), allocatable :: rays
  double precision, dimension(:, :), allocatable    :: dir, wavel


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

  print *, "a: ", hyper_a
  print *, "hyper radius: ", hyper_rad

  !CHARGIN MAH LAZOR
  call init_rays()

  !call init_optix()
  call fire_lazors(rays, dir, lobound, hibound, numrays, numoptics, par_pos, par_a, &
       hyper_pos, hyper_a, hyper_c)

  call plot_that_action(name, lobound, hibound, rays, numrays)
    
  stop

contains

  subroutine init_rays()
    logical :: check_bounds
    double precision, dimension(3,3) :: rotx, roty
    double precision :: area, gridsize, xcount = 0, ycount = 0
    logical :: done = .false.
    integer :: rayct = 1, i


    !rotation 
    rotx(1,:) = (/1.0, 0.0, 0.0/)
    rotx(2,:) = (/0.0, cos(beamrot(1)*pi/180), sin(beamrot(1)*pi/180)/)
    rotx(3,:) = (/0.0, -sin(beamrot(1)*pi/180), cos(beamrot(1)*pi/180)/)

    roty(1,:) = (/cos(beamrot(2)*pi/180), 0.0, -sin(beamrot(2)*pi/180)/)
    roty(2,:) = (/0.0, 1.0, 0.0/)
    roty(3,:) = (/sin(beamrot(2)*pi/180), 0.0, cos(beamrot(2)*pi/180)/)

    allocate(rays(100, numrays, 3))
    rays = 0.0

    allocate(dir(numrays, 3))
    !update this to calculate dir from beamcenter to primary
    !vertex. also normalize
    dir(:, 1) = 0.0
    dir(:, 2) = 0.0
    dir(:, 3) = -1.0

    if (.not. check_bounds(beamcenter, lobound, hibound)) then
       print *, "ERROR: beamcenter out of bounds"
       stop
    end if

    if (beamstyle .eq. 0) then
       area = pi * beamradius**2
       gridsize = sqrt(area/numrays) !gandalf!
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

       rays(1, rayct + 1:rayct*2, 1) = -rays(1, 1:rayct, 1)
       rays(1, rayct + 1:rayct*2, 2) = rays(1, 1:rayct, 2)

       ycount = 0
       xcount = 0
       rayct = 2*rayct + 1
       do
          if (ycount + gridsize <= sqrt(beamradius**2 - xcount**2)) then

             ycount = ycount + gridsize

             rays(1, rayct, :) = (/xcount, ycount, 0.0/)
             rayct = rayct + 1

             rays(1, rayct, :) = (/xcount, -ycount, 0.0/)
             rayct = rayct + 1
          else
             exit
          end if
       end do

       rayct = rayct - 1
       print *, "ACTUAL RAYS: ", rayct !actual rays.

    end if

    print *, "ZEROS: ", count(count(rays(1, :, :) .eq. 0, 2) .eq. 3) !lol?
    
    do i = 1, numrays
       rays(1, i, :) = matmul(rotx, rays(1, i, :))
       rays(1, i, :) = matmul(roty, rays(1, i, :))
       rays(1, i, :) = beamcenter + rays(1, i, :)
    end do

  end subroutine init_rays

end program raygun

logical function check_bounds(point, lobound, hibound) result(answer)

  implicit none

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

subroutine fire_lazors(rays, dir, lobound, hibound, numrays, numoptics, par_pos, par_a, &
     hyper_pos, hyper_a, hyper_c)

  implicit none

  integer, intent(IN)                                         :: numrays, numoptics
  integer, dimension(3), intent (IN)                          :: lobound, hibound
  double precision, dimension(3), intent(IN)                  :: par_pos, hyper_pos
  double precision, intent(IN)                                :: par_a, hyper_a, hyper_c
  double precision, dimension(100, numrays, 3), intent(INOUT) :: rays
  double precision, dimension(numrays, 3), intent(INOUT)      :: dir
  double precision, dimension(3)  :: t_arr = 0.0, prim_norm = 0.0, sec_norm = 0.0
  double precision                :: prim_a, prim_b, prim_c
  double precision                :: sec_a, sec_b, sec_c
  !double precision                :: det_a, det_b, det_c, det_norm
  double precision                :: t, a, b, c, bsquare
  integer                         :: i, bnc = 1

  
  do i = 1, numrays
     !note rederive this to include parabola location
     prim_a = (dir(i, 1)**2 + dir(i, 2)**2)/(4*par_a)
     prim_b = (2*rays(bnc, i, 1)*dir(i, 1) + 2*rays(bnc, i, 2)*dir(i, 2))/(4*par_a) - dir(i, 3)
     prim_c = (rays(bnc, i, 1)**2 + rays(bnc, i, 2)**2)/(4*par_a) - rays(bnc, i, 3)

     bsquare = prim_b**2 - 4*prim_a*prim_c

     !do a check on me.
     sec_a = ((hyper_c**2 - hyper_a**2)*(dir(i, 1)**2 + dir(i, 2)**2) - (hyper_a * dir(i, 3))**2) &
          / (hyper_a**2 * (hyper_c**2 - hyper_a**2))
     sec_b = ((hyper_c**2 - hyper_a**2)*(2*dir(i, 1)*(rays(bnc, i, 1) - hyper_pos(1)) &
          + (2 * dir(i, 2) * (rays(bnc, i, 2) - hyper_pos(2)))) + (hyper_a**2 * 2 * dir(i, 3)) &
          * (rays(bnc, i, 3) - hyper_pos(3))) / (hyper_a**2 * (hyper_c**2 - hyper_a**2))
     sec_c = ((hyper_c**2 - hyper_a**2) * (rays(bnc, i, 1)**2 - 2*rays(bnc, i, 1)*hyper_pos(1) &
          + hyper_pos(1)**2 + rays(bnc, i, 2)**2 - 2*rays(bnc, i, 2)*hyper_pos(2) &
          + hyper_pos(2)**2) + hyper_a**2 * (-rays(bnc, i, 3)**2 + 2*rays(bnc, i, 3)*hyper_pos(3) &
          - hyper_pos(3)**2 + (hyper_a**2 * (hyper_c**2 - hyper_a**2)))) / (hyper_a**2 &
          * (hyper_c**2 - hyper_a**2))

     print *, (sec_b**2 - 4*sec_a*sec_c)

     if (prim_a == 0.0) then
        if ((bsquare) >= 0.0 .and. (.true.)) then
           !go ahead
           t = -prim_c/prim_b
           t_arr = dir(i, :)*t
           rays(2, i, :) = rays(1, i, :) + t_arr
        else
           !let it go!
           print *, "ERROR, not intersect!"
        end if
     else if (bsquare == 0) then
        !right on
     else if (bsquare > 0) then
        print *, "How did we get here?"
     else if (bsquare < 0) then
        !fail
        print *, "no intersection?!"

     end if

  end do

end subroutine fire_lazors

!valgrind hates you.
subroutine plot_that_action(name, lobound, hibound, rays, numrays)
  use plplot, PI => PL_PI

  integer, intent(IN)                                  :: numrays
  integer, dimension(3), intent (IN)                   :: lobound, hibound
  character(len=40), intent(IN)                        :: name
  real(plflt), dimension(100, numrays, 3), intent(IN)  :: rays
  real(plflt), dimension(1:2)                          :: x, y
  real(plflt)                                          :: xmin, xmax, ymin, ymax, zmin, zmax
  integer                                              :: just, axis

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

  !note to self #[0x212b] is the plplot escape seq for angstrom.

  !I need black, damnit.
  call plscol0(15, 0, 0, 0)

  !name file, set plot device.  pdfcairo seems to give the best
  !bounding box. LaTeX loves it. I would take pdfcairo to prom.
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
  !call plwid(0)
  !call plline(rays(:,1,1), rays(:,1,2))
  call plline(x, y)

  call plcol0(15)
  !call plenv(zmin, zmax, ymin, ymax, just, axis)
  call plenv(-1.0, 0.4, -0.6, 0.6, just, axis)
  call pllab("Z Axis", "Y Axis", "View from X axis. #[0x212b]")
  do i = 1, numrays
     call plline(rays(1:2,i,3), rays(1:2,i,2))
  end do

  call plcol0(15)
  call plenv(zmin, zmax, xmin, xmax, just, axis)
  call pllab("Z Axis", "X Axis", "View from Y axis. #[0x212b]")
  do i = 1, numrays
     call plline(rays(1:2,i,3), rays(1:2,i,1))
  end do


  call plcol0(15)
  !call plenv(zmin, zmax, xmin, xmax, just, axis)
  call plenv(-0.6, 0.6, -0.6, 0.6, just, axis)
  call pllab("X Axis", "Y Axis", "View from detector. #[0x212b]")
  !call plssym(0.0, 1.0)
  call plpoin(rays(1,:,1), rays(1,:,2), 95)!95

  !Good news, everyone!
  call plend()

end subroutine plot_that_action
!kthxbai
