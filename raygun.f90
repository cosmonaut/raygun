program raygun

  use xml_data_config

  implicit none

  interface
!      subroutine init_optix()
       
!      end subroutine init_optix

!      subroutine fire_lazors()

!      end subroutine fire_lazors

!      subroutine plot_spot(rays, numrays, mask_ct, name, beamrot)
!        use plplot, PI => PL_PI
       
!        integer, intent(IN)                                  :: numrays
!        integer, dimension(numrays), intent(IN)              :: mask_ct
!        double precision, dimension(2), intent(IN)                    :: beamrot
!        real(plflt), dimension(100, numrays, 3), intent(IN)  :: rays
!        character(len=40), intent(IN)                        :: name  
!        double precision, dimension(:, :), allocatable       :: points
!        integer                                              :: just, axis, n = 1, good
!        double precision :: converted = 0.0
!      end subroutine plot_spot
  end interface

  double precision, parameter     :: pi = 3.1415926535897932

  integer                                           :: argc
  character(len=100)                                :: file
  double precision, dimension(:, :, :), allocatable :: rays
  integer, dimension(:), allocatable                :: mask_ct
  double precision, dimension(:), allocatable       :: wavel
  double precision, dimension(:, :), allocatable    :: dir


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

  allocate(mask_ct(numrays))
  mask_ct = 1

  allocate(wavel(numrays))
  wavel = 0.0

  !CHARGIN MAH LAZOR
  call init_rays()

  call fire_lazors(rays, dir, mask_ct, lobound, hibound, numrays, numoptics, par_pos, par_a, &
       par_rad, par_in_rad, hyper_pos, hyper_rad, hyper_a, hyper_c)

  call plot_that_action(name, lobound, hibound, rays, numrays, mask_ct, wavel, beamrot)

  call plot_spot(rays, wavel, numrays, mask_ct, name, beamrot)

  call xmlify(rays, numrays, mask_ct, name, wavel)

  deallocate(rays)
  deallocate(dir)
  deallocate(mask_ct)
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

    else if (beamstyle .eq. 1) then
       !RANDOM
       call init_random_seed()
       do i = 1, numrays
          call random_number(rays(1, i, 2))
          rays(1, i, 2) = rays(1, i, 2) - 0.5D0
          call random_number(rays(1, i, 1))
          rays(1, i, 1) = (rays(1, i, 1) - 0.5)*(abs(sqrt(beamradius**2 - rays(1, i, 2)**2))/0.5)

       end do
       rays(1, :, 3) = 0.0D0
       
    else
       print *, "ERROR, BEAMSTYLE INCORRECT"
       stop
    end if

    print *, "ZEROS: ", count(count(rays(1, :, :) .eq. 0.0, 2) .eq. 3) !lol?
    
    do i = 1, numrays
       rays(1, i, :) = matmul(rotx, rays(1, i, :))
       rays(1, i, :) = matmul(roty, rays(1, i, :))
       rays(1, i, :) = beamcenter + rays(1, i, :)
       dir(i, :) = matmul(roty, dir(i, :))
       dir(i, :) = matmul(rotx, dir(i, :))
       dir(i, :) = dir(i, :)/sqrt(dot_product(dir(i, :), dir(i, :)))

       if (rays(1, i, 1) >= 0 .and. rays(1, i, 2) >= 0) then
          wavel(i) = 7500.0
       else if (rays(1, i, 1) >= 0 .and. rays(1, i, 2) < 0) then
          wavel(i) = 5900.0
       else if (rays(1, i, 1) < 0 .and. rays(1, i, 2) < 0) then
          wavel(i) = 5700.0
       else
          wavel(i) = 4950.0
       end if
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
        print *, "POINT: ", point(:)
        print *, "ERROR, point out of bounds"
        stop
     end if
  end do

end function check_bounds

subroutine fire_lazors(rays, dir, mask_ct, lobound, hibound, numrays, numoptics, par_pos, par_a, &
     par_rad, par_in_rad, hyper_pos, hyper_rad, hyper_a, hyper_c)

  implicit none

  integer, intent(IN)                                         :: numrays, numoptics
  integer, dimension(3), intent (IN)                          :: lobound, hibound
  double precision, dimension(3), intent(IN)                  :: par_pos, hyper_pos
  double precision, intent(IN)                                :: par_a, hyper_a, hyper_c, par_rad, hyper_rad, par_in_rad
  double precision, dimension(100, numrays, 3), intent(INOUT) :: rays
  double precision, dimension(numrays, 3), intent(INOUT)      :: dir
  integer, dimension(numrays), intent(INOUT)                  :: mask_ct
  double precision, dimension(2, 20)                          :: t_pos
  logical, dimension(numrays)     :: mask
  double precision, dimension(3)  :: t_arr = 0.0, normal = 0.0
  double precision                :: prim_a, prim_b, prim_c
  double precision                :: sec_a, sec_b, sec_c
  double precision                :: det_a, det_b, det_c
  double precision                :: p_bsquare, s_bsquare, d_bsquare
  integer                         :: i, bnc = 1, t_calcd = 1
  logical                         :: switch


  t_pos = 0.0
  mask = .false.
  mask_ct = 1

  do
     do i = 1, numrays
        if (mask(i) .eqv. .false.) then
           !note rederive this to include parabola location
           prim_a = (dir(i, 1)**2 + dir(i, 2)**2)/(4.0D0*par_a)
           prim_b = (2.0D0*rays(bnc, i, 1)*dir(i, 1) + 2.0D0*rays(bnc, i, 2)*dir(i, 2))/(4.0D0*par_a) - dir(i, 3)
           prim_c = (rays(bnc, i, 1)**2 + rays(bnc, i, 2)**2)/(4.0D0*par_a) - rays(bnc, i, 3)

           p_bsquare = (prim_b**2 - 4.0D0*prim_a*prim_c)

           sec_a = ( (hyper_a**2) * (dir(i, 1)**2 + dir(i, 2)**2) - (hyper_c**2 - hyper_a**2) * (dir(i, 3))**2) &
                / ((hyper_a**2) * (hyper_c**2 - hyper_a**2))

           sec_b = ( (2.0D0*(hyper_a**2)) * (rays(bnc, i, 1)*dir(i, 1) - hyper_pos(1)*dir(i, 1) &
                + rays(bnc, i, 2)*dir(i, 2) - hyper_pos(2)*dir(i, 2)) + (2.0D0*(hyper_c**2 - hyper_a**2)) &
                * (hyper_pos(3)*dir(i, 3) - rays(bnc, i, 3)*dir(i, 3))) &
                / (hyper_a**2 * (hyper_c**2 - hyper_a**2))

           sec_c = ( (hyper_a**2) * (rays(bnc, i, 1)**2 + hyper_pos(1)**2 &
                - 2.0D0*rays(bnc, i, 1)*hyper_pos(1) + rays(bnc, i, 2)**2 + hyper_pos(2)**2 &
                - 2.0D0*rays(bnc, i, 2)*hyper_pos(2)) + (hyper_c**2 - hyper_a**2) &
                * (2.0D0*rays(bnc, i, 3)*hyper_pos(3) - rays(bnc, i, 3)**2 - hyper_pos(3)**2) &
                + (hyper_a**2)*(hyper_c**2 - hyper_a**2)) &
                / ((hyper_a**2)*(hyper_c**2 - hyper_a**2))

           s_bsquare = (sec_b**2 - 4.0D0*sec_a*sec_c)

           det_a = 0.0
           det_b = dir(i, 3)
           det_c = rays(bnc, i, 3) + 0.1

           d_bsquare = det_b**2 - 4.0D0*det_a*det_c
           
           !PRIMARY math
           if (prim_a == 0.0) then
              if (p_bsquare >= 0.0 .and. (.true.)) then
                 !go ahead
                 t_pos(:, t_calcd) = (/-prim_c/prim_b, 1.0/)
                 t_arr = dir(i, :)*t_pos(1, t_calcd)
                 if (sqrt((rays(bnc, i, 1) + t_arr(1))**2 + (rays(bnc, i, 2) + t_arr(2))**2) <= par_rad &
                      .and. sqrt((rays(bnc, i, 1) + t_arr(1))**2 + (rays(bnc, i, 2) + t_arr(2))**2) >= par_in_rad) then
                    !good
                    t_calcd = t_calcd + 1
                 else
                    !not an intersection
                    t_pos(:, t_calcd) = (/0.0, 0.0/)
                    t_arr = 0.0
                 end if

              else
                 !let it go!
                 print *, "ERROR, not intersect!"
              end if
           else if (p_bsquare == 0) then
              t_pos(:, t_calcd) = (/-prim_b/(2*prim_a), 1.0/)
              t_arr = dir(i, :)*t_pos(1, t_calcd)
              if (sqrt((rays(bnc, i, 1) + t_arr(1))**2 + (rays(bnc, i, 2) + t_arr(2))**2) <= par_rad &
                   .and. sqrt((rays(bnc, i, 1) + t_arr(1))**2 + (rays(bnc, i, 2) + t_arr(2))**2) >= par_in_rad) then
                 !good
                 t_calcd = t_calcd + 1
              else
                 t_pos(:, t_calcd) = (/0.0, 0.0/)
                 t_arr = 0.0
              end if

           else if (p_bsquare > 0) then
              !print *, "PBSQUARE > 0??????? "
              t_pos(:, t_calcd) = (/(-prim_b + sqrt(p_bsquare))/(2*prim_a), 1.0/)
              t_arr = dir(i, :)*t_pos(1, t_calcd)
              if (sqrt((rays(bnc, i, 1) + t_arr(1))**2 + (rays(bnc, i, 2) + t_arr(2))**2) <= par_rad &
                   .and. sqrt((rays(bnc, i, 1) + t_arr(1))**2 + (rays(bnc, i, 2) + t_arr(2))**2) >= par_in_rad) then    
                 t_calcd = t_calcd + 1
              else
                 t_pos(:, t_calcd) = (/0.0, 0.0/)
                 t_arr = 0.0
              end if

              t_pos(:, t_calcd) = (/(-prim_b - sqrt(p_bsquare))/(2*prim_a), 1.0/)
              t_arr = dir(i, :)*t_pos(1, t_calcd)
              if (sqrt((rays(bnc, i, 1) + t_arr(1))**2 + (rays(bnc, i, 2) + t_arr(2))**2) <= par_rad &
                   .and. sqrt((rays(bnc, i, 1) + t_arr(1))**2 + (rays(bnc, i, 2) + t_arr(2))**2) >= par_in_rad) then
                 t_calcd = t_calcd + 1
              else
                 t_pos(:, t_calcd) = (/0.0, 0.0/)
                 t_arr = 0.0
              end if
           end if

           !SECONDARY math
           if (sec_a == 0.0) then
              print *, "SECONDARY A IS 0!!!1"
              if (s_bsquare >= 0.0 .and. (.true.)) then
                 !yeah
                 t_pos(:, t_calcd) = (/ -sec_c/sec_b, 2.0/)
                 t_arr = dir(i, :)*t_pos(1, t_calcd)
                 
                 if (sqrt((rays(bnc, i, 1) + t_arr(1))**2 + (rays(bnc, i, 2) + t_arr(2))**2) <= hyper_rad &
                    .and. (rays(bnc, i, 3) + t_arr(3) >= 1.5)) then
                    !good
                    t_calcd = t_calcd + 1
                 else
                    t_pos(:, t_calcd) = (/0.0, 0.0/)
                    t_arr = 0.0
                 end if
              else
                 !how did we get here?!
                 print *, "ERROR"
                 !stop
              end if
           else if (s_bsquare == 0) then
              t_pos(:, t_calcd) = (/-sec_b/(2.0D0*sec_a), 2.0/)
              t_arr = dir(i, :)*t_pos(1, t_calcd)

              if (sqrt((rays(bnc, i, 1) + t_arr(1))**2 + (rays(bnc, i, 2) + t_arr(2))**2) <= hyper_rad &
                   .and. (rays(bnc, i, 3) + t_arr(3) >= 1.5)) then 
                 !good
                 t_calcd = t_calcd + 1
              else
                 t_pos(:, t_calcd) = (/0.0, 0.0/)
                 t_arr = 0.0
              end if
           else if (s_bsquare >= 0) then
              t_pos(:, t_calcd) = (/(-sec_b + sqrt(s_bsquare))/(2.0D0*sec_a), 2.0/)
              t_arr = dir(i, :)*t_pos(1, t_calcd)

              if (sqrt((rays(bnc, i, 1) + t_arr(1))**2 + (rays(bnc, i, 2) + t_arr(2))**2) <= hyper_rad &
                   .and. (rays(bnc, i, 3) + t_arr(3) >= 1.5)) then
                 !good
                 t_calcd = t_calcd + 1
              else
                 t_pos(:, t_calcd) = (/0.0, 0.0/)
                 t_arr = 0.0
              end if

              t_pos(:, t_calcd) = (/(-sec_b - sqrt(s_bsquare))/(2.0D0*sec_a), 2.0/)
              t_arr = dir(i, :)*t_pos(1, t_calcd)
              if (sqrt((rays(bnc, i, 1) + t_arr(1))**2 + (rays(bnc, i, 2) + t_arr(2))**2) <= hyper_rad &
                   .and. (rays(bnc, i, 3) + t_arr(3) >= 1.5)) then
                 !good
                 t_calcd = t_calcd + 1
              else
                 t_pos(:, t_calcd) = (/0.0, 0.0/)
                 t_arr = 0.0
              end if
           end if


           if (det_a == 0.0) then
              if (d_bsquare >= 0.0) then
                 t_pos(:, t_calcd) = (/ -det_c/det_b, 3.0/)
                 t_arr = dir(i, :)*t_pos(1, t_calcd)
                 !print *, "det tarr: ", t_arr
                 if (abs(rays(bnc, i, 1) + t_arr(1)) <= 0.4 .and. abs(rays(bnc, i, 2) + t_arr(2)) <= 0.4) then
                    t_calcd = t_calcd + 1
                 else
                    t_pos(:, t_calcd) = (/0.0, 0.0/)
                    t_arr = 0.0
                 end if
              else
                 print *, "ERROR"
                 !stop
              end if
           end if

           t_arr = dir(i, :) * minval(t_pos(1, :), 1, t_pos(1, :) > 0.1)
           rays(bnc + 1, i, :) = rays(bnc, i, :) + t_arr

           !Essential Logic
           if (t_pos(2, minloc(t_pos(1, :), 1, t_pos(1, :)  > 0.01)) == 1.0) then
              normal = (/(2/(4*par_a))*rays(bnc + 1, i, 1), & 
                   (2/(4*par_a))*rays(bnc + 1, i, 2), &
                   -1.0/)           
              normal = normal/(sqrt(dot_product(normal, normal)))
              
              if (dir(i, 3) > 0.0) then
                 mask(i) = .true.
                 mask_ct(i) = bnc + 1
              else
                 mask_ct = bnc + 1
              end if

              dir(i, :) = dir(i, :) - 2*normal*dot_product(dir(i, :), normal)
              dir(i, :) = dir(i, :)/(sqrt(dot_product(dir(i, :), dir(i, :))))

           else if (t_pos(2, minloc(t_pos(1, :), 1, t_pos(1, :)  > 0.01)) == 2.0) then
              !hyperboloid normal
              normal = (/(2.0D0*(rays(bnc + 1, i, 1)))/(hyper_c**2 - hyper_a**2), &
                   (2.0D0*(rays(bnc + 1, i, 2)))/(hyper_c**2 - hyper_a**2), &
                   (-2.0D0*rays(bnc + 1, i, 3) + 2.0D0*hyper_pos(3))/(hyper_a**2)/)
              normal = normal/(sqrt(dot_product(normal, normal)))
              normal(3) = -abs(normal(3))
              if (dir(i, 3) < 0.0) then
                 mask(i) = .true.
                 mask_ct(i) = bnc + 1
              else
                 mask_ct(i) = bnc + 1
              end if

              !hax(:) = dir(i, :)
              dir(i, :) = (dir(i, :) - 2*normal*dot_product(dir(i, :), normal))
              dir(i, :) = dir(i, :)/(sqrt(dot_product(dir(i, :), dir(i, :))))

           else if (t_pos(2, minloc(t_pos(1, :), 1, t_pos(1, :)  > 0.01)) == 3.0) then

              mask(i) = .true.
              mask_ct(i) = bnc + 1
              
           else
              !print *, "EPIC ERROR, NO T"
              mask(i) = .true.
              mask_ct(i) = bnc + 1
           end if
           
           t_calcd = 1
           t_pos = 0.0

           switch = .false.
        end if
     end do
     print *, "bounce", bnc
     bnc = bnc + 1
     if (switch .eqv. .true.) then
       exit
     end if
     switch = .true.
     !exit
  end do
  !print *, mask

end subroutine fire_lazors


!valgrind hates you.
subroutine plot_that_action(name, lobound, hibound, rays, numrays, mask_ct, wavel, beamrot)
  use plplot, PI => PL_PI

  integer, intent(IN)                                  :: numrays
  integer, dimension(3), intent (IN)                   :: lobound, hibound
  integer, dimension(numrays), intent(IN)              :: mask_ct
  double precision, dimension(numrays), intent(IN)     :: wavel
  character(len=40), intent(IN)                        :: name
  real(plflt), dimension(100, numrays, 3), intent(IN)  :: rays
  double precision, dimension(2), intent(IN)           :: beamrot
  real(plflt), dimension(40)                           :: x, y, xx, yy
  real(plflt)                                          :: xmin, xmax, ymin, ymax, zmin, zmax
  integer                                              :: just, axis, good = 0

  double precision, dimension(:, :), allocatable       :: points
  integer, dimension(:), allocatable                   :: color
  double precision                                     :: converted = 0.0, buffer_x = 0.0, buffer_y = 0.0
  character(len=10)                                    :: display

  xmin = lobound(1)
  ymin = lobound(2)
  xmax = hibound(1)
  ymax = hibound(2)
  zmin = lobound(3)
  zmax = hibound(3)

  just = 1
  axis = 0

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
  call plenv(-0.6, 0.6, -0.6, 0.6, just, axis)
  call pllab("X Axis (Meters)", "Y Axis (Meters)", "View of Initial Beam from Z axis")

  !call plcol0(1)
  do i = 1, numrays
     if (wavel(i) == 7500.0) then
        call plcol0(1)
     else if (wavel(i) == 5900.0) then
        call plcol0(2)
     else if (wavel(i) == 5700.0) then
        call plcol0(3)
     else if (wavel(i) == 4950.0) then
        call plcol0(9)
     else
        !wat?
     end if

     call plpoin((/rays(1, i, 1)/), (/rays(1, i, 2)/), 95)
     !do j = 1, mask_ct(i) - 1
        !call plline(rays(j:j+1,i,1), rays(j:j+1,i,2))
     !end do
  end do


  call plcol0(15)
  !call plenv(zmin, zmax, ymin, ymax, just, axis)
  call plenv(-0.3, 3.2, -0.6, 0.6, just, axis)
  call pllab("Z Axis (Meters)", "Y Axis (Meters)", "View from X Axis")
  do i = 1, numrays
     if (abs(rays(1, i, 1)) <= 0.01) then
        do j = 1, mask_ct(i) - 1
           call plline(rays(j:j+1,i,3), rays(j:j+1,i,2))
        end do
     else
        !DON'T!!!!11111oneoneoneone
     end if
  end do

  step = 1.0/16000.0

  x = 0.0
  x(1) = 2.59570709643376
  do i = 1, 40
     if (i == 1) then
        !nono
     else
        x(i) = x(i - 1) + step
     end if
     y(i) = sqrt((((x(i)-1.45)**2)/(1.1457070964337**2) - 1)*(1.55**2 - 1.1457070964337**2))
  end do

  step_par = 1.0/1900.0
  xx = (0.045**2)/(4*1.14570709)
  do i = 1, 40
     if (i == 1) then
        !non
     else
        xx(i) = xx(i - 1) + step_par
     end if
     yy(i) = sqrt(12.0*xx(i))
  end do

  call plcol0(1)
  call plline(x, y)

  y = -y
  call plline(x, y)

  call plline(xx, yy)
  yy = -yy
  call plline(xx, yy)


  call plcol0(15)
  !call plenv(zmin, zmax, xmin, xmax, just, axis)
  call plenv(-0.3, 3.2, -0.6, 0.6, just, axis)
  call pllab("Z Axis (Meters)", "X Axis (Meters)", "View from Y Axis")

  do i = 1, numrays
     if (abs(rays(1, i, 2)) <= 0.01) then
        do j = 1, mask_ct(i) - 1
           call plline(rays(j:j+1,i,3), rays(j:j+1,i,1))
        end do
     else
        !ZIP
     end if
  end do

  !#[0x212b]
  call plcol0(1)
  call plline(x, y)
  y = -y
  call plline(x, y)

  call plline(xx, yy)
  yy = -yy
  call plline(xx, yy)


  good = count(mask_ct(:) .eq. maxval(mask_ct))

  allocate(points(good, 2))
  allocate(color(good))
  !print *, size(points)
  do i = 1, numrays
     if (mask_ct(i) == maxval(mask_ct(:))) then
        points(n, :) = (/rays(maxval(mask_ct), i, 1), rays(maxval(mask_ct), i, 2)/)

        if (wavel(i) == 7500.0) then
           color(n) = 1
        else if (wavel(i) == 5900.0) then
           color(n) = 2
        else if (wavel(i) == 5700.0) then
           color(n) = 3
        else if (wavel(i) == 4950.0) then
           color(n) = 9
        else
           !wat?
        end if
        n = n + 1
     else
        !wat
     end if
  end do

  call plcol0(15)
  xmin = minval(points(:, 1))*1.0e6
  xmax = maxval(points(:, 1))*1.0e6
  ymin = minval(points(:, 2))*1.0e6
  ymax = maxval(points(:, 2))*1.0e6
  buffer_x = (maxval(points(:, 1)) - minval(points(:, 1)))*1.0e6/50
  buffer_y = (maxval(points(:, 2)) - minval(points(:, 2)))*1.0e6/50
  call plsxax(5, 5)
  call plenv(xmin - buffer_x, xmax + buffer_x, &
      ymin - buffer_y, ymax + buffer_y, 0, axis)

  converted = beamrot(2)*3600.0
  write (display, '(I4)') int(converted)

  call pllab("X Axis (#gmm)", "Y Axis (#gmm)", "Detector View:" //trim(display) &
       // " Arcseconds Off Axis")

  do i = 1, good
     call plcol0(color(i))
     call plpoin((/points(i, 1)*1.0e6/), (/points(i, 2)*1.0e6/), 95)
  end do

  print *, "Plotted"

  call plend()

end subroutine plot_that_action


!valgrind also hates you
subroutine plot_spot(rays, wavel, numrays, mask_ct, name, beamrot)
  use plplot, PI => PL_PI

  integer, intent(IN)                                  :: numrays
  integer, dimension(numrays), intent(IN)              :: mask_ct
  double precision, dimension(numrays), intent(IN)     :: wavel
  double precision, dimension(2), intent(IN)           :: beamrot
  real(plflt), dimension(100, numrays, 3), intent(IN)  :: rays
  character(len=40), intent(IN)                        :: name  
  double precision, dimension(:, :), allocatable       :: points
  integer, dimension(:), allocatable                   :: color
  integer                                              :: just, axis, n = 1, good
  double precision                                     :: converted = 0.0, buffer_x = 0.0, buffer_y = 0.0
  double precision                                     :: xmin, xmax, ymin, ymax
  character(len=10)                                    :: display


  print *, "Plotting Spot..."
  
  just = 1
  axis = 0

  !I need black, damnit.
  call plscol0(15, 0, 0, 0)

  !name file, set plot device.  pdfcairo seems to give the best
  !bounding box. LaTeX loves it. I would take pdfcairo to prom.
  call plsfnam(trim(name) // "spot" // ".pdf")
  call plsdev("pdfcairo")

  call plscolbg(255, 255, 255)
  call plinit()
  call plfont(2)
  
  print *, count(mask_ct(:) .eq. maxval(mask_ct))
  good = count(mask_ct(:) .eq. maxval(mask_ct))

  allocate(points(good, 2))
  allocate(color(good))

  do i = 1, numrays
     if (mask_ct(i) == maxval(mask_ct(:))) then
        points(n, :) = (/rays(maxval(mask_ct), i, 1), rays(maxval(mask_ct), i, 2)/)

        if (wavel(i) == 7500.0) then
           color(n) = 1
        else if (wavel(i) == 5900.0) then
           color(n) = 2
        else if (wavel(i) == 5700.0) then
           color(n) = 3
        else if (wavel(i) == 4950.0) then
           color(n) = 9
        else
           !what?
        end if
        n = n + 1
     else
        !do not want
     end if
  end do

  call plcol0(15)
  xmin = minval(points(:, 1))*1.0e6
  xmax = maxval(points(:, 1))*1.0e6
  ymin = minval(points(:, 2))*1.0e6
  ymax = maxval(points(:, 2))*1.0e6
  buffer_x = (maxval(points(:, 1)) - minval(points(:, 1)))*1.0e6/50
  buffer_y = (maxval(points(:, 2)) - minval(points(:, 2)))*1.0e6/50
  call plsxax(5, 5)
  call plenv(xmin - buffer_x, xmax + buffer_x, &
      ymin - buffer_y, ymax + buffer_y, just, axis)

  converted = beamrot(2)*3600.0
  write (display, '(I4)') int(converted)

  call pllab("X Axis (#gmm)", "Y Axis (#gmm)", "View from Detector " &
       // trim(adjustl(display)) &
       !// " 3000" &
       // " Arcseconds Off Axis")


  do i = 1, good
     call plcol0(color(i))
     call plpoin((/points(i, 1)*1.0e6/), (/points(i, 2)*1.0e6/), 95)
  end do

  print *, "Plotted Spot."
  call plend()
  
  deallocate(points)

end subroutine plot_spot

subroutine xmlify(rays, numrays, mask_ct, name, wavel)
  integer, intent(IN)                                      :: numrays
  integer, dimension(numrays), intent(IN)                  :: mask_ct
  double precision, dimension(numrays), intent(IN)         :: wavel
  double precision, dimension(100, numrays, 3), intent(IN) :: rays
  character(len=40), intent(IN)                            :: name
  double precision, dimension(:, :, :), allocatable        :: good_rays, r, o, y, g, b, v
  double precision, dimension(:), allocatable              :: good_wave
  integer                                                  :: n, good, n_good_rays, j, k, h
  character(len=20) :: num, pt, xpt, ypt, zpt

  n = maxval(mask_ct)
  good = count(mask_ct .eq. n)

  print *, "N: ", n
  print *, "good: ", good

  allocate(good_rays(n, good, 3))
  allocate(good_wave(good))

  j = 1
  do i = 1, numrays
     if (mask_ct(i) == n) then
        do k = 1, n
           good_rays(k, j, :) = rays(k, i, :)
           good_wave(j) = wavel(i)
        end do
        j = j + 1 
     end if
  end do

  print *, "done"

  print *, count((good_wave <= 7500.0 .and. good_wave > 6200.0))

  open(unit = 1, file = trim(name) // "rt3d.xml")

  print *, "file open"

  allocate(r(n, count((good_wave <= 7500.0 .and. good_wave > 6200.0)), 3)) !r
  allocate(o(n, count((good_wave <= 6200.0 .and. good_wave > 5900.0)), 3)) !o
  allocate(y(n, count((good_wave <= 5900.0 .and. good_wave > 5700.0)), 3)) !y
  allocate(g(n, count((good_wave <= 5700.0 .and. good_wave > 4950.0)), 3)) !g
  allocate(b(n, count((good_wave <= 4950.0 .and. good_wave > 4500.0)), 3)) !b
  allocate(v(n, count((good_wave <= 4500.0 .and. good_wave > 3800.0)), 3)) !v
  !print *, "hello?"

  j = 1
  do i = 1, good
     if (good_wave(i) <= 7500.0 .and. good_wave(i) > 6200.0) then
        do k = 1, n
           r(k, j, :) = good_rays(k, i, :)
        end do
        j = j + 1
     end if
  end do


  j = 1
  do i = 1, good
     if (good_wave(i) <= 6200.0 .and. good_wave(i) > 5900.0) then
        do k = 1, n
           o(k, j, :) = good_rays(k, i, :)
        end do
        j = j + 1
     end if
  end do

  j = 1
  do i = 1, good
     if (good_wave(i) <= 5900.0 .and. good_wave(i) > 5700.0) then
        do k = 1, n
           y(k, j, :) = good_rays(k, i, :)
        end do
        j = j + 1
     end if
  end do

  j = 1
  do i = 1, good
     if (good_wave(i) <= 5700.0 .and. good_wave(i) > 4950.0) then
        do k = 1, n
           g(k, j, :) = good_rays(k, i, :)
        end do
        j = j + 1
     end if
  end do


  j = 1
  do i = 1, good
     if (good_wave(i) <= 4950.0 .and. good_wave(i) > 4500.0) then
        do k = 1, n
           b(k, j, :) = good_rays(k, i, :)
        end do
        j = j + 1
     end if
  end do

  j = 1
  do i = 1, good
     if (good_wave(i) <= 4500.0 .and. good_wave(i) > 3800.0) then
        do k = 1, n
           v(k, j, :) = good_rays(k, i, :)
        end do
        j = j + 1
     end if
  end do


  write(1, *) "<raytrace><lines>"

  if (size(r, 2) /= 0) then  
     write(1, *) "<rlines>"
     do i = 1, size(r, 2)
        write(num, *) i
        !dummy = trim(trim(dummy))
        write (1, *) "<line" // trim(adjustl(num)) // ">"
        do k = 1, n
           !print *, "POINTS!"
           write(pt, *) k - 1
           write (1, *) "<point" // trim(adjustl(pt)) // ">"
           !do h = 1, 3
           
           write (xpt, '(F20.16)') r(k, i, 1)
           !print *, good_rays(i, 1)
           write (ypt, '(F20.15)') r(k, i, 2)
           write (zpt, '(F20.15)') r(k, i, 3)
           write (1, *) "<x>" // trim(adjustl(xpt)) // "</x>"
           write (1, *) "<y>" // trim(adjustl(ypt)) // "</y>"
           write (1, *) "<z>" // trim(adjustl(zpt)) // "</z>"
           
           !end do
           write (1, *) "</point" // trim(adjustl(pt)) // ">"
        end do
        
        write (1, *) "</line" // trim(adjustl(num)) // ">"
     end do
     write(1, *) "</rlines>"
  end if
  
  if (size(o, 2) /= 0) then
     write(1, *) "<olines>"
     do i = 1, size(o, 2)
        write(num, *) i
        !dummy = trim(trim(dummy))
        write (1, *) "<line" // trim(adjustl(num)) // ">"
        do k = 1, n
           !print *, "POINTS!"
           write(pt, *) k - 1
           write (1, *) "<point" // trim(adjustl(pt)) // ">"
           !do h = 1, 3
           
           write (xpt, '(F20.16)') o(k, i, 1)
           !print *, good_rays(i, 1)
           write (ypt, '(F20.15)') o(k, i, 2)
           write (zpt, '(F20.15)') o(k, i, 3)
           write (1, *) "<x>" // trim(adjustl(xpt)) // "</x>"
           write (1, *) "<y>" // trim(adjustl(ypt)) // "</y>"
           write (1, *) "<z>" // trim(adjustl(zpt)) // "</z>"
           
           !end do
           write (1, *) "</point" // trim(adjustl(pt)) // ">"
        end do
        
        write (1, *) "</line" // trim(adjustl(num)) // ">"
     end do
     write(1, *) "</olines>"
  end if

  if (size(y, 2) /= 0) then
     write(1, *) "<ylines>"
     do i = 1, size(y, 2)
        write(num, *) i
        !dummy = trim(trim(dummy))
        write (1, *) "<line" // trim(adjustl(num)) // ">"
        do k = 1, n
           !print *, "POINTS!"
           write(pt, *) k - 1
           write (1, *) "<point" // trim(adjustl(pt)) // ">"
           !do h = 1, 3
           
           write (xpt, '(F20.16)') y(k, i, 1)
           !print *, good_rays(i, 1)
           write (ypt, '(F20.15)') y(k, i, 2)
           write (zpt, '(F20.15)') y(k, i, 3)
           write (1, *) "<x>" // trim(adjustl(xpt)) // "</x>"
           write (1, *) "<y>" // trim(adjustl(ypt)) // "</y>"
           write (1, *) "<z>" // trim(adjustl(zpt)) // "</z>"
           
           !end do
           write (1, *) "</point" // trim(adjustl(pt)) // ">"
        end do
        
        write (1, *) "</line" // trim(adjustl(num)) // ">"
     end do
     write(1, *) "</ylines>"
  end if

  if (size(g, 2) /= 0) then
     write(1, *) "<glines>"
     do i = 1, size(g, 2)
        write(num, *) i
        !dummy = trim(trim(dummy))
        write (1, *) "<line" // trim(adjustl(num)) // ">"
        do k = 1, n
           !print *, "POINTS!"
           write(pt, *) k - 1
           write (1, *) "<point" // trim(adjustl(pt)) // ">"
           !do h = 1, 3
           
           write (xpt, '(F20.16)') g(k, i, 1)
           write (ypt, '(F20.15)') g(k, i, 2)
           write (zpt, '(F20.15)') g(k, i, 3)
           write (1, *) "<x>" // trim(adjustl(xpt)) // "</x>"
           write (1, *) "<y>" // trim(adjustl(ypt)) // "</y>"
           write (1, *) "<z>" // trim(adjustl(zpt)) // "</z>"
           
           !end do
           write (1, *) "</point" // trim(adjustl(pt)) // ">"
        end do
        
        write (1, *) "</line" // trim(adjustl(num)) // ">"
     end do
     write(1, *) "</glines>"
  end if

  if (size(b, 2) /= 0) then
     write(1, *) "<blines>"
     do i = 1, size(b, 2)
        write(num, *) i
        !dummy = trim(trim(dummy))
        write (1, *) "<line" // trim(adjustl(num)) // ">"
        do k = 1, n
           !print *, "POINTS!"
           write(pt, *) k - 1
           write (1, *) "<point" // trim(adjustl(pt)) // ">"
           !do h = 1, 3
           
           write (xpt, '(F20.16)') b(k, i, 1)
           write (ypt, '(F20.15)') b(k, i, 2)
           write (zpt, '(F20.15)') b(k, i, 3)
           write (1, *) "<x>" // trim(adjustl(xpt)) // "</x>"
           write (1, *) "<y>" // trim(adjustl(ypt)) // "</y>"
           write (1, *) "<z>" // trim(adjustl(zpt)) // "</z>"
           
           !end do
           write (1, *) "</point" // trim(adjustl(pt)) // ">"
        end do
        
        write (1, *) "</line" // trim(adjustl(num)) // ">"
     end do
     write(1, *) "</blines>"
  end if

  if (size(v, 2) /= 0) then
     write(1, *) "<vlines>"
     do i = 1, size(v, 2)
        write(num, *) i
        !dummy = trim(trim(dummy))
        write (1, *) "<line" // trim(adjustl(num)) // ">"
        do k = 1, n
           !print *, "POINTS!"
           write(pt, *) k - 1
           write (1, *) "<point" // trim(adjustl(pt)) // ">"
           !do h = 1, 3
           
           write (xpt, '(F20.16)') v(k, i, 1)
           write (ypt, '(F20.15)') v(k, i, 2)
           write (zpt, '(F20.15)') v(k, i, 3)
           write (1, *) "<x>" // trim(adjustl(xpt)) // "</x>"
           write (1, *) "<y>" // trim(adjustl(ypt)) // "</y>"
           write (1, *) "<z>" // trim(adjustl(zpt)) // "</z>"
           
           !end do
           write (1, *) "</point" // trim(adjustl(pt)) // ">"
        end do
        
        write (1, *) "</line" // trim(adjustl(num)) // ">"
     end do
     write(1, *) "</vlines>"
  end if

  ! print *, count((wavel <= 7500.0 .and. wavel > 6200.0)) !r
  ! print *, count((wavel <= 6200.0 .and. wavel > 5900.0)) !o
  ! print *, count((wavel <= 5900.0 .and. wavel > 5700.0)) !y
  ! print *, count((wavel <= 5700.0 .and. wavel > 4950.0)) !g
  ! print *, count((wavel <= 4950.0 .and. wavel > 4500.0)) !b
  ! print *, count((wavel <= 4500.0 .and. wavel > 3800.0)) !v




  write(1, *) "</lines></raytrace>"
  close(1)

end subroutine xmlify

SUBROUTINE init_random_seed()

  implicit none

  INTEGER                            :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  
  CALL SYSTEM_CLOCK(COUNT=clock)
  
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
  
  DEALLOCATE(seed)
END SUBROUTINE init_random_seed
!kthxbai
