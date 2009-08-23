program raygun

  use plplot, PI => PL_PI
  use xml_data_config

  implicit none

  integer                            :: argc
  character(len=100)                 :: file
  double precision                   :: test
  double precision, dimension(3000)  :: rays = 0

  test = 1.34567834563987498374
  print *, test

  argc = IARGC()

  if(argc > 1) then
     print *, "ERROR, too many arguments"
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

  

  print *, rays(1)
  call plot_that_shit()
  
  !call xml_process( 'test.xml', attribs, data, startfunc, datafunc, endfunc, 0, error )
  
  !print *,"O HAI WORLDS!"
  
contains

  subroutine plot_that_shit()
    real(plflt), dimension(1:2) :: x, y
    real(plflt)                 :: xmin, xmax, ymin, ymax, zmin, zmax, rot
    integer                     :: just, axis

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
    
    !change color index 15 to be black
    call plscol0(15, 0, 0, 0)
    
    !name file, set plot device, and rotate the plot.
    !pdfcairo seems to give the best bounding box. LaTeX loves it.
    call plsfnam("wat.pdf")
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

  end subroutine plot_that_shit

!kthxbai  
end program raygun

