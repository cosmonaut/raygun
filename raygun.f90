
program raygun

  use plplot, PI => PL_PI
  use xml_data_config
  implicit none

  integer                            :: n
  character(len=100)                 :: file
  double precision                   :: i
  double precision, dimension(300)   :: rays = 0

  i = 1.34567834563987498374
  print *, i

  n = IARGC()

  if(n > 1) then
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
  
  print *, hibound

  !call xml_process( 'test.xml', attribs, data, startfunc, datafunc, endfunc, 0, error )
  
!  print *,"O HAI WORLDS!"
  
contains

  subroutine plot_that_shit()
    real(plflt), dimension(1:2) :: x, y
    real(plflt)                 :: xmin, xmax, ymin, ymax, rot
    integer                     :: just, axis

    xmin = lobound(1)
    ymin = lobound(2)
    xmax = hibound(1)
    ymax = hibound(2)
    just = 1
    axis = 0
    
    x(1) = 0
    x(2) = 5
    y(1) = 0
    y(2) = 5

    rot = 1

    !change color index 15 to be black
    call plscol0(15, 0, 0, 0)
    
    !name file, set plot device, and rotate the plot.
    !pdfcairo seems to give the best bounding box. LaTeX loves it.
    call plsfnam("wat.pdf")
    call plsdev("pdfcairo")


    !set background, initialize.
    call plscolbg(255, 255, 255)
    call plinit()
    call plfont(2)
    !call plsfont(1, 0, 0)
    call plfontld(1)

    call plcol0(15)
    call plenv(xmin, xmax, ymin, ymax, just, axis)

    call pllab("Blah #gx", "Meow", "#[0x00c5]#[0x212b]#[0x03c8]#[0x01fa]#[0x00b0]")
    call plcol0(1)
    call plline(x, y)

    call plend()
    
  end subroutine plot_that_shit

  
end program raygun

