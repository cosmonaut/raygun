program raygun

  use plplot
  use xml_data_config
  implicit none

  integer                            :: n
  character(len=100)                 :: file
  double precision                   :: i

  i = 1.34567834563987498374D0
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

  

  !call xml_process( 'test.xml', attribs, data, startfunc, datafunc, endfunc, 0, error )
  
!  print *,"O HAI WORLDS!"
  
  stop
end program raygun

