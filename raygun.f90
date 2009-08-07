program raygun

  use plplot
  use xml_data_config
  implicit none

  integer                            :: n
  character wat*100

  n = IARGC()
  print *,"indent?"
  print *,n

  call getarg(1, wat)
  print *,wat

  
  call read_xml_file_config("test.xml")


  !call xml_process( 'test.xml', attribs, data, startfunc, datafunc, endfunc, 0, error )
  
!  print *,"O HAI WORLDS!"
  
  stop
end program raygun

