module get_config

contains

subroutine startfunc( tag, attribs, error )
  character(len=*)                 :: tag
  character(len=*), dimension(:,:) :: attribs
  logical                          :: error

  print *,tag
  !we will have to do some error/format checking...

end subroutine startfunc

subroutine datafunc( tag, data, error )
   character(len=*)               :: tag
   character(len=*), dimension(:) :: data
   logical                        :: error

   integer                        :: i

   print *,data
end subroutine datafunc
 
subroutine endfunc( tag, error )
   character(len=*)               :: tag
   logical                        :: error
   
   integer                        :: i

end subroutine endfunc

end module


program raygun

  use plplot
  use XMLPARSE
  use READ_XML_PRIMITIVES
  
  use get_config
  ! I guess we don't want this?
  !implicit none


  !see readint.f90 in the examples... I'll have to do somethign
  !similar.
  

  character(len=40), dimension(2,10) :: attribs
  character(len=80), dimension(100)  :: data
  logical                            :: error
  integer                            :: n
  character wat*100
  
  n = IARGC()
  print *,"indent?"
  print *,n

  call getarg(1, wat)
  print *,wat

  call xml_process( 'test.xml', attribs, data, startfunc, datafunc, endfunc, 0, error )
  
!  call read_xml_integer()
  
!  print *,"O HAI WORLDS!"
  
  stop
end program raygun

