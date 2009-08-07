module xml_data_config
  use READ_XML_PRIMITIVES
  use XMLPARSE
  implicit none
  integer, private :: lurep_
  logical, private :: strict_
  
  character(len=40), dimension(4)      :: optics
  logical                              :: has_optics
  integer, dimension(2,3)              :: bounds
  logical                              :: has_bounds
  
  integer                              :: numrays
  logical                              :: has_numrays

  character        :: hi*30
  logical          :: has_hi

contains

subroutine read_xml_file_config(fname, lurep, errout)
  character(len=*), intent(in)           :: fname
  integer, intent(in), optional          :: lurep
  logical, intent(out), optional         :: errout

  type(XML_PARSE)                        :: info
  logical                                :: error
  character(len=80)                      :: tag
  logical                                :: endtag
  character(len=80), dimension(1:2,1:20) :: attribs
  integer                                :: noattribs
  character(len=200), dimension(1:100)   :: data
  integer                                :: nodata

  has_hi = .false.
  

  call xml_open( info, fname, .true. )
  call xml_options( info, report_errors=.true. )
  lurep_ = 0
  if ( present(lurep) ) then
     lurep_ = lurep
     call xml_options( info, report_lun=lurep )
  endif
  strict_ = .false.
  error = .false.
  
  print *,"ERROR?"
  print *,xml_error(info)
  
  bounds(0,2) = 3
  print *,bounds(0,2)

  do
     call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
     if ( xml_error(info) ) then
        write(*,*) 'Error reading input file!'
        stop
     endif
     if ( endtag .and. noattribs .eq. 0 ) then
        if ( xml_ok(info) ) then
           cycle
        else
           exit
        endif
     endif
     select case( tag )
     case('hi')
        print *,tag
        call read_xml_line( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             hi, has_hi )
        print *, hi
        
     case default
        if ( strict_ ) then
           error = .true.
           if ( lurep_ .gt. 0 ) then
              write( lurep_, * ) 'Unknown or wrongly placed tag: ',  trim(tag)
           else
              write( *, * ) 'Unknown or wrongly placed tag: ',  trim(tag)
           endif
        endif
     end select
     
     if ( .not. xml_ok(info) ) exit
  end do


end subroutine read_xml_file_config

end module xml_data_config
