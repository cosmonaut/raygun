module xml_data_config
  use READ_XML_PRIMITIVES
  use XMLPARSE
  implicit none
  integer, private :: lurep_
  logical, private :: strict_
  
  character(len=40), dimension(4)             :: optics
  integer, dimension(2,3)                     :: bounds
  integer, dimension(:), pointer              :: lobound
  integer, dimension(:), pointer              :: hibound
  integer                                     :: numrays, numoptics, beamstyle
  double precision, dimension(:), pointer     :: beamcenter, beamrot, par_pos, hyper_pos
  double precision                            :: beamradius, par_a, hyper_a, hyper_c
  double precision                            :: hyper_rad, par_rad, par_in_rad
  character(len=40)                           :: name
  double precision                            :: grat_r, det_r, det_rad, grat_lines, grat_rad
  double precision, dimension(:), pointer     :: det_pos, grat_pos

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

  logical :: has_optics
  logical :: has_lobound
  logical :: has_hibound
  logical :: has_numrays, has_beamradius, has_beamrot, has_beamcenter, has_beamstyle
  logical :: has_name, has_par_a, has_hyper_a, has_hyper_c, has_par_pos
  logical :: has_hyper_pos, has_hyper_rad, has_numoptics, has_par_rad, has_par_in_rad
  logical :: has_det_r, has_grat_r, has_det_pos, has_grat_pos, has_grat_lines, has_det_rad
  logical :: has_grat_rad


  has_lobound = .false.
  has_hibound = .false.
  has_numrays = .false.
  has_numoptics = .false.
  has_beamstyle = .false.
  has_beamradius = .false.
  has_beamcenter = .false.
  has_beamrot = .false.
  has_name = .false.
  has_par_pos = .false.
  has_par_a = .false.
  has_par_rad = .false.
  has_par_in_rad = .false.
  has_hyper_pos = .false.
  has_hyper_rad = .false.
  has_hyper_a = .false.
  has_hyper_c = .false.
  has_det_r = .false.
  has_det_rad = .false.
  has_grat_r = .false.
  has_grat_rad = .false.
  has_det_pos = .false.
  has_grat_pos = .false.
  has_grat_lines = .false.

  call xml_open( info, fname, .true. )
  call xml_options( info, report_errors=.true. )
  lurep_ = 0
  if ( present(lurep) ) then
     lurep_ = lurep
     call xml_options( info, report_lun=lurep )
  endif
  strict_ = .false.
  error = .false.
  
  !bounds(1,1) = 3
  !print *,bounds(1,1)

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

     case('raytrace')
        print *,"READING RAYTRACE XML"

     case('name')
        call read_xml_word( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             name, has_name )

     case('lowerbound')
        call read_xml_integer_array( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             lobound, has_lobound )

     case('upperbound')
        call read_xml_integer_array( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             hibound, has_hibound )

     case ('numrays')
        call read_xml_integer( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             numrays, has_numrays )

     case ('numoptics')
        call read_xml_integer( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             numoptics, has_numoptics )

     case ('beamstyle')
        call read_xml_integer( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             beamstyle, has_beamstyle )

     case ('beamrot')
        call read_xml_double_array( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             beamrot, has_beamrot )

     case ('beamradius')
        call read_xml_double( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             beamradius, has_beamradius )
        
     case ('beamcenter')
        call read_xml_double_array( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             beamcenter, has_beamcenter )

     case ('par_pos')
        call read_xml_double_array( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             par_pos, has_par_pos )

     case ('par_a')
        call read_xml_double( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             par_a, has_par_a )

     case ('par_rad')
        call read_xml_double( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             par_rad, has_par_rad )

     case ('par_in_rad')
        call read_xml_double( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             par_in_rad, has_par_in_rad )

     case ('hyper_pos')
        call read_xml_double_array( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             hyper_pos, has_hyper_pos )

     case ('hyper_rad')
        call read_xml_double( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             hyper_rad, has_hyper_rad )

     case ('hyper_a')
        call read_xml_double( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             hyper_a, has_hyper_a )

     case ('hyper_c')
        call read_xml_double( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             hyper_c, has_hyper_c )

     case ('grat_r')
        call read_xml_double( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             grat_r, has_grat_r )

     case ('grat_lines')
        call read_xml_double( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             grat_lines, has_grat_lines )

     case ('det_r')
        call read_xml_double( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             det_r, has_det_r )

     case ('det_rad')
        call read_xml_double( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             det_rad, has_det_rad )

     case ('grat_rad')
        call read_xml_double( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             grat_rad, has_grat_rad )

     case ('grat_pos')
        call read_xml_double_array( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             grat_pos, has_grat_pos )

     case ('det_pos')
        call read_xml_double_array( &
             info, tag, endtag, attribs, noattribs, data, nodata, &
             det_pos, has_det_pos )


     case default
        print *, "DEF!"
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
