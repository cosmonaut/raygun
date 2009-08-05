program raygun

use plplot
use XMLPARSE

implicit none

character                              :: fname
type(XML_PARSE)                        :: info

fname = "test.xml"

call xml_open(info, fname, .true.)

print *,"O HAI WORLDS!"

stop
end program raygun
