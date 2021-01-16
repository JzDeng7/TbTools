program main2
  implicit none

  namelist /hopp/ a b c d e f
  open(8,'in.win','old')
  read(8,hopp)
  write(*,*) a, b, c, d, e, f
end program main2
  
