!-------------------------------------------------------------------------------
module jmk_readline_module
   implicit none
   private
   public readline,save_history_list,get_history_list
!-------------------------------------------------------------------------------
contains

  subroutine readline(line,prompt)
    implicit none
    character(len=*), intent(out) :: line
    character(len=*), intent(in)  :: prompt
    write(*,'(a,1x)',advance='no') prompt(1:len_trim(prompt))
    read(*,'(a)') line
  end subroutine readline

  subroutine show_history_list()
  end subroutine show_history_list

  subroutine save_history_list()
  end subroutine save_history_list

  subroutine get_history_list()
  end subroutine get_history_list

end module jmk_readline_module

! the test program
!program test_program
!  use jmk_readline_module
!  character(len=256):: line
! begin
!  call get_history_list()
!  write(*,*) &
!      '____________________________________________________________', &
!      '  Your input lines are now editable using the GNU', &
!      '  readline(3C) procedure.  By default, up-arrow and', &
!      '  down-arrow go thru the history lines; left and right arrow', &
!      '  keys and delete and just typing characters let you do', &
!      '  simple editing. Far more input control is available.', &
!      '  See the browser pages and man(1) pages for readline(3c).', &
!      ' ____________________________________________________________', &
!      ' Enter text and then edit it. "q" quits; "h" displays history:'
!  do
!    call readline(line,'EPSRshell>') ! read editable input line
!    if (line=='q') exit
!  enddo
!  call save_history_list()
!end program test_program
!
