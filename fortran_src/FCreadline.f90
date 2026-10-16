!-------------------------------------------------------------------------------
module jmk_readline_module
   use ISO_C_Binding
   implicit none
   private
   public readline,save_history_list,get_history_list
!-------------------------------------------------------------------------------
! The structure used to store a history entry.
  type, bind(C) :: hist_entry_c
    type(C_ptr) :: line, timestamp, data
  end type hist_entry_c

! A structure used to pass the current state of the history stuff around.
  type, bind(C) :: hist_state_c
    type(C_ptr) :: entries;     ! (hist_entry**)  Pointer to the entries themselves.
    integer(C_int) :: offset    ! The location pointer within this array.
    integer(C_int) :: length    ! Number of elements within this array.
    integer(C_int) :: size      ! Number of slots allocated to this array.
    integer(C_int) :: flags
  end type hist_state_c
!-------------------------------------------------------------------------------
  interface
    function readline_c(prompt) result(line_ptr) bind(C,name="readline")
      import
      type(C_ptr) :: line_ptr
      character(kind=C_char), intent(in) :: prompt(*)
    end function readline_c
    subroutine using_history() bind(C)
      import
    end subroutine using_history
    subroutine add_history_c(line_ptr) bind(C,name="add_history")
      import
      type(C_ptr), value, intent(in) :: line_ptr
    end subroutine add_history_c
    function history_get_history_state() result(state_ptr) bind(C)
      import
      type(C_ptr) :: state_ptr
    end function history_get_history_state
    subroutine free_c(ptr) bind(C,name="free")
      import
      type(C_ptr), value, intent(in) :: ptr
    end subroutine free_c
    function strlen(s) bind(C)
      import
      integer(C_size_t) :: strlen
      type(C_ptr), value, intent(in)  :: s
    end function strlen
  end interface
!-------------------------------------------------------------------------------
contains

  subroutine readline(line,prompt)
    implicit none
    character(len=*), intent(out) :: line
    character(len=*), intent(in)  :: prompt
  ! local
    type(C_ptr) :: line_ptr
    character(kind=C_char), pointer :: char_ptr(:)
    integer i,linenumber,ierr
    logical found
    character(len=256) test
  ! begin
    call using_history()
  ! trim to last non-blank character and append null for C, and
  ! call readline(3c) to read a line of input in edit mode
    found=.false.
    do while(.not.found)
      line_ptr = readline_c(trim(prompt)//' '//C_NULL_char)
  ! Copy the result to 'line'
      line=' '
      call C_F_pointer(line_ptr,char_ptr,(/len(line)/))
      do i=1,len(line)
        if (iachar(char_ptr(i))==0) exit
!        if (iachar(char_ptr(i)).lt.32) exit
        line(i:i)=char_ptr(i)
      end do
  ! if the "history" command is on a line by itself, show history
      if (line=='h') then 
         call show_history_list()
      else if (line(1:1).eq.'!') then
         !Try to read the line number requested
         read(line(2:len_trim(line)),*) linenumber
         line=get_history_line(linenumber,found)
         if(found) then
            write(6,'(/a,a,a)',advance='no') trim(prompt)//' ',trim(line),' ? (Press <enter> to accept.)'
            test=' '
            read(5,'(a)',iostat=ierr) test
            found=len_trim(test).eq.0
         end if
      else
         found=.true.
         if(len_trim(line).gt.0) call add_history_c(line_ptr) !Don't add blank lines to history, but return them to calling program
      end if
      if(.not.found) write(6,*) ! Blank line to fit with EPSR format
    end do
  ! free allocated result line returned from readline(3c)
    call free_c(line_ptr)
  end subroutine readline

  subroutine show_history_list()
    type(C_ptr) :: ptr
    type(hist_state_c), pointer :: state_ptr
    type(C_ptr), pointer :: entries(:)
    type(hist_entry_c), pointer :: entry
    character(kind=C_char,len=huge(0)), pointer :: line
    integer :: i
    integer(kind=8) :: line_len
    ptr = history_get_history_state()
    call C_F_pointer(ptr,state_ptr)
    call C_F_pointer(state_ptr%entries,entries,(/state_ptr%length/))
    write(6,'(/a/)') 'Current history list:'
    do i=1,size(entries,1)
      call C_F_pointer(entries(i),entry)
      call C_F_pointer(entry%line,line)
      line_len=strlen(entry%line)
!      write(*,'(I0,4A)') i,": '",line(:line_len),"'",merge('*',' ',i==state_ptr%offset)
      write(6,'(i6,1x,a)') i,line(:line_len)
    end do
    write(6,'(/a)') 'To use a particular line type !<line number>'
  end subroutine show_history_list

  function get_history_line(linenumber,found)
    character(len=256) :: get_history_line
    integer, intent(in) :: linenumber
    logical, intent(out) :: found
    type(C_ptr) :: ptr
    type(hist_state_c), pointer :: state_ptr
    type(C_ptr), pointer :: entries(:)
    type(hist_entry_c), pointer :: entry
    character(kind=C_char,len=huge(0)), pointer :: line
    integer :: i
    integer(kind=8) :: line_len
    ptr = history_get_history_state()
    call C_F_pointer(ptr,state_ptr)
    call C_F_pointer(state_ptr%entries,entries,(/state_ptr%length/))
    found=linenumber.gt.0.and.linenumber.le.size(entries,1)
    if(found) then
      i=linenumber
      call C_F_pointer(entries(i),entry)
      call C_F_pointer(entry%line,line)
      line_len=strlen(entry%line)
      get_history_line=line(1:line_len)
      get_history_line=get_history_line(1:len_trim(get_history_line))
    else
      get_history_line=''
    end if
  end function get_history_line

  subroutine current_history_list(currenthistory,ncurrenthistory)
    type(C_ptr) :: ptr
    type(hist_state_c), pointer :: state_ptr
    type(C_ptr), pointer :: entries(:)
    type(hist_entry_c), pointer :: entry
    character(kind=C_char,len=huge(0)), pointer :: line
    character(len=256), dimension(:), allocatable, intent(out) :: currenthistory
    integer, intent(out) :: ncurrenthistory
    integer :: i
    integer(kind=8) :: line_len
    ptr = history_get_history_state()
    call C_F_pointer(ptr,state_ptr)
    call C_F_pointer(state_ptr%entries,entries,(/state_ptr%length/))
    if(allocated(currenthistory)) deallocate(currenthistory)
    if(size(entries,1).gt.0) allocate(currenthistory(size(entries,1)))
    i=0
    do while(i.lt.size(entries,1))
      i=i+1
      call C_F_pointer(entries(i),entry)
      call C_F_pointer(entry%line,line)
      line_len=strlen(entry%line)
      currenthistory(i)=line(1:line_len)
    end do
    ncurrenthistory=i
!    write(6,'(a,i0)') 'current_history_list> No. of entries: ',ncurrenthistory
  end subroutine current_history_list

  subroutine save_history_list()
    type(C_ptr) :: ptr
    type(hist_state_c), pointer :: state_ptr
    type(C_ptr), pointer :: entries(:)
    type(hist_entry_c), pointer :: entry
    character(kind=C_char,len=huge(0)), pointer :: line
    character(len=256), dimension(:), allocatable :: myhistory
    character(len=256) :: testline
    integer :: i,j,ierr,nmyhistory
    integer(kind=8) :: line_len
    logical found
    ptr = history_get_history_state()
    call C_F_pointer(ptr,state_ptr)
    call C_F_pointer(state_ptr%entries,entries,(/state_ptr%length/))
    if(allocated(myhistory)) deallocate(myhistory)
    allocate(myhistory(size(entries,1)))
    nmyhistory=0
! Check each entry does not already exist in the history
    do i=1,size(entries,1)
      call C_F_pointer(entries(i),entry)
      call C_F_pointer(entry%line,line)
      line_len=strlen(entry%line)
      testline=line(:min(255,line_len))
      j=0
      found=.false.
      do while(j.lt.nmyhistory.and..not.found)
        j=j+1
        found=testline(1:len_trim(testline)).eq.myhistory(j)(1:len_trim(myhistory(j))) &
        .or.len_trim(testline).eq.0
      end do
      if(.not.found) then
         nmyhistory=nmyhistory+1
         myhistory(nmyhistory)=testline(1:len_trim(testline))
      end if
    end do
    open(10,file='FChistory.txt',status='unknown',iostat=ierr)
    if(ierr.eq.0) then
      do j=1,nmyhistory
        write(10,'(a)') myhistory(j)(1:len_trim(myhistory(j)))
      end do
    end if
    close(10)
    if(allocated(myhistory)) deallocate(myhistory)
  end subroutine save_history_list

  subroutine get_history_list()
    implicit none
    character(len=256) :: line
  ! local
    type(C_ptr) :: line_ptr
    character(kind=C_char), pointer :: char_ptr(:)
    character(len=1), target :: my_char(256)
    character(len=256), dimension(:), allocatable :: currenthistory
    character(len=256) :: myhistory
    integer :: i,j,ierr,ncurrenthistory
    logical :: found
  ! begin
    call using_history()
    !Get the current history
    call current_history_list(currenthistory,ncurrenthistory)
    open(10,file='FChistory.txt',status='unknown',iostat=ierr)
    do while(ierr.eq.0)
      read(10,'(a)',iostat=ierr) myhistory
      if(ierr.eq.0) then
! Check whether this line exists in the current history
         found=.false.
         i=0
         do while(.not.found.and.i.lt.ncurrenthistory)
             i=i+1
             found=myhistory(1:len_trim(myhistory)).eq.currenthistory(i)(1:len_trim(currenthistory(i)))
         end do
         if(.not.found) then
           line=myhistory(1:len_trim(myhistory))
           do i=1,min(len_trim(line),255)
             my_char(i)=line(i:i)
           end do
           my_char(len_trim(line)+1)=C_NULL_char
           line_ptr=c_loc(my_char(1))
           call C_F_pointer(line_ptr,char_ptr,(/len_trim(line)+1/))
           do i=1,size(char_ptr)
             char_ptr(i)=my_char(i)
           end do
           call add_history_c(line_ptr)
         end if
      end if
    end do
    close(10)
    if(allocated(currenthistory)) deallocate(currenthistory)
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
