   subroutine dump_cols(name, header, a)
      !! Write a column table to $VILMA_BENCH_DUMP/<name> for the analysis scripts.
      !!
      !! No-op unless VILMA_BENCH_DUMP names a directory, so `make check` and any
      !! plain run behave exactly as before — the dump is opt-in and costs
      !! nothing when off. `header` names the columns and is written as a leading
      !! `#` comment line, so the file is self-describing and readable with any
      !! delimited-text reader. Included by the benchmark tests that feed figures.
      character(*), intent(in) :: name, header
      real(wp),     intent(in) :: a(:,:)          ! (nrow, ncol)
      character(512) :: dir, path
      integer :: u, i, st
      call get_environment_variable('VILMA_BENCH_DUMP', dir, status=st)
      if (st /= 0 .or. len_trim(dir) == 0) return
      path = trim(dir)//'/'//name
      open(newunit=u, file=trim(path), status='replace', action='write', iostat=st)
      if (st /= 0) then
         write(*,'(3a)') '   WARNING: cannot write ', trim(path), ' (dump skipped)'
         return
      end if
      write(u,'(2a)') '# ', header
      ! Full double precision (17 significant digits). es18.10 rounded the
      ! values, which destroys structure a consumer may rely on -- the Spada disc
      ! reference is exactly piecewise linear, and rounding turned its ~40 true
      ! breakpoints into ~114 spurious ones.
      do i = 1, size(a,1)
         write(u,'(*(es24.16,1x))') a(i,:)
      end do
      close(u)
      write(*,'(3a,i0,a,i0,a)') '   dumped ', trim(path), ' (', size(a,1), ' x ', size(a,2), ')'
   end subroutine dump_cols
