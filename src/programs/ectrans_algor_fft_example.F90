program ectrans_algor_fft_example

use parkind1, only : jpim, jprb
use bluestein_mod, only : fftb_type, bluestein_init, bluestein_fft, bluestein_term

implicit none

external :: fft992, set99b

integer(kind=jpim), parameter :: n_fft992 = 12
integer(kind=jpim), parameter :: n_bluestein = 7
integer(kind=jpim), parameter :: batch = 3
real(kind=jprb), parameter :: tolerance = 500._jprb*epsilon(1._jprb)

call demo_fft992_batch(n_fft992, batch)
call demo_bluestein_batch(n_bluestein, batch)

contains

subroutine demo_fft992_batch(n, lot)

integer(kind=jpim), intent(in) :: n, lot

real(kind=jprb) :: trigs(n)
integer(kind=jpim) :: ifax(10)
logical :: lusefft992
real(kind=jprb), allocatable :: work(:,:)
real(kind=jprb), allocatable :: original(:,:)
real(kind=jprb) :: max_error

call set99b(trigs, ifax, n, lusefft992)
if (.not. lusefft992) then
  write(*,'(a,i0,a)') 'FFT992 does not support N=', n, ' in this configuration.'
  stop 1
endif

allocate(work(lot, n + 2))
allocate(original(lot, n))

call fill_batch_signal(original)
work(:,:) = 0._jprb
work(:,1:n) = original(:,:)

write(*,'(a)') 'FFT992 example'
write(*,'(a,i0,a,i0,a)') '  Supported length N=', n, ', LOT=', lot, ' batched vectors.'
write(*,'(a)') '  Layout: WORK(field,point), so INC=LOT and JUMP=1.'

call fft992(work, trigs, ifax, lot, 1_jpim, n, lot, -1_jpim)
call print_vector('  Packed spectral coefficients for field 1:', work(1,1:n+2))

call fft992(work, trigs, ifax, lot, 1_jpim, n, lot, 1_jpim)
max_error = maxval(abs(work(:,1:n) - original(:,:)))

write(*,'(a,1x,es12.4)') '  Max round-trip error:', max_error
if (max_error > tolerance) then
  write(*,'(a)') '  FFT992 round-trip check failed.'
  stop 1
endif

deallocate(original)
deallocate(work)

end subroutine demo_fft992_batch


subroutine demo_bluestein_batch(n, klot)

integer(kind=jpim), intent(in) :: n, klot

type(fftb_type) :: tb
integer(kind=jpim) :: iclen
real(kind=jprb), allocatable :: work(:,:)
real(kind=jprb), allocatable :: original(:,:)
real(kind=jprb) :: max_error

iclen = (n/2 + 1)*2

tb%ndlon = n
tb%nlat_count = 1
allocate(tb%nlats(1))
tb%nlats(1) = n
call bluestein_init(tb)

allocate(work(klot, iclen))
allocate(original(klot, n))

call fill_batch_signal(original)
work(:,:) = 0._jprb
work(:,1:n) = original(:,:)

write(*,'(a)') 'BLUESTEIN_FFT example'
write(*,'(a,i0,a,i0,a)') '  Arbitrary length N=', n, ', KLOT=', klot, ' batched vectors.'
write(*,'(a)') '  Layout: WORK(field,packed_point) with packed spectral size ICLEN=(N/2+1)*2.'

call bluestein_fft(tb, n, -1, klot, work)
work(:,1:iclen) = work(:,1:iclen) / real(n, jprb)
call print_vector('  Packed spectral coefficients for field 1:', work(1,1:iclen))

call bluestein_fft(tb, n, 1, klot, work)
max_error = maxval(abs(work(:,1:n) - original(:,:)))

write(*,'(a,1x,es12.4)') '  Max round-trip error:', max_error
if (max_error > tolerance) then
  write(*,'(a)') '  BLUESTEIN_FFT round-trip check failed.'
  stop 1
endif

call bluestein_term(tb)
deallocate(original)
deallocate(work)

end subroutine demo_bluestein_batch


subroutine fill_batch_signal(values)

real(kind=jprb), intent(out) :: values(:,:)

integer(kind=jpim) :: jf, jj, npoints
real(kind=jprb) :: angle, pi

npoints = size(values, 2)
pi = acos(-1._jprb)

do jf = 1, size(values, 1)
  do jj = 1, npoints
    angle = 2._jprb*pi*real(jj - 1, jprb)/real(npoints, jprb)
    values(jf,jj) = sin(real(jf, jprb)*angle) + 0.25_jprb*cos(real(jf + 1, jprb)*angle)
  enddo
enddo

end subroutine fill_batch_signal


subroutine print_vector(title, values)

character(len=*), intent(in) :: title
real(kind=jprb), intent(in) :: values(:)

integer(kind=jpim) :: j

write(*,'(a)') trim(title)
do j = 1, size(values)
  write(*,'(a,i0,a,1x,es18.10)') '    [', j, '] =', values(j)
enddo

end subroutine print_vector

end program ectrans_algor_fft_example