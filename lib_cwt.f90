!*|
!*|
!*|              Continuous wavelet transform library
!*|
!*|
!*|============================================================================
module lib_cwt
  use lib_spectral
  implicit none

  private cwt_constructor,cwt_transform,&
       make_wavelet, fact

  type cwt
     !private
     character(len=4), private                     :: Wt_Type
     integer, private                              :: Wt_order=4
     real   , private                              :: Wt_w0=6.0,fs
     real   , private                              :: Wt_mu=5,Wt_sigma=0.6
     real                                          :: ds,s0

     !public
     complex, dimension(:,:), allocatable          :: W
     integer                                       :: Nscales
     integer                                       :: Nsamples

   contains

     !public
     procedure :: create     => cwt_constructor
     procedure :: transform  => cwt_transform
     procedure :: scale2f    => cwt_scale2f
     procedure :: Parseval   => cwt_check_Pval

     !private
     procedure, private :: make_wavelet

  end type cwt

contains

  subroutine cwt_constructor(this,Ns,samp_freq,wt_order,wt_w0,wt_mu,wt_sigma)
    class(cwt)                       ::this
    integer                          ::Ns
    real                             ::samp_freq
    integer      ,optional           ::wt_order
    real         ,optional           ::wt_w0
    real,        optional            ::wt_mu,wt_sigma    !for Morlet
    !========================================

    this%Nsamples = Ns
    this%fs       = samp_freq

    select case (this%wt_type)
    case ('Paul')
       this%ds = 0.02d0

    case ('Morl')
       this%ds = 0.04d0

    case ('Bump')
       this%ds = 0.015d0
    end select

    this%s0 = 2.d0/this%fs
    ! this%Nscales = int( log(this%Nsamples/this%fs/this%s0/50000.)&
    !      /log(2.0)/this%ds )
    this%Nscales = 250

    print*, "Nb Scale, St_min, St_max"
    print*,this%nscales, this%scale2f(this%nscales)*0.05/340., this%scale2f(1)*0.05/340

    if (present(wt_order)) this%wt_order = wt_order
    if (present(wt_w0)   ) this%wt_w0    = wt_w0
    if (present(wt_sigma)) this%wt_sigma = wt_sigma
    if (present(wt_sigma)) this%wt_mu    = wt_mu

    if (.not.allocated(this%W)) then
       allocate (this%W(this%Nscales,this%Nsamples))
    end if

  end subroutine cwt_constructor

  !****************************************************************
  !*       Continuous Wavelet Transform
  !****************************************************************

  subroutine cwt_transform(this,sig,samp_freq,wt_type,order,w0,mu,sigma)
    class(cwt)                         ::this
    real,       dimension(:)           ::sig
    real                               ::samp_freq
    integer                            ::Nsc
    character(len=4)                   ::wt_type
    integer,     optional              ::order    !for Paul or DOG
    real,        optional              ::w0       !for Morlet
    real,        optional              ::mu,sigma !for bump

    !internal variables
    type(psd_param)                    ::fft_param
    complex, dimension(:), allocatable ::sigfft,W_hat,wavelet
    real                               ::scale,s0,ds
    integer                            ::is,iw
    !==============================================================

    !Construct the wavelet coef table
    this%wt_type = wt_type
    select case(this%wt_type)
    case ('Paul')
       if (.not.present(order)) STOP 'missing arguement in cwt_transform'
       call this%create(size(sig,1),samp_freq,wt_order=order)
    case ('Morl')
       if (.not.present(w0)) STOP 'missing arguement in cwt_transform'
       call this%create(size(sig,1),samp_freq,wt_w0=w0)
    case('Bump')
       if (.not.present(mu)) STOP 'missing arguement in cwt_transform'
       if (.not.present(sigma)) STOP 'missing arguement in cwt_transform'
       call this%create(size(sig,1),samp_freq,wt_mu=mu, wt_sigma=sigma)
    end select

    !table allocation
    allocate(wavelet(this%Nsamples))
    allocate(W_hat(this%Nsamples))

    !fft of the whole signal
    allocate(sigfft(this%nsamples/2+1))
    fft_param = psd_param (nfft=this%nsamples, norm_fft=.false.)
    call fft(sig,sigfft,fft_param)

    do is=1,this%Nscales
       !build the wavelet in the frequency domain
       !or compute the FT of a time domain wavelet
       scale = this%s0*2**((is-1)*this%ds)
       call this%make_wavelet(scale,wavelet)

       ! convolution of the scaled wavelet with the signal
       do iw=1,this%Nsamples/2
          W_hat(iw)=sigfft(this%nsamples/2-iw+2)*conjg(wavelet(iw))
       end do
       do iw=this%Nsamples/2+1,this%Nsamples
          W_hat(iw)=sigfft(iw-this%Nsamples/2)*conjg(wavelet(iw))
       end do
       !inverse Fourier transform to get the wavelet coefficient
       call ifft(W_hat,this%W(is,:),fft_param)

    end do

    !fft scaling
    this%W=this%W/this%nsamples

    if (allocated(wavelet)) deallocate(wavelet)
    if (allocated(w_hat))   deallocate(W_hat)

  end subroutine cwt_transform

  !****************************************************************
  !*       Wavelet maker
  !****************************************************************
  subroutine make_wavelet(this,scale,wavelet)
    implicit none
    class(cwt)                     ::this
    complex, dimension(:)          ::wavelet
    real                           ::scale

    !internal variables
    integer                        ::m,iw
    real                           ::omega,w0
    real, parameter                ::pi=4.0*atan(1.0)
    !===================================================

    select case (this%wt_type)
    case ('Paul')
       m=this%wt_order

       !filling the wavelet table
       wavelet(1:this%Nsamples/2)=0.d0
       do iw=this%Nsamples/2+1,this%Nsamples
          omega = 2*pi*real(iw-this%Nsamples/2)/real(this%Nsamples)*this%fs
          wavelet(iw) = cmplx(2**m/sqrt(real(m*fact(2*m-1)))*&
               (scale*omega)**m*exp(-scale*omega) , 0.0)
       end do


    case ('Morl')
       w0=this%wt_w0

       !filling the wavelet table
       wavelet(1:this%Nsamples/2)=0.d0
       do iw=this%Nsamples/2+1,this%Nsamples
          omega = 2*pi*real(iw-this%Nsamples/2)/real(this%Nsamples)*this%fs
          wavelet(iw) = cmplx(pi**(-0.25d0)*exp(-(scale*omega - w0)**2/2.d0) , 0.0)
       end do

    case ('Bump')
       w0=this%wt_w0

       !filling the wavelet table
       wavelet=cmplx(0.0,0.0)
       do iw=1,this%Nsamples
          omega = 2*pi*real(iw-this%Nsamples/2)/real(this%Nsamples)*this%fs
          if ( scale*omega >= (this%wt_mu-this%wt_sigma) .and. &
               scale*omega <= (this%wt_mu+this%wt_sigma) ) then
             wavelet(iw) = cmplx(exp(1.0 - 1.0/&
                  (1.0-(scale*omega-this%wt_mu)**2/this%wt_sigma**2)),0.d0)
          end if
       end do

    case default
       STOP 'Wavelet type unknown'
    end select

    !scaling
    wavelet = sqrt(2*pi*scale*this%fs)*wavelet


  end subroutine make_wavelet

  !****************************************************************
  !*       Check Parseval
  !****************************************************************
  function cwt_check_Pval(this) result(rms)
    class(cwt)                         ::this
    real                               ::rms

    integer                            ::is,i
    real                               ::scale
    !------------------------------------------

    rms=0.0
    do i=1,this%nsamples
       do is=1,this%nscales
          scale = this%s0*2**((is-1)*this%ds)
          rms = rms + abs(this%W(is,i))**2/scale
       end do
    end do

    rms = rms*this%ds/this%fs/real(this%nsamples)

    select case (this%wt_type)
    case ('Paul')
       rms = rms / 1.132
    case ('Morl')
       rms = rms / 0.776
    case default
       write(06,*) 'cwt_check_Pval: unknown coefficients for this wavelet'
    end select

  end function cwt_check_Pval

  !****************************************************************
  !*       factorial
  !****************************************************************
  recursive function fact(n) result(f)
    implicit none
    integer                        ::n
    integer                        ::f
    !=======================================

    if (n<=1) then
       f=1
    else
       f=n*fact(n-1)
    end if

  end function fact


  !****************************************************************
  !*       Wavelet scale to freq
  !****************************************************************
  function cwt_scale2f(this,iscale) result(freq)
    class(cwt)                         ::this
    integer                            ::iscale
    real                               ::freq

    real                               ::scale
    real, parameter                    ::pi=4.0*atan(1.0)
    !==============================================================

    scale = this%s0*2**(iscale*this%ds)

    select case (this%wt_type)
    case('Paul')
       freq=(2*this%wt_order + 1)/(4*pi*scale)
    case('Morl')
       freq=(this%wt_w0 + sqrt(2+this%wt_w0**2))/(4*pi*scale)
    case('Bump')
       freq=this%wt_mu/(2.0*pi*scale)
    case default
       STOP 'scale2f: unknwown wavelet'
    end select


  end function cwt_scale2f

end module lib_cwt
