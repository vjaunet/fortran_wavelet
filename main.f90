PROGRAM CWT_PRGM
  use lib_press_data
  use lib_cwt
  implicit none

  character(len=100)                        ::ifname,ofname,ofname_bin
  character(len=100)                        ::arg
  integer                                   ::iarg=1

  type(pressdata)                           ::micro
  type(cwt)                                 ::wavelet
  character(len=4)                          ::wvt_type='Paul'
  integer                                   ::order=4
  real                                      ::wt_w0=6.0
  real                                      ::wt_mu=5.0
  real                                      ::wt_sigma=0.6

  real                                      ::mean
  integer                                   ::imicro
  integer                                   ::nsamp_max=0,unsamp
  real, dimension(:,:)     ,allocatable     ::M_gnuplot

  !======
  real, parameter                           ::D=0.05
  real                                      ::Uj

  !======
  integer                                   ::i,is
  real, parameter                           ::pi=4.d0*atan(1.d0)
  !==============================================================


  !=======================================
  ! command line handling
  !---------------------------------------

  if (command_argument_count() < 1) then
     write(06,'(a)')"missing argument"
     call print_help();
  end if

  do while (iarg <= command_argument_count())

     call get_command_argument(iarg, arg)

     select case (trim(arg))

     case("-i")
        iarg=iarg+1
        call get_command_argument(iarg, ifname)

     case("-o")
        iarg=iarg+1
        call get_command_argument(iarg, ofname)

     case("-nsamp")
        iarg=iarg+1
        call get_command_argument(iarg, arg)
        read(arg,*)Nsamp_max

     case("-imic")
        iarg=iarg+1
        call get_command_argument(iarg, arg)
        read(arg,*)imicro

     case("-type")
        iarg=iarg+1
        call get_command_argument(iarg, wvt_type)

     case("-order")
        iarg=iarg+1
        call get_command_argument(iarg,arg)
        read(arg,*)order

     case("-w0")
        iarg=iarg+1
        call get_command_argument(iarg,arg)
        read(arg,*)wt_w0

     case("-mu")
        iarg=iarg+1
        call get_command_argument(iarg,arg)
        read(arg,*)wt_mu

     case("-sigma")
        iarg=iarg+1
        call get_command_argument(iarg,arg)
        read(arg,*)wt_sigma

     case default
        write(06,*)' Error : Unknown option ',trim(arg)
        call print_help()

     end select

     iarg=iarg+1
  end do



  !=======================================
  ! reading data
  !---------------------------------------
  call read_ccur_data(micro,ifname,Uj)

  !=======================================
  ! wavelet transform
  !---------------------------------------
  UnSamp = 10
  mean=sum(micro%p(:,imicro))/real(micro%nsamples)
  micro%p(:,imicro) = micro%p(:,imicro)-mean
  call wavelet%transform(micro%p(1:nsamp_max:UnSamp,imicro),&
       micro%fs/unsamp,wvt_type,&
       order=order,w0=wt_w0,mu=wt_mu,sigma=wt_sigma)
  call micro%destroy()

  print*,"Max Conv. time:", real(wavelet%Nsamples/(micro%fs/UnSamp)*Uj/D)
  !print*,wavelet%Nsamples,nsamp_max

  !=======================================
  ! outputing results
  !---------------------------------------
  open (unit=10, file=trim(ofname))
  write(10,'(a)')'# St, abs(C_w)'
  do is=1,wavelet%Nscales
     write(10,*)wavelet%scale2f(is)*D/Uj,sum(abs(wavelet%W(is,:)))/real(wavelet%nsamples)
  end do
  close(10)

  !in gnuplot binary format
  allocate(M_gnuplot(wavelet%Nsamples+1,wavelet%Nscales+1))
  M_gnuplot(1,1) = wavelet%Nsamples
  M_gnuplot(1,2:wavelet%Nscales+1) = (/(real(wavelet%scale2f(is)*D/Uj),is=1,wavelet%Nscales)/)
  M_gnuplot(2:wavelet%Nsamples+1,1) = (/(real((i-1)/(micro%fs/UnSamp)*Uj/D),&
       i=1,wavelet%Nsamples,1)/)
  M_gnuplot(2:,2:) = real(abs(transpose(wavelet%W(:,1:wavelet%Nsamples:1))))

  ofname_bin=ofname(1:len_trim(ofname)-3)//'bin'
  open (unit=10, file=trim(ofname_bin), action='write', access='stream',&
       form='unformatted',status='replace')
  write(10)M_gnuplot
  close(10)
  deallocate(M_gnuplot)

  !! free memory
  call micro%destroy()

contains

  subroutine read_ccur_data(pdata,filename,Uj)
    type(pressdata)                  ::pdata
    character(len=*)                 ::filename
    real                             ::Uj

    !ccur bin file data variable
    character(len=3)                 ::version
    character(len=17)                ::date_time
    integer(kind=8)                  ::n_samples,n_micros
    integer(kind=8)                  ::freq
    real                             ::M,p_amb,t_amb,p0_p,p0_s
    real                             ::T0_p,T0_s,hydro
    character(len=500)               ::comments
    real, dimension(:,:), allocatable::coord_micros
    real, dimension(:,:), allocatable::p_acou
    !=========================================================

    !read ccur bin file data
    open(unit=111, file=trim(filename), access='stream')
    read(111)version
    read(111)date_time
    read(111)n_samples,n_micros,freq
    read(111)M,p_amb,t_amb,p0_p,p0_s
    read(111)T0_p,T0_s,hydro
    read(111)comments
    allocate(coord_micros(n_micros,6))
    read(111)coord_micros
    allocate(p_acou(n_samples, n_micros))
    read(111)p_acou
    close(111)

    !! filling up the pressdata handler
    call pdata%create(int(n_samples),int(n_micros),real(freq),real(p_amb),comments)
    pdata%p=p_acou
    pdata%x=0.0
    !compute the jet flow velocity
    Uj = M*sqrt(1.4*287.05*(273.15 + t_amb))

    deallocate(p_acou)
    deallocate(coord_micros)

  end subroutine read_ccur_data

  subroutine print_help()

    write(06, '(a)')'HELP :'
    write(06, '(a)')'./cwttr -i file -o ofile -nsc 10 -type Paul -order 6 -w0 6.0'
    STOP 'execution stopped'

  end subroutine print_help


END PROGRAM CWT_PRGM
