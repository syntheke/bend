program bend
  
  !***********************************************************************
  ! Bennding of non-positive matrix to make it positive definite 
  ! 
  ! Optimised version of the bending subroutine in Sang Hong Lee's
  !   MTG2 program (https://sites.google.com/site/honglee0707/
  !
  ! Author: Gerhard Moser
  !*********************************************************************

  IMPLICIT NONE

  integer :: i,j, narg, thn, mg, n, io
  real,allocatable::GRM(:,:,:)
  real :: x1
  character(len=64),allocatable::grm_known(:)
  character(len=64)::fl1,fl2,filnam,fl3, option, cdum4

  thn=1
  
  narg=command_argument_count()
  if (narg==0.or.narg==1) then
     print*,'-p fam file -g grm file bend 1 '
     stop
  end if

  do i=1, (narg/2)
     call get_command_argument(i*2-1, option)
     call get_command_argument(i*2, filnam)
     if (option =='-p'.or.option=='-fam') then
        fl1=filnam
        print*,'ID file   : ',trim(fl1)
     else if (option=='-g'.or.option=='-grm') then
        mg=0
        fl3=filnam
        print*,'grm file  : ',trim(fl3)
     else if(option == '-threads') then
        read(filnam,*)thn
        print*,'threads   : ',thn
     else
        print*,'-p fam file -g grm file -threads 1 '
        stop
     end if
  end do

  n=0
  open (unit=38,file=trim(fl1),status='old')      !fam file
  do 
     n=n+1
     read(38,*,iostat=io)cdum4
     if (io.ne.0) exit
     !print*,trim(pedor(i,1)),trim(pedor(i,2))
  end do
  close(38)
  n=n-1
  print*,'no. ID: ',n

  allocate(GRM(1,n,n))
  GRM=0.0
  open (unit=40,file=trim(fl3),status='old')
  do
     read (40,*,iostat=io)i,j,x1 
     if (io.ne.0) exit
     GRM(1,i,j)=x1
     GRM(1,j,i)=x1
  end do
  close(40)


  call OMP_SET_NUM_THREADS(thn)

  call mtg_bend (GRM,n,n,fl3,thn)
  
contains

  subroutine mtg_bend (mbin,n,vn,grm_known,thn)
    
    !***********************************************************************
    !PCA analysis 
    !Author: Sang Hong Lee (c) 2009,2015
    !***********************************************************************
    
    implicit none
    character(len=64)::grm_known,cdum
    integer::n,i,j,k,k2,zi,zj,im,ma1,info,vn,thn,M,np,p,OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
    integer::lwork,liwork,IL,IU

    real::mbin(1,n,n)          !assumin single GRM

    double precision,allocatable::eivec(:,:),eival(:),Z(:,:),v2(:,:)
    double precision,allocatable::work(:)
    integer,allocatable::iwork(:),ISUPPZ(:)
    double precision::x1,VL,VU,v1

    !time measure **********************************************************
    INTEGER::now(8)
    CHARACTER*8 date
    CHARACTER*20 time,zone
    double precision:: timef, t1_real, t2_real, elapsed_real_secs
    double precision:: t1_cpu, t2_cpu, elapsed_cpu_secs
    !***********************************************************************

    allocate(eival(n),eivec(n,n),Z(n,n)) 
    allocate(ISUPPZ(2*n)) 

    PRINT*,'Bending GRM to be PD *********************************************'
    print*,''

    eivec=mbin(1,:,:)

    call cpu_time (t2_cpu)

    IL=1
    IU=vn

    lwork=-1
    liwork=-1
    allocate(work(1)) 
    allocate(iwork(1)) 

    call dsyevr ('V','A','U',n,eivec,n,VL,VU,IL,IU,-1.0,M,eival,Z,n,ISUPPZ,work,lwork,iwork,liwork,info)

    lwork=int(abs(work(1)))
    deallocate(work)
    allocate(work(lwork))

    liwork=int(abs(iwork(1)))
    deallocate(iwork)
    allocate(iwork(liwork))

101 continue

    call dsyevr ('V','A','U',n,eivec,n,VL,VU,IL,IU,-1.0,M,eival,Z,n,ISUPPZ,work,lwork,iwork,liwork,info)

    print*,"0 means no error >>> ",info
    print*,''

    k=0
    do i=1,n
       if (eival(i).le.0.0D0) then
          eival(i)=0.001
          k=k+1
       end if
    end do
    print*,k

    if (k==0) goto 110
    
    allocate(v2(n,thn))
!$OMP PARALLEL DO PRIVATE(i,j)
    do p=1,thn
       do i=p,n,thn
          v2(:,p) = Z(i,:) * eival(:)
          do j=1, i
             eivec(j,i) =  dot_product(v2(:,p), Z(j,:))
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO
    deallocate(v2)

    print*,'bent and check again >>>'
    goto 101

110 continue

    allocate(v2(n,thn))
!$OMP PARALLEL DO PRIVATE(i,j)
    do p=1,thn
       do i=p,n,thn
          v2(:,p) = Z(i,:) * eival(:)
          do j=1, i
             eivec(j,i) =  dot_product(v2(:,p), Z(j,:))
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO
    deallocate(v2)

    call cpu_time (t1_cpu)
    print*,'bending done - time:',real(t1_cpu-t2_cpu)

    open (UNIT=39,FILE=trim(grm_known)//'.bend',STATUS='unknown')
    do i=1,n
       do j=1,i
          if (eivec(i,j).ne.0) then
             !write(39,*)i,j,real(eivec(i,j))  !,GRM(i,j)
             write(39,*)i,j,real(eivec(j,i))  !,GRM(i,j)
          end if
       end do
    end do
    close(39)

  end subroutine


end program bend
