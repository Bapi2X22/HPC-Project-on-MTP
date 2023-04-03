subroutine force_calc(TotAtom,Box,Rcut,r,Sig,Eps,Force,PE)
 use general, only: dp,atom1,atom2
 implicit none 
 integer,intent(in) :: TotAtom
 integer :: i,ierr
 real(kind=dp),intent(in) :: r(TotAtom,3)
 real(kind=dp) :: fac2,fac6,r2,df,fc(3),dr(3),Ecut,R2cut,num_neigh(TotAtom),neighbours(TotAtom,TotAtom)
 real(kind=dp),intent(out) :: Force(TotAtom,3),PE

 real(kind=dp) :: Box,Rcut,Eps,Sig

 R2cut=Rcut*Rcut 

 ! calculating Ecut at Rcut 
 fac2=Sig*Sig/R2cut
 fac6=fac2*fac2*fac2
 Ecut=4.d0*Eps*fac6*(fac6-1)

 ! calculting force and energy 
  open(unit=20,file='verlet_list.xyz',action='read',iostat=ierr)
  open(unit=30,file='num_verlet.xyz',action='read',iostat=ierr)
  
  num_neigh=0
  neighbours=0
  
  do atom1=1,TotAtom
   read(30,*) num_neigh(atom1)
  enddo
  
  do atom1=1,TotAtom
   i=num_neigh(atom1)
   read(20,*) neighbours(atom1,1:i) 
  enddo

 PE=0.d0
 Force=0.d0
 do atom1=1,TotAtom
  do i=1,num_neigh(atom1)
     atom2=neighbours(atom1,i)
     dr=r(atom1,:)-r(atom2,:)
     dr=dr-Box*anint(dr/Box)
     r2=dot_product(dr,dr)
     if(r2<=R2cut) then          ! r2cut  is square of rcut  
       r2=1/r2 
       fac2=r2*Sig*Sig 
       fac6=fac2*fac2*fac2 
       df=48.d0*Eps*r2*fac6*(fac6-0.5d0)
       fc(1)=df*dr(1)
       fc(2)=df*dr(2)
       fc(3)=df*dr(3)
       Force(atom1,:)=Force(atom1,:)+fc(:) 
       Force(atom2,:)=Force(atom2,:)-fc(:) 
       PE=PE+4.d0*Eps*fac6*(fac6-1)-Ecut        !shifted to zero at cutoff 
     endif 
   enddo
 enddo 

!print

! writing forces 
! open(unit=200,file='md.force',action='write') 
!   write(200,*)  
! do atom1=1,TotAtom
!   write(200,"(3F15.8)") Force(atom1,:)
! enddo

end subroutine force_calc 
