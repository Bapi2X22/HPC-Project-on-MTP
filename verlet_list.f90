subroutine verlet_list(TotAtom,Rcut,Rverlet,Box,r)
  implicit none
  
  integer:: ierr,atomi,nbri,k,i
  integer, intent(in) :: TotAtom
  real(kind=8) :: dr(3)
  real(kind=8),intent(in) :: Rcut,Rverlet,Box
  real(kind=8),allocatable,intent(in) :: r(:,:)
  real(kind=8),allocatable :: nnbrs(:)
  integer,allocatable :: neighbours(:,:),vlist(:,:)
  
  allocate(nnbrs(TotAtom))
  allocate(neighbours(TotAtom,TotAtom)) 
  allocate(vlist(TotAtom,TotAtom))  
   
  neighbours=0
  nnbrs=0
  
  open(unit=2,file='cell_list.xyz',action='read',iostat=ierr)
  open(unit=3,file='num_neighbours.xyz',action='read',iostat=ierr)
  
  do atomi=1,TotAtom
   read(3,*) nnbrs(atomi)
  enddo
  
  do atomi=1,TotAtom
   k=nnbrs(atomi)
   read(2,*) neighbours(atomi,1:k) 
  enddo
  
  open(unit=4,file='verlet_list.xyz',action='write',iostat=ierr)
  open(unit=5,file='num_verlet.xyz',action='write',iostat=ierr)
  
  do atomi=1,TotAtom
   i=1
   do k=1,nnbrs(atomi)
    if(atomi/=neighbours(atomi,k)) then
     dr = r(atomi,:)-r(neighbours(atomi,k),:)
      if(dot_product(dr,dr)<Rverlet**2) then
      
       vlist(atomi,i) = neighbours(atomi,k)
       i=i+1
       
      end if
    end if
   end do
   write(4,*) vlist(atomi,1:i-1)
   write(5,*) i-1
  end do  
  
 end subroutine
