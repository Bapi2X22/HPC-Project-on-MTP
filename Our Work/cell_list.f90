subroutine cell_list(TotAtom,Rcut,Rverlet,Box,r)
  implicit none
 
  integer :: ierr, atomi,xdim,i,a,b,c,nbr,nnbrs,j,cell(3),nbri,cx,cy,cz
  integer,intent(in) :: TotAtom
  real(kind=8) :: rn,dr(3)
  real(kind=8), intent(in) :: Rcut,Rverlet,Box
 
  integer, allocatable :: hoc(:,:,:),ll(:),nbrs(:)
  real(kind=8),allocatable, intent(in) :: r(:,:)
  
  integer, allocatable :: vnbrs(:)
  
  allocate(nbrs(TotAtom))
  allocate(ll(TotAtom))

  rn = Box/int(Box/Rcut)
  xdim = Box/rn
 
  allocate(hoc(xdim,xdim,xdim))
  
  hoc = 0
  ll = 0
 
  do atomi = 1,TotAtom
  
   cx = floor(r(atomi,1)/rn) + xdim/2
   cy = floor(r(atomi,2)/rn) + xdim/2
   cz = floor(r(atomi,3)/rn) + xdim/2
   
   ll(atomi) = hoc(cx,cy,cz)
   hoc(cx,cy,cz)=atomi    
  
  end do

  open(unit=2,file='cell_list.xyz',action='write',iostat=ierr)
  open(unit=3,file='num_neighbours.xyz',action='write',iostat=ierr)
    
  do atomi=1,TotAtom
  
   cx = floor(r(atomi,1)/rn) + xdim/2
   cy = floor(r(atomi,2)/rn) + xdim/2
   cz = floor(r(atomi,3)/rn) + xdim/2
   
   nbrs = 0
   j = 1  
   do a=-1,1
    do b=-1,1
     do c=-1,1
     
        cell(1) = cx+a   
        cell(2) = cy+b   
        cell(3) = cz+c   
        
        do i=1,3
         if(cell(i)<0) then
          cell(i) = cell(i) + xdim
         else if (cell(i)>xdim) then
          cell(i) = cell(i) - xdim
         end if
        end do
         
   	nbr = hoc(cell(1),cell(2),cell(3))
   	
   	if(nbr>0) then
     	 nbrs(j) = nbr
     	 j=j+1
   	endif
   
   	do while(nbr>0)
    	 nbr = ll(nbrs(j-1))
    
    	 if(nbr>0) then
      	  nbrs(j) = nbr
      	  j=j+1
    	 end if
   	end do 
   	
     end do
    end do
   end do
   
   write(3,*) j-1
   write(2,*) nbrs(1:j-1)
   
  end do
    
 end subroutine
