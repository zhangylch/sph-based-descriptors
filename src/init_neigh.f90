module initmod
     implicit none
     integer(kind=4),parameter :: intype=4,typenum=8
     integer(kind=intype) :: interaction,max_neigh,numatom,length
     integer(kind=intype) :: nimage(3),rangebox(3)
     integer(kind=intype),allocatable :: index_rs(:,:,:),index_numrs(:,:,:,:,:)
     real(kind=typenum) :: rc,rcsq,dier
     real(kind=typenum) :: matrix(3,3),inv_matrix(3,3),rangecoor(3)
     real(kind=typenum),allocatable :: imageatom(:,:,:),shiftvalue(:,:)
end module

subroutine init_neigh(in_rc,in_dier,in_max_neigh,in_numatom,cell)
     use initmod
     implicit none
     integer(kind=intype) :: in_max_neigh,in_numatom
     real(kind=typenum) :: in_rc,in_dier
     real(kind=typenum) :: s1,s2,rlen1,rlen2
     real(kind=typenum) :: cell(3,3)
     real(kind=typenum) :: tmp(3),maxd(3),mind(3),vel1(3,3),vel2(3,3)
!f2py real(kind=typenum),intent(in) :: in_rc,in_dier,max_neigh
!f2py real(kind=typenum),intent(in) :: cell
       max_neigh=in_max_neigh
       numatom=in_numatom
       rc=in_rc
       dier=in_dier
       rcsq=rc*rc
       interaction=ceiling(rc/dier)
       matrix=cell
!Note that the fortran store the array with the column first, so the lattice parameters is the transpose of the its realistic shape
       s1=abs(matrix(1,1)*matrix(2,2))
       s2=abs(matrix(1,1)*matrix(3,3))
       rlen1=sqrt(dot_product(matrix(:,1),matrix(:,1)))
       rlen2=sqrt(dot_product(matrix(:,2),matrix(:,2)))
       tmp(1)=s1/rlen2
       tmp(2)=s1/rlen1
       tmp(3)=s2/rlen1
       nimage=ceiling(rc/abs(tmp))
       length=(2*nimage(1)+1)*(2*nimage(2)+1)*(2*nimage(3)+1)
       call inverse_matrix(matrix,inv_matrix)
       do i = 1,3
         maxd(i)=maxval(matrix(i,:))
         mind(i)=minval(matrix(i,:))
       end do
       rangecoor=maxd-mind+2.0*rc
       rangebox=ceiling(rangecoor/dier)

! allocate the array
       allocate(imageatom(3,numatom,length)) 
       allocate(index_rs(rangebox(1),rangebox(2),rangebox(3)))
       allocate(index_numrs(2,max_neigh,rangebox(1),rangebox(2),rangebox(3)))
       allocate(shiftvalue(3,length))
! obatin image 
       vec2(:,1)=-matrix(:,1)*nimage(1)
       vec2(:,2)=-matrix(:,2)*nimage(2)
       vec2(:,3)=-matrix(:,3)*nimage(3)
       l=0
       vec1(:,3)=vec2(:,3)
       do i=-nimage(3),nimage(3)
         vec1(:,2)=vec2(:,2)
         do j=-nimage(2),nimage(2)
           vec1(:,1)=vec2(:,1)
           do k=-nimage(1),nimage(1)
             l=l+1
             shiftvalue(:,l)=vec1(:,1)+vec1(:,2)+vec1(:,3)
             vec1(:,1)=vec1(:,1)+matrix(:,1)
           end do
           vec1(:,2)=vec1(:,2)+matrix(:,2)
         end do
         vec1(:,3)=vec1(:,3)+matrix(:,3)
       end do
       return
end subroutine

subroutine deallocate_all()
     use initmod
       deallocate(imageatom)
       deallocate(index_numrs)
       deallocate(index_rs)
       deallocate(shiftvalue)
end subroutine
