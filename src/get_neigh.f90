subroutine EANN_out(coor,atomindex,shift_image)
     use initmod
     implicit none
     integer(kind=intype) :: numatom,num,maxneigh
     integer(kind=intype) :: i,i1,i2,i3,j,k,l,iatom
     integer(kind=intype) :: sca(3),rangebox(3),boundary(2,3)
     real(kind=typenum) :: s1,s2,rlen1,rlen2
     real(kind=typenum) :: cutoff
     real(kind=typenum) :: oriminv(3),rangecoor(3),tmp(3),tmp1(3)
     real(kind=typenum) :: coor(3,numatom),cart(3,numatom),fcoor(3,numatom)
     real(kind=typenum) :: matrix(3,3),inv_matrix(3,3)
     real(kind=typenum) :: vec1(3,3),vec2(3,3)
!f2py integer(kind=intype),intent(in) :: numatom
!f2py real(kind=typenum),intent(in),check(0==0) :: coor,matrix
       fcoor=matmul(inv_matrix,coor)
! move all atoms to an cell which is convenient for the expansion of the image
       oriminv=cart(:,1)
       do i=2,numatom
         sca=nint(fcoor(:,i)-fcoor(:,1))
         coor(:,i)=coor(:,i)-sca(1)*scalmatrix(:,1)-sca(2)*scalmatrix(:,2)-sca(3)*scalmatrix(:,3)
         do j=1,atomdim
           if(coor(j,i)<oriminv(j)) oriminv(j)=coor(j,i)
         end do
       end do
       oriminv=oriminv-maxrc
       do i=1,numatom
         cart(:,i)=cart(:,i)-oriminv
       end do
! obatin image 
       l=0
       do i=-nimage(3),nimage(3)
         do j=-nimage(2),nimage(2)
           do k=-nimage(1),nimage(1)
             l=l+1
             do num=1,numatom
               imageatom(:,num,l)=cart(:,num)+shiftvalue(:,l)
             end do
           end do
         end do
       end do
       index_rs=0
       do num=1,length
         do i=1,numatom
           dire(:,i,num)=dire(:,i,num)-oriminv
           if(dire(1,i,num)>0d0.and.dire(1,i,num)<rangecoor(1).and.dire(2,i,num)>0d0.and.dire(2,i,num)<rangecoor(2)  &
           .and.dire(3,i,num)>0d0.and.dire(3,i,num)<rangecoor(3)) then
             sca=ceiling(dire(:,i,num)/dier)
             index_rs(sca(1),sca(2),sca(3))=index_rs(sca(1),sca(2),sca(3))+1
             index_numrs(:,index_rs(sca(1),sca(2),sca(3)),sca(1),sca(2),sca(3))=[i-1,num] 
           end if
         end do
       end do
       scutnum=0
       do iatom = 1, numatom
         sca=ceiling(cart/dier)
         ninit=(length+1)/2
         dire(:,natom,ninit)=100d0
         do i=1,3
           boundary(1,i)=max(1,sca(i)-interaction)
           boundary(2,i)=min(numrs(i),sca(i)+interaction)
         end do
         do i3=boundary(1,3),boundary(2,3)
           do i2=boundary(1,2),boundary(2,2)
             do i1=boundary(1,1),boundary(2,1)
               do i=1,index_rs(i1,i2,i3)
                 j=index_numrs(1,i,i1,i2,i3)
                 num=index_numrs(2,i,i1,i2,i3)
                 tmp1=dire(:,j,num)-cart
                 tmp=dot_product(tmp1,tmp1)
                 if(tmp<=rcsq(ntype)) then
                   atomindex(:,scutnum)=[iatom-1,j]
                   shift(:,scutnum)=[shiftvalue(1,num),shiftvalue(2,num),shiftvalue(3,num)]
                   scutnum=scutnum+1
                 end if
               end do
             end do
           end do
         end do
         dire(:,natom,ninit)=cart
       end do
     return
end subroutine
   
