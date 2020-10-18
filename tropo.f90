subroutine tropo(temp, nlon, nlat, nlev, pres, plimu, pliml, plimlex, dofill, tp, tperr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! determination of tropopause height from gridded temperature data
!
! reference:  Reichler, T., M. Dameris, R. Sausen (2003): 
!             Determining the tropopause height from gridded data, 
!             Geophys. Res. L., 30, No. 20, 2042
!
! input:    temp(nlon,nlat,nlev)    3D-temperature field
!           nlon                    # of grid points in x
!           nlat                    # of grid points in y
!           nlev                    # of vertical pressure levels
!           pres(nlev)              array of pressure levels in hPa, length = nlev
!           plimu                   upper limit for tropopause pressure in Pa, usually 45000
!           pliml                   lower limit for tropopause pressure in Pa, usually 7500
!           plimlex                 lower limit in extratropics, usually same as pliml, i.e., 7500
!           dofill                  fill undefined values with neighboring points if .true.
!
! output:   tp(nlon, nlat)          tropopause pressure in Pa, same horizontal dimension as temp
!           tperr                   # of undetermined values
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none

integer,intent(in)                        :: nlon, nlat, nlev
real,intent(in),dimension(nlon,nlat,nlev) :: temp
real,intent(in),dimension(nlev)           :: pres
real, intent(in)                          :: plimu, pliml, plimlex
logical, intent(in)                       :: dofill
real,intent(out),dimension(nlon,nlat)     :: tp
integer,intent(out)                       :: tperr

integer                                   :: i, j, invert, ifil
integer                                   :: lon, lat
real,dimension(nlev)                      :: t
real,dimension(nlev)                      :: p
real                                      :: trp

real, parameter                           :: gamma=-0.002 ! K/m

! check vertical orientation of data
if (pres(1) .gt. pres(2)) then
  invert=1
  do i=1,nlev
    p(i)=pres(nlev+1-i)*100.  ! hPa > Pa
  enddo
else
  invert=0
  do i=1,nlev
    p(i)=pres(i)*100.         ! hPa > Pa
  enddo
endif

tperr = 0
do lon=1,nlon
do lat=1,nlat
  if (invert.eq.1) then
    do i=1,nlev
      t(i)=temp(lon,lat,nlev+1-i)
    enddo
  else
    do i=1,nlev
      t(i)=temp(lon,lat,i)
    enddo
  endif
  call twmo(nlev, t, p, plimu, pliml, gamma, trp)
  if (lat.lt..15*nlat.and.trp.lt.plimlex) trp=-998.
  if (lat.gt..85*nlat.and.trp.lt.plimlex) trp=-997.
  tp(lon,lat)=trp
  if(trp.lt..0) then
    tperr = tperr+1    
  endif
end do
end do

! fill holes
if (dofill) then
  call fill(tp, nlon, nlat, ifil)
  if (ifil.ne.tperr) then
    print*, 'Inconsistent'
    stop
  endif
endif

return
end subroutine tropo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! twmo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine twmo(level, t, p, plimu, pliml, gamma, trp)

implicit none
integer,intent(in)                  :: level
real,intent(in),dimension(level)    :: t, p
real,intent(in)                     :: plimu, pliml, gamma
real,intent(out)                    :: trp

real,parameter                      :: kap=0.286
real,parameter                      :: faktor = -9.81/287.0
real,parameter                      :: deltaz = 2000.0
real,parameter                      :: ka1=kap-1.

real                                :: pmk, pm, a, b, tm, dtdp, dtdz
real                                :: ag, bg, ptph
real                                :: pm0, pmk0, dtdz0
real                                :: p2km, asum, aquer
real                                :: pmk2, pm2, a2, b2, tm2, dtdp2, dtdz2
integer                             :: icount, jj
integer                             :: j

trp=-999.0                           ! negative means not valid
do j=level,2,-1

   ! dt/dz
   pmk= .5 * (p(j-1)**kap+p(j)**kap)
   pm = pmk**(1/kap)              
   a = (t(j-1)-t(j))/(p(j-1)**kap-p(j)**kap)
   b = t(j)-(a*p(j)**kap)
   tm = a * pmk + b              
   dtdp = a * kap * (pm**ka1)
   dtdz = faktor*dtdp*pm/tm

   ! dt/dz valid?
   if (j.eq.level)    go to 999     ! no, start level, initialize first
   if (dtdz.le.gamma) go to 999     ! no, dt/dz < -2 K/km
   if (pm.gt.plimu)   go to 999     ! no, too low

   ! dtdz is valid, calculate tropopause pressure
   if (dtdz0.lt.gamma) then
      ag = (dtdz-dtdz0) / (pmk-pmk0)    
      bg = dtdz0 - (ag * pmk0)         
      ptph = exp(log((gamma-bg)/ag)/kap)
   else
      ptph = pm
   endif

   if (ptph.lt.pliml) go to 999    
   if (ptph.gt.plimu) go to 999          

   ! 2nd test: dtdz above 2 km must not exceed gamma
   p2km = ptph + deltaz*(pm/tm)*faktor          ! p at ptph + 2km
   asum = 0.0                                   ! dtdz above
   icount = 0                                   ! number of levels above

   ! test until apm < p2km
   do jj=j,2,-1

       pmk2 = .5 * (p(jj-1)**kap+p(jj)**kap)    ! p mean ^kappa
       pm2 = pmk2**(1/kap)                      ! p mean
       if(pm2.gt.ptph) go to 110                ! doesn't happen
       if(pm2.lt.p2km) go to 888                ! ptropo is valid

       a2 = (t(jj-1)-t(jj))                     ! a
       a2 = a2/(p(jj-1)**kap-p(jj)**kap)
       b2 = t(jj)-(a2*p(jj)**kap)               ! b
       tm2 = a2 * pmk2 + b2                     ! T mean
       dtdp2 = a2 * kap * (pm2**(kap-1))        ! dt/dp
       dtdz2 = faktor*dtdp2*pm2/tm2
       asum = asum+dtdz2
       icount = icount+1
       aquer = asum/float(icount)               ! dt/dz mean
  
       ! discard ptropo ?
        if (aquer.le.gamma) go to 999           ! dt/dz above < gamma

110 continue
    enddo                           ! test next level

888 continue                        ! ptph is valid
    trp = ptph
    return

999 continue                        ! continue search at next higher level
    pm0 = pm
    pmk0 = pmk
    dtdz0  = dtdz

enddo

! no tropopouse found
return
end subroutine twmo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fill(dat, ix, iy, ir)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, intent(in)                 :: ix, iy
integer, intent(out)                :: ir
real, dimension(ix,iy)              :: dat
real, dimension(4)                  :: help
integer                             :: jx, jy, icount, ipk, ic
real                                :: drop, sum

icount = 0
do jx=1,ix
do jy=1,iy
  if (loch(dat(jx,jy))) icount = icount+1
enddo
enddo
if (icount.gt.(ix*iy)/2) stop 'ERROR: Too many holes (>50%)'
ir = icount
if (icount.eq.0) return

ipk = 0
10   continue
do jx=1,ix
do jy=1,iy
     if(loch(dat(jx,jy))) then
          drop = dat(jx,jy)

        ! left edge
          if (jx.eq.1) then
          if (jy.eq.1) then
          help(1) = dat(jx,jy+1)
          help(2) = dat(jx+1,jy)
          help(3) = drop
          help(4) = drop
          go to 200
          endif
          if (jy.eq.iy) then
          help(1) = drop
          help(2) = dat(jx+1,jy)
          help(3) = dat(jx,jy-1)
          help(4) = drop
          go to 200
          endif
          help(1) = dat(jx,jy+1)
          help(2) = dat(jx+1,jy)
          help(3) = dat(jx,jy-1)
          help(4) = drop
          go to 200
          endif

          ! right edge
          if (jx.eq.ix) then
          if (jy.eq.1) then
          help(1) = dat(jx,jy+1)
          help(2) = drop
          help(3) = drop
          help(4) = dat(jx-1,jy)
          go to 200
          endif
          if (jy.eq.iy) then
          help(1) = drop
          help(2) = drop
          help(3) = dat(jx,jy-1)
          help(4) = dat(jx-1,jy)
          go to 200
          endif
          help(1) = dat(jx,jy+1)
          help(2) = drop
          help(3) = dat(jx,jy-1)
          help(4) = dat(jx-1,jy)
          go to 200
          endif

        ! bottom edge
          if (jy.eq.1) then
          help(1) = dat(jx,jy+1)
          help(2) = dat(jx+1,jy)
          help(3) = drop
          help(4) = dat(jx-1,jy)
          go to 200
          endif

          ! upper edge
          if(jy.eq.iy) then
          help(1) = drop
          help(2) = dat(jx+1,jy)
          help(3) = dat(jx,jy-1)
          help(4) = dat(jx-1,jy)
          go to 200
          endif

        ! no edge
          help(1) = dat(jx,jy+1)
          help(2) = dat(jx+1,jy)
          help(3) = dat(jx,jy-1)
          help(4) = dat(jx-1,jy)

 200      continue

          ic = 0
          sum = 0.0
          do jj=1,4
            if(.not.loch(help(jj))) then
              sum = sum+help(jj)
              ic = ic+1
            endif
          enddo

          if (ic.gt.0) then
            dat(jx,jy) = sum/float(ic)  ! fill with mean of valid
            ipk = ipk+1                 ! neighbourpoints
          endif

     endif
     if (ipk .ge. icount) return    ! until all filled
enddo
enddo
go to 10

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical function loch(x)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real, intent(in)        :: x

edge = -98.0
if (x.lt.edge) then
  loch = .true.
else
  loch = .false.
endif
return
end function loch

end
