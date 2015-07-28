        program procmss
	
	implicit none
	character, allocatable :: band(:,:)
	integer(2), allocatable :: sband(:,:,:)
	integer(2), allocatable :: oband(:,:)
	character(80) filename(4)
	character(4) suffix(4)
	integer*8 padding
	data suffix /"_B40","_B50","_B60","_B70"/
	integer ii,nr,nc,ib,als,i,j,ierr,nt
	integer igain(4)
	real lmin(4),lmax(4),qmin(4),qmax(4)
	real esun(4)
	real gain(4),offset(4)
	real xa(4),xb(4),xc(4)
	real xx
c ancillary data
	integer(2), allocatable :: ozone(:,:)        	
	integer(2), allocatable :: dem(:,:)        	
	integer(2), allocatable :: slp(:,:,:)        	
	integer(2), allocatable :: wv(:,:,:)        	
	integer(2), allocatable :: islp(:,:)        	
	integer(2), allocatable :: iwv(:,:)        	
	character(80) fileozone,filereana,filedem
	real twv(4),wvint,tslp(4),slpint
	integer it
	integer nco,nro,ncr,nrr,ntr,ncd,nrd
	
c	data llmina /-6.2, -6.0, -4.5, -4.5, -1.0,-0.35, 0.0/
c	data hlmina / -6.2, -6.0, -4.5, -4.5, -1.0,-0.35,3.2/
c        data llmaxa/ 297.5, 303.4, 235.5, 235.0, 47.70, 16.60,17.04/
c	data hlmaxa / 194.3, 202.4, 158.6, 157.5, 31.76,10.932,12.65/
	
c	data llminb / -6.2, -6.4, -5.0, -5.1, -1.0, -0.35,0.0/
c	data hlminb / -6.2, -6.4, -5.0, -5.1, -1.0, -0.35,3.2/
c        data llmaxb/ 293.7, 300.9, 234.4, 241.1, 47.57,16.54, 17.04 /
c	data hlmaxb /191.6, 196.5, 152.9, 157.4, 31.06, 10.80, 12.65/
	
	data esun /1805.42,1511.42,1241.4,908.018/
	real ah2o(4),bh2o(4)
	data ah2o /0.00251053,0.00768567,0.0249089,0.105129/
	data bh2o /0.753595,0.687902,0.609682,0.426469/
	real aoz(4)
	data aoz /-0.084081,-0.0582646,-0.00965367,0.0/
	real a1(4),b0(4),b1(4)
	data a1 /9.82973e-06,0.00716206,0.0297642,8.56296e-06/
	data b0 /0.639237,0.614214,0.690885,0.140782/
	data b1 /-0.240851,0.193631,0.564946,-0.110042/
	real tauray(4)
	data tauray /0.09817,0.04713,0.02862,0.01428/
	real m,xtv,uh2o,uoz,tgwv,tgoz,tgog,pres,xtaur,roray
	real xmus,xmuv,xtts,xttv,ttatm,satm,xphi,rotoa
	real xts,xfs,irad(4),dsol,cpi,pi,om
	real lat11,lon11,lat12,lon12,lat21,lon21,lat22,lon22
	real xlat,xlon,mlat,mlon,x,y
	integer row,col,alti,salti
	real time
	integer doy
	
! BEGIN OF HDF PARAMETER BLOCK        
	 character*4 sdsname
         integer sfstart, sfselect, sfrdata, sfendacc, sfend
	 integer sfscatt,sfwdata,sfcreate,sfginfo
         integer sd_id, sd_id2,sds_id, sds_index, status
         integer start(5), edges(5), stride(5)
         integer nstart(2), nedges(2), nstride(2)
         integer DFACC_READ,DFACC_RDWR,DFNT_CHAR8,DFACC_CREATE
	 integer DFNT_INT16
         parameter (DFACC_READ = 1)
         parameter (DFACC_RDWR = 3)
         parameter (DFACC_CREATE = 4)
         parameter (DFNT_INT16 = 22)
         parameter (DFNT_CHAR8 = 4)
	 character*80 sds_name
	 integer rank,data_type
	 integer n_attrs
	 integer dim_sizes(5)
	 integer dims(2)
         integer dim_length(5), comp_type, comp_prm(4)
! END OF HDF PARAMETER BLOCK	 
	
	cpi=atan(1.)*4.
	pi=cpi
	do i=1,4
	read(5,'(A80)') filename(i)
	enddo
	read(5,*) (lmax(i),i=1,4)
	read(5,*) (lmin(i),i=1,4)
	read(5,*) (qmax(i),i=1,4)
	read(5,*) (qmin(i),i=1,4)
	read(5,*) nr,nc
	read(5,*) xts,xfs
	read(5,*) time
	read(5,*) lat11,lon11
	read(5,*) lat12,lon12
	read(5,*) lat21,lon21
	read(5,*) lat22,lon22
        read(5,*) doy
        om=(.9856*float(doy-4))*pi/180.
        dsol=1./((1.-.01673*cos(om))**2)
        read(5,'(A80)') fileozone
	write(6,*) fileozone
        read(5,'(A80)') filereana
	write(6,*) filereana
        read(5,'(A80)') filedem
 	write(6,*) filedem
	
c reading ancillary data (ozone)	
	ii=index(fileozone," ")-1
        sd_id= sfstart(fileozone(1:ii),DFACC_RDWR)
        write(6,*) "sd_id",sd_id
        sds_index = 2
        sds_id    = sfselect(sd_id, sds_index)
        write(6,*) "sds_id", sds_id
        status= sfginfo(sds_id, sds_name, rank, dim_sizes, data_type,n_attrs)
        write(6,*) "sdsname ",sds_name
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       write(6,*) dim_sizes(1),dim_sizes(2)
       nco= dim_sizes(1)
       nro=dim_sizes(2)
       allocate (ozone(nco,nro),stat=ierr)
cc read  data
       start(1)=0
       start(2) = 0
       edges(1) = nco
       edges(2) = nro
       stride(1) = 1
       stride(2) = 1
       sds_index = 2
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,ozone)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
c close HDF file
       status = sfend(sd_id)
       write(6,*) "status sfend ",status
	
c reading ancillary data (reanalysis)	
	ii=index(filereana," ")-1
        sd_id= sfstart(filereana(1:ii),DFACC_RDWR)
        write(6,*) "sd_id",sd_id
        sds_index = 2
        sds_id    = sfselect(sd_id, sds_index)
        write(6,*) "sds_id", sds_id
        status= sfginfo(sds_id, sds_name, rank, dim_sizes, data_type,n_attrs)
        write(6,*) "sdsname ",sds_name
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       write(6,*) dim_sizes(1),dim_sizes(2),dim_sizes(3)
       ncr= dim_sizes(1)
       nrr=dim_sizes(2)
       ntr=dim_sizes(3)
       allocate (slp(ncr,nrr,ntr),stat=ierr)
       allocate (wv(ncr,nrr,ntr),stat=ierr)
cc read  slp and wv
       start(1)=0
       start(2) = 0
       start(3) = 0
       edges(1) = ncr
       edges(2) = nrr
       edges(3) = ntr
       stride(1) = 1
       stride(2) = 1
       stride(3) = 1
       sds_index = 2
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,slp)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       sds_index = 3
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,wv)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
c close HDF file
       status = sfend(sd_id)
       write(6,*) "status sfend ",status
       
       
c reading ancillary data (dem)	
	ii=index(filedem," ")-1
        sd_id= sfstart(filedem(1:ii),DFACC_RDWR)
        write(6,*) "sd_id",sd_id
        sds_index = 0
        sds_id    = sfselect(sd_id, sds_index)
        write(6,*) "sds_id", sds_id
        status= sfginfo(sds_id, sds_name, rank, dim_sizes, data_type,n_attrs)
        write(6,*) "sdsname ",sds_name
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       write(6,*) dim_sizes(1),dim_sizes(2)
       ncd= dim_sizes(1)
       nrd=dim_sizes(2)
       allocate (dem(ncd,nrd),stat=ierr)
cc read  data
       start(1)=0
       start(2) = 0
       edges(1) = ncd
       edges(2) = nrd
       stride(1) = 1
       stride(2) = 1
       sds_index = 0
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,dem)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
c close HDF file
       status = sfend(sd_id)
       write(6,*) "status sfend ",status
       
       
c	stop 
	
	
	do i=1,4
	irad(i)=esun(i)*cos(xts*cpi/180.)*dsol/cpi

	gain(i)=(lmax(i)-lmin(i))/(qmax(i)-qmin(i))
	offset(i)=lmin(i)
	enddo
	
        allocate(band(nr,nc),STAT=als)
        allocate(sband(7,nr,nc),STAT=als)
        allocate(oband(nc,nr),STAT=als)
	
	xtv=0.
	m=1./cos(xts*cpi/180.)+1/cos(xtv*cpi/180.)
c compute mean ozone, water vapor and pressure (slp) for the whole scene
        mlat=(lat11+lat12+lat21+lat22)/4.	
        mlon=(lon11+lon12+lon21+lon22)/4.
c ozone row,col
        row=(90.0-mlat)/(180./180)+1
	col=(180+mlon)/(360./288)+1
	if (row.le.1) row=1		
	if (col.le.1) col=1		
	if (row.ge.180) row=180		
	if (col.ge.288) col=288
	uoz=ozone(col,row)/1000.
	Write(6,*) "row, col , uoz ",row,col,uoz
c water vapor	
        row=(90.0-mlat)/(180./73)+1
	col=(180+mlon)/(360./144)+1
	if (row.le.1) row=1		
	if (col.le.1) col=1		
	if (row.ge.180) row=73		
	if (col.ge.288) col=144
	twv(1)=(wv(col,row,1)*0.0099999998+277.64999)/10.
	twv(2)=(wv(col,row,2)*0.0099999998+277.64999)/10.
	twv(3)=(wv(col,row,3)*0.0099999998+277.64999)/10.
	twv(4)=(wv(col,row,4)*0.0099999998+277.64999)/10.
c sea level pressure
	tslp(1)=(slp(col,row,1)*1.0+119765.)/100.
	tslp(2)=(slp(col,row,2)*1.0+119765.)/100.
	tslp(3)=(slp(col,row,3)*1.0+119765.)/100.
	tslp(4)=(slp(col,row,4)*1.0+119765.)/100.

	it=int(time/6)+1
	if (it.eq.3) it=3
	wvint=(time-(it-1)*6)*(twv(it+1)-twv(it))/6.+twv(it)
	slpint=(time-(it-1)*6)*(tslp(it+1)-tslp(it))/6.+tslp(it)
	write(6,*) "time row col wv ",time,row,col,twv(1),twv(2),twv(3),twv(4),wvint
	write(6,*) "time row col slp ",time,row,col,tslp(1),tslp(2),tslp(3),tslp(4),slpint

			
	uh2o=3.0
c	uoz=0.35
	pres=slpint/1013.
	xphi=0.
	do ib=1,4
	tgwv=exp(-ah2o(ib)*((m*uh2o)**bh2o(ib)))
	tgoz=exp(aoz(ib)*m*uoz)
	tgog=exp(-(a1(ib)*pres)*(m**(exp(-(b0(ib)+b1(ib)*pres)))))
        xtaur=tauray(ib)*pres
        xmus=cos(xts*cpi/180.)
        xmuv=cos(xtv*cpi/180.)
        call comptransray(xtaur,xmus,xtts)
        call comptransray(xtaur,xmuv,xttv)
C Compute total transmission (product downward by  upward)
       ttatm=xtts*xttv
       call local_csalbr(xtaur,satm)
       call local_chand(xphi,xmus,xmuv,xtaur,roray)
	
	ii=index(filename(ib)," ")-1
	write(6,*) "reading ",filename(ib)(1:ii)
	  open(1,file=filename(ib)(1:ii),form='UNFORMATTED',action='READ',access='DIRECT',recl=nc*nr+8)
	  read(1,rec=1) padding,((band(i,j),j=1,nc),i=1,nr)
	  close(1)
	  salti=-9999
	     do i=1,nr
	     do j=1,nc
C compute lat lon
             y=(i-1.)/(nr-1.)
	     x=(j-1.)/(nc-1.)
	     xlat=lat11*(1-x)*(1-y)+lat12*x*(1-y)+lat21*y*(1-x)+lat22*x*y
	     xlon=lon11*(1-x)*(1-y)+lon12*x*(1-y)+lon21*y*(1-x)+lon22*x*y
c compute pressure based on 1013. sea level presure 
             row=(90.0-xlat)/0.05+1
	     col=(180.0+xlon)/0.05+1
	     if (row.le.1) row=1	   
	     if (col.le.1) col=1	   
	     if (row.ge.3600) row=3600	   
	     if (col.ge.7200) col=7200	   
	     alti=dem(row,col)
	     if (alti.eq.-9999) alti=0
C update rayleigh with pressure if necessary 	     
	     if (alti.ne.salti) then
	     salti=alti
	     pres=slpint*exp(-salti/8500.)/1013.
             call comptransray(xtaur,xmus,xtts)
	     call comptransray(xtaur,xmuv,xttv)
C Compute total transmission (product downward by  upward)
             ttatm=xtts*xttv
             call local_csalbr(xtaur,satm)
             call local_chand(xphi,xmus,xmuv,xtaur,roray)
	     endif
	     rotoa=(ichar(band(i,j))*gain(ib)+offset(ib))/irad(ib)
	     xx=rotoa/(tgog*tgoz)
	     xx=(xx-roray)/(ttatm*tgwv)
	     xx=xx/(1.+satm*xx)
	     sband(ib,i,j)=int(xx*10000.)
         if ((ib .eq. 1) .and. (i .eq. 1004) .and. (j .eq. 2139)) then
           write(6,*) 'test###', sband(ib, i, j), tgog, tgoz, roray, ttatm, tgwv
         end if
	     enddo
	     enddo
	enddo 
 

!Saving in an HDF file
         sd_id= sfstart('correcteddata.hdf',DFACC_CREATE)
	 do ib=1,4
	 do i=1,nr
	 do j=1,nc
	 oband(j,i)=sband(ib,i,j)
	 enddo
	 enddo
         write(sds_name,'(A4,A4)') "SREF_",suffix(ib)
	 write(6,*) "writing band ",sds_name
	 dim_length(1)=nc
	 dim_length(2)=nr
	 comp_type=4
	 comp_prm(1)=8
	 rank=2
	 dim_sizes(1)=nc
	 dim_sizes(2)=nr
         start(1) = 0
         start(2) = 0
          edges(1) = nc
         edges(2) = nr
         stride(1) = 1
         stride(2) = 1
	 sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,oband)
	 write(6,*) "status sfwdata ",status
         status = sfendacc(sds_id)
	 write(6,*) "status sfendacc ",status
	 enddo
	 do i=1,nr
	 do j=1,nc
	 oband(j,i)=sband(2,i,j)/2.
	 enddo
	 enddo
        write(sds_name,'(A4,A4)') "SREF_","_B10"
	 write(6,*) "writing band ",sds_name
	 dim_length(1)=nc
	 dim_length(2)=nr
	 comp_type=4
	 comp_prm(1)=8
	 rank=2
	 dim_sizes(1)=nc
	 dim_sizes(2)=nr
         start(1) = 0
         start(2) = 0
          edges(1) = nc
         edges(2) = nr
         stride(1) = 1
         stride(2) = 1
	 sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,oband)
	 write(6,*) "status sfwdata ",status
         status = sfendacc(sds_id)
	 write(6,*) "status sfendacc ",status
	 
         status = sfend(sd_id) 	 
	 write(6,*) "status sfend ",status

	stop
	end
	
       subroutine comptransray(xtaur,xmus,ttray)
       real xtaur,xmus,ttray,ddiftt,ddirtt
       
       ddiftt=(2./3.+xmus)+(2./3.-xmus)*exp(-xtaur/xmus)
       ddiftt=ddiftt/((4./3.)+xtaur)-exp(-xtaur/xmus)
       ddirtt=exp(-xtaur/xmus)
       ttray=ddirtt+ddiftt
       return
       end
      subroutine local_csalbr(xtau,xalb)
      real xtau,xalb,fintexp3
      xalb=(3*xtau-fintexp3(xtau)*(4+2*xtau)+2*exp(-xtau))
      xalb=xalb/(4.+3*xtau)
      return
      end
      real function fintexp3(xtau)
      real xx,xtau,fintexp1
      xx=(exp(-xtau)*(1.-xtau)+xtau*xtau*fintexp1(xtau))/2.
      fintexp3=xx
      return
      end
      real function fintexp1(xtau)
c accuracy 2e-07... for 0<xtau<1
      real xx,a(0:5),xtau,xftau
      integer i
      data (a(i),i=0,5) /-.57721566,0.99999193,-0.24991055,
     c                  0.05519968,-0.00976004,0.00107857/
      xx=a(0)
      xftau=1.
      do i=1,5
      xftau=xftau*xtau
      xx=xx+a(i)*xftau
      enddo
      fintexp1=xx-log(xtau)
      return
      end

	subroutine local_chand (xphi,xmuv,xmus,xtau
     s			,xrray)
c input parameters: xphi,xmus,xmuv,xtau
c xphi: azimuthal difference between sun and observation (xphi=0,
c in backscattering) and expressed in degree (0.:360.)
c xmus: cosine of the sun zenith angle
c xmuv: cosine of the observation zenith angle
c xtau: molecular optical depth
c output parameter: xrray : molecular reflectance (0.:1.)
c constant : xdep: depolarization factor (0.0279)
        parameter (fac = 0.017453293)
	real xdep,pl(10)
	real fs0,fs1,fs2
	real as0(10),as1(2),as2(2)
        real xphi,xmus,xmuv,xtau,xrray,pi,phios,xcosf1,xcosf2
        real xcosf3,xbeta2,xfd,xph1,xph2,xph3,xitm, xp1, xp2, xp3
        real cfonc1,cfonc2,cfonc3,xlntau,xitot1,xitot2,xitot3
        integer i
	data (as0(i),i=1,10) /.33243832,-6.777104e-02,.16285370
     s	,1.577425e-03,-.30924818,-1.240906e-02,-.10324388
     s	,3.241678e-02,.11493334,-3.503695e-02/
	data (as1(i),i=1,2) /.19666292, -5.439061e-02/
	data (as2(i),i=1,2) /.14545937,-2.910845e-02/
C	pi=3.1415927
C	fac=pi/180.
	phios=180.-xphi
	xcosf1=1.
	xcosf2=cos(phios*fac)
	xcosf3=cos(2*phios*fac)
	xbeta2=0.5
	xdep=0.0279
	xfd=xdep/(2-xdep)
	xfd=(1-xfd)/(1+2*xfd)
	xph1=1+(3*xmus*xmus-1)*(3*xmuv*xmuv-1)*xfd/8.
	xph2=-xmus*xmuv*sqrt(1-xmus*xmus)*sqrt(1-xmuv*xmuv)
	xph2=xph2*xfd*xbeta2*1.5
	xph3=(1-xmus*xmus)*(1-xmuv*xmuv)
	xph3=xph3*xfd*xbeta2*0.375
	xitm=(1-exp(-xtau*(1/xmus+1/xmuv)))*xmus/(4*(xmus+xmuv))
	xp1=xph1*xitm
	xp2=xph2*xitm
	xp3=xph3*xitm
	xitm=(1-exp(-xtau/xmus))*(1-exp(-xtau/xmuv))
	cfonc1=xph1*xitm
	cfonc2=xph2*xitm
	cfonc3=xph3*xitm
	xlntau=log(xtau)
	pl(1)=1.
	pl(2)=xlntau
	pl(3)=xmus+xmuv
	pl(4)=xlntau*pl(3)
	pl(5)=xmus*xmuv
	pl(6)=xlntau*pl(5)
	pl(7)=xmus*xmus+xmuv*xmuv
	pl(8)=xlntau*pl(7)
	pl(9)=xmus*xmus*xmuv*xmuv
	pl(10)=xlntau*pl(9)
	fs0=0.
	do i=1,10
	fs0=fs0+pl(i)*as0(i)
	enddo
	fs1=pl(1)*as1(1)+pl(2)*as1(2)
	fs2=pl(1)*as2(1)+pl(2)*as2(2)
	xitot1=xp1+cfonc1*fs0*xmus
	xitot2=xp2+cfonc2*fs1*xmus
	xitot3=xp3+cfonc3*fs2*xmus
	xrray=xitot1*xcosf1
	xrray=xrray+xitot2*xcosf2*2
	xrray=xrray+xitot3*xcosf3*2
	xrray=xrray/xmus
	return
	end
