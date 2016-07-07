pro congrid_her,irdc
;
;** here I assume:
;   all images are all padded by zeros, not NANs!!
;   there are no NANs anywhere in the images
;   all NANs in images go to a value of zero.
;   all fits files have only one data plane in the zero-th extension
;   the SPIRE 350 um pixel scale = 10"/pix
;   the SPIRE 250 um pixel scale = 6"/pix
;
;you need to include convolve_flag.pro in your IDL path,
;or in the same directory to run this.  
;you need standard IDL astro libraries to run this.
;
;must also include in path: 
;myregrid.pro
;rebinkernel.pro
;ensure_psf_centered.pro
;
;**this code is not public. use at your own risk.
;**AMS 28/09/2010
;
;AMS 03/03/2011 
;updates: 
;* use azymuthally symmetric Gonzalo Aniano kernels
;* SPIRE beam areas
;* extended emission flux correctoins for PACS and SPIRE
;* the convolution procedure and kernel regridding to be more accurate.
;
;SR 20/04/2012
;+ adapted to my setup
;
;SR 19-3-15
; removed sub-mm data -- not needed for Kelly algorithm
;


;@~/.idl/convolution_proc/convolve.pro
;irdc = 'IRDC011.11'
;dir = '/disk1/ragan/herschel/epos/'
dir='/disk1/ragan/launhardt_tmap/'
bdir='/disk1/ragan/epos_balog/PACS_fits/'

    beam = '500' ;to which wavelength beam do you want to convolve, 350 or 500 um?  
    ;beam = '350'

    ;out dir for convolved data, if left blank then = current directory
    path_out = dir+irdc+'/conv500/'

    ;set up a file naming convention here (this is only a suggestion):
    ;Obj name + wavelength + reduction info + convol. wavelength + pixel scale
    out_fn_obj = irdc
    out_fn_info_p = 'conv'+beam 
    out_fn_info_s = 'conv'+beam ; [L2], [BM], [SC]...
    out_fn_info_sm = 'gauss_beam'+beam  

    ;set data and kernel paths
    pathp = dir+irdc+'/'             ; pacs data path
    paths = dir+irdc+'/'             ; spire data path
;    pathsub = dir+irdc+'/'           ; submm data path
    pathk = '/disk1/ragan/kernel/'

    ;set data file names to be read in
    p70_fn  = 'blue_nonan.fits'
    p100_fn = 'green_nonan.fits'
    p160_fn = 'red_nonan.fits'
    s250_fn = 'spire250.fits'
    s350_fn = 'spire350.fits'
    s500_fn = 'spire500.fits'
;    a870_fn = irdc+'_ag_fk5.fits'

    p70i = mrdfits(pathp+p70_fn,0,hp70,/silent)
    sxdelpar,hp70,'XTENSION' 
    p100i = mrdfits(pathp+p100_fn,0,hp100,/silent)
    sxdelpar,hp100,'XTENSION' 
    p160i = mrdfits(pathp+p160_fn,0,hp160,/silent)
    sxdelpar,hp160,'XTENSION'
    s250i = mrdfits(paths+s250_fn,0,hs250,/silent)
    sxdelpar,hs250,'XTENSION'
    s350i = mrdfits(paths+s350_fn ,0,hs350,/silent)
    sxdelpar,hs350,'XTENSION'
    s500i = mrdfits(paths+s500_fn ,0,hs500,/silent)  
    sxdelpar,hs500,'XTENSION'
;    a870i = mrdfits(pathsub+a870_fn,0,ha870,/silent)

;cd1=sxpar(ha870,'CD1_1') & cd2=sxpar(ha870,'CD2_2')
;sxdelpar,ha870,'CD1_1' & sxdelpar,ha870,'CD1_2' & sxdelpar,ha870,'CD2_1' & sxdelpar,ha870,'CD2_2'
;sxaddpar,ha870,'CDELT1',cd1,after='CRPIX1'
;sxaddpar,ha870,'CDELT2',cd2,after='CRPIX2'

    ;read in convolution kernel files [GA's kernels are NOT normalized]
;    if beam eq '500' then begin 
       k70  = mrdfits(pathk+'Kernel_HiRes_PACS_70_to_SPIRE_500.fits.gz',0,hk70i,/silent)
       k100 = mrdfits(pathk+'Kernel_HiRes_PACS_100_to_SPIRE_500.fits.gz',0,hk100i,/silent)
       k160 = mrdfits(pathk+'Kernel_HiRes_PACS_160_to_SPIRE_500.fits.gz',0,hk160i,/silent)
       k250 = mrdfits(pathk+'Kernel_HiRes_SPIRE_250_to_SPIRE_500.fits.gz',0,hk250i,/silent)
       k350 = mrdfits(pathk+'Kernel_HiRes_SPIRE_350_to_SPIRE_500.fits.gz',0,hk350i,/silent)
;       k870 = mrdfits(pathk+'Kernel_HiRes_LABOCA_to_SPIRE_500.fits.gz',0,hk870i,/silent)
;    endif else if beam eq '250' then begin
;       k70  = mrdfits(pathk+'Kernel_HiRes_PACS_70_to_SPIRE_250.fits.gz',0,hk70i,/silent)
;       k100 = mrdfits(pathk+'Kernel_HiRes_PACS_100_to_SPIRE_250.fits.gz',0,hk100i,/silent)
;       k160 = mrdfits(pathk+'Kernel_HiRes_PACS_160_to_SPIRE_250.fits.gz',0,hk160i,/silent)
;    endif

;help,k870

    ;***************************************************************
    ;brute-force NAN-removal, in case some are still in the images...
    nani = where(finite(p70i,/nan) eq 1L,nnan) 
    if nani[0] ne -1L then p70i[nani] = 0.0
    nani = where(finite(p100i,/nan) eq 1L,nnan) 
    if nani[0] ne -1L then p100i[nani] = 0.0
    nani = where(finite(p160i,/nan) eq 1L,nnan) 
    if nani[0] ne -1L then  p160i[nani] = 0.0
    nani = where(finite(s250i,/nan) eq 1L,nnan) 
    if nani[0] ne -1L then s250i[nani] = 0.0
    nani = where(finite(s350i,/nan) eq 1L,nnan) 
    if nani[0] ne -1L then s350i[nani] = 0.0
    nani = where(finite(s500i,/nan) eq 1L,nnan) 
    if nani[0] ne -1L then s500i[nani] = 0.0
;    nani = where(finite(a870i,/nan) eq 1L,nnan) 
;    if nani[0] ne -1L then a870i[nani] = 0.0


    ;extract astrometry from fits files
    extast,hp70,astr70
    extast,hp100,astr100
    extast,hp160,astr160
    extast,hs250,astr250
    extast,hs350,astr350
    extast,hs500,astr500
;    extast,ha870,astr870

    ;calculate pixel scales
    pixs70  = (astr70.cdelt[0]^2)^0.5*3600.   ;["/pix]
    if pixs70 eq float(3600) then $
       pixs70 = (astr70.cd[0]^2. + astr70.cd[1]^2.)^0.5*3600.;["/pix]       
    pixs100 = (astr100.cdelt[0]^2.)^0.5*3600. ;["/pix]
    if pixs100 eq float(3600) then $
       pixs100 = (astr100.cd[0]^2. + astr100.cd[1]^2.)^0.5*3600.;["/pix]
    pixs160 = (astr160.cdelt[0]^2.)^0.5*3600. ;["/pix]
    if pixs160 eq float(3600) then $
       pixs160 = (astr160.cd[0]^2. + astr160.cd[1]^2.)^0.5*3600.;["/pix]
    pixs250 = (astr250.cdelt[0]^2.)^0.5*3600. ;["/pix]
    pixs350 = (astr350.cdelt[0]^2.)^0.5*3600. ;["/pix]
    pixs500 = (astr500.cdelt[0]^2.)^0.5*3600. ;["/pix]
;    pixs870 = abs(astr870.cdelt[0])*3600. ;["/pix]

;    print,'AG image header'
;    print,astr870
;    print,'AG image pixel scale',pixs870

    ;extract astrometry, pixel scales from kernel headers. 
    extast,hk70i, astrk70
    extast,hk100i,astrk100
    extast,hk160i,astrk160
    extast,hk250i,astrk250
;    extast,hk870i,astrk870
    extast,hk350i,astrk350

;   print,'kernel header 870'
;   print,astrk870

    pixsk70 = (astrk70.cd[0]^2. + astrk70.cd[1]^2.)^0.5*3600.;["/pix]
    pixsk100 = (astrk100.cd[0]^2. + astrk100.cd[1]^2.)^0.5*3600.;["/pix]
    pixsk160 = (astrk160.cd[0]^2. + astrk160.cd[1]^2.)^0.5*3600.;["/pix]
    pixsk250 = (astrk250.cd[0]^2. + astrk250.cd[1]^2.)^0.5*3600.;["/pix]
;    pixsk870 = (astrk870.cd[0]^2. + astrk870.cd[1]^2.)^0.5*3600.;["/pix]
    if beam eq '500' then pixsk350 = (astrk350.cd[0]^2. + astrk350.cd[1]^2.)^0.5*3600.;["/pix]
    
;     print,'kernel pixel scale 870',pixsk870

; stop

    ;PACS native units: Jy/pix and the point-source -> extended source corrections
    p70i  = p70i  / (pixs70^2.) 
    p100i = p100i / (pixs100^2.) ;/ 1.151 ;[Jy/arcsec^2]
    p160i = p160i / (pixs160^2.) ;/ 1.174 ;[Jy/arcsec^2]
    ;SPIRE native units: Jy/beam and the point-source -> extended source corrections
    s250i = s250i / 423.0 ;* 0.9828 ;[Jy/arcsec^2]
    s350i = s350i / 751.0 ;* 0.9834 ;[Jy/arcsec^2] [beam area = FWHM^2*pi/(4*ln(2))]
    s500i = s500i / 1587.0 ;* 0.9710 ;[Jy/arcsec^2]
    ;Summ SCUBA & mm IRAM units: Jy/beam 
    ;beam size = 8.6" & 14.9" at 450 & 850
    ;ATLASGAL (870um) = 19.2"
;    a870i = a870i / (!Pi*19.2^2./(4.*ALOG(2)));[Jy/arcsec^2]    

    ;convolve 70 um image   
    ratkimg70 = pixsk70/pixs70
    if ratkimg70 gt 1.0 then x = 1 
    k70r = rebinkernel(k70,pixsk70,pixs70)
    ensure_psf_centered,k70r
    flag70 = 0.0
    convolve_flag,p70i,k70r,flag70,imgconv70

    astrk70r = astrk70
    astrk70r.naxis = size(k100r,/dim)
    astrk70r.cd[0] = pixs70/3600.
    astrk70r.cd[3] = pixs70/3600.
    astrk70r.crpix = astrk70r.naxis/2
    hk70ir = hk70i
    putast,hk70ir,astrk70r
    ;mwrfits,k70r,'kernel_70_'+beam+'_regrid.fits',hk70ir,/create


    ;convolve 100 um image   
    ratkimg100 = pixsk100/pixs100
    if ratkimg100 gt 1.0 then x = 1 
    k100r = rebinkernel(k100,pixsk100,pixs100)
    ensure_psf_centered,k100r
    flag100 = 0.0
    convolve_flag,p100i,k100r,flag100,imgconv100

    astrk100r = astrk100
    astrk100r.naxis = size(k100r,/dim)
    astrk100r.cd[0] = pixs100/3600.
    astrk100r.cd[3] = pixs100/3600.
    astrk100r.crpix = astrk100r.naxis/2
    hk100ir = hk100i
    putast,hk100ir,astrk100r
    ;mwrfits,k100r,'kernel_100_'+beam+'_regrid.fits',hk100ir,/create

    ;convolve 160 um image
    ratkimg160 = pixsk160/pixs160
    k160r = rebinkernel(k160,pixsk160,pixs160)
    ensure_psf_centered,k160r
    flag160 = 0.0
    convolve_flag,p160i,k160r,flag160,imgconv160
    
    astrk160r = astrk160
    astrk160r.naxis = size(k160r,/dim)
    astrk160r.cd[0] = pixs160/3600.
    astrk160r.cd[3] = pixs160/3600.
    astrk160r.crpix = astrk160r.naxis/2
    hk160ir = hk160i
    putast,hk160ir,astrk160r
    ;mwrfits,k160r,'kernel_160_'+beam+'_regrid.fits',hk160ir,/create

    ;convolve 250 um image
    ratkimg250 = pixsk250/pixs250
    k250r = rebinkernel(k250,pixsk250,pixs250)
    ensure_psf_centered,k250r
    flag250 = 0.0
    ;stop
    convolve_flag,s250i,k250r,flag250,imgconv250

    astrk250r = astrk250
    astrk250r.naxis = size(k250r,/dim)
    astrk250r.cd[0] = pixs250/3600.
    astrk250r.cd[3] = pixs250/3600.
    astrk250r.crpix = astrk250r.naxis/2
    hk250ir = hk250i
    putast,hk250ir,astrk250r
    ;mwrfits,k250r,'kernel_250_'+beam+'_regrid.fits',hk250ir,/create

    ;convolve 350 um image
    if beam eq '500' then begin
       
       ratkimg350 = pixsk350/pixs350
       k350r = rebinkernel(k350,pixsk350,pixs350)
       ensure_psf_centered,k350r
       flag350 = 0.0
       convolve_flag,s350i,k350r,flag350,imgconv350
       
       astrk350r = astrk350
       astrk350r.naxis = size(k350r,/dim)
       astrk350r.cd[0] = pixs350/3600.
       astrk350r.cd[3] = pixs350/3600.
       astrk350r.crpix = astrk350r.naxis/2
       hk350ir = hk350i
       putast,hk350ir,astrk250r
       ;mwrfits,k350r,'kernel_350_'+beam+'_regrid.fits',hk350ir,/create
       
    endif

    ;convolve 870 um image
;    ratkimg870 = pixsk870/pixs870
;    print,pixsk870,pixs870
;    k870r = rebinkernel(k870,pixsk870,pixs870)
;    ensure_psf_centered,k870r
;    flag870 = 0.0
;    convolve_flag,a870i,k870r,flag870,imgconv870

;    astrk870r = astrk870
;    astrk870r.naxis = size(k870r,/dim)
;    astrk870r.cd[0] = pixs870/3600.
;    astrk870r.cd[3] = pixs870/3600.
;    astrk870r.crpix = astrk870r.naxis/2
;    hk870ir = hk870i
;    putast,hk870ir,astrk870r
    
;    ;convolve the submm data using a Gaussian beam aproximation.
;    if beam eq '500' then beamh = (1711.0*4.*ALOG(2)/!Pi)^0.5 else $
;       if beam eq '250' then beamh = (816.0*4.*ALOG(2)/!Pi)^0.5
;    
;    fwhm_870_500 = (beamh^2. - 19.2^2.)^0.5 ;["]
;    fwhm_870_500_pix = fwhm_870_500/pixs870
;    pixpbeam870 = beamh/pixs870 ;pixels per beam
;    npixelpsf870 = long([pixpbeam870,pixpbeam870])*5
;    psf_870_500 = psf_gaussian(NPIXEL = npixelpsf870,FWHM = fwhm_870_500_pix,$
;                               /double,/normalize)
;    flag870 = 0.0
;    convolve_flag,a870i,psf_870_500,flag870,imgconv870
;    
;    ;pick a reference header image: 
;    ;combine 500 um beam w/ pixscale at 350 um = 10"/pix 
;    ;*** OR *** 
;    ;combine 350 um beam w/ pixscale at 250 = 6"/pix
    
    if beam eq '500' then begin
       imgreforg = imgconv350
       imgrefhdrorg = hs350
       out_fn_info_p = out_fn_info_p+'_bin10.fits' 
       out_fn_info_s = out_fn_info_s+'_bin10.fits' 
;       out_fn_info_sm = out_fn_info_sm+'_bin10.fits'
    endif else if beam eq '350' then begin
       imgreforg = imgconv250
       imgrefhdrorg = hs250
       out_fn_info_p = out_fn_info_p+'_bin06.fits' 
       out_fn_info_s = out_fn_info_s+'_bin06.fits' 
;       out_fn_info_sm = out_fn_info_sm+'_bin06.fits'
    endif
    size_img_ref = size(imgreforg,/dim)
    size_ref = min(size_img_ref)
    ;may need to adjust this param:
    trim_edge = 20L ;because the SPIRE images are larger than PACS.  
    xsize = floor(size_ref/2.) - trim_edge  &  ysize = xsize
    center = size_img_ref/2
    hextract,imgreforg,imgrefhdrorg,imgref,imgrefhdr,center[0]-xsize,center[0]+xsize,$
             center[1] - ysize,center[1] + ysize;,/silent
    sxdelpar,imgrefhdr,['PLANE1','PLANE2','PLANE3','PLANE4']
    history = ['AMS: Image convolution to the '+beam+' um beam','AMS: performed on '+systime()]
    sxaddhist, history,imgrefhdr
    sxaddpar,imgrefhdr,'BUNIT','Jy/arcsec^2'

    hastrom, imgconv70, hp70, p70s, hp70s, imgrefhdr, MISSING =0, INTERP = 2    
    hastrom, imgconv100, hp100, p100s, hp100s, imgrefhdr, MISSING =0, INTERP = 2
    hastrom, imgconv160, hp160, p160s, hp160s, imgrefhdr, MISSING =0, INTERP = 2

    hdr70 = imgrefhdr
    sxaddpar,hdr70,'FILTER','70'
    sxaddpar,hdr70,'INSTRUMENT','PACS'
    mwrfits,p70s,path_out+out_fn_obj+'0070_'+out_fn_info_p,hdr70,/create

    hdr100 = imgrefhdr
    sxaddpar,hdr100,'FILTER','100'
    sxaddpar,hdr100,'INSTRUMENT','PACS'
    mwrfits,p100s,path_out+out_fn_obj+'0100_'+out_fn_info_p,hdr100,/create
    
    hdr160 = imgrefhdr
    sxaddpar,hdr160,'FILTER','160'
    sxaddpar,hdr160,'INSTRUMENT','PACS'
    mwrfits,p160s,path_out+out_fn_obj+'0160_'+out_fn_info_p,hdr160,/create

    if beam eq '500' then begin

       hastrom, imgconv250, hs250, s250s, hs250s, imgrefhdr, MISSING =0, INTERP = 2
       hastrom, s500i, hs500, s500s, hs500s, imgrefhdr, MISSING =0, INTERP = 2
       
       hdr250 = imgrefhdr
       sxaddpar,hdr250,'FILTER','250'
       mwrfits,s250s,path_out+out_fn_obj+'0250_'+out_fn_info_s,hdr250,/create
       
       mwrfits,imgref,path_out+out_fn_obj+'0350_'+out_fn_info_s,imgrefhdr,/create
       
       hdr500 = imgrefhdr
       sxaddpar,hdr500,'FILTER','500'
       mwrfits,s500s,path_out+out_fn_obj+'0500_'+out_fn_info_s,hdr500,/create

    endif else if beam eq '350' then begin

       hastrom, s350i, hs350, s350s, hs350s, imgrefhdr, MISSING =0, INTERP = 2
       
       mwrfits,imgref,path_out+out_fn_obj+'0250_'+out_fn_info_s,imgrefhdr,/create
       
       hdr350 = imgrefhdr
       sxaddpar,hdr350,'FILTER','350'
       mwrfits,s350s,path_out+out_fn_obj+'0350_'+out_fn_info_s,hdr350,/create

    endif

    sxaddpar,imgrefhdr,'QTTY____','Jy/arcsec^2'

;    hastrom, imgconv870, ha870, a870s, ha870s, imgrefhdr, MISSING =0, INTERP = 2

;    hdr870 = imgrefhdr
;    sxaddpar,hdr870,'FILTER','870'
;    mwrfits,a870s,path_out+out_fn_obj+'0870_'+out_fn_info_s,hdr870,/create
    
;    stop 
  return
end
