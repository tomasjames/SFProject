pro congrid_her,irdc
;
; Butchered code from Dr. Sarah Ragan to assess regrid and convolution in IDL
;


; Define the directory containing the data
dir = '/export/home/c1158976/regrid'

; Define the directory containing the PSFs
bdir='/export/home/c1158976/regrid'

	; Define the beam
	beam = '250'

    ;out dir for convolved data, if left blank then = current directory
    path_out = dir+'/conv/'

    ;set up a file naming convention here (this is only a suggestion):
    ;Obj name + wavelength + reduction info + convol. wavelength + pixel scale
    out_fn_obj = irdc
    out_fn_info_p = 'conv'+beam 
    out_fn_info_s = 'conv'+beam ; [L2], [BM], [SC]...
    out_fn_info_sm = 'gauss_beam'+beam  

    ;set data and kernel paths
    pathp = dir+irdc+'/'             ; pacs data path
    paths = dir+irdc+'/'             ; spire data path
    ;pathk = '/disk1/ragan/kernel/'

    ;set data file names to be read in
    s250_fn = 'HorseHead.fits'

    ; Read in the data
    s250i = mrdfits(paths+s250_fn,0,hs250,/silent)
    sxdelpar,hs250,'XTENSION'

    ;read in PSF files [GA's kernels are NOT normalized]
       k250 = mrdfits('theoretical_spire_beam_model_psw_V0_2.fits',0,hk250i,/silent)


    ;***************************************************************
    ;brute-force NAN-removal, in case some are still in the images...
    nani = where(finite(s250i,/nan) eq 1L,nnan) 
    if nani[0] ne -1L then s250i[nani] = 0.0
    
    ;extract astrometry from fits files
    extast,hs250,astr250

    ;calculate pixel scales
    pixs250 = (astr250.cdelt[0]^2.)^0.5*3600. ;["/pix]

    ;extract astrometry, pixel scales from PSF headers. 
    extast,hk250i,astrk250

	;calculate pixel scales
    pixsk250 = (astrk250.cd[0]^2. + astrk250.cd[1]^2.)^0.5*3600.;["/pix]


    ;PACS native units: Jy/pix and the point-source -> extended source corrections
    s250i = s250i / 423.0 ;* 0.9828 ;[Jy/arcsec^2]   

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


    if beam eq '250' then begin

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
