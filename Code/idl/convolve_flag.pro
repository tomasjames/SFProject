pro convolve_flag,img,psf,flag,imgconv,conv = conv,flagimg = flagimg,imgflagconv=flagimgconv

  conv = convolve_sr(img,psf,/no_pad)

  flagimg = img
  flagimg[*,*] = 1.D
  lflag = where(img eq flag)
  flagimg[lflag] = 0.0
  flagimgconv = convolve_sr(flagimg,psf,/no_pad)
  imgconv = conv/flagimgconv
  imgconv[lflag] = 0.0

  return
end
