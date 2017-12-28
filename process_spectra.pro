
PRO process_spectra

  ;; example of usage:
  ;; res = qsfit('data/spec-0752-52251-0323.fits', z=0.3806, ebv=0.06846, outdir='output')
  ;; qsfit_plot, res, filename='plot\spec0752-52251-0323'

  ;; reading the data
  data = read_csv('QSO_name.csv')
  help, /struct, data
  
  ;; error handling if file for plotting doesn't exist
  for i = 0, 0 do begin
    catch, error
    if error ne 0 then begin
      catch, /cancel
      print, 'An error occured!'
      continue
    ENDIF
    
    ;; fitting the spectrum
    res = qsfit('data/'+data.field1[i], z=data.field2[i], ebv=data.field3[i], outdir='output')
    
    ;; flatten the data structure, then save to file
    flat = qsfit_flatten_results(res)
    gps, flat, out=table_result
    write_csv, 'table/'+data.field1[i]+'.txt', table_result
        
    ;; rebin and plotting, but without residual
    res.gfit.plot.(0).main.rebin = 3
    qsfit_plot, res, filename=data.field1[i], resid=0

  endfor


END