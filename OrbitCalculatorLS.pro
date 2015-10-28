pro OrbitCalculatorLS,N_ord,tim,x,y,z,master,par

; input
; N_ord: polynomial order
; tim: array of state vector times
; x,y,z: arrays of state vector positions
; 
; output
; master: matrix of polynomial coefficients that interpolate the state vectors
; par: parameter structure

;*********************************
;* Dichiarazione delle variabili *
;*********************************

t=tim

N_stat=n_elements(x)

if (N_stat le N_ord) then begin
  print,'Numero State Vector insufficienti'
  goto,eend
endif

master=dblarr(9,N_ord+1)

par={c:0.d0,freq:0.d0,lambda:0.d0,fsamp:0.d0,fsamp_hz:0.d0,prf:0.d0, $
     t_delay:0.d0,abs_tstar:0d0,t0m:0.d0,state_spac_m:0.d0, $
     ronear:0.d0,deltaro:0.d0,tstart:0.d0,deltatime:0.d0,t_mid_m:0.d0, $
     center_lat:0.d0,center_lon:0.d0,date:'',$
     n_stat: N_stat, $          ;Numero State Vector 
     t:  t,           $         ;Istanti degli State Vector
     x:  x,           $         ;  |
     y:  y,           $         ;  |  State Vector
     z:  z }                    ;  |



;parametri della scena elaborata

; ---------------------------------------

par.abs_tstar = t(0)		; start_time assoluto (master)

par.t0m       = t(0)	 	 	; tempo del primo state vector

par.state_spac_m = t(1)-t(0)	;time spacing state vectors master

; ---------------------------------------


par.t_mid_m=total(t)/n_elements(t); Tempo medio

par.tstart=par.abs_tstar-par.t_mid_m

   
; Mette l'origine nel tempo medio degli state vectors


t =t -par.t_mid_m

 master(0,*)=reform(poly_fit(t,x,N_ord)) ; the IDL poly_fit function permits to fit a number of input (x) at given positions (t) with a polynoial of N_ord degree: the output is an array containing the polynomial coefficients (N_ord+1)
 master(1,*)=reform(poly_fit(t,y,N_ord))
 master(2,*)=reform(poly_fit(t,z,N_ord))
 for i=1,N_ord do master(3,i-1)=double(i)*master(0,i)
 for i=1,N_ord do master(4,i-1)=double(i)*master(1,i)
 for i=1,N_ord do master(5,i-1)=double(i)*master(2,i)
 for i=1,N_ord-1 do master(6,i-1)=double(i)*master(3,i)
 for i=1,N_ord-1 do master(7,i-1)=double(i)*master(4,i)
 for i=1,N_ord-1 do master(8,i-1)=double(i)*master(5,i)

eend:

return
end
