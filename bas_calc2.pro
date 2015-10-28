;********************************************************************;
;The function prod_vett returns the vectorial product of two vectors;
;********************************************************************;

function prod_vett,a,b

c=a
c(*)=0.

c(0,*) =   a(1,*)*b(2,*)-a(2,*)*b(1,*)
c(1,*) = -(a(0,*)*b(2,*)-a(2,*)*b(0,*))
c(2,*) =   a(0,*)*b(1,*)-a(1,*)*b(0,*)

return, c
end

;******************************************************************************;
;; The baseline procedure evaluates the three components of the baseline.  
;;
;; Inputs:
;; master = coefficients of the polynomy (of a given degree) which 
;;          interpolates the reference orbit state vectors, in matrix form. 
;;         ( x   y   z  vx  vy  vz  ax  ay  az  )
;;         | x2  y2  z2 vx2 vy2 vz2 ax2 ay2 az2 |
;;         | .   .   .  .   .   .   .   .   .   |
;;         ( xn  yn  zn vxn vyn vzn axn ayn azn )
;;         where x,y,z are the positions, vx, vy, vz are the
;;         velocities and ax, ay and az are the accelerations. 
;; slave = coefficients of the polynomy (of a given degree) which 
;;         interpolates the orbit state vectors of the image, in
;;         matrix form. 
;; xtie,ytie,ztie = cartesian coordinates of the ground point
;;
;; Outputs:
;; b_perp = perpendicular baseline
;; b_parallel = parallel baseline
;; b_along = along track baseline
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro baseline,master,slave,xtie,ytie,ztie,b_perp,b_parall,b_along

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
;; The function sincro evaluates the zero doppler orbit time for a;
;; given ground point
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

tr=sincro(0.d0,xtie,ytie,ztie,master)
;t=sincro(0.d0,xtie,ytie,ztie,slave)     ; This function in our case can be substituted by the time of annotated state vectors.

;b=[sx(t,slave)-sx(tr,master), $            ; baseline vector
;   sy(t,slave)-sy(tr,master), $
;   sz(t,slave)-sz(tr,master)]

b=[slave(0)-sx(tr,master), $            ; baseline vector
   slave(1)-sy(tr,master), $
   slave(2)-sz(tr,master)]

; Define unit vectors
b_unit=b
va_unit=b
vpar_unit=b

b_unit(*)=0.
va_unit(*)=0.
vpar_unit(*)=0.
; Define unit vectors (end)

b_mod=sqrt(b(0,*)^2+b(1,*)^2+b(2,*)^2)

for k=0,2 do begin
   b_unit(k,*)=b(k,*)/b_mod(*)         ; baseline unit vector
endfor

vpar=[xtie-sx(tr,master), $             ; parallel vector
      ytie-sy(tr,master), $
      ztie-sz(tr,master)]

vpar_mod=sqrt(vpar(0,*)^2+vpar(1,*)^2+vpar(2,*)^2)

for k=0,2 do begin
   vpar_unit(k,*)=vpar(k,*)/vpar_mod(*)      ; parallel unit vector
endfor

va=[vx(tr,master), $                       ; along-track vector
    vy(tr,master), $
    vz(tr,master)]

va_mod=sqrt(va(0,*)^2+va(1,*)^2+va(2,*)^2)

for k=0,2 do begin
   va_unit(k,*)=va(k,*)/va_mod(*)     ; along-track unit vactor
endfor

vperp_unit=prod_vett(vpar_unit,va_unit)    ; perpendicular unit vector

b_perp   = b(0,*)*vperp_unit(0,*)+    $
           b(1,*)*vperp_unit(1,*)+    $
           b(2,*)*vperp_unit(2,*)
           
b_parall = b(0,*)*vpar_unit(0,*)+    $
           b(1,*)*vpar_unit(1,*)+    $
           b(2,*)*vpar_unit(2,*)
           
b_along  = b(0,*)*va_unit(0,*)+    $
           b(1,*)*va_unit(1,*)+    $
           b(2,*)*va_unit(2,*)
                 


return
end
