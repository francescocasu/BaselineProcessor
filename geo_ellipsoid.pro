;; This function provides the geocentric cartesian coordinates of
;; a point located at a defined time along the orbit and at a given
;; range.
;
;; Input: 
;; t = (relative) time along the reference orbit in seconds with
;; respect to the time axis of the polynomy describing the orbit. 
;; ro1 = range of the point of interest
;; master = coefficients of the polynomy (of a given degree) which 
;;          interpolates the orbit state vectors of the image, in matrix form. 
;;         ( x   y   z  vx  vy  vz  ax  ay  az  )
;;         | x2  y2  z2 vx2 vy2 vz2 ax2 ay2 az2 |
;;         | .   .   .  .   .   .   .   .   .   |
;;         ( xn  yn  zn vxn vyn vzn axn ayn azn )
;;         where x,y,z are the positions, vx, vy, vz are the
;;         velocities and ax, ay and az are the accelerations.
;; x0, y0, z0 = coordinates of a ground point close to the point of
;;              interest
;; 
;; Output:
;; [x,y,z] = vector containing the cartesian coordinates of the point
;; of interest   
  

function geo_ellipsoid,t,ro1,master,x0,y0,z0,  $
                convergence=conv, $  ;boolean output: 0 not converged, 1 converged
                maxiter=maxit,    $  ;number of iterations before exiting
                tolerance_s=tol_s,$  ;convergence tolerance thredhold 
                max_err=max_err,  $  ;output: accuracy of the retrieved position 
		par_map=par_map      ;ellipsoid axes

if (n_elements(maxit) eq 0) then maxit=100
if (n_elements(tol_s) eq 0) then tol_s=1.d-4

if (n_elements(par_map) eq 0) then begin
; WGS84
     a=6378137.000d0
     b=6356752.314d0
; GEM6
;     a = 6378.1440000d3
;     b = 6356.7590000d3
endif else begin
     a=par_map.a
     b=par_map.b
endelse


x=x0
y=y0
z=z0

aa=dblarr(3,3)
f=dblarr(1,3)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;Evaluation of the position and vlocity on the orbit in;;;;;;;;;;;
;;correspondence of time t through dedicated functions (sx, sy, sz,;;;
;;vx, vy, vz). These functions in our case can be substituted by the;;
;;;;;;;;;;;;;;;;;;;values of annotated state vectors.;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
sxm=sx(t,master)
sym=sy(t,master)
szm=sz(t,master)

vxm=vx(t,master)
vym=vy(t,master)
vzm=vz(t,master)

; Evaluation of the velocity vector square module
 
vmod=sqrt(vxm^2+vym^2+vzm^2) ; ^2 represents the power of 2

; Normalization of the velocity components

vxm = vxm/vmod
vym = vym/vmod
vzm = vzm/vmod

; 

for ii=1,maxit do begin

   ro1_temp = sqrt((x-sxm)^2+(y-sym)^2+(z-szm)^2) 
   rr_app   = sqrt(x^2+y^2+(a*z/b)^2)

; Filling the system matrix
; Note that in IDL the first dimension represts columns while the
; second dimension represents rows.
 
   aa(0,0) = (x-sxm)/ro1_temp
   aa(1,0) = (y-sym)/ro1_temp
   aa(2,0) = (z-szm)/ro1_temp

   aa(0,1) = vxm
   aa(1,1) = vym
   aa(2,1) = vzm

   aa(0,2) = x/rr_app
   aa(1,2) = y/rr_app
   aa(2,2) = (a/b)^2*z/rr_app


   f(0) = ro1_temp-ro1
   f(1) = (x-sxm)*vxm+(y-sym)*vym+(z-szm)*vzm
;   f(2) = ((x^2+y^2)/a^2)+(z^2/b^2)-1
   f(2) = rr_app-a

;The INVERT function uses the Gaussian elimination method to compute the inverse of the square matrix aa.
   aai = invert(aa,stato)  
   if (stato NE 0) then print,'Matrix inverstion problems. State =', stato

   delta = -aai##f ; the ## operator computes array elements by multiplying the rows of the first array by the columns of the second array

   x=x+delta(0)
   y=y+delta(1)
   z=z+delta(2)


   if max(total(abs(delta))) lt tol_s then goto, fine ; the function "total" computes the sum of the elements of its argument. "abs" evaluates the absolute value

endfor
conv=1

fine:
max_err=max(abs(delta))
;print,'Numero iterazioni in geo_ellipsoid :',ii-1
return, [x,y,z]

end
