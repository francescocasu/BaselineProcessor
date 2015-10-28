; *****************************************************************
; Funzione di sincronizzazione vettoriale
; funziona sia con dati scalari che vettoriali

; For a given ground point, the zero doppler orbit time is returned. 
; ******************************************************************
;; Input:
;; t0 = starting time set by default equal to zero
;; x, y, z = ground point coordinates 
;; master = coefficients of the polynomy (of a given degree) which 
;;          interpolates the orbit state vectors, in matrix form. 
;;         ( x   y   z  vx  vy  vz  ax  ay  az  )
;;         | x2  y2  z2 vx2 vy2 vz2 ax2 ay2 az2 |
;;         | .   .   .  .   .   .   .   .   .   |
;;         ( xn  yn  zn vxn vyn vzn axn ayn azn )
;;         where x,y,z are the positions, vx, vy, vz are the
;;         velocities and ax, ay and az are the accelerations.
;******************************************************************

function sincro,t0,x,y,z,master,  $
                convergence=conv, $  ;boolean output: 0 not converged, 1 converged
                maxiter=maxit,    $  ;number of iterations before exiting
                tolerance_t=tol_t,$  ;convergence tolerance threshold 
                max_err=max_err      ;output: accuracy of the retrieved position 

if (n_elements(maxit) eq 0) then maxit=100
if (n_elements(tol_t) eq 0) then tol_t=1.d-8


t=t0

for ii=1,maxit do begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;Evaluation of the position, velocity and acceleration on the;;;;;;
;;;;;orbit in correspondence of time t through dedicated functions;;;;; 
;;;;;;;;;;;;;;;;;(sx, sy, sz, vx, vy, vz, ax, ay, az).;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        sxm=sx(t,master)
        sym=sy(t,master)
        szm=sz(t,master)

        vxm=vx(t,master)
        vym=vy(t,master)
        vzm=vz(t,master)

        axm=ax(t,master)
        aym=ay(t,master)
        azm=az(t,master)


        delta=-(vxm*(x-sxm)+vym*(y-sym)+vzm*(z-szm))/ $
               (axm*(x-sxm)+aym*(y-sym)+azm*(z-szm)-vxm^2d0-vym^2d0-vzm^2d0)

        t=t+delta

        if max(abs(delta)) lt tol_t then goto, fine

endfor
conv=1

fine:
max_err=max(abs(delta))

return, t

end
