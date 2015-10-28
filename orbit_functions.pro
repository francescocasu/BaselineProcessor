; *************************************
; Function for orbit descriptions based
; on polynomial coefficient given in 
; the 'master' matrix
; *************************************

; Position

function sx, t,master
m=dblarr(15)
tmp=size(master)
grado=tmp(2)  ; maximum number of coefficient
m(0:grado-1)=master(0,*)
return, m(0)+m(1)*t+m(2)*t^2+m(3)*t^3+m(4)*t^4  $
                      +m(5)*t^5+m(6)*t^6+m(7)*t^7  $
		      +m(8)*t^8+m(9)*t^9+m(10)*t^10  $
		      +m(11)*t^11+m(12)*t^12+m(13)*t^13
end

function sy, t,master
m=dblarr(15)
tmp=size(master)
grado=tmp(2)
m(0:grado-1)=master(1,*)
return, m(0)+m(1)*t+m(2)*t^2+m(3)*t^3+m(4)*t^4  $
                      +m(5)*t^5+m(6)*t^6+m(7)*t^7  $
		      +m(8)*t^8+m(9)*t^9+m(10)*t^10  $
		      +m(11)*t^11+m(12)*t^12+m(13)*t^13
end

function sz, t,master
m=dblarr(15)
tmp=size(master)
grado=tmp(2)
m(0:grado-1)=master(2,*)
return, m(0)+m(1)*t+m(2)*t^2+m(3)*t^3+m(4)*t^4  $
                      +m(5)*t^5+m(6)*t^6+m(7)*t^7  $
		      +m(8)*t^8+m(9)*t^9+m(10)*t^10  $
		      +m(11)*t^11+m(12)*t^12+m(13)*t^13
end


; Velocity

function vx, t,master
m=dblarr(15)
tmp=size(master)
grado=tmp(2)
m(0:grado-1)=master(3,*)
return, m(0)+m(1)*t+m(2)*t^2+m(3)*t^3+m(4)*t^4  $
                      +m(5)*t^5+m(6)*t^6+m(7)*t^7  $
		      +m(8)*t^8+m(9)*t^9+m(10)*t^10  $
		      +m(11)*t^11+m(12)*t^12+m(13)*t^13
end

function vy, t,master
m=dblarr(15)
tmp=size(master)
grado=tmp(2)
m(0:grado-1)=master(4,*)
return, m(0)+m(1)*t+m(2)*t^2+m(3)*t^3+m(4)*t^4  $
                      +m(5)*t^5+m(6)*t^6+m(7)*t^7  $
		      +m(8)*t^8+m(9)*t^9+m(10)*t^10  $
		      +m(11)*t^11+m(12)*t^12+m(13)*t^13
end

function vz, t,master
m=dblarr(15)
tmp=size(master)
grado=tmp(2)
m(0:grado-1)=master(5,*)
return, m(0)+m(1)*t+m(2)*t^2+m(3)*t^3+m(4)*t^4  $
                      +m(5)*t^5+m(6)*t^6+m(7)*t^7  $
		      +m(8)*t^8+m(9)*t^9+m(10)*t^10  $
		      +m(11)*t^11+m(12)*t^12+m(13)*t^13
end



; Acceleration

function ax, t,master
m=dblarr(15)
tmp=size(master)
grado=tmp(2)
m(0:grado-1)=master(6,*)
return, m(0)+m(1)*t+m(2)*t^2+m(3)*t^3+m(4)*t^4  $
                      +m(5)*t^5+m(6)*t^6+m(7)*t^7  $
		      +m(8)*t^8+m(9)*t^9+m(10)*t^10  $
		      +m(11)*t^11+m(12)*t^12+m(13)*t^13 
end

function ay, t,master
m=dblarr(15)
tmp=size(master)
grado=tmp(2)
m(0:grado-1)=master(7,*)
return, m(0)+m(1)*t+m(2)*t^2+m(3)*t^3+m(4)*t^4  $
                      +m(5)*t^5+m(6)*t^6+m(7)*t^7  $
		      +m(8)*t^8+m(9)*t^9+m(10)*t^10  $
		      +m(11)*t^11+m(12)*t^12+m(13)*t^13
end

function az, t,master
m=dblarr(15)
tmp=size(master)
grado=tmp(2)
m(0:grado-1)=master(8,*)
return, m(0)+m(1)*t+m(2)*t^2+m(3)*t^3+m(4)*t^4  $
                      +m(5)*t^5+m(6)*t^6+m(7)*t^7  $
		      +m(8)*t^8+m(9)*t^9+m(10)*t^10  $
		      +m(11)*t^11+m(12)*t^12+m(13)*t^13
end
