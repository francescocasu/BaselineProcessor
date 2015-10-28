; please, intend this as a pesudo-code

pro OrbSave, stateVectNum,tstart,statevectime,statevectpos,stateMin=stateMin

; input
; stateVectNum: number of available state vectors
; tstart: start time of the acquisition in seconds from 00:00
; statevecttime: array of state vector times (not needed in array form, indeed), in seconds from 00:00
; statevectpos: array of state vector positions and velocity
; stateMin: minumum nuber of required state vectors

; output
; master: matrix of polynomial coefficients that interpolate the state vectors
; par: parameters structure 

if keyword_set(stateMin) eq 0 then stateMin = 5l

stateVectNumCurr = min([stateMin,stateVectNum]) ;this is not mandatory. One can choose the appropriate number of state vectors: they should be at least 5

; following 6 lines are needed to identify, among the availables, the
; state vectors that span across the considered image

tmedio=(tend-tstart)/2.d0+tstart 
indexNearest = round((tmedio-STATEVECTTIME(0))/(STATEVECTTIME(1)-STATEVECTTIME(0)))
indexStart = indexNearest - round((stateVectNumCurr-1l)/2.)
if indexStart lt 0 then indexStart = 0l
indexEnd = indexNearest + round((stateVectNumCurr-1l)/2.)
if indexEnd ge stateVectNum then indexEnd = stateVectNum-1l

svpp=STATEVECTPOS(*,indexStart:indexEnd) ; array of state vectors (not needed in array form): the important thing is to properly identify the statevectors of interest across the current image
svx=transpose(svpp(0,*))
svy=transpose(svpp(1,*))
svz=transpose(svpp(2,*))
svtp=STATEVECTTIME(indexStart:indexEnd) ;state vectors times corresponding to the previously selected state vectors
;stop
grado = leggiParam(pathgen+'/'+paramFileIn,'GRADORB') ; degree of the polynomial: read from external parameter files

OrbitCalculatorLS,grado,svtp,svx,svy,svz,master,par ;state vector interpolator


tstart = tstart - par.t_mid_m ; this line permits to refer the start time of the "image" to the reference time of the polynomial that interpolates the state vectors


; fills a parameter structure: to check if all of these values are
; needed later (they come from the S1A file). Not all of them are
; present as input in this software template

par.tstart = tstart
par.t_delay = rdel(0)           ;range delay
par.c      = c                  ;velocita' della luce
par.freq   = CARRIERFREQ        ;frequenza
par.lambda = LAMBDA             ;lunghezza d'onda
par.fsamp  = fc
par.fsamp_hz  = fc              ;frequenza di campionamento
par.prf    = prf                ;prf
par.center_lat = leggiParam(pathgen+'/'+paramFileIn,'LATCM') 
par.center_lon = leggiParam(pathgen+'/'+paramFileIn,'LONCM')
par.ronear    =  par.t_delay*par.c/2.0d0 ; near-range    [m]	
par.deltaro   =  ((1.d0/par.fsamp)*par.c/2.0d0) ; range spacing [m]	
par.deltatime = (1.d0/par.prf) ; azimuth-time spacing [s]		

par=create_struct(par,'TX_PULSE_LEN_S',tau) ;to check if needed

m=master
save,file=pathMaster+'/orbite_'+sensor+'.sav',par,m ; done by IDL but perhaps not needed

end
