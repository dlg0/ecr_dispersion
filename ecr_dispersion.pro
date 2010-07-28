;+
;	Calculate the kPer or kPar for waves
;	in the electron cyclotron range based 
;	on Brambilla.
;
;	David L Green, ORNL, 28-Jul-2010
;	greendl1@ornl.gov
;-


pro ecr_dispersion, $
	kPar = kPar, $
	kPer = kPer, $
	SI = SI 

	if (not keyword_set(kPar)) $
			and (not keyword_set(kPer)) then begin

		print, 'Input error:'
		print, 'Either kPer or kPar must be set.'
		print, 'For kPer=0 or kPar=0 just use 1e-6 or the like ;)'

	endif


	; parameters
	; ----------

	f	= 28.0d9	; Hz


	if not keyword_set ( SI ) then begin

		; Stupid Brambilla units (cgs)
		; ----------------------------

		print, ''
		print, '*** Using CGS units'
		print, '*** Ensure inputs are appropriate'
		print, ''

		ne_	= 1.0d12	; cm^-3
		b0	= 1.1d4		; Gauss = 1d4*Tesla 

		print, 'ne: ', ne_, ' [cm^-3]'
		print, 'B: ', b0, ' [Gauss]'
		print, ''

		e 	= 4.8032d-10
		Z	= -1d0
		q	= e * Z
		me	= 9.109d-28
		c	= 2.99792d10

		wce	=	q * b0 / ( me * c ) 
		wpe	= sqrt ( 4 * !pi * ne_ * q^2 / me ) 

		units = ' [cm^-1]'

	endif else begin

		; SI units
		; --------

		print, ''
		print, '*** Using SI units'
		print, '*** Ensure inputs are appropriate'
		print, ''

		ne_	= 1.0d18	; m^-3
		b0	= 1.1 		; T

		print, 'ne: ', ne_, ' [m^-3]'
		print, 'B: ', b0, ' [T]'
		print, ''

		e 	= 1.602d-19
		Z	= -1d0
		q	= e * Z
		me	= 9.109d-31	
		e0 	= 8.85d-12
		c	= 2.99792d8

		wce = q * b0 / me
		wpe	= sqrt ( ne_ * q^2 / ( me * e0 ) ) 

		units = ' [m^-1]'

	endelse

	w	= 2 * !pi * f


	; From pg 203-207 Brambilla
	; -------------------------

	X 	= wpe^2 / w^2
	u 	= wce / w

	if keyword_set(kPar) then begin

		nPar = kPar * c / w

		del1	= u^2 * ( 2 - nPar^2 )^2 + 4 * nPar^2 * ( 1-X )

		if del1 lt 0 then begin
			print, 'Error: del1 < 0'
			stop
		endif

		nPerSq_O	= 1 - nPar^2 - X $
				- X*u/2 * ( u*(1+nPar^2)-sqrt(del1) ) / ( 1 - X - u^2 )
		nPerSq_X	= 1 - nPar^2 - X $
				- X*u/2 * ( u*(1+nPar^2)+sqrt(del1) ) / ( 1 - X - u^2 )

		print, 'kPer^2_O: +/-', nPerSq_O*w^2/c^2 
		print, 'kPer^2_X: +/-', nPerSq_X*w^2/c^2

		print, 'kPer_O: +/-', sqrt(nPerSq_O*w^2/c^2), units
		print, 'kPer_X: +/-', sqrt(nPerSq_X*w^2/c^2), units

	endif

	if keyword_set(kPer) then begin

		nPer = kPer * c / w

		del2	= ( 2 * ( 1-X ) - nPer^2 )^2 + nPer^4 * ( u^2 - 1 )
		
		if del2 lt 0 then begin
			print, 'Error: del2 < 0'
			stop
		endif


		nParSq_L	= 1.0 / ( 2 * ( 1-u^2 ) ) $
			   	* ( 1-X-u^2+X*u * ( nPer^2*u + sqrt(del2) ) / (1-X) )	
	   	nParSq_R	= 1.0 / ( 2 * ( 1-u^2 ) ) $
				* ( 1-X-u^2+X*u * ( nPer^2*u - sqrt(del2) ) / (1-X) )

		print, 'nPar^2_L: +/-', nParSq_L*w^2/c^2 
		print, 'nPar^2_R: +/-', nParSq_R*w^2/c^2

		print, 'kPar_L: +/-', sqrt(nParSq_L*w^2/c^2), units
		print, 'kPar_R: +/-', sqrt(nParSq_R*w^2/c^2), units

	endif


end
