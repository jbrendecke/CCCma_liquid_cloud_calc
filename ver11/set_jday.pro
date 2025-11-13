PRO set_jday,iyr,imo,iday, iddd

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; This program calculates the Julian Day from a calendar day.
;; 
;; INPUT:
;; ------------------------------
;; iyr = 4-digit year (INTEGER)
;; imo = 2-digit month (INTEGER)
;; iday = 2-digit day (INTEGER)
;; 
;; OUTPUT:
;; ------------------------------
;; iddd = Julian Day (0-365) (INTEGER)
;;
;; This is modified from YYDDD, notice the index of iddd
;; idymon[irow,imo-1], the '-1' is due to IDL index starts at 0.
;;
;; Zhe Feng, Sep 25, 2009 @ UND
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


if (n_params() ne 4) then begin
   print,'Incorrect number of arguments. Correct usage: '
   print,'set_jday,year,mon,day, jday'
   stop
endif

idymon = INTARR(2,12)
k = 0
leap = 0

idymon[0,*] = [0,31,59,90,120,151,181,212,243,273,304,334]   ;non-leap year
idymon[1,*] = [0,31,60,91,121,152,182,213,244,274,305,335]   ;leap year

leap = iyr - ( iyr/4 )*4

;Leap Year: leap = 0
if ( leap eq 0 ) then begin
  irow = 1
;Non-leap Year: leap =! 0
endif else begin
  irow = 0
endelse

iddd = idymon[irow,imo-1] + iday
;print,iddd
END
