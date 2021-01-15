

set title "TODO--TODO"

unset mouse
set view 60,60
set xyplane 0

set title "pVector relative error [frac-of-radius]" rotate by 90
set xlabel "latitude [rad]"
set ylabel "Radius [frac-of-radius]"
# set zlabel "pVector distance error [frac-of-radius]" rotate by 90

zmap(x)=x

splot \
	  'zetaDiff.dat' u 2:4:(zmap($8)) \
		w p pt 5 ps .50 \
		lc "green" ti "2nd order" \
	, 'zetaDiff.dat' u 2:4:(zmap($6)) \
		w p pt 5 ps .50 \
		lc "red" ti "1st order" \
	;
pause -1;

