PHY2371.pdf : PHY2371.tex bratu.pdf heat_0.pdf \
	heat_1.pdf heat_xl_1.pdf heat_xl_2.pdf \
	oscillateur_energy1.pdf oscillateur_energy2.pdf \
	oscillateur_euler.pdf radio_leapfrog.pdf \
	radio.pdf error.pdf heat_3d.pdf  \
	error_withRK4.pdf \
	bvp_prob1.pdf bvp_prob2.pdf bvp_error.pdf \
	rounding_error.pdf upwind.pdf downwind.pdf compare_up_lax.pdf \
	PHY2371.f0.table PHY2371.f1.table PHY2371.f2.table
	pdflatex PHY2371.tex

pdf : PHY2371.tex
	pdflatex PHY2371.tex

bratu.pdf : bvp.py
	python bvp.py

bvp_prob1.pdf : bvp.py
	python bvp.py

bvp_prob2.pdf : bvp.py
	python bvp.py

bvp_error.pdf : bvp.py
	python bvp.py

oscillateur_energy1.pdf : oscillateur.py
	python oscillateur.py


oscillateur_energy2.pdf : oscillateur.py
	python oscillateur.py


oscillateur_euler.pdf : oscillateur.py
	python oscillateur.py

heat_0.pdf : heat.py
	python heat.py

heat_1.pdf : heat.py
	python heat.py

heat_xl_1.pdf : heat.py
	python heat.py

heat_xl_2.pdf : heat.py
	python heat.py

radio.pdf : radiodecay.py
	python radiodecay.py

radio_leapfrog.pdf : radiodecay.py
	python radiodecay.py

rounding_error.pdf : radiodecay.py
	python radiodecay.py


error.pdf : radiodecay.py
	python radiodecay.py

error_withRK4.pdf : radiodecay.py
	python radiodecay.py

upwind.pdf : advect.py
	python advect.py

downwind.pdf : advect.py
	python advect.py

compare_up_lax.pdf : advect.py
	python advect.py

PHY2371.f2.table : PHY2371.f2.gnuplot
	gnuplot PHY2371.f2.gnuplot


PHY2371.f0.table : PHY2371.f0.gnuplot
	gnuplot PHY2371.f0.gnuplot


PHY2371.f1.table : PHY2371.f1.gnuplot
	gnuplot PHY2371.f1.gnuplot


make clean: 
	rm -f *.pdf

all: PHY2371.pdf
