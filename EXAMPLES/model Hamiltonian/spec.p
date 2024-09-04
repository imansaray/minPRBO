# gnuplot script to plot the data for complex energy 
# paper

#path0="/wrk/inm4/Si/full_BSE_gen_herm/"
path0="/mnt/c/Users/inm4/OneDrive - NIST/Desktop/Publications/minPRBO/model_hamiltonian/serial/"


# To save to a png file
set terminal pngcairo enhanced


# Linestyles

# --- black ---
set style line 1 \
    linecolor rgb 'black'\
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 0.75

set style line 11 \
    linecolor rgb 'black'\
    linetype 1 linewidth 2 \
    pointtype 3 pointsize 0.75\
    dashtype 2

# --- red ---
set style line 2 \
    linecolor rgb '#ff0000'\
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 0.75

set style line 22 \
    linecolor rgb '#ff0000'\
    linetype 1 linewidth 2 \
    pointtype 3 pointsize 0.75\
    dashtype 2

# --- blue ---
set style line 3 \
    linecolor rgb '#0000ff'\
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 0.75

set style line 33 \
    linecolor rgb '#0000ff'\
    linetype 1 linewidth 2 \
    pointtype 3 pointsize 0.75\
    dashtype 2

# --- green ---
set style line 4 \
    linecolor rgb '#008000'\
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 0.75

set style line 44 \
    linecolor rgb '#008000'\
    linetype 1 linewidth 2 \
    pointtype 3 pointsize 0.75\
    dashtype 2

# --- magenta ---
set style line 5 \
    linecolor rgb '#ff00ff'\
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 0.75

set style line 55 \
    linecolor rgb '#ff00ff'\
    linetype 1 linewidth 2 \
    pointtype 3 pointsize 0.75\
    dashtype 2

# --- cyan ---
set style line 6 \
    linecolor rgb '#00ffff'\
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 0.75

set style line 66 \
    linecolor rgb '#00ffff'\
    linetype 1 linewidth 2 \
    pointtype 3 pointsize 0.75\
    dashtype 2

# --- dark cyan ---
set style line 7 \
    linecolor rgb '#008b8b'\
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 0.75

set style line 77 \
    linecolor rgb '#008b8b'\
    linetype 1 linewidth 2 \
    pointtype 3 pointsize 0.75\
    dashtype 2

# --- lime green ---
set style line 8 \
    linecolor rgb '#00ff00'\
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 0.75

set style line 88 \
    linecolor rgb '#00ff00'\
    linetype 1 linewidth 2 \
    pointtype 3 pointsize 0.75\
    dashtype 2

# --- brown ---
set style line 9 \
    linecolor rgb '#a52a2a'\
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 0.75

set style line 99 \
    linecolor rgb '#a52a2a'\
    linetype 1 linewidth 2 \
    pointtype 3 pointsize 0.75\
    dashtype 2

#---- Point styles colors and sizes ------

pt_style1="pointsize 1.0 linecolor rgb '#a52a2a'"

# --- black ---
pt_style1="pointsize 1.0 linecolor rgb 'black'"

# --- red ---
pt_style2="pointsize 1.0 linecolor rgb '#ff0000'"

# --- blue ---
pt_style3="pointsize 1.0 linecolor rgb '#0000ff'"

# --- green ---
pt_style4="pointsize 1.0 linecolor rgb '#008000'"

# --- magenta ---
pt_style5="pointsize 1.0 linecolor rgb '#ff00ff'"

# --- cyan ---
pt_style6="pointsize 1.0 linecolor rgb '#00ffff'"

# --- dark cyan ---
pt_style7="pointsize 1.0 linecolor rgb '#008b8b'"

# --- lime green ---
pt_style8="pointsize 1.0 linecolor rgb '#00ff00'"

# --- brown ---
pt_style9="pointsize 1.0 linecolor rgb '#a52a2a'"

set autoscale # scale axes automatically
unset log     # remove any log-scaling
unset label   # remove any previous labelling
#set logscale y
#set logscale x
set key top left
#unset arrow


## Fig 1: comparison of exact and approx dielectric function
# Hermitian case

# (a) e1(w)
set output 'e1_exact_herm_comp_m_1000_mlanc_200.png'
set autoscale # scale axes automatically
unset log     # remove any log-scaling
unset label   # remove any previous labelling
#set logscale y
#set logscale x
set key top right
#set xrange [0:8.0]
#set yrange [-15:35]
#set yrange [0.0:10**4]
#set xlabel "Lanczos iteration"
set xlabel "{/Symbol w}" 
set ylabel "{/Symbol e}_1({/Symbol w})" 
#set ylabel "e_2"
#set ylabel "{/Symbol e}_2({/Symbol w})" 
#set ylabel "Loss of orthogonality"
set title "Herm Lanc, m=1000, mlanc=200, spect range=22.0, eta=0.1"

# Herm
fa_1="resolvent_exact_herm_m_1000_spec_range_22.dat"
fa_2="resolvent_lanc_herm_no_m_1000_mlanc_200_spec_range_22.dat"
fa_3="resolvent_lanc_herm_full_m_1000_mlanc_200_spec_range_22.dat"

plot path0.fa_1 u 1:2 t  'Exact      ' w p pointtype 3 @pt_style2,\
     path0.fa_2 u 1:2 t  'no. orthog, niter = 200' w p pointtype 6 @pt_style3,\
     path0.fa_2 u 1:2 t  'full, niter = 200' w p pointtype 12 @pt_style4,\


# non Herm
#plot path0.fa_1 u 1:2 t  'Exact      ' w p pointtype 3 @pt_style2,\
#     path0.fa_2 u 1:2 t  'non herm., no., niter = 300' w p pointtype 6 @pt_style3,\
#     path0.fa_4 u 1:2 t  'non herm., full, niter = 300' w p pointtype 10 @pt_style5,\
#     path0.fa_3 u 1:2 t  'non herm., minPRBO {/Symbol \326}p, niter = 300' w p pointtype 12 @pt_style4,\


# (b) e2(w)
set output 'e2_exact_herm_comp_m_1000_mlanc_200.png'
set autoscale # scale axes automatically
unset log     # remove any log-scaling
unset label   # remove any previous labelling
#set logscale y
#set logscale x
set key bottom right
#set xrange [0:8.0]
#set yrange [-15:35]
#set yrange [0.0:10**4]
#set xlabel "Lanczos iteration"
set xlabel "{/Symbol w}" 
#set ylabel "{/Symbol e}_1({/Symbol w})" 
#set ylabel "e_2"
set ylabel "{/Symbol e}_2({/Symbol w})" 
#set ylabel "Loss of orthogonality"
set title "Herm Lanc, m=1000, mlanc=200, spect range=22.0, eta=0.1"


plot path0.fa_1 u 1:3 t  'Exact      ' w p pointtype 3 @pt_style2,\
     path0.fa_2 u 1:3 t  'no. orthog, niter = 200' w p pointtype 6 @pt_style3,\
     path0.fa_2 u 1:3 t  'full, niter = 200' w p pointtype 12 @pt_style4,\

# c, Loss of orthogonality 
set output 'loss_of_orthog_herm_comp_m_1000_mlanc_200.png'
set autoscale # scale axes automatically
unset log     # remove any log-scaling
#unset label   # remove any previous labelling
set logscale y
##set logscale x
set key bottom right
##set xrange [0:10]
#set yrange [1*10**-14:1]
reps=sqrt(1.11*10**-16)
set arrow from 0.0,reps to 200,reps nohead linewidth 3 linecolor rgb '#ff0000' dt 2
set label "{/Symbol \326}p" at 105.0,4.0*10**-8 tc ls 1
set xlabel "Lanczos iteration"
set ylabel "Loss of orthogonality"
set title " Loss of orthog, Herm. Lanc., m = 1000, mlanc = 200, eta = 0.1"

# Herm
fa_1="fnorm_herm_no_m_1000_mlanc_200_spec_range_22.dat"
fa_2="fnorm_herm_full_m_1000_mlanc_200_spec_range_22.dat"


# Herm
plot path0.fa_1 u 1:2 t  'no orthog.' w l linestyle 2,\
     path0.fa_2 u 1:2 t  'full orthog.' w l linestyle 3,\

## Fig 2: comparison of exact and approx dielectric function
# Non symmetric, real case

# (a) e1(w)
set output 'e1_exact_non_symm_comp_m_1000_mlanc_200.png'
set autoscale # scale axes automatically
unset log     # remove any log-scaling
unset label   # remove any previous labelling
unset arrow
#set logscale y
#set logscale x
set key bottom left
#set xrange [0:8.0]
#set yrange [-15:35]
#set yrange [0.0:10**4]
#set xlabel "Lanczos iteration"
set xlabel "{/Symbol w}" 
set ylabel "{/Symbol e}_1({/Symbol w})" 
#set ylabel "e_2"
#set ylabel "{/Symbol e}_2({/Symbol w})" 
#set ylabel "Loss of orthogonality"
set title "Non symm., m=1000, mlanc=200, spect range=22.0, eta=0.1"

# Herm
fa_1="resolvent_exact_non_symm_m_1000_spec_range_22.dat"
fa_2="resolvent_lanc_non_symm_no_m_1000_mlanc_200_spec_range_22.dat"
#fa_3="resolvent_lanc_herm_full_m_1000_mlanc_200_spec_range_22.dat"

plot path0.fa_1 u 1:2 t  'Exact      ' w p pointtype 3 @pt_style2,\
     path0.fa_2 u 1:2 t  'no. orthog, niter = 200' w p pointtype 6 @pt_style3,\
#     path0.fa_2 u 1:2 t  'full, niter = 200' w p pointtype 12 @pt_style4,\




# (b) e2(w)
set output 'e2_exact_non_symm_comp_m_1000_mlanc_200.png'
set autoscale # scale axes automatically
unset log     # remove any log-scaling
unset label   # remove any previous labelling
#set logscale y
#set logscale x
set key bottom right
#set xrange [0:8.0]
#set yrange [-15:35]
#set yrange [0.0:10**4]
#set xlabel "Lanczos iteration"
set xlabel "{/Symbol w}" 
#set ylabel "{/Symbol e}_1({/Symbol w})" 
#set ylabel "e_2"
set ylabel "{/Symbol e}_2({/Symbol w})" 
#set ylabel "Loss of orthogonality"
set title "Non symm., m=1000, mlanc=200, spect range=22.0, eta=0.1"

plot path0.fa_1 u 1:3 t  'Exact      ' w p pointtype 3 @pt_style2,\
     path0.fa_2 u 1:3 t  'no. orthog, niter = 200' w p pointtype 6 @pt_style3,\

## Fig 3: comparison of exact and approx dielectric function
# Non-Hermitian case

# (a) e1(w)
set output 'e1_exact_non_herm_comp_m_1000_mlanc_200.png'
set autoscale # scale axes automatically
unset log     # remove any log-scaling
unset label   # remove any previous labelling
#set logscale y
#set logscale x
set key top right
#set xrange [0:8.0]
#set yrange [-15:35]
#set yrange [0.0:10**4]
#set xlabel "Lanczos iteration"
set xlabel "{/Symbol w}" 
set ylabel "{/Symbol e}_1({/Symbol w})" 
#set ylabel "e_2"
#set ylabel "{/Symbol e}_2({/Symbol w})" 
#set ylabel "Loss of orthogonality"
set title "Non Herm Lanc, m=1000, mlanc=200, spect range=22.0, eta=0.1"

# Herm
fa_1="resolvent_exact_non_herm_m_1000_spec_range_22.dat"
fa_2="resolvent_lanc_non_herm_no_m_1000_mlanc_200_spec_range_22.dat"
fa_3="resolvent_lanc_non_herm_minprbo_reps_m_1000_mlanc_200_spec_range_22.dat"
fa_4="resolvent_lanc_non_herm_full_m_1000_mlanc_200_spec_range_22.dat"


# non Herm
plot path0.fa_1 u 1:2 t  'Exact      ' w p pointtype 3 @pt_style2,\
     path0.fa_2 u 1:2 t  'no. biorthog., niter = 200' w p pointtype 6 @pt_style3,\
     path0.fa_3 u 1:2 t  'minPRBO {/Symbol \326}p, niter = 200' w p pointtype 12 @pt_style4,\
     path0.fa_4 u 1:2 t  'full biorthog., niter = 200' w p pointtype 10 @pt_style5,\


# (b) e2(w)
set output 'e2_exact_non_herm_comp_m_1000_mlanc_200.png'
set autoscale # scale axes automatically
unset log     # remove any log-scaling
unset label   # remove any previous labelling
#set logscale y
#set logscale x
set key bottom right
#set xrange [0:8.0]
#set yrange [-15:35]
#set yrange [0.0:10**4]
#set xlabel "Lanczos iteration"
set xlabel "{/Symbol w}" 
#set ylabel "{/Symbol e}_1({/Symbol w})" 
#set ylabel "e_2"
set ylabel "{/Symbol e}_2({/Symbol w})" 
#set ylabel "Loss of orthogonality"
set title "Non Herm Lanc, m=1000, mlanc=200, spect range=22.0, eta=0.1"


# non Herm
plot path0.fa_1 u 1:3 t  'Exact      ' w p pointtype 3 @pt_style2,\
     path0.fa_2 u 1:3 t  'no. biorthog., niter = 200' w p pointtype 6 @pt_style3,\
     path0.fa_3 u 1:3 t  'minPRBO {/Symbol \326}p, niter = 200' w p pointtype 12 @pt_style4,\
     path0.fa_4 u 1:3 t  'full biorthog., niter = 200' w p pointtype 10 @pt_style5,\

# c, Loss of orthogonality 
set output 'loss_of_orthog_non_herm_comp_m_1000_mlanc_200.png'
set autoscale # scale axes automatically
unset log     # remove any log-scaling
#unset label   # remove any previous labelling
set logscale y
##set logscale x
set key bottom right
##set xrange [0:10]
#set yrange [1*10**-14:1]
reps=sqrt(1.11*10**-16)
set arrow from 0.0,reps to 200,reps nohead linewidth 3 linecolor rgb '#ff0000' dt 2
set label "{/Symbol \326}p" at 105.0,4.0*10**-8 tc ls 1
set xlabel "Lanczos iteration"
set ylabel "Loss of orthogonality"
set title " Loss of orthog, Non. Herm. Lanc., m = 1000, mlanc = 200, eta = 0.1"

# Non Herm
fa_1="fnorm_non_herm_no_m_1000_mlanc_200_spec_range_22.dat"
fa_2="fnorm_non_herm_minprbo_reps_m_1000_mlanc_200_spec_range_22.dat"
fa_3="fnorm_non_herm_full_m_1000_mlanc_200_spec_range_22.dat"

plot path0.fa_1 u 1:2 t  'no. biorthog.' w l linestyle 4 dt 2,\
     path0.fa_2 u 1:2 t  'minPRBO {/Symbol \326}p' w l linestyle 3 dt 4,\
     path0.fa_3 u 1:2 t  'full biorthog.' w l linestyle 2,\


## Fig 4: comparison of different flavors of PRBO, exact and approx dielectric function
# Non Herm

# (a) e1(w)
set output 'e1_exact_prbo_comp_m_1000_mlanc_200.png'
set autoscale # scale axes automatically
unset log     # remove any log-scaling
unset label   # remove any previous labelling
unset arrow
#set logscale y
#set logscale x
set key bottom left
#set xrange [0:8.0]
#set yrange [-15:35]
#set yrange [0.0:10**4]
#set xlabel "Lanczos iteration"`
set xlabel "{/Symbol w}" 
set ylabel "{/Symbol e}_1({/Symbol w})" 
#set ylabel "e_2"
#set ylabel "{/Symbol e}_2({/Symbol w})" 
#set ylabel "Loss of orthogonality"
set title "Non Herm, m=1000, mlanc=200, spect range=22.0, eta=0.1"

# Herm
fa_1="resolvent_exact_non_herm_m_1000_spec_range_22.dat"
fa_2="resolvent_lanc_non_herm_minprbo_reps_m_1000_mlanc_200_spec_range_22.dat"
fa_3="resolvent_lanc_non_herm_prbo_1_reps_m_1000_mlanc_200_spec_range_22.dat"
fa_4="resolvent_lanc_non_herm_prbo_2_reps_m_1000_mlanc_200_spec_range_22.dat"

# non Herm
plot path0.fa_1 u 1:2 t  'Exact      ' w p pointtype 3 @pt_style2,\
     path0.fa_2 u 1:2 t  'minPRBO {/Symbol \326}p, niter = 200' w p pointtype 6 @pt_style3,\
     path0.fa_3 u 1:2 t  'PRBO 1 {/Symbol \326}p, niter = 200' w p pointtype 12 @pt_style4,\
     path0.fa_4 u 1:2 t  'PRBO 2 {/Symbol \326}p, niter = 200' w p pointtype 10 @pt_style5,\




# (b) e2(w)
set output 'e2_exact_prbo_comp_m_1000_mlanc_200.png'
set autoscale # scale axes automatically
unset log     # remove any log-scaling
unset label   # remove any previous labelling
#set logscale y
#set logscale x
set key bottom right
#set xrange [0:8.0]
#set yrange [-15:35]
#set yrange [0.0:10**4]
#set xlabel "Lanczos iteration"
set xlabel "{/Symbol w}" 
#set ylabel "{/Symbol e}_1({/Symbol w})" 
#set ylabel "e_2"
set ylabel "{/Symbol e}_2({/Symbol w})" 
#set ylabel "Loss of orthogonality"
set title "Non Herm., m=1000, mlanc=200, spect range=22.0, eta=0.1"

plot path0.fa_1 u 1:3 t  'Exact      ' w p pointtype 3 @pt_style2,\
     path0.fa_2 u 1:3 t  'minPRBO {/Symbol \326}p, niter = 200' w p pointtype 6 @pt_style3,\
     path0.fa_3 u 1:3 t  'PRBO 1 {/Symbol \326}p, niter = 200' w p pointtype 12 @pt_style4,\
     path0.fa_4 u 1:3 t  'PRBO 2 {/Symbol \326}p, niter = 200' w p pointtype 10 @pt_style5,\

# c, Loss of orthogonality 
set output 'loss_of_orthog_non_herm_prbo_comp_m_1000_mlanc_200.png'
set autoscale # scale axes automatically
unset log     # remove any log-scaling
#unset label   # remove any previous labelling
set logscale y
##set logscale x
set key bottom right
##set xrange [0:10]
#set yrange [1*10**-14:1]
reps=sqrt(1.11*10**-16)
set arrow from 0.0,reps to 200,reps nohead linewidth 3 linecolor rgb '#ff0000' dt 2
set label "{/Symbol \326}p" at 105.0,4.0*10**-8 tc ls 1
set xlabel "Lanczos iteration"
set ylabel "Loss of orthogonality"
set title " Loss of orthog, Non. Herm. Lanc., m = 1000, mlanc = 200, eta = 0.1"

# Non Herm
fa_1="fnorm_non_herm_no_m_1000_mlanc_200_spec_range_22.dat"
fa_2="fnorm_non_herm_minprbo_reps_m_1000_mlanc_200_spec_range_22.dat"
fa_3="fnorm_non_herm_prbo_1_reps_m_1000_mlanc_200_spec_range_22.dat"
fa_4="fnorm_non_herm_prbo_2_reps_m_1000_mlanc_200_spec_range_22.dat"

plot path0.fa_1 u 1:2 t  'no. biorthog.' w l linestyle 4,\
     path0.fa_2 u 1:2 t  'minPRBO {/Symbol \326}p' w l linestyle 3,\
     path0.fa_3 u 1:2 t  'PRBO 1, {/Symbol \326}p' w l linestyle 5,\
     path0.fa_4 u 1:2 t  'PRBO 2, {/Symbol \326}p' w l linestyle 6,\
