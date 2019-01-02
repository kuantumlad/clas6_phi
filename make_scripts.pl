$date = $ARGV[0];
$kin = $ARGV[1];
$index = $ARGV[2];
$nfile_start = $ARGV[3];
$nfile_end = $ARGV[4];
$data_or_sim = $ARGV[5];
$momcorr_el = $ARGV[6];
$momcorr_pi = $ARGV[7];
$cutlvl = $ARGV[8];
$cutstrictness = $ARGV[9];
$kin = "kin${kin}";
$harris_parameters = "/w/hallb-scifs17exp/clas12/bclary/clas6/harris_parameters/";
$harris_codebase = "/w/hallb-scifs17exp/clas12/bclary/clas6/harris_codebase/";

print " Creaing job ${index} for cut level ${cutlvl} with strictness ${cutstrictness} start $nfile_start end $nfile_end \n";

$s_cutstrict = "";

if( $cutstrictness == -1 ){
    $s_cutstrict = "loose";
}
elsif( $cutstrictness == 0 ){
    $s_cutstrict = "nominal";
}
elsif( $cutstrictness == 1 ){
    $s_cutstrict == "tight";
}


open( outdat2, "> jsc_cutlvl${cutlvl}_strictlvl${s_cutstrict}.jsub" ) || die("Can't create jsc file");

print outdat2
"JOBNAME: ${kin}
PROJECT: clas12
MAIL: bclary\@jlab.org
OS: centos7
TRACK: debug
SINGLE_JOB: true
MEMORY: 1 GB
OTHER_FILES:
${harris_parameters}protonMeanBetaFitParameters.txt
${harris_parameters}pipMeanBetaFitParameters.txt
${harris_parameters}kpMeanBetaFitParameters.txt
${harris_parameters}protonSigmaBetaFitParameters.txt
${harris_parameters}pipSigmaBetaFitParameters.txt
${harris_parameters}kpSigmaBetaFitParameters.txt
${harris_codebase}programFiles/\*.C
${harris_codebase}\*.C\n";

for( $j = 1; $j <= 6; $j++ ){
print outdat2
"${harris_parameters}angles_s${j}.out\n"
}
for( $s = 1; $s <= 6; $s++ ){
print outdat2
"${harris_parameters}momentum2_s${s}_c0.out
${harris_parameters}momentum2_s${s}_c1.out
${harris_parameters}momentum3_s${s}_c0.out
${harris_parameters}momentum3_s${s}_c1.out\n";
}
print outdat2
"COMMAND: /work/clas12/bclary/clas6/farm/go_myphi6.sh ${index} ${nfile_start} ${nfile_end} ${data_or_sim} ${momcorr_el} ${momcorr_pi} ${cutlvl} ${cutstrictness}
TOWORK
OUTPUT_TEMPLATE: /w/hallb-scifs17exp/clas12/bclary/clas6/results/phi6_${s_cutstrict}/@OUTPUT_DATA@\n";
close outdat2;


