$size = 30;
$nodeNum = 2*$size*$size + 1;
$endTime = "0.000000000005";
$timeStep = "0.000000000001";

# voltage for 1.2V
$vin12 = "vin_1_2_1500_m.txt";
$vinAlpha = "vin_zero_1500_alpha.txt";
$vinBeta = "vin_zero_1500_beta.txt";
$vinGamma = "vin_zero_1500_gamma.txt";

# voltage for clock input
$vinCLK = "vin_clk_m.txt";

# g for the power grid
$gVal = "1000";
$gAlpha = "g_alpha.txt";
$gBeta = "g_beta.txt";
$gGamma = "g_gamma.txt";

# g for the common source amplifier
$commonGm = 0.0002;
$commonGalpha = "common_g_alpha.txt";
$commonGbeta = "common_g_beta.txt";
$commonGgamma = "common_g_gamma.txt";

# table for the nmos
$nmos_m = "nmos_ivtable_np_1.2_0.01.txt";
$nmos_alpha = "zero_gv_2_241_alpha.txt";
$nmos_beta = "zero_gv_2_241_beta.txt";
$nmos_gamma = "zero_gv_2_241_gamma.txt";
$nmos_config = "nmos_config.txt";

printf("*global_var_num\n");
printf("*2\n");
printf("*global_cor_matrix\n");
printf("*1.0 0.0\n");
printf("*0.0 1.0\n");
printf("*node_num\n");
printf("*%d\n",$nodeNum);
printf("*end_time\n");
printf("*%s\n",$endTime);
printf("*time_step\n");
printf("*%s\n",$timeStep);

printf("V1 1 0 vin_1_2_1500_m.txt vin_zero_1500_alpha.txt vin_zero_1500_beta.txt vin_zero_1500_gamma.txt\n");
printf("V2 %d 0 vin_clk_m.txt vin_zero_1500_alpha.txt vin_zero_1500_beta.txt vin_zero_1500_gamma.txt\n",$nodeNum);

$index = 1;
for($i=1;$i<=$size;$i=$i+1)
{
	for($j=1;$j<$size;$j=$j+1)
	{
		$pos = ($i-1)*$size + $j;
		$neg = $pos+1;
		printf("G%d %d %d %s %s %s %s\n",$index,$pos,$neg,$gVal,$gAlpha,$gBeta,$gGamma);
		$index = $index+1;
	}
}


for($i=1;$i<$size;$i=$i+1)
{
	for($j=1;$j<=$size;$j=$j+1)
	{		$pos = ($i-1)*$size + $j;
		$neg = $pos + $size;
		printf("G%d %d %d %s %s %s %s\n",$index,$pos,$neg,$gVal,$gAlpha,$gBeta,$gGamma);
		$index = $index+1;
	}
}


for($i=1;$i<=$size*$size;$i=$i+1)
{
	$pos = $i;
	$neg = $i + $size*$size;
	printf("G%d %d %d %s %s %s %s\n",$index,$pos,$neg,$commonGm,$commonGalpha,$commonGbeta,$commonGgamma);
	$index = $index+1;
}


$index = 1;
for($i=1;$i<=$size*$size;$i=$i+1)
{
	$pos = $i + $size*$size;
	$neg = $nodeNum;
	printf("M%d %d %d 0 %s %s %s %s %s\n",$index,$pos,$neg,$nmos_m,$nmos_alpha,$nmos_beta,$nmos_gamma,$nmos_config);
	$index = $index+1;
}



