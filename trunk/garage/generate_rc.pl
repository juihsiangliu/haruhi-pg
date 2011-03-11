$nodeNum = 100000;
$endTime = "0.0000000015";
$timeStep = "0.000000000001";

$vinM = "vin_pwl_m.txt";
$vinAlpha = "vin_zero_1500_alpha.txt";
$vinBeta = "vin_zero_1500_beta.txt";
$vinGamma = "vin_zero_1500_gamma.txt";

$gAlpha = "g_alpha.txt";
$gBeta = "g_beta.txt";
$gGamma = "g_gamma.txt";

$cAlpha = "c_alpha.txt";
$cBeta = "c_beta.txt";
$cGamma = "c_gamma.txt";

$gVal = "80";
$cVal = "0.0000000000008";

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

printf("V10 1 0 vin_pwl_m.txt vin_zero_1500_alpha.txt vin_zero_1500_beta.txt vin_zero_1500_gamma.txt\n");
for($i=1;$i<$nodeNum;$i=$i+1)
{
	printf("G%d %d %d %s %s %s %s\n",$i,$i,$i+1,$gVal,$gAlpha,$gBeta,$gGamma);
}
for($i=1;$i<=$nodeNum;$i=$i+1)
{
	printf("C%d %d 0 %s %s %s %s\n",$i,$i,$cVal,$cAlpha,$cBeta,$cGamma);
}

