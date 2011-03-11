$size = 100;

#====================================================================

$group = 4;
$part = $group / 2;
$index = 1;
$blockBase = $size*$size;
$partBase = $blockBase*2 + 2;

$nodeNum= $blockBase*$group + $group + 2;
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

$gVal = "100";
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

for($partIndex=0;$partIndex<$part;$partIndex=$partIndex+1)
{

	for($k=0;$k<2;$k=$k+1) # used to set 2 blocks
	{
		$base = $partIndex*$partBase + $k*$blockBase;

		# horizontal
		for($i=1;$i<=$size;$i=$i+1)
		{
			for($j=1;$j<$size;$j=$j+1)
			{
				$pos = $base + ($i-1)*$size + $j;
				$neg = $pos+1;
				printf("G%d %d %d %s %s %s %s\n",$index,$pos,$neg,$gVal,$gAlpha,$gBeta,$gGamma);
				$index = $index+1;
			}
		}

		# vertical
		for($i=1;$i<$size;$i=$i+1)
		{
			for($j=1;$j<=$size;$j=$j+1)
			{
				$pos = $base + ($i-1)*$size + $j;
				$neg = $pos + $size;
				printf("G%d %d %d %s %s %s %s\n",$index,$pos,$neg,$gVal,$gAlpha,$gBeta,$gGamma);
				$index = $index+1;
			}
		}
	}

	# connecting 2 blocks

	$dummy1 =  $partIndex*$partBase + 2 * ($blockBase) + 1;
	$dummy2 = $dummy1 + 1;

	$pos = $partIndex*$partBase + $blockBase;
	$neg = $dummy1;
	printf("G%d %d %d %s %s %s %s\n",$index,$pos,$neg,$gVal,$gAlpha,$gBeta,$gGamma);
	$index 	= $index+1;

	$pos = $partIndex*$partBase + 2*$blockBase;
	$neg = $dummy2;
	printf("G%d %d %d %s %s %s %s\n",$index,$pos,$neg,$gVal,$gAlpha,$gBeta,$gGamma);
	$index 	= $index+1;

	$pos = $dummy1;
	$neg = $dummy2;
	printf("G%d %d %d %s %s %s %s\n",$index,$pos,$neg,$gVal,$gAlpha,$gBeta,$gGamma);
	$index 	= $index+1;
}


# to connect the 2 parts
$dummy1 = 2*$partBase+1;
$dummy2 = $dummy1 + 1;

$pos = $partBase;
$neg = $dummy1;
printf("G%d %d %d %s %s %s %s\n",$index,$pos,$neg,$gVal,$gAlpha,$gBeta,$gGamma);
$index 	= $index+1;

$pos = 2*$partBase;
$neg = $dummy2;
printf("G%d %d %d %s %s %s %s\n",$index,$pos,$neg,$gVal,$gAlpha,$gBeta,$gGamma);
$index 	= $index+1;

$pos = $dummy1;
$neg = $dummy2;
printf("G%d %d %d %s %s %s %s\n",$index,$pos,$neg,$gVal,$gAlpha,$gBeta,$gGamma);
$index 	= $index+1;



for($i=1;$i<=$nodeNum;$i=$i+1)
{
	printf("C%d %d 0 %s %s %s %s\n",$i,$i,$cVal,$cAlpha,$cBeta,$cGamma);
}


