# vgs = x axis
# vds = y axis

$flag = 0;
$i = -1;

while( defined($line=<>))
{
	chomp($line);

#	"  *** source  0:vds =  700.0000m ***  ";
	if($line =~ /^\s* (\*)+ \s* source \s* 0:vds \s* = \s* \-? \d+.\d*\w? \s* (\*)+ \s*$/x)
	{
		@tokens = split(/\s+/,$line);
		foreach $token(@tokens)
		{
			if($token =~ /^ \-?\d+.\d*\w?  $/x)  #700.0000m
			{
				@numbers = split(/[A-Za-z]/,$token);  #700.0000
				@units = split(/[0-9.]/,$token); #_____m
				$vds = pop @numbers;
				$temp = pop @units;
				$ratio = 1;
				if($temp eq "m")
				{
					$ratio = 0.001;
				}
				$vds = $vds * $ratio;
				push(@vdsList,$vds);
				$i = $i+1;
			}
		}
		$flag = 1;   # enter the iv table block
	}

	# actual iv table block (ids)	
	if($flag == 1)
	{
#		0.          -50.9771p
#		100.00000m   -715.1975p
		if($line =~ /^ \s* \-? \d+.\d*\w? \s* \-? \d+.\d*\w? \s* $/x)
		{
			@tokens = split(/\s+/,$line);  # 100.00000m   -715.1975p
				
			$ids = pop @tokens;
			$ratio = 1;
			if($ids =~ /a$/)
		        {
				$ratio = 1e-15;
			}
			elsif($ids =~ /p$/)
			{
				$ratio = 1e-12;
			}
			elsif($ids =~ /n$/)
			{
				$ratio = 1e-9;
			}
			elsif($ids =~ /u$/)
			{
				$ratio = 1e-6;
			}
			elsif($ids =~ /m$/)
			{
				$ratio = 1e-3;
			}
			$ids = $ratio*$ids;

			$vgs = pop @tokens;
			$ratio = 1;
			if($vgs =~ /a$/)
		        {
				$ratio = 1e-15;
			}
			elsif($vgs =~ /p$/)
			{
				$ratio = 1e-12;
			}
			elsif($vgs =~ /n$/)
			{
				$ratio = 1e-9;
			}
			elsif($vgs =~ /u$/)
			{
				$ratio = 1e-6;
			}
			elsif($vgs =~ /m$/)
			{
				$ratio = 1e-3;
			}
			$vgs = $ratio*$vgs;

			if($i == 0)
			{
				push(@vgsList,$vgs);
			}
			push(@idsList,$ids);
		}
	}

#	" **** job concluded "
	if($line =~ /^\s* (\*)+ \s* job \s+ concluded \s*$/x)
	{
		$flag = 0;  # leave the iv table block
	}
}


$sizeVgsList = $#vgsList+1;
$sizeVdsList = $#vdsList+1;

print("$sizeVgsList\n");
print("$sizeVdsList\n");

$element = shift(@vgsList);
until( $element eq undef )
{
	printf("%f ",$element);
	$element = shift(@vgsList);
}
print("\n");

$element = shift(@vdsList);
until( $element eq undef )
{
	printf("%f ",$element);
	$element = shift(@vdsList);
}
print("\n");

# vgs = x axis
# vds = y axis
for($i=0;$i<$sizeVdsList;$i++)
{
	for($j=0;$j<$sizeVgsList;$j++)
	{
		printf("%g ",-1*(shift(@idsList)));
	}
	print("\n");
}
