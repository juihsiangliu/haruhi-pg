function get_error(monteName,quadName)
	dataMonte = dlmread(monteName);
	dataQuad = dlmread(quadName);
	
	meanMonteList = mean(dataMonte,1) + 10^-100;
	varMonteList = var(dataMonte,0,1) + 10^-100;

	meanQuadList = dataQuad(1,:);
	varQuadList = dataQuad(2,:);

	beginIndex = 2;
	endIndex = length(meanMonteList);

	errMean = zeros(1,endIndex-beginIndex+1);
	errVar = zeros(1,endIndex-beginIndex+1);
	for(i=beginIndex:endIndex)
		errMean(i) = abs(meanMonteList(i) - meanQuadList(i))/meanMonteList(i);
		errVar(i) = abs(varMonteList(i) - varQuadList(i));
	end

%	errMean = abs(meanMonteList - meanQuadList)./meanMonteList;
%	errVar = abs(varMonteList - varQuadList);

	mean(errMean)
	mean(errVar)

end
