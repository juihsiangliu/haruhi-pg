function plot_result(type,filename,nodeNum,endTime,y_min,y_max,outname)

	fprintf('plot_result(type,filename,nodeNum,endTime,y_min,y_max,outname)\n');
	data = dlmread(filename);
	[monteNum stepNum] = size(data);
	time = linspace(0,endTime,stepNum);
	zero = zeros(1,stepNum);

	switch(type)
		case {'quad'}
			meanList = data(1,:);
			varList = data(2,:); 
				figure
				plot(time,max(zero,meanList(1 : stepNum )), 'linewidth',2);
				hold all;
				plot(time,max(zero,meanList(1 : stepNum) + 3 * varList(1:stepNum).^0.5) ,'linewidth',2);
				plot(time,max(zero,meanList(1 : stepNum) - 3 * varList(1:stepNum).^0.5) ,'linewidth',2);
				title('Quadratic','fontsize',20);
				legend('\fontsize {16} mean','\fontsize{16} + 3 sigma','\fontsize{16} - 3 sigma');
				xlabel('time (sec)','fontsize',16);
				ylabel('voltage (v)','fontsize',16);
				axis([0 endTime y_min y_max]);
				postfix = strcat('node',num2str(nodeNum));
				postfix = strcat('.',postfix);
				outname = strcat(outname,postfix);
				saveas(gcf,strcat(outname,'.png'))
		case {'monte'}
			meanList = mean(data,1);
			varList = var(data,0,1);
				figure
				plot(time,max(zero,meanList(1 : stepNum )), 'linewidth',2);
				hold all;
				plot(time,max(zero,meanList(1 : stepNum) + 3 * varList(1:stepNum).^0.5) ,'linewidth',2);
				plot(time,max(zero,meanList(1 : stepNum) - 3 * varList(1:stepNum).^0.5) ,'linewidth',2);
				title('Monte Carlo','fontsize',20);
				legend('\fontsize {16} mean','\fontsize{16} + 3 sigma','\fontsize{16} - 3 sigma');
				xlabel('time (sec)','fontsize',16);
				ylabel('voltage (v)','fontsize',16);
				axis([0 endTime y_min y_max]);
				postfix = strcat('node',num2str(nodeNum));
				postfix = strcat('.',postfix);
				outname = strcat(outname,postfix);
				saveas(gcf,strcat(outname,'.png'))
		case {'monte-sample'}
				figure
				for(j=1:min(20,monteNum))
					plot(time,max(zero,data(j,1 : stepNum )));
					hold on
				end
				title('Monte Carlo','fontsize',20);
				xlabel('time (sec)','fontsize',16);
				ylabel('voltage (v)','fontsize',16);
				axis([0 endTime y_min y_max]);
				postfix = strcat('node',num2str(nodeNum));
				postfix = strcat('.',postfix);
				outname = strcat(outname,postfix);
				saveas(gcf,strcat(outname,'.png'))
	end
end





