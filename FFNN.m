clear
clc
%% PROCESSING
PiyasaTakasTotal = xlsread('PTF-03092018-17092019.xls');
ToplamTuketim=xlsread('GercekZamanliTuketim-03092018-17092019.xls');
YukTahminPlani = xlsread('YukTahminPlani-03092018-16092019.xls');

NGPrices=xlsread('Gaz_Referans_Fiyati_(GRF)_2018-09-03_2019-09-17.xls');
kk=1;
for i=1:1:length(NGPrices(:,2))
    NGPricesnew(kk:kk+23,1)=NGPrices(i,2);
    kk=kk+24;
end

Haftasonusilinik22=xlsread('Haftasonusilinik22.xls');

CoalPrice =  xlsread('CoalPrice.xls');
kk=1;
for i=1:1:length(CoalPrice)
    CoalPricesnew(kk:kk+23,1)=CoalPrice(i,1);
    kk=kk+24;
end
OilPrices = xlsread('OilPrices.xls');
kk=1;
for i=1:1:length(OilPrices(:,2))
    OilPricesnew(kk:kk+23,1)=OilPrices(i,2);
    kk=kk+24;
end
%% TRAINING
k=1;
j=0;
for i=24:24:336
    PiyasaTakasTotal1(:,k)=PiyasaTakasTotal(8761+j:8760+i,1);
    k=k+1;
    j=j+24;
end
%% PRICE OUTLIER
X1newSeptember = xlsread('X1newSeptember.xls');
X1new = xlsread('X1new');

satir_baslangic=1;
sutun=[1:9 10:21];
%% VALIDATION
j=24; rr=1; l=0; rrr=1; rrrr=1;
for i = 0:24:312
    net = feedforwardnet(1,'trainscg');
    [net,tr] = train(net,(X1new(satir_baslangic:8760+i,sutun))',(PiyasaTakasTotal(satir_baslangic:8760+i))');
    Feedforwardtestingpredict(:,rr)=sim(net,(X1newSeptember((i+1):(i+24),sutun))')';
    RealSeptember(:,rr) = PiyasaTakasTotal1(:,rr);
    error=0;
    for k = 1:1:24
        error=error+abs((RealSeptember(k,rr)-Feedforwardtestingpredict(k,rr))/RealSeptember(k,rr));
    end
    errorr(rr)=(error/24)*100;
    %
    error=0;
    for k = 1:1:24
        error=error+(RealSeptember(k,rrr)-Feedforwardtestingpredict(k,rrr))^2;
    end
    errorRMSE(rrr)=sqrt(error/24);
    %
    error=0;
    for k = 1:1:24
        error=error+(Feedforwardtestingpredict(k,rrrr)-RealSeptember(k,rrrr));
    end
    errorMSD(rrrr)=error/24;
    %
    Validation14days(l+1:l+24,1)=Feedforwardtestingpredict(:,rr);
    PlotActualPrices(l+1:l+24,1)=RealSeptember(:,rr);
    l=l+24;
    rr=rr+1;
    rrr=rrr+1;
    rrrr=rrrr+1;
end


% Testing Data
Errors = PlotActualPrices(1:288) - Validation14days(1:288);
MSE = mean(Errors.^2);
MAPE = mean(abs(Errors./PlotActualPrices(1:288)));
RMSE = sqrt(MSE);
ErrorMean = mean(Errors);
ErrorStd = std(Errors);
FFNNLowerBounds = [Validation14days-RMSE]
FFNNUpperBounds = [Validation14days+RMSE]

figure
x = 1:1:288';
subplot(2,2,[1 2]);
aa0 = plot(PlotActualPrices(1:288));
hold on;
aa1 = plot(Validation14days(1:288));
h1 = fill([x,fliplr(x)],[FFNNUpperBounds(1:288)',fliplr(FFNNLowerBounds(1:288)')],[0.8500 0.3250 0.0980],'facealpha',0.3);
hold off;
aa0.LineWidth = 2; aa1.LineWidth = 2;
legend('Actual Price','Predicted Price','Confidence Interval (%99)','Location','southwest','NumColumns',3);
ylabel('Price(TL/MWh)','FontSize',9);
xlabel('Hours','FontSize',9);
grid on;
title('Regular Forecast-Feedforward Neural Network','FontSize',11);


subplot(2,2,3);
plot(Errors);
title(['MSE = ' num2str(MSE, '%.2f') ', RMSE = ' num2str(RMSE, '%.2f')]);
ylabel('Errors','FontSize',9);
grid on;

subplot(2,2,4);
histfit(Errors, 50);
title(['Error Mean = ' num2str(ErrorMean, '%.2f') ', Error StD = ' num2str(ErrorStd, '%.2f')]);
% saveas(gcf,'FFNNplot','epsc')
% saveas(gcf,'FFNNplot2','png')
print(gcf,'FFNNplot2.png','-dpng','-r1000');



% print(plotperform(tr),'RegularForecastValid.png','-dpng','-r1000');

%% PLOT
% oneyearandmonths(12); %1-14 arasında parantezin içine sayı yazılabilir
% ngoilcoaleurograph(12); %1-14 arasında parantezin içine sayı yazılabilir
% laggedpricesvslaggednonprices(12); %1-14 arasında parantezin içine sayı yazılabilir
% scatterchartpricevsnonprices(12); %1-14 arasında parantezin içine sayı yazılabilir
% scatterngoilcoaleurograph(12);
% scatteroneyearandmonths(12);
%% FORECASTED FORECAST PRICE
j=24; rr=1; t=1; e=1; x=1; y=1; l=0; rrr=1; rrrr=1;
for i = 0:24:312
    if i >=24
        X1newSeptember(i+1:i+24,10) = ForecastedValidationSeptember(:,e);
        v1 = filter(ones(1,24)/24, 1, X1new(1:8760+i,10));
        X1newSeptember(i+1:i+24,12) = v1(8760+i-23:8760+i,1);
        e=e+1;
    end
    if i >=48
        X1new(8760+i-23:8760+i,10) = ForecastedValidationSeptember(:,t);
        v2 = filter(ones(1,24)/24, 1, X1new(1:8760+i,10));
        X1new(8760+i-23:8760+i,12) = v2(8760+i-23:8760+i,1);
        t=t+1;
    end
    if i>=168
        X1new(8760+i-23:8760+i,11) = ForecastedValidationSeptember(:,x);
        x=x+1;
    end
    if i>=192
        X1newSeptember((i+1):(i+24),11) = ForecastedValidationSeptember(:,y);
        y=y+1;
    end
    net = feedforwardnet(1,'trainscg');
    [net,tr] = train(net,(X1new(satir_baslangic:8760+i,sutun))',(PiyasaTakasTotal(satir_baslangic:8760+i))');
    ForecastedValidationSeptember(:,rr)=sim(net,(X1newSeptember((i+1):(i+24),sutun))')';
    RealSeptember(:,rr) = PiyasaTakasTotal1(:,rr);
    error=0;
    for k = 1:1:24
        error=error+abs((RealSeptember(k,rr)-ForecastedValidationSeptember(k,rr))/RealSeptember(k,rr));
    end
    Forecasterror(rr)=(error/24)*100;
    %
    error=0;
    for k = 1:1:24
        error=error+(RealSeptember(k,rrr)-ForecastedValidationSeptember(k,rrr))^2;
    end
    ForecasterrorRMSE(rrr)=sqrt(error/24);
    %
    error=0;
    for k = 1:1:24
        error=error+(ForecastedValidationSeptember(k,rrrr)-RealSeptember(k,rrrr));
    end
    ForecasterrorMSD(rrrr)=error/24;
    %
    Forecastforecast14days(l+1:l+24,1)=ForecastedValidationSeptember(:,rr);
    PlotActualPricesDynamicForecast(l+1:l+24,1)=RealSeptember(:,rr);
    l=l+24;
    PiyasaTakasTotal(8761+i:8760+i+24,1)=ForecastedValidationSeptember(:,rr);%Regressiondaki PTF yi değiştiriyor.
    rr=rr+1;
    rrr=rrr+1;
    rrrr=rrrr+1;
end


% Testing Data
Errors = PlotActualPrices(1:288) - Forecastforecast14days(1:288);
MSE = mean(Errors.^2);
MAPE = mean(abs(Errors./PlotActualPrices(1:288)));
RMSE = sqrt(MSE);
ErrorMean = mean(Errors);
ErrorStd = std(Errors);
DFFNNLowerBounds = [Forecastforecast14days-RMSE]
DFFNNUpperBounds = [Forecastforecast14days+RMSE]

figure
x = 1:1:288';
subplot(2,2,[1 2]);
aa2 = plot(PlotActualPrices(1:288));
hold on;
aa3 = plot(Forecastforecast14days(1:288));
h2 = fill([x,fliplr(x)],[DFFNNUpperBounds(1:288)',fliplr(DFFNNLowerBounds(1:288)')],[0.8500 0.3250 0.0980],'facealpha',0.3);
hold off;
aa2.LineWidth = 2; aa3.LineWidth = 2;
legend('Actual Price','Dynamic Forecast Price','Confidence Interval (%99)','Location','southwest','NumColumns',3);
ylabel('Price(TL/MWh)','FontSize',9);
xlabel('Hours','FontSize',9);
grid on;
title('Dynamic Forecast-Feedforward Neural Network','FontSize',11);


subplot(2,2,3);
plot(Errors);
title(['MSE = ' num2str(MSE, '%.2f') ', RMSE = ' num2str(RMSE, '%.2f')]);
ylabel('Errors','FontSize',9);
grid on;

subplot(2,2,4);
histfit(Errors, 50);
title(['Error Mean = ' num2str(ErrorMean, '%.2f') ', Error StD = ' num2str(ErrorStd, '%.2f')]);
% saveas(gcf,'FFNNDynamicplot','epsc','-depsc')
% saveas(gcf,'FFNNDynamicplot2','png','-r500')
print(gcf,'FFNNDynamicplot2.png','-dpng','-r1000');


% print(plotperform(tr),'DynamicForecastValid.png','-dpng','-r1000');

%% PLOT
% figure
% mm=plot(Validation14days(1:144),'b');
% hold on
% ll=plot(Forecastforecast14days(1:144),'g');
% hold on
% q=reshape(RealSeptember,24*14,1);
% nn=plot(q(1:144),'r');
% mm.LineWidth = 2;
% nn.LineWidth = 2;
% ll.LineWidth = 2;
% legend({'Regular Forecast Price','Dynamic Forecast Price','Real Price'},'FontSize',12,'Location','northwest','NumColumns',2);
% title('Estimated Prices vs Real Price','FontSize',12);
% xlabel('Hours','FontSize',12);
% ylabel('Price (TL/MWh)','FontSize',12);
% ylim([50 500])
