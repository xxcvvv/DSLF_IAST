clear all;close all;
% Author:jim
% Last updated:20200105

%% import setting
n=6;        % number of data set
opts = spreadsheetImportOptions("NumVariables", n);
opts.Sheet = "Sheet1";
opts.DataRange = "A1:F80";
opts.VariableNames = ["X", "Y", "X1", "Y1", "X2", "Y2"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];

%% import data
datapath=' ';        % add datapath
cd(datapath);
files=dir('*.XLSX');

%% For parameters
    file_name=files.name;
    path=[datapath,'\',file_name];
    rdata=readtable(path, opts, "UseExcel", false);
    data=table2array(rdata);
    
    for i=1:n/2       
        i
        p1=data(:,2.*i-1); 
        q1=data(:,2.*i);    
        k1=~isnan(p1);
        k2=~isnan(q1);
        p=p1(k1);
        q=q1(k2);
        save tempdata.mat p q;     
        LB=[0 0 0 0 0 0];               % lower boundary 
        UB=[20 0.01 2 20 0.01 2];       % upper boundary
        ObjectiveFunction = @simple_fitness; 
        nvars = 6;                      
        rng default;                    
        [coeff,fval]=ga(ObjectiveFunction,nvars,...
                [],[],[],[],LB,UB,[]); 
        disp('coeff'); 
        disp(coeff); 
        disp('fval');
        disp(fval); 
        clear p1 q1 k1 k2
        pm(1,i)=coeff(1);q1=num2str(pm(1,i));
        pm(2,i)=coeff(2);b1=num2str(pm(2,i));
        pm(3,i)=coeff(3);n1=num2str(pm(3,i));
        pm(4,i)=coeff(4);q2=num2str(pm(4,i));
        pm(5,i)=coeff(5);b2=num2str(pm(5,i));
        pm(6,i)=coeff(6);n2=num2str(pm(6,i));
        yf=@(c,xx) (c(1).*c(2).*xx.^c(3)./(1+c(2).*xx.^c(3)))+(c(4).*c(5).*xx.^c(6)./(1+c(5).*xx.^c(6))); 
        yfit=yf(coeff,p); 
        R=corrcoef(q,yfit);  
        r=num2str(R(1,2));
        ymax=max(yfit);
        % plot
        figure(i);
        plot(p,q,'o');
        hold on;
        plot(p,yfit); 
        xlim([0 120]);
        ylim([0 ymax+0.2]);
        legend('Êµ²â','ÄâºÏ','Location','NorthWest');
        title(opts.Sheet);
        text(90,0.7*ymax,['R:',r]);
        text(90,0.6*ymax,['q1:',q1]);text(90,0.5*ymax,['b1:',b1]);
        text(90,0.4*ymax,['n1:',n1]);text(90,0.3*ymax,['q2:',q2]);
        text(90,0.2*ymax,['b2:',b2]);text(90,0.1*ymax,['n2:',n2]);
        saveas(i,[opts.Sheet,'-',num2str(i)],'png');
        clear p q fval yf yfit nvars i q1 q2 n1 n2 b1 b2 R r ymax coeff LB UB
        delete tempdata.mat
    end
clear opts data file_name path datapath ObjectiveFunction n

%% Objective Function
function y=simple_fitness(c)
    load tempdata.mat p q; 
    xx=p; 
    yt=q; 
    yf=(c(1).*c(2).*xx.^c(3)./(1+c(2).*xx.^c(3)))+(c(4).*c(5).*xx.^c(6)./(1+c(5).*xx.^c(6))); 
    y=sum(abs(yf-yt))/length(yt); 
end