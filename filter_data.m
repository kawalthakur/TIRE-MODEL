clear 
clc
%%
%IMPORTING THE DATA TO THE DATA-BASE FOR RUN-18 AND RUN-36
f2 = {'run37.'};
f = {'run36.' 'run18.'};
for i = 1:length(f)
    t = importdata(strcat('B1464',f{i},'dat'),'\t',3);
    names  = t.textdata(2);
    nchans = size(t.data,2);
    if i == 2
        t.data(1:9000,:)=[];
    end
    clear name
    for n=1:nchans
        [name, names]=strtok(names);
        eval(strcat('data.',f{i},upper(string(name)),'=t.data(:,',num2str(n),');'));
    end
     
end
clear names nchans name file_name path_name f i 

%%
%EXTRACTING THE PRESSURE SPECIFIC DATA
c = true;
names = t.textdata(2);
nchans = size(t.data,2);
clear t 
while c==true
    press = str2num(input('Enter the value of Pressure (12/10/14 PSI):- ','s'));
    if press == 12
        for n=1:nchans
            [name, names]=strtok(names);
            eval(strcat('data.run18.',upper(string(name)),'= data.run18.',upper(string(name)),'(450:24730)'));
            eval(strcat('data.run36.',upper(string(name)),'= data.run36.',upper(string(name)),'(1:19035)'));
        end
        c = false;
    elseif press ==10
         for n=1:nchans
            [name, names]=strtok(names);
            eval(strcat('data.run18.',upper(string(name)),'=','data.run18.',upper(string(name)),'(29390:51710)'));
            eval(strcat('data.run36.',upper(string(name)),'= data.run36.',upper(string(name)),'(18932:37900)'));
         end
        c = false;
    elseif press==14
         for n=1:nchans
            [name, names]=strtok(names);
           eval(strcat('data.run18.',upper(string(name)),'=','data.run18.',upper(string(name)),'(56400:78730)'));
           eval(strcat('data.run36.',upper(string(name)),'= data.run36.',upper(string(name)),'(37890:end)'));
         end
        c = false; 
    else 
        c = true;
    end
end

FZ0 = mean(data.run18.FZ);
R0 = 0.2286;
clear names nchans c
%%
%FINDING ZERO CROSSING IN RUN-18 FOR SWEEPS OF SLIP-ANGLE  
%THIS IS FOR A PRESSURE VALUE OF 12 PSI ONLY
m = 1:length(data.run18.SA);
sp = spline(m,data.run18.SA);
z = fnzeros(sp);
z = round(z(1,:));
z = z.*repmat([1 0 1 0 0 0 0 0],1,(length(z)/8));
z(z==0)=[];
data.run18.z = unique(z);

figure
subplot(3,1,1);
sp = plot(m,data.run18.SA,'r');
hold on 
line([0 m(end)],[0 0],'color','k');
hold on 
plot(data.run18.z,zeros(length(data.run18.z)),'bo');
hold on 
grid on 
%Evaluate upto line 29 to select the suitable range 
zc_u = plot(data.run18.z,zeros(length(data.run18.z),1),'bo');

ylabel('Slip Angle');
legend([sp,zc_u],'Test data','zeros crossings');
clear zc zc_u l sp

subplot(3,1,2);
plot(data.run18.FZ,'r');
grid on 
hold on
line([0 m(end)],[-1000 -1000],'color','k');
%plot(data.run18.z,-1000*ones(length(data.run18.z),1),'bo');
ylabel('FZ(N)');

subplot(3,1,3);
plot(data.run18.IA,'r');
grid on 
hold on 
line([1 m(end)],[0 0 ],'color','k');
%plot(data.run18.z,zeros(length(data.run18.z),1),'bo');
xlabel('Counter');
ylabel('IA (deg)');

%%
%FINDING ZERO CROSSING IN RUN 36 FOR SWEEPS OF SLIP-RATIO
%THIS IS FOR PRESSURE VALUE OF 12 PSI ONLY 
data.run36.SL = ((mean(data.run36.RE)./mean(data.run36.RL)).*(1+data.run36.SR))-1;
data.run36.SL = data.run36.SL-0.01505;
m = 1:length(data.run36.SL);
sp = spline(m,data.run36.SL);
z = fnzeros(sp);
z = round(z(1,:));
c = [];
z_f = [];
for i = 1:length(z)-1
    if abs(diff([z(i) z(i+1)]))<=5
        c = [c;z(i+1)];
        if i>1
            continue;
        end
    else
        z_f = [z_f;z(i+1)];
        continue;
    end
    
    z_f = [z_f;mean(c)];
    clear c
    c = [];
end
data.run36.z = z;
clear z_f z
%Calculating the maximum valued indices
inmax = [];
for i = 1:2:length(data.run36.z)-1
    inmax = [inmax;data.run36.z(i)-105;data.run36.z(i+1)+105];
    
end
data.run36.max = inmax;
%this is the discrpancy in 10 PSI
if press == 10
    data.run36.max(find(data.run36.max==9724))=[];
    data.run36.max(find(data.run36.max==11760))=[];
end
clear inmax 

figure
subplot(4,1,1);
x = plot(m,data.run36.SL,'r');
hold on 
line([0 m(end)],[0 0],'color','k');
hold on 
sd = plot(data.run36.max,zeros(length(data.run36.max))*max(data.run36.SL),'bo');
hold on 
%{
for i=1:length(data.run36.max)
    line([data.run36.max(i) data.run36.max(i) ],[0.2 -0.3],'color','g');
    hold on 
    
end
%}
grid on 
ylabel('Slip Ratio ');
legend('Test data','Maximum Values');

subplot(4,1,2);
plot(m,data.run36.FZ,'r');
hold on 
line([0 m(end)],[-750 -750],'color','k');
hold on 
plot(data.run36.max,ones(length(data.run36.max))*(-750),'bo');
ylabel('FZ');

subplot(4,1,3);
plot(m,data.run36.IA,'r');
hold on 
line([0 m(end)],[0 0],'color','k');
hold on 
plot(data.run36.max,zeros(length(data.run36.max))*(-750),'bo');
ylabel('IA');

subplot(4,1,4);
plot(m,data.run36.SA,'r');
hold on 
line([0 m(end)],[0 0],'color','k');
hold on 
plot(data.run36.max,zeros(length(data.run36.max))*(-750),'bo');
ylabel('SA');
xlabel('Counter');

clear z_f z sp m c i ind_max temp RANGE x sd 
%%
%FILTERING THE TEST DATA FOR RUN-18 PURE SIDE SLIP  
clear data.run18.fmdata
q = 0;
data.run18.SL = ((data.run18.RE ./data.run18.RL).*(1+data.run18.SR))-1;
for i = 7:2:length(data.run18.z)
    sa = data.run18.SA(data.run18.z(i):data.run18.z(i+1));
    fy = data.run18.FY(data.run18.z(i):data.run18.z(i+1));
    fz = data.run18.FZ(data.run18.z(i):data.run18.z(i+1));
    mz = data.run18.MZ(data.run18.z(i):data.run18.z(i+1));
    mx = data.run18.MX(data.run18.z(i):data.run18.z(i+1));
    ia = data.run18.IA(data.run18.z(i):data.run18.z(i+1));
    k = data.run18.SL(data.run18.z(i):data.run18.z(i+1));
    fx = data.run18.FX(data.run18.z(i):data.run18.z(i+1));
    
    
    fz_mean = mean(data.run18.FZ(data.run18.z(i):data.run18.z(i+1)))*ones(length(data.run18.SA(data.run18.z(i):data.run18.z(i+1))),1);
    fy = fz_mean.*fy./fz;
    mz = fz_mean.*mz./fz;
    mx = fz_mean.*mx./fz;
    fx = fz_mean.*fx./fz;
   
    sp_fy = csaps(sa,fy,.1);
    sp_mz = csaps(sa,mz,.1);
    sp_mx = csaps(sa,mx,.1);
    sp_fx = csaps(sa,fx,.1);
    %{
    if i ==7
        figure
        plot(sa,fy);
        hold on 
        fnplt(sp_fy,'b',2,[-13 13]);
    end
    %}
    for sl=floor(min(sa)):1:ceil(max(sa)) %This is the range for SA ie -13 deg to 13 deg
        q=q+1;
        data.run18.fmdata(q,1)=round(mean(k));
        data.run18.fmdata(q,2)=sl;
        data.run18.fmdata(q,3)=round(mean(ia));
        data.run18.fmdata(q,4)=mean(fz);
        data.run18.fmdata(q,5)=fnval(sp_fy,sl); %Fy value at SA 
        data.run18.fmdata(q,6)=fnval(sp_fx,sl);
        data.run18.fmdata(q,7)=fnval(sp_mz,sl);
        data.run18.fmdata(q,8)=fnval(sp_mx,sl);
    end    
   
end
clear sa fy fz mz mx ia k fz_mean i imn imx q sl tmp sp_fy sp_mx sp_mz 
%{
A = [PCY1 PDY1 PDY2 PDY3 PEY1 PEY2 PEY3 PEY4 PKY1 PKY2 PKY3 PHY1...
    PHY2 PHY3 PVY1 PVY2 PVY3 PVY4]; 

data.run18.slips = unique(data.run18.fmdata(:,1))';
data.run18.nslips = length(data.run18.slips);

inx0 = find(data.run18.fmdata(:,2)==0);
inx4 = find(data.run18.fmdata(:,2)==4);

data.run18.fmdata0 = data.run18.fmdata(inx0,:);
data.run18.loads = mean(reshape(data.run18.fmdata0(:,3),data.run18.nslips,[]),1);
data.run18.nloads = length(data.run18.loads);
data.run18.fy0 = reshape(data.run18.fmdata0(:,4),data.run18.nslips,data.run18.nloads);
data.run18.mz0 = reshape(data.run18.fmdata0(:,5),data.run18.nslips,data.run18.nloads);
data.run18.INPUT = [data.run18.fmdata(inx0,3),data.run18.fmdata(inx0,1),data.run18.fmdata(inx0,7),data.run18.fmdata(inx0,2)];

%CURVE FITTING SUBROUTINE
%curve_fit(data.run18.INPUT,data.run18.fmdata(inx0,4));
pacefy = MF52_Fy_fcn(A,data.run18.INPUT);

figure
fnplt(csaps({data.run18.slips,data.run18.loads},data.run18.fy0,0.99));
colormap(white);
hold on 
plot3(data.run18.fmdata(inx0,1),data.run18.fmdata(inx0,3),data.run18.fmdata(inx0,4),'b*');
hold on 
plot3(data.run18.fmdata(inx0,1),data.run18.fmdata(inx0,3),pacefy,'ro');
%Plotting the Pacejka fit model 
figure
fnplt(csaps({data.run18.slips,data.run18.loads},data.run18.mz0,0.99));
colormap(white);

%}

%%
%FILTERING THE TEST DATA RUN-36
clear data.run36.fmdata
q = 0;
vb = 0;
shift_data = 0;
color = ['r' 'g' 'b' 'y'];
for i=1:2:length(data.run36.max)-1
    k_sa = data.run36.SL(data.run36.max(i):data.run36.max(i+1));
    k = data.run36.SL(data.run36.max(i):data.run36.max(i+1));
    sa = data.run36.SA(data.run36.max(i):data.run36.max(i+1));
    fy = data.run36.FY(data.run36.max(i):data.run36.max(i+1));
    fz = data.run36.FZ(data.run36.max(i):data.run36.max(i+1));
    mz = data.run36.MZ(data.run36.max(i):data.run36.max(i+1));
    mx = data.run36.MX(data.run36.max(i):data.run36.max(i+1));
    ia = data.run36.IA(data.run36.max(i):data.run36.max(i+1));
    fx = data.run36.FX(data.run36.max(i):data.run36.max(i+1));
    fz_mean = mean(data.run36.FZ(data.run36.max(i):data.run36.max(i+1)));
    %*ones(length(data.run36.SA(data.run36.max(i):data.run36.max(i+1))),1);
    fz_mean = round(abs(fz_mean));
    ia_mean = round(mean(ia));
    sa_mean = round(mean(sa));
    
    %Here we solve the discrepency in the FX vs SR plot 
    %This will shift the FX plot to FX= 0 at SR = 0 
  
    for iii = 1
        if fz_mean > 870 && fz_mean < 950
            if ia_mean == 0 && sa_mean==0
                k = k-0.026359;
            end
            if ia_mean == 2 && sa_mean==0
                k = k-0.027536;
            end
            if ia_mean == 4 && sa_mean==0
                k = k-0.03037;
            end
            if ia_mean == 0 && sa_mean==-3
                k = k-0.019644;
            end
            if ia_mean == 2 && sa_mean==-3
                k = k-0.022609;
            end
            if ia_mean == 4 && sa_mean==-3
                k = k-0.0270;
            end
            if ia_mean == 0 && sa_mean==-6
                k = k-0.01415;
            end
            if ia_mean == 2 && sa_mean==-6
                k = k-0.018336;
            end
            if ia_mean == 4 && sa_mean==-6
                k = k-0.024353;
            end
              
        elseif fz_mean > 600 && fz_mean < 700
            if ia_mean == 0 && sa_mean==0
                k = k-0.03303;
            end
            if ia_mean == 2 && sa_mean==0
                k = k-0.034731;
            end
            if ia_mean == 4 && sa_mean==0
                k = k-0.0393;
            end
            if ia_mean == 0 && sa_mean==-3
                k = k-0.02867;
            end
            if ia_mean == 2 && sa_mean==-3
                k = k-0.031897;
            end
            if ia_mean == 4 && sa_mean==-3
                k = k-0.0386;
            end
            if ia_mean == 0 && sa_mean==-6
                k = k-0.026708;
            end
            if ia_mean == 2 && sa_mean==-6
                k = k-0.03019;
            end
            if ia_mean == 4 && sa_mean==-6
                k = k-0.036998;
            end

        elseif fz_mean > 1000 && fz_mean < 1120
            if ia_mean == 0 && sa_mean==0
                k = k-0.017246;
            end
            if ia_mean == 2 && sa_mean==0
                k = k-0.01851;
            end
            if ia_mean == 4 && sa_mean==0
                k = k-0.020036;
            end
            if ia_mean == 0 && sa_mean==-3
                k = k-0.010182;
            end
            if ia_mean == 2 && sa_mean==-3
                k = k-0.011839;
            end
            if ia_mean == 4 && sa_mean==-3
                k = k-0.01489;
            end
            if ia_mean == 0 && sa_mean==-6
                k = k-0.00005;
            end
            if ia_mean == 2 && sa_mean==-6
                k = k-0.003597;
            end
            if ia_mean == 4 && sa_mean==-6
                k = k-0.0052107;
            end

        elseif fz_mean > 200 && fz_mean < 300
            if ia_mean == 0 && sa_mean==0
                k = k-0.046112;
            end
            if ia_mean == 2 && sa_mean==0
                k = k-0.049818;
            end
            if ia_mean == 4 && sa_mean==0
                k = k-0.06303;
            end
            if ia_mean == 0 && sa_mean==-3
                k = k-0.0436;
            end
            if ia_mean == 2 && sa_mean==-3
                k = k-0.049426;
            end
            if ia_mean == 4 && sa_mean==-3
                k = k-0.0637;
            end
            if ia_mean == 0 && sa_mean==-6
                k = k-0.042536;
            end
            if ia_mean == 2 && sa_mean==-6
                k = k-0.048117;
            end
            if ia_mean == 4 && sa_mean==-6
                k = k-0.062638;
            end

        else
            k = k-0;
         end
    end
 %}
 %{
 
    for jjj = 1
        if (ia_mean>=0 && sa_mean == -6 && fz_mean > 200 && fz_mean < 300)
            k_sa = k_sa-0.03814;
        end
        if (ia_mean>=0 && sa_mean == -3 && fz_mean > 200 && fz_mean < 300)
            k_sa = k_sa - 0.03602;
        end
        if (ia_mean>=0 && sa_mean == -6 && fz_mean > 600 && fz_mean < 700)
            k_sa = k_sa-0.03052;
        end
        if (ia_mean>=0 && sa_mean == -3 && fz_mean > 600 && fz_mean < 700)
            k_sa = k_sa -0.01744;
        end
        if (ia_mean>=0 && sa_mean == -6 && fz_mean > 870 && fz_mean < 950)
            k_sa = k_sa-0.01483;
        end
        if (ia_mean>=0 && sa_mean == -3 && fz_mean > 870 && fz_mean < 950)
            k_sa = k_sa -0.005294;
        end
        if (ia_mean>=0 && sa_mean == -6 && fz_mean > 1000 && fz_mean < 1124)
            k_sa = k_sa+0.01308;
        end
        if (ia_mean>=0 && sa_mean == -3 && fz_mean > 1000 && fz_mean < 1124)
            k_sa = k_sa +0.01166;
        end
            
    end
  %} 
    fy = -fz_mean.*fy./fz;
    mz = -fz_mean.*mz./fz;
    mx = -fz_mean.*mx./fz;
    fx = -fz_mean.*fx./fz;
    
    sp_fx = csaps(k,fx,0.9999);
    sp_mx = csaps(k,mx,0.9999);
    sp_mz = csaps(k,mz,0.9999);
    sp_fy = csaps(k_sa,fy,0.9999);
 
    %{
    for ss = [0 -3 -6]
        if (ia_mean==2 && sa_mean == ss && fz_mean > 870 && fz_mean < 950)
            vb = vb+1;
            figure(9)
            grid on
            hold on
            plot(k_sa,fy,strcat(color(vb),'.'));
            hold on
            fnplt(sp_fy,color(vb),1,[-0.218 0.218]);
            xlabel('Slip- Ratio');
            ylabel('FY');
           
        for ki = linspace(-0.05,0.1,10000)
            if fnval(fnder(sp_fy),ki)>-0.1 && fnval(fnder(sp_fy),ki)<0.1
                hold on
                %plot(ki,0,'ro');
                %text(ki,0,strcat('v = ',num2str(ki)));
                %break;
            end
        end
            
        end
    end
   %} 
    for sl = linspace(-0.2,0.2,27)
        q = q+1;
        data.run36.fmdata(q,1) = sl;        %The slip-ratio values 
        data.run36.fmdata(q,2) = round(mean(sa));
        data.run36.fmdata(q,3) = round(mean(ia)); %Inclination Angles
        data.run36.fmdata(q,4) = mean(fz);  %Normal Force
        data.run36.fmdata(q,5) = fnval(sp_fy,sl);
        data.run36.fmdata(q,6) = fnval(sp_fx,sl);   %Evaluating the Filtered values of FX 
        data.run36.fmdata(q,7) = fnval(sp_mz,sl);
        data.run36.fmdata(q,8) = fnval(sp_mx,sl);
        
    end
    
end
clear k sa fy fz mz mx ia fx fz_mean sp_fx sl q 
%%

data.run36.slips = unique(data.run36.fmdata(:,1))';
data.run36.nslips = length(data.run36.slips);
inx0 = find(data.run36.fmdata(:,3)==0 & data.run36.fmdata(:,2)==0);
%inx1 = find(data.run36.fmdata(:,2)==0 & data.run36.fmdata(:,5)==-6);
data.run36.fmdata0 = data.run36.fmdata(inx0,:);
%data.run36.fmdata1 = data.run36.fmdata(inx1,:);
data.run36.loads = mean(reshape(data.run36.fmdata0(:,4),data.run36.nslips,[]),1);
data.run36.nloads = length(data.run36.loads);


%data.run36.loads1 = mean(reshape(data.run36.fmdata1(:,3),data.run36.nslips,[]),1);
%data.run36.nloads1 = length(data.run36.loads1);

data.run36.fx0 = reshape(data.run36.fmdata0(:,6),data.run36.nslips,data.run36.nloads);
%data.run36.fx1 = reshape(data.run36.fmdata1(:,4),data.run36.nslips,data.run36.nloads1);
figure
%fnplt(csaps({data.run36.slips,data.run36.loads},data.run36.fx0,0.99999));
%hold on 
%fnplt(csaps({data.run36.slips,data.run36.loads1},data.run36.fx1,0.99999));
%colormap(white);

%%
load('COEFF.mat','V');
save('F_DATA.mat','data','V');



