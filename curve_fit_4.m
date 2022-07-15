clear
%%
clear 
%Fitting Routine ----- FY-SA-FZ ------ @IA = 0,1,2,3,4 Pure Side slip 
fit_curve = 1; %Set this to do the curve fitting 
for ii=1
    clc
    load('F_DATA.mat');
    global  V
    %The filtered data from previous run is imported
    %Coefficients files containing values is imported
    %Declaring coefficients as the global coefficients
    for i=1:length(V)
        eval(strcat('global'," ",V{i,1}));
        eval(strcat(V{i,1},'=',num2str(V{i,2}),';'));
    end
    clear i
    ind0 = find(data.run18.fmdata(:,3)==0); %To segregate all the index values for data points with IA = 0
    ind1 = find(data.run18.fmdata(:,3)==1);
    ind4 = find(data.run18.fmdata(:,3)==4);
    ind3 = find(data.run18.fmdata(:,3)==3);
    ind2 = find(data.run18.fmdata(:,3)==2);
    
    
    fz= data.run18.fmdata(:,4);
    ia = data.run18.fmdata(:,3);
    sa = data.run18.fmdata(:,2);
    k = data.run18.fmdata(:,1);
    fy = data.run18.fmdata(:,5);
    %This will goin the lsqcurve fit for curve fitting 
    INPUT = [fz,sa,k,ia];
    d = fy;
    A_OLD = [];
    
    for i=1:18
        s = V{i,1};
        eval(strcat('A_OLD = [A_OLD,',s,'];'));
    end
    clear fy sa fz ia k s A RESNORM
    if fit_curve ==1
        options = optimset('MaxFunEvals',20000,'MaxIter',20000,'Display','final','TolX',1e-7,'TolFun',1e-7);
        fig2 = figure('MenuBar','none','Name',['Pacejka_97 fy_fitting Results'],'Position',[2 2 1600 1180],'NumberTitle','off');
        for k=1:5
            [A,RESNORM(k),RESIDUAL,EXITFLAG]=lsqcurvefit('mf_fy_ps',A_OLD,INPUT,d,[],[],options);
            AA(:,k) = A;
            for n = 1:18
                subplot(6,3,n);
                bar([AA(n,:)],'group');
                title(['A(',num2str(n),')','=',V{n,1}],'FontSize',8);
            end
            %Updating the new coefficients
            for n= 1:18
                eval(['A_OLD(' num2str(n) ') = ' num2str(A(n)) ' -1*eps*rand;']); % bootstrap
            end
            %set(fig1,'Name',filename,'\tFree Rolling Lateral Force\t\tIteration: ',num2str(k),'\tRESNORM:,' num2str(RESNORM(k)));
            drawnow
        end
        clear   sa fy fz ia k RESIDUAL RESSNORM
        for i=1:18
            V{i,2} = A(i);
        end
        A_OLD = A;
    end
    
    %Now we plot the fitting curves for IA = 1 and IA = 4 side by side
    q =0;
    fy_fitted =[];
    fy_filtered = [];
    for l = [0 1 2 3 4]
        q = q+1;
        figure
        %subplot(1,2,q)
        eval(strcat('fz= data.run18.fmdata(ind',num2str(l),',4);'));
        eval(strcat('ia= data.run18.fmdata(ind',num2str(l),',3);'));
        eval(strcat('sa= data.run18.fmdata(ind',num2str(l),',2);'));
        eval(strcat('k= data.run18.fmdata(ind',num2str(l),',1);'));
        eval(strcat('fy= data.run18.fmdata(ind',num2str(l),',5);'));
        INPUT = [fz,sa,k,ia];
        
        pacefy = mf_fy_ps(A_OLD,INPUT);   %Evaluating The Pacejka Equation with the fitted coefficients 
        eval(strcat('slips = unique(data.run18.fmdata(ind',num2str(l),',2));'));
        nslips = length(slips);
        eval(strcat('loads  = mean(reshape(data.run18.fmdata(ind',num2str(l),',4),nslips,[]),1);'));
        nloads = length(loads);
        eval(strcat('fy0 = reshape(data.run18.fmdata(ind',num2str(l),',5),nslips,nloads);'));
        fnplt(csaps({slips loads},fy0,0.99));
        %colormap(white);
        alpha 0.3
        shading interp
        hold on 
        eval(strcat('plot3(data.run18.fmdata(ind',num2str(l),',2),data.run18.fmdata(ind',num2str(l),',4),data.run18.fmdata(ind',num2str(l),',5),''k.'');'));
        hold on 
        eval(strcat('plot3(data.run18.fmdata(ind',num2str(l),',2),data.run18.fmdata(ind',num2str(l),',4),pacefy,''ro'');'));
        xlabel('SA');
        ylabel('FZ');
        zlabel('FY');
        title(strcat('IA = ',num2str(l)));
        fy_fitted = [fy_fitted;pacefy];
        eval(strcat('fy_filtered = [fy_filtered;data.run18.fmdata(ind',num2str(l),',5)];'))
    end
    
    SEOE = [];    %Here SEOE refers to the Standard Error of Estimate
    %SEOE should approach 0, 0 means the best fit, value in between 0-2
    %can be considered as a good fit
    R2 = [];
    for i = 1:27:216
        SEOE = [SEOE;sqrt(sum((fy_filtered(i:i+26)-fy_fitted(i:i+26)).^2))./length(fy_fitted(i:i+26))];
        %R2 estimation
        sigma_ref = std(fy_filtered(i:i+26));
        sigma_fit = std(fy_fitted(i:i+26));
        mean_ref = mean(fy_filtered(i:i+26));
        mean_fit = mean(fy_fitted(i:i+26));
        num = sum((fy_filtered(i:i+26)-mean_ref*ones(27,1)).*(fy_fitted(i:i+26)-mean_fit*ones(27,1)));
        den = 26*sigma_ref*sigma_fit;
        R2 = [R2;(num/den)];
    end
    SEOE
    R2
    
    save('F_DATA.mat','V','data');
    %clear 
    
end
%%
%Fitting Routine------MZ-SA-FZ--------@IA = 0,1,2,3,4 Pure side slip
fit_curve = 0; %Set this to do the curve fitting 
iterations = 2; %Define the number of iterations to be performed, minimum 2 
for ii=1
    load('F_DATA.mat');
    global  V
    for i=1:length(V)
        eval(strcat('global'," ",V{i,1}));
        eval(strcat(V{i,1},'=',num2str(V{i,2}),';'));
    end
    clear i 
    %FIRST STAGE  
        q1 = [19:22 25:28 31:32 35:38 40:41];
        q2 = [23 29:30 33:34 39 42:43];
        %Finding the index of required data, here we trim the data beyound the
        %Slip Angle Range(-10,10) 
        ind = find( data.run18.fmdata(:,2)<10 & data.run18.fmdata(:,2)>-10);
        ind0 = find(data.run18.fmdata(:,3)==0 & data.run18.fmdata(:,2)<10 & data.run18.fmdata(:,2)>-10);
        ind2 = find(data.run18.fmdata(:,3)==2 & data.run18.fmdata(:,2)<10 & data.run18.fmdata(:,2)>-10);
        ind4 = find(data.run18.fmdata(:,3)==4 & data.run18.fmdata(:,2)<10 & data.run18.fmdata(:,2)>-10);
        ind3 = find(data.run18.fmdata(:,3)==3 & data.run18.fmdata(:,2)<10 & data.run18.fmdata(:,2)>-10);
        ind1 = find(data.run18.fmdata(:,3)==1 & data.run18.fmdata(:,2)<10 & data.run18.fmdata(:,2)>-10);

        %Preparing Data for Stage- 1 curve fitting 
        fz= data.run18.fmdata(ind0,4);
        ia = data.run18.fmdata(ind0,3);
        sa = data.run18.fmdata(ind0,2);
        k = data.run18.fmdata(ind0,1);
        mz = data.run18.fmdata(ind0,7);
        INPUT = [fz,sa,k,ia];
        d = mz;
        B_OLD = [];

        for i=q1
            s = V{i,1};
            eval(strcat('B_OLD = [B_OLD,',s,'];'));
        end

        clear mz sa fz ia k s
        clear B RESNORM
        if fit_curve==1
            options = optimset('MaxFunEvals',20000,'MaxIter',20000,'Display','final','TolX',1e-7,'TolFun',1e-7);
            fig2 = figure('MenuBar','none','Name',['Pacejka_97 MZ_fitting(STAGE-I)'],'Position',[2 2 1600 1180],'NumberTitle','off');
            for k=1:iterations
                [B,RESNORM(k),RESIDUAL,EXITFLAG,OUTPUT]=lsqcurvefit('mf_mz_ps1',B_OLD,INPUT,d,[],[],options);
                BB(:,k) = B;
                for n = 1:16
                    subplot(4,4,n);
                    bar([BB(n,:)],'group');
                    title(['B(',num2str(n),')','=',V{q1(n),1}],'FontSize',8);
                end
                %Updating the new coefficients
                for n= 1:16
                    eval(['B_OLD(' num2str(n) ') = ' num2str(B(n)) ' -1*eps*rand;']); % bootstrap
                end
                %set(fig1,'Name',filename,'\tFree Rolling Lateral Force\t\tIteration: ',num2str(k),'\tRESNORM:,' num2str(RESNORM(k)));
                drawnow
            end
            clear   sa fy fz ia k RESIDUAL RESSNORM
            for i=1:16
                V{q1(i),2} = B(i);
            end
            save('12PSI.mat','V','data');
        end
        clear V data
        load('F_DATA.mat');
    
    %SECOND STAGE
        fz= data.run18.fmdata(ind,4);
        ia = data.run18.fmdata(ind,3);
        sa = data.run18.fmdata(ind,2);
        k = data.run18.fmdata(ind,1);
        mz = data.run18.fmdata(ind,7);
        INPUT = [fz,sa,k,ia];
        d = mz;
        B_OLD = [];
        for i=q2
            s = V{i,1};
            eval(strcat('B_OLD = [B_OLD,',s,'];'));
        end
        clear mz sa fz ia k s BB RESNORM
        if fit_curve==1
            options = optimset('MaxFunEvals',20000,'MaxIter',20000,'Display','final','TolX',1e-7,'TolFun',1e-7);
            fig2 = figure('MenuBar','none','Name',['Pacejka_97 MZ(STAGE-II)'],'Position',[2 2 1600 1180],'NumberTitle','off');
            for k=1:iterations
                [B,RESNORM(k),RESIDUAL,EXITFLAG,OUTPUT]=lsqcurvefit('mf_mz_ps2',B_OLD,INPUT,d,[],[],options);
                BB(:,k) = B;
                for n = 1:8
                    subplot(4,2,n);
                    bar([BB(n,:)],'group');
                    title(['B(',num2str(n),')','=',V{q2(n),1}],'FontSize',8);
                end
                %Updating the new coefficients
                for n= 1:8
                    eval(['B_OLD(' num2str(n) ') = ' num2str(B(n)) ' -1*eps*rand;']); % bootstrap
                end
                %set(fig1,'Name',filename,'\tFree Rolling Lateral Force\t\tIteration: ',num2str(k),'\tRESNORM:,' num2str(RESNORM(k)));
                drawnow
            end
            clear   sa fy fz ia k RESIDUAL RESSNORM
            for i=1:8
                
                V{q2(i),2} = B(i);
            end
            B_OLD = B;
        end
        
        %Now we plot the fitting curves for IA = 1 and IA = 4 side by side
        figure
        q =0;
        mz_fitted = [];
        mz_filtered = [];
        for l = [0 4]
            q = q+1;
            %subplot(1,2,q)
            figure
            eval(strcat('fz= data.run18.fmdata(ind',num2str(l),',4);'));
            eval(strcat('ia= data.run18.fmdata(ind',num2str(l),',3);'));
            eval(strcat('sa= data.run18.fmdata(ind',num2str(l),',2);'));
            eval(strcat('k= data.run18.fmdata(ind',num2str(l),',1);'));
            eval(strcat('mz= data.run18.fmdata(ind',num2str(l),',7);'));
            INPUT = [fz,sa,k,ia];
            
            pacemz = mf_mz_ps2(B_OLD,INPUT);   %Evaluating The Pacejka Equation with the fitted coefficients
            eval(strcat('slips = unique(data.run18.fmdata(ind',num2str(l),',2));'));
            nslips = length(slips);
            eval(strcat('loads  = mean(reshape(data.run18.fmdata(ind',num2str(l),',4),nslips,[]),1);'));
            nloads = length(loads);
            eval(strcat('mz0 = reshape(data.run18.fmdata(ind',num2str(l),',7),nslips,nloads);'));
            fnplt(csaps({slips loads},mz0,0.99));
            alpha 0.3 
            shading interp
            %colormap(white);
            hold on
            eval(strcat('plot3(data.run18.fmdata(ind',num2str(l),',2),data.run18.fmdata(ind',num2str(l),',4),data.run18.fmdata(ind',num2str(l),',7),''k.'');'));
            hold on
            eval(strcat('plot3(data.run18.fmdata(ind',num2str(l),',2),data.run18.fmdata(ind',num2str(l),',4),pacemz,''ro'');'));
            xlabel('SA');
            ylabel('FZ');
            zlabel('MZ');
            title(strcat('IA = ',num2str(l)));
            mz_fitted = [mz_fitted;pacemz];
            eval(strcat('mz_filtered = [mz_filtered;data.run18.fmdata(ind',num2str(l),',7)];'));
        end
        SEOE = [];    %Here SEOE refers to the Standard Error of Estimate 
        %SEOE should approach 0, 0 means the best fit, value in between 0-2
        %can be considered as a good fit 
        R2 = [];
        for i = 1:19:190
            SEOE = [SEOE;sqrt(sum((mz_filtered(i:i+18)-mz_fitted(i:i+18)).^2))./length(mz_fitted(i:i+18))];
            %R2 estimation 
            sigma_ref = std(mz_filtered(i:i+18));
            sigma_fit = std(mz_fitted(i:i+18));
            mean_ref = mean(mz_filtered(i:i+18));
            mean_fit = mean(mz_fitted(i:i+18));
            num = sum((mz_filtered(i:i+18)-mean_ref*ones(19,1)).*(mz_fitted(i:i+18)-mean_fit*ones(19,1)));
            den = 18*sigma_ref*sigma_fit;
            R2 = [R2;(num/den)];
        end
        SEOE
        R2
        %{
        [slips_mesh, loads_mesh] = meshgrid(slips,loads);
        indf = find(error_per<20);
        es = mean(error_per(indf));
        es
        e = abs(error_per(1:(nslips*5)));
        e = reshape(e,nslips,nloads);
        figure
        surf(slips_mesh,loads_mesh,e');
        shading interp
        %}
        save('F_DATA.mat','V','data');
        clear
end
%%
%Fitting Routine ------MX-SA-FZ--------@IA = 0,1,2,3,4 Pure side slip
fit_curve = 0; %Set this to do the curve fitting 
for ii=1
    load('F_DATA.mat');
        global  V
    for i=1:length(V)
        eval(strcat('global'," ",V{i,1}));
        eval(strcat(V{i,1},'=',num2str(V{i,2}),';'));
    end
    clear i 
    ind0 = find(data.run18.fmdata(:,3)==0);
    ind1 = find(data.run18.fmdata(:,3)==1);
    ind4 = find(data.run18.fmdata(:,3)==4);
    ind2 = find(data.run18.fmdata(:,3)==2);
    ind3 = find(data.run18.fmdata(:,3)==3);
    fz =  data.run18.fmdata(:,4);
    ia = data.run18.fmdata(:,3);
    sa = data.run18.fmdata(:,2);
    k = data.run18.fmdata(:,1);
    mx = data.run18.fmdata(:,8);
    %This will go into the curve fitting routine
    INPUT = [fz,sa,k,ia];
    d = mx;
    C_OLD = [];
    q = [85:95];
    for i = q
        s = V{i,1};
        eval(strcat('C_OLD = [C_OLD,',s,'];'));
    end
    clear mx sa fz ia k s
    clear C RESNORM CC
    if fit_curve ==1
        options = optimset('MaxFunEvals',20000,'MaxIter',20000,'Display','final','TolX',1e-7,'TolFun',1e-7);
        fig4 = figure('MenuBar','none','Name',['Pacejka_97 fy_fitting Results'],'Position',[2 2 1600 1180],'NumberTitle','off');
        
        for k=1:5
            [C,RESNORM(k),RESIDUAL,EXITFLAG]=lsqcurvefit('mf_mx_pcs',C_OLD,INPUT,d,[],[],options);
            CC(:,k) = C;
            for n = 1:11
                subplot(4,3,n);
                bar([CC(n,:)],'group');
                title(['C(',num2str(n),')','=',V{q(n),1}],'FontSize',8);
            end
            %Updating the new coefficients
            for n= 1:11
                eval(['C_OLD(' num2str(n) ') = ' num2str(C(n)) ' -1*eps*rand;']); % bootstrap
            end
            %set(fig1,'Name',filename,'\tFree Rolling Lateral Force\t\tIteration: ',num2str(k),'\tRESNORM:,' num2str(RESNORM(k)));
            drawnow
        end
        clear   sa mx fz ia k RESIDUAL RESSNORM
        for i=1:11
            V{q(i),2} = C(i);
        end
        C_OLD = C;
    end
    %Now we plot the fitting curves for IA = 1 and IA = 4 side by side
    figure
    q =0;
    mx_fitted = [];
    mx_filtered = [];
    for l = [0 1 2 4]
        q = q+1;
        %subplot(1,2,q)
        figure
        eval(strcat('fz= data.run18.fmdata(ind',num2str(l),',4);'));
        eval(strcat('ia= data.run18.fmdata(ind',num2str(l),',3);'));
        eval(strcat('sa= data.run18.fmdata(ind',num2str(l),',2);'));
        eval(strcat('k= data.run18.fmdata(ind',num2str(l),',1);'));
        eval(strcat('mx= data.run18.fmdata(ind',num2str(l),',8);'));
        INPUT = [fz,sa,k,ia];
        
        pacemx = mf_mx_pcs(C_OLD,INPUT);   %Evaluating The Pacejka Equation with the fitted coefficients 
        eval(strcat('slips = unique(data.run18.fmdata(ind',num2str(l),',2));'));
        nslips = length(slips);
        eval(strcat('loads  = mean(reshape(data.run18.fmdata(ind',num2str(l),',4),nslips,[]),1);'));
        nloads = length(loads);
        eval(strcat('mx0 = reshape(data.run18.fmdata(ind',num2str(l),',8),nslips,nloads);'));
        fnplt(csaps({slips loads},mx0,0.99));
        %colormap(white);
        alpha 0.3
        shading interp
        hold on 
        eval(strcat('plot3(data.run18.fmdata(ind',num2str(l),',2),data.run18.fmdata(ind',num2str(l),',4),data.run18.fmdata(ind',num2str(l),',8),''k.'');'));
        hold on 
        eval(strcat('plot3(data.run18.fmdata(ind',num2str(l),',2),data.run18.fmdata(ind',num2str(l),',4),pacemx,''ro'');'));
        xlabel('SA');
        ylabel('FZ');
        zlabel('MX');
        title(strcat('IA = ',num2str(l)));
        mx_fitted = [mx_fitted;pacemx];
        eval(strcat('mx_filtered = [mx_filtered;data.run18.fmdata(ind',num2str(l),',8)];'));
    end
    SEOE = [];    %Here SEOE refers to the Standard Error of Estimate
    %SEOE should approach 0, 0 means the best fit, value in between 0-2
    %can be considered as a good fit
    R2 = [];
    for i = 1:27:540
        SEOE = [SEOE;sqrt(sum((mx_filtered(i:i+26)-mx_fitted(i:i+26)).^2))./length(mx_fitted(i:i+26))];
        %R2 estimation
        sigma_ref = std(mx_filtered(i:i+26));
        sigma_fit = std(mx_fitted(i:i+26));
        mean_ref = mean(mx_filtered(i:i+26));
        mean_fit = mean(mx_fitted(i:i+26));
        num = sum((mx_filtered(i:i+26)-mean_ref*ones(27,1)).*(mx_fitted(i:i+26)-mean_fit*ones(27,1)));
        den = 26*sigma_ref*sigma_fit;
        R2 = [R2;(num/den)];
    end
    SEOE
    R2
    save('F_DATA.mat','V','data');
      
end
%%
%Fitting Routine --------FX-SR-FZ---------@IA = 0,2,4 Pure Longitudinal
%Slip
fit_curve = 0; %Set this to do the curve fitting 
for ii=1
   
    load('F_DATA.mat');
    global V
    for i=1:length(V)
        eval(strcat('global'," ",V{i,1}));
        eval(strcat(V{i,1},'=',num2str(V{i,2}),';'));
    end
    clear i 
    q = 48:62;
    %Seperating the index for data points with specific IA value at SA= 0
    %deg (Since we are dealing only with Pure Longitudinal Slip)
    alp0 = find(data.run36.fmdata(:,2)==0);
    ind0 = find(data.run36.fmdata(:,3)==0 & data.run36.fmdata(:,2)==0);
    ind2 = find(data.run36.fmdata(:,3)==2 & data.run36.fmdata(:,2)==0);
    ind4 = find(data.run36.fmdata(:,3)==4 & data.run36.fmdata(:,2)==0);
    %Averaging the values of Fz for all the IA and SA
    ind_fz1 = find(data.run36.fmdata(:,4)<-850 & data.run36.fmdata(:,4)>-900);
    data.run36.fmdata(ind_fz1,4) = mean(data.run36.fmdata(ind_fz1,4));
    
    ind_fz2 = find(data.run36.fmdata(:,4)<-650 & data.run36.fmdata(:,4)>-685);
    data.run36.fmdata(ind_fz2,4) = mean(data.run36.fmdata(ind_fz2,4));
    
    ind_fz3 = find(data.run36.fmdata(:,4)<-1090 & data.run36.fmdata(:,4)>-1130);
    data.run36.fmdata(ind_fz3,4) = mean(data.run36.fmdata(ind_fz3,4));
    
    ind_fz4 = find(data.run36.fmdata(:,4)<-200 & data.run36.fmdata(:,4)>-230);
    data.run36.fmdata(ind_fz4,4) = mean(data.run36.fmdata(ind_fz4,4));
    
    
    fz = data.run36.fmdata(alp0,4);
    sa = data.run36.fmdata(alp0,2);
    k = data.run36.fmdata(alp0,1);
    ia = data.run36.fmdata(alp0,3);
    fx = data.run36.fmdata(alp0,6);
    INPUT = [fz,sa,k,ia];
    d = fx;
    A_OLD = [];
    for i=1:15
        s = V{q(i),1};
        eval(strcat('A_OLD = [A_OLD,',s,'];'));
    end
    clear fx sa fz ia k s A RESNORM AA
    if fit_curve==1
        options = optimset('MaxFunEvals',20000,'MaxIter',20000,'Display','final','TolX',1e-7,'TolFun',1e-7);
        fig2 = figure('MenuBar','none','Name',['Pacejka_97 fx_fitting Results'],'Position',[2 2 1600 1180],'NumberTitle','off');
        for k=1:5
            [A,RESNORM(k),RESIDUAL,EXITFLAG]=lsqcurvefit('mf_fx_ps',A_OLD,INPUT,d,[],[],options);
            AA(:,k) = A;
            for n = 1:15
                subplot(5,3,n);
                bar([AA(n,:)],'group');
                title(['A(',num2str(n),')','=',V{q(n),1}],'FontSize',8);
            end
            %Updating the new coefficients
            for n= 1:15
                eval(['A_OLD(' num2str(n) ') = ' num2str(A(n)) ' -1*eps*rand;']); % bootstrap
            end
            %set(fig1,'Name',filename,'\tFree Rolling Lateral Force\t\tIteration: ',num2str(k),'\tRESNORM:,' num2str(RESNORM(k)));
            drawnow
        end
        clear sa fx fz ia k RESIDUAL RESSNORM
        
        for i=1:15
            V{q(i),2} = A(i);
        end
        A_OLD =A;
    end

    %Now we plot the fitting curves for IA = 1 and IA = 4 side by side
    figure
    q =0;
    fx_fitted = [];
    fx_filtered = [];
    for l = [0 2 4]
        q = q+1;
        %subplot(1,2,q)
        figure
        eval(strcat('fz= data.run36.fmdata(ind',num2str(l),',4);'));
        eval(strcat('ia= data.run36.fmdata(ind',num2str(l),',3);'));
        eval(strcat('sa= data.run36.fmdata(ind',num2str(l),',2);'));
        eval(strcat('k= data.run36.fmdata(ind',num2str(l),',1);'));
        eval(strcat('fx= data.run36.fmdata(ind',num2str(l),',6);'));
        INPUT = [fz,sa,k,ia];
        
        pacefx = mf_fx_ps(A_OLD,INPUT);   %Evaluating The Pacejka Equation with the fitted coefficients 
        eval(strcat('slips = unique(data.run36.fmdata(ind',num2str(l),',1));'));
        nslips = length(slips);
        eval(strcat('loads  = mean(reshape(data.run36.fmdata(ind',num2str(l),',4),nslips,[]),1);'));
        nloads = length(loads);
        eval(strcat('fx0 = reshape(data.run36.fmdata(ind',num2str(l),',6),nslips,nloads);'));
        fnplt(csaps({slips loads},fx0,0.999999));
        alpha 0.3
        shading interp
        hold on 
        eval(strcat('plot3(data.run36.fmdata(ind',num2str(l),',1),data.run36.fmdata(ind',num2str(l),',4),data.run36.fmdata(ind',num2str(l),',6),''k.'');'));
        hold on 
        eval(strcat('plot3(data.run36.fmdata(ind',num2str(l),',1),data.run36.fmdata(ind',num2str(l),',4),pacefx,''ro'');'));
        xlabel('SR');
        ylabel('FZ');
        zlabel('FX');
        title(strcat('IA = ',num2str(l)));
        fx_fitted = [fx_fitted;pacefx];
        eval(strcat('fx_filtered = [fx_filtered;data.run36.fmdata(ind',num2str(l),',6)];'));
    end
    SEOE = [];    %Here SEOE refers to the Standard Error of Estimate
    %SEOE should approach 0, 0 means the best fit, value in between 0-2
    %can be considered as a good fit
    R2 = [];
    for i = 1:27:324
        SEOE = [SEOE;sqrt(sum((fx_filtered(i:i+26)-fx_fitted(i:i+26)).^2))./length(fx_fitted(i:i+26))];
        %R2 estimation
        sigma_ref = std(fx_filtered(i:i+26));
        sigma_fit = std(fx_fitted(i:i+26));
        mean_ref = mean(fx_filtered(i:i+26));
        mean_fit = mean(fx_fitted(i:i+26));
        num = sum((fx_filtered(i:i+26)-mean_ref*ones(27,1)).*(fx_fitted(i:i+26)-mean_fit*ones(27,1)));
        den = 26*sigma_ref*sigma_fit;
        R2 = [R2;(num/den)];
    end
    R2
    SEOE
    save('F_DATA.mat','V','data');
    %clear 
end
%%
clear
%Fitting Routine ---------FY-SR-SA----------@const FZ and constant IA =
%0,2,4
fit_curve = [0 0 0]; %Set this to do the curve fitting [STAGE 1 , STAGE 2 & STAGE 3 ] 
for ii=1
    load('F_DATA.mat');
    global V
    for i=1:length(V)
        eval(strcat('global'," ",V{i,1}));
        eval(strcat(V{i,1},'=',num2str(V{i,2}),';'));
    end
    clear i
    
    %For STAGE 1 we have to take data-set with IA = 0 and varying FZ the
    %plot of FY-SA-SR
    %Index values for IA = 0
  
    ind0 = find(data.run36.fmdata(:,3)==0 );  %data for Stage 2
    
    ind2 = find( data.run36.fmdata(:,4)<-870 & data.run36.fmdata(:,4)>-950);
    data.run36.fmdata(ind2,4) = mean(data.run36.fmdata(ind2,4));
    ind2 = find( data.run36.fmdata(:,4)<-1000 & data.run36.fmdata(:,4)>-1120);
    data.run36.fmdata(ind2,4) = mean(data.run36.fmdata(ind2,4));
    ind2 = find( data.run36.fmdata(:,4)<-200 & data.run36.fmdata(:,4)>-300);
    data.run36.fmdata(ind2,4) = mean(data.run36.fmdata(ind2,4));
   
    ind2 = find( data.run36.fmdata(:,4)<-650 & data.run36.fmdata(:,4)>-750); %Nominal Load 
    data.run36.fmdata(ind2,4)=-700;
    ind1 = find(data.run36.fmdata(:,3)==0 & data.run36.fmdata(:,4)<-650 & data.run36.fmdata(:,4)>-750); % data for STAGE 1
    loads = unique(data.run36.fmdata(:,4));
    for jj = [0 2 4]
        for w = 1:length(loads)
            %inf1 = find(data.run36.fmdata(:,3)==0 & data.run36.fmdata(:,4)==loads(i));
            eval(strcat('inf',num2str(w),num2str(jj),' = find(data.run36.fmdata(:,3)==',num2str(jj),' & data.run36.fmdata(:,4)==loads(',num2str(w),'));'));
        end
    end
    
    sa = data.run36.fmdata(ind1,2);
    k = data.run36.fmdata(ind1,1);
    fz = data.run36.fmdata(ind1,4);
    ia = data.run36.fmdata(ind1,3);
    fy = data.run36.fmdata(ind1,5);
    INPUT = [fz,sa,k,ia];
    d = fy;
    q1 = [70 71 72 74 76 77 79 82 83 84];
    q2 = [75 78 80];
    q3 = [73 81];
    
    R_OLD = [];
    for i=1:length(q1)
        s = V{q1(i),1};
        eval(strcat('R_OLD = [R_OLD,',s,'];'));
    end
    clear fy sa fz ia k  R RESNORM RR
    if fit_curve(1)==1
        options = optimset('MaxFunEvals',20000,'MaxIter',20000,'Display','final','TolX',1e-7,'TolFun',1e-7);
        fig2 = figure('MenuBar','none','Name',['Pacejka_97 FY_combined_slip'],'Position',[2 2 1600 1180],'NumberTitle','off');
        for k=1:10
            [R,RESNORM(k),RESIDUAL,EXITFLAG]=lsqcurvefit('mf_fy_cs1',R_OLD,INPUT,d,[],[],options);
            RR(:,k) = R;
            for n = 1:10
                subplot(3,4,n);
                bar([RR(n,:)],'group');
                title(['R(',num2str(n),')','=',V{q1(n),1}],'FontSize',8);
            end
            %Updating the new coefficients
            for n= 1:10
                eval(['R_OLD(' num2str(n) ') = ' num2str(R(n)) ' -1*eps*rand;']); % bootstrap
            end
            %set(fig1,'Name',filename,'\tFree Rolling Lateral Force\t\tIteration: ',num2str(k),'\tRESNORM:,' num2str(RESNORM(k)));
            drawnow
        end
        clear  RESIDUAL RESSNORM
        for i=1:10
            V{q1(i),2} = R(i);
        end
        R_OLD = R;
    end
    %}
    %STAGE 2
    sa = data.run36.fmdata(ind0,2);
    k = data.run36.fmdata(ind0,1);
    fz = data.run36.fmdata(ind0,4);
    ia = data.run36.fmdata(ind0,3);
    fy = data.run36.fmdata(ind0,5);
    INPUT = [fz,sa,k,ia];
    d = fy;
    clear fy sa fz ia k  R RESNORM RR
    R_OLD = [];
    for i=1:length(q2)
        s = V{q2(i),1};
        eval(strcat('R_OLD = [R_OLD,',s,'];'));
    end
    if fit_curve(2)==1
        options = optimset('MaxFunEvals',20000,'MaxIter',20000,'Display','final','TolX',1e-7,'TolFun',1e-7,'Algorithm','levenberg-marquardt');
        fig2 = figure('MenuBar','none','Name',['Pacejka_97 FY_combined_slip_STAGE_2'],'Position',[2 2 1600 1180],'NumberTitle','off');
        for k=1:10
            [R,RESNORM(k),RESIDUAL,EXITFLAG]=lsqcurvefit('mf_fy_cs2',R_OLD,INPUT,d,[],[],options);
            RR(:,k) = R;
            for n = 1:3
                subplot(1,3,n);
                bar([RR(n,:)],'group');
                title(['R(',num2str(n),')','=',V{q2(n),1}],'FontSize',8);
            end
            %Updating the new coefficients
            for n= 1:3
                eval(['R_OLD(' num2str(n) ') = ' num2str(R(n)) ' -1*eps*rand;']); % bootstrap
            end
            %set(fig1,'Name',filename,'\tFree Rolling Lateral Force\t\tIteration: ',num2str(k),'\tRESNORM:,' num2str(RESNORM(k)));
            drawnow
        end
        clear sa fx fz ia k RESIDUAL RESSNORM
        for i=1:3
            V{q2(i),2} = R(i);
        end
        R_OLD = R;
    end
    
    %STAGE 3
    sa = data.run36.fmdata(ind0,2);
    k = data.run36.fmdata(ind0,1);
    fz = data.run36.fmdata(ind0,4);
    ia = data.run36.fmdata(ind0,3);
    fy = data.run36.fmdata(ind0,5);
    INPUT = [fz,sa,k,ia];
    d = fy;
    clear fy sa fz ia k  R RESNORM RR
    R_OLD = [];
    for i=1:length(q3)
        s = V{q3(i),1};
        eval(strcat('R_OLD = [R_OLD,',s,'];'));
    end
    if fit_curve(3)==1
        options = optimset('MaxFunEvals',20000,'MaxIter',20000,'Display','final','TolX',1e-7,'TolFun',1e-7,'Algorithm','levenberg-marquardt');
        fig2 = figure('MenuBar','none','Name',['Pacejka_97 FY_combined_slip_STAGE_2'],'Position',[2 2 1600 1180],'NumberTitle','off');
        for k=1:10
            [R,RESNORM(k),RESIDUAL,EXITFLAG]=lsqcurvefit('mf_fy_cs3',R_OLD,INPUT,d,[],[],options);
            RR(:,k) = R;
            for n = 1:2
                subplot(1,2,n);
                bar([RR(n,:)],'group');
                title(['R(',num2str(n),')','=',V{q3(n),1}],'FontSize',8);
            end
            %Updating the new coefficients
            for n= 1:2
                eval(['R_OLD(' num2str(n) ') = ' num2str(R(n)) ' -1*eps*rand;']); % bootstrap
            end
            %set(fig1,'Name',filename,'\tFree Rolling Lateral Force\t\tIteration: ',num2str(k),'\tRESNORM:,' num2str(RESNORM(k)));
            drawnow
        end
        clear sa fx fz ia k RESIDUAL RESSNORM
        for i=1:2
            V{q3(i),2} = R(i);
        end
        R_OLD = R;
    end
    
    %Plotting the results
    figure(1)
    q = 0 ;
    for j = [0]
        for i = [ 3 ]
            q = q+1;
            %subplot(3,3,q);
            figure(q);
            eval(strcat('fz = data.run36.fmdata(inf',num2str(i),num2str(j),',4);'))
            eval(strcat('ia = data.run36.fmdata(inf',num2str(i),num2str(j),',3);'))
            eval(strcat('sa = data.run36.fmdata(inf',num2str(i),num2str(j),',2);'))
            eval(strcat('k =  data.run36.fmdata(inf',num2str(i),num2str(j),',1);'))
            INPUT = [fz,sa,k,ia];
            pacefy = mf_fy_cs3(R_OLD,INPUT);
            eval(strcat('slip_ratio = unique(data.run36.fmdata(inf',num2str(i),num2str(j),',1));'))
            nslip_ratio = length(slip_ratio);
            eval(strcat(' slip_angle = unique(data.run36.fmdata(inf',num2str(i),num2str(j),",2));"));
            slip_angle = sort(slip_angle,'descend');
            nslip_angle = length(slip_angle);
            eval(strcat('fy0 = reshape(data.run36.fmdata(inf',num2str(i),num2str(j),',5),nslip_ratio,nslip_angle);'));
            eval(strcat('fnplt(csaps({slip_ratio,slip_angle},fy0,0.999999));'));
            alpha 0.3
            shading interp
            hold on
            eval(strcat('plot3(data.run36.fmdata(inf',num2str(i),num2str(j),',1),data.run36.fmdata(inf',num2str(i),num2str(j),',2),data.run36.fmdata(inf',num2str(i),num2str(j),',5),''k.'');'))
            hold on
            eval(strcat(' plot3(data.run36.fmdata(inf',num2str(i),num2str(j),',1),data.run36.fmdata(inf',num2str(i),num2str(j),',2),pacefy,''ro'');'))
            xlabel('SR');
            ylabel('SA');
            zlabel('FY');
            title(strcat('FZ = ',num2str(-round(loads(i))),' & IA = ',num2str(j)));
        end
    end
    save('F_DATA.mat','V','data');
%clear
end
%%
%Fitting Routine -----------FX-SR-SA----------@ FZ and  IA = 0,2,4
fit_curve = [0 0]; %Set this to do the curve fitting [STAGE 1 and STAGE 2]
for ii=1
    load('F_DATA.mat');
    global V
    for i=1:length(V)
        eval(strcat('global'," ",V{i,1}));
        eval(strcat(V{i,1},'=',num2str(V{i,2}),';'));
    end
    clear i 
    %For STAGE 1 we have to take data-set with IA = 0 and varying FZ the
    %plot of FY-SA-SR
    %Index values for IA = 0 
    loads = unique(data.run36.fmdata(:,4));
    ind0 = find(data.run36.fmdata(:,3)==0 );
    indl = find(data.run36.fmdata(:,3)==0 & data.run36.fmdata(:,4)==loads(1));
    for jj = [0 2 4]
        for w = 1:length(loads)
            %inf1 = find(data.run36.fmdata(:,3)==0 & data.run36.fmdata(:,4)==loads(i));
            eval(strcat('inf',num2str(w),num2str(jj),' = find(data.run36.fmdata(:,3)==',num2str(jj),' & data.run36.fmdata(:,4)==loads(',num2str(w),'));'));
        end
    end
    
    sa = data.run36.fmdata(ind0,2);
    k = data.run36.fmdata(ind0,1);
    fz = data.run36.fmdata(ind0,4);
    ia = data.run36.fmdata(ind0,3);
    fx = data.run36.fmdata(ind0,6);
    INPUT = [fz,sa,k,ia];
    d = fx;
    q1 = [63:68];
    q2 = [69];
    R_OLD = [];
    for i=1:6
        s = V{q1(i),1};
        eval(strcat('R_OLD = [R_OLD,',s,'];'));
    end
    clear fy sa fz ia k  R RESNORM RR
    if fit_curve(1)==1
        options = optimset('MaxFunEvals',80000,'MaxIter',80000,'Display','final','TolX',1e-7,'TolFun',1e-7);
        fig2 = figure('MenuBar','none','Name',['Pacejka_97 FX_combined_slip'],'Position',[2 2 1600 1180],'NumberTitle','off');
        for k=1:5
            [R,RESNORM(k),RESIDUAL,EXITFLAG]=lsqcurvefit('mf_fx_cs1',R_OLD,INPUT,d,[],[],options);
            RR(:,k) = R;
            for n = 1:6
                subplot(5,3,n);
                bar([RR(n,:)],'group');
                title(['R(',num2str(n),')','=',V{q1(n),1}],'FontSize',8);
            end
            %Updating the new coefficients
            for n= 1:6
                eval(['R_OLD(' num2str(n) ') = ' num2str(R(n)) ' -1*eps*rand;']); % bootstrap
            end
            %set(fig1,'Name',filename,'\tFree Rolling Lateral Force\t\tIteration: ',num2str(k),'\tRESNORM:,' num2str(RESNORM(k)));
            drawnow
        end
        clear  RESIDUAL RESSNORM
        for i=1:6
            V{q1(i),2} = R(i);
        end
    end
    %}
    %STAGE 2 
    sa = data.run36.fmdata(:,2);
    k = data.run36.fmdata(:,1);
    fz = data.run36.fmdata(:,4);
    ia = data.run36.fmdata(:,3);
    fx = data.run36.fmdata(:,6);
    INPUT = [fz,sa,k,ia];
    d = fx;
    clear fy sa fz ia k  R RESNORM RR
    R_OLD = [];
    for i=1
        s = V{q2(i),1};
        eval(strcat('R_OLD = [R_OLD,',s,'];'));
    end
    if fit_curve(2)==1
        options = optimset('MaxFunEvals',20000,'MaxIter',20000,'Display','final','TolX',1e-7,'TolFun',1e-7);
        fig2 = figure('MenuBar','none','Name',['Pacejka_97 FX_combined_slip_STAGE_2'],'Position',[2 2 1600 1180],'NumberTitle','off');
        for k=1:5
            [R,RESNORM(k),RESIDUAL,EXITFLAG]=lsqcurvefit('mf_fx_cs2',R_OLD,INPUT,d,[],[],options);
            RR(:,k) = R;
            for n = 1
                %subplot(1,2,n);
                bar([RR(n,:)],'group');
                title(['R(',num2str(n),')','=',V{q2(n),1}],'FontSize',8);
            end
            %Updating the new coefficients
            for n= 1
                eval(['R_OLD(' num2str(n) ') = ' num2str(R(n)) ' -1*eps*rand;']); % bootstrap
            end
            %set(fig1,'Name',filename,'\tFree Rolling Lateral Force\t\tIteration: ',num2str(k),'\tRESNORM:,' num2str(RESNORM(k)));
            drawnow
        end
        clear sa fx fz ia k RESIDUAL RESSNORM
        for i=1
            V{q2(i),2} = R(i);
        end
        R_OLD = R;
    end
    
    %Plotting the results 
    q = 1 ;
    for j = [0 2 4]
        for i = [1 2 3 4]
            q = q+1;
            %subplot(3,3,q);
            figure(q);
            eval(strcat('fz = data.run36.fmdata(inf',num2str(i),num2str(j),',4);'))
            eval(strcat('ia = data.run36.fmdata(inf',num2str(i),num2str(j),',3);'))
            eval(strcat('sa = data.run36.fmdata(inf',num2str(i),num2str(j),',2);'))
            eval(strcat('k =  data.run36.fmdata(inf',num2str(i),num2str(j),',1);'))
            INPUT = [fz,sa,k,ia];
            pacefx = mf_fx_cs2(R_OLD,INPUT);
            eval(strcat('slip_ratio = unique(data.run36.fmdata(inf',num2str(i),num2str(j),',1));'))
            nslip_ratio = length(slip_ratio);
            eval(strcat(' slip_angle = unique(data.run36.fmdata(inf',num2str(i),num2str(j),",2));"));
            slip_angle = sort(slip_angle,'descend');
            nslip_angle = length(slip_angle);
            eval(strcat('fx0 = reshape(data.run36.fmdata(inf',num2str(i),num2str(j),',6),nslip_ratio,nslip_angle);'));
            eval(strcat('fnplt(csaps({slip_ratio,slip_angle},fx0,0.999999));'));
            colormap(white);
            hold on
            eval(strcat('plot3(data.run36.fmdata(inf',num2str(i),num2str(j),',1),data.run36.fmdata(inf',num2str(i),num2str(j),',2),data.run36.fmdata(inf',num2str(i),num2str(j),',6),''k.'');'))
            hold on
            eval(strcat(' plot3(data.run36.fmdata(inf',num2str(i),num2str(j),',1),data.run36.fmdata(inf',num2str(i),num2str(j),',2),pacefx,''ro'');'))
            xlabel('SR');
            ylabel('SA');
            zlabel('FX');
            title(strcat('FZ = ',num2str(-round(loads(i))),' & IA = ',num2str(j)));
        end
    end
    
    save('F_DATA.mat','V','data');
    clear 
end
%%
%Fitting Routine -----------------MZ-SA-SR--------@const FZ and IA= 0,2,4
fit_curve = [0 0];
for ii = 1
    load('F_DATA.mat');
    global V
    for i=1:length(V)
        eval(strcat('global'," ",V{i,1}));
        eval(strcat(V{i,1},'=',num2str(V{i,2}),';'));
    end
    clear i 
    %For STAGE 1 we have to take data-set with IA = 0 and varying FZ the
    %plot of FY-SA-SR
    %Index values for IA = 0 
    loads = unique(data.run36.fmdata(:,4));
    ind0 = find(data.run36.fmdata(:,3)==0 );
    indl = find(data.run36.fmdata(:,3)==0 & data.run36.fmdata(:,4)==loads(1));
    for jj = [0 2 4]
        for w = 1:length(loads)
            %inf1 = find(data.run36.fmdata(:,3)==0 & data.run36.fmdata(:,4)==loads(i));
            eval(strcat('inf',num2str(w),num2str(jj),' = find(data.run36.fmdata(:,3)==',num2str(jj),' & data.run36.fmdata(:,4)==loads(',num2str(w),'));'));
        end
    end
    
    sa = data.run36.fmdata(ind0,2);
    k = data.run36.fmdata(ind0,1);
    fz = data.run36.fmdata(ind0,4);
    ia = data.run36.fmdata(ind0,3);
    mz = data.run36.fmdata(ind0,7);
    INPUT = [fz,sa,k,ia];
    d = mz;
    q1 = [44:47];
    S_OLD = [];
    for i=1:4
        s = V{q1(i),1};
        eval(strcat('S_OLD = [S_OLD,',s,'];'));
    end
    clear mz sa fz ia k  S RESNORM SS
    if fit_curve(1)==1
        options = optimset('MaxFunEvals',80000,'MaxIter',80000,'Display','final','TolX',1e-7,'TolFun',1e-7);
        fig2 = figure('MenuBar','none','Name',['Pacejka_97 FX_combined_slip'],'Position',[2 2 1600 1180],'NumberTitle','off');
        for k=1:5
            [S,RESNORM(k),RESIDUAL,EXITFLAG]=lsqcurvefit('mf_mz_cs',S_OLD,INPUT,d,[],[],options);
            SS(:,k) = S;
            for n = 1:4
                subplot(2,2,n);
                bar([SS(n,:)],'group');
                title(['S(',num2str(n),')','=',V{q1(n),1}],'FontSize',8);
            end
            %Updating the new coefficients
            for n= 1:9
                eval(['S_OLD(' num2str(n) ') = ' num2str(S(n)) ' -1*eps*rand;']); % bootstrap
            end
            %set(fig1,'Name',filename,'\tFree Rolling Lateral Force\t\tIteration: ',num2str(k),'\tRESNORM:,' num2str(RESNORM(k)));
            drawnow
        end
        clear  RESIDUAL RESSNORM
        for i=1:4
            V{q1(i),2} = S(i);
        end
        S_OLD = S;
    end
    
    %Plotting the results 
    q = 1 ;
    for j = [0 2 4]
        for i = [1 2 3 4]
            q = q+1;
            %subplot(3,3,q);
            figure(q);
            eval(strcat('fz = data.run36.fmdata(inf',num2str(i),num2str(j),',4);'))
            eval(strcat('ia = data.run36.fmdata(inf',num2str(i),num2str(j),',3);'))
            eval(strcat('sa = data.run36.fmdata(inf',num2str(i),num2str(j),',2);'))
            eval(strcat('k =  data.run36.fmdata(inf',num2str(i),num2str(j),',1);'))
            INPUT = [fz,sa,k,ia];
            pacemz = mf_mz_cs(S_OLD,INPUT);
            eval(strcat('slip_ratio = unique(data.run36.fmdata(inf',num2str(i),num2str(j),',1));'))
            nslip_ratio = length(slip_ratio);
            eval(strcat(' slip_angle = unique(data.run36.fmdata(inf',num2str(i),num2str(j),",2));"));
            slip_angle = sort(slip_angle,'descend');
            nslip_angle = length(slip_angle);
            eval(strcat('mz0 = reshape(data.run36.fmdata(inf',num2str(i),num2str(j),',7),nslip_ratio,nslip_angle);'));
            eval(strcat('fnplt(csaps({slip_ratio,slip_angle},mz0,0.999999));'));
            colormap(white);
            hold on
            eval(strcat('plot3(data.run36.fmdata(inf',num2str(i),num2str(j),',1),data.run36.fmdata(inf',num2str(i),num2str(j),',2),data.run36.fmdata(inf',num2str(i),num2str(j),',7),''k.'');'))
            hold on
            eval(strcat(' plot3(data.run36.fmdata(inf',num2str(i),num2str(j),',1),data.run36.fmdata(inf',num2str(i),num2str(j),',2),pacemz,''ro'');'))
            xlabel('SR');
            ylabel('SA');
            zlabel('MZ');
            title(strcat('FZ = ',num2str(-round(loads(i))),' & IA = ',num2str(j)));
        end
    end
    
    save('F_DATA.mat','V','data');
    clear 
end
