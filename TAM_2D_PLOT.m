%% 清屏

clc;
clear;
close all; 

%% 参数设置

time_zero = 2015;   % 时间零点t0，可以自行设置也可以在line(32,34)读入时设置

spacing = 0.0352360817477097;                                               % 1 pixel = 0.0352360817477097 μm

filebox = ['C:\Users\eatsleep\Desktop\20231107\20231113\fe2o3-s2'];
file = search_folder(filebox,'dat');
n_max = length(file);

%% 给出最强点的坐标和需要删除的边框，根据得到的结果略微调整

xl = 3;   xr = 55; yd = 34; yu = 218; % xleft,xright,ydown,yup;
xc = (xl+xr)/2; yc = (yd+yu)/2; % xcenter,ycenter;
px = 5; py = 6; msdt =1.9; msdp = 101;

msd = zeros(n_max,1);                                                       % 均方位移σ^2
delay_time = zeros(n_max,1);                                                % 时间延时t
delta_R = zeros(n_max,1);
temp_x = cell(n_max,1);
temp_y = cell(n_max,1);
temp_guass = cell(n_max,1);                                                 % 预分配一个元胞数组拟合高斯函数
tempcell =cell(1,1);
%% 数据读入

figure
% pause(10);
for ii = 1:n_max
    %% 输入文件名，即延时时间
    
    temp_file = cell2mat(file(ii));
    st_position = regexp(temp_file,'\\');
    dot_position = regexp(temp_file,'\.');
    decay_raw = str2double(temp_file((st_position(end)+1):dot_position-1)); % 读入的时间为四位数
    
    %% 如果读取了时间零点取消注释
%     if (ii==1)
%         time_zero = decay_raw;
%     end
    
    %% 计算延时时间坐标
    decay = (decay_raw-time_zero)/15;                                       % 单位ps
    time_tit = [num2str(decay) 'ps'];                                       % 标题为延时时间
    
    %% 读入强度矩阵
    
    delta_OD = importdata(temp_file);
    [xlength,ylength] = size(delta_OD);                                     % x为竖轴(横坐标)，y为横轴(纵坐标)
    
    %% 如果没有多留左上角的一行一列，请注释掉
    if (xlength==65)
      delta_OD = delta_OD(2:xlength,2:ylength);
      xlength = xlength-1;
      ylength = ylength-1;
    end
    
    %% 给出intensity矩阵
    
    intensity = zeros(xlength*ylength,3); 
    for yy = 1:ylength
      for xx = 1:xlength              
          intensity((yy-1).*xlength+xx,1) = xx.*spacing;
          intensity((yy-1).*xlength+xx,2) = yy.*spacing;
          intensity((yy-1).*xlength+xx,3) = delta_OD(xx,yy);
      end
    end

    %% 画图
    intensity_max = max(intensity(xlength*yd:xlength*yu,3));
    intensity_min = min(intensity(xlength*yd:xlength*yu,3));
    
     xx_bar = reshape(intensity(:,1),xlength,ylength)-xc*spacing;            % x轴的矩阵，y固定，x随竖轴递增
     yy_bar = reshape(intensity(:,2),xlength,ylength)-yc*spacing;            % y轴的矩阵，x固定，y随横轴递减
%     xx_bar = reshape(intensity(:,1),xlength,ylength);            % x轴的矩阵，y固定，x随竖轴递增
%     yy_bar = reshape(intensity(:,2),xlength,ylength);            % y轴的矩阵，x固定，y随横轴递减
    
    %% 调整画图方式，表示进行修改
     subplot(px,py,ii);
    
    %% 在中心点附近寻找最值
    dc = 5;    % 在中心点14范围内寻找最值
    delta_R(ii) = delta_OD(xc,yc);
    sq_OD = delta_OD(xc-dc:xc+dc,yc-dc:yc+dc);
    c_pos = [xc-dc-1 yc-dc-1];
    maxr = -realmax; minr = realmax;
    for nn = 1:2*dc+1
        for mm = 1:2*dc+1
            if sq_OD(mm,nn)>maxr
                maxr = sq_OD(mm,nn);
                max_pos = [mm nn];
            end
            if sq_OD(mm,nn)<minr
                minr = sq_OD(mm,nn);
                min_pos = [mm nn];
            end
        end
    end
    if abs(minr)>abs(maxr)
        delta_R(ii) = minr;
        c_pos = min_pos+c_pos;
    else
        delta_R(ii) = maxr;
        c_pos = max_pos+c_pos;
    end

%     dc = dc*4;
%     sq_OD = delta_OD(c_pos(1)-dc:c_pos(1)+dc,c_pos(2)-dc:c_pos(2)+dc);
%     xx_c = xx_bar(c_pos(1)-dc:c_pos(1)+dc,c_pos(2)-dc:c_pos(2)+dc)-c_pos(1).*spacing;
%     yy_c = yy_bar(c_pos(1)-dc:c_pos(1)+dc,c_pos(2)-dc:c_pos(2)+dc)-c_pos(2).*spacing; 
%     pcolor(xx_c,yy_c,sq_OD);
    
    %% 设置一个时间延时集合
    tttime = [0,1,10,20,100,250,500];
    [temptt,idtt] = find(tttime == decay);
    idtt = 1;
    if (idtt)
%         subplot(1,length(tttime),idtt);
    pcolor(xx_bar(xl:xr,yd:yu),yy_bar(xl:xr,yd:yu),delta_OD(xl:xr,yd:yu));
    shading interp
    set(gcf,'Colormap',turbo);
    axis equal
    xlim([(xl-xc)*spacing (xr-xc)*spacing]);
    ylim([(yd-yc)*spacing (yu-yc)*spacing]);
%     xlim([-dc*spacing dc*spacing]);
%     ylim([-dc*spacing dc*spacing]);
    
%     xlabel('X position (μm)','FontName','Arial','FontSize',12);             % 设置x轴标签内容和字体
%     ylabel('Y position (μm)','FontName','Arial','FontSize',12);             % 设置y轴标签内容和字体
    title(time_tit,'FontName','Arial','FontSize',6);
%     set(gca, 'Fontname', 'Arial', 'Fontsize', 12);                          % 设置xy轴的字体类型和字号大小的
    
    set(gca,'CLim',[intensity_min intensity_max]);                          % 设置强度坐标
    
    end
    %% 动画制作可能使用的函数    
%     pause(0.1);

    %% 导入需要进行guass拟合的函数

    temp_guass{ii} = delta_OD(xl:xr,yd:yu);
    temp_x{ii} = xx_bar(xl:xr,yd:yu);
    temp_y{ii} = yy_bar(xl:xr,yd:yu);
    
    delay_time(ii) = decay;                                                 % 读取对应时间坐标

end

print([filebox,'\2Dplot.tif'], '-dtiffn','-r600');



%     int_guass = cell2mat(temp_guass);  
%     xxfit = cell2mat(temp_x);
%     yyfit = cell2mat(temp_y);
%% 进行高斯拟合并画图

% figure
for ii = 1:n_max
    
    tempcell{1} = temp_guass{ii};
    zfit = cell2mat(tempcell);
    tempcell{1} = temp_x{ii};
    xfit = cell2mat(tempcell);
    tempcell{1} = temp_y{ii};
    yfit = cell2mat(tempcell);
    
    zfit = reshape(zfit,[],1);
    xfit = reshape(xfit,[],1);
    yfit = reshape(yfit,[],1);


    %% guass拟合(使用origin里的guass）

    A = 0.001;  xc = 0; yc = 0; z0 = 0; w = 0.5;                                    % 初始参数，需要一定的调整
    Ac = 0.0001; k = 0; u = 0; lam = 0.593681576065333; b = 0 ;
%     ini_0 = [A xc yc z0 w Ac k u lam];                                                    % 参数顺序
%     ini_0 = [A xc yc z0 w]; 
    ini_0 = [A xc yc z0 w];
    guass_fun =@(p,q) ...
        p(4)+(p(1)./(p(5).*sqrt(pi./2))).*exp(-2.*((q(:,1)-p(2)).^2+(q(:,2)-p(3)).^2)./p(5).^2);
%         +p(6).*cos(p(7).*sqrt((q(:,1)-p(2)).^2+(q(:,2)-p(3)).^2)+p(8)).*cos(-p(9).*sqrt((q(:,1)-p(2)).^2+(q(:,2)-p(3)).^2));     % 拟合函数为y = y0 + (A/w*sqrt(pi/2))*exp(-2((x-xc)/w)^2)
%     options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');% 非线性最小二乘法使用levenberg-marquardt算法
    para = lsqcurvefit(guass_fun,ini_0,[xfit,yfit],zfit,[-20,-1,-1,-0.01,0.1],[20,1,1,0.01,4]);% lsqcurvefit(函数，参数初始值，横坐标，纵坐标）返回拟合参数值               
    
    %% 数据处理，画图并获取msd
    A = para(1); xc = para(2); yc = para(3); z0 = para(4); w = para(5);                    % 获取对应值
    delta_R(ii) = delta_R(ii) - z0;
    zfit = zfit-z0; para(4) = 0;                      % 减去基线
    zzplot = guass_fun(para,[xfit,yfit]);
    zzplot = reshape(zzplot,[xr-xl+1,yu-yd+1]);
    zfit = reshape(zfit,[xr-xl+1,yu-yd+1]);
    xfit = reshape(xfit,[xr-xl+1,yu-yd+1]);
    yfit = reshape(yfit,[xr-xl+1,yu-yd+1]);
    onefit(ii,:) = zfit((xl+xr)./2,:);
    oneplot(ii,:) = zzplot((xl+xr)./2,:);
%     figure
%     plot(yfit((xl+xr)./2,:),oneplot,'r-',yfit((xl+xr)./2,:),onefit,'b.');
    


    %% 导出减去基线后的二维图像数据，并添加横纵坐标
%     temp_file = cell2mat(file(ii));
%     dot_position = regexp(temp_file,'\.');
%     out_write = [0,yfit(1,:);xfit(:,1),zfit];
%     file_out = [temp_file(1:dot_position(end)),'csv'];
%     writematrix(out_write,file_out);

    %% 画拟合图区域

%     figure 
%     subplot(2,1,1); 
%     surf(xfit,yfit,500*zfit);
%     shading interp
%     set(gcf,'Colormap',turbo);
%     axis equal
%     subplot(2,1,2);
%     surf(xfit,yfit,500*zzplot);
%     shading interp
%     set(gcf,'Colormap',turbo);
%     axis equal

    msd(ii) = (w./2).^2;                                                    % 单位为μm^2
end

%% 画出x=0处的y坐标拟合图

figure

for ii = 1:n_max
    if (ii)
%     if (mod(ii,2) == 1)
% floor(n_max./2)+1
        subplot(5,6,ii);
        plot(yfit((xl+xr)./2,:),oneplot(ii,:),'-',yfit((xl+xr)./2,:),onefit(ii,:),'.');
        hold on
    end
end


%% 画出msd随时间的曲线

figure

%% 对start和stop位置进行选取

sr = 1; so = n_max;
for ii = 1:n_max
    if delay_time(ii)<msdt
        sr = ii;
    end
    if delay_time(ii)<msdp
        so = ii;
    end
end

%% 线性拟合y=k1*x+k2，画图

k_para = polyfit(delay_time(sr:so,1),msd(sr:so,1),1);
plot(delay_time,msd,'o');    hold on
plot(delay_time(sr:so,1),polyval(k_para,delay_time(sr:so,1)),'-');
% xlim([0 50]);

mobilty = k_para(1)*193400                                                  % 得到最后的迁移率
xlabel('Time delay (ps)','FontName','Arial','FontSize',12);             % 设置x轴标签内容和字体
ylabel('σ^2 (μm^2)','FontName','Arial','FontSize',12);             % 设置y轴标签内容和字体
% xlim([delay_time(1),delay_time(end)]);
title('MSD','FontName','Arial','FontSize',14);
text(10,1,['mobility=',num2str(mobilty),'cm^2V^-^1s^-^1']);
% text(0,0,['w=',num2str(w*1000),'um']);
print([filebox,'\MSD.tif'], '-dtiffn','-r600');

%% 画图ΔR/R图

figure
plot(delay_time,delta_R);
title('ΔR/R','FontName','Arial','FontSize',16);
xlabel('Time delay (ps)','FontName','Arial','FontSize',12);             % 设置x轴标签内容和字体
ylabel('ΔR/R (OD)','FontName','Arial','FontSize',12);             % 设置y轴标签内容和字体
% xlim([delay_time(1),delay_time(end)]);
print([filebox,'\deltaR.tif'], '-dtiffn','-r600');
% close all; 