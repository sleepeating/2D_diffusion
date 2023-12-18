%% 清屏

clc;
clear;
close all; 

%% 参数设置

time_zero = 2020;   % 时间零点t0，可以自行设置也可以在line(32,34)读入时设置

spacing = 0.0352360817477097;                                               % 1 pixel = 0.0352360817477097 μm

filebox = 'C:\Users\eatsleep\Desktop\20231107\Si';
file = search_folder(filebox,'dat');
n_max = length(file);

%% 给出最强点的坐标和需要删除的边框，根据得到的结果略微调整

xl = 3;   xr = 55; yd = 44; yu = 192; % xleft,xright,ydown,yup;
xc = (xl+xr)/2; yc = (yd+yu)/2; % xcenter,ycenter;
px = 5; py = 6; msdt =1; msdp = 41;

msd = zeros(n_max,1);                                                       % 均方位移σ^2
delay_time = zeros(n_max,1);                                                % 时间延时t
delta_R = zeros(n_max,1);
temp_guass = cell(n_max,1);                                                 % 预分配一个元胞数组拟合高斯函数

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
%     if (ii == 1)
        delta_OD = importdata(temp_file);
        delta_OD0 = delta_OD;
        [xlength,ylength] = size(delta_OD);                                     % x为竖轴(横坐标)，y为横轴(纵坐标)
%     else
%         if (ii == 4)
%             delta_OD = importdata(temp_file)+delta_OD0;
%         else
%             delta_OD = importdata(temp_file)-delta_OD0;
%         end
%     end
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
    
    %% 动画制作可能使用的函数    
%     pause(0.1);

    %% 导入需要进行guass拟合的函数

    sum(delta_OD(xc-2:xc+2,yd:yu),1);
    temp_guass{ii} = delta_OD(xc,yd:yu);
    delay_time(ii) = decay;                                                 % 读取对应时间坐标

end

print([filebox,'\2Dplot.tif'], '-dtiffn','-r600');

%% 创建高斯分布函数

int_guass = cell2mat(temp_guass);                                           % 转换为可分析数组
y_position = spacing*(yd-yc):spacing:spacing*(yu-yc);                       % 创建横坐标

%% 进行高斯拟合并画图

figure
for ii = 1:n_max

    %% guass拟合(使用origin里的guass）

    A = 0.001;  xc = 0; y0 = 0; w = 0.7;                                    % 初始参数，需要一定的调整
    ini_0 = [A xc y0 w];                                                    % 参数顺序
    
    guass_fun =@(p,x) ...
        p(3)+(p(1)./(p(4).*sqrt(pi./2))).*exp(-2.*((x-p(2))./p(4)).^2);     % 拟合函数为y = y0 + (A/w*sqrt(pi/2))*exp(-2((x-xc)/w)^2)
    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');% 非线性最小二乘法使用levenberg-marquardt算法
    para = lsqcurvefit(guass_fun,ini_0,y_position,int_guass(ii,:),[-10,-1,-0.005,0.1],[10,1,0.005,10],options);% lsqcurvefit(函数，参数初始值，横坐标，纵坐标）返回拟合参数值               
    
    %% 数据处理，画图并获取msd
    A = para(1); xc = para(2); y0 =para(3); w = para(4);                    % 获取对应值
    delta_R(ii) = delta_R(ii) - y0;
    int_guass(ii,:) = int_guass(ii,:)-y0; para(3) = 0;                      % 减去基线
    plot(y_position,int_guass(ii,:),'.');    hold on
    trueguass(ii,:) = guass_fun(para,y_position);
    plot(y_position,guass_fun(para,y_position),'-');    hold on
    
    msd(ii) = (w./2).^2;                                                    % 单位为μm^2
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
xlim([delay_time(1),delay_time(end)]);
title('MSD','FontName','Arial','FontSize',14);
text(10,1,['mobility=',num2str(mobilty),'cm^2V^-^1s^-^1']);
print([filebox,'\MSD.tif'], '-dtiffn','-r600');

%% 画图ΔR/R图

figure
plot(delay_time,delta_R);
title('ΔR/R','FontName','Arial','FontSize',16);
xlabel('Time delay (ps)','FontName','Arial','FontSize',12);             % 设置x轴标签内容和字体
ylabel('ΔR/R (OD)','FontName','Arial','FontSize',12);             % 设置y轴标签内容和字体
xlim([delay_time(1),delay_time(end)]);
print([filebox,'\deltaR.tif'], '-dtiffn','-r600');