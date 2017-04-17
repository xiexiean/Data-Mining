clear ; close all; clc
s=[1 2 3 7 8 9 10 11 12 13 14 15 17 18 21 23 24 25 26 27 28]; %标量属性
n=[4 5 6 16 19 20 22];    %数值属性
name={'surgery','Age','Hospital Number ','rectal temperature','pulse ','respiratory rate','temperature of extremities','peripheral pulse ','mucous membranes','capillary refill time','pain','peristalsis','abdominal distension','nasogastric tube','nasogastric reflux','nasogastric reflux PH','rectal examination - feces','abdomen','packed cell volume','total protein','abdominocentesis appearance','abdomcentesis total protein','outcome','surgical lesion','type of lesion1','type of lesion2','type of lesion3','cp_data'};
[data]=xlsread('data'); %加载数据

%对标称属性，给出每个可能取值的频数
Data=cell(1,28);
for k=s      
    c=tabulate(data(:,k));  
    Data{k}=c;
end
A=Data{1,2}(:,2)~=0;
Data{1,2}=Data{1,2}(A,:); 
A=Data{1,3}(:,2)~=0;
Data{1,3}=Data{1,3}(A,:);

for k=s
    fprintf('%s:\n',name{k});
    disp(Data{1,k}(:,1:2));
end

%求最大值
max=max(data);     
%求最小值
min=min(data);      
%求均值
mean=nanmean(data);     
%求中位数
median=nanmedian(data);    
%求上四分位数
Q1=prctile(data,25);    
%求下四分位数
Q3=prctile(data,75);     
%缺失值
empty=sum(isnan(data)); 
%最大、最小、均值、中位数、四分位数及缺失值的个数
num=[max;min;mean;median;Q1;Q3;empty];  
%num(:,s)=[];

for k=n
    fprintf('%s:\n',name{k});
    fprintf('最大值:%4.4f\n最小值:%4.4f\n均值:%4.4f\n中位数:%4.4f\n上四分位数:%4.4f\n下四分位数:%4.4f\n缺失值:%4.4f\n',num(1,k),num(2,k),num(3,k),num(4,k),num(5,k),num(6,k),num(7,k));  
end


for k=n
    figure(1);
    subplot(1,2,1),hist(data(:,k));
    title([name{k},'：直方图']);
    subplot(1,2,2),qqplot(data(:,k));
    title([name{k},'：qq图']);
    hold on;
    pause;
    hold off;
end

    figure(1);
    subplot(2,4,1),qqplot(data(:,4));
    title([name{4},'：qq图']);
    subplot(2,4,2),qqplot(data(:,5));
    title([name{5},'：qq图']);
    subplot(2,4,3),qqplot(data(:,6));
    title([name{6},'：qq图']);   
    subplot(2,4,4),qqplot(data(:,16));
    title([name{16},'：qq图']);
    subplot(2,4,5),qqplot(data(:,19));
    title([name{19},'：qq图']);
    subplot(2,4,6),qqplot(data(:,20));
    title([name{20},'：qq图']);
    subplot(2,4,7),qqplot(data(:,22));
    title([name{22},'：qq图']);
    hold on;
    pause;
    hold off;

close(figure(gcf));
%绘制盒图
% for k=n
%     figure(2);
%     boxplot(data(:,k));
%     title([name{k} ,'：盒图']);
%     hold on;
%     pause;
%     hold off;
% end
    figure(2);
    subplot(2,4,1),boxplot(data(:,4));
    title([name{4},'：盒图']);
    subplot(2,4,2),boxplot(data(:,5));
    title([name{5},'：盒图']);
    subplot(2,4,3),boxplot(data(:,6));
    title([name{6},'：盒图']);   
    subplot(2,4,4),boxplot(data(:,16));
    title([name{16},'：盒图']);
    subplot(2,4,5),boxplot(data(:,19));
    title([name{19},'：盒图']);
    subplot(2,4,6),boxplot(data(:,20));
    title([name{20},'：盒图']);
    subplot(2,4,7),boxplot(data(:,22));
    title([name{22},'：盒图']);
    hold on;
    pause;
    hold off;
close(figure(gcf));

%剔除缺失值 
nan=data;
nan(any(isnan(nan)'),:)=[];
for k=n
    figure(3);
    subplot(2,2,1),hist(data(:,k));
    title([name{k} ,'：剔除之前的直方图']);
    subplot(2,2,2),hist(nan(:,k));
    title([name{k} ,'：剔除之后的直方图']);
    subplot(2,2,3),qqplot(data(:,k));
    title([name{k} ,'：剔除之前的qq图']);
    subplot(2,2,4),qqplot(nan(:,k));
    title([name{k} ,'：剔除之后的qq图']);     
    hold on;
    pause;
    hold off;
end
close(figure(gcf));
xlswrite('nonan.xls',nan,'sheet1','A1');
%最高频值 
hf=data;
for k=n
    h=mode(hf(:,k));
    hf1=hf(:,k);
    hf1(isnan(hf1)) = h;
    figure(4);
    subplot(2,2,1),hist(data(:,k));
    title([name{k} ,'：修改之前的直方图']);
    subplot(2,2,2),hist(hf1);
    title([name{k} ,'：修改之后的直方图']);
    subplot(2,2,3),qqplot(data(:,k));
    title([name{k} ,'：修改之前的qq图']);
    subplot(2,2,4),qqplot(hf1);
    title([name{k} ,'：修改之后的qq图']);     
    hold on;
    pause;
    hold off;
end
close(figure(gcf));
%相关性
analytic_mat=data;
[m, n1] = size(analytic_mat); 
ATTRIBUTE_L = 1;
ATTRIBUTE_H = 28;
standard_line = analytic_mat(11, ATTRIBUTE_L: ATTRIBUTE_H);
% 获得相关性矩阵，函数见下
%CORRELATION_MAT_ATTRIBUTE 计算属性之间的相关性矩阵。
% 就是行和列都是属性，value(i,j)是属性i和属性j的相关性。它是对称矩阵哟。
COR_SIZE = ATTRIBUTE_H - ATTRIBUTE_L + 1; % 相关性矩阵的大小
cor_mat = -ones(COR_SIZE, COR_SIZE); % 初始化相关性矩阵，由于要取最大相关性，初始为最小值（-1）
for i = ATTRIBUTE_L: ATTRIBUTE_H - 1
    for j = i + 1: ATTRIBUTE_H
        merge = [analytic_mat(:, i), analytic_mat(:, j)]; % 将待求相关系数的两列并起来
        [NaN_line, ~] = find(isnan(merge) == 1);
        merge(NaN_line, :) = []; % 删掉含有NaN的行以便正确求解相关系数
        cor_indx = i - ATTRIBUTE_L + 1;
        cor_indy = j - ATTRIBUTE_L + 1; % 相关性矩阵下标
        cor_mat(cor_indx, cor_indy) = corr(merge(:, 1), merge(:, 2)); % merge的两列即去除NaN的两属性，求相关系数
        cor_mat(cor_indy, cor_indx) = cor_mat(cor_indx, cor_indy); % 对称矩阵
    end
end
cor_mat(isnan(cor_mat)) = -1;
cor_size = size(cor_mat, 1); % 矩阵大小，正常情况下是方阵
for i = 1: m
    for j = ATTRIBUTE_L: ATTRIBUTE_H
        if(isnan(analytic_mat(i, j)))
            [~, index] = sort(cor_mat(j - ATTRIBUTE_L + 1, :));
            index_list = fliplr(index); % sort升序，fliplr翻转，变成降序，得到参考的属性优先度列表
            flag = 0; % 标识是否补全成功
            for k = 1: cor_size
                ref_attr = index_list(k); % 用于补全参考的属性
                if(~isnan(analytic_mat(i, ref_attr)))
                    analytic_mat(i, j) = standard_line(j - ATTRIBUTE_L + 1) / standard_line(ref_attr) * ...
                        analytic_mat(i, ref_attr + ATTRIBUTE_L - 1); % 按比例补全（这不是最好的方法）
                    flag = 1;
                    break;
                end
            end
            if(flag == 0)
                disp(['Insert fail at row ', num2str(i), ' col ', num2str(j)]);
                return ;
            end
        end
    end
end
hf=analytic_mat;
for k=n
    figure(5);
    subplot(2,2,1),hist(data(:,k));
    title([name{k} ,'：修改之前的直方图']);
    subplot(2,2,2),hist(hf(:,k));
    title([name{k} ,'：修改之后的直方图']);
    subplot(2,2,3),qqplot(data(:,k));
    title([name{k} ,'：修改之前的qq图']);
    subplot(2,2,4),qqplot(hf(:,k));
    title([name{k} ,'：修改之后的qq图']);     
    hold on;
    pause;
    hold off;
end
close(figure(gcf));
xlswrite('cor.xls',analytic_mat,'sheet1','A1');

%相似性
%CORRELATION_MAT_ATTRIBUTE 计算属性之间的相关性矩阵。
%   就是行和列都是属性，value(i,j)是属性i和属性j的相关性。它是对称矩阵哟。
analytic_mat=data;
COR_SIZE = ATTRIBUTE_H - ATTRIBUTE_L + 1; % 相关性矩阵的大小

cor_mat = -ones(COR_SIZE, COR_SIZE); % 初始化相关性矩阵，由于要取最大相关性，初始为最小值（-1）
for i = ATTRIBUTE_L: ATTRIBUTE_H - 1
    for j = i + 1: ATTRIBUTE_H
        merge = [analytic_mat(:, i), analytic_mat(:, j)]; % 将待求相关系数的两列并起来
        [NaN_line, ~] = find(isnan(merge) == 1);
        merge(NaN_line, :) = []; % 删掉含有NaN的行以便正确求解相关系数
        cor_indx = i - ATTRIBUTE_L + 1;
        cor_indy = j - ATTRIBUTE_L + 1; % 相关性矩阵下标
        cor_mat(cor_indx, cor_indy) = corr(merge(:, 1), merge(:, 2)); % merge的两列即去除NaN的两属性，求相关系数
        cor_mat(cor_indy, cor_indx) = cor_mat(cor_indx, cor_indy); % 对称矩阵
    end
end
cor_mat(isnan(cor_mat)) = -1;
%SIMILARITY_MAT_SAMPLE 计算各个样本间的相似性。与correlation_mat_attribute类似。
%   此处相似性实际上求了欧几里得距离，因此越小越相似。其余与上述函数类似，不做多余的注释。
SIM_SIZE = size(analytic_mat, 1); % 相似矩阵大小，与analytic_mat样本数一致
sim_mat = ones(SIM_SIZE, SIM_SIZE) * 999; % 初始化为最大距离
for i = 1: SIM_SIZE - 1
    for j = i + 1: SIM_SIZE
        merge = [analytic_mat(i, ATTRIBUTE_L: ATTRIBUTE_H)', ...
            analytic_mat(j, ATTRIBUTE_L: ATTRIBUTE_H)']; % 将两行样本转置合并为属性数x2的矩阵
        [NaN_line, ~] = find(isnan(merge) == 1);
        merge(NaN_line, :) = [];     
        sim_mat(i, j) = norm(merge(:, 1) - merge(:, 2)); % 两样本的欧几里得距离
        sim_mat(j, i) = sim_mat(i, j); % 对称矩阵
    end
end
        sim_size = size(sim_mat, 1); % 矩阵大小，正常情况下是方阵
        for i = 1: m
            for j = ATTRIBUTE_L: ATTRIBUTE_H
                if(isnan(analytic_mat(i, j)))
                    [~, index_list] = sort(sim_mat(i, :));
                    flag = 0; % 标识是否补全成功
                    for k = 1: sim_size
                        ref_samp = index_list(k); % 用于补全参考的属性
                        if(~isnan(analytic_mat(ref_samp, j)))
                            analytic_mat(i, j) = analytic_mat(ref_samp, j); % 原样填上，补全
                            flag = 1;
                            break;
                        end
                    end
                    if(flag == 0)
                        disp(['Insert fail at row ', num2str(i), ' col ', num2str(j)]);
                        return ;
                    end
                end
            end
        end
hf=analytic_mat;
for k=n
    figure(6);
    subplot(2,2,1),hist(data(:,k));
    title([name{k} ,'：修改之前的直方图']);
    subplot(2,2,2),hist(hf(:,k));
    title([name{k} ,'：修改之后的直方图']);
    subplot(2,2,3),qqplot(data(:,k));
    title([name{k} ,'：修改之前的qq图']);
    subplot(2,2,4),qqplot(hf(:,k));
    title([name{k} ,'：修改之后的qq图']);     
    hold on;
    pause;
    hold off;
end
close(figure(gcf));
xlswrite('sim.xls',analytic_mat,'sheet1','A1');





